"""
Consumption footprint with a REX LCI extension and EF 3.1
==========================================================
    footprint = C . F . diag(1/x) . (I - A)^-1 . y_country

C: EF 3.1 characterisation factors (categories x LCI flows), F: extension from
rex_env_extension.py / rex_calibrate.py, A = T diag(1/x), x = T 1 + Y 1.
The output vector (I - A)^-1 y is solved iteratively (GMRES) against T stored
once as a float32 memmap, so the 30807^2 Leontief inverse is never formed.
Covers industry emissions only (direct household emissions are not extended).

Step 1, characterisation factors (in the env of the flows' Brightway project):
    python rex_footprint.py cfs --project default --prefix "EF v3.1" --out cfs_ecoinvent.csv
    <werkbank python> rex_footprint.py cfs --project "bafu test" --prefix sentier --out cfs_bafu.csv
Step 2, footprint (AB env, needs h5py):
    python rex_footprint.py footprint --rex-dir <REX3_2022> --t-mat <T_REX3.mat> --cache <dir>
        --run <out dir> --cfs cfs_ecoinvent.csv --country Switzerland --label ecoinvent
"""

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

# sentier EF v3.1 names -> ecoinvent EF v3.1 names (used as common category names)
SENTIER_TO_EI = {
    "Acidification": "acidification",
    "Climate change": "climate change",
    "Climate change-Biogenic": "climate change: biogenic",
    "Climate change-Fossil": "climate change: fossil",
    "Climate change-Land use and land use change": "climate change: land use and land use change",
    "EF-particulate Matter": "particulate matter formation",
    "Eutrophication marine": "eutrophication: marine",
    "Eutrophication, freshwater": "eutrophication: freshwater",
    "Eutrophication, terrestrial": "eutrophication: terrestrial",
    "Human toxicity, cancer": "human toxicity: carcinogenic",
    "Human toxicity, cancer_organics": "human toxicity: carcinogenic, organics",
    "Human toxicity, cancer_inorganics": "human toxicity: carcinogenic, inorganics",
    "Human toxicity, non-cancer": "human toxicity: non-carcinogenic",
    "Human toxicity, non-cancer_organics": "human toxicity: non-carcinogenic, organics",
    "Human toxicity, non-cancer_inorganics": "human toxicity: non-carcinogenic, inorganics",
    "Ionising radiation, human health": "ionising radiation: human health",
    "Land use": "land use",
    "Ozone depletion": "ozone depletion",
    "Photochemical ozone formation - human health": "photochemical oxidant formation: human health",
    "Resource use, fossils": "energy resources: non-renewable",
    "Resource use, minerals and metals": "material resources: metals/minerals",
    "Water use": "water use",
    "Ecotoxicity, freshwater": "ecotoxicity: freshwater",
    "Ecotoxicity, freshwater_organics": "ecotoxicity: freshwater, organics",
    "Ecotoxicity, freshwater_inorganics": "ecotoxicity: freshwater, inorganics",
}


# ---------------------------------------------------------------------------
def export_cfs(project, prefix, out):
    """Write (database, code, category, unit, cf) for all methods starting with prefix."""
    import bw2data as bd
    bd.projects.set_current(project)
    rows = []
    id_to_key = None if hasattr(bd, "get_node") else {v: k for k, v in bd.mapping.items()}  # bw2data 3
    for m in bd.methods:
        if m[0] != prefix and not (prefix == "sentier" and m[0] == "sentier"):
            continue
        cat = m[1] if prefix != "sentier" else SENTIER_TO_EI[m[2]]
        unit = bd.methods[m].get("unit")
        for entry in bd.Method(m).load():
            flow, cf = entry[0], entry[1]
            if isinstance(flow, int):
                key = id_to_key[flow] if id_to_key is not None else bd.get_node(id=flow).key
            else:
                key = tuple(flow)
            rows.append((key[0], key[1], cat, unit, float(cf)))
    df = pd.DataFrame(rows, columns=["database", "code", "category", "unit", "cf"])
    df.to_csv(out, index=False)
    print(f"{df.category.nunique()} categories, {len(df)} CFs -> {out}")


# ---------------------------------------------------------------------------
def flow_keys(run):
    """(database, code) per row of F_REX, reconstructed from the cached pools."""
    run = Path(run)
    pools = sorted(run.glob("pools_*.json"))[0]
    flows = json.load(open(pools))
    E = np.load(pools.with_suffix(".npz"))["E"]
    keep = np.abs(E[:, :, -1]).sum(axis=1) > 0
    keys = [tuple(f["key"]) for f, k in zip(flows, keep) if k]
    names = pd.read_csv(run / "F_REX_flows.csv").name.tolist()
    assert [f["name"] for f, k in zip(flows, keep) if k] == names, "pools and F_REX_flows.csv differ"
    return keys


def t_memmap(t_mat, cache):
    """REX T as float32 memmap (T orientation: rows = supplying sectors), built once."""
    import h5py
    path = Path(cache) / "T_REX3_f32.npy"
    if path.exists():
        return np.load(path, mmap_mode="r")
    with h5py.File(t_mat) as f:
        d = f[next(k for k in f if not k.startswith("#"))]      # stored transposed (MATLAB)
        n = d.shape[0]
        T = np.lib.format.open_memmap(path, mode="w+", dtype=np.float32, shape=(n, n))
        step = 1000
        for k0 in range(0, n, step):                             # h5py column k = T row k
            T[k0:k0 + step] = d[:, k0:k0 + step].T.astype(np.float32)
        T.flush()
    return np.load(path, mmap_mode="r")


def matvec(T, w, block=2048):
    out = np.empty(T.shape[0])
    w32 = w.astype(np.float32)
    for r0 in range(0, T.shape[0], block):
        out[r0:r0 + block] = T[r0:r0 + block] @ w32
    return out


def row_sums(T, block=2048):
    return np.concatenate([np.asarray(T[r0:r0 + block], dtype=np.float64).sum(axis=1)
                           for r0 in range(0, T.shape[0], block)])


def solve_output(T, x, y, tol=1e-6):
    """x_y = (I - T diag(1/x))^-1 y by GMRES (float32 matvec limits precision to ~1e-7)."""
    from scipy.sparse.linalg import LinearOperator, gmres
    inv_x = np.where(x > 0, 1 / np.where(x > 0, x, 1), 0.0)
    n = len(x)
    op = LinearOperator((n, n), matvec=lambda v: v - matvec(T, np.ravel(v) * inv_x), dtype=float)
    it = [0]
    sol, info = gmres(op, y, x0=y.copy(), rtol=tol, restart=50, maxiter=20,
                      callback=lambda r: it.__setitem__(0, it[0] + 1), callback_type="pr_norm")
    resid = np.linalg.norm(op.matvec(sol) - y) / np.linalg.norm(y)
    print(f"GMRES: {it[0]} iterations, relative residual {resid:.1e}")
    if info != 0 and resid > 1e-5:
        raise RuntimeError(f"GMRES did not converge (info={info}, residual {resid:.1e})")
    return sol


def footprint(args):
    import scipy.sparse as sp
    from rex_env_extension import N_SECTORS, read_cell, find_file, load_q_rex
    import h5py

    run = Path(args.run)
    meta = json.load(open(run / "meta.json"))
    countries = meta["countries"]
    n = len(countries) * N_SECTORS
    c_idx = countries.index(args.country)

    T = t_memmap(args.t_mat, args.cache)
    with h5py.File(find_file(args.rex_dir, "Y_REX")) as f:
        Y = f[next(k for k in f if not k.startswith("#"))][()]   # (countries, n) = Y'
    x_path = Path(args.cache) / "x_REX3.npy"
    if x_path.exists():
        x = np.load(x_path)
    else:
        x = row_sums(T) + Y.sum(axis=0)
        np.save(x_path, x)
    y = Y[c_idx].astype(float)
    xy_path = Path(args.cache) / f"x_{args.country}.npy"
    if xy_path.exists():
        xy = np.load(xy_path)
    else:
        xy = solve_output(T, x, y)
        np.save(xy_path, xy)

    F = sp.load_npz(run / args.f_file).tocsc()
    keys = flow_keys(run)
    cfs = pd.read_csv(args.cfs)
    cats = list(dict.fromkeys(cfs.category))
    units = cfs.drop_duplicates("category").set_index("category").unit
    row = {k: i for i, k in enumerate(keys)}
    cf = cfs[[(d, c) in row for d, c in zip(cfs.database, cfs.code)]]
    C = sp.csr_matrix((cf.cf, ([cats.index(c) for c in cf.category],
                               [row[(d, c)] for d, c in zip(cf.database, cf.code)])),
                      shape=(len(cats), len(keys)))

    scale = np.where(x > 0, xy / np.where(x > 0, x, 1), 0.0)   # output for the country per unit output
    G = (C @ F) @ sp.diags(scale)                               # impacts per producing country-sector
    G = np.asarray(G.todense()) if sp.issparse(G) else G
    by_country = G.reshape(len(cats), len(countries), N_SECTORS).sum(axis=2)

    Q = load_q_rex(args.rex_dir, len(meta["q_labels"]), n)
    ghg_rex = (Q[1] * scale).sum()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    total = pd.DataFrame({"unit": units.reindex(cats).values, "total": by_country.sum(axis=1)}, index=cats)
    per_country = pd.DataFrame(by_country.T, index=countries, columns=cats)
    total.to_csv(out / f"footprint_{args.country}_{args.label}_total.csv", index_label="category")
    per_country.to_csv(out / f"footprint_{args.country}_{args.label}_by_producer.csv", index_label="producer")
    by_sector = G.reshape(len(cats), len(countries), N_SECTORS).sum(axis=1)
    pd.DataFrame(by_sector.T, index=meta["sectors"], columns=cats).to_csv(
        out / f"footprint_{args.country}_{args.label}_by_sector.csv", index_label="sector")
    print(f"REX GHG footprint (Q_REX, check): {ghg_rex / 1e9:.1f} Mt CO2-eq")
    print(total.to_string(float_format="{:.4g}".format))
    return total, per_country


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("cfs")
    a.add_argument("--project", required=True)
    a.add_argument("--prefix", required=True, help='"EF v3.1" (ecoinvent) or "sentier" (BAFU project)')
    a.add_argument("--out", required=True)
    b = sub.add_parser("footprint")
    b.add_argument("--rex-dir", required=True)
    b.add_argument("--t-mat", required=True, help="T_REX3.mat (extracted)")
    b.add_argument("--cache", required=True, help="folder for the float32 T memmap and output vectors")
    b.add_argument("--run", required=True, help="extension output folder")
    b.add_argument("--f-file", default="F_REX_calibrated.npz")
    b.add_argument("--cfs", required=True)
    b.add_argument("--country", default="Switzerland")
    b.add_argument("--label", required=True)
    b.add_argument("--out", required=True)
    args = ap.parse_args()
    if args.cmd == "cfs":
        export_cfs(args.project, args.prefix, args.out)
    else:
        footprint(args)


if __name__ == "__main__":
    main()
