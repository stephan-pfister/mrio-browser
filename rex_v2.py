"""
REX LCI extensions v2: agriculture merged, calibration beyond air, fixes, sensitivity
====================================================================================
Starts from the raw extensions of rex_env_extension.py (F_REX.npz) and applies, in order:

1. agriculture: crop sectors (0-7) from HESTIA on-farm emissions (rex_hestia.py), livestock
   sectors (8-14) from Agribalyse farm subsystems (rex_agribalyse.py). Their pollutants
   (NH3, N2O, CH4, CO2, NOx, NMVOC, CO, SO2, PM2.5 to air; nitrate and phosphorus to water)
   replace the same pollutants of the LCI extension in those sectors; all other flows
   (land, water, resources, metals, ...) stay from ecoinvent / BAFU.
2. air: pigs and poultry are calibrated fully to EXIOBASE (their GHG is too small to anchor
   NH3), then all sectors as in rex_calibrate.py (cells off by more than --threshold).
3. beyond air, per country x sector against REX's own extensions, same threshold:
   - land occupation vs REX land (transformations scaled along, then "from" = "to" per cell,
     since the per-flow outlier cap breaks their area balance)
   - fresh water consumption (withdrawals - releases, without salt water, turbine water and
     releases to the ocean) vs REX blue water; turbine water is removed on both sides and a
     negative balance is closed (releases dropped where REX has blue water, else cut to zero net)
   - extracted mass (energy carriers via heating value) vs REX material; extraction in sectors
     without REX material is dropped (REX books it on the extracting sectors)
4. EF factors: unspecified chromium gets the Cr(III) factors (as BAFU's own mapping does).

Variants for sensitivity (all for both sources):
    base         all of the above, threshold 3, global ratios, 90th-percentile cap
    thr1.5       threshold 1.5 (air and beyond air)
    regional     regional ratios with global fallback (needs <run>_regional)
    nocap        no outlier cap (needs <run>_nocap)
    ghg_livestock  pigs/poultry with the GHG anchor and the 3x rule only
    no_agri      without step 1
    no_beyond    without step 3

    python rex_v2.py --rex <REX coinvent folder> --cache <rex_cache> [--variants base thr1.5 ...]
Writes F_REX_v2_<variant>.npz into each run folder, footprints into <rex>/footprint_CH/v2 and
the sensitivity table sensitivity_EF31.csv there.
"""

import argparse
import json
import pickle
import re
import subprocess
import sys
from argparse import Namespace
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

import rex_agribalyse as ra
import rex_hestia as rh
from rex_calibrate import calibrate, load_exiobase, scale_rows
from rex_env_extension import N_SECTORS, anchor_vectors, load_q_rex, norm_categories
from rex_footprint import footprint

SOURCES = {"ecoinvent": ("out_ecoinvent_final", "cfs_ecoinvent.csv"),
           "bafu": ("out_bafu_final", "cfs_bafu.csv")}
VARIANTS = {
    "base": {},
    "thr1.5": {"threshold": 1.5},
    "regional": {"run_suffix": "_regional"},
    "nocap": {"run_suffix": "_nocap"},
    "ghg_livestock": {"full_livestock": False},
    "no_agri": {"agri": False},
    "no_beyond": {"beyond": False},
}
FULL_LIVESTOCK = [9, 10]                     # pigs, poultry
# canonical pollutant, compartment group, names in ecoinvent / EF 3.1 / HESTIA / Agribalyse (lower case)
SPECIES = [
    ("NH3", "air", r"^(ammonia|nh3)$"),
    ("N2O", "air", r"^(dinitrogen monoxide|nitrous oxide|n2o)$"),
    ("CH4", "air", r"^(methane, (non-fossil|biogenic)|methane \(biogenic\)|ch4)$"),
    ("CH4 fossil", "air", r"^(methane, fossil|methane \(fossil\))$"),
    ("CO2", "air", r"^(carbon dioxide, fossil|carbon dioxide \(fossil\)|co2)$"),
    ("NOx", "air", r"^(nitrogen oxides|nox)$"),
    ("NMVOC", "air", r"^(nmvoc|non-methane volatile organic compounds)"),
    ("CO", "air", r"^(carbon monoxide(, fossil| \(fossil\))?|co)$"),
    ("SO2", "air", r"^(sulfur dioxide|sulphur dioxide|so2)$"),
    ("PM2.5", "air", r"^(particulates?, < 2\.5 um|particles \(pm2\.5\)|pm25)$"),
    ("NO3", "water", r"^(nitrate|no3)$"),
    ("P", "water", r"^(phosphate|phosphorus|p)$"),
]
PREFERRED = {"air": r"non-urban|low population|rural|unspecified$", "water": r"surface water|fresh water|river"}
ENERGY_CARRIERS = [  # name pattern, MJ per kg (heating value), kg per m3
    (r"natural gas|gas, natural", 47.0, 0.75), (r"crude oil|oil, crude", 42.6, 0.0),
    (r"hard coal|coal, hard", 27.0, 0.0), (r"lignite|brown coal|coal, brown", 10.0, 0.0),
    (r"peat", 9.0, 0.0), (r"uranium", 560000.0, 0.0)]
CR_UNSPEC, CR_III = r"^chromium$", r"^chromium,? ?(iii|\(iii\))$"


def run_flows(run):
    """Flow dicts (with categories and key) of the rows of F_REX in a run folder."""
    pools = sorted(Path(run).glob("pools_*.json"))[0]
    flows = json.load(open(pools))
    E = np.load(pools.with_suffix(".npz"))["E"]
    keep = np.abs(E[:, :, -1]).sum(axis=1) > 0
    out = [{**f, "categories": tuple(f["categories"]), "key": tuple(f["key"])} for f, k in zip(flows, keep) if k]
    assert [f["name"] for f in out] == pd.read_csv(Path(run) / "F_REX_flows.csv").name.tolist()
    return out


def canon(name, group):
    return next((sp_ for sp_, g, pat in SPECIES if g == group and re.search(pat, name.lower())), None)


def group_of(compartment):
    c = compartment.lower()
    return "air" if "air" in c else "water" if "water" in c and "resource" not in c else None


def target_rows(flows):
    """Row in the run for each canonical pollutant (preferred sub-compartment first)."""
    rows, all_rows = {}, {}
    for i, f in enumerate(flows):
        g = group_of(f["compartment"]) if not f["compartment"].startswith("natural resource") else None
        sp_ = canon(f["name"], g) if g else None
        if sp_ is None:
            continue
        all_rows.setdefault(sp_, []).append(i)
        if sp_ not in rows or (re.search(PREFERRED[g], f["compartment"], re.I)
                               and not re.search(PREFERRED[g], flows[rows[sp_]]["compartment"], re.I)):
            rows[sp_] = i
    return rows, all_rows


def agri_blocks(hestia_dir, agb_zip, agb_cache, countries, q_ghg):
    """[(F, labels, sectors)]: HESTIA crops and Agribalyse livestock, labels = (group, name)."""
    df = rh.tidy(rh.extract(Path(hestia_dir) / "Cycle", Path(hestia_dir) / "onfarm_emissions.csv"))
    Fh, fl_h, _ = rh.build(rh.ratios(df), countries, q_ghg)
    lab_h = [("air" if f.endswith("to air") else "water", f.split(" to ")[0]) for f in fl_h]
    R, _ = ra.livestock_ratios(agb_zip, agb_cache)
    Fa = ra.scale(R, len(countries), q_ghg)
    lab_a = [(group_of(c), n) for n, c in zip(R.index.get_level_values("name"), R.index.get_level_values("category"))]
    return [(Fh, lab_h, list(range(8))), (Fa, lab_a, [int(s) for s in R.columns])]


def merge_agri(F, flows, blocks, n_countries):
    """Replace the agricultural pollutants of F by the HESTIA / Agribalyse ones."""
    rows, all_rows = target_rows(flows)
    F = F.tocoo()
    log = []
    for Fs, labels, sectors in blocks:
        species = [canon(n, g) if g else None for g, n in labels]
        mapped = [(j, rows[s]) for j, s in enumerate(species) if s in rows]
        log += [{"source_flow": f"{g}: {n}", "species": s, "mapped": s in rows}
                for (g, n), s in zip(labels, species)]
        cols = np.array([c * N_SECTORS + s for c in range(n_countries) for s in sectors])
        drop_rows = np.array(sorted({r for s in set(species) if s in all_rows for r in all_rows[s]}))
        drop = np.isin(F.row, drop_rows) & np.isin(F.col, cols)
        F = sp.coo_matrix((F.data[~drop], (F.row[~drop], F.col[~drop])), shape=F.shape)
        M = sp.csr_matrix((np.ones(len(mapped)), ([r for _, r in mapped], [j for j, _ in mapped])),
                          shape=(F.shape[0], Fs.shape[0]))
        add = (M @ Fs.tocsr()).tocoo()
        keep_cols = np.isin(add.col, cols)
        F = sp.coo_matrix((np.r_[F.data, add.data[keep_cols]],
                           (np.r_[F.row, add.row[keep_cols]], np.r_[F.col, add.col[keep_cols]])), shape=F.shape)
    return F.tocsr(), pd.DataFrame(log).drop_duplicates()


def calibrate_air(F, flows, meta, regions, exio, ghg, threshold, full_livestock):
    cflows = pd.DataFrame({"name": [f["name"] for f in flows], "compartment": [f["compartment"] for f in flows]})
    cells = []
    if full_livestock:
        only = {p: np.where(np.isin(np.arange(N_SECTORS), FULL_LIVESTOCK), x, 0.0) for p, x in exio.items()}
        F, c = calibrate(F, cflows, meta["country_region"], regions, only, ghg, threshold=1.0001)
        cells.append(c.assign(step="pigs/poultry full"))
    F, c = calibrate(F, cflows, meta["country_region"], regions, exio, ghg, threshold=threshold)
    cells.append(c.assign(step=f"air {threshold}x"))
    return F, pd.concat(cells)


def to_kg(flow):
    """kg per unit for energy-carrier resources in MJ or m3 (0 for anything else)."""
    if norm_categories(flow["categories"])[0] != "natural resource":
        return 0.0
    name, unit = flow["name"].lower(), flow["unit"].lower()
    for pat, per_mj, per_m3 in ENERGY_CARRIERS:
        if re.search(pat, name):
            return 1 / per_mj if unit == "megajoule" else per_m3 if unit in ("cubic meter", "sm3") else 0.0
    return 0.0


def calibrate_anchors(F, flows, Q, threshold):
    """Land, water and material flows per country x sector against REX land, blue water, material."""
    occ, w_in, w_out, mat = anchor_vectors(flows)
    land = np.array([norm_categories(f["categories"]) == ("natural resource", "land") for f in flows], float)
    t_from = land * np.array([bool(re.match(r"(transformation, )?from ", f["name"], re.I)) for f in flows])
    t_to = land * np.array([bool(re.match(r"(transformation, )?to ", f["name"], re.I)) for f in flows])
    turb = np.array([w and "turbine" in f["name"].lower() for w, f in zip(w_in, flows)], float)
    salt = np.array([bool(re.search(r"salt|ocean|sole", f["name"], re.I)) for f in flows])
    w_in_nt = w_in * (1 - turb) * ~salt                        # fresh water withdrawals
    w_out = w_out * np.array(["ocean" not in f["compartment"] for f in flows])
    F = F.tocsr()
    v = lambda vec: np.asarray(F.T @ vec).ravel()
    lt = np.log10(threshold)
    off = lambda ext, q: (ext > 0) & (q > 0) & (np.abs(np.log10(np.where((ext > 0) & (q > 0), q / np.where(ext > 0, ext, 1), 1))) > lt)
    stats = {}

    # land: occupation (m2a) vs REX land (km2); occupation and transformation scaled together
    q_land = Q[5] * 1e6
    o = v(occ)
    k = np.where(off(o, q_land), q_land / np.where(o > 0, o, 1), 1.0)
    F = scale_rows(F, np.where(land > 0)[0], k)
    stats["land cells scaled"] = int((k != 1).sum())
    # the per-flow outlier cap breaks the area balance of transformations: "from" = "to" per cell
    to, fr = v(t_to), v(t_from)
    F = scale_rows(F, np.where(t_from > 0)[0], np.where((fr > 0) & (to > 0), to / np.where(fr > 0, fr, 1), 1.0))

    # water: turbine water out on both sides, then consumption vs REX blue water (Mm3)
    q_w = Q[3] * 1e6
    tu, rel, win = v(turb), v(w_out), v(w_in_nt)
    a = np.where(rel > 0, np.maximum(rel - tu, 0) / np.where(rel > 0, rel, 1), 1.0)
    net = win - rel * a
    stats["water cells negative"] = int(((net < 0) & (q_w > 0)).sum())
    # negative balance: releases dropped where REX has blue water, else cut back to a zero balance
    ra = rel * a
    b = np.where(net >= 0, 1.0, np.where(q_w > 0, 0.0, win / np.where(ra > 0, ra, 1)))
    net = win - ra * b
    kw = np.where(off(net, q_w), q_w / np.where(net > 0, net, 1), 1.0)
    F = scale_rows(F, np.where(turb > 0)[0], np.zeros(F.shape[1]))
    F = scale_rows(F, np.where(w_in_nt > 0)[0], kw)
    F = scale_rows(F, np.where(w_out > 0)[0], a * b * kw)
    stats["water cells scaled"] = int((kw != 1).sum())

    # materials: extracted mass (kg; energy carriers via heating value) vs REX material (kt).
    # REX books extraction on the extracting sectors: where REX has no material, LCI extraction is dropped.
    mass = mat + np.array([to_kg(f) for f in flows])
    q_m = Q[0] * 1e6
    m = v(mass)
    km = np.where(off(m, q_m), q_m / np.where(m > 0, m, 1), 1.0)
    km = np.where((q_m <= 0) & (m > 0), 0.0, km)
    F = scale_rows(F, np.where(mass > 0)[0], km)
    stats["material cells scaled"] = int(((km != 1) & (km != 0)).sum())
    stats["material cells dropped"] = int((km == 0).sum())
    return F, stats


def fix_cfs(cfs_path, flows, out_path):
    """Unspecified chromium gets the factors of Cr(III) in the same compartment."""
    cfs = pd.read_csv(cfs_path)
    by_key = {f["key"]: f for f in flows}
    cr3 = {f["compartment"]: f["key"] for f in flows if re.search(CR_III, f["name"].lower())}
    cr3_top = {f["compartment"].split("/")[0]: f["key"] for f in flows if re.search(CR_III, f["name"].lower())}
    n = 0
    for key, f in by_key.items():
        if not re.search(CR_UNSPEC, f["name"].lower()):
            continue
        src = cr3.get(f["compartment"]) or cr3_top.get(f["compartment"].split("/")[0])
        if src is None:
            continue
        cfs = cfs[~((cfs.database == key[0]) & (cfs.code == key[1]))]
        new = cfs[(cfs.database == src[0]) & (cfs.code == src[1])].assign(code=key[1])
        cfs = pd.concat([cfs, new])
        n += 1
    cfs.to_csv(out_path, index=False)
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rex", required=True, help="REX coinvent folder (runs, REX3_2022, EXIOBASE_2022)")
    ap.add_argument("--cache", required=True, help="rex_cache: T memmap, output vectors, CF tables")
    ap.add_argument("--hestia", required=True, help="HESTIA release folder (Cycle/, onfarm_emissions.csv)")
    ap.add_argument("--agribalyse", required=True, help="Agribalyse 3.2 openLCA JSON-LD export (zip)")
    ap.add_argument("--agribalyse-cache", required=True)
    ap.add_argument("--variants", nargs="+", default=list(VARIANTS))
    ap.add_argument("--country", default="Switzerland")
    args = ap.parse_args()

    rex, cache = Path(args.rex), Path(args.cache)
    out = rex / "footprint_CH" / "v2"
    out.mkdir(parents=True, exist_ok=True)
    dicts = pickle.load(open(rex / "sector_matching_dicts.pickle", "rb"))
    blocks, logs, labels = None, [], []

    for source, (run_name, cf_name) in SOURCES.items():
        for variant in args.variants:
            opt = {"threshold": 3.0, "agri": True, "beyond": True, "full_livestock": True, "run_suffix": "",
                   **VARIANTS[variant]}
            run = rex / (run_name.replace("_final", opt["run_suffix"]) if opt["run_suffix"] else run_name)
            if not (run / "F_REX.npz").exists():
                print(f"skip {source} {variant}: {run} has no F_REX.npz")
                continue
            meta = json.load(open(run / "meta.json"))
            flows = run_flows(run)
            n_c = len(meta["countries"])
            Q = load_q_rex(rex / "REX3_2022", len(meta["q_labels"]), n_c * N_SECTORS)
            regions, exio = load_exiobase(rex / "EXIOBASE_2022" / "F.txt", dicts["exio_name_dict"], meta["exio_codes"])
            F = sp.load_npz(run / "F_REX.npz").tocsr()
            if opt["agri"]:
                if blocks is None:
                    blocks = agri_blocks(args.hestia, args.agribalyse, args.agribalyse_cache, meta["countries"], Q[1])
                F, maplog = merge_agri(F, flows, blocks, n_c)
                if variant == "base":
                    maplog.assign(source=source).to_csv(out / f"agri_flow_mapping_{source}.csv", index=False)
            F, cells = calibrate_air(F, flows, meta, regions, exio, Q[1], opt["threshold"], opt["full_livestock"])
            stats = {}
            if opt["beyond"]:
                F, stats = calibrate_anchors(F, flows, Q, opt["threshold"])
            sp.save_npz(run / f"F_REX_v2_{variant}.npz", F)
            cfs = cache / cf_name.replace(".csv", "_v2.csv")
            n_cr = fix_cfs(cache / cf_name, flows, cfs)
            logs.append({"source": source, "variant": variant, "air cells changed": len(cells),
                         "Cr flows refactored": n_cr, **stats})
            print(source, variant, logs[-1])
            label = f"{source}_v2_{variant}"
            footprint(Namespace(run=str(run), f_file=f"F_REX_v2_{variant}.npz", cfs=str(cfs), label=label,
                                out=str(out), rex_dir=str(rex / "REX3_2022"), t_mat=str(cache / "T_REX3.mat"),
                                cache=str(cache), country=args.country))
            labels.append(label)

    pd.DataFrame(logs).to_csv(out / "v2_run_log.csv", index=False)
    tot = {lab: pd.read_csv(out / f"footprint_{args.country}_{lab}_total.csv", index_col=0)["total"] for lab in labels}
    tab = pd.DataFrame(tot)
    v1 = rex / "footprint_CH"
    for s in SOURCES:
        f = v1 / f"footprint_{args.country}_{s}_total.csv"
        if f.exists():
            tab[f"{s}_v1"] = pd.read_csv(f, index_col=0)["total"]
    tab.to_csv(out / "sensitivity_EF31.csv", index_label="category")
    base = [l for l in labels if l.endswith("_base")]
    if len(base) == 2:
        subprocess.run([sys.executable, str(Path(__file__).parent / "compare_rankings.py"), "--dir", str(out),
                        "--country", args.country, "--approaches", *base,
                        "--out", str(out / "ranking_comparison_v2.xlsx")], check=True)


if __name__ == "__main__":
    main()
