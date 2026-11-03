"""
REX livestock-sector extension from Agribalyse 3.2 farm datasets
================================================================
Agribalyse livestock datasets ("... at farm (gate)") carry their emissions in on-farm
sub-processes (enteric fermentation, manure management, housing), not as direct
exchanges. The on-farm inventory of each livestock product is therefore the biosphere of
its farm subsystem:

    farm nodes: livestock processes plus enteric / manure / housing / grazing sub-processes
    g = B_farm (I - A_farm)^-1 e_product      feed, energy and background are cut off

Processes "Copied from Ecoinvent" are never part of the farm subsystem; inputs that
"Adapted from WFLDB" farm sub-processes carry are included like Agribalyse's own.

Scaling as in rex_env_extension.py with the GHG anchor and one global pool per sector:

    R[f, s] = sum_j g(f, j) / sum_j ghg(j)      j: livestock products of sector s (per unit)
    F[f, c, s] = R[f, s] x Q_REX[GHG, c, s], then calibrated to EXIOBASE air emissions of the
    livestock sectors (rex_calibrate.calibrate; flow split from Agribalyse, levels from EXIOBASE)

ghg(j) leaves out manure CH4 (exchange descriptions "feces / storage", or WFLDB manure
management processes): REX livestock GHG follows EXIOBASE, whose livestock CH4 is enteric.
Manure CH4 is still an emission flow of the extension.

    python rex_agribalyse.py --zip <openLCA JSON-LD export> --cache <folder>
        --rex-dir <REX3_2022> --exio <EXIOBASE F.txt> --map-dir <mapping folder>
        --meta <extension meta.json> --compare <run> ... --out <out dir>
"""

import argparse
import json
import pickle
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from rex_calibrate import POLLUTANTS, calibrate, load_exiobase
from rex_env_extension import N_SECTORS, load_q_rex
from rex_hestia import sector_totals

LIVESTOCK_SECTORS = range(8, 15)             # cattle, pigs, poultry, meat animals nec, animal products nec, raw milk, wool
SECTOR_RULES = [                             # first match wins
    (14, r"\bwool\b"),
    (12, r"\beggs?\b|honey"),
    (8, r"\bcalf\b|calves|heifer|cull cow|\bbull\b|steer|beef|suckler|young cattle|cattle"),
    (13, r"\bmilk\b"),
    (9, r"\bpigs?\b|\bsows?\b|piglet|pork|swine"),
    (10, r"broiler|chicken|poultry|turkey|\bduck|guinea fowl|goose|\bhens?\b"),
    (11, r"sheep|lamb|goat|\bewe|rabbit|horse"),
]
FARM_NODE = re.compile(r"enteric|manure management|manure, |slurry|housing|grazing|excretion|at farm", re.I)
NOT_FARM = re.compile(r"feed|fodder|forage|silage|\bhay\b|concentrate|ration|basket|Copied from Ecoinvent"
                      r"|processed in|at plant|at slaughter|at supermarket|at consumer|packaging", re.I)
MASS = {"kg": 1.0, "g": 1e-3, "mg": 1e-6, "µg": 1e-9, "ug": 1e-9, "ng": 1e-12, "t": 1e3, "ton": 1e3}
MANURE_CH4 = re.compile(r"feces|faeces|manure|excret|stockage|storage", re.I)
GWP = [(r"^dinitrogen monoxide$|^nitrous oxide$", 273.0), (r"^methane, fossil", 29.8),
       (r"^methane(, biogenic|, non-fossil)?$", 27.0), (r"^carbon dioxide, fossil", 1.0)]


def sector_of(name):
    for s, pat in SECTOR_RULES:
        if re.search(pat, name, re.I):
            return s
    return None


def extract(zip_path, cache):
    """Processes, technosphere inputs and biosphere exchanges per unit reference product."""
    if Path(cache).exists():
        return pickle.load(open(cache, "rb"))
    z = zipfile.ZipFile(zip_path)
    procs, tech, bio = [], [], []
    for n in z.namelist():
        if not (n.startswith("processes/") and n.endswith(".json")):
            continue
        p = json.loads(z.read(n).decode("utf-8-sig"))
        ex = p.get("exchanges", [])
        ref = next((e for e in ex if e.get("isQuantitativeReference")), None)
        if ref is None or not ref.get("amount"):
            continue
        procs.append({"id": p["@id"], "name": p["name"], "location": (p.get("location") or {}).get("name"),
                      "ref_flow": ref["flow"]["@id"]})
        for e in ex:
            if e is ref:
                continue
            f, amount = e["flow"], e.get("amount", 0) / ref["amount"]
            if f.get("flowType") == "ELEMENTARY_FLOW":
                unit = (e.get("unit") or {}).get("name", "")
                bio.append({"process": p["@id"], "name": f["name"], "category": f.get("category", ""),
                            "unit": "kg" if unit in MASS else unit,
                            "amount": amount * MASS.get(unit, 1.0) * (-1 if e.get("isInput") else 1),
                            "description": (e.get("description") or "")[:200]})
            elif e.get("isInput"):
                tech.append({"process": p["@id"], "flow": f["@id"], "amount": amount})
    data = pd.DataFrame(procs), pd.DataFrame(tech), pd.DataFrame(bio)
    pickle.dump(data, open(cache, "wb"))
    return data


def farm_inventories(procs, tech, bio):
    """On-farm inventory (flows x livestock products) per unit product, and the product table."""
    procs = procs.assign(sector=procs.name.map(sector_of))
    is_farm = procs.name.str.contains(FARM_NODE) & ~procs.name.str.contains(NOT_FARM)
    targets = procs[is_farm & procs.sector.notna() & ~procs.name.str.contains(r"enteric|manure|housing", case=False)]
    nodes = procs[is_farm].reset_index(drop=True)
    idx = {pid: i for i, pid in enumerate(nodes.id)}
    provider = dict(zip(nodes.ref_flow, nodes.id))          # product flow -> farm process
    t = tech[tech.process.isin(idx) & tech.flow.isin(provider)]
    A = sp.csc_matrix((t.amount, ([idx[provider[f]] for f in t.flow], [idx[p] for p in t.process])),
                      shape=(len(nodes), len(nodes)))
    b = bio[bio.process.isin(idx)].copy()
    # manure CH4 is kept as an emission but not counted in the anchor (REX/EXIOBASE livestock GHG is enteric)
    manure_proc = set(nodes.id[nodes.name.str.contains(r"manure management", case=False)])
    b["part"] = np.where(b.name.str.match(r"methane", case=False)
                         & (b.description.str.contains(MANURE_CH4) | b.process.isin(manure_proc)), "manure", "")
    keys = b[["name", "category", "unit", "part"]].drop_duplicates().reset_index(drop=True)
    kidx = {tuple(k): i for i, k in enumerate(keys.itertuples(index=False))}
    B = sp.csr_matrix((b.amount, ([kidx[k] for k in zip(b.name, b.category, b.unit, b.part)],
                                  [idx[p] for p in b.process])), shape=(len(keys), len(nodes)))
    D = np.zeros((len(nodes), len(targets)))
    D[[idx[p] for p in targets.id], np.arange(len(targets))] = 1.0
    X = spla.spsolve((sp.identity(len(nodes), format="csc") - A), D)
    G = B @ (X.reshape(len(nodes), -1))
    return pd.DataFrame(G, index=pd.MultiIndex.from_frame(keys), columns=targets.id), targets


def ghg(inv, manure_ch4=False):
    """IPCC 2021 GWP100 per product; manure CH4 only counted if manure_ch4."""
    names = inv.index.get_level_values("name").str.lower()
    air = (inv.index.get_level_values("category").str.contains("air", case=False)
           & (manure_ch4 | (inv.index.get_level_values("part") != "manure")))
    w = np.zeros(len(inv))
    for pat, gwp in GWP:
        w = np.where((w == 0) & air & names.str.match(pat), gwp, w)
    return (inv.mul(w, axis=0)).sum(axis=0)


def livestock_ratios(zip_path, cache, manure_ch4=False):
    """R (flows x livestock sectors, per kg CO2-eq anchor GHG) and the products behind it."""
    procs, tech, bio = extract(zip_path, Path(cache) / "agribalyse_parsed_v2.pkl")
    inv, targets = farm_inventories(procs, tech, bio)
    g = ghg(inv, manure_ch4)
    keep = g > 0
    inv, targets, g = inv.loc[:, keep], targets[keep.values], g[keep]
    sector = targets.set_index("id").sector.astype(int)
    R = inv.T.groupby(sector).sum().T.div(g.groupby(sector).sum(), axis=1)
    return R, targets.assign(ghg_kg_co2eq=g.values)


def scale(R, n_countries, q_ghg):
    """F (flows x countries*sectors): one global ratio per sector times REX GHG."""
    rows, cols, vals = [], [], []
    for s in R.columns:
        r = R[s].to_numpy()
        nz = np.nonzero(r)[0]
        for ci in range(n_countries):
            q = q_ghg[ci * N_SECTORS + s]
            if q > 0:
                rows += list(nz); cols += [ci * N_SECTORS + s] * len(nz); vals += list(r[nz] * q)
    return sp.csr_matrix((vals, (rows, cols)), shape=(len(R), n_countries * N_SECTORS))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--zip", required=True, help="openLCA JSON-LD export of Agribalyse 3.2")
    ap.add_argument("--cache", required=True, help="folder for the parsed export")
    ap.add_argument("--rex-dir", required=True)
    ap.add_argument("--exio", required=True)
    ap.add_argument("--map-dir", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--compare", nargs="*", default=[])
    ap.add_argument("--manure-ch4-in-anchor", action="store_true",
                    help="count manure CH4 in the anchor GHG (default: enteric CH4 only)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    meta = json.load(open(args.meta))
    countries = meta["countries"]
    R, targets = livestock_ratios(args.zip, args.cache, args.manure_ch4_in_anchor)
    targets.to_csv(out / "livestock_products.csv", index=False)
    R.to_csv(out / "agribalyse_ratios.csv")
    Q = load_q_rex(args.rex_dir, len(meta["q_labels"]), len(countries) * N_SECTORS)
    F = scale(R, len(countries), Q[1])
    sp.save_npz(out / "F_REX_agribalyse_livestock.npz", F)
    flows = R.index.to_frame(index=False)
    flows.to_csv(out / "F_REX_agribalyse_livestock_flows.csv", index=False)

    dicts = pickle.load(open(Path(args.map_dir) / "sector_matching_dicts.pickle", "rb"))
    regions, exio = load_exiobase(args.exio, dicts["exio_name_dict"], meta["exio_codes"])
    secs = list(LIVESTOCK_SECTORS)
    air = flows.category.str.contains("air", case=False)

    # calibration as in rex_calibrate.py, restricted to the livestock sectors this extension covers
    exio_ls = {p: np.where(np.isin(np.arange(N_SECTORS), secs), x, 0.0) for p, x in exio.items()}
    cflows = pd.DataFrame({"name": flows.name, "compartment": np.where(air, "air", "other")})
    Fc, cells = calibrate(F, cflows, meta["country_region"], regions, exio_ls, Q[1])
    sp.save_npz(out / "F_REX_agribalyse_livestock_calibrated.npz", Fc)
    cells.assign(sector=[meta["sectors"][s] for s in cells.sector] if len(cells) else []).to_csv(
        out / "calibration_cells.csv", index=False)

    long = []
    for p, (_, pat) in POLLUTANTS.items():
        sel = np.where(air & flows.name.str.match(pat, case=False))[0]
        res = {"EXIOBASE": exio[p][:, secs].sum(axis=0),
               "Agribalyse": sector_totals(F, sel, len(countries), secs),
               "Agribalyse calibrated": sector_totals(Fc, sel, len(countries), secs)}
        for run in args.compare:
            fl = pd.read_csv(Path(run) / "F_REX_flows.csv")
            Fr = sp.load_npz(Path(run) / "F_REX.npz").tocsr()
            sel = np.where(fl.compartment.str.startswith("air") & fl.name.str.match(pat, case=False))[0]
            res[Path(run).name.replace("out_", "").replace("_final", "") + " raw"] = sector_totals(
                Fr, sel, len(countries), secs)
        for src, v in res.items():
            for s, val in zip(secs, v):
                long.append({"pollutant": p, "source": src, "sector": meta["sectors"][s], "kg": val})
    long = pd.DataFrame(long)
    long.to_csv(out / "livestock_sector_comparison.csv", index=False)
    tot = long.groupby(["pollutant", "source"]).kg.sum().unstack("source").loc[list(POLLUTANTS)] / 1e9
    tot.to_csv(out / "livestock_sector_totals_Mt.csv")
    pd.set_option("display.width", 200)
    print(f"{len(targets)} livestock products:", targets.sector.astype(int).map(
        lambda s: meta["sectors"][s]).value_counts().to_dict())
    print("livestock sectors 8-14, global totals [Mt]:")
    print(tot.to_string(float_format="{:.3g}".format))


if __name__ == "__main__":
    main()
