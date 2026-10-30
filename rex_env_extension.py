"""
REX3 environmental extension from ecoinvent
============================================
Extends the REX3 satellite (Q_REX) with all ecoinvent elementary flows by
scaling sector-level ecoinvent flow ratios with an existing REX extension
("anchor"):

    F[f, c, s] = R[f, s, region(c)] * Q_REX[anchor(f), c, s]
    R[f, s, r] = sum_j flow_f,j / sum_j anchor_j     (j: ecoinvent activities
                                                     of sector s in region r,
                                                     annual = per unit * PV)

Anchors by compartment:
    air, other           -> Climate change impacts [kg CO2-eq]  (IPCC 2021 GWP100)
    water                -> Blue water consumption [Mio. m3]   (water in - water to water)
    soil, resource land  -> Total land use area EXIOBASE3 [km2] (occupation m2a)
    resource in water/air-> Blue water consumption [Mio. m3]
    other resources (kg) -> Total material footprint [kt]

Ratio fallback: regional pool -> global pool -> global pool against GHG.
Regional pools only contain activities whose ecoinvent location maps to a
single Exiobase region; GLO, RoW, RER etc. only enter the global pool.

Mapping files (Peng & Pfister 2023, JIE): sector_matching_dicts.pickle,
S2_ecoinvent_region_TO_exiobase_region.xlsx.

Run with the Activity Browser environment (bw2data/bw2calc 1.x + h5py):
    python rex_env_extension.py --rex-dir <REX3_2022> --map-dir <mapping folder> --out <dir>
"""

import argparse
import difflib
import json
import pickle
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

N_SECTORS = 163

# REX extension rows (Labels_Q_REX3) used as anchors, with the factor that
# converts the ecoinvent anchor unit into the REX unit
ANCHORS = {
    "ghg": (1, 1.0),          # kg CO2-eq
    "water": (3, 1e-6),       # m3 -> Mio. m3
    "land": (5, 1e-6),        # m2a -> km2
    "material": (0, 1e-6),    # kg -> kt
}
GHG_METHOD = ("IPCC 2021", "climate change", "global warming potential (GWP100)")

# Resource flows left out of the material anchor (kg, but not extracted material)
MATERIAL_EXCLUDE = ("Carbon dioxide", "Oxygen", "Nitrogen, in air")

# Several ISIC codes map to many Exiobase sectors (e.g. 3510 -> all electricity
# technologies). Activities are assigned by the first matching keyword on
# "name | reference product"; unmatched ones go to the default (or stay split).
KEYWORDS = [
    # electricity (3510)
    ("i40.12", r"voltage transformation|transmission"),
    ("i40.13", r"distribution"),
    ("i40.11.c", r"nuclear"),
    ("i40.11.d", r"hydro|reservoir|run-of-river|pumped storage"),
    ("i40.11.e", r"wind"),
    ("i40.11.h", r"photovoltaic|\bpv\b"),
    ("i40.11.i", r"solar thermal|solar tower|parabolic|concentrat"),
    ("i40.11.j", r"tidal|wave|ocean"),
    ("i40.11.k", r"geothermal"),
    ("i40.11.g", r"biomass|wood|biogas|bagasse|straw|landfill gas|sewage|waste|biomethane"),
    ("i40.11.a", r"coal|lignite|peat"),
    ("i40.11.f", r"\boil\b|diesel|petroleum|refinery gas"),
    ("i40.11.b", r"natural gas|\bgas\b|methane"),
    # basic non-ferrous metals (2420) and ores (0729)
    ("i27.41", r"gold|silver|platinum|palladium|rhodium|iridium|ruthenium|osmium|precious"),
    ("i27.42", r"alumin"),
    ("i27.43", r"\blead\b|zinc|\btin\b"),
    ("i27.44", r"copper"),
    ("i13.20.11", r"copper"),
    ("i13.20.12", r"nickel"),
    ("i13.20.13", r"bauxite|alumin"),
    ("i13.20.14", r"gold|silver|platinum|palladium|precious"),
    ("i13.20.15", r"\blead\b|zinc|\btin\b"),
    # quarrying, gas, food, chemicals, crops
    ("i14.2", r"sand|clay|gravel|kaolin|bentonite"),
    ("i11.c", r"liquef|\blng\b|regasif"),
    ("i15.a", r"cattle|beef|\bcow\b|veal"),
    ("i15.b", r"\bpig|pork|swine"),
    ("i15.c", r"chicken|poultry|broiler|turkey|duck"),
    ("i15.g", r"\brice\b"),
    ("i24.1", r"poly|pvc|resin|plastic|nylon|\babs\b|\bpet\b"),
    ("i24.2", r"urea|ammonium|ammonia|nitrogen fertili"),
    ("i24.3", r"phosph|potass|potash|\bp\b fertili|\bk\b fertili"),
    ("i01.b", r"wheat"),
    ("i01.e", r"rape|soybean|sunflower|oil ?seed|groundnut|peanut|linseed|mustard"),
    ("i23.3", r"nuclear|radioactive"),
]
DEFAULTS = ["i40.11.l", "i27.45", "i13.20.16", "i14.1", "i11.b", "i15.d", "i15.i",
            "i24.4", "i01.c", "i37", "i22"]
# waste treatment (3811/3812): treatment type x material
WASTE_TREATMENT = [("i90.1", r"incinerat|open burning"), ("i90.6", r"landfill|open dump"), ("i90.4", r"compost"),
                   ("i90.3", r"anaerobic|digestion|biogas"), ("i90.5", r"wastewater|waste water|sewage")]
WASTE_MATERIAL = [("a", r"food|biowaste|bio-waste|manure"), ("b", r"paper|cardboard"),
                  ("c", r"plastic|poly"), ("e", r"textile"), ("f", r"wood"),
                  ("g", r"\boil\b|hazardous|solvent|spent"), ("d", r"metal|inert|glass|concrete|mineral")]

# REX country -> Exiobase region where the name does not match the S2 table
COUNTRY_ALIASES = {
    "Czech Republic": "CZ", "UK": "GB", "USA": "US", "South Korea": "KR", "Taiwan": "TW",
    "Afghanistan": "WA", "Bhutan": "WA", "Brunei": "WA", "North Korea": "WA", "Fiji": "WA",
    "French Polynesia": "WA", "Laos": "WA", "Macao SAR": "WA", "Maldives": "WA",
    "New Caledonia": "WA", "Samoa": "WA", "Vanuatu": "WA", "Former USSR": "WE",
    "Antigua": "WL", "Aruba": "WL", "Bahamas": "WL", "Barbados": "WL", "Belize": "WL",
    "Bermuda": "WL", "Bolivia": "WL", "British Virgin Islands": "WL", "Cayman Islands": "WL",
    "Guyana": "WL", "Netherlands Antilles": "WL", "Suriname": "WL", "Greenland": "WE",
    "Andorra": "WE", "Liechtenstein": "WE", "Monaco": "WE", "Moldova": "WE",
    "San Marino": "WE", "TFYR Macedonia": "WE",
    "Burkina Faso": "WF", "Burundi": "WF", "Cape Verde": "WF", "Central African Republic": "WF",
    "Chad": "WF", "Cote dIvoire": "WF", "DR Congo": "WF", "Djibouti": "WF", "Gambia": "WF",
    "Guinea": "WF", "Lesotho": "WF", "Liberia": "WF", "Malawi": "WF", "Mali": "WF",
    "Mauritania": "WF", "Rwanda": "WF", "Sao Tome and Principe": "WF", "Seychelles": "WF",
    "Sierra Leone": "WF", "Somalia": "WF", "Swaziland": "WF", "Uganda": "WF", "Tanzania": "WF",
    "Iran": "WM", "Gaza Strip": "WM", "Syria": "WM", "UAE": "WM",
}


# ---------------------------------------------------------------------------
# REX inputs
# ---------------------------------------------------------------------------
def read_cell(path, key):
    """Read a MATLAB v7.3 cell array of strings."""
    with h5py.File(path) as f:
        return ["".join(map(chr, f[r][()].ravel())) for r in f[key][()].ravel()]


def find_file(folder, stem):
    hits = sorted(Path(folder).rglob(f"{stem}*.mat"))
    if not hits:
        raise FileNotFoundError(f"{stem}*.mat not found in {folder}")
    return hits[0]


def load_rex_labels(label_dir):
    countries = read_cell(find_file(label_dir, "Labels_Countries_REX3"), "Labels_Countries_REX3")
    sectors = read_cell(find_file(label_dir, "Labels_Sectors_REX3"), "Labels_Sectors_REX3")[:N_SECTORS]
    q_labels = read_cell(find_file(label_dir, "Labels_Q_REX3"), "Labels_Q_REX3")
    return countries, sectors, q_labels


def load_q_rex(rex_dir, n_ext, n_cols):
    """Load Q_REX as (extensions, country x sector), whatever the stored orientation."""
    path = find_file(rex_dir, "Q_REX")
    try:
        with h5py.File(path) as f:
            key = next(k for k in f if not k.startswith("#"))
            q = f[key][()]
    except OSError:  # pre-v7.3 mat file
        from scipy.io import loadmat
        m = loadmat(path)
        q = next(v for k, v in m.items() if not k.startswith("__"))
    q = np.asarray(q, dtype=float)
    if q.shape == (n_cols, n_ext):
        q = q.T
    if q.shape != (n_ext, n_cols):
        raise ValueError(f"Q_REX has shape {q.shape}, expected ({n_ext}, {n_cols})")
    return q


# ---------------------------------------------------------------------------
# Mappings
# ---------------------------------------------------------------------------
def _norm(s):
    return re.sub(r"\s*\(\d+(\.\d+)?\)\s*$", "", s.strip().strip('"')).strip().lower()


def rex_sector_codes(sectors, exio_name_dict):
    """Exiobase code for each REX sector (REX sector names have ISIC suffixes)."""
    by_name = {_norm(v): k for k, v in exio_name_dict.items()}
    codes = [by_name.get(_norm(s)) for s in sectors]
    left = [k for k in exio_name_dict if k not in codes]
    for i, c in enumerate(codes):
        if c is None:
            names = {_norm(exio_name_dict[k]): k for k in left}
            best = difflib.get_close_matches(_norm(sectors[i]), list(names), 1, 0.0)[0]
            codes[i] = names[best]
            left.remove(codes[i])
    assert len(set(codes)) == N_SECTORS
    return codes


def refine(cands, text):
    """Narrow multiple Exiobase candidates using activity / product keywords."""
    if len(cands) == 1:
        return cands
    text = text.lower()
    if any(c.startswith("i90.") for c in cands):
        if re.search(r"nuclear|radioactive", text) and "i23.3" in cands:
            return ["i23.3"]
        for prefix, pat in WASTE_TREATMENT:
            if re.search(pat, text):
                sub = [c for c in cands if c.startswith(prefix)]
                for letter, mpat in WASTE_MATERIAL:
                    if re.search(mpat, text) and f"{prefix}.{letter}" in sub:
                        return [f"{prefix}.{letter}"]
                if prefix == "i90.5" and "i90.5.b" in sub:
                    return ["i90.5.b"]
                return sub or cands
        return cands
    for code, pat in KEYWORDS:
        if code in cands and re.search(pat, text):
            return [code]
    default = [c for c in DEFAULTS if c in cands]
    return default[:1] or cands


def isic_sectors(isic, isic_to_exio, code_index, text=""):
    """REX sector indices for an ecoinvent ISIC code, coarsening if unmatched."""
    isic = isic.split(":")[0].strip()
    for n in (len(isic), 3, 2):
        key = isic[:n]
        hits = isic_to_exio.get(key)
        if hits is None:
            hits = sorted({e for k, v in isic_to_exio.items() if k.startswith(key) for e in v})
        if hits:
            hits = refine(list(hits), text)
            return [code_index[h] for h in hits if h in code_index]
    return []


def location_regions(s2):
    """ecoinvent location -> Exiobase region, single-region locations only."""
    out = {}
    for loc, reg in zip(s2["ecoinvent location"], s2["exiobase region"]):
        if isinstance(loc, str) and isinstance(reg, str) and "," not in reg:
            out[loc.strip()] = reg.strip()
    return out


def country_regions(countries, s2):
    by_name = {n.strip(): r for n, r in zip(s2["ecoinvent location - full name"], s2["exiobase region"])
               if isinstance(r, str) and "," not in r}
    regions = []
    for c in countries:
        r = COUNTRY_ALIASES.get(c, by_name.get(c))
        if r is None:
            raise KeyError(f"No Exiobase region for REX country '{c}'")
        regions.append(r)
    return regions


# ---------------------------------------------------------------------------
# ecoinvent pools
# ---------------------------------------------------------------------------
def flow_anchor(categories, name, unit):
    comp = categories[0]
    sub = categories[1] if len(categories) > 1 else ""
    if comp == "water":
        return "water"
    if comp == "soil":
        return "land"
    if comp == "natural resource":
        if sub == "land":
            return "land"
        if name.startswith("Water") and unit == "cubic meter":
            return "water"
        return "material"
    return "ghg"


def ecoinvent_pools(project, database, isic_to_exio, code_index, loc_region, regions):
    """Annual flow totals and anchor totals per (sector, pool).

    Pools: one per Exiobase region (index as in `regions`) + global (last).
    Returns flows (list of biosphere keys + metadata), E (n_flows, n_sectors, n_pools),
    A (n_anchors, n_sectors, n_pools), n_acts (n_sectors, n_pools).
    """
    import bw2data as bd
    import bw2calc as bc
    from bw2data.backends.peewee import ActivityDataset, ExchangeDataset

    bd.projects.set_current(project)
    db = bd.Database(database)

    pv = {}
    for ex in ExchangeDataset.select().where(
        (ExchangeDataset.output_database == database) & (ExchangeDataset.type == "production")
    ):
        pv[(database, ex.output_code)] = ex.data.get("production volume") or 0.0

    lca = bc.LCA({next(iter(db)): 1}, GHG_METHOD)
    lca.load_lci_data()
    lca.load_lcia_data()
    B = lca.biosphere_matrix.tocsc()
    diag = np.abs(lca.technosphere_matrix.diagonal())  # waste treatments: -1
    cf = np.asarray(lca.characterization_matrix.diagonal()).ravel()

    rev_bio = {v: k for k, v in lca.biosphere_dict.items()}
    flows = []
    for i in range(B.shape[0]):
        fl = bd.get_activity(rev_bio[i])
        cats = tuple(fl.get("categories") or ("",))
        flows.append({"key": rev_bio[i], "code": rev_bio[i][1], "name": fl["name"],
                      "compartment": "/".join(cats), "unit": fl.get("unit", ""),
                      "anchor": flow_anchor(cats, fl["name"], fl.get("unit", ""))})
    names = [f["name"] for f in flows]
    cats0 = [f["compartment"].split("/") for f in flows]
    occ = np.array([c[0] == "natural resource" and len(c) > 1 and c[1] == "land"
                    and n.startswith("Occupation") for c, n in zip(cats0, names)], float)
    # water consumption = withdrawals - releases to water (ecoinvent water balance)
    w_in = np.array([f["anchor"] == "water" and c[0] == "natural resource" and n != "Water, in air"
                     for c, n, f in zip(cats0, names, flows)], float)
    w_out = np.array([c[0] == "water" and n == "Water" for c, n in zip(cats0, names)], float)
    mat = np.array([f["anchor"] == "material" and f["unit"] == "kilogram"
                    and not n.startswith(MATERIAL_EXCLUDE) for f, n in zip(flows, names)], float)

    n_pool = len(regions) + 1
    reg_index = {r: i for i, r in enumerate(regions)}
    E = np.zeros((B.shape[0], N_SECTORS, n_pool))
    A = np.zeros((len(ANCHORS), N_SECTORS, n_pool))
    n_acts = np.zeros((N_SECTORS, n_pool), int)
    n_pos = np.zeros((len(ANCHORS), N_SECTORS, n_pool), int)  # activities with anchor > 0
    anchor_names = list(ANCHORS)

    for rec in ActivityDataset.select().where(ActivityDataset.database == database):
        d = rec.data
        if d.get("activity type") != "ordinary transforming activity":
            continue
        isic = dict(d.get("classifications") or []).get("ISIC rev.4 ecoinvent")
        text = f"{d.get('name', '')} | {d.get('reference product', '')}"
        secs = isic_sectors(isic, isic_to_exio, code_index, text) if isic else []
        scale = pv.get((database, rec.code), 0.0)
        j = lca.activity_dict[(database, rec.code)]
        if not secs or scale <= 0 or diag[j] == 0:
            continue
        col = B[:, j].toarray().ravel() * (scale / diag[j] / len(secs))
        anchors = {
            "ghg": cf @ col,
            "water": max(w_in @ col - w_out @ col, 0.0),
            "land": occ @ col,
            "material": mat @ col,
        }
        a_vec = np.array([anchors[k] * ANCHORS[k][1] for k in anchor_names])
        pools = [n_pool - 1]
        r = loc_region.get(rec.location)
        if r in reg_index:
            pools.append(reg_index[r])
        for s in secs:
            for p in pools:
                E[:, s, p] += col
                A[:, s, p] += a_vec
                n_pos[:, s, p] += a_vec > 0
                n_acts[s, p] += 1
    return flows, E, A, n_acts, n_pos


def ratios(flows, E, A, n_acts, n_pos, min_acts=1, min_share=0.5):
    """R[f, s, pool] and the fallback level used per (anchor, sector, pool).

    A pool/anchor is usable if its total is > 0, it has >= min_acts activities and
    >= min_share of them have a positive anchor (avoids exploding ratios on
    near-zero anchors). Levels: 0 regional, 1 global, 2 regional GHG, 3 global GHG, 4 none.
    """
    anchor_names = list(ANCHORS)
    a_idx = np.array([anchor_names.index(f["anchor"]) for f in flows])
    glob = E.shape[2] - 1
    ok = (A > 0) & (n_acts >= min_acts) & (n_pos >= min_share * np.maximum(n_acts, 1))
    R = np.zeros_like(E)
    level = np.full((len(ANCHORS), N_SECTORS, E.shape[2]), 4, np.int8)
    for a in range(len(ANCHORS)):
        fa = a_idx == a
        for s in range(N_SECTORS):
            for p in range(E.shape[2]):
                chain = [(a, p, 0), (a, glob, 1), (0, p, 2), (0, glob, 3)]
                for anc, pool, lv in chain:
                    if ok[anc, s, pool]:
                        R[fa, s, p] = E[fa, s, pool] / A[anc, s, pool]
                        level[a, s, p] = lv if not (pool == glob and lv == 0) else 1
                        break
    return R, level, a_idx


def build_extension(R, level, a_idx, Q, country_pool):
    """F (n_flows, n_countries * n_sectors), sparse.

    Ratios that fell back to GHG (level 2/3) are scaled with REX GHG instead of
    the compartment anchor.
    """
    q_rows = np.array([v[0] for v in ANCHORS.values()])
    ghg_row = ANCHORS["ghg"][0]
    blocks = []
    for c, p in enumerate(country_pool):
        cols = slice(c * N_SECTORS, (c + 1) * N_SECTORS)
        q_anchor = np.where(level[:, :, p] >= 2, Q[ghg_row, cols], Q[q_rows, cols])  # (anchors, sectors)
        blocks.append(sp.csr_matrix(R[:, :, p] * q_anchor[a_idx]))
    return sp.hstack(blocks, format="csr")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rex-dir", help="folder with Q_REX*.mat (e.g. REX3_2022); omit to only build ratios")
    ap.add_argument("--map-dir", required=True, help="folder with mapping files and REX3_Labels")
    ap.add_argument("--label-dir", help="REX3_Labels folder (default: inside --map-dir)")
    ap.add_argument("--project", default="default")
    ap.add_argument("--database", default="cutoff391")
    ap.add_argument("--min-acts", type=int, default=1, help="min. activities for a regional ratio")
    ap.add_argument("--min-share", type=float, default=0.5,
                    help="min. share of activities with a positive anchor to use that anchor")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    map_dir, out = Path(args.map_dir), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    label_dir = Path(args.label_dir) if args.label_dir else map_dir
    if not list(label_dir.rglob("Labels_Q_REX3*.mat")) and (map_dir / "REX3_Labels.zip").exists():
        import zipfile
        zipfile.ZipFile(map_dir / "REX3_Labels.zip").extractall(out)
        label_dir = out

    countries, sectors, q_labels = load_rex_labels(label_dir)
    dicts = pickle.load(open(map_dir / "sector_matching_dicts.pickle", "rb"))
    s2 = pd.read_excel(map_dir / "S2_ecoinvent_region_TO_exiobase_region.xlsx", sheet_name="region_matching")

    codes = rex_sector_codes(sectors, dicts["exio_name_dict"])
    code_index = {c: i for i, c in enumerate(codes)}
    c_regions = country_regions(countries, s2)
    regions = sorted(set(c_regions))
    country_pool = [regions.index(r) for r in c_regions]

    cache = out / f"ecoinvent_pools_{args.database}.pkl"
    if cache.exists():
        flows, E, A, n_acts, n_pos = pickle.load(open(cache, "rb"))
    else:
        pools = ecoinvent_pools(args.project, args.database, dicts["isic_to_exio"],
                                code_index, location_regions(s2), regions)
        pickle.dump(pools, open(cache, "wb"))
        flows, E, A, n_acts, n_pos = pools

    keep = np.abs(E[:, :, -1]).sum(axis=1) > 0           # flows present in some sector
    flows = [f for f, k in zip(flows, keep) if k]
    E = E[keep]
    R, level, a_idx = ratios(flows, E, A, n_acts, n_pos, args.min_acts, args.min_share)

    pd.DataFrame([{k: v for k, v in f.items() if k != "key"} for f in flows]).assign(
        rex_anchor=[q_labels[ANCHORS[f["anchor"]][0]] for f in flows]
    ).to_csv(out / "F_ecoinvent_REX_flows.csv", index_label="row")
    np.savez_compressed(out / "ratios.npz", R=R, level=level, A=A, n_acts=n_acts)

    pools = regions + ["GLO"]
    cov = []
    for s in range(N_SECTORS):
        row = {"sector": sectors[s], "exio_code": codes[s], "acts_global": n_acts[s, -1],
               "regions_with_acts": int((n_acts[s, :-1] > 0).sum())}
        for a, k in enumerate(ANCHORS):
            lv = level[a, s, country_pool]
            row[f"{k}_regional_countries"] = int((lv == 0).sum())
            row[f"{k}_ghg_fallback_countries"] = int(((lv == 2) | (lv == 3)).sum())
            row[f"{k}_no_ratio"] = bool((lv == 4).all())
        cov.append(row)
    pd.DataFrame(cov).to_csv(out / "sector_coverage.csv", index=False)
    json.dump({"countries": countries, "country_region": c_regions, "pools": pools,
               "sectors": sectors, "exio_codes": codes, "anchors": {k: q_labels[v[0]] for k, v in ANCHORS.items()},
               "ghg_method": GHG_METHOD, "database": args.database},
              open(out / "meta.json", "w"), indent=1)

    print(f"{len(flows)} flows, sectors without any ecoinvent activity: {int((n_acts[:, -1] == 0).sum())}")
    if args.rex_dir is None:
        print("no --rex-dir given: ratios written, F not built")
        return None

    Q = load_q_rex(args.rex_dir, len(q_labels), len(countries) * N_SECTORS)
    F = build_extension(R, level, a_idx, Q, country_pool)
    sp.save_npz(out / "F_ecoinvent_REX.npz", F)
    print(f"F: {F.shape[0]} flows x {F.shape[1]} country-sectors, {F.nnz} nonzeros")
    return F


if __name__ == "__main__":
    main()
