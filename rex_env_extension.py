"""
REX3 environmental extension from ecoinvent (or BAFU LCI)
==========================================================
Extends the REX3 satellite (Q_REX) with all LCI elementary flows by scaling
sector-level flow ratios with an existing REX extension ("anchor"):

    F[f, c, s] = R[f, s, region(c)] * Q_REX[anchor(f), c, s]
    R[f, s, r] = sum_j w_j flow_f,j / sum_j w_j anchor_j   (j: activities of sector s
                                                            in region r)
    w_j = production volume (ecoinvent, default) or 1 / direct GHG_j (equal weights;
          always used for BAFU, which has no production volumes)

Anchors by compartment:
    air, other           -> Climate change impacts [kg CO2-eq]  (IPCC 2021 GWP100)
    water                -> Blue water consumption [Mio. m3]   (water in - water to water)
    soil, resource land  -> Total land use area EXIOBASE3 [km2] (occupation m2a)
    resource in water/air-> Blue water consumption [Mio. m3]
    other resources (kg) -> Total material footprint [kt]

Outliers: within each pool, activity intensities (flow per unit of anchor) are
winsorised at the 90th percentile (--cap-quantile) when >= 5 activities emit the flow.

Ratio fallback: regional pool -> global pool -> global pool against GHG.
Regional pools only contain activities whose ecoinvent location maps to a
single Exiobase region; GLO, RoW, RER etc. only enter the global pool.

Mapping files (Peng & Pfister 2023, JIE): sector_matching_dicts.pickle,
S2_ecoinvent_region_TO_exiobase_region.xlsx.

BAFU datasets have no ISIC code: BAFU category -> EcoSpold01 category -> most
frequent ISIC in ecoinvent 3 (bridge written by every ecoinvent run).

Run ecoinvent with the Activity Browser environment (bw2data/bw2calc 1.x + h5py):
    python rex_env_extension.py --rex-dir <REX3_2022> --map-dir <mapping folder> --out <dir>
BAFU pools need Brightway 2.5 (e.g. the Werkbank venv), then F in the AB env:
    <werkbank python> rex_env_extension.py --source bafu --project "bafu test" --database bafu-2026
        --meta <ecoinvent out>/meta.json --bridge <ecoinvent out>/ecospold1_isic_bridge.csv
        --map-dir <mapping folder> --out <bafu out>
    python rex_env_extension.py --source bafu --database bafu-2026 --meta ... --map-dir ...
        --rex-dir <REX3_2022> --out <bafu out>        (reuses the cached pools)
"""

import argparse
import difflib
import json
import pickle
import re
from pathlib import Path

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
BAFU_GHG_METHOD = ("sentier", "EF v3.1", "Climate change")

# BAFU (ecoinvent v2 structure) datasets that are mixes or markets
BAFU_SKIP = r"regional storage|market|mix\b|delivered via network|kbob recommendation"
# BAFU main categories without an EcoSpold01 counterpart in ecoinvent 3:
# -> EcoSpold01 main category, or "isic:<code>" (sector then refined by keywords)
BAFU_CATEGORY_ALIAS = {
    "heat": "isic:3530", "heating": "building components", "electricity by fuel": "isic:3510",
    "wood": "wooden materials", "construction": "construction processes",
    "agricultural": "agricultural production", "water": "water supply", "landfill": "isic:3821",
    "incineration": "isic:3821", "nuclear waste": "isic:3822", "construction waste": "isic:3821",
    "electronics waste": "isic:3821", "transport waste": "isic:3821", "waste": "isic:3821",
    "computers & network": "electronics", "paper+ board": "paper & cardboard",
    "cardboard": "paper & cardboard", "minerals": "construction materials",
    "wastewater treatment": "isic:3700", "compressed air": "mechanical engineering",
    "mechanical": "mechanical engineering", "pipeline": "isic:4930", "recycling": "isic:3830",
    "flooring": "building components", "building processes": "construction processes",
    "power plants": "isic:3510",
}
# BAFU dataset names that override the (coarse) category -> ISIC bridge
BAFU_NAME_RULES = [
    (r"\b(rhodium|platinum|palladium|gold|silver|iridium|ruthenium)\b.*(refinery|primary|at plant)", "2420"),
    (r"at mine|\bore\b|concentrate|beneficiation", "0729"),
    (r"formwork", "4390"),
    (r"^(?!disposal)(?!.*li-ion).*\b(steel|pig iron|sinter)\b", "2410"),
    (r"^transport, aircraft", "5110"),
    (r"^transport, .*(transoceanic|freight ship|tanker)", "5012"),
    (r"^transport, .*barge", "5022"),
    (r"^transport, .*(train|rail)", "4912"),
    (r"^transport, .*(lorry|truck|passenger car|bus|coach|van|motor)", "4923"),
]
# BAFU "fuels / <sub>": EcoSpold01 main category by fuel
BAFU_FUELS = {"natural gas": "natural gas", "uranium": "nuclear power", "oil": "oil",
              "biomass": "biomass", "wood": "wood energy", "hard coal": "hard coal",
              "lignite": "lignite", "synthetic": "isic:1920", "hydrogen": "isic:2011"}

# Resource flows left out of the material anchor (kg, but not extracted material)
MATERIAL_EXCLUDE = ("carbon dioxide", "oxygen", "nitrogen, in air")  # lower case

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
    import h5py
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
    import h5py
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
# Flows and pools
# ---------------------------------------------------------------------------
def norm_categories(cats):
    """ecoinvent-style (compartment, subcompartment) for ecoinvent and EF 3.x flows."""
    cats = tuple(c for c in cats if c) or ("",)
    top = cats[0].lower()
    if top == "emissions" and len(cats) > 1:
        comp = cats[1].lower().replace("emissions to ", "")
        return (comp, cats[2] if len(cats) > 2 else "")
    if top == "resources" and len(cats) > 1:
        sub = {"resources from ground": "in ground", "resources from water": "in water",
               "resources from air": "in air", "resources from biosphere": "biotic"}
        return ("natural resource", sub.get(cats[1].lower(), cats[1]))
    if top == "land use":
        return ("natural resource", "land")
    return (top, cats[1] if len(cats) > 1 else "")


def flow_anchor(categories, name, unit):
    comp, sub = norm_categories(categories)
    if comp == "water":
        return "water"
    if comp == "soil":
        return "land"
    if comp == "natural resource":
        if sub == "land":
            return "land"
        if "water" in name.lower() and unit == "cubic meter":
            return "water"
        return "material"
    return "ghg"


def flow_meta(key, node):
    cats = tuple(node.get("categories") or ("",))
    unit = node.get("unit", "")
    return {"key": key, "name": node["name"], "categories": cats,
            "compartment": "/".join(norm_categories(cats)).rstrip("/"), "unit": unit,
            "anchor": flow_anchor(cats, node["name"], unit)}


def anchor_vectors(flows):
    """Row vectors that turn a biosphere column into the non-GHG anchors."""
    comp = [norm_categories(f["categories"]) for f in flows]
    occ = np.array([c == ("natural resource", "land") and f["unit"] == "square meter-year"
                    for c, f in zip(comp, flows)], float)
    # water consumption = withdrawals - releases to water (water balance)
    w_in = np.array([f["anchor"] == "water" and c[0] == "natural resource" and c[1] != "in air"
                     and f["name"] != "Water, in air" for c, f in zip(comp, flows)], float)
    w_out = np.array([c[0] == "water" and f["name"] == "Water" for c, f in zip(comp, flows)], float)
    mat = np.array([f["anchor"] == "material" and f["unit"] == "kilogram"
                    and not f["name"].lower().startswith(MATERIAL_EXCLUDE) for f in flows], float)
    return occ, w_in, w_out, mat


def cap_intensities(X, a, q, min_n):
    """Winsorise activity intensities (flow per unit of its anchor) at the pool quantile q.

    X: (activities, flows) inventories, a: (activities, flows) anchor of each flow's
    group. Only flows emitted by >= min_n activities with a positive anchor are capped;
    activities without a positive anchor are left as they are.
    """
    pos = a > 0
    r = np.where(pos, X / np.where(pos, a, 1.0), 0.0)
    nz = pos & (r > 0)
    if not (nz.sum(axis=0) >= min_n).any():
        return X
    caps = np.nanquantile(np.where(nz, r, np.nan)[:, nz.sum(axis=0) >= min_n], q, axis=0)
    cap = np.full(X.shape[1], np.inf)
    cap[nz.sum(axis=0) >= min_n] = caps
    return np.where(nz & (r > cap), cap * a, X)


def accumulate(records, B, diag, cf, flows, loc_region, regions, weighting, cap_q=0.9, cap_min=5):
    """Pool activities per (sector, region pool + global).

    records: iterable of (column, location, sector indices, production volume).
    weighting 'pv': annual flows (per unit * production volume);
    weighting 'ghg': each activity normalised to 1 kg CO2-eq direct GHG (equal weights);
    weighting 'none': per-unit inventories (1 unit of each reference product).
    cap_q: outlier cap per pool, see cap_intensities (None or >= 1 disables).
    Returns E (n_flows, n_sectors, n_pools), A (n_anchors, ...), n_acts, n_pos.
    """
    occ, w_in, w_out, mat = anchor_vectors(flows)
    n_pool = len(regions) + 1
    reg_index = {r: i for i, r in enumerate(regions)}
    E = np.zeros((B.shape[0], N_SECTORS, n_pool))
    A = np.zeros((len(ANCHORS), N_SECTORS, n_pool))
    n_acts = np.zeros((N_SECTORS, n_pool), int)
    n_pos = np.zeros((len(ANCHORS), N_SECTORS, n_pool), int)  # activities with anchor > 0
    factors = np.array([v[1] for v in ANCHORS.values()])
    a_of_flow = np.array([list(ANCHORS).index(f["anchor"]) for f in flows])
    cols, a_vecs, members = [], [], {}

    for j, location, secs, pv in records:
        if not secs or diag[j] == 0:
            continue
        col = B[:, j].toarray().ravel() / diag[j]
        if not col.any():  # no direct exchanges (pure technosphere datasets)
            continue
        ghg = cf @ col
        if weighting == "pv":
            scale = pv
        elif weighting == "ghg":
            scale = 1 / ghg if ghg > 0 else 0.0
        else:  # none: per-unit inventories summed as they are
            scale = 1.0
        if scale <= 0:
            continue
        col *= scale / len(secs)
        a_vec = factors * np.array([cf @ col, max(w_in @ col - w_out @ col, 0.0), occ @ col, mat @ col])
        pools = [n_pool - 1]
        r = loc_region.get(location)
        if r in reg_index:
            pools.append(reg_index[r])
        k = len(cols)
        cols.append(sp.csr_matrix(col))
        a_vecs.append(a_vec)
        for s in secs:
            for p in pools:
                members.setdefault((s, p), []).append(k)

    a_vecs = np.array(a_vecs)
    for (s, p), ks in members.items():
        X = sp.vstack([cols[k] for k in ks]).toarray()
        a = a_vecs[ks]
        if cap_q is not None and cap_q < 1 and len(ks) >= cap_min:
            X = cap_intensities(X, a[:, a_of_flow], cap_q, cap_min)
        E[:, s, p] = X.sum(axis=0)
        A[:, s, p] = a.sum(axis=0)
        n_pos[:, s, p] = (a > 0).sum(axis=0)
        n_acts[s, p] = len(ks)
    return E, A, n_acts, n_pos


def ecoinvent_pools(project, database, isic_to_exio, code_index, loc_region, regions, weighting="pv",
                    cap=(0.9, 5)):
    """Pools from ecoinvent 3.x (Brightway 2, bw2calc 1.x): ISIC classification,
    production volumes, ordinary transforming activities only."""
    import bw2data as bd
    import bw2calc as bc
    from bw2data.backends.peewee import ActivityDataset, ExchangeDataset

    bd.projects.set_current(project)
    pv = {ex.output_code: ex.data.get("production volume") or 0.0
          for ex in ExchangeDataset.select().where(
              (ExchangeDataset.output_database == database) & (ExchangeDataset.type == "production"))}

    lca = bc.LCA({next(iter(bd.Database(database))): 1}, GHG_METHOD)
    lca.load_lci_data()
    lca.load_lcia_data()
    B = lca.biosphere_matrix.tocsc()
    diag = np.abs(lca.technosphere_matrix.diagonal())  # waste treatments: -1
    cf = np.asarray(lca.characterization_matrix.diagonal()).ravel()
    rev_bio = {v: k for k, v in lca.biosphere_dict.items()}
    flows = [flow_meta(list(rev_bio[i]), bd.get_activity(rev_bio[i])) for i in range(B.shape[0])]

    def records():
        for rec in ActivityDataset.select().where(ActivityDataset.database == database):
            d = rec.data
            if d.get("activity type") != "ordinary transforming activity":
                continue
            isic = dict(d.get("classifications") or []).get("ISIC rev.4 ecoinvent")
            text = f"{d.get('name', '')} | {d.get('reference product', '')}"
            secs = isic_sectors(isic, isic_to_exio, code_index, text) if isic else []
            yield lca.activity_dict[(database, rec.code)], rec.location, secs, pv.get(rec.code, 0.0)

    return (flows, *accumulate(records(), B, diag, cf, flows, loc_region, regions, weighting, *cap))


def ecospold1_bridge(project, database):
    """EcoSpold01 category -> ISIC counts from ecoinvent 3 classifications."""
    import bw2data as bd
    bd.projects.set_current(project)
    rows = []
    for a in bd.Database(database):
        c = dict(a.get("classifications") or [])
        if "EcoSpold01Categories" in c and "ISIC rev.4 ecoinvent" in c:
            rows.append((c["EcoSpold01Categories"], c["ISIC rev.4 ecoinvent"].split(":")[0]))
    return pd.DataFrame(rows, columns=["ecospold1", "isic"]).value_counts().rename("n").reset_index()


def bafu_pools(project, database, bridge, isic_to_exio, code_index, loc_region, regions, weighting="none",
               cap=(0.9, 5)):
    """Pools from the BAFU LCI database (Brightway 2.5 / bw2calc 2.x, EF 3.1 flows).

    No ISIC and no production volumes: sector via BAFU category -> EcoSpold01 ->
    most frequent ecoinvent ISIC (+ keyword refinement); activities unweighted (per
    unit of product) or equal per kg CO2-eq. Skips obsolete datasets, mixes/markets
    and infrastructure (unit 'unit').
    """
    import bw2data as bd
    import bw2calc as bc

    bd.projects.set_current(project)
    db = bd.Database(database)
    top = bridge.sort_values("n", ascending=False).drop_duplicates("ecospold1")
    cat_isic = dict(zip(top.ecospold1, top.isic))
    main = (bridge.assign(main=bridge.ecospold1.str.split("/").str[0])
            .groupby(["main", "isic"]).n.sum().reset_index()
            .sort_values("n", ascending=False).drop_duplicates("main"))
    main_isic = dict(zip(main["main"], main.isic))

    def isic_of(act):
        rule = next((isic for pat, isic in BAFU_NAME_RULES if re.search(pat, act["name"], re.I)), None)
        if rule:
            return rule
        m = re.search(r"BAFU category: (.*?)\. Source", act.get("comment") or "")
        if not m:
            return None
        cat, _, sub = m.group(1).partition(" / ")
        cat, sub = cat.strip(), sub.split(chr(92))[0].strip()
        if cat == "fuels":
            alias = next((v for k, v in BAFU_FUELS.items() if m.group(1).lower().find(k) > 0), None)
        else:
            alias = BAFU_CATEGORY_ALIAS.get(cat)
        if alias and alias.startswith("isic:"):
            return alias[5:]
        if alias:
            return cat_isic.get(f"{alias}/{sub}") or main_isic.get(alias)
        return cat_isic.get(f"{cat}/{sub}") or main_isic.get(cat)

    lca = bc.LCA({next(iter(db)).id: 1}, BAFU_GHG_METHOD)
    lca.load_lci_data()
    lca.load_lcia_data()
    B = lca.biosphere_matrix.tocsc()
    diag = np.abs(lca.technosphere_matrix.diagonal())
    cf = np.asarray(lca.characterization_matrix.diagonal()).ravel()
    rev_bio = {v: k for k, v in lca.dicts.biosphere.items()}
    flows = []
    for i in range(B.shape[0]):
        node = bd.get_node(id=rev_bio[i])
        flows.append(flow_meta([node["database"], node["code"]], node))

    stats = {"used": 0, "skipped": 0, "no_isic": 0}

    def records():
        for act in db:
            text = f"{act['name']} | {act.get('reference product', '')}"
            if (act.get("sector") == "obsolete" or act.get("unit") == "unit"
                    or re.search(BAFU_SKIP, text, re.I)):
                stats["skipped"] += 1
                continue
            isic = isic_of(act)
            if isic is None:
                stats["no_isic"] += 1
                continue
            stats["used"] += 1
            secs = isic_sectors(isic, isic_to_exio, code_index, text)
            yield lca.dicts.activity[act.id], act.get("location"), secs, 1.0

    pools = (flows, *accumulate(records(), B, diag, cf, flows, loc_region, regions, weighting, *cap))
    print("BAFU activities:", stats)
    return pools


def save_pools(path, pools):
    flows, E, A, n_acts, n_pos = pools
    np.savez_compressed(path.with_suffix(".npz"), E=E, A=A, n_acts=n_acts, n_pos=n_pos)
    json.dump([{**f, "categories": list(f["categories"])} for f in flows],
              open(path.with_suffix(".json"), "w"))


def load_pools(path):
    z = np.load(path.with_suffix(".npz"))
    return json.load(open(path.with_suffix(".json"))), z["E"], z["A"], z["n_acts"], z["n_pos"]


# ---------------------------------------------------------------------------
# Ratios and extension
# ---------------------------------------------------------------------------
def ratios(flows, E, A, n_acts, n_pos, min_acts=1, min_share=0.5, global_only=False):
    """R[f, s, pool] and the fallback level used per (anchor, sector, pool).

    A pool/anchor is usable if its total is > 0, it has >= min_acts activities and
    >= min_share of them have a positive anchor (avoids exploding ratios on
    near-zero anchors). Levels: 0 regional, 1 global, 2 regional GHG, 3 global GHG, 4 none.
    """
    anchor_names = list(ANCHORS)
    a_idx = np.array([anchor_names.index(f["anchor"]) for f in flows])
    glob = E.shape[2] - 1
    ok = (A > 0) & (n_acts >= min_acts) & (n_pos >= min_share * np.maximum(n_acts, 1))
    if global_only:
        ok[:, :, :glob] = False
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
    ap.add_argument("--meta", help="meta.json of an earlier run, replaces the REX label files")
    ap.add_argument("--source", choices=["ecoinvent", "bafu"], default="ecoinvent")
    ap.add_argument("--project", default="default")
    ap.add_argument("--database", default="cutoff391")
    ap.add_argument("--weighting", choices=["pv", "ghg", "none"],
                    help="pv: production volume (ecoinvent default), ghg: equal weight per kg CO2-eq, "
                         "none: per-unit inventories (BAFU default; BAFU has no production volumes)")
    ap.add_argument("--bridge", help="ecospold1_isic_bridge.csv (BAFU only; written by ecoinvent runs)")
    ap.add_argument("--min-acts", type=int, default=3, help="min. activities in a pool for its ratio")
    ap.add_argument("--min-share", type=float, default=0.5,
                    help="min. share of activities with a positive anchor to use that anchor")
    ap.add_argument("--cap-quantile", type=float, default=0.9,
                    help="winsorise activity intensities per pool at this quantile (1 disables)")
    ap.add_argument("--cap-min", type=int, default=5,
                    help="min. emitting activities in a pool before the cap applies")
    ap.add_argument("--global-only", action="store_true",
                    help="use global sector ratios for all countries (no regional pools)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    map_dir, out = Path(args.map_dir), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    dicts = pickle.load(open(map_dir / "sector_matching_dicts.pickle", "rb"))
    s2 = pd.read_excel(map_dir / "S2_ecoinvent_region_TO_exiobase_region.xlsx", sheet_name="region_matching")

    if args.meta:
        meta = json.load(open(args.meta))
        countries, sectors, q_labels = meta["countries"], meta["sectors"], meta["q_labels"]
    else:
        label_dir = Path(args.label_dir) if args.label_dir else map_dir
        if not list(label_dir.rglob("Labels_Q_REX3*.mat")) and (map_dir / "REX3_Labels.zip").exists():
            import zipfile
            zipfile.ZipFile(map_dir / "REX3_Labels.zip").extractall(out)
            label_dir = out
        countries, sectors, q_labels = load_rex_labels(label_dir)

    codes = rex_sector_codes(sectors, dicts["exio_name_dict"])
    code_index = {c: i for i, c in enumerate(codes)}
    c_regions = country_regions(countries, s2)
    regions = sorted(set(c_regions))
    country_pool = [regions.index(r) for r in c_regions]

    weighting = args.weighting or ("none" if args.source == "bafu" else "pv")
    if args.source == "bafu" and weighting == "pv":
        ap.error("BAFU has no production volumes: use --weighting none or ghg")
    cap = (args.cap_quantile, args.cap_min)
    cache = out / f"pools_{args.database}_{weighting}_cap{args.cap_quantile:g}"
    if cache.with_suffix(".npz").exists():
        flows, E, A, n_acts, n_pos = load_pools(cache)
    else:
        if args.source == "bafu":
            bridge = pd.read_csv(args.bridge, dtype={"isic": str})
            pools = bafu_pools(args.project, args.database, bridge, dicts["isic_to_exio"],
                               code_index, location_regions(s2), regions, weighting, cap)
        else:
            pools = ecoinvent_pools(args.project, args.database, dicts["isic_to_exio"],
                                    code_index, location_regions(s2), regions, weighting, cap)
            ecospold1_bridge(args.project, args.database).to_csv(out / "ecospold1_isic_bridge.csv", index=False)
        save_pools(cache, pools)
        flows, E, A, n_acts, n_pos = pools

    keep = np.abs(E[:, :, -1]).sum(axis=1) > 0           # flows present in some sector
    flows = [f for f, k in zip(flows, keep) if k]
    E = E[keep]
    R, level, a_idx = ratios(flows, E, A, n_acts, n_pos, args.min_acts, args.min_share, args.global_only)

    pd.DataFrame([{k: v for k, v in f.items() if k not in ("key", "categories")} for f in flows]).assign(
        rex_anchor=[q_labels[ANCHORS[f["anchor"]][0]] for f in flows]
    ).to_csv(out / "F_REX_flows.csv", index_label="row")
    np.savez_compressed(out / "ratios.npz", R=R, level=level, A=A, n_acts=n_acts)

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
    json.dump({"countries": countries, "country_region": c_regions, "pools": regions + ["GLO"],
               "sectors": sectors, "exio_codes": codes, "q_labels": q_labels,
               "anchors": {k: q_labels[v[0]] for k, v in ANCHORS.items()},
               "source": args.source, "database": args.database, "weighting": weighting},
              open(out / "meta.json", "w"), indent=1)

    print(f"{len(flows)} flows, sectors without any activity: {int((n_acts[:, -1] == 0).sum())}")
    if args.rex_dir is None:
        print("no --rex-dir given: ratios written, F not built")
        return None

    Q = load_q_rex(args.rex_dir, len(q_labels), len(countries) * N_SECTORS)
    F = build_extension(R, level, a_idx, Q, country_pool)
    sp.save_npz(out / "F_REX.npz", F)
    print(f"F: {F.shape[0]} flows x {F.shape[1]} country-sectors, {F.nnz} nonzeros")
    return F


if __name__ == "__main__":
    main()
