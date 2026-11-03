"""
REX crop-sector extension from HESTIA on-farm emissions
=======================================================
HESTIA aggregated cycles (CC-BY 4.0, s3://hestia-aggregated-data) give national and world
average crop production per hectare with on-farm emissions. Background emissions (input
production, methodTier "background") are left out: REX covers upstream production itself,
and they are likely ecoinvent-derived.

Same scaling as rex_env_extension.py, with the GHG anchor and per-country pools:

    R[f, s, p] = sum_j e(f, j) / sum_j ghg(j)      j: HESTIA cycles of crop sector s in pool p
                                                    (each cycle = 1 ha, equal weight)
    F[f, c, s] = R[f, s, pool(c)] x Q_REX[GHG, c, s]  pool: country if HESTIA has it, else world

ghg(j) is IPCC 2021 GWP100 of N2O, CH4 and fossil CO2 (fuel, urea, lime). Land-use change,
stock changes, organic soils and natural vegetation burning are excluded from both the
anchor and the emissions, as EXIOBASE crop sectors do not carry them.

    python rex_hestia.py --hestia <release folder with Cycle/> --rex-dir <REX3_2022>
        --exio <EXIOBASE F.txt> --map-dir <mapping folder> --meta <extension meta.json>
        --compare <out_ecoinvent_final> <out_bafu_final> --out <out dir>
"""

import argparse
import json
import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from rex_calibrate import POLLUTANTS, load_exiobase
from rex_env_extension import N_SECTORS, load_q_rex

CROP_SECTORS = {
    0: ["riceGrainInHuskFlooded"],
    1: ["wheatGrain", "durumWheatGrain"],
    2: ["maizeGrain", "barleyGrainWhole", "oatGrain", "oatGrainInHusk", "ryeGrain", "sorghumGrain",
        "fingerMilletGrainInHusk", "pearlMilletGrainInHusk", "quinoaGrain"],
    4: ["rapeseedSeedWhole", "soybeanSeedWhole", "sunflowerSeedWhole", "mustardSeedWhole",
        "sesameSeedWhole", "safflowerSeedWhole", "castorSeed", "peanutInShell", "hempSeedWhole",
        "oilPalmFruit", "oliveFruit"],
    5: ["sugarBeetRoot", "sugarcaneStalk"],
    6: ["cottonLint", "juteStalk"],
    7: ["coffeeSeedDehulled", "cocoaSeedWhole"],
}
# processed products belong to food-processing sectors, not crop cultivation
PROCESSED = re.compile(r"Meal$|Oil(Crude|Refined|Virgin)$|Flour|Bran$|Husk$|Refined$|DistillersGrain|^wine$")
COUNTRY_ALIASES = {"united-states": "USA", "united-kingdom": "UK", "turkey-country": "Turkey",
                   "cote-divoire": "Cote dIvoire", "viet-nam": "Viet Nam"}
EXCLUDE = re.compile(r"LandUseChange|ManagementChange|StockChange|NaturalVegetationBurning|OrganicSoil")
GWP = {"n2O": 273.0, "ch4": 27.0, "ch4_fossil": 29.8, "co2": 1.0}
SPECIES = {"so2": "SOx", "nox": "NOx", "nh3": "NH3", "pm25": "PM2.5", "nmvoc": "NMVOC", "co": "CO",
           "ch4": "CH4", "n2O": "N2O"}


def crop_sector(product):
    if PROCESSED.search(product):
        return None
    for s, products in CROP_SECTORS.items():
        if product in products:
            return s
    return 3                                    # vegetables, fruit, nuts (incl. pulses, roots)


def extract(cycle_dir, cache):
    """On-farm emissions per cycle (kg per ha), one row per cycle x emission term."""
    if Path(cache).exists():
        return pd.read_csv(cache)
    rows = []
    for f in sorted(Path(cycle_dir).glob("*.jsonld")):
        m = re.match(r"(.+?)-(.+)-(\d{4})-(\d{4})-\d{8}$", f.stem)
        product, country, start = m.group(1), m.group(2), int(m.group(3))
        c = json.load(open(f, encoding="utf-8"))
        for e in c.get("emissions", []):
            if e.get("methodTier") == "background" or not e.get("value"):
                continue
            rows.append({"cycle": f.stem, "product": product, "country": country, "start": start,
                         "term": e["term"]["@id"], "units": e["term"].get("units"), "value": e["value"][0]})
    df = pd.DataFrame(rows)
    df.to_csv(cache, index=False)
    return df


def tidy(df):
    """Latest period per product x country, crop sector, species and compartment, GHG per term."""
    df = df[df.start == df.groupby(["product", "country"]).start.transform("max")].copy()
    df["sector"] = df["product"].map(crop_sector)
    df = df[df.sector.notna() & ~df.term.str.contains(EXCLUDE)].copy()
    df["sector"] = df.sector.astype(int)
    parts = df.term.str.extract(r"^(?P<species>.+?)To(?P<compartment>[A-Z][A-Za-z]*?)(?=[A-Z]|$)")
    df = df.join(parts).dropna(subset=["species"])           # soil fluxes without a compartment
    fossil_co2 = df.term.str.match(r"co2ToAir(FuelCombustion|UreaHydrolysis|LimeNeutralisation)")
    gwp = np.select([df.species == "n2O", (df.species == "ch4") & df.term.str.contains("FuelCombustion"),
                     df.species == "ch4", fossil_co2], [GWP["n2O"], GWP["ch4_fossil"], GWP["ch4"], GWP["co2"]], 0)
    df["ghg"] = df.value * gwp
    df["flow"] = df.species + " to " + df.compartment.str.lower()
    return df


def ratios(df):
    """R[flow, sector, pool] as a long table; pools are HESTIA country slugs plus 'world'."""
    per_cycle = df.groupby(["sector", "country", "cycle", "flow"]).value.sum().unstack("flow", fill_value=0)
    ghg = df.groupby(["sector", "country", "cycle"]).ghg.sum()
    num = per_cycle.groupby(["sector", "country"]).sum()
    den = ghg.groupby(["sector", "country"]).sum()
    # world pool: all national cycles, plus world cycles of crops without national data
    world_cycles = df[(df.country != "world")
                      | ~df["product"].isin(df.loc[df.country != "world", "product"])]
    wnum = world_cycles.groupby(["sector", "flow"]).value.sum().unstack("flow", fill_value=0)
    wden = world_cycles.groupby("sector").ghg.sum()
    R = num.div(den, axis=0)
    Rw = wnum.div(wden, axis=0)
    Rw.index = pd.MultiIndex.from_arrays([Rw.index, ["world"] * len(Rw)], names=["sector", "country"])
    return pd.concat([R.drop(index="world", level="country", errors="ignore"), Rw]).fillna(0)


def build(R, countries, Q_ghg):
    """F (flows x 189*163) for the crop sectors, and the pool used per country x sector."""
    flows = list(R.columns)
    slug = {c: c for c in R.index.get_level_values("country")}
    name_to_slug = {}
    for s in slug:
        name_to_slug[COUNTRY_ALIASES.get(s, s.replace("-", " ").title())] = s
    rows, cols, vals, used = [], [], [], []
    for ci, c in enumerate(countries):
        for s in sorted(CROP_SECTORS) + [3]:
            key = (s, name_to_slug.get(c, "world"))
            if key not in R.index:
                key = (s, "world")
            if key not in R.index:
                continue
            q = Q_ghg[ci * N_SECTORS + s]
            if q <= 0:
                continue
            r = R.loc[key].to_numpy()
            nz = np.nonzero(r)[0]
            rows += list(nz)
            cols += [ci * N_SECTORS + s] * len(nz)
            vals += list(r[nz] * q)
            used.append({"country": c, "sector": s, "pool": key[1]})
    F = sp.csr_matrix((vals, (rows, cols)), shape=(len(flows), len(countries) * N_SECTORS))
    return F, flows, pd.DataFrame(used)


def sector_totals(F, select, n_countries, sectors=range(8)):
    """Per sector (default: crop sectors 0..7) sum of the selected rows over all countries, kg."""
    x = np.asarray(F[select].sum(axis=0)).ravel().reshape(n_countries, N_SECTORS)
    return x[:, list(sectors)].sum(axis=0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hestia", required=True, help="HESTIA release folder with Cycle/*.jsonld")
    ap.add_argument("--rex-dir", required=True)
    ap.add_argument("--exio", required=True, help="EXIOBASE satellite F.txt")
    ap.add_argument("--map-dir", required=True)
    ap.add_argument("--meta", required=True, help="meta.json of an existing extension run")
    ap.add_argument("--compare", nargs="*", default=[], help="extension run folders to compare (F_REX.npz)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    meta = json.load(open(args.meta))
    countries = meta["countries"]
    df = tidy(extract(Path(args.hestia) / "Cycle", Path(args.hestia) / "onfarm_emissions.csv"))
    R = ratios(df)
    R.to_csv(out / "hestia_ratios.csv")
    Q = load_q_rex(args.rex_dir, len(meta["q_labels"]), len(countries) * N_SECTORS)
    F, flows, used = build(R, countries, Q[1])
    sp.save_npz(out / "F_REX_hestia_crops.npz", F)
    pd.DataFrame({"flow": flows}).to_csv(out / "F_REX_hestia_crops_flows.csv", index=False)
    used.to_csv(out / "pools_used.csv", index=False)

    dicts = pickle.load(open(Path(args.map_dir) / "sector_matching_dicts.pickle", "rb"))
    _, exio = load_exiobase(args.exio, dicts["exio_name_dict"], meta["exio_codes"])
    sectors = meta["sectors"][:8]
    table = []
    for p in POLLUTANTS:
        sp_name = [k for k, v in SPECIES.items() if v == p]
        row = {"pollutant": p, "EXIOBASE": exio[p][:, :8].sum(axis=0)}
        if sp_name:
            row["HESTIA"] = sector_totals(F, [i for i, f in enumerate(flows) if f == f"{sp_name[0]} to air"],
                                        len(countries))
        for run in args.compare:
            fl = pd.read_csv(Path(run) / "F_REX_flows.csv")
            Fr = sp.load_npz(Path(run) / "F_REX.npz").tocsr()
            sel = np.where(fl.compartment.str.startswith("air") & fl.name.str.match(POLLUTANTS[p][1], case=False))[0]
            row[Path(run).name.replace("out_", "").replace("_final", "") + " raw"] = sector_totals(Fr, sel, len(countries))
        table.append(row)
    long = []
    for row in table:
        for k, v in row.items():
            if k != "pollutant":
                for s, val in enumerate(v):
                    long.append({"pollutant": row["pollutant"], "source": k, "sector": sectors[s], "kg": val})
    long = pd.DataFrame(long)
    long.to_csv(out / "crop_sector_comparison.csv", index=False)
    tot = long.groupby(["pollutant", "source"]).kg.sum().unstack("source") / 1e9
    tot = tot.loc[[p for p in POLLUTANTS if p in tot.index]]
    tot.to_csv(out / "crop_sector_totals_Mt.csv")
    pd.set_option("display.width", 200)
    print("crop sectors 0-7, global totals [Mt]:")
    print(tot.to_string(float_format="{:.3g}".format))
    print("\npools:", used.pool.eq("world").map({True: "world", False: "country"}).value_counts().to_dict())


if __name__ == "__main__":
    main()
