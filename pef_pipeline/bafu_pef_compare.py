"""EF 3.1 results of the BAFU food/feed datasets: BAFU as published vs. a PEF version.

Correction of the Brightway import, applied first (column BAFU; BAFU_raw is without it):
  fix  : unspecified chromium emitted to soil and water is Cr(III), as in BAFU's own EF mapping

PEF layers on top of it (each reported on its own and combined):
  pest : pesticide active ingredients emitted to soil are re-routed
         90 % soil / 9 % air / 1 % fresh water (EF method, agricultural modelling)
  hm   : negative heavy-metal balances to soil are set to zero
  water: water use is characterised with the country factor of each dataset's location
         (JRC EF 3.1 CF file; regions without a country factor keep the global factor)

All layers change only biosphere exchanges or their factors, so their effect is linear:
delta score = sum_j s_j * sum_f dB_fj * CF_f, with s the supply vector of the unchanged BAFU
system. No copy of bafu-2026 is written.

    python bafu_pef_compare.py ["<path to EF-LCIAMethod_CF(EF-v3.1).xlsx>"]
"""

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

import bw2calc as bc
import bw2data as bd
import numpy as np
import openpyxl
from bw2data.backends import ActivityDataset as AD, ExchangeDataset as ED

PROJECT, DB, BIO = "bafu test", "bafu-2026", "ef-3.1-biosphere"
HERE = Path(__file__).parent
OUT = HERE.parent / "reports" / "bafu_pef_compare.csv"  # licensed data: stays local
JRC_CF = Path.home() / "Downloads" / "EF-LCIAMethod_CF(EF-v3.1).xlsx"
CR_REMAP = ("Emissions to soil", "Emissions to water")  # compartments where chromium is Cr(III)
CR3_NAMES = ("chromium(3+)", "chromium (iii)")
SPLIT = {"soil": 0.90, "air": 0.09, "water": 0.01}
AIR = ("Emissions", "Emissions to air", "Emissions to non-urban air or from high stacks")
WATER = ("Emissions", "Emissions to water", "Emissions to fresh water")
AG_SOIL, NONAG_SOIL = "Emissions to agricultural soil", "Emissions to non-agricultural soil"
HEAVY_METAL = re.compile(r"^(cadmium|copper|zinc|lead|nickel|chromium|mercury|arsenic|cobalt|molybdenum|selenium)", re.I)
NOT_PESTICIDE = re.compile(r"^(ammoni|nitrat|nitrogen|phosph|carbon|water|oils|sulfur|sulphur|calcium|potassium"
                           r"|magnesium|sodium|chlorid|fluorid|iron|alumin|boron|manganese|silver|tin|barium"
                           r"|strontium|vanadium|titanium|antimony|thallium|beryllium|uranium|radium|heat)", re.I)
NON_FOOD = re.compile(r"seed|seedling|fibres|machinery|diesel|application|drying|grading|freezing|canning|mineral"
                      r"|biowaste|vinasse|store|husked|dung|housing|parlour|tunnel|greenhouse|silo|shed|oil mill"
                      r"|sugar refinery", re.I)


def food_products():
    acts = [a for a in bd.Database(DB) if a.get("sector") == "agriculture" and a["unit"] == "kilogram"
            and not NON_FOOD.search(a["name"])]
    return sorted(acts, key=lambda a: (a["name"], a["location"]))


def water_factors(path):
    """JRC EF 3.1 water-use factors: {(flow name lower, compartment): {location or None: CF}}."""
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    out = defaultdict(dict)
    for r in wb["lciamethods_CF"].iter_rows(min_row=2, values_only=True):
        if r[3] == "Water use":
            out[(r[1].lower(), r[8])][r[5]] = r[4]
    wb.close()
    return out


def emission_changes(water_cf):
    """Per activity id: {flow id: delta amount} for each layer, plus water-use score deltas (water)."""
    flows = {f.key: f for f in bd.Database(BIO)}
    by_name = defaultdict(dict)
    for f in flows.values():
        cat = tuple(f.get("categories") or ())
        if cat in (AIR, WATER):
            by_name[f["name"].lower()].setdefault(cat, f.id)
    cr3 = {tuple(f["categories"]): f.id for f in flows.values() if f["name"].lower() in CR3_NAMES}
    water = {k: water_cf[(f["name"].lower(), f["categories"][-1])] for k, f in flows.items()
             if f.get("categories") and (f["name"].lower(), f["categories"][-1]) in water_cf}
    ag_ids = {a.id for a in bd.Database(DB) if a.get("sector") == "agriculture"}
    acts = {(a.database, a.code): (a.id, a.location)
            for a in AD.select(AD.id, AD.database, AD.code, AD.location).where(AD.database == DB)}

    layers = {k: defaultdict(lambda: defaultdict(float)) for k in ("fix", "pest", "hm", "water")}
    regional = defaultdict(float)
    pesticides, missing, no_cr3 = set(), set(), set()
    for exc in ED.select().where((ED.output_database == DB) & (ED.type == "biosphere") & (ED.input_database == BIO)):
        key = (exc.input_database, exc.input_code)
        f = flows[key]
        (j, loc), v = acts[(exc.output_database, exc.output_code)], exc.data["amount"]
        cat, fid, name = tuple(f.get("categories") or ()), f.id, f["name"]

        if key in water and water[key].get(loc) is not None:
            regional[j] += v * (water[key][loc] - water[key][None])
        if name.lower() == "chromium" and len(cat) > 1 and cat[1] in CR_REMAP:
            if cat in cr3:
                layers["fix"][j][fid] -= v
                layers["fix"][j][cr3[cat]] += v
                fid = cr3[cat]
            else:
                no_cr3.add(cat[-1])

        if not cat or cat[-1] not in (AG_SOIL, NONAG_SOIL):
            continue
        nonag = cat[-1] == NONAG_SOIL
        if HEAVY_METAL.match(name):
            if v < 0 and not nonag:
                layers["hm"][j][fid] -= v
        elif not NOT_PESTICIDE.match(name) and (not nonag or j in ag_ids):
            targets = by_name[name.lower()]
            if AIR not in targets or WATER not in targets:
                missing.add(name)
                continue
            pesticides.add(name)
            layers["pest"][j][f.id] += (SPLIT["soil"] - 1) * v
            layers["pest"][j][targets[AIR]] += SPLIT["air"] * v
            layers["pest"][j][targets[WATER]] += SPLIT["water"] * v
    return layers, regional, pesticides, missing, no_cr3


def main(jrc_cf=JRC_CF):
    bd.projects.set_current(PROJECT)
    methods = [m for m in bd.methods if m[:2] == ("sentier", "EF v3.1") and "_" not in m[2]
               and m[2] not in ("Climate change-Biogenic", "Climate change-Fossil",
                                "Climate change-Land use and land use change")]
    products = food_products()
    layers, regional, pesticides, missing, no_cr3 = emission_changes(water_factors(jrc_cf))
    print(f"{len(products)} food/feed products, {len(methods)} EF categories, {len(pesticides)} pesticides re-routed"
          f" ({len(missing)} without air/water flow: {sorted(missing)[:10]})")
    print(f"activities changed: fix {len(layers['fix'])}; pest {len(layers['pest'])}, hm {len(layers['hm'])},"
          f" water {len(regional)}; chromium kept (no Cr(III) flow): {sorted(no_cr3)}")

    cfs = {m: {(k if isinstance(k, int) else bd.get_node(database=k[0], code=k[1]).id): v
               for k, v in bd.Method(m).load()} for m in methods}

    lca = bc.LCA({products[0].id: 1}, methods[0])
    lca.lci()
    col = lca.dicts.activity
    # delta vectors over activity columns, per layer and method
    delta = {}
    for layer, changes in layers.items():
        for m in methods:
            d = np.zeros(len(col))
            for j, flows in changes.items():
                if j in col:
                    d[col[j]] = sum(dv * cfs[m].get(f, 0.0) for f, dv in flows.items())
            if layer == "water" and m[2] == "Water use":
                for j, ds in regional.items():
                    if j in col:
                        d[col[j]] += ds
            delta[layer, m] = d
    cvec = {}
    for m in methods:
        lca.switch_method(m)
        cvec[m] = lca.characterization_matrix.diagonal()

    rows = []
    for act in products:
        lca.redo_lci({act.id: 1})
        g = lca.biosphere_matrix @ lca.supply_array
        for m in methods:
            raw = float(cvec[m] @ g)
            base = raw + float(delta["fix", m] @ lca.supply_array)
            dp = float(delta["pest", m] @ lca.supply_array)
            dh = float(delta["hm", m] @ lca.supply_array)
            dw = float(delta["water", m] @ lca.supply_array)
            rows.append({"product": act["name"], "location": act["location"], "category": m[2],
                         "unit": bd.methods[m]["unit"], "BAFU_raw": raw, "BAFU": base, "PEF_pest": base + dp,
                         "PEF_hm": base + dh, "PEF_water": base + dw, "PEF": base + dp + dh + dw})

    OUT.parent.mkdir(exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print(f"\nPEF / BAFU ratio across {len(products)} products (median, min, max) -> {OUT}")
    for m in methods:
        r = [x["PEF"] / x["BAFU"] for x in rows if x["category"] == m[2] and x["BAFU"]]
        changed = sum(abs(v - 1) > 0.01 for v in r)
        print(f"  {m[2]:<45} {np.median(r):8.3f} {min(r):9.3f} {max(r):9.3f}   >1% change: {changed}")


if __name__ == "__main__":
    main(*sys.argv[1:])
