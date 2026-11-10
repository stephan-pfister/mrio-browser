"""EF 3.1 results of the BAFU food/feed datasets: BAFU as published vs. a PEF version.

Correction of the Brightway import, applied first (column BAFU; BAFU_raw is without it):
  fix  : unspecified chromium emitted to soil and water is Cr(III), as in BAFU's own EF mapping

PEF layers on top of it (each reported on its own and combined), Recommendation (EU) 2021/2279 Annex I:
  pest : pesticide active ingredients emitted to soil are re-routed 90 % agricultural soil / 9 % air /
         1 % fresh water (default of 4.4.1.4; a PEFCR may set otherwise, the dairy PEFCR keeps 100 % soil)
  hm   : negative heavy-metal balances to agricultural soil are set to zero, a partial version of 4.4.1.6
         approach (a), under which crop uptake is not modelled at all
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
import uuid
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
AG_SOIL_CAT = ("Emissions", "Emissions to soil", AG_SOIL)
HEAVY_METAL = re.compile(r"^(cadmium|copper|zinc|lead|nickel|chromium|mercury|arsenic|cobalt|molybdenum|selenium)", re.I)
NOT_PESTICIDE = re.compile(r"^(ammoni|nitrat|nitrogen|phosph|carbon|water|oils|sulfur|sulphur|calcium|potassium"
                           r"|magnesium|sodium|chlorides?\b|fluorides?\b|iron|alumin|boron|manganese|silver|tin"
                           r"|barium|strontium|vanadium|titanium|antimony|thallium|beryllium|uranium|radium|heat"
                           r"|silicon|chlorine|bromine|bromide|iodide|iodine)", re.I)
PESTICIDE_ALWAYS = re.compile(r"^sodium methylarsonate$", re.I)  # herbicide caught by the "sodium" prefix
# fungicide (sulfur) and potato haulm desiccant (sulfuric acid) when emitted by agricultural datasets
PESTICIDE_IN_AGRICULTURE = re.compile(r"^(sulfur|sulfuric acid)$", re.I)
# agricultural datasets per kg that are not food or feed: services, inputs, fibres, sowing seed, obsolete
NOT_FOOD = re.compile(r"^(xx |agricultural machinery|application|average mineral fertiliser|biowaste|canning"
                      r"|diesel|dried roughage|drying|freezing|grain drying|maize drying|potato grading"
                      r"|solid manure|slurry|vinasse)|fibres|seedling"
                      r"|\b(barley|clover|grass|maize|pea|potato|rice|rye|sugar beet|wheat) seed\b", re.I)
PRODUCTS = HERE.parent / "reports" / "bafu_food_products.csv"


def food_products():
    """Food and feed products of bafu-2026 (agriculture sector, per kg); the selection goes to PRODUCTS."""
    candidates = sorted((a for a in bd.Database(DB) if a.get("sector") == "agriculture" and a["unit"] == "kilogram"),
                        key=lambda a: (a["name"], a["location"]))
    PRODUCTS.parent.mkdir(exist_ok=True)
    with open(PRODUCTS, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["product", "location", "included"])
        w.writerows([a["name"], a["location"], not NOT_FOOD.search(a["name"])] for a in candidates)
    return [a for a in candidates if not NOT_FOOD.search(a["name"])]


def is_pesticide(name, agricultural):
    """Soil-emitted substance that is not a heavy metal: an active ingredient unless it is inorganic."""
    if PESTICIDE_ALWAYS.match(name) or (agricultural and PESTICIDE_IN_AGRICULTURE.match(name)):
        return True
    return not NOT_PESTICIDE.match(name)


def water_factors(path):
    """JRC EF 3.1 water-use factors: {(flow name lower, compartment): {location or None: CF}}."""
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    out = defaultdict(dict)
    for r in wb["lciamethods_CF"].iter_rows(min_row=2, values_only=True):
        if r[3] == "Water use":
            out[(r[1].lower(), r[8])][r[5]] = r[4]
    wb.close()
    return out


def uuid_version(code):
    try:
        return uuid.UUID(code).version
    except ValueError:
        return None


def pick(candidates, ref):
    """One flow among EF duplicates (same name and compartment), chosen deterministically: same CAS
    number as the soil flow, then the same UUID generation (legacy ILCD vs. newer), the same spelling."""
    return min(candidates, key=lambda c: (c.get("CAS number") != ref.get("CAS number"),
                                          uuid_version(c["code"]) != uuid_version(ref["code"]),
                                          c["name"] != ref["name"], c["code"]))


def flow_lookups(water_cf):
    """EF 3.1 flows by key; pesticide target flows by name and compartment; the Cr(III) flow per
    compartment; JRC water-use factors {location or None: CF} per water flow key."""
    flows = {f.key: f for f in bd.Database(BIO)}
    by_name = defaultdict(lambda: defaultdict(list))
    for f in flows.values():
        cat = tuple(f.get("categories") or ())
        if cat in (AIR, WATER, AG_SOIL_CAT):
            by_name[f["name"].lower()][cat].append(f)
    cr3 = {}
    for f in sorted(flows.values(), key=lambda f: f["code"]):
        if f["name"].lower() in CR3_NAMES:
            cr3.setdefault(tuple(f["categories"]), f)
    water = {k: water_cf[(f["name"].lower(), f["categories"][-1])] for k, f in flows.items()
             if f.get("categories") and (f["name"].lower(), f["categories"][-1]) in water_cf}
    return flows, by_name, cr3, water


def emission_changes(water_cf):
    """Per activity id: {flow id: delta amount} for each layer, plus water-use score deltas (water)."""
    flows, by_name, cr3, water = flow_lookups(water_cf)
    ag_ids = {a.id for a in bd.Database(DB) if a.get("sector") == "agriculture"}
    acts = {(a.database, a.code): (a.id, a.location)
            for a in AD.select(AD.id, AD.database, AD.code, AD.location).where(AD.database == DB)}

    layers = {k: defaultdict(lambda: defaultdict(float)) for k in ("fix", "pest", "hm", "water")}
    regional, hm_net = defaultdict(list), defaultdict(float)  # regional: [(flow id, amount, country CF)]
    pesticides, missing, no_cr3, duplicates = set(), set(), set(), {}
    for exc in ED.select().where((ED.output_database == DB) & (ED.type == "biosphere") & (ED.input_database == BIO)):
        key = (exc.input_database, exc.input_code)
        f = flows[key]
        (j, loc), v = acts[(exc.output_database, exc.output_code)], exc.data["amount"]
        cat, fid, name = tuple(f.get("categories") or ()), f.id, f["name"]

        if key in water and water[key].get(loc) is not None:
            regional[j].append((fid, v, water[key][loc]))
        if name.lower() == "chromium" and len(cat) > 1 and cat[1] in CR_REMAP:
            if cat in cr3:
                layers["fix"][j][fid] -= v
                layers["fix"][j][cr3[cat].id] += v
                fid = cr3[cat].id
            else:
                no_cr3.add(cat[-1])

        if not cat or cat[-1] not in (AG_SOIL, NONAG_SOIL):
            continue
        nonag = cat[-1] == NONAG_SOIL
        if HEAVY_METAL.match(name):
            if not nonag:
                hm_net[(j, fid)] += v
        elif is_pesticide(name, j in ag_ids) and (not nonag or j in ag_ids):
            # field-applied active ingredient: 90 % agricultural soil (also when BAFU used
            # non-agricultural soil), 9 % air, 1 % fresh water
            targets = by_name[name.lower()]
            soil = (pick(targets[AG_SOIL_CAT], f) if AG_SOIL_CAT in targets else None) if nonag else f
            if AIR not in targets or WATER not in targets or soil is None:
                missing.add(name)
                continue
            pesticides.add(name)
            air, wat = pick(targets[AIR], f), pick(targets[WATER], f)
            if len(targets[AIR]) > 1 or len(targets[WATER]) > 1:
                duplicates[(name, f["code"])] = (air["code"], wat["code"])
            layers["pest"][j][f.id] -= v
            layers["pest"][j][soil.id] += SPLIT["soil"] * v
            layers["pest"][j][air.id] += SPLIT["air"] * v
            layers["pest"][j][wat.id] += SPLIT["water"] * v
    for (j, fid), net in hm_net.items():  # per dataset and flow: a negative net balance becomes zero
        if net < 0:
            layers["hm"][j][fid] -= net
    return layers, regional, pesticides, missing, no_cr3, duplicates


def main(jrc_cf=JRC_CF):
    bd.projects.set_current(PROJECT)
    methods = [m for m in bd.methods if m[:2] == ("sentier", "EF v3.1") and "_" not in m[2]
               and m[2] not in ("Climate change-Biogenic", "Climate change-Fossil",
                                "Climate change-Land use and land use change")]
    products = food_products()
    layers, regional, pesticides, missing, no_cr3, duplicates = emission_changes(water_factors(jrc_cf))
    print(f"{len(products)} food/feed products, {len(methods)} EF categories, {len(pesticides)} pesticides re-routed"
          f" ({len(missing)} without air/water flow: {sorted(missing)[:10]})")
    for (name, soil), (air, wat) in sorted(duplicates.items()):
        print(f"  EF duplicates for {name} (soil flow {soil[:8]}): air {air[:8]}, water {wat[:8]}")
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
            if layer == "water" and m[2] == "Water use":  # country factor instead of the method's factor
                for j, entries in regional.items():
                    if j in col:
                        d[col[j]] += sum(v * (cf - cfs[m].get(f, 0.0)) for f, v, cf in entries)
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

    print(f"\nPEF / BAFU across {len(products)} products (median, min, max where BAFU > 0) -> {OUT}")
    for m in methods:
        x = [y for y in rows if y["category"] == m[2] and y["BAFU"]]
        r = [y["PEF"] / y["BAFU"] for y in x if y["BAFU"] > 0]
        changed = sum(abs(y["PEF"] - y["BAFU"]) > 0.01 * abs(y["BAFU"]) for y in x)
        flips = sum((y["PEF"] < 0) != (y["BAFU"] < 0) for y in x)
        stats = f"{np.median(r):8.3f} {min(r):9.3f} {max(r):9.3f}" if r else f"{'':>28}"
        print(f"  {m[2]:<45} {stats}   >1% change: {changed}, sign change: {flips}")


if __name__ == "__main__":
    main(*sys.argv[1:])
