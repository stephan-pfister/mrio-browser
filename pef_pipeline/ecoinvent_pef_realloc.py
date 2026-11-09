"""PEF re-allocation of ecoinvent 3.9.1 cut-off dairy processes, EF 3.1 results for the markets of
agricultural, food and beverage products (ISIC 01, 10, 11; this includes seed, field services and fibres).

ecoinvent allocates multi-output activities economically and keeps, for every co-product
dataset, its production volume (PV) and properties (price, dry mass). With
AF_i = key_i * PV_i / sum(key * PV), the unallocated inventory is x_i * PV_i / AF_i for every
co-product i, so a different allocation key only rescales each co-product column:

    x'_i = x_i * AF_new_i / AF_old_i

PEF rules applied (Recommendation (EU) 2021/2279, Annex I 4.5.1; dairy PEFCR, February 2025, 5.8.3):
  idf       : dairy farm, AF_milk = 1 - 6.04 * M_cattle / M_milk (masses from PV)
  dry mass  : dairy plant co-products (cheese/whey, yogurt/cream, butter/buttermilk/skimmed milk)

The columns are rescaled in the technosphere and biosphere matrices of one LCA; ecoinvent
itself is not changed or copied.
"""

import csv
import re
from collections import defaultdict
from pathlib import Path

import bw2calc as bc
import bw2data as bd
import numpy as np
from bw2data.backends import ActivityDataset as AD, ExchangeDataset as ED

PROJECT, DB = "salton sea pef", "cutoff391"  # Brightway 2.5 copy of "salton sea"
HERE = Path(__file__).parent
OUT = HERE.parent / "reports" / "ecoinvent_pef_realloc.csv"  # licensed data: stays local
IDF_COEFF = 6.04
RULES = {
    "milk production, from cow": "idf",
    "cheese production, soft, from cow milk": "dry mass",
    "yogurt production, from cow milk": "dry mass",
    "butter production, from cream, from cow milk": "dry mass",
    "butter production, from cow milk": "dry mass",
}
# EF 3.1 category names in ecoinvent -> names in ef31_nw.csv
CATEGORIES = {
    "acidification": "Acidification",
    "climate change": "Climate change",
    "ecotoxicity: freshwater": "Ecotoxicity, freshwater",
    "energy resources: non-renewable": "Resource use, fossils",
    "eutrophication: freshwater": "Eutrophication, freshwater",
    "eutrophication: marine": "Eutrophication marine",
    "eutrophication: terrestrial": "Eutrophication, terrestrial",
    "human toxicity: carcinogenic": "Human toxicity, cancer",
    "human toxicity: non-carcinogenic": "Human toxicity, non-cancer",
    "ionising radiation: human health": "Ionising radiation, human health",
    "land use": "Land use",
    "material resources: metals/minerals": "Resource use, minerals and metals",
    "ozone depletion": "Ozone depletion",
    "particulate matter formation": "EF-particulate Matter",
    "photochemical oxidant formation: human health": "Photochemical ozone formation - human health",
    "water use": "Water use",
}
FOOD_ISIC = re.compile(r"^(01[0-9]{2}|10[0-9]{2}|11[0-9]{2})")


def co_product_groups():
    """{(name, location): [{id, product, p, pv, price, dry mass}]} for the RULES activities."""
    groups = defaultdict(list)
    for a in AD.select().where((AD.database == DB) & (AD.name.in_(list(RULES)))):
        prod = ED.get((ED.output_database == DB) & (ED.output_code == a.code) & (ED.type == "production"))
        props = prod.data.get("properties") or {}
        groups[(a.name, a.location)].append({
            "id": a.id, "product": a.product, "p": prod.data["amount"],
            "pv": prod.data.get("production volume"),
            "price": props.get("price", {}).get("amount"),
            "dry mass": props.get("dry mass", {}).get("amount")})
    return groups


def normalise(values):
    v = np.array(values, dtype=float)
    return v / v.sum()


def shares(members, key):
    """Allocation shares by a property per unit of product times the production volume."""
    for m in members:
        if not m.get(key) or not m.get("pv") or m[key] < 0 or m["pv"] < 0:
            raise ValueError(f"{m['product']}: missing or negative '{key}' or production volume")
    return normalise([m[key] * m["pv"] for m in members])


def allocation_changes(groups):
    """Per dataset id: factor k = AF_new / AF_old, plus a log for the report."""
    k, log = {}, []
    for (name, loc), members in sorted(groups.items()):
        af_old = shares(members, "price")
        if RULES[name] == "idf":
            milk = next(m for m in members if m["product"] == "cow milk")
            others = [m for m in members if m is not milk]
            bmr = sum(m["pv"] for m in others) / milk["pv"]
            af_milk = 1 - IDF_COEFF * bmr
            it = iter(normalise([m["pv"] for m in others]) * (1 - af_milk))  # live animals by mass
            af_new = np.array([af_milk if m is milk else next(it) for m in members])
        else:
            af_new = shares(members, RULES[name])
        if af_new.min() < 0 or af_new.max() > 1 or abs(af_new.sum() - 1) > 1e-9:
            raise ValueError(f"{name} {loc}: invalid PEF allocation factors {af_new}")
        for m, o, n in zip(members, af_old, af_new):
            k[m["id"]] = n / o
            log.append({"activity": name, "location": loc, "product": m["product"], "rule": RULES[name],
                        "AF_economic": o, "AF_PEF": n, "factor": n / o})
    return k, log


def check_economic(groups):
    """Confirm ecoinvent's allocation is the price key: x_i * PV_i / AF_i equal for all co-products."""
    worst = 0.0
    for members in groups.values():
        af = shares(members, "price")
        inputs = []
        for m, a in zip(members, af):
            code = AD.get(AD.id == m["id"]).code
            amounts = {e.input_code: e.data["amount"] for e in
                       ED.select().where((ED.output_database == DB) & (ED.output_code == code) & (ED.type == "technosphere"))}
            inputs.append({c: x * m["pv"] / a for c, x in amounts.items()})
        common = set.intersection(*(set(i) for i in inputs))
        for c in common:
            vals = np.array([i[c] for i in inputs])
            if abs(vals).max() > 0:
                worst = max(worst, float(np.ptp(vals) / abs(vals).max()))
    return worst


def food_markets():
    acts = []
    for a in AD.select().where((AD.database == DB) & (AD.name.startswith("market for"))):
        isic = [c[1] for c in a.data.get("classifications", []) if c[0].startswith("ISIC")]
        if isic and FOOD_ISIC.match(isic[0]):
            acts.append(a)
    return acts


def rescale_columns(matrix, cols, factors, production=None):
    """Multiply columns by factors; ``production`` = [(row, amount)] per column keeps the production
    entry (self-consumption in the same cell is rescaled)."""
    m = matrix.tocsc(copy=True)
    for j, f in zip(cols, factors):
        start, end = m.indptr[j], m.indptr[j + 1]
        m.data[start:end] *= f
    if production is not None:
        for j, f, (i, p) in zip(cols, factors, production):
            m[i, j] = p - f * (p - matrix[i, j])
    return m.tocsr()


def main():
    bd.projects.set_current(PROJECT)
    methods = {CATEGORIES[m[1]]: m for m in bd.methods if m[0] == "EF v3.1" and m[1] in CATEGORIES}
    assert len(methods) == 16, sorted(methods)
    groups = co_product_groups()
    deviation = check_economic(groups)
    if deviation > 0.01:
        raise ValueError(f"ecoinvent allocation is not the price key (deviation {deviation:.2%})")
    print(f"{len(groups)} multi-output dairy activities; ecoinvent allocation = price key "
          f"(max deviation {deviation:.2%})")
    k, log = allocation_changes(groups)
    for row in log:
        print(f"  {row['activity'][:40]:<40} {row['location']:<6} {row['product'][:38]:<38} "
              f"AF {row['AF_economic']:.3f} -> {row['AF_PEF']:.3f}  x{row['factor']:.2f}")

    markets = food_markets()
    fus = [(a.id, a.name, a.location, a.product) for a in markets] + \
          [(m["id"], name, loc, m["product"]) for (name, loc), ms in groups.items() for m in ms]
    print(f"{len(markets)} agricultural, food and beverage markets (ISIC 01, 10, 11) + {len(k)} re-allocated datasets")

    lca = bc.LCA({fus[0][0]: 1}, methods["Climate change"])
    lca.lci()
    cvec = {}
    for c, m in methods.items():
        lca.switch_method(m)
        cvec[c] = lca.characterization_matrix.diagonal()

    def run():
        out = {}
        for fu in fus:
            lca.redo_lci({fu[0]: 1})
            g = lca.biosphere_matrix @ lca.supply_array
            out[fu[0]] = {c: float(v @ g) for c, v in cvec.items()}
        return out

    base = run()
    cols = [lca.dicts.activity[j] for j in k]
    f = [k[j] for j in k]
    prod = [(lca.dicts.product[j], next(m["p"] for ms in groups.values() for m in ms if m["id"] == j)) for j in k]
    lca.technosphere_matrix = rescale_columns(lca.technosphere_matrix, cols, f, production=prod)
    lca.biosphere_matrix = rescale_columns(lca.biosphere_matrix, cols, f)
    pef = run()

    rows = [{"activity": n, "location": loc, "product": p, "category": c,
             "unit": bd.methods[methods[c]]["unit"], "ecoinvent": base[i][c], "PEF": pef[i][c]}
            for i, n, loc, p in fus for c in methods]
    OUT.parent.mkdir(exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(OUT.with_name("ecoinvent_pef_realloc_factors.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(log[0]))
        w.writeheader()
        w.writerows(log)

    cc = [(r["PEF"] / r["ecoinvent"], r) for r in rows if r["category"] == "Climate change" and r["ecoinvent"]]
    changed = [x for x in cc if abs(x[0] - 1) > 0.01]
    print(f"\nClimate change: {len(changed)} of {len(cc)} products change by more than 1 % -> {OUT}")
    for ratio, r in sorted(changed, key=lambda x: -np.inf if x[0] < 0 else -abs(np.log(x[0])))[:25]:
        print(f"  {r['activity'][:45]:<45} {r['location']:<6} {r['product'][:35]:<35} "
              f"{r['ecoinvent']:8.3f} -> {r['PEF']:8.3f}  x{ratio:.2f}")


if __name__ == "__main__":
    main()
