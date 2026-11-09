"""PEF pipeline on the 'bafu test' project: Swiss dairy chain on bafu-2026 + EF 3.1.

Writes (or replaces) the foreground database ``pef-dairy-foreground`` in the
project, applies PEF allocation and the CFF, and computes the 16 EF 3.1
categories plus the single score for each co-product. Foreground amounts are
illustrative per-cow / per-batch values; replace them with primary data.
"""

import csv
from pathlib import Path

import bw2calc as bc
import bw2data as bd
import wurst

from pef_allocation import (add_cff_exchanges, cff_energy_disposal, cff_material,
                            fpcm, single_score, split_multioutput)

PROJECT = "bafu test"
BG, BIO, FG = "bafu-2026", "ef-3.1-biosphere", "pef-dairy-foreground"
HERE = Path(__file__).parent
OUT = HERE.parent / "reports" / "pef_bafu_results.csv"  # reports/ stays local (licensed data)

# bafu-2026 codes
GRASS = "7521559b-1ba5-3113-88de-95a1effa5532"      # Grass from meadow intensive IP, at field
SILAGE = "37b896c5-d2f1-3bf9-8e90-60ca7d4dbf14"     # Silage maize IP, at farm
WHEAT = "a469599a-f6b9-3b5a-b90b-eb512ef684f6"      # Wheat IP, at feed mill
MAIZE = "dd8175b2-cdc5-38ed-83d8-c4bbc93f8b71"      # Grain maize IP, at feed mill
SOY = "996170a3-0a18-3951-8a77-4e615f8af698"        # Soybean meal, at oil mill, BR
DIESEL = "dc92a1d0-abfa-33b3-821f-f82f741ea40e"     # Diesel, burned in agricultural machine
ELEC_LV = "ba889c1a-3fae-3489-94ef-97a4b6f552ef"    # Electricity, low voltage, supply mix GO, CH
ELEC_MV = "baf52276-8f3a-3683-b0df-77e80f0748dc"    # Electricity, medium voltage, supply mix GO, CH
HEAT = "2fcda415-e506-38d9-9c3e-8ce2c2e7c8f2"       # Heat, natural gas, at industrial furnace 1MW
WATER = "1a2f0384-d903-4f9d-8757-1e78d84366ab"      # Tap water, at user
PET_V = "c7af8832-9551-38fd-a628-363e6fc17c8e"      # PET granulate, bottle grade, RER
PET_REC = "e30a7ffc-9ec5-4676-a6e9-73dc77cb0d7a"    # xx Recycling PET, RER: obsolete and empty (no inputs
                                                    # or emissions), the only PET recycling dataset in bafu-2026
PET_INC = "43268043-3dcf-45a5-ac35-45ea38d1412c"    # Disposal, PET, to municipal incineration
PET_LF = "e8828bdf-675f-4add-933f-66d2ed424341"     # Disposal, plastics, mixture, to sanitary landfill

# ef-3.1-biosphere codes (emissions to non-urban air)
CH4_BIO = "fe0acd60-3ddc-11dd-a8eb-0050c2490048"
N2O = "08a91e70-3ddc-11dd-94c6-0050c2490048"
NH3 = "08a91e70-3ddc-11dd-96af-0050c2490048"


def tech(code, amount):
    act = bd.get_node(database=BG, code=code)
    return {"type": "technosphere", "amount": amount, "input": (BG, code), "name": act["name"],
            "product": act["reference product"], "unit": act["unit"], "location": act["location"]}


def link(code, product, amount):
    return {"type": "technosphere", "amount": amount, "input": (FG, code), "name": product,
            "product": product, "unit": "kilogram", "location": "CH"}


def bio(code, amount):
    return {"type": "biosphere", "amount": amount, "input": (BIO, code)}


def prod(product, amount, role=None, **props):
    return {"type": "production", "product": product, "amount": amount, "unit": "kilogram",
            "name": product, "location": "CH", "properties": props, "role": role}


def ds(code, name, rule, exchanges, **kw):
    return dict(database=FG, code=code, name=name, location="CH", unit="kilogram", type="process",
                pef_rule=rule, exchanges=exchanges,
                **{"reference product": exchanges[0]["product"]}, **kw)


def foreground():
    return [
        ds("dairy_farm", "dairy farm, per cow and year", "dairy_farm", [
            prod("raw milk, FPCM", fpcm(7000, fat_pct=4.1, protein_pct=3.4), role="milk"),
            prod("live cattle", 260),  # culled cow + calves sold, live weight
            tech(GRASS, 4500), tech(SILAGE, 1200),  # no current hay dataset (only "xx" ones)
            tech(WHEAT, 300), tech(MAIZE, 300), tech(SOY, 150),
            tech(DIESEL, 80), tech(ELEC_LV, 600), tech(WATER, 25000),
            bio(CH4_BIO, 140), bio(N2O, 1.5), bio(NH3, 25)]),
        ds("slaughterhouse", "cattle slaughtering, per t live weight", "slaughterhouse", [
            prod("fresh meat and edible offal", 490, price=3.0),
            prod("hides and skins", 70, price=0.3),
            prod("food grade bones and fat", 80, price=0.2),
            prod("category 1/2 material", 360, price=0.0),
            link("dairy_farm", "live cattle", 1000),
            tech(ELEC_MV, 90), tech(HEAT, 400), tech(WATER, 1500)]),
        ds("cheese_plant", "cheese production", "physical", [
            prod("cheese", 1, dry_matter=0.60),
            prod("whey", 8, dry_matter=0.065),
            link("dairy_farm", "raw milk, FPCM", 10),
            tech(ELEC_MV, 1.5), tech(HEAT, 4), tech(WATER, 15)], pef_property="dry_matter"),
        ds("cheese_packed", "cheese, packed in PET tray", "physical", [
            prod("cheese, packed", 1),
            link("cheese_plant", "cheese", 1)]),
    ]


def add_packaging(data):
    """20 g PET tray via the CFF. A = 0.5 is the PEF value when Part C of Annex II has no
    application- or material-specific A (check Part C for PET trays); R1 = 0 is the dairy PEFCR
    default for PET cups. R2, R3 and Qs/Qp are placeholders to replace with Part C values.
    Recycling (Erec, ErecEoL) is burden-free because bafu-2026 has no PET recycling inventory."""
    lookup = {(BG, c): {k: a[k] for k in ("name", "unit", "location")} | {"product": a["reference product"]}
              for c in (PET_V, PET_REC, PET_INC, PET_LF, HEAT, ELEC_MV)
              for a in [bd.get_node(database=BG, code=c)]}
    terms = cff_material(A=0.5, R1=0.0, R2=0.4, Qsin_Qp=0.9, Qsout_Qp=0.9,
                         Ev=(BG, PET_V), Erec=(BG, PET_REC), ErecEoL=(BG, PET_REC))
    terms += cff_energy_disposal(B=0, R2=0.4, R3=0.5, LHV=22.0, X_heat=0.3, X_elec=0.15,
                                 Eer=(BG, PET_INC), ESE_heat=(BG, HEAT), ESE_elec=(BG, ELEC_MV),
                                 ED=(BG, PET_LF))
    packed = next(d for d in data if d["code"] == "cheese_packed")
    add_cff_exchanges(packed, mass=0.02, terms=terms, lookup=lookup)


def load_nw():
    with open(HERE / "ef31_nw.csv", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return ({r["category"]: float(r["normalisation_per_person"]) for r in rows},
            {r["category"]: float(r["weighting_pct"]) / 100 for r in rows})


def main():
    bd.projects.set_current(PROJECT)
    if FG in bd.databases:
        del bd.databases[FG]

    data, log = split_multioutput(foreground())
    add_packaging(data)
    for d in data:
        d.pop("pef_rule", None), d.pop("pef_property", None)
    wurst.write_brightway2_database(data, FG)

    print("\nAllocation factors applied:")
    for row in log:
        print(f"  {row['dataset']:<38} {row['product']:<30} {row['rule']:<15} {row['factor']:.3f}")

    nf, wf = load_nw()
    methods = [("sentier", "EF v3.1", c) for c in nf]
    fus = ["cheese_packed", "cheese_plant|whey", "slaughterhouse|fresh meat and edible offal",
           "dairy_farm|raw milk, FPCM", "dairy_farm|live cattle"]
    acts = [bd.get_node(database=FG, code=c) for c in fus]

    lca = bc.LCA({acts[0]: 1}, methods[0])
    lca.lci()
    results = {}
    for act in acts:
        lca.redo_lci({act.id: 1})
        scores = {}
        for m in methods:
            lca.switch_method(m)
            lca.lcia()
            scores[m[2]] = lca.score
        results[act["reference product"]] = scores

    OUT.parent.mkdir(exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["category", "unit"] + list(results))
        for c in nf:
            w.writerow([c, bd.methods[("sentier", "EF v3.1", c)]["unit"]] + [results[p][c] for p in results])
        w.writerow(["Single score", "mPt"] + [1000 * sum(single_score(r, nf, wf).values()) for r in results.values()])
    print(f"\nPer kg product (full table: {OUT}):")
    print(f"  {'product':<30} {'kg CO2 eq':>10} {'m3 depr.':>9} {'mPt':>8}")
    for p, r in results.items():
        print(f"  {p:<30} {r['Climate change']:>10.3f} {r['Water use']:>9.3f} "
              f"{1000 * sum(single_score(r, nf, wf).values()):>8.4f}")


if __name__ == "__main__":
    main()
