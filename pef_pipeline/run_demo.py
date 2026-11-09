"""End-to-end demo of the PEF pipeline in an isolated Brightway project.

Foreground primary data (multi-output, wurst format) -> PEF allocation + CFF
-> written as single-output Brightway database -> LCIA.
All numbers are illustrative; replace the background with EF datasets /
ecoinvent and the method with EF 3.1 for real studies (see README).
"""

import bw2calc as bc
import bw2data as bd
import wurst

from pef_allocation import (add_cff_exchanges, cff_energy_disposal, cff_material,
                            fpcm, split_multioutput)

PROJECT = "pef_pipeline_demo"
BIO, BG, FG = "demo-biosphere", "demo-background", "demo-foreground"
METHOD = ("demo", "climate change", "GWP100")


def build_background():
    bio = {
        (BIO, "co2"): {"name": "Carbon dioxide, fossil", "unit": "kilogram", "type": "emission", "categories": ("air",)},
        (BIO, "ch4"): {"name": "Methane, biogenic", "unit": "kilogram", "type": "emission", "categories": ("air",)},
        (BIO, "n2o"): {"name": "Dinitrogen monoxide", "unit": "kilogram", "type": "emission", "categories": ("air",)},
    }
    bd.Database(BIO).write(bio)

    def proc(code, name, unit, co2):
        return (BG, code), {"name": name, "reference product": name, "unit": unit, "location": "GLO",
                            "type": "process", "exchanges": [
                                {"input": (BG, code), "amount": 1, "type": "production"},
                                {"input": (BIO, "co2"), "amount": co2, "type": "biosphere"}]}

    bd.Database(BG).write(dict([
        proc("feed", "compound feed", "kilogram", 0.6),
        proc("diesel", "diesel, burned", "kilogram", 3.2),
        proc("elec", "electricity, country mix", "kilowatt hour", 0.12),
        proc("pet_v", "PET granulate, virgin", "kilogram", 2.2),
        proc("pet_r", "PET granulate, recycled", "kilogram", 0.5),
        proc("pet_rec_eol", "PET recycling at end of life", "kilogram", 0.3),
        proc("pet_inc", "PET incineration with energy recovery", "kilogram", 2.0),
        proc("pet_landfill", "PET landfill", "kilogram", 0.05),
        proc("heat", "heat, natural gas", "megajoule", 0.07),
    ]))
    bd.Method(METHOD).write([((BIO, "co2"), 1), ((BIO, "ch4"), 27.0), ((BIO, "n2o"), 273)])


def tech(code, amount):
    """Background technosphere exchange in wurst format (linked by key)."""
    act = bd.get_node(database=BG, code=code)
    return {"type": "technosphere", "amount": amount, "input": (BG, code), "name": act["name"],
            "product": act["reference product"], "unit": act["unit"], "location": act["location"]}


def link(code, product, amount):
    """Foreground exchange to one co-product of a (multi-output) foreground dataset."""
    return {"type": "technosphere", "amount": amount, "input": (FG, code), "name": product,
            "product": product, "unit": "kilogram", "location": "CH"}


def bio(code, amount):
    return {"type": "biosphere", "amount": amount, "input": (BIO, code)}


def prod(product, amount, unit="kilogram", **props):
    role = props.pop("role", None)
    return {"type": "production", "product": product, "amount": amount, "unit": unit,
            "name": product, "location": "CH", "properties": props, "role": role}


def foreground():
    """Primary data as it comes from a farm / plant questionnaire."""
    milk = fpcm(raw_milk_kg=7000, fat_pct=4.1, protein_pct=3.4)
    ds = lambda code, name, rule, exchanges, **kw: dict(
        database=FG, code=code, name=name, location="CH", unit="kilogram", type="process",
        **{"reference product": exchanges[0]["product"]}, pef_rule=rule, exchanges=exchanges, **kw)
    return [
        ds("dairy_farm", "dairy farm, per year and cow", "dairy_farm", [
            prod("raw milk, FPCM", milk, role="milk"),
            prod("live cattle", 300),
            tech("feed", 5000), tech("diesel", 60), tech("elec", 500),
            bio("ch4", 120), bio("n2o", 3)]),
        ds("slaughterhouse", "cattle slaughtering", "slaughterhouse", [
            prod("fresh meat and edible offal", 490, price=3.0),
            prod("hides and skins", 70, price=0.3),
            prod("food grade bones and fat", 80, price=0.2),
            prod("category 1/2 material", 360, price=0.0),
            link("dairy_farm", "live cattle", 1000),
            tech("elec", 90), tech("heat", 400)]),
        ds("cheese_plant", "cheese production", "physical", [
            prod("cheese", 1, dry_matter=0.60),
            prod("whey", 8, dry_matter=0.065),
            link("dairy_farm", "raw milk, FPCM", 10),
            tech("elec", 1.5), tech("heat", 4)], pef_property="dry_matter"),
        ds("cheese_packed", "cheese, packed in PET tray", "physical", [
            prod("cheese, packed", 1),
            link("cheese_plant", "cheese", 1)]),
    ]


def main():
    if PROJECT in bd.projects:
        bd.projects.delete_project(PROJECT, delete_dir=True)
    bd.projects.set_current(PROJECT)
    build_background()

    data, log = split_multioutput(foreground())
    print("Allocation factors applied:")
    for row in log:
        print(f"  {row['dataset']:<32} {row['product']:<30} {row['rule']:<15} {row['factor']:.3f}")

    # Packaging with the Circular Footprint Formula (PET tray, 20 g; parameters illustrative,
    # take A, R1, R2, R3, Qs/Qp, B from Part C of Annex II of the PEF method)
    lookup = {(BG, c): {k: a[k] for k in ("name", "unit", "location")} | {"product": a["reference product"]}
              for c in ("pet_v", "pet_r", "pet_rec_eol", "pet_inc", "pet_landfill", "heat", "elec")
              for a in [bd.get_node(database=BG, code=c)]}
    packed = next(ds for ds in data if ds["code"] == "cheese_packed")
    terms = cff_material(A=0.5, R1=0.3, R2=0.4, Qsin_Qp=0.9, Qsout_Qp=0.9,
                         Ev=(BG, "pet_v"), Erec=(BG, "pet_r"), ErecEoL=(BG, "pet_rec_eol"))
    terms += cff_energy_disposal(B=0, R2=0.4, R3=0.4, LHV=22.0, X_heat=0.3, X_elec=0.15,
                                 Eer=(BG, "pet_inc"), ESE_heat=(BG, "heat"), ESE_elec=(BG, "elec"),
                                 ED=(BG, "pet_landfill"))
    add_cff_exchanges(packed, mass=0.02, terms=terms, lookup=lookup)

    for ds in data:  # strip pipeline-only keys before writing
        ds.pop("pef_rule", None), ds.pop("pef_property", None)
    wurst.write_brightway2_database(data, FG)

    print("\nClimate change, kg CO2-eq per functional unit:")
    for code in ("cheese_packed", "cheese_plant|whey", "slaughterhouse|fresh meat and edible offal",
                 "dairy_farm|raw milk, FPCM"):
        act = bd.get_node(database=FG, code=code)
        lca = bc.LCA({act: 1}, METHOD)
        lca.lci(), lca.lcia()
        print(f"  1 {act['unit']:<9} {act['reference product']:<30} {lca.score:8.3f}")


if __name__ == "__main__":
    main()
