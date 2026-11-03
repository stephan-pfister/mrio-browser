"""Standalone HESTIA crop database on BAFU, with PEF allocation (Brightway 2.5 + wurst).

Builds ``hestia-2026-04-30`` in the ``bafu test`` project from the HESTIA aggregated
release (CC-BY 4.0, s3://hestia-aggregated-data): one dataset per crop (or processed crop)
and country, latest period, linked only to ``bafu-2026`` and ``ef-3.1-biosphere`` -- no
ecoinvent anywhere in the system.

Per HESTIA cycle (functional unit: 1 ha, or the processed amount):

- products with an economic value share become production exchanges; multi-output cycles
  (grain + straw, oil + meal, rice milling fractions) are split with
  ``pef_allocation.split_multioutput`` under the economic rule ("crop", LEAP / PEF Annex IV)
- inputs are linked to bafu-2026 (fertilisers per nutrient, fuels as supply only, grid
  electricity, machinery, pesticides, lime, tap water, processing chemicals); seed and crop
  inputs of processing cycles are linked to HESTIA datasets; manure and compost carry no
  burden (residues, PEF cut-off)
- on-farm emissions are mapped to EF 3.1 flows; HESTIA background emissions are dropped
  (BAFU supplies the upstream chains)
- pesticides follow the PEF default: 90 % of the applied active ingredient to agricultural
  soil, 9 % to air, 1 % to fresh water (HESTIA's own pesticide fate is not used)
- land occupation (site + unused duration) as Arable / Permanent crops, irrigation water as
  ground, river or lake water withdrawal (taken as consumed)

Inputs and emissions that cannot be mapped are logged to ``reports/hestia_unmapped.csv``.

    python hestia_db.py [--cycles <HESTIA release>/Cycle] [--score]
"""

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

import bw2calc as bc
import bw2data as bd
import wurst

from pef_allocation import single_score, split_multioutput

PROJECT = "bafu test"
BG, BIO, FG = "bafu-2026", "ef-3.1-biosphere", "hestia-2026-04-30"
HERE = Path(__file__).parent
CYCLES = Path("C:/Users/Pfister/rex_cache/hestia/HESTIA_2026-04-30/Cycle")
REPORTS = HERE.parent / "reports"  # stays local (numbers derived from BAFU)

ISO2 = {"albania": "AL", "argentina": "AR", "australia": "AU", "bangladesh": "BD", "bolivia": "BO",
        "brazil": "BR", "cambodia": "KH", "canada": "CA", "china": "CN", "colombia": "CO",
        "costa-rica": "CR", "cote-divoire": "CI", "cyprus": "CY", "dominican-republic": "DO",
        "ecuador": "EC", "france": "FR", "germany": "DE", "ghana": "GH", "hungary": "HU", "india": "IN",
        "indonesia": "ID", "iran": "IR", "ireland": "IE", "israel": "IL", "italy": "IT", "kenya": "KE",
        "latvia": "LV", "malaysia": "MY", "mexico": "MX", "myanmar": "MM", "nepal": "NP",
        "netherlands": "NL", "pakistan": "PK", "peru": "PE", "portugal": "PT", "romania": "RO",
        "russia": "RU", "south-africa": "ZA", "spain": "ES", "sweden": "SE", "tanzania": "TZ",
        "thailand": "TH", "turkey-country": "TR", "ukraine": "UA", "united-kingdom": "GB",
        "united-states": "US", "viet-nam": "VN", "world": "GLO"}
PROCESSED = re.compile(r"Meal$|Oil(Crude|Refined|Virgin)$|Flour|Bran$|Husk$|Refined|DistillersGrain|^wine$")
PERMANENT = re.compile(r"^(apple|avocado|banana|blueberry|cocoa|coffee|grapes|guava|hazelnut|jujube|kiwifruit"
                       r"|lemon|mandarin|mango|olive|orange|peach|pistachio|almond|oilPalm)", re.I)

# HESTIA input term -> (bafu-2026 name, locations in order of preference, amount factor)
N_AN = ("Ammonium nitrate, as N, at regional storehouse", ["RER"], 1.0)
P_TSP = ("Triple superphosphate, as P2O5, at regional storehouse", ["RER"], 1.0)
K_KCL = ("Potassium chloride, as K2O, at regional storehouse", ["RER"], 1.0)
DAP_N = ("Diammonium phosphate, as N, at regional storehouse", ["RER"], 1.0)
DAP_P = ("Diammonium phosphate, as P2O5, at regional storehouse", ["RER"], 1.0)
LIME = ("Limestone, milled, loose, at plant", ["CH"], 1.0)
LFO = ("Light fuel oil, at regional storage", ["RER"], 1.0)
INPUTS = {
    "inorganicNitrogenFertiliserUnspecifiedKgN": N_AN, "ammoniumNitrateKgN": N_AN, "npkBlendKgN": N_AN,
    "npBlendKgN": N_AN,
    "ureaKgN": ("Urea, as N, at regional storehouse", ["RER"], 1.0),
    "ureaAmmoniumNitrateKgN": N_AN,
    "ammoniumSulphateKgN": ("Ammonium sulphate, as N, at regional storehouse", ["RER"], 1.0),
    "ammoniumSulphateNitrateKgN": ("Ammonium sulphate, as N, at regional storehouse", ["RER"], 1.0),
    "pkBlendKgP2O5": P_TSP, "pkBlendKgK2O": K_KCL,
    "calciumAmmoniumNitrateKgN": ("Calcium ammonium nitrate, as N, at regional storehouse", ["RER"], 1.0),
    "calciumNitrateKgN": ("Calcium nitrate, as N, at regional storehouse", ["RER"], 1.0),
    "potassiumNitrateKgN": ("Potassium nitrate, as N, at regional storehouse", ["RER"], 1.0),
    "diammoniumPhosphateKgN": DAP_N, "monoammoniumPhosphateKgN": DAP_N,
    "anhydrousAmmoniaKgN": ("Ammonia, liquid, at regional storehouse", ["RER"], 17 / 14),
    "inorganicPhosphorusFertiliserUnspecifiedKgP2O5": P_TSP, "npkBlendKgP2O5": P_TSP,
    "npBlendKgP2O5": P_TSP, "tripleSuperPhosphateKgP2O5": P_TSP,
    "singleSuperPhosphateKgP2O5": ("Single superphosphate, as P2O5, at regional storehouse", ["RER"], 1.0),
    "diammoniumPhosphateKgP2O5": DAP_P, "monoammoniumPhosphateKgP2O5": DAP_P,
    "phosphateRockKgP2O5": ("Phosphate rock, as P2O5, beneficiated, dry, at plant", ["MA"], 1.0),
    "inorganicPotassiumFertiliserUnspecifiedKgK2O": K_KCL, "npkBlendKgK2O": K_KCL,
    "potassiumChlorideKgK2O": K_KCL,
    "potassiumSulphateKgK2O": ("Potassium sulphate, as K2O, at regional storehouse", ["RER"], 1.0),
    "potassiumNitrateKgK2O": ("Potassium nitrate, as K2O, at regional storehouse", ["RER"], 1.0),
    "lime": LIME, "agriculturalLimestone": LIME,
    "dolomiticLime": ("Dolomite, at plant", ["RER"], 1.0),
    "burntLime": ("Lime, hydrated, loose, at plant", ["CH"], 74 / 56),       # CaO -> Ca(OH)2
    "calciumHydroxideProcessingAid": ("Lime, hydrated, loose, at plant", ["CH"], 1.0),
    "diesel": ("Diesel, at regional storage", ["RER"], 1.0),
    "petrol": ("Petrol, at regional storage", ["RER"], 1.0),
    "liquefiedPetroleumGasesLpg": ("Propane, butane, at refinery", ["RER"], 1.0),
    "propane": ("Propane, butane, at refinery", ["RER"], 1.0),
    "heatingOil": LFO, "fuelOil": LFO, "heavyGasOil": LFO,
    "lubricant": ("Lubricating oil, at plant", ["RER"], 1.0),
    "naturalGas": ("Natural gas, high pressure, at consumer", ["RoW", "GLO"], 36.0),  # m3 -> MJ (LHV)
    "machineryInfrastructureDepreciatedAmountPerCycle": (
        "Agricultural machinery, general, production", ["CH"], 1.0),
    "steelUnspecifiedDepreciatedAmountPerCycle": ("Steel, low-alloyed, at plant", ["RER"], 1.0),
    "plasticTypeUnspecified": ("Polyethylene, LDPE, granulate, at plant", ["RER"], 1.0),
    "plasticTypeUnspecifiedDepreciatedAmountPerCycle": ("Polyethylene, LDPE, granulate, at plant", ["RER"], 1.0),
    "waterMunicipalTap": ("Tap water, at user", ["RER"], 1000.0),                 # m3 -> kg
    "sodiumHydroxideProcessingAid": ("Sodium hydroxide, 50% in H2O, production mix, at plant", ["RER"], 1.0),
    "CAS-110-54-3": ("Hexane, at plant", ["RER"], 1.0),
    "phosphoricAcidProcessingAid": ("Phosphoric acid, industrial grade, 85% in H2O, at plant", ["RER"], 1.0),
    "sulphuricAcidProcessingAid": ("Sulphuric acid, liquid, at plant", ["RER"], 1.0),
    "liquidAmmoniaProcessingAid": ("Ammonia, liquid, at regional storehouse", ["RER"], 1.0),
    "inorganicMagnesiumFertiliserUnspecifiedKgMg": ("Magnesium sulphate, at plant", ["RER"], 120.4 / 24.3),
    "boron": ("Borax, anhydrous, powder, at plant", ["RER"], 201.2 / 43.2),
}
MIN_MASS_SHARE = 0.005  # co-products below this mass relative to the main product are dropped
PESTICIDE = ("Pesticide unspecified, at regional storehouse", ["RER"], 1.0)
ELECTRICITY = re.compile(r"^electricityGrid")
NO_BURDEN = re.compile(r"manure|compost|slurry|excreta|organicFertiliser|digestate|residue", re.I)
WATER_RESOURCE = {"waterPumpedGroundwater": "Ground Water", "waterRiverStream": "river water",
                  "waterLakeOrReservoir": "lake water", "waterSourceUnspecified": "Water"}

AIR = ("Emissions", "Emissions to air", "Emissions to non-urban air or from high stacks")
AIR_UNSPEC = ("Emissions", "Emissions to air", "Emissions to air, unspecified")
FRESH = ("Emissions", "Emissions to water", "Emissions to fresh water")
AGSOIL = ("Emissions", "Emissions to soil", "Emissions to agricultural soil")
LUC = re.compile(r"LandUseChange|ManagementChange|StockChange|NaturalVegetationBurning|OrganicSoil")
PEF_PESTICIDE = [(AGSOIL, 0.90), (AIR, 0.09), (FRESH, 0.01)]


def emission_flow(term):
    """EF 3.1 (name, compartment) for an on-farm HESTIA emission term, or None."""
    sp, comp = re.match(r"^(.+?)To([A-Z][a-z]+)", term).groups() if re.match(r"^.+?To[A-Z]", term) else (term, "")
    fuel = "FuelCombustion" in term
    luc = bool(LUC.search(term))
    if comp == "Air":
        name = {"nh3": "Ammonia", "n2O": "Nitrous Oxide", "nox": "Nitrogen Oxides", "so2": "Sulfur Dioxide",
                "nmvoc": "Non-methane Volatile Organic Compounds", "pm25": "Particles (PM2.5)",
                "h2S": "Hydrogen Sulfide",
                "ch4": "Methane (land Use Change)" if luc else "Methane (fossil)" if fuel else "Methane (biogenic)",
                "co": "Carbon Monoxide (land Use Change)" if luc else
                      "Carbon Monoxide (fossil)" if fuel else "Carbon Monoxide (biogenic)",
                "co2": "Carbon Dioxide (land Use Change)" if luc else "Carbon Dioxide (fossil)"}.get(sp)
        return (name, AIR) if name else None
    if sp == "no3" and comp in ("Groundwater", "Surface"):
        return ("Nitrate", FRESH)
    if sp == "p" and comp in ("Drainage", "Groundwater", "Surface"):
        return ("Phosphorus", FRESH)
    return None


def load_cycles(cycle_dir, cache):
    """Compact copy of the latest-period cycles (the release files carry large covariance matrices)."""
    if cache.exists():
        return json.load(open(cache, encoding="utf-8"))
    latest = {}
    for f in sorted(Path(cycle_dir).glob("*.jsonld")):
        product, country, start = re.match(r"(.+?)-(.+)-(\d{4})-\d{4}-\d{8}$", f.stem).groups()
        if int(start) >= latest.get((product, country), (0, None))[0]:
            latest[(product, country)] = (int(start), f)
    out = []
    for (product, country), (_, f) in sorted(latest.items()):
        c = json.load(open(f, encoding="utf-8"))
        term = lambda x: {"id": x["term"]["@id"], "name": x["term"].get("name", x["term"]["@id"]),
                          "units": x["term"].get("units"), "type": x["term"].get("termType"),
                          "value": (x.get("value") or [0])[0]}
        out.append({"id": f.stem, "product": product, "country": country,
                    "period": f"{c.get('startDate', '')[:4]}-{c.get('endDate', '')[:4]}",
                    "siteArea": c.get("siteArea", 1), "siteDuration": c.get("siteDuration") or 365,
                    "siteUnusedDuration": c.get("siteUnusedDuration") or 0,
                    "products": [dict(term(p), share=p.get("economicValueShare"), primary=p.get("primary"))
                                 for p in c.get("products", [])],
                    "inputs": [term(i) for i in c.get("inputs", [])],
                    "emissions": [term(e) for e in c.get("emissions", [])
                                  if e.get("methodTier") != "background" and not e["term"]["@id"].startswith("pesticide")]})
    json.dump(out, open(cache, "w", encoding="utf-8"))
    return out


class Lookup:
    """bafu-2026 activities by (name, location) and EF 3.1 flows by (name, compartment) or (CAS, compartment)."""

    def __init__(self):
        self.bg = {(a["name"], a["location"]): a for a in bd.Database(BG)}
        self.elec = {a["location"]: a for (n, _), a in self.bg.items() if n == "Electricity mix"}
        self.bio, self.cas = {}, {}
        for f in bd.Database(BIO):
            cat = tuple(f.get("categories") or ())
            self.bio.setdefault((f["name"].lower(), cat), f["code"])
            if f.get("CAS number"):
                self.cas.setdefault((f["CAS number"].lstrip("0"), cat), f["code"])

    def tech(self, name, locations, amount):
        for loc in locations:
            a = self.bg.get((name, loc))
            if a is not None:
                return {"type": "technosphere", "amount": amount, "input": (BG, a["code"]), "name": a["name"],
                        "product": a.get("reference product"), "unit": a["unit"], "location": a["location"]}
        raise KeyError(f"bafu-2026 has no '{name}' in {locations}")

    def flow(self, name, *cats):
        for cat in cats:
            code = self.bio.get((name.lower(), cat))
            if code:
                return code
        return None


def build(cycles, lk):
    """wurst-format datasets (multi-output where HESTIA gives several valued products) and an unmapped log."""
    by_product = defaultdict(dict)
    for c in cycles:
        by_product[c["product"]][c["country"]] = c
    unmapped = defaultdict(lambda: [0, 0.0, ""])
    data = []

    def miss(term, value, units):
        unmapped[term][0] += 1
        unmapped[term][1] += value
        unmapped[term][2] = units

    for c in cycles:
        loc = ISO2.get(c["country"], "GLO")
        code = c["id"]
        processing = bool(PROCESSED.search(c["product"]))
        products = [p for p in c["products"] if p["share"]] or [p for p in c["products"] if p["primary"]]
        main = next((p for p in products if p["primary"]), products[0])
        # aggregated cycles can carry intercrops with near-zero averaged amounts but a value share;
        # per kg they would get absurd burdens, so they are dropped (main product keeps the burden)
        products = [p for p in products if p is main or p["value"] >= MIN_MASS_SHARE * main["value"]]
        name = f"{main['name']}, {c['country']}, HESTIA {c['period']}"
        ex = [{"type": "production", "product": p["name"], "amount": p["value"], "unit": "kilogram",
               "name": name, "location": loc, "input": (FG, code),
               "properties": {"price": (p["share"] or 100) / p["value"] if p["value"] else 0}}
              for p in products if p["value"]]
        for i in c["inputs"]:
            t, v = i["id"], i["value"]
            if not v:
                continue
            if t in INPUTS:
                ex.append(lk.tech(*INPUTS[t][:2], v * INPUTS[t][2]))
            elif ELECTRICITY.match(t):
                a = lk.elec.get(loc) or lk.elec["GLO"]
                ex.append({"type": "technosphere", "amount": v, "input": (BG, a["code"]), "name": a["name"],
                           "product": a.get("reference product"), "unit": a["unit"], "location": a["location"]})
            elif i["type"] == "pesticideAI":
                ex.append(lk.tech(*PESTICIDE[:2], v))
                cas = t[4:].lstrip("0") if t.startswith("CAS-") else None
                flows = [(lk.cas.get((cas, comp)), share) for comp, share in PEF_PESTICIDE] if cas else []
                if flows and all(f for f, _ in flows):
                    ex += [{"type": "biosphere", "amount": v * s, "input": (BIO, f)} for f, s in flows]
                else:
                    miss(f"{t} (emission)", v, i["units"])
            elif t in WATER_RESOURCE:
                f = lk.flow(WATER_RESOURCE[t], ("Resources", "Resources from water",
                                                "Renewable material resources from water"))
                ex.append({"type": "biosphere", "amount": v, "input": (BIO, f)})
            elif t == "seed" and not processing:
                ex.append({"type": "technosphere", "amount": v, "input": (FG, code), "product": main["name"],
                           "name": name, "unit": "kilogram", "location": loc})
            elif i["type"] == "crop" and t in by_product:
                src = by_product[t].get(c["country"]) or by_product[t].get("world")
                if src is None:
                    miss(t, v, i["units"])
                    continue
                src_main = next((p for p in src["products"] if p["primary"]), src["products"][0])
                ex.append({"type": "technosphere", "amount": v, "input": (FG, src["id"]),
                           "product": src_main["name"], "name": src_main["name"], "unit": "kilogram",
                           "location": ISO2.get(src["country"], "GLO")})
            elif NO_BURDEN.search(t) or i["type"] in ("organicFertiliser", "excreta"):
                continue
            else:
                miss(t, v, i["units"])
        for e in c["emissions"]:
            fl = emission_flow(e["id"])
            code_f = lk.flow(fl[0], fl[1], AIR_UNSPEC) if fl else None
            if code_f:
                ex.append({"type": "biosphere", "amount": e["value"], "input": (BIO, code_f)})
            elif e["value"]:
                miss(e["id"], e["value"], e["units"])
        if not processing:
            irrigated = any(i["id"] in WATER_RESOURCE and i["value"] for i in c["inputs"])
            land = ("Permanent Crops, Non-irrigated, Intensive" if PERMANENT.match(c["product"])
                    else "Arable, Irrigated" if irrigated else "Arable")
            m2a = 1e4 * c["siteArea"] * (c["siteDuration"] + c["siteUnusedDuration"]) / 365
            ex.append({"type": "biosphere", "amount": m2a,
                       "input": (BIO, lk.flow(land, ("Land use", "Land occupation")))})
        data.append({"database": FG, "code": code, "name": name, "location": loc, "unit": "kilogram",
                     "type": "process", "reference product": main["name"], "exchanges": ex,
                     "pef_rule": "economic" if processing else "crop",
                     "comment": f"HESTIA aggregated cycle {code} (CC-BY 4.0); on-farm emissions, "
                                f"background from bafu-2026; PEF allocation, pesticide 90/9/1."})
    return data, unmapped


def score(data):
    """EF 3.1 (16 categories) and single score per kg for every dataset's reference product."""
    with open(HERE / "ef31_nw.csv", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    nf = {r["category"]: float(r["normalisation_per_person"]) for r in rows}
    wf = {r["category"]: float(r["weighting_pct"]) / 100 for r in rows}
    methods = [("sentier", "EF v3.1", c) for c in nf]
    acts = [bd.get_node(database=FG, code=d["code"]) for d in data]
    lca = bc.LCA({acts[0]: 1}, methods[0])
    lca.lci()
    out = []
    for act in acts:
        lca.redo_lci({act.id: 1})
        s = {}
        for m in methods:
            lca.switch_method(m)
            lca.lcia()
            s[m[2]] = lca.score
        out.append({"dataset": act["name"], "product": act["reference product"], **s,
                    "single score mPt": 1000 * sum(single_score(s, nf, wf).values())})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cycles", default=str(CYCLES), help="HESTIA release Cycle folder")
    ap.add_argument("--score", action="store_true", help="also score every product with EF 3.1")
    args = ap.parse_args()

    cycles = load_cycles(args.cycles, Path(args.cycles).parent / "cycles_latest_compact.json")
    bd.projects.set_current(PROJECT)
    if FG in bd.databases:
        del bd.databases[FG]
    lk = Lookup()
    data, unmapped = build(cycles, lk)
    data, log = split_multioutput(data)
    for d in data:
        d.pop("pef_rule", None), d.pop("pef_property", None)
        for e in d["exchanges"]:
            e.pop("properties", None)
    wurst.write_brightway2_database(data, FG)

    REPORTS.mkdir(exist_ok=True)
    with open(REPORTS / "hestia_unmapped.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["term", "cycles", "summed amount", "unit"])
        for t, (n, v, u) in sorted(unmapped.items(), key=lambda x: -x[1][0]):
            w.writerow([t, n, v, u])
    with open(REPORTS / "hestia_allocation.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["dataset", "product", "rule", "factor"])
        w.writeheader()
        w.writerows(log)
    print(f"{FG}: {len(data)} datasets from {len(cycles)} HESTIA cycles; "
          f"{len({r['dataset'] for r in log})} multi-output cycles split; "
          f"{len(unmapped)} unmapped terms (reports/hestia_unmapped.csv)")

    if args.score:
        rows = score(data)
        with open(REPORTS / "hestia_ef31_scores.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"scores: reports/hestia_ef31_scores.csv ({len(rows)} products)")


if __name__ == "__main__":
    main()
