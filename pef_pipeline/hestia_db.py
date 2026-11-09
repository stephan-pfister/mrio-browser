"""Standalone HESTIA crop database on BAFU, with PEF allocation (Brightway 2.5 + wurst).

Builds ``hestia-2026-04-30`` in the ``bafu test`` project from the HESTIA aggregated
release (CC-BY 4.0, s3://hestia-aggregated-data): one dataset per crop (or processed crop)
and country, latest period, linked only to ``bafu-2026`` and ``ef-3.1-biosphere`` -- no
ecoinvent anywhere in the system.

Per HESTIA cycle (functional unit: 1 ha, or the processed amount):

- products with an economic value share become production exchanges; multi-output cycles
  (grain + straw, oil + meal, rice milling fractions) are split with
  ``pef_allocation.split_multioutput`` under the economic rule ("crop", LEAP / PEF Annex I);
  co-products under 0.5 % of the main product's mass and 2 % of the value are dropped
- inputs are linked to bafu-2026 (fertilisers per nutrient, fuels as supply only, grid
  electricity, machinery, pesticides, lime, tap water, processing chemicals); seed, and crop
  or crude-oil inputs of processing cycles, are linked to HESTIA datasets; manure and compost
  carry no burden (residues, PEF cut-off)
- on-farm emissions are mapped to EF 3.1 flows, including P lost with eroded soil. N lost with
  eroded soil is left out: PEF, ecoinvent and Agribalyse count N to water as nitrate leaching,
  and soil-bound N scored as total N would dominate marine eutrophication (cotton, millets).
  HESTIA background emissions are dropped (BAFU supplies the upstream chains). Fuel-combustion
  totals that duplicate per-fuel entries are dropped, CO2 from biofuels is biogenic, and carbon
  gains from land-use change are not credited (PEF: no credit for carbon storage)
- pesticides follow the PEF default: 90 % of the applied active ingredient to agricultural
  soil, 9 % to air, 1 % to fresh water (HESTIA's own pesticide fate is not used). Unidentified
  active ingredients (unspecified, classes) are emitted as the cycle's identified mix where that
  carries most of the mass, else as the mean mix over all cycles; bulk, inert and processing
  substances HESTIA types as ingredients (oils, kaolin, chlorine, sulphur, acids...) are left out
- land occupation (site + unused duration) as Arable / Permanent crops, irrigation water as
  ground, river or lake water withdrawal (canal water as river water; taken as consumed)

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
PROCESSED = re.compile(r"Meal$|Oil(Crude|Refined|Virgin)$|Flour|Bran$|(?<!In)Husk$|Refined|DistillersGrain|^wine$")
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
# co-products below both thresholds (mass relative to the main product, value share in %) are
# dropped: aggregation artefacts such as intercrops with near-zero averaged amounts
MIN_MASS_SHARE, MIN_VALUE_SHARE = 0.005, 2.0
PESTICIDE = ("Pesticide unspecified, at regional storehouse", ["RER"], 1.0)
ELECTRICITY = re.compile(r"^electricityGrid")
NO_BURDEN = re.compile(r"manure|compost|slurry|excreta|organicFertiliser|digestate|residue", re.I)
WATER_RESOURCE = {"waterPumpedGroundwater": "Ground Water", "waterRiverStream": "river water",
                  "waterLakeOrReservoir": "lake water", "waterIrrigationCanal": "river water",
                  "waterSourceUnspecified": "Water"}
# inputs HESTIA types as pesticide ingredients that are bulk, inert or processing substances
NOT_PESTICIDE = re.compile(r"chlorine|thiosul|sodium sul|zinc sul|kaolin|paraffin|mineral oil|petroleum oil"
                           r"|vegetable oil|buttermilk|zeolite|boric|borax|adjuvant|surfactant|sulphur|sulfur"
                           r"|hydroxide|phosphoric|sulphuric|sulfuric|nitrogen|benzoic|ethylene glycol", re.I)
# unidentified active ingredients: the only ones the proxy stands in for
UNIDENTIFIED_AI = re.compile(r"unspecifiedai$|compoundsunspecified$", re.I)
BIOFUEL = re.compile(r"biodiesel|ethanol|biogas|vegetable oil|wood", re.I)

AIR = ("Emissions", "Emissions to air", "Emissions to non-urban air or from high stacks")
AIR_GROUND = ("Emissions", "Emissions to air", "Emissions to non-urban air close to ground")
AIR_UNSPEC = ("Emissions", "Emissions to air", "Emissions to air, unspecified")
FRESH = ("Emissions", "Emissions to water", "Emissions to fresh water")
WATER_UNSPEC = ("Emissions", "Emissions to water", "Emissions to water, unspecified")
AGSOIL = ("Emissions", "Emissions to soil", "Emissions to agricultural soil")
SOIL_UNSPEC = ("Emissions", "Emissions to soil", "Emissions to soil, unspecified")
LUC = re.compile(r"LandUseChange|ManagementChange|StockChange|NaturalVegetationBurning|OrganicSoil")
# PEF default split of the applied active ingredient, with compartments to fall back on
PEF_PESTICIDE = [((AGSOIL, SOIL_UNSPEC), 0.90), ((AIR, AIR_GROUND, AIR_UNSPEC), 0.09), ((FRESH, WATER_UNSPEC), 0.01)]


def emission_flow(term, biofuel=False):
    """EF 3.1 [(name, compartment), ...] candidates for an on-farm HESTIA emission term, or []."""
    if term == "pErosionSoilFlux":                       # P reaching water with eroded soil
        return [("Phosphorus", FRESH), ("Phosphorus", WATER_UNSPEC)]
    sp, comp = re.match(r"^(.+?)To([A-Z][a-z]+)", term).groups() if re.match(r"^.+?To[A-Z]", term) else (term, "")
    fuel = "FuelCombustion" in term and not biofuel
    luc = bool(LUC.search(term))
    if comp == "Air":
        if sp == "pm25":                                 # field and tractor emissions, close to ground
            return [("Particles (PM0.2 - PM2.5)", AIR_GROUND), ("Particles (PM2.5)", AIR), ("Particles (PM2.5)", AIR_UNSPEC)]
        name = {"nh3": "Ammonia", "n2O": "Nitrous Oxide", "nox": "Nitrogen Oxides", "so2": "Sulfur Dioxide",
                "nmvoc": "Non-methane Volatile Organic Compounds", "h2S": "Hydrogen Sulfide",
                "ch4": "Methane (land Use Change)" if luc else "Methane (fossil)" if fuel else "Methane (biogenic)",
                "co": "Carbon Monoxide (land Use Change)" if luc else
                      "Carbon Monoxide (fossil)" if fuel else "Carbon Monoxide (biogenic)",
                "co2": "Carbon Dioxide (land Use Change)" if luc else
                       "Carbon Dioxide (fossil)" if fuel or "FuelCombustion" not in term else "Carbon Dioxide (biogenic)"
                }.get(sp)
        return [(name, AIR), (name, AIR_UNSPEC)] if name else []
    if sp == "no3" and comp in ("Groundwater", "Surface"):
        return [("Nitrate", FRESH), ("Nitrate", WATER_UNSPEC)]
    if sp == "p" and comp in ("Drainage", "Groundwater", "Surface"):
        return [("Phosphorus", FRESH), ("Phosphorus", WATER_UNSPEC)]
    return []


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
                          "value": (x.get("value") or [0])[0],
                          "inputs": [i.get("@id") for i in x.get("inputs") or []]}
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
        for f in sorted(bd.Database(BIO), key=lambda f: f["code"]):     # deterministic where names repeat
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


def pesticide_flows(term, lk):
    """[(EF flow code, share)] of the PEF 90/9/1 split for an active ingredient given by CAS, or []."""
    if not term.startswith("CAS-"):
        return []
    cas = term[4:].lstrip("0")
    flows = [(next((lk.cas[(cas, c)] for c in comps if (cas, c) in lk.cas), None), share)
             for comps, share in PEF_PESTICIDE]
    return flows if all(f for f, _ in flows) else []


def is_pesticide(i, processing):
    """Inputs HESTIA types as pesticide ingredients that are real crop-protection active ingredients."""
    return (i["type"] == "pesticideAI" and not processing
            and not NOT_PESTICIDE.search(i["id"]) and not NOT_PESTICIDE.search(i["name"]))


def farm_emissions(c):
    """On-farm emissions without double counts: where HESTIA gives fuel-combustion emissions per
    fuel input and also as an unlinked total from another model, the unlinked total is dropped."""
    linked = {e["id"] for e in c["emissions"] if "FuelCombustion" in e["id"] and e.get("inputs")}
    return [e for e in c["emissions"] if e.get("inputs") or e["id"] not in linked]


def build(cycles, lk):
    """wurst-format datasets (multi-output where HESTIA gives several valued products) and an unmapped log."""
    by_product = defaultdict(dict)
    for c in cycles:
        by_product[c["product"]][c["country"]] = c
    unmapped = defaultdict(lambda: [0, 0.0, ""])
    data = []
    # global mix of identified active ingredients: mean of the cycles' mass shares, so that no
    # single cycle (e.g. an implausible dose) dominates it
    global_mix, n_mix = defaultdict(float), 0
    for c in cycles:
        proc = bool(PROCESSED.search(c["product"]))
        ident = {i["id"]: i["value"] for i in c["inputs"]
                 if i["value"] and is_pesticide(i, proc) and pesticide_flows(i["id"], lk)}
        if ident:
            n_mix += 1
            for t, v in ident.items():
                global_mix[t] += v / sum(ident.values())
    global_mix = {t: v / n_mix for t, v in global_mix.items()}

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
        # aggregated cycles can carry intercrops with near-zero averaged amounts and a small value
        # share; per kg they would get absurd burdens, so they are dropped (the remaining products
        # share their burden)
        products = [p for p in products if p is main or p["value"] >= MIN_MASS_SHARE * main["value"]
                    or (p["share"] or 0) >= MIN_VALUE_SHARE]
        name = f"{main['name']}, {c['country']}, HESTIA {c['period']}"
        ex = [{"type": "production", "product": p["name"], "amount": p["value"], "unit": "kilogram",
               "name": name, "location": loc, "input": (FG, code),
               "properties": {"price": (p["share"] or 100) / p["value"] if p["value"] else 0}}
              for p in products if p["value"]]
        known, unknown = defaultdict(float), 0.0
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
                if not is_pesticide(i, processing):
                    miss(f"{t} (not a crop-protection ingredient, left out)", v, i["units"])
                    continue
                ex.append(lk.tech(*PESTICIDE[:2], v))
                flows = pesticide_flows(t, lk)
                if flows:
                    ex += [{"type": "biosphere", "amount": v * s, "input": (BIO, f)} for f, s in flows]
                    known[t] += v
                elif UNIDENTIFIED_AI.search(t):
                    unknown += v
                    miss(f"{t} (unidentified, emitted via proxy)", v, i["units"])
                else:
                    miss(f"{t} (no EF flow, emission left out)", v, i["units"])
            elif t in WATER_RESOURCE:
                f = lk.flow(WATER_RESOURCE[t], ("Resources", "Resources from water",
                                                "Renewable material resources from water"))
                ex.append({"type": "biosphere", "amount": v, "input": (BIO, f)})
            elif t == "seed" and not processing:
                ex.append({"type": "technosphere", "amount": v, "input": (FG, code), "product": main["name"],
                           "name": name, "unit": "kilogram", "location": loc})
            elif i["type"] in ("crop", "processedFood") and t in by_product:
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
        if unknown:
            # proxy for unidentified active ingredients: the cycle's identified mix where it carries at
            # least as much mass as the unidentified part, else the global mix
            mix = known if sum(known.values()) >= unknown else global_mix
            total = sum(mix.values())
            for t, kg in mix.items():
                ex += [{"type": "biosphere", "amount": unknown * kg / total * s, "input": (BIO, f),
                        "comment": "PEF pesticide proxy"} for f, s in pesticide_flows(t, lk)]
        for e in farm_emissions(c):
            if not e["value"]:
                continue
            if LUC.search(e["id"]) and e["value"] < 0:        # PEF: no credit for carbon storage
                miss(f"{e['id']} (carbon gain, not credited)", e["value"], e["units"])
                continue
            if e["id"] == "nErosionSoilFlux":
                miss(f"{e['id']} (soil-bound N, left out: nitrate leaching only)", e["value"], e["units"])
                continue
            biofuel = any(BIOFUEL.search(x or "") for x in e.get("inputs") or [])
            code_f = next((f for f in (lk.flow(n, comp) for n, comp in emission_flow(e["id"], biofuel)) if f), None)
            if code_f:
                ex.append({"type": "biosphere", "amount": e["value"], "input": (BIO, code_f)})
            else:
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

    cycles = load_cycles(args.cycles, Path(args.cycles).parent / "cycles_latest_compact_v2.json")
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
