"""PEF multifunctionality handling for Brightway foreground databases, using wurst.

Pipeline
--------
1. extract   : wurst.extract_brightway2_databases([fg_db])  -> list of dicts
2. allocate  : split every multi-output dataset into single-output datasets
               according to the PEF rule named in ``ds["pef_rule"]``
3. cff       : packaging / end-of-life exchanges built with the Circular
               Footprint Formula (linear in datasets -> plain exchanges)
4. write     : wurst.write_brightway2_database(data, f"{fg_db} (PEF)")
5. lcia      : bw2calc with the EF 3.x method + normalisation & weighting

Foreground convention (wurst format): a multi-output dataset has several
``production`` exchanges; each carries ``properties`` such as
``{"mass": kg, "price": EUR/unit, "energy": MJ}`` and the dataset names its
rule in ``pef_rule``.  Rules follow the EF Recommendation (EU) 2021/2279,
Annex I, 4.5 and the agricultural modelling requirements (Annex IV).
"""

import copy

IDF_DAIRY_COEFF = 6.04  # IDF (2015) biophysical allocation, used by the PEF dairy rules


def fpcm(raw_milk_kg, fat_pct, protein_pct):
    """Fat and protein corrected milk (IDF 2015), fat/protein in %."""
    return raw_milk_kg * (0.1226 * fat_pct + 0.0776 * protein_pct + 0.2534)


def _property_shares(products, prop):
    values = [p["amount"] * p["properties"][prop] for p in products]
    total = sum(values)
    if total <= 0:
        raise ValueError(f"Sum of property '{prop}' is zero")
    return [v / total for v in values]


def rule_physical(ds, products, prop="mass"):
    """Allocation on an underlying physical relationship (PEF hierarchy step 2)."""
    return _property_shares(products, ds.get("pef_property", prop))


def rule_economic(ds, products):
    """Economic allocation (hierarchy step 3; slaughterhouse defaults, LEAP crops)."""
    return _property_shares(products, "price")


def rule_dairy_farm(ds, products):
    """Biophysical milk/meat split at farm gate: AF_milk = 1 - 6.04 * M_meat / M_milk.

    The milk output must be expressed in kg FPCM and flagged ``"role": "milk"``;
    all other outputs (live animals sold, kg live weight) share the remainder by mass.
    """
    milk = [p for p in products if p.get("role") == "milk"]
    if len(milk) != 1:
        raise ValueError(f"{ds['name']}: exactly one production exchange needs role 'milk'")
    meat = [p for p in products if p is not milk[0]]
    m_meat = sum(p["amount"] for p in meat)
    af_milk = 1 - IDF_DAIRY_COEFF * m_meat / milk[0]["amount"]
    shares = []
    for p in products:
        if p is milk[0]:
            shares.append(af_milk)
        else:
            shares.append((1 - af_milk) * p["amount"] / m_meat)
    return shares


RULES = {
    "physical": rule_physical,
    "economic": rule_economic,
    "dairy_farm": rule_dairy_farm,
    "slaughterhouse": rule_economic,  # use the PEF default mass fractions and prices
    "crop": rule_economic,            # crop / residue (e.g. straw), LEAP guidelines
}


def allocation_factors(ds):
    products = [e for e in ds["exchanges"] if e["type"] == "production"]
    shares = RULES[ds.get("pef_rule", "physical")](ds, products)
    if abs(sum(shares) - 1) > 1e-9 or min(shares) < 0:
        raise ValueError(f"{ds['name']}: invalid allocation factors {shares}")
    return [(p, s) for p, s in zip(products, shares)]


def split_multioutput(data):
    """Replace multi-output datasets by allocated single-output datasets.

    Returns the new list and a log of the factors applied (for the PEF report).
    Internal consumers are relinked to the split dataset matching their product.
    """
    out, log, relink = [], [], {}
    for ds in data:
        products = [e for e in ds["exchanges"] if e["type"] == "production"]
        if len(products) < 2:
            out.append(ds)
            continue
        for prod, share in allocation_factors(ds):
            new = copy.deepcopy(ds)
            new["code"] = f"{ds['code']}|{prod['product']}"
            new["name"] = f"{ds['name']}, PEF-allocated to {prod['product']}"
            new["reference product"] = prod["product"]
            new["unit"] = prod["unit"]
            new["exchanges"] = [copy.deepcopy(prod)] + [
                dict(e, amount=e["amount"] * share)
                for e in ds["exchanges"]
                if e["type"] != "production"
            ]
            new["exchanges"][0]["name"] = new["name"]
            new["exchanges"][0]["input"] = (ds["database"], new["code"])
            new["allocation"] = {"rule": ds.get("pef_rule", "physical"), "factor": share}
            relink[((ds["database"], ds["code"]), prod["product"])] = (ds["database"], new["code"])
            log.append({"dataset": ds["name"], "product": prod["product"],
                        "rule": new["allocation"]["rule"], "factor": share})
            out.append(new)

    for ds in out:
        for e in ds["exchanges"]:
            if e["type"] == "technosphere" and "input" in e:
                key = (tuple(e["input"]), e.get("product"))
                if key in relink:
                    e["input"] = relink[key]
    return out, log


def cff_material(A, R1, R2, Qsin_Qp, Qsout_Qp, Ev, Erec, ErecEoL, Ev_star=None):
    """Material part of the Circular Footprint Formula as technosphere exchanges.

    (1-R1)Ev + R1(A Erec + (1-A) Ev Qsin/Qp) + (1-A) R2 (ErecEoL - Ev* Qsout/Qp)

    E* arguments are Brightway keys (database, code); returns ``[(key, coeff), ...]``
    per unit of material, ready to multiply by the material mass in the product.
    """
    Ev_star = Ev_star or Ev
    return [
        (Ev, (1 - R1) + R1 * (1 - A) * Qsin_Qp),
        (Erec, R1 * A),
        (ErecEoL, (1 - A) * R2),
        (Ev_star, -(1 - A) * R2 * Qsout_Qp),
    ]


def cff_energy_disposal(B, R2, R3, LHV, X_heat, X_elec, Eer, ESE_heat, ESE_elec, ED):
    """Energy recovery and disposal parts of the CFF, per unit of material.

    (1-B) R3 (Eer - LHV X_heat ESE_heat - LHV X_elec ESE_elec) + (1-R2-R3) ED
    """
    return [
        (Eer, (1 - B) * R3),
        (ESE_heat, -(1 - B) * R3 * LHV * X_heat),
        (ESE_elec, -(1 - B) * R3 * LHV * X_elec),
        (ED, 1 - R2 - R3),
    ]


def add_cff_exchanges(ds, mass, terms, lookup):
    """Append CFF terms (from the two functions above) to a wurst dataset.

    ``lookup`` maps a Brightway key to a dict with name/product/unit/location.
    Terms with the same dataset (e.g. Ev and Ev*) are merged.
    """
    merged = {}
    for key, coeff in terms:
        merged[key] = merged.get(key, 0) + coeff * mass
    for key, amount in merged.items():
        if amount:
            ds["exchanges"].append(dict(lookup[key], type="technosphere",
                                        input=key, amount=amount, comment="CFF"))


def single_score(scores, normalisation, weighting):
    """EF single score: sum_c score_c / NF_c * WF_c (NF/WF from the EF reference package)."""
    return {c: scores[c] / normalisation[c] * weighting[c] for c in scores}
