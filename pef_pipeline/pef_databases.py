"""PEF versions of bafu-2026 and the HESTIA database as Brightway databases, and portable files.

Builds in the ``bafu test`` project (Brightway 2.5):

- ``ef-3.1-water-regional``: country copies of the EF 3.1 water flows the databases use (biosphere)
- ``bafu-2026-pef``: bafu-2026 with the rules of bafu_pef_compare.py written into its exchanges
    fix    unspecified chromium to soil and water as Cr(III), as in BAFU's own EF mapping
    pest   pesticide active ingredients 90 % agricultural soil / 9 % air / 1 % fresh water
           (Recommendation (EU) 2021/2279 Annex I 4.4.1.4; also BAFU's emissions to
           non-agricultural soil in agricultural datasets)
    hm     negative net heavy-metal balances to agricultural soil removed (partial 4.4.1.6 (a))
    water  water flows of datasets in a country with a JRC factor become that country's flow
- ``hestia-2026-04-30-pef``: hestia-2026-04-30 linked to bafu-2026-pef, water flows by country
  (its pesticides already follow the PEF split)
- methods ``("EF v3.1 PEF", <category>)``: the 25 sentier EF v3.1 methods; water use adds the JRC
  country factors of the regional flows. For the original databases they equal sentier EF v3.1.

``--verify`` checks that LCA on bafu-2026-pef reproduces the PEF column of bafu_pef_compare.py and
summarises the change from hestia-2026-04-30 to hestia-2026-04-30-pef. ``--export`` writes the
databases and methods in Brightway's write format to reports/pef_bw/ with a manifest;
import_pef_databases.py loads them into any project (Brightway 2.5 or legacy 2). The build writes
the change log and CHANGES.md there too.

The files contain modified BAFU data: keep them local (BAFU terms of use; CHANGES.md and the change
log are the documentation of changes the terms require).

    python pef_databases.py [--build] [--verify] [--export] [--jrc <EF-LCIAMethod_CF(EF-v3.1).xlsx>]
"""

import argparse
import csv
import datetime
import gzip
import json
from collections import defaultdict
from pathlib import Path

import bw2calc as bc
import bw2data as bd
import numpy as np

import bafu_pef_compare as cmp
from export_hestia_db import portable_db
from pef_allocation import scaled

PROJECT, BIO, BAFU = cmp.PROJECT, cmp.BIO, cmp.DB
RESIDUAL, HESTIA = "bafu-2026-residual", "hestia-2026-04-30"
WATER_DB, BAFU_PEF, HESTIA_PEF = "ef-3.1-water-regional", "bafu-2026-pef", "hestia-2026-04-30-pef"
SOURCE_METHODS, PEF_METHODS = ("sentier", "EF v3.1"), "EF v3.1 PEF"
RELINK = {BAFU: BAFU_PEF, HESTIA: HESTIA_PEF}
ORDER = [BIO, WATER_DB, RESIDUAL, BAFU_PEF, HESTIA_PEF]  # import order
OUT = Path(__file__).parent.parent / "reports" / "pef_bw"
CHANGES = ["database", "dataset", "location", "code", "rule", "flow", "compartment", "new flow",
           "new compartment", "amount", "new amount"]
RULES = {
    "fix": "unspecified chromium to soil and water mapped to Cr(III), as in BAFU's own EF mapping "
           "(correction of the Brightway import, not a PEF rule)",
    "pest": "pesticide active ingredient split 90 % agricultural soil / 9 % air / 1 % fresh water "
            "(PEF Annex I 4.4.1.4); emissions to non-agricultural soil of agricultural datasets included",
    "hm": "negative net heavy-metal balance to agricultural soil set to zero (partial PEF Annex I 4.4.1.6 (a))",
    "water": "water flow replaced by the flow of the dataset's country, characterised with the JRC "
             "EF 3.1 country factor (PEF: water use at country level)",
}


class Builder:
    """Applies the rules exchange by exchange and records every change."""

    def __init__(self, jrc_cf):
        self.flows, self.by_name, self.cr3, self.water = cmp.flow_lookups(cmp.water_factors(jrc_cf))
        self.regional = {}  # regional flow key -> (EF 3.1 flow key, location)
        self.log = []

    def describe(self, key):
        if key is None:
            return "", ""
        if key[0] == WATER_DB:
            base, loc = self.regional[key]
            f = self.flows[base]
            return f"{f['name']}, {loc}", f["categories"][-1]
        f = self.flows[key]
        return f["name"], (f.get("categories") or ("",))[-1]

    def record(self, db, ds, rule, key, new_key, amount, new_amount):
        self.log.append(dict(zip(CHANGES, [db, ds["name"], ds.get("location"), ds["code"], rule,
                                           *self.describe(key), *self.describe(new_key), amount, new_amount])))

    def regionalise(self, db, ds, e):
        """A water flow of a dataset in a country with a JRC factor becomes that country's flow."""
        key, loc = tuple(e["input"]), ds.get("location")
        if key not in self.water or self.water[key].get(loc) is None:
            return e
        new = (WATER_DB, f"{key[1]}_{loc}")
        self.regional[new] = (key, loc)
        self.record(db, ds, "water", key, new, e["amount"], e["amount"])
        return dict(e, input=new)

    def pesticide_targets(self, f, nonag):
        """[(flow, share)] of the 90/9/1 split, or None where the EF biosphere lacks a target flow."""
        targets = self.by_name[f["name"].lower()]
        soil = (cmp.pick(targets[cmp.AG_SOIL_CAT], f) if cmp.AG_SOIL_CAT in targets else None) if nonag else f
        if soil is None or cmp.AIR not in targets or cmp.WATER not in targets:
            return None
        return [(soil, cmp.SPLIT["soil"]), (cmp.pick(targets[cmp.AIR], f), cmp.SPLIT["air"]),
                (cmp.pick(targets[cmp.WATER], f), cmp.SPLIT["water"])]

    def bafu_exchanges(self, ds):
        agricultural = ds.get("sector") == "agriculture"
        out, metals = [], defaultdict(list)
        for e in ds["exchanges"]:
            key = tuple(e["input"])
            if e["type"] != "biosphere" or key[0] != BIO:
                out.append(e)
                continue
            f = self.flows[key]
            cat, name = tuple(f.get("categories") or ()), f["name"]
            if name.lower() == "chromium" and len(cat) > 1 and cat[1] in cmp.CR_REMAP and cat in self.cr3:
                new = self.cr3[cat].key
                self.record(BAFU_PEF, ds, "fix", key, new, e["amount"], e["amount"])
                e, key = dict(e, input=new), new
            if key in self.water:
                out.append(self.regionalise(BAFU_PEF, ds, e))
                continue
            comp = cat[-1] if cat else ""
            if comp == cmp.AG_SOIL and cmp.HEAVY_METAL.match(name):
                metals[key].append(e)
                continue
            if (comp == cmp.AG_SOIL or (comp == cmp.NONAG_SOIL and agricultural)) \
                    and not cmp.HEAVY_METAL.match(name) and cmp.is_pesticide(name, agricultural):
                split = self.pesticide_targets(f, comp == cmp.NONAG_SOIL)
                if split:
                    for target, share in split:
                        out.append(dict(scaled(e, share), input=target.key))
                        self.record(BAFU_PEF, ds, "pest", key, target.key, e["amount"], e["amount"] * share)
                    continue
            out.append(e)
        for key, group in metals.items():  # per dataset and flow: a negative net balance becomes zero
            net = sum(x["amount"] for x in group)
            if net < 0:
                self.record(BAFU_PEF, ds, "hm", key, None, net, 0.0)
            else:
                out.extend(group)
        return out

    def hestia_exchanges(self, ds):
        return [self.regionalise(HESTIA_PEF, ds, e) if e["type"] == "biosphere" else e for e in ds["exchanges"]]

    def water_database(self):
        data = {}
        for new, (key, loc) in sorted(self.regional.items()):
            f = self.flows[key]
            data[new] = {"database": WATER_DB, "code": new[1], "name": f"{f['name']}, {loc}",
                         "categories": tuple(f["categories"]), "unit": f["unit"], "type": f.get("type", "emission"),
                         "location": loc, "CAS number": f.get("CAS number"),
                         "comment": f"EF 3.1 flow {key[1]} for datasets in {loc}; water use is characterised "
                                    f"with the JRC EF 3.1 country factor", "exchanges": []}
        return data


def pef_copy(source, target, exchanges):
    """{(target, code): dataset}: exchanges(ds) gives the new exchanges; BAFU and HESTIA links are renamed."""
    data = {}
    for (_, code), ds in bd.Database(source).load().items():
        ds = {k: v for k, v in ds.items() if k not in ("id", "output")}
        new = []
        for e in exchanges(ds):
            db, c = e["input"]
            new.append(dict({k: v for k, v in e.items() if k != "output"}, input=(RELINK.get(db, db), c)))
        data[(target, code)] = dict(ds, database=target, exchanges=new)
    return data


def write(name, data):
    if name in bd.databases:
        del bd.databases[name]
    bd.Database(name).write(data)
    print(f"{name}: {len(data)} datasets")


def method_keys():
    keys = {}
    for db in (BIO, WATER_DB):
        keys.update({a.id: a.key for a in bd.Database(db)})
    return keys


def write_methods(builder):
    keys = method_keys()
    for m in [m for m in bd.methods if m[:2] == SOURCE_METHODS]:
        cfs = [(keys[r[0]] if isinstance(r[0], int) else tuple(r[0]), r[1]) for r in bd.Method(m).load()]
        meta = {k: v for k, v in bd.methods[m].items() if k == "unit"}
        meta["description"] = "EF v3.1 (sentier)"
        if m[2] == "Water use":
            cfs += [(new, builder.water[key][loc]) for new, (key, loc) in sorted(builder.regional.items())]
            meta["description"] += f"; plus the JRC EF 3.1 country factors for the flows of {WATER_DB}"
        method = bd.Method((PEF_METHODS, m[2]))
        if method.name in bd.methods:
            method.deregister()
        method.register(**meta)
        method.write(cfs)
    print(f"{len([m for m in bd.methods if m[0] == PEF_METHODS])} methods {PEF_METHODS}")


def write_log(builder, out):
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "pef_changes.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=CHANGES)
        w.writeheader()
        w.writerows(builder.log)
    print("written", out / "pef_changes.csv")


def write_documentation(out):
    """CHANGES.md from the change log; a pesticide exchange has one log row per target compartment,
    counted by its agricultural-soil row."""
    with open(out / "pef_changes.csv", encoding="utf-8") as f:
        log = list(csv.DictReader(f))
    exchanges, datasets = defaultdict(int), defaultdict(set)
    for r in log:
        group = (r["database"], r["rule"])
        exchanges[group] += r["rule"] != "pest" or r["new compartment"] == cmp.AG_SOIL
        datasets[group].add(r["code"])
    n_regional = len(bd.Database(WATER_DB))
    lines = [f"# PEF versions of BAFU:2026 and HESTIA ({datetime.date.today()})", "",
             "Modified data. `bafu-2026-pef` is derived from the Life Cycle Inventory database of the Swiss "
             "Federal Administration, BAFU:2026; `hestia-2026-04-30-pef` from the HESTIA aggregated release "
             "2026-04-30 (CC BY 4.0). Neither BAFU nor HESTIA has reviewed or endorsed these changes. Cite the "
             "original sources. BAFU's terms of use do not allow distributing the BAFU data separately: keep "
             "these files local and pass them on only under terms equivalent to BAFU's.", "",
             "## Changes", "", "Every changed exchange is listed in `pef_changes.csv` (dataset, flow, old and new "
             "amount and flow).", "", "| Database | Rule | What changed | Datasets | Exchanges |", "|---|---|---|---|---|"]
    for (db, rule), n in sorted(exchanges.items()):
        lines.append(f"| `{db}` | {rule} | {RULES[rule]} | {len(datasets[(db, rule)])} | {n} |")
    lines += ["", f"`{HESTIA_PEF}` also links to `{BAFU_PEF}` instead of `{BAFU}`. Its pesticides already "
              "follow the PEF split (built by hestia_db.py).", "",
              f"`{WATER_DB}` holds {n_regional} country copies of EF 3.1 water flows. The methods "
              f"`(\"{PEF_METHODS}\", <category>)` are the 25 sentier EF v3.1 methods; water use adds the JRC "
              "country factors of these flows. Use them for all databases: for the original ones they equal "
              "sentier EF v3.1.", "",
              "PEF category rules can override these defaults (the dairy PEFCR keeps pesticides 100 % in soil).", "",
              "## Import", "", "`python import_pef_databases.py <project> --files <this folder>` (Brightway 2.5 "
              "or legacy Brightway 2)."]
    (out / "CHANGES.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("written", out / "CHANGES.md")


def build(jrc_cf, out):
    builder = Builder(jrc_cf)
    bafu = pef_copy(BAFU, BAFU_PEF, builder.bafu_exchanges)
    hestia = pef_copy(HESTIA, HESTIA_PEF, builder.hestia_exchanges)
    write(WATER_DB, builder.water_database())
    write(BAFU_PEF, bafu)
    write(HESTIA_PEF, hestia)
    write_methods(builder)
    write_log(builder, out)
    write_documentation(out)


def scores(acts, methods):
    """{node id: {category: score}} for all nodes with one LCA (its matrices span the systems of the
    first and the last node)."""
    lca = bc.LCA({acts[0].id: 1, acts[-1].id: 1}, methods[0])
    lca.lci()
    cvec = {}
    for m in methods:
        lca.switch_method(m)
        cvec[m[-1]] = lca.characterization_matrix.diagonal()
    out = {}
    for a in acts:
        lca.redo_lci({a.id: 1})
        g = lca.biosphere_matrix @ lca.supply_array
        out[a.id] = {c: float(v @ g) for c, v in cvec.items()}
    return out


def verify():
    pef = [m for m in bd.methods if m[0] == PEF_METHODS]
    with open(cmp.OUT, encoding="utf-8") as f:
        ref = {(r["product"], r["location"], r["category"]): float(r["PEF"]) for r in csv.DictReader(f)}
    cats = sorted({c for _, _, c in ref})
    products = cmp.food_products()
    acts = [bd.get_node(database=BAFU_PEF, code=a["code"]) for a in products]
    got = scores(acts, [m for m in pef if m[1] in cats])
    print(f"\n{BAFU_PEF} vs. the PEF column of bafu_pef_compare.py, {len(acts)} products (largest relative difference):")
    for c in cats:
        dev = [abs(got[a.id][c] - ref[(a["name"], a["location"], c)]) / max(abs(ref[(a["name"], a["location"], c)]), 1e-30)
               for a in acts if ref[(a["name"], a["location"], c)]]
        print(f"  {c:<45} {max(dev):.1e}")

    h_old = list(bd.Database(HESTIA))
    h_new = [bd.get_node(database=HESTIA_PEF, code=a["code"]) for a in h_old]
    both = scores(h_old + h_new, [m for m in pef if m[1] in cats])
    print(f"\n{HESTIA_PEF} / {HESTIA}, {len(h_old)} datasets (median, 5th-95th percentile, changed > 1 %):")
    for c in cats:
        r = np.array([both[n.id][c] / both[o.id][c] for o, n in zip(h_old, h_new) if both[o.id][c] > 0])
        print(f"  {c:<45} {np.median(r):6.3f}  {np.percentile(r, 5):6.3f}-{np.percentile(r, 95):<7.3f}"
              f" {int((abs(r - 1) > 0.01).sum()):4d} of {len(r)}")


def portable_methods():
    keys = method_keys()
    return [{"name": list(m), "metadata": {k: v for k, v in bd.methods[m].items() if k in ("unit", "description")},
             "cfs": [[list(keys[r[0]] if isinstance(r[0], int) else r[0]), r[1]] for r in bd.Method(m).load()]}
            for m in bd.methods if m[0] == PEF_METHODS]


def export(out):
    out.mkdir(parents=True, exist_ok=True)
    for db in ORDER:
        with gzip.open(out / f"{db}.json.gz", "wt", encoding="utf-8") as f:
            json.dump(portable_db(db), f)
        print("written", out / f"{db}.json.gz")
    with gzip.open(out / "ef-3.1-pef-methods.json.gz", "wt", encoding="utf-8") as f:
        json.dump(portable_methods(), f)
    manifest = {"databases": ORDER, "methods": "ef-3.1-pef-methods.json.gz", "created": str(datetime.date.today()),
                "project": PROJECT, "documentation": "CHANGES.md"}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8")
    print("written", out / "manifest.json")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--export", action="store_true")
    ap.add_argument("--docs", action="store_true", help="rewrite CHANGES.md from the change log")
    ap.add_argument("--jrc", default=str(cmp.JRC_CF))
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()
    bd.projects.set_current(PROJECT)
    if args.build:
        build(args.jrc, Path(args.out))
    if args.verify:
        verify()
    if args.docs:
        write_documentation(Path(args.out))
    if args.export:
        export(Path(args.out))


if __name__ == "__main__":
    main()
