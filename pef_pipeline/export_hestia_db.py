"""Make the HESTIA database usable outside ``bafu test``: own project, portable files, Excel.

1. Project ``hestia-bafu-pef`` (Brightway 2.5, Activity Browser 3): a copy of ``bafu test`` reduced to
   ``hestia-2026-04-30``, ``bafu-2026``, ``bafu-2026-residual``, ``ef-3.1-biosphere`` and the 25 EF 3.1
   methods.
2. ``<out>/<database>.json.gz`` and ``<out>/ef-3.1-methods.json.gz``: the same in Brightway's write
   format -- datasets and exchanges keyed by (database, code), characterisation factors keyed by flow
   (database, code). ``import_hestia_db.py`` writes them into any project, in Brightway 2.5 or in
   legacy Brightway 2 (Activity Browser 2). BW2Packages are not used: the 2.5 ones store factors by
   internal ids and do not open in legacy Brightway.
3. ``<out>/lci-hestia-2026-04-30.xlsx``: the HESTIA database as a Brightway Excel sheet, for reading.

    python export_hestia_db.py [--out <folder>]

The files contain BAFU data: keep them local (BAFU terms of use).
"""

import argparse
import gzip
import json
from pathlib import Path

import bw2data as bd
import bw2io as bi

SOURCE, PROJECT = "bafu test", "hestia-bafu-pef"
KEEP = ["ef-3.1-biosphere", "bafu-2026-residual", "bafu-2026", "hestia-2026-04-30"]  # import order
METHODS = ("sentier", "EF v3.1")
OUT = Path(__file__).parent.parent / "reports" / "hestia_bw"


def dump(obj, path):
    with gzip.open(path, "wt", encoding="utf-8") as f:
        json.dump(obj, f)
    print("written", path)


def portable_db(name):
    """{'database|code': dataset}; legacy Brightway knows only the type 'process'."""
    out = {}
    for (db, code), ds in bd.Database(name).load().items():
        ds = {k: v for k, v in ds.items() if k not in ("id", "output")}
        if ds.get("type") == "processwithreferenceproduct":
            ds["type"] = "process"
        ds["exchanges"] = [{k: v for k, v in e.items() if k != "output"} for e in ds.get("exchanges", [])]
        out[f"{db}|{code}"] = ds
    return out


def portable_methods():
    keys = {a.id: a.key for db in KEEP for a in bd.Database(db)}
    return [{"name": list(m), "metadata": {k: v for k, v in bd.methods[m].items()
                                           if k in ("unit", "description", "abbreviation")},
             "cfs": [[list(keys[i] if isinstance(i, int) else i), cf] for i, cf in bd.Method(m).load()]}
            for m in bd.methods if m[:2] == METHODS]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    bd.projects.set_current(SOURCE)
    if PROJECT in bd.projects:
        bd.projects.delete_project(PROJECT, delete_dir=True)
    bd.projects.copy_project(PROJECT, switch=True)
    for db in list(bd.databases):
        if db not in KEEP:
            del bd.databases[db]
    for m in list(bd.methods):
        if m[:2] != METHODS:
            bd.Method(m).deregister()
    print(f"project {PROJECT}: {sorted(bd.databases)}, {len(bd.methods)} methods")

    for db in KEEP:
        dump(portable_db(db), out / f"{db}.json.gz")
    dump(portable_methods(), out / "ef-3.1-methods.json.gz")
    print("written", bi.export.write_lci_excel(KEEP[-1], dirpath=str(out)))


if __name__ == "__main__":
    main()
