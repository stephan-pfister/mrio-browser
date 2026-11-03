"""Write the portable HESTIA files (from export_hestia_db.py) into a Brightway project.

Works with Brightway 2.5 (Activity Browser 3) and legacy Brightway 2 (Activity Browser 2): it only
uses ``Database.write`` and ``Method.write``. Existing databases and methods of the same name in the
target project are replaced.

    python import_hestia_db.py <project> [--files <folder>]
"""

import argparse
import gzip
import json
from pathlib import Path

import bw2data as bd

ORDER = ["ef-3.1-biosphere", "bafu-2026-residual", "bafu-2026", "hestia-2026-04-30"]
FILES = Path(__file__).parent.parent / "reports" / "hestia_bw"


def load(path):
    with gzip.open(path, "rt", encoding="utf-8") as f:
        return json.load(f)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project")
    ap.add_argument("--files", default=str(FILES))
    args = ap.parse_args()
    folder = Path(args.files)

    bd.projects.set_current(args.project)
    for name in ORDER:
        data = {}
        for key, ds in load(folder / f"{name}.json.gz").items():
            if ds.get("categories"):
                ds["categories"] = tuple(ds["categories"])
            else:
                ds.pop("categories", None)  # legacy search index cannot handle None
            for e in ds["exchanges"]:
                e["input"] = tuple(e["input"])
            data[tuple(key.split("|", 1))] = ds
        if name in bd.databases:
            del bd.databases[name]
        bd.Database(name).register()
        bd.Database(name).write(data)
        print(f"{name}: {len(data)} datasets")
    methods = load(folder / "ef-3.1-methods.json.gz")
    for m in methods:
        method = bd.Method(tuple(m["name"]))
        if method.name in bd.methods:
            method.deregister()
        method.register(**m["metadata"])
        method.write([(tuple(k), cf) for k, cf in m["cfs"]])
    print(f"{len(methods)} EF 3.1 methods")


if __name__ == "__main__":
    main()
