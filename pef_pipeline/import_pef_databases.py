"""Write the portable PEF files (from pef_databases.py --export) into a Brightway project.

Reads ``manifest.json`` in the folder: the databases in import order and the methods file. Works
with Brightway 2.5 (Activity Browser 3) and legacy Brightway 2 (Activity Browser 2): it only uses
``Database.write`` and ``Method.write``. Databases and methods of the same name are replaced.

    python import_pef_databases.py <project> [--files <folder>]
"""

import argparse
import gzip
import json
from pathlib import Path

import bw2data as bd

FILES = Path(__file__).parent.parent / "reports" / "pef_bw"


def load(path):
    with gzip.open(path, "rt", encoding="utf-8") as f:
        return json.load(f)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project")
    ap.add_argument("--files", default=str(FILES))
    args = ap.parse_args()
    folder = Path(args.files)
    manifest = json.loads((folder / "manifest.json").read_text(encoding="utf-8"))

    bd.projects.set_current(args.project)
    for name in manifest["databases"]:
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
    methods = load(folder / manifest["methods"])
    for m in methods:
        method = bd.Method(tuple(m["name"]))
        if method.name in bd.methods:
            method.deregister()
        method.register(**m["metadata"])
        method.write([(tuple(k), cf) for k, cf in m["cfs"]])
    print(f"{len(methods)} methods")


if __name__ == "__main__":
    main()
