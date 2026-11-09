"""Cross-check our Brightway EF 3.1 baseline against BAFU's published LCIA results.

Reads the BAFU column of reports/bafu_pef_compare.csv (written by bafu_pef_compare.py) and
the release workbook "BAFU-2026 v1 LCIA Results_corrected.xlsx" (second sheet, header in
row 2, products keyed "<name> - <location>"), and reports the ratio per product and category.

    python bafu_workbook_check.py "<path to the workbook>"
"""

import sys
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

HERE = Path(__file__).parent
OURS = HERE.parent / "reports" / "bafu_pef_compare.csv"
OUT = HERE.parent / "reports" / "bafu_workbook_check.csv"  # licensed data: stays local
# our category name -> workbook column header (without the unit)
COLUMNS = {
    "Acidification": "Acidification",
    "Climate change": "Climate change",
    "Ecotoxicity, freshwater": "Ecotoxicity freshwater",
    "Eutrophication, freshwater": "Eutrophication freshwater",
    "Eutrophication marine": "Eutrophication marine",
    "Eutrophication, terrestrial": "Eutrophication terrestrial",
    "Human toxicity, cancer": "Human toxicity cancer",
    "Human toxicity, non-cancer": "Human toxicity non-cancer",
    "Ionising radiation, human health": "Ionising radiation (human health)",
    "Land use": "Land use",
    "Ozone depletion": "Ozone depletion",
    "EF-particulate Matter": "Particulate matter",
    "Photochemical ozone formation - human health": "Photochemical ozone formation (human health)",
    "Resource use, fossils": "Resource use fossils",
    "Resource use, minerals and metals": "Resource use minerals and metals",
    "Water use": "Water use",
}


def read_workbook(path):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    rows = wb[wb.sheetnames[1]].iter_rows(values_only=True)
    banner, header = next(rows), next(rows)
    start = banner.index("EF 3.1")  # the Ecological Scarcity "Land use [UBP]" comes earlier
    index = {}
    for i, h in enumerate(header):
        if h and i >= start:
            index.setdefault(h.split("[")[0].strip(), i)
    cols = {c: index[w] for c, w in COLUMNS.items()}
    data = {r[0]: {c: r[i] for c, i in cols.items()} for r in rows if r[0]}
    wb.close()
    return data


def main(path):
    published = read_workbook(path)
    ours = pd.read_csv(OURS, keep_default_na=False)  # "NA" is Namibia, not a missing value
    ours["key"] = ours["product"] + " - " + ours["location"]
    missing = sorted(set(ours.key) - set(published))
    print(f"{ours.key.nunique()} products, {len(missing)} not in the workbook: {missing[:5]}")

    ours = ours[ours.key.isin(published)].copy()
    ours["published"] = [published[k][c] for k, c in zip(ours.key, ours.category)]
    ours = ours[abs(ours.published) > 1e-9 * ours.groupby("category").published.transform(lambda p: abs(p).max())]
    ours["ratio"] = ours.BAFU / ours.published
    ours["ratio_raw"] = ours.BAFU_raw / ours.published
    ours[["product", "location", "category", "unit", "BAFU_raw", "BAFU", "published", "ratio_raw", "ratio"]] \
        .to_csv(OUT, index=False)

    print(f"\nours / published per category: products within 1 % / 10 %, median; "
          f"without -> with the import correction -> {OUT}")
    for c in COLUMNS:
        x = ours[ours.category == c]
        stats = [(int((abs(r - 1) <= .01).sum()), int((abs(r - 1) <= .1).sum()), r.median())
                 for r in (x.ratio_raw, x.ratio)]
        print(f"  {c:<45} " + "  ->  ".join(f"{a:3d} / {b:3d}  {m:6.3f}" for a, b, m in stats) + f"   of {len(x)}")
    worst = ours.assign(dev=lambda d: np.where(d.ratio < 0, np.inf, abs(np.log(abs(d.ratio)))))
    print("\nLargest deviations (sign changes first):")
    print(worst.sort_values("dev", ascending=False).head(15)[["key", "category", "BAFU", "published", "ratio"]]
          .to_string(index=False, float_format=lambda v: f"{v:.4g}"))


if __name__ == "__main__":
    main(sys.argv[1])
