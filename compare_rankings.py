"""
Ranking comparison of footprint approaches across EF 3.1 categories
====================================================================
For each EF category and each pair of approaches (footprints from rex_footprint.py),
how similarly are producer countries and producer sectors ranked?

    spearman   rank correlation over all countries (or sectors)
    top10      how many of approach A's top 10 are in approach B's top 10
    top1       whether both name the same largest contributor
    ratio      total B / total A

    python compare_rankings.py --dir <footprint folder> --country Switzerland
        --approaches ecoinvent_raw ecoinvent bafu_raw bafu --out rankings.xlsx
"""

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def load(d, country, label, by):
    return pd.read_csv(Path(d) / f"footprint_{country}_{label}_by_{by}.csv", index_col=0)


def compare(a, b):
    """Rank agreement of two contribution vectors (same index)."""
    rho = spearmanr(a, b).statistic if a.std() > 0 and b.std() > 0 else np.nan
    ta, tb = set(a.nlargest(10).index), set(b.nlargest(10).index)
    return {"spearman": rho, "top10": len(ta & tb), "top1": a.idxmax() == b.idxmax(),
            "ratio": b.sum() / a.sum() if a.sum() else np.nan}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--country", default="Switzerland")
    ap.add_argument("--approaches", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows, tops = [], {}
    for by in ("producer", "sector"):
        tabs = {l: load(args.dir, args.country, l, by) for l in args.approaches}
        cats = list(tabs[args.approaches[0]].columns)
        for cat in cats:
            for a, b in combinations(args.approaches, 2):
                rows.append({"by": by, "category": cat, "A": a, "B": b,
                             **compare(tabs[a][cat], tabs[b][cat])})
            for l in args.approaches:
                top = tabs[l][cat].nlargest(3)
                tops[(by, cat, l)] = "; ".join(f"{i[:32]} {v / tabs[l][cat].sum():.0%}" for i, v in top.items())
    res = pd.DataFrame(rows)
    top = pd.Series(tops).unstack()
    top.index.names = ["by", "category"]

    summary = (res.groupby(["by", "A", "B"])
               .agg(spearman_median=("spearman", "median"), top10_median=("top10", "median"),
                    top1_share=("top1", "mean"), ratio_median=("ratio", "median"))
               .reset_index())
    with pd.ExcelWriter(args.out) as w:
        summary.to_excel(w, sheet_name="summary", index=False)
        for by in ("producer", "sector"):
            r = res[res.by == by]
            r.pivot_table(index="category", columns=["A", "B"], values="spearman").to_excel(w, sheet_name=f"spearman {by}")
            r.pivot_table(index="category", columns=["A", "B"], values="top10").to_excel(w, sheet_name=f"top10 {by}")
        top.to_excel(w, sheet_name="top 3 per approach")
        res.to_excel(w, sheet_name="all", index=False)
    pd.set_option("display.width", 250)
    print(summary.to_string(index=False, float_format="{:.2f}".format))
    return res, top


if __name__ == "__main__":
    main()
