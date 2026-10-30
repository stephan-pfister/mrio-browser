"""
Compare REX extension variants built by rex_env_extension.py
=============================================================
Per sector and key pollutant: emission per kg CO2-eq (global pool ratio against
the GHG anchor), summed over subcompartments. Flow names differ between
ecoinvent biosphere3 and EF 3.1, so pollutants are matched by name patterns.

    python compare_extensions.py <out_dir_1> <out_dir_2> ... --labels a b ... --csv result.csv
If the runs contain F_REX.npz (built with --rex-dir), global totals are compared too.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

# pollutant -> (compartment prefix, name regex, case-insensitive; ecoinvent 3 and EF 3.1 names)
POLLUTANTS = {
    "SO2 [air]": ("air", r"^sulfur dioxide$"),
    "NOx [air]": ("air", r"^nitrogen oxides$"),
    "NH3 [air]": ("air", r"^ammonia$"),
    "PM2.5 [air]": ("air", r"^(particulates?|particulate matter), < 2\.5 um$|^particles \(pm2\.5\)$"),
    "NMVOC [air]": ("air", r"^nmvoc|^non-methane volatile"),
    "CO [air]": ("air", r"^carbon monoxide(, fossil| \(fossil\))$"),
    "CH4 fossil [air]": ("air", r"^methane(, fossil| \(fossil\))$"),
    "N2O [air]": ("air", r"^dinitrogen monoxide$|^nitrous oxide$"),
    "Hg [air]": ("air", r"^mercury( ii|\(2\+\))?$"),
    "Cd [air]": ("air", r"^cadmium( ii|\(2\+\))?$"),
    "Pb [air]": ("air", r"^lead( ii|\(2\+\))?$"),
    "Phosphate [water]": ("water", r"^phosphate$"),
    "Nitrate [water]": ("water", r"^nitrate$"),
    "COD [water]": ("water", r"^cod|^chemical oxygen demand"),
}


def select(flows, comp, pat):
    return (flows.compartment.str.startswith(comp) & flows.name.str.match(pat, case=False)).values


def pollutant_ratios(out):
    """DataFrame sector x pollutant: kg per kg CO2-eq (global pool, GHG anchor)."""
    meta = json.load(open(out / "meta.json"))
    z = np.load(out / "ratios.npz")
    flows = pd.read_csv(out / "F_REX_flows.csv")
    pools = np.load(sorted(out.glob("pools_*.npz"))[0])
    E, A = pools["E"], pools["A"]
    keep = np.abs(E[:, :, -1]).sum(axis=1) > 0
    E = E[keep][:, :, -1]
    ghg = A[0, :, -1]
    res = {}
    for label, (comp, pat) in POLLUTANTS.items():
        tot = E[select(flows, comp, pat)].sum(axis=0)
        res[label] = np.where(ghg > 0, tot / np.where(ghg > 0, ghg, 1), np.nan)
    df = pd.DataFrame(res, index=meta["sectors"])
    df.insert(0, "n_acts", z["n_acts"][:, -1])
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--labels", nargs="+")
    ap.add_argument("--csv")
    args = ap.parse_args()
    labels = args.labels or [Path(r).name for r in args.runs]

    tables = {l: pollutant_ratios(Path(r)) for l, r in zip(labels, args.runs)}
    long = pd.concat({l: t.drop(columns="n_acts").stack(future_stack=True) for l, t in tables.items()}, axis=1)
    long.index.names = ["sector", "pollutant"]
    acts = pd.concat({l: t.n_acts for l, t in tables.items()}, axis=1)
    if args.csv:
        long.to_csv(args.csv)

    print("sectors with activities:", {l: int((a > 0).sum()) for l, a in acts.items()})
    pairs = [(labels[i], labels[j]) for i in range(len(labels)) for j in range(i + 1, len(labels))]
    for base, other in pairs:
        both = long[[base, other]].replace(0, np.nan).dropna()
        both = both[(both > 0).all(axis=1)]
        lr = np.log10(both[other] / both[base])
        print(f"\n{other} vs {base}: {len(both)} sector-pollutant pairs, "
              f"{int((acts[base] > 0).mul(acts[other] > 0).sum())} common sectors")
        print(f"  median ratio {10 ** lr.median():.2f}, within x2: {(lr.abs() < np.log10(2)).mean():.0%}, "
              f"within x10: {(lr.abs() < 1).mean():.0%}")
        print("  median ratio by pollutant:")
        print((10 ** lr.groupby(level="pollutant").median()).round(2).to_string())

    # global totals when F_REX.npz exists in all runs
    if all((Path(r) / "F_REX.npz").exists() for r in args.runs):
        import scipy.sparse as sp
        tot = {}
        for l, r in zip(labels, args.runs):
            F = sp.load_npz(Path(r) / "F_REX.npz")
            flows = pd.read_csv(Path(r) / "F_REX_flows.csv")
            s = np.asarray(F.sum(axis=1)).ravel()
            tot[l] = {p: s[select(flows, c, pat)].sum() for p, (c, pat) in POLLUTANTS.items()}
        print("\nglobal totals (REX-scaled):")
        print(pd.DataFrame(tot).to_string(float_format="{:.3g}".format))
    return long


if __name__ == "__main__":
    main()
