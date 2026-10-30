"""
Calibrate a REX LCI extension against EXIOBASE air-emission accounts
=====================================================================
rex_env_extension.py scales LCI flow ratios with REX anchors; where the LCI
datasets of a sector are not representative (e.g. Norilsk PGM refining for
"precious metals"), pollutant totals can be far off. EXIOBASE 3 has sector x
region accounts for the main air pollutants, used here as the reference:

    per pollutant P, EXIOBASE region r, sector s:
        lci[r, s]  = sum of all LCI flows of P (all air compartments / species)
                     over the REX countries in r
        exio[r, s] = sum of the EXIOBASE "P - ... - air" stressors
        k = exio / lci
    |log10 k| > log10(threshold)  -> all flows of P in (r, s) scaled by k
                                     (country split and flow split kept)
    lci = 0 and exio > 0 (fill)   -> exio added to the main LCI flow of P,
                                     split over countries by REX GHG
    lci > 0 and exio = 0          -> kept, flagged
    other air flows of the cell (not CO2) -> scaled down by the median factor over the
                                     compared pollutants if below 1/threshold

    python rex_calibrate.py --run <out dir with F_REX.npz> --exio <EXIOBASE satellite F.txt>
        --map-dir <mapping folder> --rex-dir <REX3_2022> [--threshold 3] [--no-fill]
Writes F_REX_calibrated.npz, calibration_cells.csv and calibration_totals.csv into --run.
"""

import argparse
import json
import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from rex_env_extension import N_SECTORS, load_q_rex, rex_sector_codes

# EXIOBASE stressor prefix(es) -> LCI flow names (air, case-insensitive; ecoinvent 3 and EF 3.1)
METAL = r"( ii| iii| iv| vi| ion|\(\d\+\))?$"
POLLUTANTS = {
    "SOx": (["SOx"], r"^(sulfur|sulphur) (dioxide|oxides)$"),
    "NOx": (["NOx", "NOX"], r"^nitrogen oxides$"),
    "NH3": (["NH3"], r"^ammonia$"),
    "PM2.5": (["PM2.5"], r"^(particulates?|particulate matter), < 2\.5 um$|^particles \(pm2\.5\)$"),
    "NMVOC": (["NMVOC"], r"^nmvoc|^non-methane volatile"),
    "CO": (["CO"], r"^carbon monoxide"),
    "CH4": (["CH4"], r"^methane(, fossil|, non-fossil|, biogenic| \(fossil\)| \(biogenic\))?$"),
    "N2O": (["N2O"], r"^dinitrogen monoxide$|^nitrous oxide$"),
    "As": (["As"], r"^arsenic" + METAL),
    "Cd": (["Cd"], r"^cadmium" + METAL),
    "Cr": (["Cr"], r"^chromium" + METAL),
    "Cu": (["Cu"], r"^copper" + METAL),
    "Hg": (["Hg"], r"^mercury" + METAL),
    "Ni": (["Ni"], r"^nickel" + METAL),
    "Pb": (["Pb"], r"^lead" + METAL),
    "Se": (["Se"], r"^selenium" + METAL),
    "Zn": (["Zn"], r"^zinc" + METAL),
}


def load_exiobase(path, exio_name_dict, codes):
    """EXIOBASE air pollutants as {pollutant: (regions, n_sectors) array in REX sector order}."""
    raw = pd.read_csv(path, sep="\t", header=[0, 1], index_col=0)
    regions = list(dict.fromkeys(raw.columns.get_level_values(0)))
    names = list(dict.fromkeys(raw.columns.get_level_values(1)))
    order = [codes.index(c) for c in rex_sector_codes(names, exio_name_dict)]  # exio -> REX index
    out = {}
    for p, (prefixes, _) in POLLUTANTS.items():
        rows = [i for i in raw.index if i.endswith("- air") and i.split(" - ")[0] in prefixes]
        x = raw.loc[rows].sum(axis=0).to_numpy().reshape(len(regions), len(names))
        arr = np.zeros((len(regions), N_SECTORS))
        arr[:, order] = x
        out[p] = arr
    return regions, out


def scale_rows(F, rows, col_scale):
    """F with the given rows multiplied column-wise by col_scale."""
    Fs = F[rows] @ sp.diags(col_scale)
    S = sp.csr_matrix((np.ones(len(rows)), (rows, np.arange(len(rows)))), shape=(F.shape[0], len(rows)))
    return (F + S @ (Fs - F[rows])).tocsr()


def calibrate(F, flows, country_region, regions, exio, ghg, threshold=3.0, fill=True, propagate=True):
    """Return calibrated F (csr) and a DataFrame of changed / flagged cells.

    propagate: per (region, sector) cell, the median log factor exio/lci over all
    compared pollutants is applied to the other air flows of that cell (except CO2,
    which carries the REX GHG anchor) if it is below 1/threshold.
    """
    n_c = len(country_region)
    reg_idx = np.array([regions.index(r) for r in country_region])
    F = F.tocsr()
    cells = []
    logk = []                       # per pollutant: log10(exio/lci) where both > 0, else nan
    covered = np.zeros(F.shape[0], bool)
    for p, (_, pat) in POLLUTANTS.items():
        sel = np.where(flows.compartment.str.startswith("air")
                       & flows.name.str.match(pat, case=False))[0]
        covered[sel] = True
        x = exio[p]
        lci_cs = (np.asarray(F[sel].sum(axis=0)).ravel().reshape(n_c, N_SECTORS) if len(sel)
                  else np.zeros((n_c, N_SECTORS)))
        lci = np.zeros_like(x)
        np.add.at(lci, reg_idx, lci_cs)

        k = np.ones_like(x)
        both = (lci > 0) & (x > 0)
        ratio = np.where(both, x / np.where(both, lci, 1), 1.0)
        logk.append(np.where(both, np.log10(ratio), np.nan))
        out = both & (np.abs(np.log10(ratio)) > np.log10(threshold))
        k[out] = ratio[out]
        add = (lci <= 0) & (x > 0) & fill & (len(sel) > 0)
        flag = (lci > 0) & (x <= 0)

        if out.any() and len(sel):
            F = scale_rows(F, sel, k[reg_idx].ravel())          # per country x sector column
        if add.any():
            main = sel[np.argmax(np.abs(np.asarray(F[sel].sum(axis=1)).ravel()))]
            share = ghg.reshape(n_c, N_SECTORS).clip(min=0)
            reg_tot = np.zeros_like(x)
            np.add.at(reg_tot, reg_idx, share)
            n_in = np.bincount(reg_idx, minlength=len(regions))[:, None]
            w = np.where(reg_tot[reg_idx] > 0, share / np.where(reg_tot[reg_idx] > 0, reg_tot[reg_idx], 1),
                         1.0 / n_in[reg_idx])
            vals = np.where(add[reg_idx], x[reg_idx] * w, 0.0).ravel()
            nz = np.nonzero(vals)[0]
            F = (F + sp.csr_matrix((vals[nz], (np.full(len(nz), main), nz)), shape=F.shape)).tocsr()

        for action, mask in (("scaled", out), ("filled", add), ("flag: EXIOBASE zero", flag)):
            for r, s in zip(*np.nonzero(mask)):
                cells.append({"pollutant": p, "region": regions[r], "sector": s, "lci_kg": lci[r, s],
                              "exiobase_kg": x[r, s], "factor": x[r, s] / lci[r, s] if lci[r, s] > 0 else np.nan,
                              "action": action})

    if propagate:
        with np.errstate(all="ignore"):
            med = np.nanmedian(np.array(logk), axis=0)          # (regions, sectors)
        # downwards only: too low LCI totals mean missing activities, not a scalable profile
        g = np.where(np.isfinite(med) & (med < -np.log10(threshold)), 10 ** np.nan_to_num(med), 1.0)
        other = np.where(flows.compartment.str.startswith("air").values & ~covered
                         & ~flows.name.str.match(r"^carbon dioxide", case=False).values)[0]
        if (g != 1).any() and len(other):
            F = scale_rows(F, other, g[reg_idx].ravel())
        for r, s in zip(*np.nonzero(g != 1)):
            cells.append({"pollutant": "other air flows", "region": regions[r], "sector": s,
                          "lci_kg": np.nan, "exiobase_kg": np.nan, "factor": g[r, s],
                          "action": "propagated"})
    return F, pd.DataFrame(cells)


def pollutant_totals(F, flows, exio):
    tot = np.asarray(F.sum(axis=1)).ravel()
    return {p: tot[(flows.compartment.str.startswith("air")
                    & flows.name.str.match(pat, case=False)).values].sum() for p, (_, pat) in POLLUTANTS.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True, help="output folder of rex_env_extension.py (with F_REX.npz)")
    ap.add_argument("--exio", required=True, help="EXIOBASE satellite F.txt (ixi, same year)")
    ap.add_argument("--map-dir", required=True, help="folder with sector_matching_dicts.pickle")
    ap.add_argument("--rex-dir", required=True, help="folder with Q_REX*.mat (GHG split for filled cells)")
    ap.add_argument("--threshold", type=float, default=3.0, help="scale cells off by more than this factor")
    ap.add_argument("--no-fill", action="store_true", help="do not fill cells without LCI emissions")
    ap.add_argument("--no-propagate", action="store_true",
                    help="do not apply the cell's median factor to the other air flows")
    args = ap.parse_args()

    run = Path(args.run)
    meta = json.load(open(run / "meta.json"))
    flows = pd.read_csv(run / "F_REX_flows.csv")
    F = sp.load_npz(run / "F_REX.npz")
    dicts = pickle.load(open(Path(args.map_dir) / "sector_matching_dicts.pickle", "rb"))
    regions, exio = load_exiobase(args.exio, dicts["exio_name_dict"], meta["exio_codes"])
    missing = set(meta["country_region"]) - set(regions)
    if missing:
        raise KeyError(f"REX regions not in EXIOBASE: {missing}")
    Q = load_q_rex(args.rex_dir, len(meta["q_labels"]), F.shape[1])

    Fc, cells = calibrate(F, flows, meta["country_region"], regions, exio, Q[1],
                          args.threshold, not args.no_fill, not args.no_propagate)
    sp.save_npz(run / "F_REX_calibrated.npz", Fc)
    cells.assign(sector=[meta["sectors"][s] for s in cells.sector] if len(cells) else []).to_csv(
        run / "calibration_cells.csv", index=False)
    totals = pd.DataFrame({"EXIOBASE": {p: x.sum() for p, x in exio.items()},
                           "LCI": pollutant_totals(F, flows, exio),
                           "calibrated": pollutant_totals(Fc, flows, exio)}) / 1e9
    totals.to_csv(run / "calibration_totals.csv", index_label="pollutant [Mt]")
    print("global totals [Mt] (EXIOBASE industries only, no households):")
    print(totals.to_string(float_format="{:.4g}".format))
    print("\ncells:", cells.action.value_counts().to_dict() if len(cells) else {})


if __name__ == "__main__":
    main()
