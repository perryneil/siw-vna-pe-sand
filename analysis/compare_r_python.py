#!/usr/bin/env python3
"""
compare_r_python.py
===================

Cross-checks siw_analysis.R against siw_analysis.py.

1. Runs both scripts and diffs the printed reports line by line.
2. Recomputes every underlying quantity at full double precision with the
   Python functions, reads the same quantities dumped by the R script
   (--dump), and reports the largest absolute and relative difference for
   each.

Usage (from the repository root):
    python analysis/compare_r_python.py
    python analysis/compare_r_python.py --rscript "C:/Program Files/R/R-4.4.1/bin/Rscript.exe"

Requires: numpy, scipy, and R (Rscript on PATH or given with --rscript).
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import siw_analysis as S  # noqa: E402


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def num(rows, key):
    return np.array([float(r[key]) if r[key] not in ("NA", "") else np.nan for r in rows])


def python_reference(data_dir):
    D = S.load_all(data_dir)
    freq, air, pe, sand = D["freq"], D["air"], D["pe"], D["sand"]
    out = {}

    raw_p = np.array([stats.ttest_ind(sand[:, k], pe[:, k], equal_var=False).pvalue
                      for k in range(len(freq))])
    raw_p = np.nan_to_num(raw_p, nan=1.0)
    holm, bh = S.holm_bonferroni(raw_p), S.benjamini_hochberg(raw_p)
    d_all = np.array([S.cohens_d(sand[:, k], pe[:, k])[0] for k in range(len(freq))])
    out["pointwise"] = dict(raw_p=raw_p, holm=holm, bh=bh,
                            diff=sand.mean(0) - pe.mean(0), d=d_all)

    t2 = {k: [] for k in ("air", "sand_mean", "sand_sd", "pe_mean", "pe_sd", "diff", "d")}
    for f0 in S.FEATURES.values():
        k = int(np.argmin(np.abs(freq - f0)))
        s, p_ = sand[:, k], pe[:, k]
        t2["air"].append(air[k])
        t2["sand_mean"].append(s.mean()); t2["sand_sd"].append(s.std(ddof=1))
        t2["pe_mean"].append(p_.mean()); t2["pe_sd"].append(p_.std(ddof=1))
        t2["diff"].append(s.mean() - p_.mean()); t2["d"].append(S.cohens_d(s, p_)[0])
    out["table2"] = {k: np.array(v) for k, v in t2.items()}

    notch = {k: [] for k in ("fr", "depth", "bw", "q_loaded", "resolved_points")}
    for band in S.BANDS.values():
        for tr in [air] + list(pe) + list(sand):
            r = S.notch_parameters(freq, tr, band)
            for key in notch:
                notch[key].append(r[key])
    out["notch"] = {k: np.array(v, float) for k, v in notch.items()}

    sh = {"disp_hz": [], "sd_hz": [], "p": []}
    for fname in ("F2", "F4"):
        band = S.BANDS[fname]
        fr_pe = np.array([S.notch_parameters(freq, t, band)["fr"] for t in pe])
        fr_sa = np.array([S.notch_parameters(freq, t, band)["fr"] for t in sand])
        sh["disp_hz"].append(fr_pe.mean() - fr_sa.mean())
        sh["sd_hz"].append(np.std(fr_pe.mean() - fr_sa, ddof=1))
        sh["p"].append(stats.ttest_ind(fr_sa, fr_pe, equal_var=False).pvalue)
    out["shift"] = {k: np.array(v) for k, v in sh.items()}

    grad = np.gradient(air, freq / 1e6)
    out["gradient"] = {"grad": np.array([abs(grad[int(np.argmin(np.abs(freq - f0)))])
                                         for f0 in (2.341e9, 3.468e9, 4.596e9, 5.496e9)])}

    f3 = (freq >= S.BANDS["F3"][0]) & (freq <= S.BANDS["F3"][1])
    pe_min = np.array([t[f3].min() for t in pe]); sa_min = np.array([t[f3].min() for t in sand])
    b56 = (freq >= 5e9) & (freq <= 6e9)
    out["scalars"] = {"value": np.array([
        stats.ttest_ind(sa_min, pe_min, equal_var=False).pvalue,
        np.median(air[b56]), np.percentile(air[b56], 10)])}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rscript", default="Rscript")
    ap.add_argument("--data", default=S.DEFAULT_DATA)
    args = ap.parse_args()

    tmp = tempfile.mkdtemp(prefix="siw_cmp_")
    r_txt = os.path.join(tmp, "r.txt")
    py_txt = os.path.join(tmp, "py.txt")
    dump = os.path.join(tmp, "dump")

    subprocess.run([args.rscript, os.path.join(HERE, "siw_analysis.R"), "--data", args.data,
                    "--out", r_txt, "--dump", dump], check=True, stdout=subprocess.DEVNULL)
    subprocess.run([sys.executable, os.path.join(HERE, "siw_analysis.py"), "--data", args.data,
                    "--out", py_txt], check=True, stdout=subprocess.DEVNULL)

    r_lines = open(r_txt).read().splitlines()
    p_lines = open(py_txt).read().splitlines()
    diffs = [(i + 1, a, b) for i, (a, b) in enumerate(zip(p_lines, r_lines)) if a != b]
    print("=" * 78)
    print("1. Printed reports")
    print("=" * 78)
    print(f"Python lines: {len(p_lines)}   R lines: {len(r_lines)}   differing lines: {len(diffs)}")
    for i, a, b in diffs:
        print(f"  line {i}\n    py: {a}\n    R : {b}")
    print("  -> IDENTICAL" if not diffs and len(r_lines) == len(p_lines) else "  -> DIFFERENT")
    print()

    ref = python_reference(args.data)
    print("=" * 78)
    print("2. Full-precision values (Python reference vs R --dump)")
    print("=" * 78)
    print(f"{'quantity':<28}{'n':>5}{'max |abs diff|':>18}{'max rel diff':>16}")
    worst = 0.0
    for table, cols in ref.items():
        rows = read_csv(os.path.join(dump, f"{table}.csv"))
        for col, py in cols.items():
            r = num(rows, col)
            both_nan = np.isnan(py) & np.isnan(r)
            if np.any(np.isnan(py) ^ np.isnan(r)):
                print(f"{table + '.' + col:<28}{len(py):>5}   NaN pattern differs!")
                worst = np.inf
                continue
            a = np.abs(py - r)[~both_nan]
            rel = a / np.maximum(np.abs(py[~both_nan]), 1e-300)
            worst = max(worst, rel.max() if rel.size else 0.0)
            print(f"{table + '.' + col:<28}{len(py):>5}{a.max():18.3e}{rel.max():16.3e}")
    print()
    print(f"Largest relative difference over all quantities: {worst:.3e}")


if __name__ == "__main__":
    main()
