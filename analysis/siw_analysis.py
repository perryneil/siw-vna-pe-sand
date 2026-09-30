#!/usr/bin/env python3
"""
siw_analysis.py
===============

Reproduces every statistic reported in

    "Higher-order microwave resonances discriminate polyethylene
     microplastic from sand"

directly from the eleven raw LibreVNA CSV exports.

Inputs
------
A directory containing the eleven exports (default: ../data/vna, relative to
this script):

    ses06062026_bat0001_s11s21_empty.csv     air baseline
    ses06062026_s11s21_PE_0001.csv ... _0005.csv
    ses06062026_s11s21_sand_0001.csv ... _0005.csv

Each file is a LibreVNA export with 501 linearly spaced points from
1 MHz to 6 GHz and columns for S11 and S21 magnitude (linear), phase,
real and imaginary parts.

Usage
-----
    python analysis/siw_analysis.py                      # from the repository root
    python analysis/siw_analysis.py --out results.txt    # also save a transcript
    python analysis/siw_analysis.py --data /path/to/csvs

Requires: numpy, scipy   (pip install numpy scipy)

Author: P.N.J. Fernandez
"""

import argparse
import glob
import os
import sys
import warnings
from itertools import combinations

import numpy as np
from scipy import stats

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------

# Fixed analysis frequencies (Hz). These are exact sample points of the sweep.
FEATURES = {
    "F1 min":    2.3646e9,
    "F2 flank1": 3.4204e9,
    "F2 flank2": 3.4444e9,
    "F2 flank3": 3.4684e9,
    "F2 min":    3.5044e9,
    "F3 dip":    4.6322e9,
    "F4 null":   5.5081e9,
}

# Search bands used to locate each notch minimum per trial (Hz).
BANDS = {
    "F1": (2.20e9, 2.50e9),
    "F2": (3.40e9, 3.60e9),
    "F3": (4.50e9, 4.75e9),
    "F4": (5.40e9, 5.62e9),
}

ALPHA = 0.05

DEFAULT_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "vna")

# At the 1 MHz sweep edge the five PE and five sand traces are nearly identical,
# which makes SciPy warn about precision loss. The warning does not affect any
# reported value (that point has p = 1 after correction), so it is silenced.
warnings.filterwarnings("ignore", message="Precision loss occurred")


# --------------------------------------------------------------------------
# I/O
# --------------------------------------------------------------------------

def load_export(path):
    """Return (freq_Hz, s21_dB, s21_phase_deg, s21_real, s21_imag)."""
    with open(path, "r") as fh:
        header = fh.readline().rstrip("\n").rstrip(",").split(",")
    col = {name.strip(): i for i, name in enumerate(header)}

    data = np.genfromtxt(path, delimiter=",", skip_header=1)

    f = data[:, col["Frequency"]]
    mag = data[:, col["S21_Magnitude (linear)"]]
    with np.errstate(divide="ignore"):
        s21_db = 20.0 * np.log10(mag)
    phase = data[:, col["S21_Phase"]]
    real = data[:, col["S21_Real"]]
    imag = data[:, col["S21_Imaginary"]]
    return f, s21_db, phase, real, imag


def load_all(data_dir):
    air_path = os.path.join(data_dir, "ses06062026_bat0001_s11s21_empty.csv")
    pe_paths = sorted(glob.glob(os.path.join(data_dir, "*_PE_*.csv")))
    sand_paths = sorted(glob.glob(os.path.join(data_dir, "*_sand_*.csv")))

    if not os.path.exists(air_path):
        sys.exit(f"Air baseline not found: {air_path}")
    if len(pe_paths) != 5 or len(sand_paths) != 5:
        sys.exit(f"Expected 5 PE and 5 sand files, found "
                 f"{len(pe_paths)} and {len(sand_paths)} in {data_dir}")

    freq, air_db, air_phase, air_re, air_im = load_export(air_path)
    pe = np.array([load_export(p)[1] for p in pe_paths])
    sand = np.array([load_export(p)[1] for p in sand_paths])

    # every export must share the same frequency grid
    for p in pe_paths + sand_paths:
        assert np.allclose(load_export(p)[0], freq), f"grid mismatch in {p}"

    return dict(freq=freq, air=air_db, air_phase=air_phase,
                air_re=air_re, air_im=air_im,
                pe=pe, sand=sand,
                pe_files=[os.path.basename(p) for p in pe_paths],
                sand_files=[os.path.basename(p) for p in sand_paths])


# --------------------------------------------------------------------------
# Multiplicity corrections
# --------------------------------------------------------------------------

def holm_bonferroni(p):
    """Holm step-down adjusted p-values (family-wise error rate)."""
    p = np.asarray(p, float)
    n = p.size
    order = np.argsort(p)
    adj = np.empty(n)
    running = 0.0
    for rank, idx in enumerate(order):
        val = (n - rank) * p[idx]
        running = max(running, val)
        adj[idx] = min(running, 1.0)
    return adj


def benjamini_hochberg(p):
    """Benjamini-Hochberg adjusted p-values (false discovery rate)."""
    p = np.asarray(p, float)
    n = p.size
    order = np.argsort(p)
    adj = np.empty(n)
    running = 1.0
    for rank in range(n - 1, -1, -1):
        idx = order[rank]
        val = n * p[idx] / (rank + 1)
        running = min(running, val)
        adj[idx] = min(running, 1.0)
    return adj


def cohens_d(a, b):
    """Standardised effect size with pooled SD (sample SD, ddof=1)."""
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * np.var(a, ddof=1) +
                  (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    return abs(np.mean(a) - np.mean(b)) / sp, sp


# --------------------------------------------------------------------------
# Notch geometry
# --------------------------------------------------------------------------

def notch_parameters(freq, trace_db, band):
    """
    Locate the notch minimum inside `band` and characterise it.

    The 3 dB reference level is taken 3 dB below the mean of the two
    shoulder maxima bounding the notch, following the convention stated
    in the Methods. Band edges are found by linear interpolation of the
    crossings, so the returned width is not quantised to the sample grid.

    Returns dict with fr, depth, bw, q_loaded, resolved_points.
    """
    lo, hi = band
    sel = (freq >= lo) & (freq <= hi)
    idx = np.where(sel)[0]
    local = trace_db[idx]
    imin_local = int(np.argmin(local))
    imin = idx[imin_local]

    # walk outwards to the first local maximum on each side
    i = imin
    while i - 1 >= 0 and trace_db[i - 1] > trace_db[i]:
        i -= 1
    left_shoulder = i
    i = imin
    while i + 1 < len(trace_db) and trace_db[i + 1] > trace_db[i]:
        i += 1
    right_shoulder = i

    shoulder_mean = 0.5 * (trace_db[left_shoulder] + trace_db[right_shoulder])
    depth = shoulder_mean - trace_db[imin]
    level = shoulder_mean - 3.0

    def crossing(start, step):
        j = start
        while 0 <= j + step < len(trace_db) and trace_db[j] < level:
            j += step
        # linear interpolation between j and j-step
        j0, j1 = j - step, j
        y0, y1 = trace_db[j0], trace_db[j1]
        if y1 == y0:
            return freq[j1]
        t = (level - y0) / (y1 - y0)
        return freq[j0] + t * (freq[j1] - freq[j0])

    if depth <= 3.0:
        bw = np.nan
        q = np.nan
        npts = 0
    else:
        f_lo = crossing(imin, -1)
        f_hi = crossing(imin, +1)
        bw = f_hi - f_lo
        q = freq[imin] / bw if bw > 0 else np.nan
        npts = int(np.sum(trace_db[left_shoulder:right_shoulder + 1] < level))

    return dict(fr=freq[imin], idx=imin, depth=depth, bw=bw,
                q_loaded=q, resolved_points=npts,
                shoulder_mean=shoulder_mean)


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------

def fmt_p(p):
    if p >= 0.001:
        return f"{p:.3f}"
    return f"{p:.2e}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=DEFAULT_DATA, help="directory holding the 11 CSV exports")
    ap.add_argument("--out", default=None, help="write the report to this file as well")
    args = ap.parse_args()

    D = load_all(args.data)
    freq, air, pe, sand = D["freq"], D["air"], D["pe"], D["sand"]

    lines = []
    def say(s=""):
        lines.append(s)
        print(s)

    step = float(np.median(np.diff(freq)))
    say("=" * 78)
    say("SIW-VNA microplastic discrimination: full statistical reproduction")
    say("=" * 78)
    say(f"Points per sweep      : {len(freq)}")
    say(f"Span                  : {freq[0]/1e6:.3f} MHz to {freq[-1]/1e9:.3f} GHz")
    say(f"Median point spacing  : {step/1e6:.3f} MHz")
    say(f"PE trials             : {len(pe)}  {D['pe_files']}")
    say(f"Sand trials           : {len(sand)}  {D['sand_files']}")
    say("")

    # ------------------------------------------------------------------
    # 1. Point-by-point Welch tests over the whole sweep, with corrections
    # ------------------------------------------------------------------
    say("-" * 78)
    say("1. Point-wise Welch t-tests across all frequency points")
    say("-" * 78)

    raw_p = np.empty(len(freq))
    for k in range(len(freq)):
        raw_p[k] = stats.ttest_ind(sand[:, k], pe[:, k], equal_var=False).pvalue
    raw_p = np.nan_to_num(raw_p, nan=1.0)

    holm = holm_bonferroni(raw_p)
    bh = benjamini_hochberg(raw_p)

    say(f"Raw p < 0.05                 : {int(np.sum(raw_p < ALPHA))} points")
    say(f"Holm-Bonferroni p_adj < 0.05 : {int(np.sum(holm < ALPHA))} points")
    say(f"Benjamini-Hochberg q < 0.05  : {int(np.sum(bh < ALPHA))} points")
    say(f"Benjamini-Hochberg q < 0.01  : {int(np.sum(bh < 0.01))} points")
    say("")
    say("Holm survivors (p_adj < 0.05):")
    say(f"  {'Freq (GHz)':>12} {'sand-PE (dB)':>13} {'d':>6} {'raw p':>10} {'Holm':>8} {'BH q':>8}")
    for k in np.where(holm < ALPHA)[0]:
        d, _ = cohens_d(sand[:, k], pe[:, k])
        say(f"  {freq[k]/1e9:12.4f} {np.mean(sand[:,k])-np.mean(pe[:,k]):13.3f} "
            f"{d:6.2f} {fmt_p(raw_p[k]):>10} {holm[k]:8.3f} {bh[k]:8.4f}")
    say("")

    # ------------------------------------------------------------------
    # 2. Table 2: replicate statistics at the analysed features
    # ------------------------------------------------------------------
    say("-" * 78)
    say("2. Table 2 - replicate statistics by spectral feature")
    say("-" * 78)
    say(f"{'Feature':<11}{'GHz':>9}{'air':>9}{'sand mean+-SD':>19}{'PE mean+-SD':>19}"
        f"{'diff':>8}{'d':>6}{'raw p':>10}{'Holm':>9}{'BH q':>9}")
    for name, f0 in FEATURES.items():
        k = int(np.argmin(np.abs(freq - f0)))
        s, p_ = sand[:, k], pe[:, k]
        d, _ = cohens_d(s, p_)
        diff = np.mean(s) - np.mean(p_)
        say(f"{name:<11}{freq[k]/1e9:9.4f}{air[k]:9.2f}"
            f"{np.mean(s):11.2f} +-{np.std(s, ddof=1):5.2f}"
            f"{np.mean(p_):11.2f} +-{np.std(p_, ddof=1):5.2f}"
            f"{diff:+8.2f}{d:6.2f}{fmt_p(raw_p[k]):>10}"
            f"{holm[k]:9.3f}{bh[k]:9.4f}")
    say("")
    say("Note: the F3 row above is the point-wise test at 4.6322 GHz. The paper")
    say("      also reports a per-trial dip-depth test for F3, computed below.")
    say("")

    # F3 dip-depth test (a single derived test, outside the 501-point family)
    pe_depth = np.array([notch_parameters(freq, t, BANDS["F3"])["depth"] for t in pe])
    sand_depth = np.array([notch_parameters(freq, t, BANDS["F3"])["depth"] for t in sand])
    # dip depth relative to the trace minimum value itself
    pe_min = np.array([t[np.argmin(np.where((freq >= BANDS['F3'][0]) & (freq <= BANDS['F3'][1]),
                                            t, np.inf))] for t in pe])
    sand_min = np.array([t[np.argmin(np.where((freq >= BANDS['F3'][0]) & (freq <= BANDS['F3'][1]),
                                              t, np.inf))] for t in sand])
    t_dip = stats.ttest_ind(sand_min, pe_min, equal_var=False)
    say(f"F3 per-trial minimum depth: sand {np.mean(sand_min):.2f} +- {np.std(sand_min, ddof=1):.2f} dB, "
        f"PE {np.mean(pe_min):.2f} +- {np.std(pe_min, ddof=1):.2f} dB, "
        f"Welch p = {t_dip.pvalue:.3f}")
    say("")

    # ------------------------------------------------------------------
    # 3. Table 3: loaded quality factors and notch parameters
    # ------------------------------------------------------------------
    say("-" * 78)
    say("3. Table 3 - loaded quality factors and notch parameters")
    say("-" * 78)
    say(f"{'Feature':<8}{'Load':<7}{'fr (GHz)':>12}{'depth (dB)':>14}{'Q_L':>14}{'pts <3dB':>10}")
    for fname, band in BANDS.items():
        a = notch_parameters(freq, air, band)
        say(f"{fname:<8}{'air':<7}{a['fr']/1e9:12.4f}{a['depth']:14.2f}"
            f"{a['q_loaded']:14.1f}{a['resolved_points']:10d}")
        for label, group in (("PE", pe), ("sand", sand)):
            res = [notch_parameters(freq, t, band) for t in group]
            fr = np.array([r["fr"] for r in res]) / 1e9
            dep = np.array([r["depth"] for r in res])
            q = np.array([r["q_loaded"] for r in res])
            say(f"{'':<8}{label:<7}{np.mean(fr):8.4f} +-{np.std(fr, ddof=1):4.4f}"
                f"{np.mean(dep):9.2f} +-{np.std(dep, ddof=1):5.2f}"
                f"{np.nanmean(q):9.1f} +-{np.nanstd(q, ddof=1):5.1f}"
                f"{'':>10}")
    say("")

    # ------------------------------------------------------------------
    # 4. Displacement of the loaded minimum
    # ------------------------------------------------------------------
    say("-" * 78)
    say("4. Displacement of the loaded transmission minimum")
    say("-" * 78)
    for fname in ("F2", "F4"):
        band = BANDS[fname]
        fr_pe = np.array([notch_parameters(freq, t, band)["fr"] for t in pe])
        fr_sand = np.array([notch_parameters(freq, t, band)["fr"] for t in sand])
        shift = fr_pe.mean() - fr_sand  # sand below PE is positive
        t_shift = stats.ttest_ind(fr_sand, fr_pe, equal_var=False)
        say(f"{fname}: PE minima  {np.unique(fr_pe/1e9)} GHz")
        say(f"{fname}: sand minima {np.unique(fr_sand/1e9)} GHz")
        say(f"{fname}: mean sand-below-PE displacement = "
            f"{np.mean(fr_pe) - np.mean(fr_sand):.3e} Hz "
            f"({(np.mean(fr_pe)-np.mean(fr_sand))/1e6:.1f} +- {np.std(shift, ddof=1)/1e6:.1f} MHz), "
            f"Welch p = {t_shift.pvalue:.4f}")
        say(f"{fname}: one sweep step = {step/1e6:.3f} MHz -- displacement is "
            f"{(np.mean(fr_pe)-np.mean(fr_sand))/step:.2f} steps")
    say("")

    # ------------------------------------------------------------------
    # 5. Flank steepness on the air baseline
    # ------------------------------------------------------------------
    say("-" * 78)
    say("5. Flank steepness |dS21/df| on the air baseline")
    say("-" * 78)
    grad = np.gradient(air, freq / 1e6)  # dB per MHz
    for label, f0 in (("F1", 2.341e9), ("F2", 3.468e9), ("F3", 4.596e9), ("F4", 5.496e9)):
        k = int(np.argmin(np.abs(freq - f0)))
        say(f"{label} flank at {freq[k]/1e9:.4f} GHz: |dS21/df| = {abs(grad[k]):.3f} dB/MHz")
    say("")

    # ------------------------------------------------------------------
    # 6. The F4 transmission zero: phase diagnostic
    # ------------------------------------------------------------------
    say("-" * 78)
    say("6. F4 diagnostic - is the deep air null a transmission zero?")
    say("-" * 78)
    k = int(np.argmin(np.abs(freq - 5.5081e9)))
    say(f"Air |S21| at {freq[k]/1e9:.4f} GHz : {air[k]:.2f} dB")
    say(f"  neighbour below ({freq[k-1]/1e9:.4f} GHz): {air[k-1]:.2f} dB")
    say(f"  neighbour above ({freq[k+1]/1e9:.4f} GHz): {air[k+1]:.2f} dB")
    ph = D["air_phase"]
    say(f"Phase: {ph[k-2]:.1f} -> {ph[k-1]:.1f} -> {ph[k]:.1f} -> {ph[k+1]:.1f} -> {ph[k+2]:.1f} deg")
    say(f"  reversal across the null = {abs(ph[k+1] - ph[k-1]):.1f} deg")
    say(f"Real/imag at the null: Re = {D['air_re'][k]:.3e}, Im = {D['air_im'][k]:.3e}")
    band = (freq >= 5.0e9) & (freq <= 6.0e9)
    say(f"5-6 GHz band median = {np.median(air[band]):.1f} dB, "
        f"10th percentile = {np.percentile(air[band], 10):.1f} dB "
        f"(no noise-floor plateau)")
    say("")

    # ------------------------------------------------------------------
    # 7. Class separation and leave-one-out classification
    # ------------------------------------------------------------------
    say("-" * 78)
    say("7. Class separation, LOO nearest-centroid, permutation test")
    say("-" * 78)
    for name, f0 in (("F2 flank 3.4684 GHz", 3.4684e9), ("F4 null 5.5081 GHz", 5.5081e9)):
        k = int(np.argmin(np.abs(freq - f0)))
        s, p_ = sand[:, k], pe[:, k]
        gap = max(s.min(), p_.min()) - min(s.max(), p_.max())
        say(f"{name}")
        say(f"  PE range   : {p_.min():.2f} to {p_.max():.2f} dB")
        say(f"  sand range : {s.min():.2f} to {s.max():.2f} dB")
        say(f"  gap between class ranges : {abs(gap):.2f} dB "
            f"({'separated' if gap > 0 else 'OVERLAP'})")

        # leave-one-out nearest centroid on this single feature
        X = np.concatenate([s, p_])
        y = np.array([0] * len(s) + [1] * len(p_))
        correct = 0
        for i in range(len(X)):
            mask = np.ones(len(X), bool)
            mask[i] = False
            c0 = X[mask & (y == 0)].mean()
            c1 = X[mask & (y == 1)].mean()
            pred = 0 if abs(X[i] - c0) < abs(X[i] - c1) else 1
            correct += int(pred == y[i])
        say(f"  LOO nearest-centroid accuracy : {correct}/{len(X)}")
    n_splits = len(list(combinations(range(10), 5)))
    say(f"Permutation p-value for complete separation with n = 5 per class: "
        f"1/C(10,5) = 1/{n_splits} = {1/n_splits:.2e}")
    say("")
    say("=" * 78)
    say("Reminder: the two features above were selected using these same ten")
    say("measurements, so the leave-one-out figure is optimistically biased and")
    say("is a restatement of the observed separation, not an out-of-sample estimate.")
    say("=" * 78)

    if args.out:
        with open(args.out, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print(f"\n[written to {args.out}]")


if __name__ == "__main__":
    main()
