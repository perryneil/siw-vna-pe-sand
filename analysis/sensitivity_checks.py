#!/usr/bin/env python3
"""
sensitivity_checks.py
=====================

Calibration-integrity, measurement-order and robustness checks reported in the
section "Calibration integrity and measurement-order effects" of

    "Higher-order microwave resonances discriminate polyethylene
     microplastic from sand"

It complements analysis/siw_analysis.py (which reproduces Tables 4 and 5) and
recomputes, from the raw files:

  1. the S-parameters present in the exports;
  2. the passivity check (|S11| > 1) on every export;
  3. VSWR at the four notches (from |S11|);
  4. the port error terms implied by the SOLT calibration file;
  5. drift with measurement order (air, sand_0001-0005, then PE_0001-0005);
  6. the feature differences adjusted for that drift;
  7. robustness to the step change between PE_0003 and PE_0004;
  8. exact permutation limits, confidence intervals for Cohen's d and the
     exact test for the one-step displacement of the minimum;
  9. the simulated resonances and their 3 dB loaded Q from the CST export.

Usage (from the repository root):
    python analysis/sensitivity_checks.py
    python analysis/sensitivity_checks.py --out results_sensitivity.txt

Compare the transcript with analysis/expected_output_sensitivity.txt.
Requires: numpy, scipy.
"""

import argparse
import datetime
import json
import os
import sys
import warnings
from itertools import combinations

import numpy as np
from scipy import stats
from scipy.optimize import brentq
from scipy.signal import find_peaks

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "..")
VNA = os.path.join(REPO, "data", "vna")
CAL = os.path.join(REPO, "data", "calibration", "SOLT_12_1_00M-6_00G_501pt.cal")
SIM = os.path.join(REPO, "data", "simulation", "S21-air-cst-v2.txt")
AIR = "ses06062026_bat0001_s11s21_empty.csv"
ORDER = [("sand", i) for i in range(1, 6)] + [("PE", i) for i in range(1, 6)]   # measurement order
NOTCHES = {"F1": 2.3646, "F2": 3.5044, "F3": 4.6322, "F4": 5.5081}            # GHz
TZ = datetime.timezone(datetime.timedelta(hours=8))

warnings.filterwarnings("ignore", message="Precision loss occurred")

OUT = []


def say(*a):
    line = " ".join(str(x) for x in a)
    print(line)
    OUT.append(line)


def section(title):
    say("")
    say("-" * 78)
    say(title)
    say("-" * 78)


def load(path):
    head = open(path).readline().strip().rstrip(",").split(",")
    col = {n.strip(): i for i, n in enumerate(head)}
    d = np.genfromtxt(path, delimiter=",", skip_header=1)
    s21 = d[:, col["S21_Magnitude (linear)"]]
    return head, dict(f=d[:, col["Frequency"]] / 1e9, s21lin=s21, s21=20 * np.log10(s21),
                      s11lin=d[:, col["S11_Magnitude (linear)"]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", help="also write the transcript to this file")
    args = ap.parse_args()

    head, air = load(os.path.join(VNA, AIR))
    tr = [load(os.path.join(VNA, f"ses06062026_s11s21_{m}_000{i}.csv"))[1] for m, i in ORDER]
    f = air["f"]
    Y = np.array([t["s21"] for t in tr])
    S, P = Y[:5], Y[5:]
    lab = np.array([0] * 5 + [1] * 5)
    idx = lambda g: int(np.argmin(abs(f - g)))

    say("=" * 78)
    say("SIW-VNA: calibration, measurement-order and robustness checks")
    say("=" * 78)
    say("Measurement order: air baseline, sand_0001-0005, PE_0001-0005")

    section("1. S-parameters present in the exports")
    say("columns found for:", sorted({h.split("_")[0] for h in head if h.startswith("S")}))

    section("2. Passivity check: |S11| > 1 is impossible for a passive device")
    for name, t in [("air", air)] + [(f"{m}_{i}", x) for (m, i), x in zip(ORDER, tr)]:
        bad = t["s11lin"] > 1
        say(f"{name:7s} |S11| > 1 at {bad.sum():2d} points ({f[bad].min():.3f}-{f[bad].max():.3f} GHz); "
            f"max |S11| = {20*np.log10(t['s11lin'].max()):+.1f} dB")

    section("3. VSWR at the notch frequencies (from |S11|)")
    for k, g in NOTCHES.items():
        i = idx(g)
        vs = [(1 + t["s11lin"][i]) / (1 - t["s11lin"][i]) if t["s11lin"][i] < 1 else np.inf for t in [air] + tr]
        vs = np.array(vs)
        if np.all(np.isfinite(vs)):
            say(f"{k} {f[i]:.4f} GHz: air {vs[0]:.2f}; sand SD {vs[1:6].std(ddof=1):.2f} (range {np.ptp(vs[1:6]):.2f}); "
                f"PE SD {vs[6:].std(ddof=1):.2f} (range {np.ptp(vs[6:]):.2f})")
        else:
            say(f"{k} {f[i]:.4f} GHz: |S11| = {air['s11lin'][i]:.2f} on the air trace; VSWR undefined")

    section("4. Port error terms implied by the SOLT calibration file")
    cal = json.load(open(CAL))
    say("calibration:", cal.get("type"), "ports", cal.get("ports"), "| kit:", cal["calkit"]["Description"])
    stamps = [m["data"]["timestamp"] for m in cal["measurements"]]
    say("standards measured", datetime.datetime.fromtimestamp(min(stamps), TZ).strftime("%Y-%m-%d %H:%M:%S"),
        "to", datetime.datetime.fromtimestamp(max(stamps), TZ).strftime("%H:%M:%S"), "(UTC+8)")
    R = {}
    for m in cal["measurements"]:
        pts = m["data"]["points"]; fc = np.array([p["frequency"] for p in pts]) / 1e9
        if m["type"] != "Through":
            R[(m["type"], m["data"]["port"])] = np.array([p["real"] + 1j * p["imag"] for p in pts])
    for port in (1, 2):
        Sh, Op, Lo = R[("Short", port)], R[("Open", port)], R[("Load", port)]
        a, b = Sh - Lo, Op - Lo
        e11 = (a + b) / (b - a)
        bad = np.abs(e11) > 1
        say(f"port {port}: |e11| > 1 at {bad.sum()} points"
            + (f" ({fc[bad].min():.3f}-{fc[bad].max():.3f} GHz)" if bad.any() else ""))
        for g in (3.0, 3.5, 4.0, 4.3, 4.6, 4.8, 5.0, 5.5):
            j = np.argmin(abs(fc - g))
            say(f"   {fc[j]:.3f} GHz  |e11| = {abs(e11[j]):.2f}  angle((O-L)/(S-L)) = {np.degrees(np.angle(b[j]/a[j])):+5.0f} deg")
        dev = np.abs(180 - np.abs(np.degrees(np.angle(b / a))))
        say(f"   max |e11| = {np.abs(e11).max():.2f}; max departure of the open-short phase difference from 180 deg = {dev.max():.0f} deg")
        say("   at the reported features: " + "; ".join(
            f"{g:.4f} GHz |e11| = {abs(e11[np.argmin(abs(fc - g))]):.2f}"
            for g in (2.3646, 3.4204, 3.4444, 3.4684, 4.6322, 5.5081)))

    section("5. Drift with measurement order (mean S21 minus air, off-resonance bands)")
    bands = {"0.1-1": (0.1, 1), "1-2": (1, 2), "2-3.3": (2, 3.3)}
    X = np.array([[np.mean((y - air["s21"])[(f >= lo) & (f < hi)]) for lo, hi in bands.values()] for y in Y])
    for j, k in enumerate(bands):
        r, p = stats.spearmanr(np.arange(1, 11), X[:, j])
        say(f"{k:>6s} GHz: offsets in order {np.round(X[:, j], 3)}; Spearman rho = {r:+.2f}, p = {p:.3f}")

    section("6. Feature differences adjusted for drift (value ~ material + drift index)")
    say("Material and order are confounded (all sand before all PE), so the adjustment is a")
    say("sensitivity check. It is repeated for three definitions of the drift index, each the mean")
    say("of S21 - air over a band without discriminating features.")
    for lo_, hi_ in ((0.1, 2.0), (1.0, 2.0), (1.0, 3.3)):
        drift = np.array([np.mean((y - air["s21"])[(f >= lo_) & (f < hi_)]) for y in Y])
        say(f"drift index over {lo_}-{hi_} GHz: correlation with material r = {np.corrcoef(lab, drift)[0, 1]:+.2f}")
        for g in (2.3646, 3.4204, 3.4444, 3.4684, 4.6322, 5.5081):
            i = idx(g); v = Y[:, i]
            A = np.c_[np.ones(10), lab, drift]; beta = np.linalg.lstsq(A, v, rcond=None)[0]
            e = v - A @ beta; cov = (e @ e / 7) * np.linalg.inv(A.T @ A); t = beta[1] / np.sqrt(cov[1, 1])
            raw = S[:, i].mean() - P[:, i].mean()
            say(f"   {f[i]:.4f} GHz: sand - PE raw {raw:+.2f} dB, adjusted {-beta[1]:+.2f} dB "
                f"({-beta[1]/raw:5.0%} of raw), p = {2*stats.t.sf(abs(t), 7):.3g}")
    say("Trend within the sand trials (measured first), slope per trial and Spearman rho:")
    for g in (3.4204, 3.4684, 5.5081):
        i = idx(g); sl = np.polyfit(np.arange(5), S[:, i], 1)[0]
        say(f"   {f[i]:.4f} GHz: slope {sl:+.3f} dB per trial, rho = {stats.spearmanr(np.arange(5), S[:, i])[0]:+.2f}; "
            f"sand - PE = {S[:, i].mean()-P[:, i].mean():+.2f} dB")

    section("7. Robustness to the step between PE_0003 and PE_0004")
    for g in (2.3646, 3.4204, 3.4444, 3.4684, 4.6322, 5.5081):
        i = idx(g)
        for name, sl in (("PE 1-3", slice(0, 3)), ("PE 4-5", slice(3, 5))):
            a, b = S[:, i], P[sl, i]
            sep = a.max() < b.min() or b.max() < a.min()
            say(f"{f[i]:.4f} GHz sand vs {name}: sand - PE {a.mean()-b.mean():+.2f} dB, "
                f"Welch p {stats.ttest_ind(a, b, equal_var=False)[1]:.2g}, ranges separated: {sep}")

    section("8. Small-sample inference")
    for g in (3.4684, 5.5081):
        i = idx(g); x = np.r_[S[:, i], P[:, i]]; obs = abs(S[:, i].mean() - P[:, i].mean()); c = n = 0
        for comb in combinations(range(10), 5):
            sel = list(comb); n += 1
            c += abs(x[sel].mean() - np.delete(x, sel).mean()) >= obs - 1e-12
        say(f"{f[i]:.4f} GHz: exact two-sided permutation p = {c}/{n} = {c/n:.4f} "
            f"(Holm threshold for 501 tests = {0.05/501:.1e})")

    def d_ci(d, n1=5, n2=5):
        se = np.sqrt(1 / n1 + 1 / n2); tt = d / se; df = n1 + n2 - 2
        lo = brentq(lambda nc: stats.nct.sf(tt, df, nc) - 0.025, -50, 50)
        hi = brentq(lambda nc: stats.nct.sf(tt, df, nc) - 0.975, -50, 60)
        return lo * se, hi * se

    for g in (3.4204, 3.4444, 3.4684, 3.4804, 5.5081):
        i = idx(g); a, b = S[:, i], P[:, i]
        d = abs(a.mean() - b.mean()) / np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
        lo, hi = d_ci(d)
        say(f"{f[i]:.4f} GHz: Cohen's d = {d:.2f}, 95% CI {lo:.2f} to {hi:.2f} (noncentral t)")
    for k, (lo_, hi_) in {"F2": (3.40, 3.60), "F4": (5.40, 5.62)}.items():
        m = (f >= lo_) & (f <= hi_)
        f_air = f[m][np.argmin(air["s21"][m])]          # air-baseline minimum frequency
        low = lambda X: int(np.sum(f[m][np.argmin(X[:, m], axis=1)] < f_air - 0.5 * np.median(np.diff(f))))
        ns, npe = low(S), low(P)
        p = stats.fisher_exact([[ns, 5 - ns], [npe, 5 - npe]])[1]
        say(f"{k}: minimum below the air minimum ({f_air:.4f} GHz) in {ns}/5 sand and {npe}/5 PE trials; "
            f"Fisher exact two-sided p = {p:.3f}")

    section("9. Simulated resonances (CST export) and 3 dB loaded Q")
    cst = np.loadtxt(SIM, comments="#"); fs, ss = cst[:, 0], cst[:, 1]
    for i in find_peaks(ss, prominence=3)[0]:
        lvl = ss[i] - 3; l = r = i
        while l > 0 and ss[l] > lvl:
            l -= 1
        while r < len(ss) - 1 and ss[r] > lvl:
            r += 1
        fl = np.interp(lvl, [ss[l], ss[l + 1]], [fs[l], fs[l + 1]])
        fr = np.interp(lvl, [ss[r], ss[r - 1]], [fs[r], fs[r - 1]])
        say(f"peak {fs[i]:.4f} GHz, {ss[i]:.2f} dB, 3 dB bandwidth {(fr-fl)*1e3:.0f} MHz, Q_L = {fs[i]/(fr-fl):.1f}")

    if args.out:
        with open(args.out, "w") as fh:
            fh.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    sys.exit(main())
