#!/usr/bin/env python3
"""
uncertainty_budget.py
=====================

Uncertainty budget for the sand - PE differences in |S21| reported in table 7 of

    "Higher-frequency transmission features of a microwave cavity discriminate
     polyethylene microplastic from sand"

The measurand is the difference between the class means of |S21| (dB) at a common
sampled frequency, Delta = mean(sand) - mean(PE). Following JCGM 100:2008 (GUM):

  u_A  Type A, reload variability:   sqrt(s_sand^2/5 + s_PE^2/5), with
       Welch-Satterthwaite effective degrees of freedom.
  u_D  Type B, measurement-order drift. Material and order are confounded, so the
       drift adjustment cannot be applied as a correction. The largest change in
       Delta produced by any of the three drift-index definitions (the same
       ordinary least-squares model as sensitivity_checks.py) is taken as the
       half-width a of a rectangular distribution: u_D = a / sqrt(3).
  u_c  combined standard uncertainty, sqrt(u_A^2 + u_D^2); U = 2 u_c.

The frequency-grid contribution (rectangular over one 11.998 MHz step,
u = step / sqrt(12)) applies to feature positions, not to Delta at a common
sample frequency, and is printed separately. Calibration error is not evaluated
(see the manuscript, section 2.5).

Usage:  python analysis/uncertainty_budget.py      (from the repository root)
Requires: numpy, scipy
"""

import os
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
VNA = os.path.join(HERE, "..", "data", "vna")
AIR = "ses06062026_bat0001_s11s21_empty.csv"
ORDER = [("sand", i) for i in range(1, 6)] + [("PE", i) for i in range(1, 6)]
FEATURES = {"F1 min": 2.3646, "F2 flank": 3.4204, "F2 flank ": 3.4444,
            "F2 flank  ": 3.4684, "F4 null": 5.5081}
DRIFT_BANDS = ((0.1, 2.0), (1.0, 2.0), (1.0, 3.3))


def load(path):
    head = open(path).readline().strip().rstrip(",").split(",")
    col = {n.strip(): i for i, n in enumerate(head)}
    d = np.genfromtxt(path, delimiter=",", skip_header=1)
    return d[:, col["Frequency"]] / 1e9, 20 * np.log10(d[:, col["S21_Magnitude (linear)"]])


def main():
    f, air = load(os.path.join(VNA, AIR))
    Y = np.array([load(os.path.join(VNA, f"ses06062026_s11s21_{m}_000{i}.csv"))[1] for m, i in ORDER])
    S, P = Y[:5], Y[5:]
    lab = np.array([0] * 5 + [1] * 5)
    step = (f[1] - f[0]) * 1e3
    print(f"sample step {step:.3f} MHz; frequency-grid standard uncertainty "
          f"(rectangular over one step) {step/np.sqrt(12):.2f} MHz\n")
    print(f"{'feature':12s}{'f (GHz)':>9s}{'Delta':>8s}{'u_A':>7s}{'nu_eff':>8s}{'a':>7s}"
          f"{'u_D':>7s}{'u_c':>7s}{'U(k=2)':>8s}{'|D|/U':>7s}")
    for name, g in FEATURES.items():
        i = int(np.argmin(abs(f - g)))
        s, p = S[:, i], P[:, i]
        delta = s.mean() - p.mean()
        vs, vp = s.var(ddof=1) / 5, p.var(ddof=1) / 5
        uA = np.sqrt(vs + vp)
        nu = (vs + vp) ** 2 / (vs ** 2 / 4 + vp ** 2 / 4)
        adj = []
        for lo, hi in DRIFT_BANDS:
            drift = np.array([np.mean((y - air)[(f >= lo) & (f < hi)]) for y in Y])
            A = np.c_[np.ones(10), lab, drift]
            beta = np.linalg.lstsq(A, Y[:, i], rcond=None)[0]
            adj.append(-beta[1])
        a = max(abs(delta - x) for x in adj)
        uD = a / np.sqrt(3)
        uc = np.hypot(uA, uD)
        print(f"{name.strip():12s}{f[i]:9.4f}{delta:+8.2f}{uA:7.3f}{nu:8.1f}{a:7.2f}"
              f"{uD:7.3f}{uc:7.3f}{2*uc:8.2f}{abs(delta)/(2*uc):7.1f}")
    print("\nWith k = t_0.975(7) = %.2f instead of 2, every |Delta|/U for the F2 flank and F4 "
          "remains above 1." % stats.t.ppf(0.975, 7))


if __name__ == "__main__":
    main()
