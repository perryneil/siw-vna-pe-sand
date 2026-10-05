#!/usr/bin/env python3
"""
Figure 3 - air-baseline |S21| from the CST simulation and from the LibreVNA
measurement, overlaid.

Each trace is drawn on its own frequency axis, exactly as exported: no
resampling, alignment or interpolation. Annotated:
  - the three simulated resonance peaks below 5 GHz (blue dots),
  - the measured maxima near 2.00, 2.65 and 4.50 GHz (open triangles),
  - the measured minima F1-F4 (red triangles), located with the same search
    bands as analysis/siw_analysis.py.

Inputs:
  data/simulation/S21-air-cst-v2.txt            CST export (GHz, dB)
  data/vna/ses06062026_bat0001_s11s21_empty.csv  measured air baseline

Usage (from the repository root):
    python figures/make_figure3.py        -> figures/output/Figure_3.png
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

from common import AIR_FILE, REPO, VNA_DIR, load_export, out_path

CST_FILE = os.path.join(REPO, "data", "simulation", "S21-air-cst-v2.txt")

# Search bands (GHz) for the measured minima, as in analysis/siw_analysis.py
MINIMA_BANDS = {"F1": (2.20, 2.50), "F2": (3.40, 3.60), "F3": (4.50, 4.75), "F4": (5.40, 5.62)}
# Search bands (GHz) for the measured maxima discussed in the Results
MAXIMA_BANDS = [(1.90, 2.10), (2.55, 2.75), (4.40, 4.60)]

BLUE, GREY, RED = "#1f77b4", "#2b2b2b", "#d62728"

# Label positions (GHz, dB) for the F labels, chosen to keep text clear of the traces
F_LABEL_XY = {"F1": (2.38, -24.5), "F2": (3.21, -40.5), "F3": (4.63, -51.5), "F4": (4.85, -62.0)}
SIM_LABEL_XY = [(1.90, -19.5), (2.98, -13.5), (4.72, -10.5)]


def local_extreme(f, y, lo, hi, fn):
    m = np.where((f >= lo) & (f <= hi))[0]
    return m[fn(y[m])]


def main():
    cst = np.loadtxt(CST_FILE, comments="#")
    f_sim, s_sim = cst[:, 0], cst[:, 1]
    air = load_export(os.path.join(VNA_DIR, AIR_FILE))
    f_meas, s_meas = air["f"] / 1e9, air["db"]

    sim_peaks, _ = find_peaks(s_sim, prominence=3)
    sim_peaks = [i for i in sim_peaks if f_sim[i] < 5.0]
    maxima = [local_extreme(f_meas, s_meas, lo, hi, np.argmax) for lo, hi in MAXIMA_BANDS]
    minima = {k: local_extreme(f_meas, s_meas, lo, hi, np.argmin) for k, (lo, hi) in MINIMA_BANDS.items()}

    plt.rcParams.update({"font.size": 10})
    fig, ax = plt.subplots(figsize=(7.1, 3.8))
    ax.grid(color="0.88", lw=0.6)
    ax.plot(f_sim, s_sim, color=BLUE, lw=1.8, label="CST simulation", zorder=3)
    ax.plot(f_meas, s_meas, color=GREY, lw=1.0, label="LibreVNA measurement", zorder=2)

    for i, (tx, ty) in zip(sim_peaks, SIM_LABEL_XY):
        ax.plot(f_sim[i], s_sim[i], "o", color=BLUE, ms=5, zorder=4)
        ax.text(tx, ty, f"{f_sim[i]:.4f} GHz", color=BLUE, fontsize=9)

    for i in maxima:
        ax.plot(f_meas[i], s_meas[i], "^", mfc="none", mec="0.35", ms=6, zorder=4)
    ax.text(1.6, 3.5, "△  measured maxima: " + ", ".join(f"{f_meas[i]:.2f}" for i in maxima) + " GHz",
            color="0.35", fontsize=9, va="center")

    for name, i in minima.items():
        ax.plot(f_meas[i], s_meas[i], "v", color=RED, ms=6, zorder=5)
        tx, ty = F_LABEL_XY[name]
        ax.annotate(f"{name}  {f_meas[i]:.3f} GHz", xy=(f_meas[i], s_meas[i]), xytext=(tx, ty),
                    color="#b0202a", fontsize=9,
                    arrowprops=dict(arrowstyle="-", color="#b0202a", lw=0.8, shrinkB=4))

    ax.set_xlim(1.5, 6.0)
    ax.set_ylim(-78, 8)
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("|$S_{21}$| (dB)")
    ax.legend(loc="lower left", framealpha=0.95)

    fig.savefig(out_path("Figure_3.png"), dpi=300, bbox_inches="tight")
    print("written", out_path("Figure_3.png"))
    print("simulated peaks (GHz):", [round(f_sim[i], 4) for i in sim_peaks])
    print("measured maxima (GHz):", [round(f_meas[i], 4) for i in maxima])
    print("measured minima (GHz):", {k: round(f_meas[i], 4) for k, i in minima.items()})


if __name__ == "__main__":
    main()
