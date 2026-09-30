#!/usr/bin/env python3
"""
Figure 5 - diagnostic of the F4 feature in the air-baseline trace.

(a) |S21| across 5.40-5.62 GHz: the -67.2 dB minimum is a single sample of an
    unresolved transmission zero; the dotted line is the 5-6 GHz band median.
(b) S21 phase over the same span, showing the reversal across the null.
Markers are the 11.998 MHz sample points.

Usage (from the repository root):
    python figures/make_figure5.py        -> figures/output/Figure_5.png
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import AIR_FILE, VNA_DIR, load_export, out_path

F_LO, F_HI = 5.40e9, 5.62e9
F_NULL = 5.5081e9
RED = "#b5401f"


def main():
    air = load_export(os.path.join(VNA_DIR, AIR_FILE))
    f, db, ph = air["f"], air["db"], air["phase"]
    m = (f >= F_LO - 1e6) & (f <= F_HI)
    k = int(np.argmin(np.abs(f - F_NULL)))
    band = (f >= 5.0e9) & (f <= 6.0e9)
    median = np.median(db[band])

    plt.rcParams.update({
        "font.size": 9, "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
    })
    fig, (a, b) = plt.subplots(2, 1, figsize=(3.4, 3.8), sharex=True,
                               gridspec_kw=dict(hspace=0.08, height_ratios=[1, 0.95]))
    marker = dict(marker="o", ms=3.2, mfc="white", mec="#1a1a1a", mew=0.9)

    a.plot(f[m] / 1e9, db[m], color="#1a1a1a", lw=0.9, **marker)
    a.axhline(median, color="0.55", ls=":", lw=0.9)
    a.text(5.403, median - 1.8, f"5–6 GHz band median ({median:.1f} dB)".replace("-", "−"),
           fontsize=7, color="0.4", va="top")
    a.annotate(f"single sample\n{db[k]:.1f} dB".replace("-", "−"),
               xy=(f[k] / 1e9, db[k]), xytext=(5.543, -57),
               fontsize=7.5, color="0.2",
               arrowprops=dict(arrowstyle="-", color="0.3", lw=0.7))
    a.set_ylabel("|S21| (dB)")
    a.set_ylim(-72, -17)
    a.set_yticks(np.arange(-70, -19, 10))
    a.text(0.02, 0.03, "(a)", transform=a.transAxes, fontweight="bold")
    a.grid(color="0.88", lw=0.5)

    b.plot(f[m] / 1e9, ph[m], color="#1a1a1a", lw=0.9, **marker)
    b.axvline(f[k] / 1e9, color=RED, ls="--", lw=0.9)
    b.text(5.538, -30, "≈180° reversal", color=RED, fontsize=7.5, va="center")
    b.set_ylabel("S21 phase (deg)")
    b.set_ylim(-100, 100)
    b.set_yticks([-90, -45, 0, 45, 90])
    b.set_xlabel("Frequency (GHz)")
    b.text(0.02, 0.97, "(b)", transform=b.transAxes, fontweight="bold", va="top")
    b.grid(color="0.88", lw=0.5)
    b.set_xlim(5.389, 5.628)

    fig.savefig(out_path("Figure_5.png"), dpi=600, bbox_inches="tight")
    print("written", out_path("Figure_5.png"))


if __name__ == "__main__":
    main()
