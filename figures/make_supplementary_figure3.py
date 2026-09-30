#!/usr/bin/env python3
"""
Supplementary figure 3 - baseline-corrected Raman spectrum of the PE fraction.

The spectrum is baseline-corrected by asymmetric least squares (same settings
as characterisation/raman_band_analysis.py), normalised to its strongest band,
and every peak with a prominence of at least six times the noise is labelled.
The 1418 cm-1 band of the orthorhombic crystalline phase is starred. Dotted
lines mark where polypropylene (PP), polystyrene (PS) and calcite (CaCO3)
bands would appear.

Usage (from the repository root):
    python figures/make_supplementary_figure3.py
        -> figures/output/Supplementary_figure_3.png
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

from common import RAMAN_DIR, out_path

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "characterisation"))
from raman_band_analysis import (PROMINENCE_SIGMA, als_baseline,  # noqa: E402
                                 load_spectrum, noise_level)

MARKERS = [("PP 809", 809), ("PS 1001", 1001), ("CaCO$_3$ 1086", 1086), ("PS 1602", 1602)]
CRYSTALLINE = 1418
LABEL_MIN = 0.05          # label peaks above 5% of the strongest band
LABEL_SHIFT = {1419: -22, 1460: 22}   # cm-1, keeps crowded labels apart
BLUE = "#2a78d6"


def main():
    x, y = load_spectrum(os.path.join(RAMAN_DIR, "Php-MP-Raman-PE.txt"))
    c = y - als_baseline(y)
    noise = noise_level(x, c)
    c_norm = c / c.max()
    peaks, _ = find_peaks(c, prominence=PROMINENCE_SIGMA * noise, distance=3)

    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                         "font.size": 9})
    fig, ax = plt.subplots(figsize=(6.3, 2.75))
    for label, pos in MARKERS:
        ax.axvline(pos, color="0.55", ls=":", lw=0.8)
        ax.text(pos - 6, 1.2, label, rotation=90, ha="right", va="top", fontsize=7, color="0.45")
    ax.plot(x, c_norm, color=BLUE, lw=0.9)

    # label the main peaks; merge a peak that sits within 5 cm-1 of a stronger one
    labelled = []
    for i in sorted(peaks, key=lambda j: -c_norm[j]):
        if c_norm[i] < LABEL_MIN or any(abs(x[i] - x[j]) < 5 for j in labelled):
            continue
        labelled.append(i)
        pos = int(round(x[i]))
        star = "*" if abs(x[i] - CRYSTALLINE) <= 4 else ""
        ax.text(x[i] + LABEL_SHIFT.get(pos, 0), c_norm[i] + 0.05, f"{pos}{star}",
                ha="center", va="bottom", fontsize=7)

    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(-0.04, 1.25)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_xlabel("Raman shift (cm$^{-1}$)")
    ax.set_ylabel("Normalised intensity")
    ax.spines[["top", "right"]].set_visible(False)

    fig.savefig(out_path("Supplementary_figure_3.png"), dpi=300, bbox_inches="tight")
    print("written", out_path("Supplementary_figure_3.png"))


if __name__ == "__main__":
    main()
