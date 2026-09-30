#!/usr/bin/env python3
"""
Figure 4 - replicated transmission-domain perturbations for the air baseline,
PE microplastic and sand.

(a) full-band mean |S21| with +/-1 SD bands (n = 5 per loaded material),
    with the two repeatable discriminating regions highlighted;
(b) the F2 flank; (c) the F4 feature.

Usage (from the repository root):
    python figures/make_figure4.py        -> figures/output/Figure_4.png
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import load_all, out_path

C_AIR, C_PE, C_SAND = "#404040", "#1f77b4", "#d62728"
HIGHLIGHTS = [(3.42, 3.51), (5.46, 5.555)]   # GHz, discriminating regions


def band(ax, f, group, colour, label, lw=2.0):
    m, s = group.mean(axis=0), group.std(axis=0, ddof=1)
    ax.fill_between(f, m - s, m + s, color=colour, alpha=0.25, lw=0)
    ax.plot(f, m, color=colour, lw=lw, label=label)


def main():
    air, pe, sand = load_all()
    f = air["f"] / 1e9

    plt.rcParams.update({"font.size": 12})
    fig = plt.figure(figsize=(11.17, 7.95))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.62], hspace=0.32, wspace=0.22)

    ax = fig.add_subplot(gs[0, :])
    for lo, hi in HIGHLIGHTS:
        ax.axvspan(lo, hi, color="#fff3c4", alpha=0.9, lw=0, zorder=0)
    ax.plot(f, air["db"], color=C_AIR, lw=2.0, label="Air (baseline)")
    band(ax, f, pe, C_PE, "PE microplastic (mean, n=5)")
    band(ax, f, sand, C_SAND, "Sand (mean, n=5)")
    ax.set_xlim(1, 6)
    ax.set_ylim(-75, 2)
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("|S21| (dB)")
    ax.set_title("(a) Replicated transmission spectra with ±1 SD bands", loc="left")
    ax.legend(loc="lower left")

    zooms = [
        (gs[1, 0], (3.30, 3.58), "(b) Second-mode flank ≈ 3.45 GHz", "upper right"),
        (gs[1, 1], (5.30, 5.69), "(c) High-order feature ≈ 5.5 GHz", "lower left"),
    ]
    for spec, (lo, hi), title, loc in zooms:
        a = fig.add_subplot(spec)
        m = (f >= lo) & (f <= hi)
        a.plot(f[m], air["db"][m], color=C_AIR, lw=2.0, label="Air")
        band(a, f[m], pe[:, m], C_PE, "PE")
        band(a, f[m], sand[:, m], C_SAND, "Sand")
        a.set_xlabel("Frequency (GHz)")
        a.set_ylabel("|S21| (dB)")
        a.set_title(title, loc="left")
        a.legend(loc=loc)

    fig.savefig(out_path("Figure_4.png"), dpi=200, bbox_inches="tight")
    print("written", out_path("Figure_4.png"))


if __name__ == "__main__":
    main()
