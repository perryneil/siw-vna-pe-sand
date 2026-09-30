#!/usr/bin/env python3
"""
raman_band_analysis.py
======================

Band analysis of the Raman spectra used to verify the identity of the
polyethylene (PE) fraction in

    "Higher-order microwave resonances discriminate polyethylene
     microplastic from sand"

For the PE spectrum the script:

  1. removes the fluorescence background by asymmetric least squares (ALS),
  2. estimates the noise over the featureless 1680-1760 cm-1 region,
  3. finds peaks with a prominence of at least six times that noise,
  4. matches them (within +/- 8 cm-1) to literature bands of PE, PP and PS,
  5. checks marker bands of the other polymers and of calcite: whether a
     peak is detected there (same 6-sigma criterion), and the baseline-
     corrected intensity at that position as a % of the strongest band.

A marker counts as present only if a peak is detected. The raw intensity is
shown for context: in the PE spectrum the intensity near 1086 cm-1 is the
broad amorphous v(C-C) shoulder of PE itself (PE band near 1080 cm-1), not a
sharp calcite band.

Inputs
------
Php-MP-Raman-PE.txt (two tab-separated columns: Raman shift in cm-1,
intensity), in ../data/raman relative to this script by default.

Usage
-----
    python characterisation/raman_band_analysis.py
    python characterisation/raman_band_analysis.py --data /path/to/raman --out raman.txt

Requires: numpy, scipy
"""

import argparse
import os

import numpy as np
from scipy import sparse
from scipy.signal import find_peaks
from scipy.sparse.linalg import spsolve

DEFAULT_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "raman")

# Literature Raman bands (cm-1) inside the measured 660-1760 cm-1 window.
REFERENCE_BANDS = {
    "PE": [1062, 1080, 1130, 1170, 1295, 1370, 1416, 1440, 1460],
    "PP": [809, 841, 900, 941, 973, 998, 1040, 1103, 1152, 1168, 1220, 1256,
           1330, 1360, 1435, 1458],
    "PS": [795, 1001, 1031, 1155, 1183, 1197, 1330, 1450, 1583, 1602],
}

# Marker bands used to look for contamination by other materials.
MARKERS = {
    "PS 1001": 1001, "PS 1602": 1602,
    "PP 809": 809, "PP 841": 841, "PP 973": 973,
    "calcite 1086": 1086,
}

MATCH_TOL = 8.0            # cm-1, band assignment tolerance
NOISE_WINDOW = (1680, 1760)
PROMINENCE_SIGMA = 6.0


def load_spectrum(path):
    """Return (shift_cm1, intensity) sorted by increasing Raman shift."""
    d = np.loadtxt(path, comments="#")
    d = d[np.argsort(d[:, 0])]
    return d[:, 0], d[:, 1]


def als_baseline(y, lam=1e5, p=0.01, n_iter=15):
    """Asymmetric least-squares baseline (Eilers and Boelens, 2005)."""
    n = len(y)
    D = sparse.diags([1.0, -2.0, 1.0], [0, -1, -2], shape=(n, n - 2), dtype=float)
    D = lam * D.dot(D.T)
    w = np.ones(n)
    for _ in range(n_iter):
        W = sparse.spdiags(w, 0, n, n)
        z = spsolve(sparse.csc_matrix(W + D), w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z


def noise_level(x, c):
    """Noise from the point-to-point differences in the featureless window."""
    m = (x > NOISE_WINDOW[0]) & (x < NOISE_WINDOW[1])
    return np.std(np.diff(c[m])) / np.sqrt(2)


def marker_signal(x, c, centre, half_width=4.0):
    m = (x > centre - half_width) & (x < centre + half_width)
    return c[m].max()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=DEFAULT_DATA, help="folder holding Php-MP-Raman-PE.txt")
    ap.add_argument("--out", default=None, help="write the report to this file as well")
    args = ap.parse_args()

    lines = []

    def say(s=""):
        lines.append(s)
        print(s)

    spectra = {}
    for polymer in ("PE",):
        x, y = load_spectrum(os.path.join(args.data, f"Php-MP-Raman-{polymer}.txt"))
        c = y - als_baseline(y)
        noise = noise_level(x, c)
        peaks, _ = find_peaks(c, prominence=PROMINENCE_SIGMA * noise, distance=3)
        spectra[polymer] = (x, c, noise, peaks)
        strongest = c[peaks].max()

        say("=" * 72)
        say(f"{polymer}: {len(x)} points, {x.min():.1f}-{x.max():.1f} cm-1, "
            f"mean step {np.mean(np.diff(x)):.3f} cm-1, noise {noise:.1f} counts")
        say("=" * 72)
        say(f"{'shift':>8} {'rel. int.':>10} {'SNR':>8}   literature match")
        for i in peaks:
            tags = []
            for ref, bands in REFERENCE_BANDS.items():
                nearest = min(bands, key=lambda b: abs(b - x[i]))
                if abs(nearest - x[i]) <= MATCH_TOL:
                    tags.append(f"{ref} {nearest}")
            say(f"{x[i]:8.1f} {c[i] / strongest:10.3f} {c[i] / noise:8.1f}   {', '.join(tags)}")

        own = REFERENCE_BANDS[polymer]
        found = [b for b in own if np.any(np.abs(x[peaks] - b) <= MATCH_TOL)]
        missing = [b for b in own if b not in found]
        unassigned = [round(x[i]) for i in peaks
                      if np.all(np.abs(np.array(own) - x[i]) > MATCH_TOL)]
        say(f"Reference {polymer} bands found: {len(found)}/{len(own)}   missing: {missing}")
        say(f"Peaks not assigned to {polymer}: {unassigned}")
        say()

    say("=" * 72)
    say("Marker bands of other materials")
    say("  peak detected = a 6-sigma peak within +/- 4 cm-1 of the marker;")
    say("  intensity = baseline-corrected signal there, % of the strongest band")
    say("=" * 72)
    for polymer, (x, c, noise, peaks) in spectra.items():
        top = c.max()
        say(f"{polymer}:")
        for name, centre in MARKERS.items():
            if name.split()[0] == polymer:
                continue  # a polymer's own bands are not contamination markers
            v = marker_signal(x, c, centre)
            hit = np.any(np.abs(x[peaks] - centre) <= 4.0)
            say(f"  {name:<13} peak detected: {'YES' if hit else 'no ':<4}"
                f"  intensity {100 * v / top:5.1f} %")
    say()

    if args.out:
        with open(args.out, "w") as fh:
            fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
