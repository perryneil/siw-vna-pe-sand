# SIW–VNA discrimination of polyethylene microplastic from sand

Data and analysis code for the article:

> **Higher-order microwave resonances discriminate polyethylene microplastic from sand**
> Perry Neil J. Fernandez, Romeric F. Pobre, Yuen-Wa Ho
> Submitted to the *South African Journal of Science* (2026).

A two-port substrate integrated waveguide (SIW) cavity was measured with an open-source
LibreVNA vector network analyser. The measurements are an air baseline plus five independent
fill-and-reload trials each of polyethylene (PE) microplastic and beach sand. Each sweep has
501 points from 1 MHz to 6 GHz. The code in this repository recomputes the statistics and
characterisation results in the article, and redraws the data figures, directly from the raw
files.

## Repository layout

```
data/
  vna/        11 raw LibreVNA exports (air baseline, PE_0001–0005, sand_0001–0005)
  raman/      Raman spectrum of the PE fraction (Php-MP-Raman-PE.txt)
  simulation/ CST simulation of the air-filled cavity (S21-air-cst-v2.txt)
analysis/
  siw_analysis.py          all transmission statistics (Tables 2 and 3, Results)
  expected_output.txt      the output of that script, for comparison
  siw_analysis.R           independent R version of siw_analysis.py
  expected_output_R.txt    the output of the R version
  compare_r_python.py      runs both versions and compares their results
  R_vs_Python_comparison.md   outcome of that comparison
characterisation/
  raman_band_analysis.py   Raman band assignment for the PE fraction
  expected_output.txt      the output of that script, for comparison
figures/
  common.py                shared loading helpers
  make_figure3.py          Figure 3
  make_figure4.py          Figure 4
  make_figure5.py          Figure 5
  make_supplementary_figure1.py   Supplementary figure 1
  output/                  figures as produced by the scripts
requirements.txt
LICENSE
```

## Setup

The code was run with Python 3.11, NumPy 2.4.4, SciPy 1.17.1 and Matplotlib 3.10.9.

```bash
pip install -r requirements.txt
```

## Running

Run every command from the repository root.

```bash
python analysis/siw_analysis.py                     # transmission statistics
python characterisation/raman_band_analysis.py      # Raman band analysis
python figures/make_figure3.py                      # -> figures/output/Figure_3.png
python figures/make_figure4.py                      # -> figures/output/Figure_4.png
python figures/make_figure5.py                      # -> figures/output/Figure_5.png
python figures/make_supplementary_figure1.py        # -> figures/output/Supplementary_figure_1.png
```

The two analysis scripts accept `--out FILE` to save a transcript. They also accept `--data DIR`
if the data are stored elsewhere. Compare your transcript with the `expected_output.txt` in the
same folder.

## R version of the statistics

`analysis/siw_analysis.R` reimplements `analysis/siw_analysis.py` in base R (R 4.0 or later, no
add-on packages). It uses the same feature frequencies, search bands, notch convention and tests,
and prints the same report. Run it from the repository root:

```bash
Rscript analysis/siw_analysis.R                     # same options: --data DIR, --out FILE
python analysis/compare_r_python.py                 # run both versions and compare
```

If `Rscript` is not on the PATH, pass its location to the comparison script, for example
`--rscript "C:/Program Files/R/R-4.4.1/bin/Rscript.exe"`.

The R report is identical, line for line, to `analysis/expected_output.txt`. Compared at full
double precision, the largest relative difference between the two versions over every
p-value, effect size, notch parameter and summary statistic is 3.6 × 10⁻¹². This comes from
floating-point rounding. Two implementation details differ:

- Welch's t-test is written out explicitly instead of calling `t.test()`, so that degenerate
  frequency points give p = 1, as in the Python script, instead of an error. On these data it
  agrees with `t.test(var.equal = FALSE)` at all 501 points.
- The Holm and Benjamini–Hochberg adjustments use `p.adjust()`, which applies the same
  formulas as the hand-written Python functions.

Details are in `analysis/R_vs_Python_comparison.md`.

## What each script reproduces

**`analysis/siw_analysis.py`** reads the eleven VNA exports and prints, in order:

1. Point-wise Welch t-tests at all 501 frequencies, with Holm–Bonferroni (family-wise error
   rate) and Benjamini–Hochberg (false discovery rate) adjustment. This gives the seven
   Holm-surviving points: 1.0088 GHz, 3.4204–3.4804 GHz and 5.5081 GHz.
2. **Table 2**: mean ± SD for sand and PE at each analysed feature, the difference, Cohen's d,
   raw and adjusted p. It also gives the F3 per-trial dip-depth test.
3. **Table 3**: notch frequency, depth and loaded quality factor Q_L for air, PE and sand.
4. The one-sample-step displacement of the loaded minimum at F2 and F4.
5. The air-baseline flank slope |dS21/df| at F1–F4.
6. The F4 diagnostic: the −67.2 dB single sample, its neighbours, the phase reversal across
   the null, and the 5–6 GHz band median and tenth percentile.
7. Class ranges and gaps, the leave-one-out nearest-centroid result, and the permutation
   p-value for complete separation, 1/C(10,5).

**`characterisation/raman_band_analysis.py`** processes the PE spectrum in five steps:

1. Removes the background by asymmetric least squares.
2. Estimates the noise over 1680–1760 cm⁻¹.
3. Detects peaks with a prominence of at least six times the noise.
4. Assigns the peaks to literature PE, PP and PS bands (±8 cm⁻¹).
5. Checks the PP, PS and calcite marker positions for any detected peak.

The PE spectrum matches 8 of 9 reference PE bands, including the 1418 cm⁻¹ band of the
orthorhombic crystalline phase. No peak is detected at any PP, PS or calcite marker. The
intensity near 1086 cm⁻¹ is the amorphous ν(C–C) shoulder of PE itself, not a calcite band.
Because a single spectrum was measured, this confirms the identity of the probed material,
not the purity of the whole batch.

The particle-size and shape values reported in the article (Dv10, Dv50, Dv90 and circularity)
are taken directly from the instrument report of the automated image analysis. They are not
computed by code in this repository.

## Data format

Each file in `data/vna/` is a LibreVNA CSV export. It has a header row and 501 rows at 11.998 MHz
spacing from 1 MHz to 6 GHz. The columns include `Frequency` (Hz),
`S21_Magnitude (linear)`, `S21_Phase` (degrees), `S21_Real` and `S21_Imaginary`, and the
matching S11 columns. The scripts convert magnitude to dB as 20·log10|S21|. The file names
are kept exactly as exported, because the scripts identify the PE and sand trials by the
`_PE_` and `_sand_` patterns.

`data/raman/Php-MP-Raman-PE.txt` has two tab-separated columns: Raman shift (cm⁻¹) and
intensity (counts). The spectrum covers 661–1760 cm⁻¹. It was acquired with 785 nm excitation
through a 10× objective.

`data/simulation/S21-air-cst-v2.txt` is the CST Studio Suite export of the simulated air-filled
cavity: two columns, frequency (GHz) and |S21| (dB), from 0.1 to 6 GHz. Figure 3 plots it and the
measured air baseline each on its own frequency axis, without resampling or alignment.

## Conventions made explicit by the code

- **Notch shoulders and bandwidth.** From each notch minimum the code walks outward to the
  first local maximum on each side. The 3 dB level is set 3 dB below the mean of those two
  shoulders. Band edges are found by linear interpolation, so the bandwidth is not quantised to
  the sample grid.
- **Effect size.** Cohen's d uses the pooled sample standard deviation (ddof = 1).
- **Feature selection.** The two discriminating features were selected on the same ten
  measurements used for the leave-one-out result. That result therefore restates the observed
  separation; it is not an estimate of out-of-sample accuracy.

## Licence

The code is released under the MIT Licence (see `LICENSE`).

## Contact

Perry Neil J. Fernandez, University of the Philippines Visayas — pjfernandez@up.edu.ph
