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
analysis/
  siw_analysis.py          all transmission statistics (Tables 2 and 3, Results)
  expected_output.txt      the output of that script, for comparison
characterisation/
  raman_band_analysis.py   Raman band assignment for the PE fraction
  expected_output.txt      the output of that script, for comparison
figures/
  common.py                shared loading helpers
  make_figure4.py          Figure 4
  make_figure5.py          Figure 5
  make_supplementary_figure3.py   Supplementary figure 3
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
python figures/make_figure4.py                      # -> figures/output/Figure_4.png
python figures/make_figure5.py                      # -> figures/output/Figure_5.png
python figures/make_supplementary_figure3.py        # -> figures/output/Supplementary_figure_3.png
```

The two analysis scripts accept `--out FILE` to save a transcript. They also accept `--data DIR`
if the data are stored elsewhere. Compare your transcript with the `expected_output.txt` in the
same folder.

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
