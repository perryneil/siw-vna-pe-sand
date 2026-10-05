# R port of `siw_analysis.py` — comparison with the Python original

Prepared 3 October 2026.

## Files

| File | Purpose |
|---|---|
| `siw_analysis.R` | Base-R port of `siw_analysis.py` (no add-on packages). The arguments are the same: `--data` and `--out`. `--dump DIR` also writes full-precision CSVs for checking. |
| `expected_output_R.txt` | Transcript of `Rscript analysis/siw_analysis.R`. |
| `compare_r_python.py` | Runs both scripts, diffs the reports and compares every underlying number at full double precision. |

To run it from the repository root:

```
Rscript analysis/siw_analysis.R
python analysis/compare_r_python.py          # needs Rscript on PATH, or --rscript "<path to Rscript.exe>"
```

## Environment used for the check

R 4.3.3 (base packages only); Python 3.13, NumPy 2.5.3, SciPy 1.18.1. The Python script
also reproduced `expected_output.txt` exactly under these newer versions (the pinned
versions are NumPy 2.4.4, SciPy 1.17.1). Line endings were not counted.

## Result 1 — printed reports

The R transcript is **identical line for line** (112 of 112 lines) to the Python output
and to `expected_output.txt`. This includes Table 2, Table 3, the F2/F4 displacement,
flank steepness, the F4 phase diagnostic, and the LOO and permutation results.

## Result 2 — full-precision values

| Quantity | n | max abs diff | max rel diff |
|---|---:|---:|---:|
| Point-wise Welch raw p | 501 | 1.2e-13 | 7.2e-13 |
| Holm-adjusted p | 501 | 1.8e-13 | 7.2e-13 |
| Benjamini–Hochberg q | 501 | 9.0e-14 | 6.2e-13 |
| Sand − PE mean difference | 501 | 1.4e-14 | 3.6e-12 |
| Cohen's d | 501 | 5.1e-13 | 3.6e-12 |
| Table 2 (air, means, SDs, diff, d) | 7 each | ≤ 5.8e-14 | ≤ 3.8e-14 |
| Notch f_r, resolved points | 44 | 0 | 0 |
| Notch depth (dB) / bandwidth (Hz) / Q_L | 44 | 5.0e-14 / 4.8e-7 / 7.1e-14 | ≤ 3.8e-15 |
| F2/F4 displacement (Hz) / SD (Hz) / Welch p | 2 each | 0 / 4.7e-9 / 4.5e-17 | ≤ 2.8e-15 |
| Flank gradient | 4 | 4.0e-15 | 2.8e-15 |
| F3 dip p, 5–6 GHz median, 10th percentile | 3 | 3.6e-14 | 4.4e-15 |

The largest relative difference is **3.6 × 10⁻¹²**. This is floating-point rounding: R and
NumPy add up means and variances in different orders. It is about nine orders of magnitude
below the last digit printed in the paper. The discrete results agree exactly: the notch
frequencies, the points below the 3 dB level, the count of significant points (273 / 7 / 183 / 52),
the Holm survivors and the LOO accuracies.

## How each test was translated

| Python | R | Equivalence |
|---|---|---|
| `scipy.stats.ttest_ind(..., equal_var=False)` | `welch_test()` (explicit Welch–Satterthwaite) | Agrees with `t.test(var.equal = FALSE)` at all 501 points (≤ 4e-16). An explicit function is used so that degenerate points give NaN → p = 1 as in SciPy, instead of `t.test()` stopping with an error. |
| `holm_bonferroni()` (hand-written) | `p.adjust(method = "holm")` | Same step-down formula; agrees to 7e-13. |
| `benjamini_hochberg()` (hand-written) | `p.adjust(method = "BH")` | Same step-up formula; agrees to 6e-13. |
| `cohens_d()` (pooled SD, ddof = 1) | `cohens_d()` using `var()` | Same estimator. |
| `np.std(ddof=1)`, `np.nanmean/nanstd` | `sd()`, `mean(na.rm = TRUE)` | Same estimator. |
| `notch_parameters()` | `notch_parameters()` | Line-for-line, with indices shifted to start at 1. Same shoulder walk, 3 dB level and interpolated crossings. |
| `np.gradient(y, x)` | `np_gradient()` | Same second-order non-uniform central difference, first-order at the ends. |
| `np.percentile(..., 10)` | `quantile(type = 7)` | Both use linear interpolation. |
| `np.argmin` (first minimum on ties) | `which.min` (first minimum on ties) | Same tie rule. |
| LOO nearest centroid, 1/C(10,5) | same loop, `choose(10, 5)` | Same. |

## Conclusion

The R implementation independently reproduces every statistic produced by `siw_analysis.py`,
and therefore every number in Tables 2 and 3 and in the text that the script supports.
