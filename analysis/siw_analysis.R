#!/usr/bin/env Rscript
# =============================================================================
# siw_analysis.R
# =============================================================================
#
# R port of siw_analysis.py. Reproduces every statistic reported in
#
#     "Higher-order microwave resonances discriminate polyethylene
#      microplastic from sand"
#
# directly from the eleven raw LibreVNA CSV exports, using the same
# feature frequencies, search bands, notch convention and tests as the
# Python script. The printed report is laid out line-for-line like the
# Python transcript (expected_output.txt) so the two can be diffed.
#
# Inputs
# ------
# A directory containing the eleven exports (default: ../data/vna, relative
# to this script):
#
#     ses06062026_bat0001_s11s21_empty.csv     air baseline
#     ses06062026_s11s21_PE_0001.csv ... _0005.csv
#     ses06062026_s11s21_sand_0001.csv ... _0005.csv
#
# Usage
# -----
#     Rscript analysis/siw_analysis.R                     # from the repository root
#     Rscript analysis/siw_analysis.R --out results.txt   # also save a transcript
#     Rscript analysis/siw_analysis.R --data /path/to/csvs
#     Rscript analysis/siw_analysis.R --dump dump_dir     # full-precision CSVs
#                                                         # (used for the Python
#                                                         #  comparison)
#
# Requires: base R only (>= 4.0). No add-on packages.
#
# Notes on the translation
# ------------------------
# * Welch's t-test is computed explicitly (welch_test below) rather than with
#   stats::t.test(): t.test() stops with "data are essentially constant" on
#   near-identical groups, whereas SciPy returns NaN, which the Python script
#   then sets to p = 1. welch_test() does the same. On these data it agrees
#   with t.test(var.equal = FALSE) at all 501 points to within 4e-16.
# * Holm and Benjamini-Hochberg adjustments use stats::p.adjust(), which is
#   algebraically identical to the hand-written holm_bonferroni() and
#   benjamini_hochberg() in the Python script.
# * numpy.gradient (second-order central differences, first-order at the
#   ends) and numpy.percentile (linear interpolation = R quantile type 7)
#   are reproduced exactly.
# * R indexes from 1; every index in the notch walk is shifted accordingly.
#
# Author: P.N.J. Fernandez
# =============================================================================

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Fixed analysis frequencies (Hz). These are exact sample points of the sweep.
FEATURES <- c(
  "F1 min"    = 2.3646e9,
  "F2 flank1" = 3.4204e9,
  "F2 flank2" = 3.4444e9,
  "F2 flank3" = 3.4684e9,
  "F2 min"    = 3.5044e9,
  "F3 dip"    = 4.6322e9,
  "F4 null"   = 5.5081e9
)

# Search bands used to locate each notch minimum per trial (Hz).
BANDS <- list(
  F1 = c(2.20e9, 2.50e9),
  F2 = c(3.40e9, 3.60e9),
  F3 = c(4.50e9, 4.75e9),
  F4 = c(5.40e9, 5.62e9)
)

ALPHA <- 0.05

script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f) == 1) dirname(normalizePath(f)) else getwd()
}
DEFAULT_DATA <- file.path(script_dir(), "..", "data", "vna")

# -----------------------------------------------------------------------------
# I/O
# -----------------------------------------------------------------------------

load_export <- function(path) {
  # Return list(f, s21_db, phase, re, im). Columns are found by header name,
  # because the air and loaded exports have different column sets.
  d <- utils::read.csv(path, check.names = FALSE, strip.white = TRUE)
  d <- d[, names(d) != "", drop = FALSE]          # trailing comma -> empty column
  mag <- d[["S21_Magnitude (linear)"]]
  list(
    f     = as.numeric(d[["Frequency"]]),
    db    = 20 * log10(as.numeric(mag)),
    phase = as.numeric(d[["S21_Phase"]]),
    re    = as.numeric(d[["S21_Real"]]),
    im    = as.numeric(d[["S21_Imaginary"]])
  )
}

load_all <- function(data_dir) {
  air_path   <- file.path(data_dir, "ses06062026_bat0001_s11s21_empty.csv")
  pe_paths   <- sort(Sys.glob(file.path(data_dir, "*_PE_*.csv")))
  sand_paths <- sort(Sys.glob(file.path(data_dir, "*_sand_*.csv")))

  if (!file.exists(air_path)) stop("Air baseline not found: ", air_path)
  if (length(pe_paths) != 5 || length(sand_paths) != 5)
    stop(sprintf("Expected 5 PE and 5 sand files, found %d and %d in %s",
                 length(pe_paths), length(sand_paths), data_dir))

  air  <- load_export(air_path)
  pe_l <- lapply(pe_paths, load_export)
  sa_l <- lapply(sand_paths, load_export)

  # every export must share the same frequency grid
  for (x in c(pe_l, sa_l))
    if (!isTRUE(all.equal(x$f, air$f))) stop("grid mismatch")

  list(freq = air$f, air = air$db, air_phase = air$phase,
       air_re = air$re, air_im = air$im,
       pe   = do.call(rbind, lapply(pe_l, `[[`, "db")),   # 5 x 501
       sand = do.call(rbind, lapply(sa_l, `[[`, "db")),
       pe_files = basename(pe_paths), sand_files = basename(sand_paths))
}

# -----------------------------------------------------------------------------
# Statistics
# -----------------------------------------------------------------------------

welch_test <- function(a, b) {
  # Two-sided Welch (unequal-variance) t-test; NaN when undefined, as SciPy.
  na <- length(a); nb <- length(b)
  va <- stats::var(a) / na; vb <- stats::var(b) / nb
  t  <- (mean(a) - mean(b)) / sqrt(va + vb)
  df <- (va + vb)^2 / (va^2 / (na - 1) + vb^2 / (nb - 1))
  p  <- suppressWarnings(2 * stats::pt(-abs(t), df))
  list(statistic = t, df = df, pvalue = p)
}

cohens_d <- function(a, b) {
  # Standardised effect size with pooled SD (sample SD, n - 1).
  na <- length(a); nb <- length(b)
  sp <- sqrt(((na - 1) * stats::var(a) + (nb - 1) * stats::var(b)) / (na + nb - 2))
  c(d = abs(mean(a) - mean(b)) / sp, sp = sp)
}

np_gradient <- function(y, x) {
  # numpy.gradient(y, x): 2nd-order central differences (non-uniform
  # spacing allowed) in the interior, 1st-order one-sided at the ends.
  n <- length(y); g <- numeric(n)
  hs <- x[3:n] - x[2:(n - 1)]          # forward step
  hd <- x[2:(n - 1)] - x[1:(n - 2)]    # backward step
  g[2:(n - 1)] <- (hd^2 * y[3:n] + (hs^2 - hd^2) * y[2:(n - 1)] - hs^2 * y[1:(n - 2)]) /
                  (hs * hd * (hd + hs))
  g[1] <- (y[2] - y[1]) / (x[2] - x[1])
  g[n] <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])
  g
}

# -----------------------------------------------------------------------------
# Notch geometry
# -----------------------------------------------------------------------------

notch_parameters <- function(freq, trace_db, band) {
  # Locate the notch minimum inside `band` and characterise it.
  # The 3 dB reference level is 3 dB below the mean of the two shoulder
  # maxima bounding the notch; band edges are linearly interpolated, so the
  # width is not quantised to the sample grid.
  n   <- length(trace_db)
  idx <- which(freq >= band[1] & freq <= band[2])
  imin <- idx[which.min(trace_db[idx])]

  # walk outwards to the first local maximum on each side
  i <- imin
  while (i - 1 >= 1 && trace_db[i - 1] > trace_db[i]) i <- i - 1
  left_shoulder <- i
  i <- imin
  while (i + 1 <= n && trace_db[i + 1] > trace_db[i]) i <- i + 1
  right_shoulder <- i

  shoulder_mean <- 0.5 * (trace_db[left_shoulder] + trace_db[right_shoulder])
  depth <- shoulder_mean - trace_db[imin]
  level <- shoulder_mean - 3.0

  crossing <- function(start, step) {
    j <- start
    while (j + step >= 1 && j + step <= n && trace_db[j] < level) j <- j + step
    j0 <- j - step; j1 <- j
    y0 <- trace_db[j0]; y1 <- trace_db[j1]
    if (y1 == y0) return(freq[j1])
    t <- (level - y0) / (y1 - y0)
    freq[j0] + t * (freq[j1] - freq[j0])
  }

  if (depth <= 3.0) {
    bw <- NA_real_; q <- NA_real_; npts <- 0L
  } else {
    f_lo <- crossing(imin, -1)
    f_hi <- crossing(imin, +1)
    bw <- f_hi - f_lo
    q  <- if (bw > 0) freq[imin] / bw else NA_real_
    npts <- as.integer(sum(trace_db[left_shoulder:right_shoulder] < level))
  }
  list(fr = freq[imin], idx = imin, depth = depth, bw = bw,
       q_loaded = q, resolved_points = npts, shoulder_mean = shoulder_mean)
}

# -----------------------------------------------------------------------------
# Formatting helpers (mirror Python's f-strings / numpy printing)
# -----------------------------------------------------------------------------

fmt_p <- function(p) if (p >= 0.001) sprintf("%.3f", p) else sprintf("%.2e", p)

py_list <- function(x) paste0("[", paste0("'", x, "'", collapse = ", "), "]")

np_array <- function(x) {
  # numpy default float-array printing: shortest digits (max 8) shared by all
  digs <- vapply(x, function(v) {
    s <- sub("0+$", "", sprintf("%.8f", v)); nchar(sub("^[^.]*\\.", "", s))
  }, integer(1))
  paste0("[", paste(sprintf(paste0("%.", max(digs), "f"), x), collapse = " "), "]")
}

nan_sd <- function(x) stats::sd(x[!is.na(x)])

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

parse_args <- function(argv) {
  opt <- list(data = DEFAULT_DATA, out = NULL, dump = NULL)
  i <- 1
  while (i <= length(argv)) {
    key <- sub("^--", "", argv[i])
    if (!key %in% names(opt) || i == length(argv))
      stop("usage: Rscript siw_analysis.R [--data DIR] [--out FILE] [--dump DIR]")
    opt[[key]] <- argv[i + 1]; i <- i + 2
  }
  opt
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  D <- load_all(args$data)
  freq <- D$freq; air <- D$air; pe <- D$pe; sand <- D$sand
  nf <- length(freq)

  lines <- character(0)
  say <- function(s = "") { lines <<- c(lines, s); cat(s, "\n", sep = "") }
  rule <- function(ch) strrep(ch, 78)

  step <- stats::median(diff(freq))
  say(rule("="))
  say("SIW-VNA microplastic discrimination: full statistical reproduction")
  say(rule("="))
  say(sprintf("Points per sweep      : %d", nf))
  say(sprintf("Span                  : %.3f MHz to %.3f GHz", freq[1] / 1e6, freq[nf] / 1e9))
  say(sprintf("Median point spacing  : %.3f MHz", step / 1e6))
  say(sprintf("PE trials             : %d  %s", nrow(pe), py_list(D$pe_files)))
  say(sprintf("Sand trials           : %d  %s", nrow(sand), py_list(D$sand_files)))
  say()

  # ---------------------------------------------------------------------------
  # 1. Point-by-point Welch tests over the whole sweep, with corrections
  # ---------------------------------------------------------------------------
  say(rule("-"))
  say("1. Point-wise Welch t-tests across all frequency points")
  say(rule("-"))

  raw_p <- vapply(seq_len(nf), function(k) welch_test(sand[, k], pe[, k])$pvalue, numeric(1))
  raw_p[is.nan(raw_p) | is.na(raw_p)] <- 1.0

  holm <- stats::p.adjust(raw_p, method = "holm")
  bh   <- stats::p.adjust(raw_p, method = "BH")

  say(sprintf("Raw p < 0.05                 : %d points", sum(raw_p < ALPHA)))
  say(sprintf("Holm-Bonferroni p_adj < 0.05 : %d points", sum(holm < ALPHA)))
  say(sprintf("Benjamini-Hochberg q < 0.05  : %d points", sum(bh < ALPHA)))
  say(sprintf("Benjamini-Hochberg q < 0.01  : %d points", sum(bh < 0.01)))
  say()
  say("Holm survivors (p_adj < 0.05):")
  say(sprintf("  %12s %13s %6s %10s %8s %8s",
              "Freq (GHz)", "sand-PE (dB)", "d", "raw p", "Holm", "BH q"))
  for (k in which(holm < ALPHA)) {
    d <- cohens_d(sand[, k], pe[, k])[["d"]]
    say(sprintf("  %12.4f %13.3f %6.2f %10s %8.3f %8.4f",
                freq[k] / 1e9, mean(sand[, k]) - mean(pe[, k]), d,
                fmt_p(raw_p[k]), holm[k], bh[k]))
  }
  say()

  # ---------------------------------------------------------------------------
  # 2. Table 2: replicate statistics at the analysed features
  # ---------------------------------------------------------------------------
  say(rule("-"))
  say("2. Table 2 - replicate statistics by spectral feature")
  say(rule("-"))
  say(sprintf("%-11s%9s%9s%19s%19s%8s%6s%10s%9s%9s",
              "Feature", "GHz", "air", "sand mean+-SD", "PE mean+-SD",
              "diff", "d", "raw p", "Holm", "BH q"))
  table2 <- NULL
  for (name in names(FEATURES)) {
    k <- which.min(abs(freq - FEATURES[[name]]))
    s <- sand[, k]; p_ <- pe[, k]
    d <- cohens_d(s, p_)[["d"]]
    df_ <- mean(s) - mean(p_)
    say(sprintf("%-11s%9.4f%9.2f%11.2f +-%5.2f%11.2f +-%5.2f%+8.2f%6.2f%10s%9.3f%9.4f",
                name, freq[k] / 1e9, air[k],
                mean(s), stats::sd(s), mean(p_), stats::sd(p_),
                df_, d, fmt_p(raw_p[k]), holm[k], bh[k]))
    table2 <- rbind(table2, data.frame(
      feature = name, freq = freq[k], air = air[k],
      sand_mean = mean(s), sand_sd = stats::sd(s),
      pe_mean = mean(p_), pe_sd = stats::sd(p_),
      diff = df_, d = d, raw_p = raw_p[k], holm = holm[k], bh = bh[k]))
  }
  say()
  say("Note: the F3 row above is the point-wise test at 4.6322 GHz. The paper")
  say("      also reports a per-trial dip-depth test for F3, computed below.")
  say()

  # F3 dip-depth test (a single derived test, outside the 501-point family)
  f3 <- freq >= BANDS$F3[1] & freq <= BANDS$F3[2]
  pe_min   <- apply(pe,   1, function(t) min(t[f3]))
  sand_min <- apply(sand, 1, function(t) min(t[f3]))
  t_dip <- welch_test(sand_min, pe_min)
  say(sprintf("F3 per-trial minimum depth: sand %.2f +- %.2f dB, PE %.2f +- %.2f dB, Welch p = %.3f",
              mean(sand_min), stats::sd(sand_min), mean(pe_min), stats::sd(pe_min),
              t_dip$pvalue))
  say()

  # ---------------------------------------------------------------------------
  # 3. Table 3: loaded quality factors and notch parameters
  # ---------------------------------------------------------------------------
  say(rule("-"))
  say("3. Table 3 - loaded quality factors and notch parameters")
  say(rule("-"))
  say(sprintf("%-8s%-7s%12s%14s%14s%10s", "Feature", "Load", "fr (GHz)",
              "depth (dB)", "Q_L", "pts <3dB"))
  notch_tab <- NULL
  add_notch <- function(fname, load, trial, r)
    notch_tab <<- rbind(notch_tab, data.frame(
      feature = fname, load = load, trial = trial, fr = r$fr, depth = r$depth,
      bw = r$bw, q_loaded = r$q_loaded, resolved_points = r$resolved_points))
  for (fname in names(BANDS)) {
    band <- BANDS[[fname]]
    a <- notch_parameters(freq, air, band)
    add_notch(fname, "air", 0, a)
    say(sprintf("%-8s%-7s%12.4f%14.2f%14.1f%10d", fname, "air",
                a$fr / 1e9, a$depth, a$q_loaded, a$resolved_points))
    for (label in c("PE", "sand")) {
      group <- if (label == "PE") pe else sand
      res <- lapply(seq_len(nrow(group)), function(i) notch_parameters(freq, group[i, ], band))
      for (i in seq_along(res)) add_notch(fname, label, i, res[[i]])
      fr  <- vapply(res, `[[`, numeric(1), "fr") / 1e9
      dep <- vapply(res, `[[`, numeric(1), "depth")
      q   <- vapply(res, `[[`, numeric(1), "q_loaded")
      say(sprintf("%-8s%-7s%8.4f +-%4.4f%9.2f +-%5.2f%9.1f +-%5.1f%10s",
                  "", label, mean(fr), stats::sd(fr), mean(dep), stats::sd(dep),
                  mean(q, na.rm = TRUE), nan_sd(q), ""))
    }
  }
  say()

  # ---------------------------------------------------------------------------
  # 4. Displacement of the loaded minimum
  # ---------------------------------------------------------------------------
  say(rule("-"))
  say("4. Displacement of the loaded transmission minimum")
  say(rule("-"))
  shift_tab <- NULL
  for (fname in c("F2", "F4")) {
    band <- BANDS[[fname]]
    fr_pe   <- apply(pe,   1, function(t) notch_parameters(freq, t, band)$fr)
    fr_sand <- apply(sand, 1, function(t) notch_parameters(freq, t, band)$fr)
    shift <- mean(fr_pe) - fr_sand            # sand below PE is positive
    t_shift <- welch_test(fr_sand, fr_pe)
    disp <- mean(fr_pe) - mean(fr_sand)
    say(sprintf("%s: PE minima  %s GHz", fname, np_array(sort(unique(fr_pe / 1e9)))))
    say(sprintf("%s: sand minima %s GHz", fname, np_array(sort(unique(fr_sand / 1e9)))))
    say(sprintf("%s: mean sand-below-PE displacement = %.3e Hz (%.1f +- %.1f MHz), Welch p = %.4f",
                fname, disp, disp / 1e6, stats::sd(shift) / 1e6, t_shift$pvalue))
    say(sprintf("%s: one sweep step = %.3f MHz -- displacement is %.2f steps",
                fname, step / 1e6, disp / step))
    shift_tab <- rbind(shift_tab, data.frame(
      feature = fname, disp_hz = disp, sd_hz = stats::sd(shift), p = t_shift$pvalue))
  }
  say()

  # ---------------------------------------------------------------------------
  # 5. Flank steepness on the air baseline
  # ---------------------------------------------------------------------------
  say(rule("-"))
  say("5. Flank steepness |dS21/df| on the air baseline")
  say(rule("-"))
  grad <- np_gradient(air, freq / 1e6)        # dB per MHz
  flanks <- c(F1 = 2.341e9, F2 = 3.468e9, F3 = 4.596e9, F4 = 5.496e9)
  grad_tab <- NULL
  for (label in names(flanks)) {
    k <- which.min(abs(freq - flanks[[label]]))
    say(sprintf("%s flank at %.4f GHz: |dS21/df| = %.3f dB/MHz", label, freq[k] / 1e9, abs(grad[k])))
    grad_tab <- rbind(grad_tab, data.frame(feature = label, freq = freq[k], grad = abs(grad[k])))
  }
  say()

  # ---------------------------------------------------------------------------
  # 6. The F4 transmission zero: phase diagnostic
  # ---------------------------------------------------------------------------
  say(rule("-"))
  say("6. F4 diagnostic - is the deep air null a transmission zero?")
  say(rule("-"))
  k <- which.min(abs(freq - 5.5081e9))
  say(sprintf("Air |S21| at %.4f GHz : %.2f dB", freq[k] / 1e9, air[k]))
  say(sprintf("  neighbour below (%.4f GHz): %.2f dB", freq[k - 1] / 1e9, air[k - 1]))
  say(sprintf("  neighbour above (%.4f GHz): %.2f dB", freq[k + 1] / 1e9, air[k + 1]))
  ph <- D$air_phase
  say(sprintf("Phase: %.1f -> %.1f -> %.1f -> %.1f -> %.1f deg",
              ph[k - 2], ph[k - 1], ph[k], ph[k + 1], ph[k + 2]))
  say(sprintf("  reversal across the null = %.1f deg", abs(ph[k + 1] - ph[k - 1])))
  say(sprintf("Real/imag at the null: Re = %.3e, Im = %.3e", D$air_re[k], D$air_im[k]))
  b56 <- freq >= 5.0e9 & freq <= 6.0e9
  med56 <- stats::median(air[b56])
  p10   <- stats::quantile(air[b56], 0.10, type = 7, names = FALSE)
  say(sprintf("5-6 GHz band median = %.1f dB, 10th percentile = %.1f dB (no noise-floor plateau)",
              med56, p10))
  say()

  # ---------------------------------------------------------------------------
  # 7. Class separation and leave-one-out classification
  # ---------------------------------------------------------------------------
  say(rule("-"))
  say("7. Class separation, LOO nearest-centroid, permutation test")
  say(rule("-"))
  sep_tab <- NULL
  for (item in list(list("F2 flank 3.4684 GHz", 3.4684e9), list("F4 null 5.5081 GHz", 5.5081e9))) {
    k <- which.min(abs(freq - item[[2]]))
    s <- sand[, k]; p_ <- pe[, k]
    gap <- max(min(s), min(p_)) - min(max(s), max(p_))
    say(item[[1]])
    say(sprintf("  PE range   : %.2f to %.2f dB", min(p_), max(p_)))
    say(sprintf("  sand range : %.2f to %.2f dB", min(s), max(s)))
    say(sprintf("  gap between class ranges : %.2f dB (%s)", abs(gap),
                if (gap > 0) "separated" else "OVERLAP"))

    # leave-one-out nearest centroid on this single feature
    X <- c(s, p_)
    y <- c(rep(0, length(s)), rep(1, length(p_)))
    correct <- 0
    for (i in seq_along(X)) {
      mask <- rep(TRUE, length(X)); mask[i] <- FALSE
      c0 <- mean(X[mask & y == 0]); c1 <- mean(X[mask & y == 1])
      pred <- if (abs(X[i] - c0) < abs(X[i] - c1)) 0 else 1
      correct <- correct + as.integer(pred == y[i])
    }
    say(sprintf("  LOO nearest-centroid accuracy : %d/%d", correct, length(X)))
    sep_tab <- rbind(sep_tab, data.frame(feature = item[[1]], gap = abs(gap), loo = correct))
  }
  n_splits <- choose(10, 5)
  say(sprintf("Permutation p-value for complete separation with n = 5 per class: 1/C(10,5) = 1/%d = %.2e",
              as.integer(n_splits), 1 / n_splits))
  say()
  say(rule("="))
  say("Reminder: the two features above were selected using these same ten")
  say("measurements, so the leave-one-out figure is optimistically biased and")
  say("is a restatement of the observed separation, not an out-of-sample estimate.")
  say(rule("="))

  if (!is.null(args$out)) {
    writeLines(lines, args$out)
    cat(sprintf("\n[written to %s]\n", args$out))
  }

  # Optional full-precision dump for numerical comparison with the Python script
  if (!is.null(args$dump)) {
    dir.create(args$dump, showWarnings = FALSE, recursive = TRUE)
    w <- function(x, f) utils::write.csv(x, file.path(args$dump, f), row.names = FALSE)
    d_all <- vapply(seq_len(nf), function(k) cohens_d(sand[, k], pe[, k])[["d"]], numeric(1))
    w(data.frame(freq = freq, raw_p = raw_p, holm = holm, bh = bh,
                 diff = colMeans(sand) - colMeans(pe), d = d_all), "pointwise.csv")
    w(table2, "table2.csv")
    w(notch_tab, "notch.csv")
    w(shift_tab, "shift.csv")
    w(grad_tab, "gradient.csv")
    w(sep_tab, "separation.csv")
    w(data.frame(name = c("f3_dip_p", "band56_median", "band56_p10"),
                 value = c(t_dip$pvalue, med56, p10)), "scalars.csv")
  }
  invisible(NULL)
}

main()
