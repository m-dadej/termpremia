# Survey-augmented dynamics against the Philadelphia Fed SPF
# =========================================================================
#
# Fits the US model with survey-disciplined P-dynamics and reports what the
# discipline buys. Lives in analysis/ rather than in the test suite or a
# vignette because it downloads data that cannot be redistributed: the
# Philadelphia Fed asserts copyright over the SPF with all rights reserved and
# limits use to research purposes. The package therefore fetches it and never
# ships it, and nothing that has to build on CRAN is allowed to depend on it.
#
# Results on the data as of September 2026 (763 month-ends to 2024-12):
#
#                       survey RMSE   rho     mean TP   sd TP   sd E[r] 10y
#     OLS                     -      0.9912    129.7    139.3       100
#     short horizons only    27.9    0.9953    135.4     83.7       160
#     long horizon only      42.1    0.9920    137.1    120.8        96
#     both                   28.7    0.9906    134.8    106.1       100
#
# THE FINDING. Short-horizon surveys alone make the shifting-endpoint problem
# WORSE, not better. Fitting forecasts that reach only four quarters pushes the
# VAR to be more persistent, and the expected short rate ten years out then
# swings more than it did under OLS. Anchoring the far end needs a far-end
# forecast: BILL10, asked once a year since 1992, is the only free one, and 35
# observations of it do more for the endpoint than 700 quarterly ones.
#
# Usage: Rscript analysis/validate-spf.R
# =========================================================================

N_FACTORS <- 5L
HORIZON   <- 120L          # months; the 10-year tenor
SD_GRID   <- c(0.001, 0.002, 0.003, 0.005, 0.01)
MAKE_PLOT <- TRUE

# --- load the package ----------------------------------------------------

find_pkg_root <- function() {
  candidates <- c(".", "..", "termpremia", file.path("..", "termpremia"))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    script_dir <- dirname(sub("^--file=", "", file_arg[1]))
    candidates <- c(file.path(script_dir, ".."), script_dir, candidates)
  }
  for (path in candidates) {
    desc <- file.path(path, "DESCRIPTION")
    if (file.exists(desc) &&
        any(grepl("^Package:\\s*termpremia", readLines(desc, warn = FALSE)))) {
      return(normalizePath(path))
    }
  }
  NULL
}

if (requireNamespace("termpremia", quietly = TRUE)) {
  library(termpremia)
} else {
  root <- find_pkg_root()
  if (is.null(root)) {
    stop("Could not locate the termpremia package source.", call. = FALSE)
  }
  pkgload::load_all(root, quiet = TRUE)
}

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("This script needs readxl: install.packages(\"readxl\").", call. = FALSE)
}

# --- fetch the surveys ---------------------------------------------------

cat("Downloading SPF forecasts (not redistributed; Philadelphia Fed terms)\n")
short <- spf_tbill()        # 1-4 quarters ahead, quarterly average of the 3m bill
long  <- spf_bill10()       # 10-year average of the 3m bill, Q1 surveys only

cat(sprintf("  TBILL : %4d forecasts, %s to %s\n", nrow(short),
            format(min(short$date)), format(max(short$date))))
cat(sprintf("  BILL10: %4d forecasts, %s to %s\n\n", nrow(long),
            format(min(long$date)), format(max(long$date))))

# Integrity check on the column indexing. TBILL1 is not a forecast: it is the
# previous quarter's realised rate as known to respondents. If our reading of
# the layout is right it must match the bundled bill series almost exactly, and
# if it does not, every horizon below is off by one and nothing downstream
# means anything.
check_layout <- function() {
  url <- paste0("https://www.philadelphiafed.org/-/media/frbp/assets/",
                "surveys-and-data/survey-of-professional-forecasters/",
                "data-files/files/mean_tbill_level.xlsx")
  dest <- tempfile(fileext = ".xlsx")
  on.exit(unlink(dest), add = TRUE)
  utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  raw <- as.data.frame(readxl::read_excel(dest))

  prev <- suppressWarnings(as.numeric(as.character(raw$TBILL1)))
  qstart <- as.Date(sprintf("%d-%02d-01", as.integer(raw$YEAR),
                            (as.integer(raw$QUARTER) - 1L) * 3L + 1L))
  prev_q <- as.Date(cut(qstart - 1L, "quarter"))

  tb <- tbill_3m
  tb$q <- as.Date(cut(tb$date, "quarter"))
  qavg <- aggregate(value ~ q, data = tb, FUN = mean)

  target <- qavg$value[match(prev_q, qavg$q)]
  ok <- is.finite(prev) & is.finite(target)
  sqrt(mean((prev[ok] - target[ok])^2)) * 100
}

layout_bp <- check_layout()
cat(sprintf("Layout check: TBILL1 vs realised previous quarter, %.1f bp RMSE\n",
            layout_bp))
if (layout_bp > 10) {
  stop("TBILL1 does not match the previous quarter's realised bill rate. ",
       "The SPF column layout has probably changed; every horizon below ",
       "would be wrong.", call. = FALSE)
}
cat("  (this is what pins the horizon indexing; anything above ~10bp is a bug)\n\n")

# --- fit -----------------------------------------------------------------

panel <- yield_panel(gsw_monthly, units = "percent", maturity_unit = "months",
                     instrument = "government", issuer = "US")

fits <- list(ols = suppressWarnings(atsm(panel, n_factors = N_FACTORS)))
specs <- list(short = short, long = long, both = rbind(short, long))

for (nm in names(specs)) {
  fits[[nm]] <- suppressWarnings(suppressMessages(
    atsm(panel, n_factors = N_FACTORS, p_dynamics = "survey",
         survey = specs[[nm]])
  ))
}

print(fits$both)

# --- the diagnostic that matters -----------------------------------------

# Expected short rate h months ahead, at every date. The shifting-endpoint
# problem is this quantity swinging around at long h.
expected_short_rate_ahead <- function(f, h) {
  ex <- f$factors
  for (p in seq_len(h)) {
    ex <- ex %*% t(f$pars$phi) + rep(f$pars$mu, each = nrow(ex))
  }
  (f$pars$delta0 + drop(ex %*% f$pars$delta1)) * 12 * 1e4
}

j <- match(HORIZON, fits$ols$maturities)

cat(sprintf("\n=== %d-year decomposition and endpoint behaviour, bp ===\n",
            HORIZON %/% 12))
cat(sprintf("%-8s %9s %8s %8s %8s %11s %11s\n",
            "spec", "svy RMSE", "rho", "mean TP", "sd TP",
            "sd E[r] 5y", "sd E[r] 10y"))
for (nm in names(fits)) {
  f <- fits[[nm]]
  rmse <- if (is.null(f$survey)) NA_real_ else f$survey$rmse_bp[["survey"]]
  cat(sprintf("%-8s %9.1f %8.4f %8.1f %8.1f %11.0f %11.0f\n", nm, rmse,
              f$spectral_radius[["real_world"]],
              mean(f$term_premium[, j]) * 1e4,
              stats::sd(f$term_premium[, j]) * 1e4,
              stats::sd(expected_short_rate_ahead(f, 60L)),
              stats::sd(expected_short_rate_ahead(f, 120L))))
}

cat("\nRead the last column, not the first. Fitting short-horizon surveys\n")
cat("better is not the same as anchoring the far end, and on its own it\n")
cat("moves the far end the wrong way.\n")

# --- how well is the long-horizon anchor itself matched? -----------------

cat("\n=== fit to BILL10 specifically ===\n")
for (nm in names(fits)) {
  f <- fits[[nm]]
  tg <- suppressMessages(survey_targets(long, panel$dates, survey_control()))
  d <- survey_design(tg)
  implied <- survey_implied(f$pars$mu, f$pars$phi,
                            f$factors[d$rows, , drop = FALSE],
                            f$pars$delta0, f$pars$delta1, d)
  cat(sprintf("%-8s RMSE %5.0f bp   model mean %5.0f vs survey %5.0f bp\n", nm,
              sqrt(mean((tg$value - implied)^2)) * 12 * 1e4,
              mean(implied) * 12 * 1e4, mean(tg$value) * 12 * 1e4))
}

# --- sensitivity to how much the surveys are trusted ---------------------

cat("\n=== sensitivity to survey_control(sd) ===\n")
cat(sprintf("%-8s %9s %8s %8s %11s\n", "sd (bp)", "svy RMSE", "rho", "sd TP",
            "sd E[r] 10y"))
for (s in SD_GRID) {
  f <- suppressWarnings(suppressMessages(
    atsm(panel, n_factors = N_FACTORS, p_dynamics = "survey",
         survey = specs$both, survey_control = survey_control(sd = s))
  ))
  cat(sprintf("%-8.0f %9.1f %8.4f %8.1f %11.0f\n", s * 1e4,
              f$survey$rmse_bp[["survey"]], f$spectral_radius[["real_world"]],
              stats::sd(f$term_premium[, j]) * 1e4,
              stats::sd(expected_short_rate_ahead(f, 120L))))
}

cat("\nThe model cannot get below about 27bp of disagreement with the SPF\n")
cat("however hard it is pushed, which is why the default sd is 30bp: it is\n")
cat("roughly self-consistent with the disagreement actually achieved.\n")

# --- plot ----------------------------------------------------------------

if (MAKE_PLOT) {
  op <- graphics::par(mfrow = c(2, 1), mar = c(3, 4, 2, 1))

  d <- fits$ols$dates
  tp_ols <- fits$ols$term_premium[, j] * 1e4
  tp_svy <- fits$both$term_premium[, j] * 1e4

  plot(d, tp_ols, type = "l", lwd = 1.4, col = "firebrick", xlab = "",
       ylab = "bp", main = sprintf("US %dy term premium", HORIZON %/% 12),
       ylim = range(tp_ols, tp_svy))
  graphics::lines(d, tp_svy, lwd = 1.4, col = "darkgreen")
  graphics::abline(h = 0, col = "grey85", lty = 2)
  graphics::legend("topright", c("OLS dynamics", "survey-disciplined"),
                   col = c("firebrick", "darkgreen"), lwd = 1.4, bty = "n",
                   cex = 0.8)

  e_ols <- expected_short_rate_ahead(fits$ols, 120L)
  e_short <- expected_short_rate_ahead(fits$short, 120L)
  e_both <- expected_short_rate_ahead(fits$both, 120L)

  plot(d, e_ols, type = "l", lwd = 1.4, col = "firebrick", xlab = "",
       ylab = "bp", main = "Expected short rate 10 years ahead",
       ylim = range(e_ols, e_short, e_both))
  graphics::lines(d, e_short, lwd = 1.4, col = "grey55")
  graphics::lines(d, e_both, lwd = 1.4, col = "darkgreen")
  graphics::legend("topright",
                   c("OLS", "short-horizon surveys only", "short + BILL10"),
                   col = c("firebrick", "grey55", "darkgreen"), lwd = 1.4,
                   bty = "n", cex = 0.8)

  graphics::par(op)
}

invisible(fits)
