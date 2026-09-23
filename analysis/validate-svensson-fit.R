# Rebuilding a dense curve from a handful of tenors with svensson_fit()
# =========================================================================
#
# Samples a known monthly curve at eight standard tenors -- 3m, 6m, 1y, 2y,
# 3y, 5y, 7y, 10y -- rebuilds the monthly grid from them with svensson_fit(),
# and asks how far the ACM decomposition moves compared with fitting the
# original curve. Lives in analysis/ because the UK half reads a Bank of
# England archive.
#
# TWO CURVES, AND WHY THE SECOND ONE IS THE TEST. The US curve is GSW, which
# IS a Svensson curve, so rebuilding it with Svensson is flattering: the model
# can recover the truth exactly. The Bank of England curve is a variable
# roughness spline, not Svensson, so rebuilding it is a fair test of what
# happens to a user's own data. The UK sample is the one unbroken run in which
# the Bank publishes the curve all the way down to one month, 1997:03-2012:12
# (T = 190), so that the true short end is an observation, not the Bank's own
# extrapolation.
#
# The UK comparison uses FOUR factors because five are not estimable on the
# original curve itself over 190 months: the risk-adjusted dynamics come out
# explosive (spectral radius 1.12) before any rebuilding is involved, and
# atsm() says so. Four fit it to 1.9bp.
#
# Results as of September 2026. The 10-year term premium from the rebuilt
# curve, against the same model fitted to the original curve, with the true
# one-month yield supplied as `short_rate`:
#
#                          corr, changes   mean |diff|   max |diff|
#   US  Svensson                1.0000         0.1bp         1bp
#   US  Nelson-Siegel           0.9694         5.8bp        50bp
#   UK  Svensson                0.9987         4.4bp        13bp
#   UK  Nelson-Siegel           0.9954         5.0bp        15bp
#
#   curve error, bp           between tenors     one month (extrapolated)
#                              RMSE     max         RMSE     max
#   US  Svensson                0.0     0.7          1.2     20
#   US  Nelson-Siegel           2.1    20.8         15.2    259
#   UK  Svensson                0.3     8.9         11.9     57
#   UK  Nelson-Siegel           2.1    14.1          9.5     42
#
# WHAT THIS SETTLES.
#
# 1. Estimating the decays date by date is unstable, as feared -- the first
#    decay doubles or halves month to month in 30% of US months and 48% of
#    UK ones under Svensson -- but the instability does not reach the
#    premium. The fitted CURVE is well determined even where its parameters
#    are not, and the curve is all atsm() sees. No smoothing across dates is
#    needed.
#
# 2. Svensson is the better choice whenever there are tenors enough for it.
#    On the fair test it fits the curve seven times more closely between
#    tenors than Nelson-Siegel and lands nearer on the premium too, which is
#    why svensson_fit() defaults to it from seven tenors up.
#
# 3. The extrapolated one-month yield is poor in individual months, 20-260bp
#    off at worst. Supplying the observed short rate instead changed the
#    10-year premium's distance from the original by under 1bp on average in
#    both tests -- the long end is not where that error lands -- but the
#    one-month yield is what every excess return is measured against, so
#    atsm() warns when an extrapolated one is used.
#
# 4. The hold-out check is conservative. Dropping a tenor leaves the refit one
#    tenor short of the real fit, and on the UK curve the true error between
#    tenors was a fifth to three quarters of what the check reported. It
#    ranks the models correctly.
#
# 5. The decay search must start from several places. Refining only the best
#    grid cell, the first version of svensson_fit() missed basins lying
#    between grid points; refining the best three local minima halved the UK
#    premium error (mean 8.6bp to 4.4bp, max 26bp to 13bp) and cut the worst
#    UK one-month miss from 152bp to 57bp.
#
# Usage: Rscript analysis/validate-svensson-fit.R
# =========================================================================

pkgload::load_all(".", quiet = TRUE)

TENORS  <- c(3, 6, 12, 24, 36, 60, 84, 120)
BOE_ZIP <- file.path("data-raw", ".cache", "glcnominalmonthedata.zip")

bp <- function(x) round(x * 1e4, 1)

compare <- function(label, truth, rebuilt, n_factors, short_rate = NULL) {
  base <- suppressWarnings(atsm(truth, n_factors = n_factors,
                                short_rate = short_rate,
                                short_rate_units = "decimal"))
  fit <- suppressWarnings(atsm(rebuilt, n_factors = n_factors,
                               short_rate = short_rate,
                               short_rate_units = "decimal"))
  a <- base$term_premium[, "120"]
  b <- fit$term_premium[, "120"]
  data.frame(
    case = label,
    corr_levels = round(cor(a, b), 4),
    corr_changes = round(cor(diff(a), diff(b)), 4),
    mean_abs_bp = bp(mean(abs(a - b))),
    max_abs_bp = bp(max(abs(a - b))),
    rho_q = round(fit$spectral_radius[["risk_adjusted"]], 4)
  )
}

round_trip <- function(name, truth, n_factors) {
  y <- curve_matrix(truth)
  sparse <- yield_panel(y[, as.character(TENORS)], units = "decimal",
                        maturity_unit = "months")
  sr <- data.frame(date = truth$dates, value = y[, "1"])
  # Between the tenors only: 1m and 2m are extrapolation, reported apart.
  off <- setdiff(min(TENORS):max(TENORS), TENORS)

  cat("\n=== ", name, ": ", length(truth$dates), " months, ",
      format(min(truth$dates)), " to ", format(max(truth$dates)), " ===\n",
      sep = "")

  rows <- list()
  for (model in c("svensson", "nelson_siegel")) {
    t0 <- Sys.time()
    fit <- svensson_fit(sparse, yield_type = "zero", model = model)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    dense <- predict(fit, maturities = 1:120)
    err <- curve_matrix(dense) - y

    p <- coef(fit)
    jumps <- mean(abs(diff(log(p$tau1))) > log(2))

    cat(sprintf("\n%s (%.0fs with the hold-out check)\n", model, secs))
    cat(sprintf("  between tenors: RMSE %.1fbp, max %.1fbp\n",
                bp(sqrt(mean(err[, off]^2))), bp(max(abs(err[, off])))))
    cat(sprintf("  one-month yield, extrapolated: RMSE %.1fbp, max %.1fbp\n",
                bp(sqrt(mean(err[, "1"]^2))), bp(max(abs(err[, "1"])))))
    cat(sprintf("  tau1 doubles or halves month to month: %.0f%% of months\n",
                100 * jumps))

    # Does the hold-out check predict the true error between tenors? Compare
    # each held-out tenor's refit error with the actual error of the full
    # fit at the maturities in the gaps either side of it.
    chk <- fit$curves[[1]]$check
    for (i in seq_len(nrow(chk))) {
      h <- chk$maturity[i]
      k <- match(h, TENORS)
      gap <- off[off > TENORS[k - 1] & off < TENORS[k + 1]]
      cat(sprintf("  hold-out %3dm: RMSE %5.1fbp | true error in (%d, %d)m: ",
                  h, bp(chk$rmse[i]), TENORS[k - 1], TENORS[k + 1]))
      cat(sprintf("RMSE %5.1fbp\n",
                  bp(sqrt(mean(err[, as.character(gap)]^2)))))
    }

    rows[[length(rows) + 1]] <- compare(paste(model, "/ observed short rate"),
                                        truth, dense, n_factors, sr)
    rows[[length(rows) + 1]] <- compare(paste(model, "/ extrapolated 1m"),
                                        truth, dense, n_factors)
  }
  cat("\n10-year term premium against the fit on the original curve:\n")
  print(do.call(rbind, rows), row.names = FALSE)
}

# --- US: GSW, a Svensson curve ------------------------------------------

us <- yield_panel(gsw_monthly, units = "percent", maturity_unit = "months",
                  issuer = "US")
round_trip("US, GSW (Svensson truth: flattering)", us, n_factors = 5)

# --- UK: Bank of England spline -----------------------------------------

boe <- boe_yield_curve("nominal", cache = BOE_ZIP)
boe <- boe[boe$date >= as.Date("1997-03-01") &
           boe$date <= as.Date("2012-12-31"), ]
stopifnot(!any(boe$extrapolated), length(unique(boe$date)) == 190L)
uk <- yield_panel(boe, units = "percent", maturity_unit = "months",
                  issuer = "GB")
round_trip("UK, Bank of England (spline truth: the fair test)", uk,
           n_factors = 4)
