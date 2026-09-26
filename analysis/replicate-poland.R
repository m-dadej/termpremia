# Replicate the yieldcartography.com Polish ACM term premium
# =========================================================================
#
# Reproduces the 10-year Polish term premium published at
# https://yieldcartography.com/term-premia/ using this package's ACM
# implementation.
#
# Result on the reference data (1,132 weekly observations, 2005-2026):
#
#     correlation, levels          0.9997
#     correlation, weekly changes  0.9978
#     fit to their curve           0.30 bp RMSE
#     level offset                 +21.2 bp (ours above theirs)
#
# METHOD. Their site describes weekly Friday closes but monthly excess
# returns, which cannot both describe the estimation. The `pca_*` columns in
# their curve export settle it: every populated row is a month-end. So the
# parameters are estimated on MONTHLY data and the decomposition is then
# EVALUATED weekly -- valid because, once the parameters are known, the
# decomposition at any date needs only the factors on that date. That is what
# predict() does here.
#
# DATA. Two CSVs exported from yieldcartography.com:
#   - curve time series (full panel): needs `date` and the NSS parameters
#     beta0_pct..beta3_pct, tau1_years, tau2_years
#   - 10y term premium: needs `date`, `acm_bp`, `brw_bp`
#
# Their data is released for academic and personal use, which is a use
# restriction rather than a redistribution grant. It is therefore NOT bundled
# with this package, and this script reads it from local paths you supply.
#
# Usage: Rscript analysis/replicate-poland.R
# =========================================================================

# --- configuration -------------------------------------------------------

DATA_DIR    <- "C:/Users/Mateusz/Documents/R/test_package"
CURVE_CSV   <- file.path(DATA_DIR, "yc_pl_curve_timeseries.csv")
TP_CSV      <- file.path(DATA_DIR, "yc_pl_tp_10y_acm_brw.csv")

N_FACTORS   <- 5L        # their stated specification
MATURITIES  <- 1:120     # months; 120 = the 10y horizon they publish
HORIZON     <- 120L
MAKE_PLOT   <- TRUE
library(termpremia)


# --- read ----------------------------------------------------------------

curve <- read.csv(CURVE_CSV, stringsAsFactors = FALSE)
theirs <- read.csv(TP_CSV, stringsAsFactors = FALSE)

# Exported as M/D/YYYY, which as.Date() will not parse without a format.
curve$date <- as.Date(curve$date, format = "%m/%d/%Y")
theirs$date <- as.Date(theirs$date, format = "%m/%d/%Y")

needed <- c("beta0_pct", "beta1_pct", "beta2_pct", "beta3_pct",
            "tau1_years", "tau2_years")
stopifnot(all(needed %in% names(curve)))

cat(sprintf("curve : %d obs, %s to %s\n", nrow(curve),
            format(min(curve$date)), format(max(curve$date))))
cat(sprintf("theirs: %d obs\n\n", nrow(theirs)))

# --- evaluate the NSS curve on a monthly maturity grid -------------------

# Their export publishes only 2y/5y/10y zero rates, but ACM needs a dense grid
# including the short end, so the curve is rebuilt from the NSS parameters.
params <- data.frame(
  beta0 = curve$beta0_pct, beta1 = curve$beta1_pct,
  beta2 = curve$beta2_pct, beta3 = curve$beta3_pct,
  tau1  = curve$tau1_years, tau2 = curve$tau2_years
)

yields_pct <- svensson_curve(params, maturity = MATURITIES / 12)
rownames(yields_pct) <- as.character(curve$date)
colnames(yields_pct) <- as.character(MATURITIES)

stopifnot(all(is.finite(yields_pct)))

# Sanity check: our reconstruction must reproduce their own published tenors.
# If this fails, nothing downstream is trustworthy.
for (chk in list(c(24, "y2_pct"), c(60, "y5_pct"), c(120, "y10_pct"))) {
  m <- as.integer(chk[1]); col <- chk[2]
  if (!col %in% names(curve)) next
  worst <- max(abs(yields_pct[, m] - curve[[col]]))
  cat(sprintf("NSS check %3sm vs %-8s : max abs diff %.2e\n", chk[1], col, worst))
  stopifnot(worst < 1e-3)
}
cat("\n")

# --- build panels --------------------------------------------------------

weekly <- yield_panel(
  yields_pct, units = "percent", maturity_unit = "months",
  curve_name = "PL", instrument = "government", issuer = "PL"
)

# Month-ends: the last available observation within each calendar month.
eom_rows <- sort(unname(
  tapply(seq_len(nrow(curve)), format(curve$date, "%Y-%m"), max)
))
monthly <- yield_panel(
  yields_pct[eom_rows, , drop = FALSE],
  units = "percent", maturity_unit = "months",
  curve_name = "PL", instrument = "government", issuer = "PL"
)

cat(sprintf("panels: %d monthly (estimation), %d weekly (evaluation)\n\n",
            length(monthly$dates), length(weekly$dates)))

# --- estimate monthly, evaluate weekly -----------------------------------

fit <- atsm(monthly, n_factors = N_FACTORS, maturities = MATURITIES)
print(fit)

# Poland sits close to the edge of stationarity: 261 monthly observations is a
# short sample for a near-unit-root VAR, and the risk-adjusted spectral radius
# comes out around 1.001. At a 10-year horizon that compounds to only ~1.17, so
# it does not corrupt these results, but longer maturities would be unreliable.
# It is also why the site applies Bauer-Rudebusch-Wu bias correction to Poland
# and not to the US.
cat(sprintf("\nspectral radius: risk-adjusted %.4f, real-world %.4f\n",
            fit$spectral_radius[["risk_adjusted"]],
            fit$spectral_radius[["real_world"]]))

# The same model with the small-sample bias in the factor VAR corrected. This
# is the second series they publish, and the reason this script fits both: 261
# monthly observations of a near-unit-root process is exactly the case the
# correction exists for.
fit_brw <- atsm(monthly, n_factors = N_FACTORS, maturities = MATURITIES,
                p_dynamics = "brw")
cat("\n")
print(fit_brw)

pred <- predict(fit, weekly)
ours <- pred[pred$maturity == HORIZON, c("date", "term_premium")]

pred_brw <- predict(fit_brw, weekly)
ours_brw <- pred_brw[pred_brw$maturity == HORIZON, c("date", "term_premium")]
names(ours_brw)[2] <- "term_premium_brw"

# --- compare -------------------------------------------------------------

cmp <- merge(merge(ours, ours_brw, by = "date"),
             theirs[, c("date", "acm_bp", "brw_bp")], by = "date")
cmp <- cmp[complete.cases(cmp), ]          # their series has one missing row
cmp <- cmp[order(cmp$date), ]
cmp$ours_bp <- cmp$term_premium * 1e4      # decimals -> basis points
cmp$ours_brw_bp <- cmp$term_premium_brw * 1e4

report <- function(label, mine, theirs_bp) {
  cat(sprintf(
    "%-22s corr_levels %.4f  corr_changes %.4f  mean %6.1f vs %6.1f  offset %+6.1f bp\n",
    label, cor(mine, theirs_bp), cor(diff(mine), diff(theirs_bp)),
    mean(mine), mean(theirs_bp), mean(mine - theirs_bp)
  ))
}

cat(sprintf("\n=== 10y term premium, %d weekly observations ===\n", nrow(cmp)))
cat("-- our OLS estimate --\n")
report("vs their ACM", cmp$ours_bp, cmp$acm_bp)
report("vs their BRW", cmp$ours_bp, cmp$brw_bp)
cat("-- our bias-corrected estimate --\n")
report("vs their ACM", cmp$ours_brw_bp, cmp$acm_bp)
report("vs their BRW", cmp$ours_brw_bp, cmp$brw_bp)

cat("\nLevels and changes are reported separately on purpose. Implementations\n")
cat("of the same model routinely agree closely on direction while differing\n")
cat("substantially in level; see Cohen, Hordahl & Xia (BIS QR, Sept 2018).\n")
cat("Our offset traces to the real-world VAR, the least-identified part of\n")
cat("the model, not to the pricing side.\n")

# --- plot ----------------------------------------------------------------

if (MAKE_PLOT) {
  op <- par(mfrow = c(2, 1), mar = c(3, 4, 2, 1))

  plot(cmp$date, cmp$acm_bp, type = "l", lwd = 2, col = "grey40",
       xlab = "", ylab = "bp", main = "Polish 10y term premium",
       ylim = range(cmp$acm_bp, cmp$brw_bp, cmp$ours_bp, cmp$ours_brw_bp))
  lines(cmp$date, cmp$brw_bp, lwd = 2, col = "grey70")
  lines(cmp$date, cmp$ours_bp, lwd = 1.5, col = "firebrick")
  lines(cmp$date, cmp$ours_brw_bp, lwd = 1.5, col = "steelblue")
  abline(h = 0, col = "grey80", lty = 2)
  legend("topleft",
         c("yieldcartography ACM", "yieldcartography BRW",
           "termpremia OLS", "termpremia BRW"),
         col = c("grey40", "grey70", "firebrick", "steelblue"),
         lwd = c(2, 2, 1.5, 1.5), bty = "n", cex = 0.8)

  plot(cmp$date, cmp$ours_bp - cmp$acm_bp, type = "l", col = "firebrick",
       xlab = "", ylab = "bp", main = "Difference (ours - theirs)",
       ylim = range(cmp$ours_bp - cmp$acm_bp, cmp$ours_brw_bp - cmp$brw_bp))
  lines(cmp$date, cmp$ours_brw_bp - cmp$brw_bp, col = "steelblue")
  abline(h = 0, col = "grey80", lty = 2)
  legend("topleft", c("OLS vs their ACM", "BRW vs their BRW"),
         col = c("firebrick", "steelblue"), lwd = 1, bty = "n", cex = 0.8)

  par(op)
}

invisible(list(fit = fit, comparison = cmp))
