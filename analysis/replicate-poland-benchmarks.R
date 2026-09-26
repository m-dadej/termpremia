# Replicate the yieldcartography.com Polish ACM term premium from benchmarks
# =========================================================================
#
# Estimates the Polish 10-year ACM term premium from five sparse rates --
# 3m, 6m, 2y, 5y, 10y, monthly averages -- and compares it with the series
# published at https://yieldcartography.com/term-premia/. The companion
# script, replicate-poland.R, does the same from the site's own zero-coupon
# curve. This one asks how close you get without that curve.
#
# Result, 10-year term premium against the site's ACM averaged by calendar
# month, 2005:01-2026:08:
#
#                 months   corr, levels   corr, changes   mean gap   sd of gap
#   all months       260       0.993           0.943        -0.9bp      13.0bp
#   2005-2015        132       0.990           0.934        -2.6bp      11.9bp
#   2016-2026        128       0.994           0.947        +0.9bp      13.9bp
#
#   largest gap +51bp, in 2022:08
#
# For scale, the site's own curve averaged the same way reaches 0.995 on
# changes, and at month-end, as in replicate-poland.R, 0.998.
#
# DATA. data/polgovs.xlsx holds monthly averages of daily rates, in percent,
# dated the first of each month:
#   - 1month, 3month, 6month: money-market rates. They sit 9, 20 and 25bp
#     above the site's government curve on average, so they behave like
#     interbank rates, not bills.
#   - 1year, 2year, 5year, 7year, 10year: coupon benchmark yields. 7year
#     starts only in 2018:08.
#   - 1year_es: not identified; unused.
# The site's term premium export (needs `date`, `acm_bp`) is read from a local
# path, as in replicate-poland.R, and used only for the comparison.
#
# WHY THESE FIVE RATES.
#
# 1. 1month is dropped. The money-market rates are not what the site's short
#    end measures, and the shortest one pulls the fitted curve's short end
#    furthest towards them.
# 2. 1year is dropped. It is the noisiest series against the site's curve:
#    10bp below its 1-year zero rate on average, with an 18bp standard
#    deviation, against 6-8bp for 2y, 5y and 10y.
# 3. 7year is not used. Adding it in 2018 would change how the curve is
#    fitted, and with it the factors, part of the way through the sample.
# 4. Five rates support Nelson-Siegel (four parameters) but not Svensson
#    (six). On all seven tenors Svensson bent towards the money-market rates
#    -- median first decay 0.27 years, against 1.16 in the site's fits -- and
#    matched worst of everything tried.
#
# THE BENCHMARKS ARE COUPON YIELDS, DECLARED "zero" KNOWINGLY. ACM prices
# zero-coupon bonds. Fitting the 2y-10y yields properly as annual-coupon par
# yields, with the money-market rates made continuously compounded, matched
# the site less well (0.921 on changes, against 0.943), perhaps because the
# spreadsheet's quoting conventions are not the ones that fit assumed. For
# scale: in a simulation on the site's curve, treating par yields as zero
# yields moved this premium by +7.5bp on average, with a 10bp standard
# deviation.
#
# WHAT WAS TRIED AND DID NOT HELP, full sample:
#
#                                                 corr, changes   sd of gap   mean gap
#   Svensson, all seven rates (the default)            0.835        14.0bp    +13.2bp
#   Nelson-Siegel, all seven rates                     0.882        16.6bp     -5.8bp
#   Nelson-Siegel, all but 1year                       0.908        17.2bp     -4.3bp
#   2y-10y fitted as par yields (above)                0.921        13.8bp     -8.9bp
#   money-market rates at a quarter of the weight   0.875-0.889   15-17bp   -11 to -37bp
#   1y-10y only, short end extrapolated                0.939         9.6bp    +62.8bp
#
# A linear map from the rates onto the site's curve, trained on one half of
# the sample and tested on the other, got the short end of the curve within
# 10-12bp but was not stable on the premium: 0.971 on changes in one half,
# 0.903 in the other.
#
# WHAT LIMITS THE MATCH.
#
# 1. The short end. The site fits its curve to 6-18 government bonds and
#    extrapolates 1-12 months from those near maturity. Nothing in the
#    spreadsheet measures that. The rebuilt curve is within about 10bp of the
#    site's from 13 to 120 months, but about 30bp away from 1 to 12 months.
# 2. Monthly averages. On the site's own curve, estimating on monthly averages
#    instead of month-ends moves the premium by 12bp (standard deviation) and
#    caps the match at 0.995 on changes.
# 3. The implementation offset. On the site's own curve at month-end this
#    package runs 21bp above the site. The small mean gap here is offsetting
#    errors -- about +21bp implementation, +3bp averaging, -25bp data -- not
#    a level match.
# 4. The model sits at the edge of stability on this sample: risk-adjusted
#    spectral radius about 1.00. Small differences in the inputs move the
#    10-year premium noticeably.
#
# The tenor choice was made by comparing with the site, so the match above is
# somewhat flattering. It holds about equally in both halves of the sample.
#
# BRW is not fitted here. From these rates the bias-corrected premium matched
# the site's BRW series at only 0.2-0.5 on monthly changes: the correction
# makes the factor dynamics more persistent and amplifies input differences.
#
# What would close the gap: month-end prices, coupons and maturities of the
# individual fixed-rate bonds (the site uses BondSpot closing fixings), with a
# curve fitted to prices; failing that, month-end benchmark yields instead of
# monthly averages.
#
# Usage: Rscript analysis/replicate-poland-benchmarks.R
# =========================================================================

# --- configuration -------------------------------------------------------

XLSX        <- "data/polgovs.xlsx"
DATA_DIR    <- "C:/Users/Mateusz/Documents/R/test_package"
TP_CSV      <- file.path(DATA_DIR, "yc_pl_tp_10y_acm_brw.csv")

# Spreadsheet columns used, with their maturities in months.
TENORS      <- c("3month" = 3, "6month" = 6, "2year" = 24, "5year" = 60,
                 "10year" = 120)
START       <- as.Date("2005-01-01")   # the site's sample starts here
N_FACTORS   <- 5L                      # the site's stated specification
MATURITIES  <- 1:120
HORIZON     <- 120L
MAKE_PLOT   <- TRUE
library(termpremia)


# --- read ----------------------------------------------------------------

raw <- as.data.frame(readxl::read_excel(XLSX))
stopifnot("date" %in% names(raw), all(names(TENORS) %in% names(raw)))

x <- data.frame(date = as.Date(raw$date))
for (col in names(TENORS)) x[[col]] <- as.numeric(raw[[col]])
x <- x[x$date >= START & complete.cases(x), ]
x <- x[order(x$date), ]

# ACM forms one-month holding returns, so the months must be consecutive.
month_index <- as.integer(format(x$date, "%Y")) * 12L +
  as.integer(format(x$date, "%m"))
stopifnot(all(diff(month_index) == 1L))

cat(sprintf("rates: %d months, %s to %s\n\n", nrow(x),
            format(min(x$date), "%Y-%m"), format(max(x$date), "%Y-%m")))

# --- rebuild the monthly maturity grid -----------------------------------

sparse <- as.matrix(x[, names(TENORS)])
dimnames(sparse) <- list(as.character(x$date), TENORS)

# "other": money-market rates at the short end, government benchmarks beyond.
sparse_panel <- yield_panel(
  sparse, units = "percent", maturity_unit = "months",
  curve_name = "PL", instrument = "other", issuer = "PL"
)

# No hold-out check: Nelson-Siegel needs five rates, so none can be dropped.
curve_fit <- svensson_fit(sparse_panel, yield_type = "zero",
                          model = "nelson_siegel", check = FALSE)
print(curve_fit)

dense <- predict(curve_fit, maturities = MATURITIES)

# --- estimate ------------------------------------------------------------

# The one-month short rate is the fitted curve's extrapolation below 3m, and
# atsm() warns about that. Here it is deliberate: the site's short end is an
# extrapolation too, from government bonds, and the observed 1month rate is an
# interbank rate. Supplying that rate through `short_rate` instead changed
# little: 0.936 on changes against 0.943, mean gap +9.4bp against -0.9bp.
fit <- withCallingHandlers(
  atsm(dense, n_factors = N_FACTORS, maturities = MATURITIES),
  warning = function(w) {
    if (grepl("The short rate is the one-month yield", conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  }
)
cat("\n")
print(fit)
cat(sprintf("\nspectral radius: risk-adjusted %.4f, real-world %.4f\n",
            fit$spectral_radius[["risk_adjusted"]],
            fit$spectral_radius[["real_world"]]))

pred <- predict(fit, dense)
ours <- pred[pred$maturity == HORIZON, c("date", "term_premium")]

# --- compare -------------------------------------------------------------

theirs <- read.csv(TP_CSV, stringsAsFactors = FALSE)
theirs$date <- as.Date(theirs$date, format = "%m/%d/%Y")

# The spreadsheet holds monthly averages, so the site's weekly series is
# averaged by calendar month too. With the parameters fixed, the premium is
# affine in the yields, so the premium of an averaged curve is the average
# premium: this is the like-for-like target, and the site's month-end value
# is not. The formula interface drops the one missing week.
theirs$month <- as.Date(format(theirs$date, "%Y-%m-01"))
theirs_m <- aggregate(acm_bp ~ month, data = theirs, FUN = mean)
names(theirs_m)[1] <- "date"

cmp <- merge(ours, theirs_m, by = "date")
cmp <- cmp[order(cmp$date), ]
cmp$ours_bp <- cmp$term_premium * 1e4      # decimals -> basis points
cmp$gap_bp <- cmp$ours_bp - cmp$acm_bp

report <- function(label, keep) {
  a <- cmp$ours_bp[keep]
  b <- cmp$acm_bp[keep]
  cat(sprintf("%-11s %6d   %.3f    %.3f   %+6.1f bp  %5.1f bp\n",
              label, sum(keep), cor(a, b), cor(diff(a), diff(b)),
              mean(a - b), sd(a - b)))
}

second_half <- cmp$date >= as.Date("2016-01-01")
cat("\n=== 10y term premium vs the site's ACM, monthly averages ===\n")
cat("period      months   corr_lvl  corr_chg   mean gap    sd gap\n")
report("all", rep(TRUE, nrow(cmp)))
report("2005-2015", !second_half)
report(paste0("2016-", format(max(cmp$date), "%Y")), second_half)

worst <- which.max(abs(cmp$gap_bp))
cat(sprintf("largest gap: %+.0f bp in %s\n", cmp$gap_bp[worst],
            format(cmp$date[worst], "%Y-%m")))

cat("\nThe mean gap is small by coincidence of offsetting errors, not because\n")
cat("the levels agree; see the header. On the site's own curve the same model\n")
cat("runs about 21bp above the site.\n")

# --- plot ----------------------------------------------------------------

if (MAKE_PLOT) {
  op <- par(mfrow = c(2, 1), mar = c(3, 4, 2, 1))

  plot(cmp$date, cmp$acm_bp, type = "l", lwd = 2, col = "grey40",
       xlab = "", ylab = "bp",
       main = "Polish 10y term premium, monthly averages",
       ylim = range(cmp$acm_bp, cmp$ours_bp))
  lines(cmp$date, cmp$ours_bp, lwd = 1.5, col = "firebrick")
  abline(h = 0, col = "grey80", lty = 2)
  legend("topleft",
         c("yieldcartography ACM", "termpremia, 3m 6m 2y 5y 10y"),
         col = c("grey40", "firebrick"), lwd = c(2, 1.5), bty = "n", cex = 0.8)

  plot(cmp$date, cmp$gap_bp, type = "l", col = "firebrick",
       xlab = "", ylab = "bp", main = "Difference (ours - theirs)")
  abline(h = 0, col = "grey80", lty = 2)

  par(op)
}

invisible(list(fit = fit, curve_fit = curve_fit, comparison = cmp))
