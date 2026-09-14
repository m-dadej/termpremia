# Build the bundled GSW nominal Treasury zero-coupon yield snapshot.
#
# Source: Gurkaynak, Sack and Wright (2007), "The U.S. Treasury yield curve:
# 1961 to the present", maintained by the Federal Reserve Board as series
# feds200628. A Federal Reserve Board staff research product; no copyright is
# asserted on the page, and it is not an official statistical release.
#
# Why this script exists: the published file cannot be used directly.
#
#   1. Nine lines of preamble precede the header.
#   2. The published SVENY columns run SVENY01..SVENY30 -- ANNUAL steps,
#      starting at one year. ACM needs a monthly maturity grid including the
#      short end, so the curve must be re-evaluated from the Svensson
#      parameters rather than read off.
#   3. Rows before 1980 are Nelson-Siegel, encoded as beta3 = 0 with a tau2
#      sentinel of -999.99.
#   4. Early in the sample only a few maturities were actually traded. The
#      parameters still define the whole curve, so evaluating them at 10 years
#      in 1961 produces an EXTRAPOLATION beyond the longest traded bond. We
#      record that explicitly rather than pass it off as observed data.
#
# Run with: source("data-raw/gsw.R")

pkgload::load_all(".", quiet = TRUE)

SOURCE_URL <- "https://www.federalreserve.gov/data/yield-curve-tables/feds200628.csv"
SNAPSHOT_END <- as.Date("2024-12-31")  # frozen, so the bundled data is stable
MATURITIES_M <- 1:120                  # months; ACM's grid

# --- download -----------------------------------------------------------

# Cached on disk rather than in tempdir() so that re-running this script, or
# running it alongside data-raw/acm.R, does not re-download 10MB each time.
cache_dir <- file.path("data-raw", ".cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
cache <- file.path(cache_dir, "feds200628.csv")

if (!file.exists(cache)) {
  message("Downloading ", SOURCE_URL)
  utils::download.file(SOURCE_URL, cache, mode = "wb", quiet = TRUE)
}

raw <- utils::read.csv(cache, skip = 9, stringsAsFactors = FALSE)
raw$Date <- as.Date(raw$Date)

# --- end-of-month sampling ----------------------------------------------

raw <- raw[raw$Date <= SNAPSHOT_END, , drop = FALSE]
raw <- raw[order(raw$Date), , drop = FALSE]

# Drop dates with no fitted curve BEFORE picking month ends. Doing it the other
# way round silently loses an entire month whenever the last row of that month
# happens to carry no parameters -- which costs 20 months against the NY Fed's
# monthly series, concentrated in March and May.
has_params <- stats::complete.cases(raw[, c("BETA0", "BETA1", "BETA2", "TAU1")])
raw <- raw[has_params, , drop = FALSE]

ym <- format(raw$Date, "%Y-%m")
eom_idx <- tapply(seq_len(nrow(raw)), ym, max)
eom <- raw[sort(unname(eom_idx)), , drop = FALSE]

# --- longest maturity actually published, per date ----------------------

sveny_cols <- sprintf("SVENY%02d", 1:30)
published <- as.matrix(eom[, sveny_cols])
max_published_m <- 12 * apply(published, 1, function(r) {
  ok <- which(is.finite(r))
  if (!length(ok)) NA_integer_ else max(ok)
})

# --- evaluate the curve on a monthly grid -------------------------------

params <- data.frame(
  beta0 = eom$BETA0, beta1 = eom$BETA1, beta2 = eom$BETA2,
  beta3 = eom$BETA3, tau1 = eom$TAU1,  tau2 = eom$TAU2
)

yields_pct <- svensson_curve(params, maturity = MATURITIES_M / 12)

# --- sanity check: do we reproduce GSW's own published yields? ----------

# Only tenors that are actually on the grid can be checked -- asking for 20y
# against a 120-month grid silently compares nothing and "passes".
check_years <- c(1, 2, 5, 10)
stopifnot(all(check_years * 12 <= max(MATURITIES_M)))

for (yr in check_years) {
  col <- which(MATURITIES_M == yr * 12)
  stopifnot(length(col) == 1L)

  mine <- yields_pct[, col]
  theirs <- eom[[sprintf("SVENY%02d", yr)]]
  ok <- is.finite(mine) & is.finite(theirs)

  stopifnot(sum(ok) > 0)  # refuse to pass vacuously
  worst <- max(abs(mine[ok] - theirs[ok]))
  stopifnot(worst < 1e-3)

  message(sprintf("  %2dy: %5d overlapping obs, max abs diff %.2e",
                  yr, sum(ok), worst))
}
message("Svensson evaluator matches published SVENY at all checked tenors.")

# --- assemble -----------------------------------------------------------

gsw_monthly <- data.frame(
  date = rep(eom$Date, times = length(MATURITIES_M)),
  maturity = rep(MATURITIES_M, each = nrow(eom)),
  yield = as.vector(yields_pct),
  extrapolated = rep(MATURITIES_M, each = nrow(eom)) >
    rep(max_published_m, times = length(MATURITIES_M)),
  stringsAsFactors = FALSE
)
gsw_monthly <- gsw_monthly[order(gsw_monthly$date, gsw_monthly$maturity), ]
rownames(gsw_monthly) <- NULL

attr(gsw_monthly, "source_url") <- SOURCE_URL
attr(gsw_monthly, "retrieved") <- Sys.Date()
attr(gsw_monthly, "snapshot_end") <- SNAPSHOT_END

message(
  "gsw_monthly: ", format(nrow(gsw_monthly), big.mark = ","), " rows, ",
  length(unique(gsw_monthly$date)), " dates (",
  format(min(gsw_monthly$date)), " to ", format(max(gsw_monthly$date)), "), ",
  round(100 * mean(gsw_monthly$extrapolated), 1), "% extrapolated."
)

usethis::use_data(gsw_monthly, overwrite = TRUE, compress = "xz")
