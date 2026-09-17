# Build the bundled 3-month Treasury bill snapshot.
#
# Source: Federal Reserve Board H.15 Selected Interest Rates, distributed as
# FRED series TB3MS (3-Month Treasury Bill Secondary Market Rate, monthly
# average, 1934 onwards). A Federal Reserve Board statistical release; no
# copyright is asserted, the same provenance as the GSW curve in data-raw/gsw.R.
#
# Why bundle it: the short rate is a modelling choice with real consequences.
# GSW exclude Treasury bills and every security with under three months to
# maturity, so their fitted short end is extrapolated and they caution against
# it. Against this bill series, the GSW fitted one-month yield averages about
# 28bp HIGHER, with a standard deviation of 127bp. ACM's own stated inputs
# include H.15, so an actual bill rate is not an exotic alternative -- it is
# closer to what the published model uses.
#
# Run with: source("data-raw/tbill.R")

SOURCE_URL <- "https://fred.stlouisfed.org/graph/fredgraph.csv?id=TB3MS"
SNAPSHOT_END <- as.Date("2024-12-31")   # aligned with data-raw/gsw.R

cache_dir <- file.path("data-raw", ".cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
cache <- file.path(cache_dir, "TB3MS.csv")

if (!file.exists(cache)) {
  message("Downloading ", SOURCE_URL)
  utils::download.file(SOURCE_URL, cache, mode = "wb", quiet = TRUE)
}

raw <- utils::read.csv(cache, stringsAsFactors = FALSE)
names(raw) <- c("date", "value")
raw$date <- as.Date(raw$date)
raw$value <- suppressWarnings(as.numeric(raw$value))

raw <- raw[!is.na(raw$value) & raw$date <= SNAPSHOT_END, , drop = FALSE]

# FRED stamps monthly observations to the first of the month. The yield panels
# in this package are month-END, so key on the calendar month rather than the
# day, and re-stamp to the panel's own month-end dates.
load("data/gsw_monthly.rda")
panel_dates <- sort(unique(gsw_monthly$date))

idx <- match(format(panel_dates, "%Y-%m"), format(raw$date, "%Y-%m"))
if (anyNA(idx)) {
  stop("TB3MS does not cover ", sum(is.na(idx)), " panel month(s), first: ",
       format(panel_dates[which(is.na(idx))[1]]))
}

tbill_3m <- data.frame(
  date = panel_dates,
  value = raw$value[idx],
  stringsAsFactors = FALSE
)

attr(tbill_3m, "source_url") <- SOURCE_URL
attr(tbill_3m, "retrieved") <- Sys.Date()

message(sprintf(
  "tbill_3m: %d obs, %s to %s, mean %.2f%%",
  nrow(tbill_3m), format(min(tbill_3m$date)), format(max(tbill_3m$date)),
  mean(tbill_3m$value)
))

usethis::use_data(tbill_3m, overwrite = TRUE, compress = "xz")
