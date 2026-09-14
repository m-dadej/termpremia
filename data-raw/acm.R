# Build the bundled New York Fed ACM reference snapshot.
#
# Source: Federal Reserve Bank of New York, "Treasury Term Premia", the
# published output of the five-factor Adrian, Crump and Moench (2013) model.
#
# LICENSING. The New York Fed's Terms of Use permit copying and distributing
# their content provided that: (a) source identifiers and copyright notices are
# retained, (b) modifications are clearly labelled, and (c) the same permissions
# are passed on to downstream users. See LICENSE.note at the package root, which
# records the attribution string and states that this dataset is governed by the
# New York Fed's terms rather than by the package's MIT licence.
#
# MODIFICATIONS MADE HERE, declared per (b) above:
#   1. Only the "ACM Monthly" sheet is retained; the daily sheet is dropped.
#   2. The series are reshaped from wide (one column per tenor) to long.
#   3. The sample is truncated at SNAPSHOT_END so the bundled data is stable.
# No values are altered.
#
# This snapshot exists to validate our own ACM implementation against the
# published numbers. It is a reference, not an input to estimation.
#
# `readxl` is used here only; it is deliberately NOT a package dependency,
# because data-raw/ is excluded from the build.
#
# Run with: source("data-raw/acm.R")

SOURCE_URL <- paste0(
  "https://www.newyorkfed.org/medialibrary/media/research/",
  "data_indicators/ACMTermPremium.xls"
)
SNAPSHOT_END <- as.Date("2024-12-31")  # aligned with data-raw/gsw.R

stopifnot(requireNamespace("readxl", quietly = TRUE))

# --- download -----------------------------------------------------------

cache_dir <- file.path("data-raw", ".cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
cache <- file.path(cache_dir, "ACMTermPremium.xls")

if (!file.exists(cache)) {
  message("Downloading ", SOURCE_URL)
  utils::download.file(SOURCE_URL, cache, mode = "wb", quiet = TRUE)
}

wide <- as.data.frame(readxl::read_excel(cache, sheet = "ACM Monthly"))

# --- reshape ------------------------------------------------------------

# Published dates look like "30-Jun-1961"
wide$DATE <- as.Date(wide$DATE, format = "%d-%b-%Y")
stopifnot(!anyNA(wide$DATE))

wide <- wide[wide$DATE <= SNAPSHOT_END, , drop = FALSE]
wide <- wide[order(wide$DATE), , drop = FALSE]

# ACMY = model-fitted yield, ACMTP = term premium,
# ACMRNY = risk-neutral (expected average short rate) yield.
series <- c(fitted = "ACMY", term_premium = "ACMTP", risk_neutral = "ACMRNY")
tenors <- 1:10

for (s in series) {
  missing <- setdiff(sprintf("%s%02d", s, tenors), names(wide))
  if (length(missing)) {
    stop("Expected column(s) absent from the workbook: ",
         paste(missing, collapse = ", "))
  }
}

pull <- function(prefix) {
  cols <- sprintf("%s%02d", prefix, tenors)
  as.vector(as.matrix(wide[, cols]))
}

acm_published <- data.frame(
  date = rep(wide$DATE, times = length(tenors)),
  maturity = rep(tenors * 12L, each = nrow(wide)),
  fitted = pull("ACMY"),
  term_premium = pull("ACMTP"),
  risk_neutral = pull("ACMRNY"),
  stringsAsFactors = FALSE
)
acm_published <- acm_published[order(acm_published$date, acm_published$maturity), ]
rownames(acm_published) <- NULL

# --- internal consistency of the published series -----------------------

# The decomposition must add up in the source data. If it does not, our
# understanding of the columns is wrong and nothing downstream can be trusted.
gap <- with(acm_published, fitted - (risk_neutral + term_premium))
worst <- max(abs(gap), na.rm = TRUE)
message(sprintf("Published decomposition residual: max |fitted - (rn + tp)| = %.2e",
                worst))
stopifnot(worst < 1e-6)

attr(acm_published, "source_url") <- SOURCE_URL
attr(acm_published, "retrieved") <- Sys.Date()
attr(acm_published, "snapshot_end") <- SNAPSHOT_END
attr(acm_published, "attribution") <- paste0(
  "© ", format(Sys.Date(), "%Y"), " Federal Reserve Bank of New York. ",
  "Content from the New York Fed subject to the Terms of Use at newyorkfed.org."
)

message(
  "acm_published: ", format(nrow(acm_published), big.mark = ","), " rows, ",
  length(unique(acm_published$date)), " dates (",
  format(min(acm_published$date)), " to ", format(max(acm_published$date)), ")."
)

usethis::use_data(acm_published, overwrite = TRUE, compress = "xz")
