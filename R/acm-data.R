# The New York Fed's published ACM series, fetched rather than bundled.
#
# Unlike the Federal Reserve Board series behind gsw_monthly, this one IS
# withheld for licensing reasons: the New York Fed's Terms of Use are not a
# licence CRAN recognises, they can change without notice, and they would have
# to pass to every downstream user of the package. Downloading it on the
# user's behalf puts them in a direct relationship with those terms instead.

#' Download the New York Fed's published ACM term premium estimates
#'
#' Fetches the Federal Reserve Bank of New York's published decomposition of
#' US Treasury yields into an expected-average-short-rate component and a term
#' premium, from the five-factor model of Adrian, Crump and Moench (2013), and
#' reshapes its monthly sheet to one row per date and maturity.
#'
#' This is a **reference for validation**, not an input to estimation. It is
#' what this package's own ACM implementation is checked against; see
#' `vignette("validation")`.
#'
#' @section Licensing:
#' Downloaded, never redistributed. The data are the New York Fed's and are
#' governed by its Terms of Use, <https://www.newyorkfed.org/privacy/termsofuse>,
#' not by this package's MIT licence. Those terms permit copying and
#' redistribution subject to attribution, clear labelling of any
#' modifications, and passing the same terms on. The attribution string is
#' attached to the result as its `"attribution"` attribute. See `LICENSE.note`.
#'
#' Modifications made here: only the monthly sheet is read, the series are
#' reshaped from wide to long, and the sample is optionally truncated at `end`.
#' No values are altered.
#'
#' @section The published series are revised:
#' The New York Fed re-estimates the model over the full sample at each update,
#' so historical values move from one release to the next: by up to about 2bp
#' on the 10-year term premium between two releases in 2026. Pass `end`, and
#' keep a copy through `cache`, when a comparison has to be reproducible.
#'
#' @section The decomposition is exact:
#' `fitted == risk_neutral + term_premium` holds to machine precision in the
#' published data. It is checked on every read: if it fails, the columns are
#' not what this function assumes, and it stops rather than return them.
#'
#' @param end Optional last date to keep, inclusive.
#' @param url Location of the New York Fed's workbook.
#' @param cache Optional path to a downloaded copy. When it exists it is read
#'   instead of downloading; when it does not, the download is saved there.
#' @param quiet Suppress the download progress bar.
#'
#' @return A data frame with one row per month-end and maturity:
#'   \describe{
#'     \item{date}{Month-end observation date.}
#'     \item{maturity}{Maturity in months: 12, 24, ..., 120.}
#'     \item{fitted}{Model-fitted zero-coupon yield, in percent (`ACMY`).}
#'     \item{term_premium}{Term premium, in percent (`ACMTP`).}
#'     \item{risk_neutral}{Risk-neutral yield, that is the expected average
#'       future short rate over the life of the bond, in percent (`ACMRNY`).}
#'   }
#'   with attributes `source_url` and `attribution`.
#'
#' @source
#' Federal Reserve Bank of New York, Treasury Term Premia,
#' <https://www.newyorkfed.org/research/data_indicators/term-premia-tabs>.
#'
#' @references
#' Adrian, T., R. K. Crump and E. Moench (2013). "Pricing the term structure
#' with linear regressions." *Journal of Financial Economics* 110(1), 110-138.
#'
#' @seealso [atsm()], [gsw_monthly]
#'
#' @examplesIf interactive()
#' acm <- acm_published(end = "2024-12-31")
#'
#' # The published decomposition adds up exactly
#' with(acm, max(abs(fitted - (risk_neutral + term_premium))))
#'
#' # Against this package's estimate on the same sample
#' panel <- yield_panel(gsw_monthly, units = "percent",
#'                      maturity_unit = "months", issuer = "US")
#' fit <- atsm(panel, n_factors = 5)
#' cmp <- merge(term_premium(fit, maturity = 120),
#'              acm[acm$maturity == 120, c("date", "term_premium")])
#' cor(cmp$value, cmp$term_premium)
#'
#' @export
acm_published <- function(end = NULL,
                          url = paste0(
                            "https://www.newyorkfed.org/medialibrary/media/",
                            "research/data_indicators/ACMTermPremium.xls"
                          ),
                          cache = NULL,
                          quiet = TRUE) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("`acm_published()` needs the readxl package: ",
         "install.packages(\"readxl\").", call. = FALSE)
  }

  path <- fetch_to_file(url, cache, quiet,
                        what = "the New York Fed's ACM term premia",
                        fileext = ".xls")

  wide <- tryCatch(
    as.data.frame(readxl::read_excel(path, sheet = "ACM Monthly")),
    error = function(e) {
      stop("Could not read the sheet \"ACM Monthly\" from\n  ", path, "\n",
           conditionMessage(e),
           "\nThe download may have returned an error page instead of the ",
           "workbook, or the New York Fed may have changed its layout.",
           call. = FALSE)
    }
  )

  out <- acm_reshape(wide, end)
  attr(out, "source_url") <- url
  # The copyright sign is built from its code point: R code must be ASCII.
  attr(out, "attribution") <- paste0(
    intToUtf8(169L), " ", format(Sys.Date(), "%Y"),
    " Federal Reserve Bank of New York. ",
    "Content from the New York Fed subject to the Terms of Use at ",
    "newyorkfed.org."
  )
  out
}


#' Reshape the New York Fed's monthly ACM sheet to long form
#'
#' @param wide The sheet as read: `DATE`, then `ACMY01`-`ACMY10`,
#'   `ACMTP01`-`ACMTP10` and `ACMRNY01`-`ACMRNY10`, in percent.
#' @param end Optional last date to keep.
#' @return A data frame ordered by date, then maturity.
#' @keywords internal
#' @noRd
acm_reshape <- function(wide, end = NULL) {
  tenors <- 1:10
  series <- c(fitted = "ACMY", term_premium = "ACMTP",
              risk_neutral = "ACMRNY")

  needed <- c("DATE", as.vector(outer(series, sprintf("%02d", tenors),
                                      paste0)))
  missing_cols <- setdiff(needed, names(wide))
  if (length(missing_cols)) {
    stop("The ACM sheet is missing column(s): ",
         paste(utils::head(missing_cols, 6L), collapse = ", "),
         if (length(missing_cols) > 6L) ", ...",
         ". The file layout may have changed.", call. = FALSE)
  }

  date <- parse_acm_date(wide$DATE)
  if (anyNA(date)) {
    bad <- utils::head(as.character(wide$DATE[is.na(date)]), 3L)
    stop("Could not read ", sum(is.na(date)), " date(s) in the ACM sheet, ",
         "e.g. ", paste0("\"", bad, "\"", collapse = ", "),
         ". The file layout may have changed.", call. = FALSE)
  }

  keep <- rep(TRUE, length(date))
  if (!is.null(end)) {
    end <- as.Date(end)
    if (length(end) != 1L || is.na(end)) {
      stop("`end` must be a single date.", call. = FALSE)
    }
    keep <- date <= end
  }
  if (!any(keep)) {
    stop("The ACM sheet has no observations on or before `end`.",
         call. = FALSE)
  }

  ord <- order(date[keep])
  wide <- wide[keep, , drop = FALSE][ord, , drop = FALSE]
  date <- date[keep][ord]

  # One date per row of the sheet, ten maturities per date: transposing the
  # date-by-tenor block reads it out date-major.
  pull <- function(prefix) {
    cols <- sprintf("%s%02d", prefix, tenors)
    block <- matrix(
      unlist(lapply(cols, function(col) {
        suppressWarnings(as.numeric(wide[[col]]))
      }), use.names = FALSE),
      nrow = nrow(wide)
    )
    as.vector(t(block))
  }

  out <- data.frame(
    date = rep(date, each = length(tenors)),
    maturity = rep(tenors * 12L, times = length(date)),
    fitted = pull(series[["fitted"]]),
    term_premium = pull(series[["term_premium"]]),
    risk_neutral = pull(series[["risk_neutral"]])
  )

  # If the components do not add up, the columns are not what they are
  # assumed to be, and nothing compared against them can be trusted.
  gap <- with(out, fitted - (risk_neutral + term_premium))
  gap <- abs(gap[is.finite(gap)])
  if (length(gap) && max(gap) > 1e-6) {
    stop("The published components do not add up: fitted differs from ",
         "risk_neutral + term_premium by up to ", format(max(gap), digits = 3),
         " percentage points. The columns may have been relabelled.",
         call. = FALSE)
  }

  out
}


#' Read the dates in the New York Fed's ACM sheet
#'
#' The sheet stores dates as text such as `"30-Jun-1961"`. Month names are
#' matched against `month.abb`, which is always English, rather than parsed
#' with `%b`, which follows the locale and returns `NA` for every row under a
#' non-English `LC_TIME`. Dates that are already dates, or ISO strings, are
#' accepted too, in case the layout changes.
#'
#' @param x The `DATE` column as read.
#' @return A `Date` vector, `NA` where a value cannot be read.
#' @keywords internal
#' @noRd
parse_acm_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x, tz = "UTC"))

  x <- trimws(as.character(x))
  out <- as.Date(x, format = "%Y-%m-%d")

  parts <- regmatches(x, regexec("^([0-9]{1,2})-([A-Za-z]{3})-([0-9]{4})$", x))
  dmy <- lengths(parts) == 4L
  if (any(dmy)) {
    p <- do.call(rbind, parts[dmy])
    month <- match(tolower(p[, 3L]), tolower(month.abb))
    out[dmy] <- as.Date(sprintf("%s-%02d-%02d", p[, 4L], month,
                                as.integer(p[, 2L])),
                        format = "%Y-%m-%d")
  }
  out
}
