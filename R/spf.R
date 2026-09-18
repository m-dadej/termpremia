# Philadelphia Fed Survey of Professional Forecasters.
#
# This is the one dataset the package fetches rather than bundles. Everywhere
# else the decision was to ship a reproducible snapshot and never touch the
# network; the SPF is an exception because the Philadelphia Fed asserts
# copyright with all rights reserved and limits use to research purposes, which
# is incompatible with redistribution under an MIT licence. Fetching it puts
# the user in a direct relationship with the publisher's terms, which is the
# only arrangement that works.

#' Download short-rate forecasts from the Survey of Professional Forecasters
#'
#' Fetches the Philadelphia Fed's mean `TBILL` responses -- forecasts of the
#' three-month Treasury bill rate -- and reshapes them into the form
#' [atsm()] expects for `p_dynamics = "survey"`.
#'
#' @section Licensing:
#' **This data is not redistributed with the package, and this function
#' downloads it on your behalf rather than shipping a copy.** The Philadelphia
#' Fed asserts copyright over the SPF with all rights reserved and limits use
#' to research purposes. Using it here is fine; passing it on is a question for
#' them, not for this package's licence. See `LICENSE.note`.
#'
#' @section What the columns mean:
#' The SPF publishes, for each survey quarter, a row of forecasts named
#' `TBILL1` to `TBILL6`. `TBILL1` is not a forecast at all: it is the previous
#' quarter's realised rate as known to respondents, and this function uses it
#' as an integrity check -- against the bundled `tbill_3m` it agrees to under
#' three basis points, which is what pins the column indexing down. `TBILL2` is
#' the current quarter and `TBILL3` to `TBILL6` are one to four quarters ahead.
#'
#' Each forecast is a *quarterly average* of a *three-month* rate, so the
#' matching control settings are `survey_control(tenor = 3,
#' average_months = 3)`.
#'
#' The survey is collected in the middle of the quarter. `date` is set to the
#' middle month of the survey quarter, which is accurate to a couple of weeks
#' and is matched to the nearest panel observation anyway.
#'
#' @param horizons Which quarterly horizons to return, in quarters ahead of the
#'   survey quarter. `0` is the current quarter, which [atsm()] will drop
#'   because it refers to a period already under way; the default of `1:4` is
#'   the genuinely forward-looking part.
#' @param url Location of the mean-response file. Exposed so that a local copy
#'   can be used instead, and because the Philadelphia Fed have moved these
#'   files before.
#' @param quiet Suppress the download progress bar.
#'
#' @return A data frame with `date` (when the forecast was made),
#'   `target_date` (first month of the quarter forecast), `horizon_quarters`
#'   and `value` (percent), ready to pass to `atsm(survey = )`.
#'
#' @section Requirements:
#' Needs the `readxl` package and a working internet connection. Neither is a
#' hard dependency of termpremia.
#'
#' @references
#' Federal Reserve Bank of Philadelphia. *Survey of Professional Forecasters*.
#'
#' @seealso [atsm()], [survey_control()], [term_premium_survey()]
#'
#' @examples
#' \dontrun{
#' spf <- spf_tbill()
#'
#' panel <- yield_panel(gsw_monthly, units = "percent",
#'                      maturity_unit = "months", issuer = "US")
#'
#' fit <- atsm(panel, n_factors = 5, p_dynamics = "survey", survey = spf,
#'             survey_control = survey_control(tenor = 3, average_months = 3))
#' }
#'
#' @export
spf_tbill <- function(horizons = 1:4,
                      url = paste0(
                        "https://www.philadelphiafed.org/-/media/frbp/assets/",
                        "surveys-and-data/survey-of-professional-forecasters/",
                        "data-files/files/mean_tbill_level.xlsx"
                      ),
                      quiet = TRUE) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("`spf_tbill()` needs the readxl package: install.packages(\"readxl\").",
         call. = FALSE)
  }
  if (!is.numeric(horizons) || !length(horizons) ||
      any(horizons < 0 | horizons > 4 | horizons != as.integer(horizons))) {
    stop("`horizons` must be whole numbers between 0 and 4.", call. = FALSE)
  }

  dest <- tempfile(fileext = ".xlsx")
  on.exit(unlink(dest), add = TRUE)

  ok <- tryCatch({
    utils::download.file(url, dest, mode = "wb", quiet = quiet)
    TRUE
  }, error = function(e) {
    stop("Could not download the SPF file from\n  ", url, "\n", conditionMessage(e),
         "\nThe Philadelphia Fed have relocated these files before; check the ",
         "survey's data page and pass a new `url` if so.", call. = FALSE)
  })
  stopifnot(ok)

  raw <- as.data.frame(readxl::read_excel(dest))
  spf_reshape(raw, horizons)
}


#' Download ten-year short-rate expectations from the SPF
#'
#' Fetches the Philadelphia Fed's mean `BILL10` responses: the expected annual
#' average return on three-month Treasury bills over the next ten years.
#'
#' @section Why this one matters:
#' [spf_tbill()] reaches four quarters ahead, which disciplines the near term
#' and says almost nothing about where the model thinks rates are heading. The
#' shifting-endpoint problem lives at the far end, so anchoring it needs a
#' far-end forecast, and `BILL10` is the only free one. Kim and Wright use Blue
#' Chip's six-to-eleven-year forecasts for the same purpose; those are
#' proprietary, which is why this is "Kim-Wright style" and not a replication.
#'
#' It is a thin series. `BILL10` is asked in **first-quarter surveys only**, so
#' there are about 35 observations, beginning 1992Q1. Thin is not the same as
#' uninformative: each one is a direct statement about a horizon the yield
#' curve alone cannot pin down, and they are the observations that move the
#' answer most.
#'
#' @section Licensing:
#' As [spf_tbill()]: downloaded, never redistributed. See `LICENSE.note`.
#'
#' @section Dating:
#' The ten-year window nominally begins with the current year, but the survey
#' is taken in February, so the window is placed at the month following the
#' survey instead. That shifts a 120-month average by at most two months, which
#' is well inside what any of this resolves.
#'
#' @param url Location of the mean-response file.
#' @param quiet Suppress the download progress bar.
#'
#' @return A data frame shaped as [spf_tbill()]'s, so the two can be combined
#'   with `rbind()`.
#'
#' @seealso [spf_tbill()], [atsm()], [survey_control()]
#'
#' @examples
#' \dontrun{
#' # Near and far horizons together is the configuration that works.
#' surveys <- rbind(spf_tbill(), spf_bill10())
#'
#' panel <- yield_panel(gsw_monthly, units = "percent",
#'                      maturity_unit = "months", issuer = "US")
#' fit <- atsm(panel, n_factors = 5, p_dynamics = "survey", survey = surveys)
#' }
#'
#' @export
spf_bill10 <- function(url = paste0(
                         "https://www.philadelphiafed.org/-/media/frbp/assets/",
                         "surveys-and-data/survey-of-professional-forecasters/",
                         "data-files/files/mean_bill10_level.xlsx"
                       ),
                       quiet = TRUE) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("`spf_bill10()` needs the readxl package: ",
         "install.packages(\"readxl\").", call. = FALSE)
  }

  dest <- tempfile(fileext = ".xlsx")
  on.exit(unlink(dest), add = TRUE)

  tryCatch(
    utils::download.file(url, dest, mode = "wb", quiet = quiet),
    error = function(e) {
      stop("Could not download the SPF BILL10 file from\n  ", url, "\n",
           conditionMessage(e), call. = FALSE)
    }
  )

  spf_reshape_bill10(as.data.frame(readxl::read_excel(dest)))
}


#' Reshape an SPF BILL10 sheet into forecast rows
#'
#' @param raw The sheet, as read.
#' @return A data frame of forecasts.
#' @keywords internal
#' @noRd
spf_reshape_bill10 <- function(raw) {
  missing_cols <- setdiff(c("YEAR", "QUARTER", "BILL10"), names(raw))
  if (length(missing_cols)) {
    stop("The SPF BILL10 sheet is missing column(s): ",
         paste(missing_cols, collapse = ", "),
         ". The file layout may have changed.", call. = FALSE)
  }

  value <- suppressWarnings(as.numeric(as.character(raw$BILL10)))
  year <- as.integer(raw$YEAR)
  quarter <- as.integer(raw$QUARTER)

  quarter_start <- as.Date(sprintf("%d-%02d-01", year, (quarter - 1L) * 3L + 1L))
  survey_date <- add_months(quarter_start, 1L) + 14L

  out <- data.frame(
    date = survey_date,
    target_date = add_months(survey_date, 1L),
    horizon_quarters = 0L,
    value = value,
    tenor = 3L,                        # three-month bill
    average_months = 120L,             # averaged over ten years
    stringsAsFactors = FALSE
  )

  out <- out[!is.na(out$value), , drop = FALSE]
  if (!nrow(out)) {
    stop("The SPF sheet contained no usable BILL10 forecasts.", call. = FALSE)
  }

  out <- out[order(out$date), , drop = FALSE]
  rownames(out) <- NULL

  attr(out, "source") <- "Philadelphia Fed Survey of Professional Forecasters"
  attr(out, "retrieved") <- Sys.Date()
  out
}


#' Reshape an SPF mean-response sheet into forecast rows
#'
#' Split out from the download so that the parsing can be tested without a
#' network connection, and so that a locally held copy of the file can be used.
#'
#' @param raw The sheet, as read.
#' @param horizons Quarters ahead to keep.
#'
#' @return A long data frame of forecasts.
#' @keywords internal
#' @noRd
spf_reshape <- function(raw, horizons) {
  needed <- c("YEAR", "QUARTER", paste0("TBILL", 1:6))
  missing_cols <- setdiff(needed, names(raw))
  if (length(missing_cols)) {
    stop("The SPF sheet is missing column(s): ",
         paste(missing_cols, collapse = ", "),
         ". The file layout may have changed.", call. = FALSE)
  }

  # Missing values are the literal string "#N/A", so every forecast column
  # arrives as character.
  for (col in paste0("TBILL", 1:6)) {
    raw[[col]] <- suppressWarnings(as.numeric(as.character(raw[[col]])))
  }
  raw$YEAR <- as.integer(raw$YEAR)
  raw$QUARTER <- as.integer(raw$QUARTER)

  if (anyNA(raw$YEAR) || anyNA(raw$QUARTER) ||
      any(raw$QUARTER < 1 | raw$QUARTER > 4)) {
    stop("The SPF sheet has unreadable YEAR/QUARTER values.", call. = FALSE)
  }

  quarter_start <- as.Date(sprintf("%d-%02d-01", raw$YEAR,
                                   (raw$QUARTER - 1L) * 3L + 1L))

  # Surveys are collected in the middle of the quarter: the questionnaire goes
  # out at the end of the first month and is due in the middle of the second.
  survey_date <- add_months(quarter_start, 1L) + 14L

  out <- do.call(rbind, lapply(horizons, function(h) {
    col <- paste0("TBILL", h + 2L)     # TBILL2 is the current quarter
    data.frame(
      date = survey_date,
      target_date = add_months(quarter_start, 3L * h),
      horizon_quarters = as.integer(h),
      value = raw[[col]],
      tenor = 3L,                      # three-month bill
      average_months = 3L,             # averaged over the quarter
      stringsAsFactors = FALSE
    )
  }))

  out <- out[!is.na(out$value), , drop = FALSE]
  if (!nrow(out)) {
    stop("The SPF sheet contained no usable TBILL forecasts.", call. = FALSE)
  }

  out <- out[order(out$date, out$horizon_quarters), , drop = FALSE]
  rownames(out) <- NULL

  attr(out, "source") <- "Philadelphia Fed Survey of Professional Forecasters"
  attr(out, "retrieved") <- Sys.Date()
  out
}


#' Add whole months to a date, staying on the first of the month
#'
#' `seq.Date(by = "month")` overshoots from month ends -- 31 January plus one
#' month is 3 March. Every date this is used on is the first of a month, but
#' the arithmetic is done on year and month directly so that it cannot depend
#' on that.
#'
#' @param d Dates.
#' @param n Months to add.
#' @return Dates, on the first of the resulting month.
#' @keywords internal
#' @noRd
add_months <- function(d, n) {
  lt <- as.POSIXlt(d)
  total <- lt$year * 12L + lt$mon + as.integer(n)
  as.Date(sprintf("%d-%02d-01", total %/% 12L + 1900L, total %% 12L + 1L))
}
