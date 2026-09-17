#' Model-free term premium from survey expectations
#'
#' The term premium without a model: the observed yield minus a survey's
#' expected average short rate over the same horizon.
#'
#' \deqn{\textrm{term premium}^{(n)}_t = y^{(n)}_t - \textrm{survey expected
#' average short rate over } n}
#'
#' This is crude -- surveys are infrequent, cover few horizons, and forecasters
#' may not represent the marginal investor -- but it depends on no model at all,
#' which makes it the natural sanity check on a fitted decomposition. It is what
#' exposed a roughly 100 basis point discrepancy in published euro area ACM
#' estimates in the BIS's 2018 review.
#'
#' @section Frequency, and why nothing is interpolated by default:
#' The relevant US series, the Philadelphia Fed's `BILL10` (expected annual
#' average return on 3-month Treasury bills over the current and next nine
#' years), is collected in **first-quarter surveys only**. It is an annual
#' series. A monthly model-free term premium at the 10-year horizon cannot be
#' constructed from it, which is why the BIS plot it as points rather than a
#' line.
#'
#' This function therefore returns one row per matched *survey* date and never
#' silently fills the gaps. Set `interpolate = TRUE` only if you have a specific
#' reason, and know that the result between survey dates is invention, not data.
#'
#' @section Data sources:
#' None is bundled. The Philadelphia Fed asserts copyright with all rights
#' reserved and limits use to research purposes, which is incompatible with
#' redistribution under this package's licence. Free options include Philadelphia
#' Fed SPF `BILL10` and `TBILL`, the FOMC Summary of Economic Projections
#' longer-run federal funds rate, and the New York Fed's Survey of Primary
#' Dealers. For the euro area the BIS use Consensus Economics, which is
#' proprietary.
#'
#' @param panel A [yield_panel].
#' @param survey A data frame with a `date` column and a `value` column giving
#'   the expected average short rate over `maturity` months.
#' @param maturity Horizon in months. Must be on the panel's grid and must match
#'   the horizon the survey actually asks about.
#' @param curve Curve to use. Defaults to the first.
#' @param survey_units `"auto"`, `"percent"` or `"decimal"`, interpreted as for
#'   [yield_panel()].
#' @param max_gap Maximum number of days between a survey date and the panel
#'   observation matched to it. Survey dates with no observation inside this
#'   window are dropped and reported.
#' @param interpolate If `TRUE`, linearly interpolate the survey onto every
#'   panel date. Off by default, and reported loudly when used.
#'
#' @return A data frame with `date` (the panel observation date), `survey_date`,
#'   `maturity`, `yield`, `expected_short_rate` and `term_premium`. Yields are
#'   annualised decimals.
#'
#' @examples
#' panel <- yield_panel(gsw_monthly, units = "percent",
#'                      maturity_unit = "months", issuer = "US")
#'
#' # FABRICATED survey data, for illustration only -- this is not a forecast of
#' # anything. Real series are named under "Data sources" above.
#' fake_survey <- data.frame(
#'   date = seq(as.Date("2000-02-15"), as.Date("2024-02-15"), by = "year"),
#'   value = 3.5 + sin(seq(0, 6, length.out = 25)) / 2
#' )
#'
#' tp <- term_premium_survey(panel, fake_survey, maturity = 120,
#'                           survey_units = "percent")
#' head(tp)
#'
#' @export
term_premium_survey <- function(panel,
                                survey,
                                maturity = 120L,
                                curve = NULL,
                                survey_units = c("auto", "percent", "decimal"),
                                max_gap = 45L,
                                interpolate = FALSE) {
  survey_units <- match.arg(survey_units)

  if (!inherits(panel, "yield_panel")) {
    stop("`panel` must be a yield_panel; see ?yield_panel.", call. = FALSE)
  }
  if (!is.data.frame(survey)) {
    stop("`survey` must be a data frame with `date` and `value` columns.",
         call. = FALSE)
  }
  missing_cols <- setdiff(c("date", "value"), names(survey))
  if (length(missing_cols)) {
    stop("`survey` is missing column(s): ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  curve <- curve %||% panel$meta$curve[1L]
  grid <- panel$maturities[[curve]]
  if (!maturity %in% grid) {
    stop("Maturity ", maturity, " is not on the panel's grid.", call. = FALSE)
  }

  y <- curve_matrix(panel, curve)[, match(maturity, grid)]

  s_date <- as_date_strict(survey$date, arg = "survey$date")
  s_val <- as_yield_decimal(survey$value, survey_units, arg = "survey$value")

  keep <- !is.na(s_val)
  s_date <- s_date[keep]
  s_val <- s_val[keep]

  ord <- order(s_date)
  s_date <- s_date[ord]
  s_val <- s_val[ord]

  if (interpolate) {
    return(survey_interpolated(panel$dates, y, s_date, s_val, maturity))
  }

  # Match each survey date to the nearest panel observation. Surveys are run on
  # their own calendar and will rarely land exactly on a month end.
  idx <- vapply(s_date, function(d) {
    gaps <- abs(as.numeric(panel$dates - d))
    if (min(gaps) > max_gap) NA_integer_ else which.min(gaps)
  }, integer(1))

  unmatched <- sum(is.na(idx))
  if (unmatched) {
    message(
      unmatched, " of ", length(idx), " survey date(s) had no panel ",
      "observation within ", max_gap, " days and were dropped."
    )
  }
  if (all(is.na(idx))) {
    stop("No survey date could be matched to a panel observation. Check that ",
         "the survey and panel cover overlapping periods.", call. = FALSE)
  }

  ok <- !is.na(idx)
  out <- data.frame(
    date = panel$dates[idx[ok]],
    survey_date = s_date[ok],
    maturity = maturity,
    yield = unname(y[idx[ok]]),
    expected_short_rate = s_val[ok],
    stringsAsFactors = FALSE
  )
  out$term_premium <- out$yield - out$expected_short_rate
  out
}

#' @noRd
survey_interpolated <- function(dates, y, s_date, s_val, maturity) {
  if (length(s_date) < 2L) {
    stop("At least two survey observations are needed to interpolate.",
         call. = FALSE)
  }

  message(
    "Interpolating ", length(s_date), " survey observations onto ",
    length(dates), " panel dates. Values between survey dates are ",
    "constructed, not observed."
  )

  inside <- dates >= min(s_date) & dates <= max(s_date)
  est <- stats::approx(s_date, s_val, xout = dates[inside])$y

  out <- data.frame(
    date = dates[inside],
    survey_date = as.Date(NA),
    maturity = maturity,
    yield = unname(y[inside]),
    expected_short_rate = est,
    stringsAsFactors = FALSE
  )
  out$term_premium <- out$yield - out$expected_short_rate
  out
}
