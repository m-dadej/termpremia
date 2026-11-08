#' US Treasury zero-coupon yields, month-end, 1961-2024
#'
#' A frozen snapshot of the Gurkaynak, Sack and Wright (GSW) nominal Treasury
#' zero-coupon curve, evaluated on a monthly maturity grid and sampled at
#' month ends. This is the input used to reproduce the New York Fed's published
#' ACM term premium estimates, which [acm_published()] downloads.
#'
#' @section Why this is not simply the published file:
#' The Federal Reserve Board publishes fitted yields as `SVENY01`-`SVENY30`,
#' at **annual** steps starting at one year. The ACM model needs a monthly
#' maturity grid, so this dataset is produced by evaluating the published
#' Svensson parameters (`BETA0`-`BETA3`, `TAU1`, `TAU2`) at every maturity from
#' 1 to 120 months via [svensson_curve()]. The build script verifies that doing
#' so reproduces the published `SVENY` columns to within their rounding
#' precision (5e-05) at the 1, 2, 5 and 10 year tenors.
#'
#' @section Extrapolation:
#' The Svensson parameters define the curve at every maturity, including
#' maturities longer than any bond actually trading on that date. Early in the
#' sample, ten-year yields are therefore *extrapolations*, not observations:
#' `SVENY10` is absent for 122 of the 763 month-ends. The `extrapolated` column
#' flags these, and about 4.8% of rows are affected. Treat them as model output,
#' and consider excluding them when the distinction matters.
#'
#' @format A data frame with 91,560 rows and 4 columns:
#' \describe{
#'   \item{date}{Month-end observation date. The last date in each calendar
#'     month for which a curve was fitted, which is not always the last
#'     calendar day.}
#'   \item{maturity}{Maturity in months, 1 to 120.}
#'   \item{yield}{Continuously compounded zero-coupon yield, **in percent**.
#'     Pass `units = "percent"` to [yield_panel()].}
#'   \item{extrapolated}{`TRUE` where the maturity exceeds the longest maturity
#'     the Board published for that date.}
#' }
#'
#' @source
#' Federal Reserve Board series feds200628,
#' <https://www.federalreserve.gov/data/nominal-yield-curve.htm>. A staff
#' research product, not an official statistical release, and subject to
#' revision without notice. Built by `data-raw/gsw.R`.
#'
#' @references
#' Gurkaynak, R. S., B. Sack and J. H. Wright (2007). "The U.S. Treasury yield
#' curve: 1961 to the present." *Journal of Monetary Economics* 54(8),
#' 2291-2304.
#'
#' @examples
#' p <- yield_panel(gsw_monthly, units = "percent", maturity_unit = "months",
#'                  instrument = "government", issuer = "US")
#' p
#'
#' # Excluding extrapolated long maturities early in the sample
#' observed <- gsw_monthly[!gsw_monthly$extrapolated, ]
#' range(observed$date[observed$maturity == 120])
"gsw_monthly"


#' Three-month US Treasury bill rate, month-end aligned, 1961-2024
#'
#' The 3-month Treasury bill secondary market rate, as an alternative short rate
#' for [atsm()].
#'
#' @section Why an alternative short rate matters:
#' Gurkaynak, Sack and Wright exclude Treasury bills, and every security with
#' under three months to maturity, when fitting their curve. Their short end is
#' therefore an extrapolation, and they caution against it. Measured against
#' this series, the GSW fitted one-month yield sits about **28bp higher on
#' average**, with a standard deviation of 127bp.
#'
#' That matters because the short rate enters excess returns directly as
#' \eqn{-r_t}, so it moves the estimated prices of risk and hence the level of
#' the term premium. Substituting this series for the fitted one-month yield
#' shifts the estimated 10-year US term premium by roughly 27bp. See
#' `vignette("validation")`.
#'
#' A caveat: this is a *three*-month rate standing in for the model's
#' *one*-month period, which is internally inconsistent and degrades the fit to
#' the curve. It is a diagnostic for bounding the short rate's influence, not a
#' drop-in improvement.
#'
#' @format A data frame with 763 rows and 2 columns:
#' \describe{
#'   \item{date}{Month-end date, aligned to [gsw_monthly].}
#'   \item{value}{Secondary market rate, monthly average, in percent.}
#' }
#'
#' @source Federal Reserve Board H.15,
#'   <https://www.federalreserve.gov/releases/h15/>, as distributed by FRED
#'   (series `TB3MS`). Built by `data-raw/tbill.R`.
#'
#' @examples
#' panel <- yield_panel(gsw_monthly, units = "percent",
#'                      maturity_unit = "months", issuer = "US")
#'
#' # The fitted short end runs well above actual bill rates
#' fitted_1m <- gsw_monthly$yield[gsw_monthly$maturity == 1]
#' mean(fitted_1m - tbill_3m$value)
"tbill_3m"
