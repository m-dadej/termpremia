#' Nelson-Siegel and Svensson zero-coupon yield curves
#'
#' Evaluate the Nelson-Siegel (1987) or Svensson (1994) parametric yield curve at
#' arbitrary maturities. These are the parameterisations published by
#' Gurkaynak, Sack and Wright (2007) for the US Treasury curve, and by several
#' other central banks for their own markets.
#'
#' The Svensson zero-coupon yield at maturity `n` (in years) is
#'
#' \deqn{y(n) = \beta_0 + \beta_1 f(n, \tau_1) + \beta_2 [f(n, \tau_1) -
#'   e^{-n/\tau_1}] + \beta_3 [f(n, \tau_2) - e^{-n/\tau_2}]}
#'
#' where \eqn{f(n, \tau) = (1 - e^{-n/\tau}) / (n/\tau)}. Setting
#' \eqn{\beta_3 = 0} recovers Nelson-Siegel.
#'
#' Before 1980 Gurkaynak, Sack and Wright fitted only the Nelson-Siegel form,
#' and encode this with `beta3 = 0` alongside a **sentinel** `tau2` of
#' `-999.99` (not `0`, and not `NA`). Because the fourth term is skipped
#' entirely whenever `beta3` is zero, such rows evaluate correctly rather than
#' producing `NaN`. A non-zero `beta3` combined with a non-positive `tau2` is an
#' error, not silently ignored: that combination means the inputs are wrong.
#'
#' The yields returned are continuously compounded and carry the same units as
#' the supplied parameters. Gurkaynak, Sack and Wright publish parameters in
#' percent, so the result is in percent unless the inputs were rescaled.
#'
#' @param maturity Numeric vector of maturities, **in years**. Must be
#'   non-negative and finite.
#' @param beta0,beta1,beta2,beta3 Level, slope and curvature parameters.
#'   `beta3` defaults to `0`, giving Nelson-Siegel.
#' @param tau1,tau2 Decay parameters, in years. `tau1` must be positive.
#'   `tau2` may be `0` or `NA` when `beta3` is zero.
#'
#' @return A numeric vector of zero-coupon yields, the same length as
#'   `maturity`.
#'
#' @references
#' Gurkaynak, R. S., B. Sack and J. H. Wright (2007). "The U.S. Treasury yield
#' curve: 1961 to the present." *Journal of Monetary Economics* 54(8),
#' 2291-2304.
#'
#' Svensson, L. E. O. (1994). "Estimating and interpreting forward interest
#' rates: Sweden 1992-1994." NBER Working Paper 4871.
#'
#' @examples
#' # A upward-sloping curve, maturities of 1 to 10 years
#' svensson_yield(1:10, beta0 = 4, beta1 = -2, beta2 = 1, beta3 = 0.5,
#'                tau1 = 1.5, tau2 = 10)
#'
#' # Nelson-Siegel: beta3 and tau2 omitted
#' svensson_yield(1:10, beta0 = 4, beta1 = -2, beta2 = 1, tau1 = 1.5)
#'
#' @export
svensson_yield <- function(maturity, beta0, beta1, beta2, beta3 = 0,
                           tau1, tau2 = NA_real_) {
  maturity <- as.numeric(maturity)

  if (anyNA(maturity)) {
    stop("`maturity` must not contain NA.", call. = FALSE)
  }
  if (any(maturity < 0)) {
    stop("`maturity` must be non-negative.", call. = FALSE)
  }
  if (any(!is.finite(maturity))) {
    stop("`maturity` must be finite.", call. = FALSE)
  }

  scalars <- list(beta0 = beta0, beta1 = beta1, beta2 = beta2,
                  beta3 = beta3, tau1 = tau1)
  for (nm in names(scalars)) {
    if (length(scalars[[nm]]) != 1L) {
      stop("`", nm, "` must be a single value.", call. = FALSE)
    }
  }

  # A missing parameter set (a date with no fitted curve) yields all-NA rather
  # than an error: GSW files contain such rows.
  if (anyNA(c(beta0, beta1, beta2, tau1))) {
    return(rep(NA_real_, length(maturity)))
  }

  if (!is.finite(tau1) || tau1 <= 0) {
    stop("`tau1` must be a positive, finite number.", call. = FALSE)
  }

  beta3 <- if (is.na(beta3)) 0 else beta3

  # Nelson-Siegel is Svensson with the fourth term switched off. GSW encode this
  # as beta3 = 0 together with tau2 = 0, which would otherwise divide by zero.
  use_second_hump <- beta3 != 0
  if (use_second_hump && (is.na(tau2) || !is.finite(tau2) || tau2 <= 0)) {
    stop(
      "`tau2` must be a positive, finite number when `beta3` is non-zero.",
      call. = FALSE
    )
  }

  out <- beta0 +
    beta1 * ns_level_factor(maturity, tau1) +
    beta2 * ns_hump_factor(maturity, tau1)

  if (use_second_hump) {
    out <- out + beta3 * ns_hump_factor(maturity, tau2)
  }

  out
}

#' Nelson-Siegel basis functions
#'
#' `ns_level_factor()` is \eqn{(1 - e^{-n/\tau}) / (n/\tau)} and
#' `ns_hump_factor()` subtracts \eqn{e^{-n/\tau}} from it. Both are continuous
#' at `n = 0`, where the level factor tends to 1 and the hump factor to 0; the
#' limits are imposed explicitly because the expressions are 0/0 there.
#'
#' @param maturity Numeric vector of maturities in years.
#' @param tau Positive decay parameter in years.
#' @return A numeric vector the same length as `maturity`.
#' @noRd
ns_level_factor <- function(maturity, tau) {
  x <- maturity / tau
  out <- (1 - exp(-x)) / x

  # The expression is 0/0 at x = 0 and loses precision catastrophically for
  # very small x, where the limit and the second-order expansion are exact to
  # machine precision anyway.
  tiny <- x < 1e-8
  out[tiny] <- 1 - x[tiny] / 2

  out
}

#' @rdname ns_level_factor
#' @noRd
ns_hump_factor <- function(maturity, tau) {
  ns_level_factor(maturity, tau) - exp(-maturity / tau)
}

#' Evaluate a panel of Nelson-Siegel or Svensson curves
#'
#' Vectorised counterpart to [svensson_yield()]: evaluates one curve per row of
#' `params` on a shared maturity grid, which is how a published parameter file
#' such as Gurkaynak, Sack and Wright's is turned into a yield panel.
#'
#' @param params A data frame or matrix with one row per date and columns
#'   `beta0`, `beta1`, `beta2`, `beta3`, `tau1`, `tau2`. Column names are
#'   matched case-insensitively. Rows whose parameters are missing produce a row
#'   of `NA`.
#' @param maturity Numeric vector of maturities, **in years**.
#'
#' @return A numeric matrix with `nrow(params)` rows and `length(maturity)`
#'   columns. Row names are taken from `params` if present; column names are the
#'   maturities.
#'
#' @examples
#' params <- data.frame(
#'   beta0 = c(4, 4.1), beta1 = c(-2, -1.9), beta2 = c(1, 1.1),
#'   beta3 = c(0.5, 0.4), tau1 = c(1.5, 1.5), tau2 = c(10, 10)
#' )
#' svensson_curve(params, maturity = c(1, 5, 10))
#'
#' @export
svensson_curve <- function(params, maturity) {
  params <- as.data.frame(params)
  needed <- c("beta0", "beta1", "beta2", "beta3", "tau1", "tau2")

  idx <- match(needed, tolower(names(params)))
  if (anyNA(idx)) {
    stop(
      "`params` is missing required column(s): ",
      paste(needed[is.na(idx)], collapse = ", "),
      call. = FALSE
    )
  }
  params <- params[, idx, drop = FALSE]
  names(params) <- needed

  out <- vapply(
    seq_len(nrow(params)),
    function(i) {
      svensson_yield(
        maturity = maturity,
        beta0 = params$beta0[i], beta1 = params$beta1[i],
        beta2 = params$beta2[i], beta3 = params$beta3[i],
        tau1  = params$tau1[i],  tau2  = params$tau2[i]
      )
    },
    numeric(length(maturity))
  )

  out <- if (length(maturity) == 1L) matrix(out, ncol = 1L) else t(out)
  dimnames(out) <- list(rownames(params), format(maturity))
  out
}
