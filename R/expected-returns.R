# Expected one-period excess returns.
#
# The quantity the estimator is actually built around. ACM's step 2 regresses
# realised excess returns on innovations and lagged factors, and step 3 turns
# the result into prices of risk; the expected excess return is that machinery
# read forwards rather than backwards. It is also the cleanest way to see what
# the prices of risk are claiming: a term premium is an accumulated expected
# excess return, so if the premium looks implausible this is where to look.

#' Expected one-period excess bond returns
#'
#' The return an `n`-month bond is expected to earn over the next month in
#' excess of the one-month rate, given the factors today.
#'
#' @section The formula:
#' \deqn{E_t\left[rx^{(n-1)}_{t+1}\right] =
#'   B_{n-1}'(\lambda_0 + \lambda_1 X_t)
#'   - \tfrac12\left(B_{n-1}'\Sigma B_{n-1} + \sigma^2\right)}
#'
#' The first term is compensation for risk: the bond's exposure to the factors
#' priced by the prices of risk. The second is a Jensen term, the convexity
#' adjustment that makes a log return differ from the log of an expected
#' return; `expected_excess_return()` can return it separately.
#'
#' Note the loading is \eqn{B_{n-1}}, not \eqn{B_n}. A bond bought with `n`
#' months to run is sold a month later with `n - 1`, so it is the shorter
#' bond's sensitivity that determines the risk taken. A one-month bond held
#' one month is riskless and its expected excess return is exactly zero, which
#' the recursion delivers through \eqn{B_0 = 0} and which is asserted in the
#' tests.
#'
#' @section Loadings:
#' `expected_excess_return_loadings()` gives \eqn{B_{n-1}'\lambda_1}, scaled
#' by each factor's standard deviation, so an entry reads as the effect on the
#' expected excess return of a one standard deviation move in that factor. It
#' is the compact answer to "which factor is driving the premium".
#'
#' @param object An `atsm_fit` or `atsm_real_fit`.
#' @param maturity Maturities in months. Defaults to every fitted maturity.
#' @param convexity When `TRUE`, the returned frame carries the risk
#'   compensation and the convexity adjustment as separate components as well
#'   as their sum.
#' @param curve For a joint fit, `"nominal"` or `"real"`.
#' @param ... Unused.
#'
#' @return `expected_excess_return()` returns a data frame with `date`,
#'   `maturity` and `value`, in monthly decimal units -- these are one-month
#'   returns, not annualised rates, so multiplying by 12 would misstate them.
#'   With `convexity = TRUE` there is also a `component` column.
#'   `expected_excess_return_loadings()` returns a data frame with `maturity`,
#'   `factor` and `value`.
#'
#' @seealso [term_premium()], [forward_rate()]
#'
#' @examples
#' panel <- yield_panel(gsw_monthly, units = "percent",
#'                      maturity_unit = "months", issuer = "US")
#' fit <- atsm(panel, n_factors = 5)
#'
#' er <- expected_excess_return(fit, maturity = 120)
#' mean(er$value) * 1e4        # basis points per month
#'
#' expected_excess_return_loadings(fit, maturity = c(24, 120))
#'
#' @name expected-excess-returns
NULL

#' @rdname expected-excess-returns
#' @export
expected_excess_return <- function(object, maturity = NULL,
                                   convexity = FALSE, ...) {
  UseMethod("expected_excess_return")
}

#' @export
expected_excess_return.atsm_fit <- function(object, maturity = NULL,
                                            convexity = FALSE, ...) {
  er_frame(object, object$maturities, object$recursion$p$b,
           object$pars, maturity, convexity, object$pars$sigma2)
}

#' @rdname expected-excess-returns
#' @export
expected_excess_return.atsm_real_fit <- function(object, maturity = NULL,
                                                 convexity = FALSE,
                                                 curve = c("nominal", "real"),
                                                 ...) {
  curve <- match.arg(curve)
  grid <- if (curve == "nominal") object$maturities else object$real_maturities
  b <- if (curve == "nominal") {
    object$recursion$q$b
  } else {
    # The real curve's risk exposure runs through B + pi1: one period of
    # inflation enters the payoff, so it is carried into the next step.
    sweep(object$recursion$q$b_real, 2L, object$pars$pi1, "+")
  }
  # ACMY put measurement error on returns in the likelihood rather than a
  # sigma2 term in the recursion, so there is none to add here.
  er_frame(object, grid, b, object$pars, maturity, convexity, 0)
}

#' @keywords internal
#' @noRd
er_frame <- function(object, grid, b, pars, maturity, convexity, sigma2) {
  maturity <- maturity %||% grid
  bad <- setdiff(maturity, grid)
  if (length(bad)) {
    stop("Maturity(ies) not fitted: ", paste(bad, collapse = ", "),
         ". Available: ", min(grid), "-", max(grid), " months.",
         call. = FALSE)
  }

  parts <- er_parts(maturity, b, pars, object$factors, sigma2)
  n <- length(object$dates)

  if (!convexity) {
    return(data.frame(
      date = rep(object$dates, times = length(maturity)),
      maturity = rep(maturity, each = n),
      value = as.vector(parts$risk + parts$convexity),
      stringsAsFactors = FALSE
    ))
  }

  pieces <- list(expected = parts$risk + parts$convexity,
                 risk_compensation = parts$risk,
                 convexity = parts$convexity)
  data.frame(
    date = rep(object$dates, times = length(maturity) * length(pieces)),
    maturity = rep(rep(maturity, each = n), times = length(pieces)),
    component = rep(names(pieces), each = n * length(maturity)),
    value = unlist(lapply(pieces, as.vector), use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

#' Risk compensation and convexity, by maturity and date
#'
#' @param maturity Maturities in months.
#' @param b Loading matrix indexed by maturity.
#' @param pars Parameter list.
#' @param x Factor matrix.
#' @param sigma2 Return pricing error variance, or 0.
#' @return A list of two `T x length(maturity)` matrices.
#' @keywords internal
#' @noRd
er_parts <- function(maturity, b, pars, x, sigma2) {
  n_t <- nrow(x)
  risk <- matrix(0, n_t, length(maturity))
  conv <- matrix(0, n_t, length(maturity))

  # lambda_t = lambda0 + lambda1 X_t, one row per date.
  price <- sweep(x %*% t(pars$lambda1), 2L, pars$lambda0, "+")

  for (j in seq_along(maturity)) {
    n <- maturity[j]
    if (n <= 1L) next          # a one-month bond held one month is riskless

    bp <- b[n - 1L, ]
    risk[, j] <- drop(price %*% bp)
    conv[, j] <- -0.5 * (drop(crossprod(bp, pars$sigma %*% bp)) + sigma2)
  }

  dimnames(risk) <- dimnames(conv) <-
    list(rownames(x), as.character(maturity))
  list(risk = risk, convexity = conv)
}

#' @rdname expected-excess-returns
#' @export
expected_excess_return_loadings <- function(object, maturity = NULL, ...) {
  UseMethod("expected_excess_return_loadings")
}

#' @export
expected_excess_return_loadings.atsm_fit <- function(object, maturity = NULL,
                                                     ...) {
  er_loadings(object, object$maturities, object$recursion$p$b, maturity)
}

#' @export
expected_excess_return_loadings.atsm_real_fit <- function(object,
                                                          maturity = NULL,
                                                          curve = c("nominal",
                                                                    "real"),
                                                          ...) {
  curve <- match.arg(curve)
  grid <- if (curve == "nominal") object$maturities else object$real_maturities
  b <- if (curve == "nominal") {
    object$recursion$q$b
  } else {
    sweep(object$recursion$q$b_real, 2L, object$pars$pi1, "+")
  }
  er_loadings(object, grid, b, maturity)
}

#' @keywords internal
#' @noRd
er_loadings <- function(object, grid, b, maturity) {
  maturity <- maturity %||% grid
  bad <- setdiff(maturity, grid)
  if (length(bad)) {
    stop("Maturity(ies) not fitted: ", paste(bad, collapse = ", "),
         ". Available: ", min(grid), "-", max(grid), " months.",
         call. = FALSE)
  }

  # Scaled by each factor's own standard deviation, so the entries are
  # comparable across factors and read as "effect of a one sd move".
  sd_x <- apply(object$factors, 2L, stats::sd)
  nms <- colnames(object$factors) %||% paste0("factor", seq_along(sd_x))

  value <- matrix(0, length(maturity), length(sd_x))
  for (j in seq_along(maturity)) {
    n <- maturity[j]
    if (n <= 1L) next
    value[j, ] <- drop(b[n - 1L, ] %*% object$pars$lambda1) * sd_x
  }

  data.frame(
    maturity = rep(maturity, times = length(sd_x)),
    factor = rep(nms, each = length(maturity)),
    value = as.vector(value),
    stringsAsFactors = FALSE
  )
}
