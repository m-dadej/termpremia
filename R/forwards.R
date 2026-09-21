# Forward rates and their decomposition.
#
# Why this exists as its own file. A forward rate is a difference of two log
# prices, and every component of a fitted decomposition is affine in the
# factors, so the forward version of each component is the same difference of
# the same coefficients. That makes forwards nearly free once the recursion
# has run -- and they are what the literature actually argues about. Abrahams,
# Adrian, Crump, Moench and Yu state their headline results for the FIVE-TO-TEN
# YEAR FORWARD breakeven, not for the ten-year spot rate, and a spot
# comparison against a forward claim is not the test it looks like.

#' Forward rates and their decomposition
#'
#' The `start`-to-`end` forward rate implied by a fitted model, split the same
#' way the spot decomposition is.
#'
#' @section How a forward decomposes:
#' A log bond price is affine in the factors, \eqn{p^{(n)}_t = A_n + B_n'X_t},
#' so the rate locked in between `start` and `end` months is
#'
#' \deqn{f_t = -\frac{(A_{m_2} - A_{m_1}) + (B_{m_2} - B_{m_1})'X_t}
#'   {m_2 - m_1}}
#'
#' Every component of the spot decomposition is affine in the same way, so its
#' forward counterpart is that same difference applied to that component's own
#' coefficients. The forward term premium is therefore the difference of the
#' fitted and risk-neutral forwards, exactly as at the spot level, and the
#' identities that hold for spot rates hold here too. The tests assert that.
#'
#' @section Why forwards rather than spot rates:
#' A ten-year spot rate mixes the next five years with the five after that. A
#' five-to-ten year forward removes the near end, which is where policy is
#' known and expectations are sharp, and isolates the part of the curve where
#' risk premia dominate. That is why the source papers state their conclusions
#' in forward terms, and why `analysis/replicate-acmy.R` compares forwards.
#'
#' @param object An `atsm_fit` or `atsm_real_fit`.
#' @param start,end Forward window ends, **in months**, with
#'   `start < end` and both on the fitted grid. `start = 60, end = 120` is the
#'   five-to-ten year forward.
#' @param ... Unused.
#'
#' @return A data frame with `date`, `start`, `end`, `component` and `value`,
#'   in annualised decimals. The components are `fitted`, `risk_neutral` and
#'   `term_premium` for a single-curve fit; a joint fit adds the real curve,
#'   the breakeven, expected inflation, the inflation risk premium and, when
#'   fitted with one, the liquidity components.
#'
#' @seealso [forward_curve()] for the whole forward curve on one date,
#'   [term_premium()], [atsm_real()]
#'
#' @examples
#' panel <- yield_panel(gsw_monthly, units = "percent",
#'                      maturity_unit = "months", issuer = "US")
#' fit <- atsm(panel, n_factors = 5)
#'
#' # The five-to-ten year forward and its decomposition.
#' fwd <- forward_rate(fit, start = 60, end = 120)
#' head(fwd)
#'
#' @export
forward_rate <- function(object, start, end, ...) {
  UseMethod("forward_rate")
}

#' @export
forward_rate.atsm_fit <- function(object, start = 60L, end = 120L, ...) {
  w <- check_forward_window(start, end, object$maturities)

  fitted <- forward_from(object$recursion$p, object$factors, w)
  rn <- forward_from(object$recursion$q, object$factors, w)

  forward_frame(object$dates, w, list(
    fitted = fitted,
    risk_neutral = rn,
    term_premium = fitted - rn
  ))
}

#' @export
forward_rate.atsm_real_fit <- function(object, start = 60L, end = 120L, ...) {
  w_n <- check_forward_window(start, end, object$maturities, "nominal")
  w_r <- check_forward_window(start, end, object$real_maturities, "real")

  co_q <- object$recursion$q
  co_p <- object$recursion$p
  x <- object$factors

  fit_n <- forward_from(list(a = co_q$a, b = co_q$b), x, w_n)
  rn_n <- forward_from(list(a = co_p$a, b = co_p$b), x, w_n)
  fit_r <- forward_from(list(a = co_q$a_real, b = co_q$b_real), x, w_r)
  rn_r <- forward_from(list(a = co_p$a_real, b = co_p$b_real), x, w_r)

  out <- list(
    fitted = fit_n,
    risk_neutral = rn_n,
    term_premium = fit_n - rn_n,
    fitted_real = fit_r,
    risk_neutral_real = rn_r,
    term_premium_real = fit_r - rn_r,
    breakeven = fit_n - fit_r,
    expected_inflation = rn_n - rn_r,
    inflation_risk_premium = (fit_n - rn_n) - (fit_r - rn_r)
  )

  if (object$has_liquidity) {
    j <- object$pca$liq_index
    out$liquidity_nominal <- forward_liquidity(co_q$b, x, w_n, j)
    out$liquidity_real <- forward_liquidity(co_q$b_real, x, w_r, j)
    # Section 2.2 reports inflation measures net of the liquidity factor's
    # own contribution to each curve.
    out$breakeven_liquidity_adjusted <- out$breakeven -
      (out$liquidity_nominal - out$liquidity_real)
  }

  forward_frame(object$dates, w_n, out)
}

#' The forward rate implied by one set of recursion coefficients
#'
#' @param rec A list with `a` and `b` indexed by maturity in months.
#' @param x Factor matrix.
#' @param w Output of `check_forward_window()`.
#' @return A vector of annualised decimal forward rates, one per date.
#' @keywords internal
#' @noRd
forward_from <- function(rec, x, w) {
  d_a <- rec$a[w$end] - rec$a[w$start]
  d_b <- rec$b[w$end, ] - rec$b[w$start, ]
  -(d_a + drop(x %*% d_b)) / w$span * 12
}

#' The liquidity factor's contribution to a forward rate
#'
#' The same difference of coefficients, restricted to the liquidity column.
#'
#' @keywords internal
#' @noRd
forward_liquidity <- function(b, x, w, liq_index) {
  d_b <- b[w$end, liq_index] - b[w$start, liq_index]
  -(d_b * x[, liq_index]) / w$span * 12
}

#' @keywords internal
#' @noRd
check_forward_window <- function(start, end, grid, what = NULL) {
  label <- if (is.null(what)) "" else paste0(" on the ", what, " grid")

  for (nm in c("start", "end")) {
    v <- get(nm)
    if (length(v) != 1L || !is.numeric(v) || !is.finite(v) ||
        v != as.integer(v) || v < 1) {
      stop("`", nm, "` must be a single whole number of months, at least 1.",
           call. = FALSE)
    }
  }
  start <- as.integer(start)
  end <- as.integer(end)

  if (start >= end) {
    stop("`start` (", start, ") must be before `end` (", end,
         "): a forward rate covers the period between them.", call. = FALSE)
  }
  bad <- setdiff(c(start, end), grid)
  if (length(bad)) {
    stop("Maturity(ies) not fitted", label, ": ", paste(bad, collapse = ", "),
         ". Available: ", min(grid), "-", max(grid), " months.", call. = FALSE)
  }

  list(start = start, end = end, span = end - start)
}

#' @keywords internal
#' @noRd
forward_frame <- function(dates, w, components) {
  n <- length(dates)
  data.frame(
    date = rep(dates, times = length(components)),
    start = w$start,
    end = w$end,
    component = rep(names(components), each = n),
    value = unlist(components, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}


# The forward curve on one date -------------------------------------------

#' The one-month forward curve on a single date
#'
#' The cross-section that [forward_rate()] takes a slice of: the rate locked
#' in for each single month ahead, from the observed curve and from the
#' model's fitted and risk-neutral curves.
#'
#' @section Reading it:
#' The gap between the fitted and risk-neutral forward curves is the forward
#' term premium by maturity. Where they run together, the model is saying the
#' curve at that horizon is expectations; where they separate, it is premium.
#' The observed series is there to show the fit, and differs from `fitted` by
#' the pricing error alone.
#'
#' @param object An `atsm_fit` or `atsm_real_fit`.
#' @param date A date in the estimation sample. Defaults to the last one.
#' @param ... Unused.
#'
#' @return A data frame with `maturity`, `component` and `value`. Maturities
#'   run from the second fitted month, since a one-month forward needs a
#'   predecessor.
#'
#' @seealso [forward_rate()]
#'
#' @examples
#' panel <- yield_panel(gsw_monthly, units = "percent",
#'                      maturity_unit = "months", issuer = "US")
#' fit <- atsm(panel, n_factors = 5)
#' fc <- forward_curve(fit)
#' head(fc)
#'
#' @export
forward_curve <- function(object, date = NULL, ...) {
  UseMethod("forward_curve")
}

#' @export
forward_curve.atsm_fit <- function(object, date = NULL, ...) {
  i <- resolve_fit_date(object, date)
  mats <- object$maturities

  forward_curve_frame(mats, list(
    observed = one_month_forwards(object$observed[i, ], mats),
    fitted = one_month_forwards(object$fitted[i, ], mats),
    risk_neutral = one_month_forwards(object$risk_neutral[i, ], mats)
  ))
}

#' @export
forward_curve.atsm_real_fit <- function(object, date = NULL, ...) {
  i <- resolve_fit_date(object, date)
  mats <- object$maturities
  mats_r <- object$real_maturities

  nominal <- forward_curve_frame(mats, list(
    observed = one_month_forwards(object$observed[i, ], mats),
    fitted = one_month_forwards(object$fitted[i, ], mats),
    risk_neutral = one_month_forwards(object$risk_neutral[i, ], mats)
  ))
  real <- forward_curve_frame(mats_r, list(
    observed_real = one_month_forwards(object$observed_real[i, ], mats_r),
    fitted_real = one_month_forwards(object$fitted_real[i, ], mats_r),
    risk_neutral_real = one_month_forwards(object$risk_neutral_real[i, ],
                                           mats_r)
  ))
  out <- rbind(nominal, real)
  rownames(out) <- NULL
  out
}

#' One-month forward rates implied by a spot curve
#'
#' \eqn{f(n-1 \to n) = n\,y_n - (n-1)\,y_{n-1}}, which is the difference of
#' the two log prices. Only defined where the predecessor maturity is on the
#' grid, so a gap in the grid leaves an `NA` rather than a silently wrong
#' number spanning it.
#'
#' @param y A vector of annualised yields for one date.
#' @param maturities Their maturities in months.
#' @return A vector the same length, `NA` where no predecessor exists.
#' @keywords internal
#' @noRd
one_month_forwards <- function(y, maturities) {
  prev <- match(maturities - 1L, maturities)
  out <- rep(NA_real_, length(maturities))
  ok <- !is.na(prev)
  out[ok] <- maturities[ok] * y[ok] - (maturities[ok] - 1L) * y[prev[ok]]
  out
}

#' @keywords internal
#' @noRd
forward_curve_frame <- function(maturities, components) {
  n <- length(maturities)
  out <- data.frame(
    maturity = rep(maturities, times = length(components)),
    component = rep(names(components), each = n),
    value = unlist(components, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  out[!is.na(out$value), , drop = FALSE]
}

#' @keywords internal
#' @noRd
resolve_fit_date <- function(object, date) {
  if (is.null(date)) return(length(object$dates))

  date <- as.Date(date)
  i <- match(date, object$dates)
  if (is.na(i)) {
    stop("No fitted observation on ", format(date), ". The sample runs ",
         format(min(object$dates)), " to ", format(max(object$dates)), ".",
         call. = FALSE)
  }
  i
}
