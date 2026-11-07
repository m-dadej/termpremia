# Fitting Nelson-Siegel and Svensson curves to a handful of observed zero or
# par yields, so that a curve published only at standard tenors can be put on
# the monthly maturity grid ACM needs. svensson_curve() is the evaluation half;
# this file is its inverse.
#
# Units. The fit runs on yields in BASIS POINTS internally, so that sums of
# squared errors are of order one rather than 1e-9 and the optimiser's
# relative tolerances mean something. Parameters are converted back to
# decimals, the package's canonical unit, on the way out.

#' Fit Nelson-Siegel or Svensson curves to sparse zero-coupon or par yields
#'
#' Fits a Nelson-Siegel or Svensson zero-coupon curve to the yields observed on
#' every date of a [yield_panel] -- zero-coupon yields, or the par yields of
#' coupon bonds -- so that a curve published at a handful of tenors can be
#' evaluated on the dense monthly grid that [atsm()] needs. The fitted curves
#' are evaluated with `predict()`, which returns a new `yield_panel` of
#' zero-coupon yields that marks extrapolated maturities.
#'
#' @section Why this is needed:
#' ACM forms one-month holding returns, which need yields at maturities `n` and
#' `n - 1` months, and it extracts its factors from the whole cross-section. A
#' curve quoted at 3m, 6m, 1y, 2y, 5y and 10y supports neither. Every
#' published application bridges the gap with a parametric curve: the New
#' York Fed's own input is the Gurkaynak-Sack-Wright Svensson curve evaluated
#' monthly, which is what [gsw_monthly] is. If you hold published parameters,
#' [svensson_curve()] evaluates them directly and this function is not needed.
#'
#' @section What is estimated:
#' On each date the curve is fitted by least squares on yields, all tenors
#' weighted equally. Given the decay parameters the betas are linear, so they
#' are profiled out and only the decays are searched: a grid over `tau_range`,
#' then bounded local refinements from the most promising grid points. Each
#' date is fitted independently. The default range, 0.1 to 30 years, is the
#' range Gurkaynak, Sack and Wright's own published decays occupy.
#'
#' The decays are not well identified. A fifth of the published GSW Svensson
#' fits have two nearly equal decays with large offsetting betas, and here the
#' first decay doubles or halves from one month to the next in 30% (US) to 48%
#' (UK) of months. That matters less than it sounds: what `atsm()` consumes is
#' the fitted curve, which is well determined even when its parameters are
#' not. Rebuilding the Bank of England curve from eight tenors moved its
#' 10-year ACM term premium by 4bp on average, with a correlation of 0.999 in
#' monthly changes; see `analysis/validate-svensson-fit.R`.
#'
#' The sum of squares has several near-equal basins, and the best can sit in
#' a valley narrower than the grid. The grid therefore only proposes starting
#' points: the best three local minima are each refined and the best result
#' kept.
#'
#' @section Choosing the model:
#' `model = "auto"` fits Svensson when every date has at least seven observed
#' tenors, and Nelson-Siegel when every date has at least five -- one more
#' than each model's number of free parameters, so that there is always a
#' residual to fit. With fewer, only fixed decays are possible. The choice is
#' made once per curve, not per date, because a curve that switched model
#' part-way through its sample would carry an artificial break into the
#' factors.
#'
#' @section Fixed decays, and the factor limit they impose:
#' Supplying `tau` fixes the decays on every date, as Diebold and Li (2006) do
#' with Nelson-Siegel. The fit is then linear and fast, but every fitted curve
#' is a combination of the same three (Nelson-Siegel) or four (Svensson)
#' loadings. The dense panel therefore varies along exactly that many
#' directions, and a five-factor model cannot be estimated on it: its extra
#' principal components are numerically zero. [atsm()] refuses such a panel
#' with an explicit error. Estimated decays do not have this limit.
#'
#' @section Extrapolation:
#' Maturities outside the observed range on a given date are extrapolated.
#' `predict()` flags them, and [atsm()] warns when its default short rate --
#' the curve's one-month yield -- is one of them. From a shortest tenor of
#' three months, the rebuilt one-month yield missed the true one by 1 to 15bp
#' on average across the US and UK curves and both models, and by 20 to 260bp
#' in the worst single month. Every excess return is measured against that
#' yield, so supply an observed bill or policy rate through `short_rate`
#' instead.
#'
#' @section The hold-out check:
#' Comparing a fitted curve with the yields it was fitted to measures nothing
#' useful. With `check = TRUE` each interior tenor is instead dropped in turn,
#' every date is refitted without it, and the refitted curve's prediction at
#' that tenor is compared with the observation. That measures how well the
#' curve interpolates *between* your tenors on your data. It says nothing
#' about extrapolation beyond them.
#'
#' Read it as an upper bound. The refit works from one tenor fewer than the
#' real fit, so it overstates the error: on the Bank of England curve the true
#' error between tenors was a fifth to three quarters of what the check
#' reported.
#'
#' @section Par yields:
#' Yields quoted at standard tenors are very often par yields -- benchmark
#' government bond yields, swap rates, US constant-maturity Treasury yields
#' and most vendor generic series -- which are not zero-coupon yields and
#' cannot be told apart from them by inspection. Treating one as the other
#' biases the fitted curve, and the bias grows with maturity and with the
#' slope of the curve. `yield_type` therefore has no default.
#'
#' With `yield_type = "par"` the curve is fitted to the par yields themselves.
#' On every date the zero-coupon curve is chosen so that bonds paying
#' `coupon_frequency` coupons a year, priced off it, have par yields as close
#' as possible to the observed ones, in least squares. A par bond maturing in
#' `n` years pays coupons at `n`, `n - 1/f`, `n - 2/f` and so on back to the
#' first date after today, that first period being short when `n` is not a
#' whole number of coupon periods, and its par yield is the coupon that prices
#' it at par. A tenor no longer than one coupon period pays a single cash
#' flow, so its par yield is the simple money-market yield: money-market rates
#' and bills at the short end of a coupon curve need no separate treatment.
#'
#' Two things are not modelled. Day counts: maturities and accrual periods
#' are measured in months divided by twelve, so a money-market rate quoted on
#' an actual/360 basis should be scaled by 365/360 first. And the coupon
#' effect: a benchmark bond's yield to maturity is a par yield only when the
#' bond trades at par. Away from par the error grows with the gap between its
#' coupon and its yield; for a bond near par it is a small fraction of the
#' error made by treating the yield as zero-coupon.
#'
#' The decay search is the zero-yield one, run on par yields converted to
#' zero yields through the current curve, and its result is then refined
#' against the par yields directly. A par fit therefore takes several times
#' as long as a zero fit -- about five times for Svensson -- and the hold-out
#' check multiplies that again. The `rmse` in [coef()] and the hold-out check
#' measure errors in par yields, the quantity fitted.
#'
#' Rebuilt this way from eight semi-annual par yields, the US and UK curves
#' come back within 0.0 and 0.3bp of the originals between tenors (Svensson),
#' and the 10-year ACM term premium within 0.0 and 4.9bp on average, about
#' what fitting their true zero yields achieves. Taking the same par yields
#' for zero yields misses the curves by 14 and 7bp and moves the US premium
#' by up to 92bp. See `analysis/validate-svensson-fit.R`.
#'
#' @param panel A [yield_panel]. Every curve in it is fitted.
#' @param yield_type What the yields are: `"zero"` for zero-coupon yields, or
#'   `"par"` for the par yields of coupon bonds. Must be supplied; see the
#'   section on par yields. One value for every curve, or a named vector keyed
#'   by curve name.
#' @param coupon_frequency Coupons a year paid by the bonds whose par yields
#'   are supplied: 1, 2, 4 or 12. Required for par yields, which are quoted at
#'   this compounding frequency, so it cannot be inferred from them: 1 for
#'   annual coupons (Bunds, most central European government bonds), 2 for
#'   semi-annual (US Treasuries, gilts). Ignored for zero-coupon yields. One
#'   value, or a named vector keyed by curve name.
#' @param model `"auto"`, `"svensson"` or `"nelson_siegel"`.
#' @param tau Optional fixed decays, in years: one value for Nelson-Siegel,
#'   two for Svensson. `NULL`, the default, estimates them on every date.
#' @param tau_range Bounds for estimated decays, in years.
#' @param check Run the hold-out check described above. It refits every date
#'   once per interior tenor, so it multiplies the running time.
#'
#' @return An object of class `svensson_fit`. Use `predict()` to evaluate it
#'   on a maturity grid, `coef()` for the parameters, and `print()` for fit
#'   and hold-out diagnostics.
#'
#' @references
#' Diebold, F. X. and C. Li (2006). "Forecasting the term structure of
#' government bond yields." *Journal of Econometrics* 130(2), 337-364.
#'
#' Nelson, C. R. and A. F. Siegel (1987). "Parsimonious modeling of yield
#' curves." *Journal of Business* 60(4), 473-489.
#'
#' Svensson, L. E. O. (1994). "Estimating and interpreting forward interest
#' rates: Sweden 1992-1994." NBER Working Paper 4871.
#'
#' @seealso [predict.svensson_fit()], [svensson_curve()], [yield_panel()]
#'
#' @examples
#' # The US curve observed only at eight standard tenors, over ten years
#' tenors <- c(3, 6, 12, 24, 36, 60, 84, 120)
#' sparse <- gsw_monthly[gsw_monthly$maturity %in% tenors &
#'                       gsw_monthly$date >= as.Date("2015-01-01"), ]
#' panel <- yield_panel(sparse, units = "percent", maturity_unit = "months",
#'                      issuer = "US")
#'
#' fit <- svensson_fit(panel, yield_type = "zero")
#' fit
#'
#' # Back on a monthly grid, extrapolation below three months flagged
#' dense <- predict(fit, maturities = 1:120)
#' dense
#'
#' \donttest{
#' # Benchmark bond yields are par yields. Here, money-market rates at 3m and
#' # 6m and annual-coupon par yields from 1y, as on most European curves.
#' params <- data.frame(beta0 = 0.05, beta1 = -0.02, beta2 = 0.01,
#'                      beta3 = 0, tau1 = 2, tau2 = NA)
#' tenors <- c(3, 6, 12, 24, 36, 60, 84, 120)
#' zero_at <- function(m) drop(svensson_curve(params, m / 12))
#' par <- vapply(tenors, function(n) {
#'   pay <- rev(seq(n, 1e-9, by = -12))           # coupon dates, in months
#'   d <- exp(-zero_at(pay) * pay / 12)
#'   (1 - d[length(d)]) / sum(diff(c(0, pay)) / 12 * d)
#' }, numeric(1))
#' bench <- yield_panel(matrix(par, 1, dimnames = list("2024-01-31", tenors)),
#'                      units = "decimal", maturity_unit = "months")
#'
#' fit_par <- svensson_fit(bench, yield_type = "par", coupon_frequency = 1,
#'                         check = FALSE)
#' # Zero-coupon yields again: the 10-year point matches the true curve
#' curve_matrix(predict(fit_par, maturities = 120)) - zero_at(120)
#' }
#'
#' @export
svensson_fit <- function(panel,
                         yield_type,
                         coupon_frequency = NULL,
                         model = c("auto", "svensson", "nelson_siegel"),
                         tau = NULL,
                         tau_range = c(0.1, 30),
                         check = TRUE) {
  if (!inherits(panel, "yield_panel")) {
    stop("`panel` must be a yield_panel; see ?yield_panel.", call. = FALSE)
  }
  if (missing(yield_type)) {
    stop("Declare what the yields are: `yield_type = \"zero\"` for ",
         "zero-coupon yields, or `yield_type = \"par\"` for par yields. ",
         "Yields published at a handful of standard tenors are often PAR ",
         "yields (benchmark bond yields, swap rates, US constant-maturity ",
         "Treasury yields, most vendor generic series), and par and ",
         "zero-coupon yields cannot be told apart from the numbers. See ",
         "?svensson_fit.", call. = FALSE)
  }
  curves <- stats::setNames(nm = panel$meta$curve)
  yield_type <- check_yield_type(yield_type, curves)
  freq <- check_coupon_frequency(coupon_frequency, yield_type, curves)
  model <- match.arg(model)
  tau_range <- check_tau_range(tau_range)
  tau <- check_fixed_tau(tau, model)

  fits <- lapply(curves, function(nm) {
    fit_one_curve(curve_matrix(panel, nm), panel$maturities[[nm]], nm,
                  model, tau, tau_range, check, freq[[nm]])
  })

  for (nm in curves) {
    dropped <- fits[[nm]]$dropped
    if (length(dropped)) {
      warning(
        length(dropped), " date(s) on curve '", nm, "' have fewer than ",
        fits[[nm]]$min_obs, " observed maturities and were left unfitted, ",
        "the first on ", format(panel$dates[dropped[1L]]), ". Their rows ",
        "are missing in predict(), and atsm() refuses a panel with missing ",
        "rows: restrict the panel to dates with enough tenors, or fix the ",
        "decays through `tau`, which needs fewer.", call. = FALSE
      )
    }
  }

  structure(
    list(
      curves = fits,
      dates = panel$dates,
      meta = panel$meta,
      frequency = panel$frequency,
      tau_range = tau_range
    ),
    class = "svensson_fit"
  )
}


# Argument checks ---------------------------------------------------------

#' One value per curve, from a scalar, a vector in curve order, or a vector
#' named by curve
#'
#' @param need Curves that must have an entry when `v` is named. The others
#'   come back `NA`.
#' @keywords internal
#' @noRd
per_curve <- function(v, curves, arg, need = curves) {
  if (!is.null(names(v))) {
    missing_curves <- setdiff(need, names(v))
    if (length(missing_curves)) {
      stop("`", arg, "` is named but has no entry for curve(s): ",
           paste(missing_curves, collapse = ", "), call. = FALSE)
    }
    return(stats::setNames(unname(v[match(curves, names(v))]), curves))
  }
  if (length(v) == 1L) return(stats::setNames(rep(v, length(curves)), curves))
  if (length(v) == length(curves)) return(stats::setNames(v, curves))
  stop("`", arg, "` must be length 1, length ", length(curves), ", or a named ",
       "vector keyed by curve name.", call. = FALSE)
}

#' @keywords internal
#' @noRd
check_yield_type <- function(yield_type, curves) {
  msg <- "`yield_type` must be \"zero\" or \"par\"."
  if (!is.character(yield_type) || !length(yield_type)) {
    stop(msg, call. = FALSE)
  }
  out <- per_curve(yield_type, curves, "yield_type")
  if (anyNA(out) || !all(out %in% c("zero", "par"))) stop(msg, call. = FALSE)
  out
}

#' Coupon frequency of each par curve
#'
#' @return A list named by curve: the frequency for par curves, `NULL` for
#'   zero curves, which is what the fitting code takes to mean zero yields.
#' @keywords internal
#' @noRd
check_coupon_frequency <- function(coupon_frequency, yield_type, curves) {
  out <- stats::setNames(vector("list", length(curves)), curves)
  par_curves <- curves[yield_type == "par"]

  if (!length(par_curves)) {
    if (!is.null(coupon_frequency)) {
      warning("`coupon_frequency` was supplied but no curve has ",
              "`yield_type = \"par\"`, so it was ignored.", call. = FALSE)
    }
    return(out)
  }
  if (is.null(coupon_frequency)) {
    stop("Par yields need `coupon_frequency`, the number of coupons a year ",
         "paid by the bonds they describe: 1 for annual coupons (Bunds, most ",
         "central European government bonds), 2 for semi-annual (US ",
         "Treasuries, gilts). It sets the compounding of a par yield, so it ",
         "cannot be inferred from the numbers.", call. = FALSE)
  }

  msg <- "`coupon_frequency` must be 1, 2, 4 or 12 coupons a year."
  if (!is.numeric(coupon_frequency)) stop(msg, call. = FALSE)
  f <- per_curve(coupon_frequency, curves, "coupon_frequency",
                 need = par_curves)[par_curves]
  if (anyNA(f) || !all(f %in% c(1, 2, 4, 12))) stop(msg, call. = FALSE)

  out[par_curves] <- as.list(as.numeric(f))
  out
}

#' @keywords internal
#' @noRd
check_tau_range <- function(tau_range) {
  if (!is.numeric(tau_range) || length(tau_range) != 2L ||
      anyNA(tau_range) || any(!is.finite(tau_range)) ||
      tau_range[1L] <= 0 || tau_range[1L] >= tau_range[2L]) {
    stop("`tau_range` must be two increasing, positive numbers of years.",
         call. = FALSE)
  }
  as.numeric(tau_range)
}

#' @keywords internal
#' @noRd
check_fixed_tau <- function(tau, model) {
  if (is.null(tau)) return(NULL)

  if (!is.numeric(tau) || !length(tau) %in% 1:2 || anyNA(tau) ||
      any(!is.finite(tau)) || any(tau <= 0)) {
    stop("`tau` must be one positive number of years (Nelson-Siegel) or two ",
         "(Svensson).", call. = FALSE)
  }
  implied <- if (length(tau) == 1L) "nelson_siegel" else "svensson"
  if (model != "auto" && model != implied) {
    stop("`model = \"", model, "\"` needs ", if (model == "svensson") "two" else
         "one", " decay(s) in `tau`, but ", length(tau), " were supplied.",
         call. = FALSE)
  }
  # Equal decays make the two hump loadings identical, so the fourth beta is
  # not identified by any amount of data.
  if (length(tau) == 2L && isTRUE(all.equal(tau[1L], tau[2L]))) {
    stop("The two Svensson decays in `tau` must differ; equal decays make the ",
         "two hump terms identical.", call. = FALSE)
  }
  as.numeric(tau)
}


# Fitting -----------------------------------------------------------------

#' Fit one curve of a panel on every date
#'
#' @param y `T x N` decimal yields, `NA` where unobserved.
#' @param maturities Maturities in months, matching the columns of `y`.
#' @param freq Coupons a year when `y` holds par yields; `NULL` for zero
#'   yields.
#' @return A list describing the fit; see the fields assembled at the end.
#' @keywords internal
#' @noRd
fit_one_curve <- function(y, maturities, name, model, tau, tau_range, check,
                          freq = NULL) {
  n_obs <- rowSums(!is.na(y))
  model <- resolve_curve_model(model, tau, n_obs, name)
  min_obs <- curve_min_obs(model, tau)

  params <- curve_params(y * 1e4, maturities / 12, model, tau, tau_range,
                         min_obs, freq)
  # Dates with no observations at all are not failures -- a TIPS curve that
  # starts decades after its nominal partner has hundreds -- only dates that
  # had yields but too few of them.
  dropped <- which(n_obs > 0L & is.na(params[, "beta0"]))

  observed_range <- function(f) {
    apply(y, 1L, function(r) {
      if (all(is.na(r))) NA_real_ else f(maturities[!is.na(r)])
    })
  }

  list(
    model = model,
    yield_type = if (is.null(freq)) "zero" else "par",
    coupon_frequency = freq,
    tau = tau,
    params = params,
    maturities = maturities,
    n_obs = n_obs,
    shortest = observed_range(min),
    longest = observed_range(max),
    min_obs = min_obs,
    dropped = dropped,
    check = if (isTRUE(check)) {
      holdout_check(y, maturities, model, tau, tau_range, min_obs, freq)
    }
  )
}

#' @keywords internal
#' @noRd
resolve_curve_model <- function(model, tau, n_obs, name) {
  if (!is.null(tau)) {
    return(if (length(tau) == 1L) "nelson_siegel" else "svensson")
  }
  if (model != "auto") return(model)

  have <- n_obs[n_obs > 0]
  if (!length(have)) {
    stop("Curve '", name, "' has no observed yields.", call. = FALSE)
  }
  fewest <- min(have)
  if (fewest >= 7L) return("svensson")
  if (fewest >= 5L) return("nelson_siegel")

  stop("Curve '", name, "' has as few as ", fewest, " observed maturities on ",
       "some dates. With estimated decays, Nelson-Siegel needs at least 5 and ",
       "Svensson at least 7 -- one more than their free parameters. Fix the ",
       "decays through `tau`, which needs 4 (Nelson-Siegel) or 5 (Svensson), ",
       "or drop the sparse dates.", call. = FALSE)
}

#' Observations needed for at least one residual degree of freedom
#' @keywords internal
#' @noRd
curve_min_obs <- function(model, tau) {
  n_beta <- if (model == "svensson") 4L else 3L
  n_tau <- if (!is.null(tau)) 0L else if (model == "svensson") 2L else 1L
  n_beta + n_tau + 1L
}

#' Loadings of the Nelson-Siegel or Svensson curve
#'
#' The same loadings as [svensson_yield()], written out directly because the
#' decay search builds this matrix tens of thousands of times per curve. The
#' zero-maturity limit that `ns_level_factor()` guards against cannot arise:
#' observed maturities are strictly positive and decays at most `tau_range`.
#'
#' @param m Maturities in years, all positive.
#' @param tau1,tau2 Decays in years; `tau2 = NULL` gives Nelson-Siegel.
#' @keywords internal
#' @noRd
curve_design <- function(m, tau1, tau2 = NULL) {
  x1 <- m / tau1
  e1 <- exp(-x1)
  level <- (1 - e1) / x1
  if (is.null(tau2)) return(cbind(1, level, level - e1))

  x2 <- m / tau2
  e2 <- exp(-x2)
  cbind(1, level, level - e1, (1 - e2) / x2 - e2)
}

#' Sum of squared residuals of each row of `y` projected on `x`
#'
#' Residuals are formed explicitly rather than as `y'y - y'QQ'y`: on curves
#' that fit almost exactly the difference cancels catastrophically, and the
#' grid search would then pick its minimum from rounding noise. `.lm.fit()`
#' pivots out aliased columns, which two nearly equal Svensson decays produce,
#' and it is a single compiled call -- this function is evaluated tens of
#' thousands of times per curve, so that matters.
#'
#' @param y A `G x N` matrix, one date per row, or a vector for one date.
#' @keywords internal
#' @noRd
pattern_ssr <- function(y, x) {
  y <- if (is.matrix(y)) t(y) else matrix(y, ncol = 1L)
  colSums(stats::.lm.fit(x, y)$residuals^2)
}

#' Fit every date of the panel
#'
#' Dates are grouped by which maturities they observe, so that within a group
#' the loadings are shared and the grid search runs on all its dates at once.
#'
#' @param y `T x N` yields in basis points.
#' @param m Maturities in years.
#' @param freq Coupons a year when `y` holds par yields; `NULL` for zero
#'   yields.
#' @return A `T x 7` matrix: `beta0`-`beta3` and `rmse` in decimals, `tau1`
#'   and `tau2` in years. Unfittable dates are `NA`.
#' @keywords internal
#' @noRd
curve_params <- function(y, m, model, tau, tau_range, min_obs, freq = NULL) {
  out <- matrix(NA_real_, nrow(y), 7L, dimnames = list(
    rownames(y), c("beta0", "beta1", "beta2", "beta3", "tau1", "tau2", "rmse")
  ))

  obs <- !is.na(y)
  key <- apply(obs, 1L, function(r) paste(which(r), collapse = ","))

  for (k in unique(key)) {
    rows <- which(key == k)
    cols <- which(obs[rows[1L], ])
    if (length(cols) < min_obs) next
    out[rows, ] <- if (is.null(freq)) {
      fit_pattern(y[rows, cols, drop = FALSE], m[cols], model, tau, tau_range)
    } else {
      fit_pattern_par(y[rows, cols, drop = FALSE], m[cols], model, tau,
                      tau_range, freq)
    }
  }
  out
}

#' Fit a group of dates that share observed maturities
#' @keywords internal
#' @noRd
fit_pattern <- function(y, m, model, tau, tau_range) {
  sv <- model == "svensson"

  taus <- if (!is.null(tau)) {
    matrix(tau, nrow(y), length(tau), byrow = TRUE)
  } else if (sv) {
    search_svensson(y, m, tau_range)
  } else {
    matrix(search_nelson_siegel(y, m, tau_range), ncol = 1L)
  }

  out <- matrix(NA_real_, nrow(y), 7L)
  for (i in seq_len(nrow(y))) {
    x <- curve_design(m, taus[i, 1L], if (sv) taus[i, 2L])
    q <- qr(x)
    beta <- qr.coef(q, y[i, ])
    # An aliased column only arises with two nearly equal decays, where the
    # remaining columns already reproduce the curve; drop its contribution.
    beta[is.na(beta)] <- 0
    resid <- y[i, ] - drop(x %*% beta)
    if (!sv) beta <- c(beta, 0)
    out[i, ] <- c(beta / 1e4, taus[i, 1L],
                  if (sv) taus[i, 2L] else NA_real_,
                  sqrt(mean(resid^2)) / 1e4)
  }
  out
}

# The decay search. The sum of squares is very flat in the decays, with
# several basins that fit the tenors almost equally well, and the true
# minimum can sit in a valley narrower than the grid spacing: on exact
# Nelson-Siegel data the grid's best cell was a basin with a 0.1bp fit while
# the exact one lay between two grid points that each scored worse. So the
# grid only proposes starting points; the best few LOCAL minima are each
# refined, and the best refinement wins.

#' Grid of decays, log-spaced over `tau_range`
#' @keywords internal
#' @noRd
tau_grid <- function(tau_range, n) {
  exp(seq(log(tau_range[1L]), log(tau_range[2L]), length.out = n))
}

#' Log-decay bounds of the grid cells either side of point `b`
#' @keywords internal
#' @noRd
tau_cell <- function(grid, b) {
  log(grid[c(max(b - 1L, 1L), min(b + 1L, length(grid)))])
}

#' Decay search for Nelson-Siegel: grid, then Brent from each local minimum
#' @keywords internal
#' @noRd
search_nelson_siegel <- function(y, m, tau_range, n_grid = 60L,
                                 n_start = 3L) {
  grid <- tau_grid(tau_range, n_grid)
  ssr <- vapply(grid, function(t) pattern_ssr(y, curve_design(m, t)),
                numeric(nrow(y)))
  if (!is.matrix(ssr)) ssr <- matrix(ssr, nrow = 1L)

  pad <- cbind(Inf, ssr, Inf)
  is_min <- ssr <= pad[, seq_len(n_grid), drop = FALSE] &
    ssr <= pad[, seq_len(n_grid) + 2L, drop = FALSE]

  vapply(seq_len(nrow(y)), function(i) {
    starts <- which(is_min[i, ])
    starts <- starts[order(ssr[i, starts])][seq_len(min(n_start,
                                                        length(starts)))]
    f <- function(lt) pattern_ssr(y[i, ], curve_design(m, exp(lt)))

    best_tau <- grid[starts[1L]]
    best_value <- ssr[i, starts[1L]]
    for (b in starts) {
      o <- stats::optimize(f, tau_cell(grid, b))
      if (o$objective < best_value) {
        best_tau <- exp(o$minimum)
        best_value <- o$objective
      }
    }
    best_tau
  }, numeric(1))
}

#' Decay search for Svensson: 2-D grid, then a bounded refinement from each
#' of the best local minima
#'
#' Both orderings of the decays are searched. They are not interchangeable,
#' because the first decay also shapes the slope loading. The diagonal, where
#' the two hump loadings coincide, is excluded.
#'
#' Dates are processed in blocks so that the `dates x grid x grid` array of
#' sums of squares stays small on long daily panels.
#'
#' @keywords internal
#' @noRd
search_svensson <- function(y, m, tau_range, n_grid = 30L, n_start = 3L,
                            block = 500L) {
  if (nrow(y) > block) {
    blocks <- split(seq_len(nrow(y)), ceiling(seq_len(nrow(y)) / block))
    return(do.call(rbind, lapply(blocks, function(r) {
      search_svensson(y[r, , drop = FALSE], m, tau_range, n_grid, n_start,
                      block)
    })))
  }

  grid <- tau_grid(tau_range, n_grid)
  n <- nrow(y)
  g <- seq_len(n_grid)

  ssr <- array(Inf, c(n, n_grid, n_grid))
  for (i in g) {
    for (j in g) {
      if (i == j) next
      ssr[, i, j] <- pattern_ssr(y, curve_design(m, grid[i], grid[j]))
    }
  }

  # A cell is a local minimum when none of its eight neighbours is lower.
  pad <- array(Inf, c(n, n_grid + 2L, n_grid + 2L))
  pad[, g + 1L, g + 1L] <- ssr
  is_min <- is.finite(ssr)
  for (di in -1:1) {
    for (dj in -1:1) {
      if (di == 0L && dj == 0L) next
      is_min <- is_min & ssr <= pad[, g + 1L + di, g + 1L + dj, drop = FALSE]
    }
  }

  out <- matrix(NA_real_, n, 2L)
  for (k in seq_len(n)) {
    s <- ssr[k, , ]
    starts <- which(is_min[k, , ], arr.ind = TRUE)
    starts <- starts[order(s[starts])[seq_len(min(n_start, nrow(starts)))], ,
                     drop = FALSE]
    f <- function(p) pattern_ssr(y[k, ], curve_design(m, exp(p[1L]),
                                                     exp(p[2L])))

    best_value <- s[starts[1L, , drop = FALSE]]
    best_tau <- grid[starts[1L, ]]
    for (r in seq_len(nrow(starts))) {
      li <- tau_cell(grid, starts[r, 1L])
      lj <- tau_cell(grid, starts[r, 2L])
      # The objective is a sum of squares in bp^2, so a relative tolerance of
      # about 1e-6 is far below anything visible in a fitted yield.
      o <- stats::optim(log(grid[starts[r, ]]), f, method = "L-BFGS-B",
                        lower = c(li[1L], lj[1L]), upper = c(li[2L], lj[2L]),
                        control = list(factr = 1e10))
      if (o$value < best_value) {
        best_value <- o$value
        best_tau <- exp(o$par)
      }
    }
    out[k, ] <- best_tau
  }
  out
}

#' Evaluate fitted parameters at maturities given in years
#' @keywords internal
#' @noRd
curve_evaluate <- function(params, m) {
  svensson_curve(as.data.frame(params[, c("beta0", "beta1", "beta2", "beta3",
                                          "tau1", "tau2"), drop = FALSE]), m)
}


# Par yields --------------------------------------------------------------
#
# A par yield is not a point on the zero curve but a summary of it: the coupon
# that prices a bond at par, which depends on every discount factor up to the
# bond's maturity. The curve is therefore fitted by pricing those bonds off
# it. Given the decays the betas are no longer linear, because discount
# factors are exponential in them, but they are very nearly so -- a par yield
# sits a smooth few tens of basis points at most from the zero yield of the
# same maturity -- and Gauss-Newton converges in a few steps from the fit that
# treats par yields as zero yields.
#
# The decay search stays the zero-yield one. Its grid is what protects the fit
# from the narrow basins described above, and the grid needs linear betas to
# be cheap. So it is run on PSEUDO-ZERO yields: the observed par yields plus
# the current curve's zero-minus-par spread at each tenor, which are exact
# zero yields once the curve is right. The search is repeated with the spread
# of the latest curve until the spread stops moving, which on smooth curves
# takes two rounds, and the decays it settles on are then refined against the
# par yields themselves.

#' Cash flows of a par bond at each tenor
#'
#' Coupon dates run back from maturity in steps of one coupon period to the
#' first date after today; the first period is short when the maturity is not
#' a whole number of periods. Computed in months, where a coupon period is a
#' whole number, so that coupon dates shared by several tenors coincide
#' exactly.
#'
#' @param m Maturities in years.
#' @param freq Coupons a year.
#' @return A list: `times`, the distinct cash-flow times in years; `w`, an
#'   `N x K` matrix holding each tenor's accrual fraction at each of its
#'   cash-flow times and zero elsewhere; and `final`, the column of each
#'   tenor's maturity.
#' @keywords internal
#' @noRd
par_schedule <- function(m, freq) {
  period <- 12 / freq
  flows <- lapply(m * 12, function(n) {
    # The tolerance keeps a maturity that is a whole number of periods, up to
    # rounding, from growing a spurious stub of a few nanoseconds.
    k <- ceiling(n / period - 1e-6) - 1
    rev(n - (0:k) * period)
  })
  key <- function(x) round(x, 6)
  times <- sort(unique(key(unlist(flows))))

  w <- matrix(0, length(m), length(times))
  final <- integer(length(m))
  for (i in seq_along(flows)) {
    j <- match(key(flows[[i]]), times)
    w[i, j] <- diff(c(0, flows[[i]])) / 12
    final[i] <- j[length(j)]
  }
  list(times = times / 12, w = w, final = final)
}

#' Par yields from zero-coupon yields at the cash-flow times
#'
#' @param z Continuously compounded decimal zero yields at `sched$times`: a
#'   vector, or a matrix with one date per row.
#' @return A matrix of decimal par yields, one row per date and one column per
#'   tenor of `sched`.
#' @keywords internal
#' @noRd
par_price <- function(z, sched) {
  if (!is.matrix(z)) z <- matrix(z, nrow = 1L)
  d <- exp(-sweep(z, 2L, sched$times, "*"))
  unname((1 - d[, sched$final, drop = FALSE]) / (d %*% t(sched$w)))
}

#' Par yields of fitted parameters, one row per date
#' @keywords internal
#' @noRd
par_evaluate <- function(params, sched) {
  par_price(curve_evaluate(params, sched$times), sched)
}

#' Betas that fit par yields given the decays, by Gauss-Newton
#'
#' The default start is the linear fit that treats the par yields as zero
#' yields, which is already within a few tens of basis points; `start`, the
#' betas of nearby decays, is used instead when it fits better. Steps that
#' would raise the sum of squares are halved, and an aliased column -- only
#' possible with two nearly equal Svensson decays -- gets no step, as in
#' `fit_pattern()`.
#'
#' @param obs Decimal par yields at tenors `m` (years).
#' @param taus One decay (Nelson-Siegel) or two (Svensson), in years.
#' @param start Optional starting betas.
#' @return A list with `beta`, `resid` (decimal) and `ssr` (basis points
#'   squared).
#' @keywords internal
#' @noRd
par_betas <- function(obs, m, sched, taus, start = NULL, max_iter = 20L) {
  tau2 <- if (length(taus) == 2L) taus[2L]
  xc <- curve_design(sched$times, taus[1L], tau2)
  tt <- sched$times

  price <- function(beta) {
    d <- exp(-drop(xc %*% beta) * tt)
    a <- drop(sched$w %*% d)
    p <- (1 - d[sched$final]) / a
    list(beta = beta, d = d, a = a, p = p, ssr = sum((obs - p)^2))
  }

  beta <- qr.coef(qr(curve_design(m, taus[1L], tau2)), obs)
  beta[is.na(beta)] <- 0
  cur <- price(beta)
  if (!is.null(start)) {
    warm <- price(start)
    if (is.finite(warm$ssr) && (!is.finite(cur$ssr) || warm$ssr < cur$ssr)) {
      cur <- warm
    }
  }

  for (it in seq_len(max_iter)) {
    if (!is.finite(cur$ssr)) break
    # d p / d beta, from p = (1 - d_n) / a and d d_t / d beta = -t d_t x_t.
    dd <- -(tt * cur$d) * xc
    jac <- (-dd[sched$final, , drop = FALSE] -
              cur$p * (sched$w %*% dd)) / cur$a
    step <- qr.coef(qr(jac), obs - cur$p)
    step[is.na(step)] <- 0

    h <- 1
    new <- price(cur$beta + step)
    while ((!is.finite(new$ssr) || new$ssr > cur$ssr) && h > 1e-4) {
      h <- h / 2
      new <- price(cur$beta + h * step)
    }
    if (!is.finite(new$ssr) || new$ssr > cur$ssr) break
    cur <- new
    # Converged once the step no longer moves any fitted par yield by more
    # than 1e-6bp. Measured on yields rather than betas, because two nearly
    # equal Svensson decays leave the betas free to drift along a direction
    # the curve does not see.
    if (max(abs(h * (jac %*% step))) < 1e-10) break
  }

  resid <- obs - cur$p
  list(beta = cur$beta, resid = resid, ssr = sum((resid * 1e4)^2))
}

#' Fit one date's par yields at given decays, optionally refining the decays
#'
#' With `refine = TRUE` the decays move within the basin the search chose, by
#' at most `width` in log terms, to where the par yields themselves fit best.
#' The proposal is kept if nothing does better.
#'
#' @return A list with `row`, the date's seven parameters as in
#'   `curve_params()`, and `ssr` in basis points squared.
#' @keywords internal
#' @noRd
par_refine <- function(obs, m, sched, taus, refine, tau_range, width = 0.2) {
  best <- par_betas(obs, m, sched, taus)
  best_taus <- taus

  if (refine) {
    lt <- log(taus)
    lower <- pmax(lt - width, log(tau_range[1L]))
    upper <- pmin(lt + width, log(tau_range[2L]))
    # Each evaluation starts Gauss-Newton from the betas of the last, whose
    # decays the optimiser has only just left, which roughly halves the
    # iterations it takes.
    last <- best$beta
    f <- function(l) {
      s <- par_betas(obs, m, sched, exp(l), start = last)
      if (is.finite(s$ssr)) last <<- s$beta
      if (is.finite(s$ssr)) s$ssr else 1e30
    }
    cand <- if (length(lt) == 2L) {
      stats::optim(lt, f, method = "L-BFGS-B", lower = lower, upper = upper,
                   control = list(factr = 1e10))$par
    } else {
      stats::optimize(f, c(lower, upper))$minimum
    }
    refined <- par_betas(obs, m, sched, exp(cand), start = last)
    if (refined$ssr < best$ssr) {
      best <- refined
      best_taus <- exp(cand)
    }
  }

  sv <- length(taus) == 2L
  list(
    row = c(if (sv) best$beta else c(best$beta, 0), best_taus[1L],
            if (sv) best_taus[2L] else NA_real_, sqrt(mean(best$resid^2))),
    ssr = best$ssr
  )
}

#' Fit a group of dates that share observed maturities, from par yields
#'
#' @param y `G x N` par yields in basis points.
#' @param m Maturities in years.
#' @param tol Change in the pseudo-zero yields, in basis points, below which
#'   another round of the decay search would not move the result.
#' @keywords internal
#' @noRd
fit_pattern_par <- function(y, m, model, tau, tau_range, freq,
                            max_rounds = 3L, tol = 0.01) {
  sv <- model == "svensson"
  sched <- par_schedule(m, freq)
  obs <- y / 1e4
  n <- nrow(y)

  out <- matrix(NA_real_, n, 7L, dimnames = list(NULL, c(
    "beta0", "beta1", "beta2", "beta3", "tau1", "tau2", "rmse"
  )))
  best_ssr <- rep(Inf, n)
  keep_better <- function(i, r) {
    if (r$ssr < best_ssr[i]) {
      out[i, ] <<- r$row
      best_ssr[i] <<- r$ssr
    }
  }

  # Fixed decays leave nothing to search: the betas are directly the least
  # squares fit to the par yields.
  if (!is.null(tau)) {
    for (i in seq_len(n)) {
      keep_better(i, par_refine(obs[i, ], m, sched, tau, FALSE, tau_range))
    }
    return(out)
  }

  # The rounds only have to choose the basin, so the betas are solved at the
  # proposed decays and the decays are left alone; refining them against the
  # par yields happens once, from the best proposal, at the end.
  pseudo <- y   # the first round treats the par yields as zero yields
  for (round in seq_len(max_rounds)) {
    taus <- if (sv) {
      search_svensson(pseudo, m, tau_range)
    } else {
      matrix(search_nelson_siegel(pseudo, m, tau_range), ncol = 1L)
    }
    for (i in seq_len(n)) {
      keep_better(i, par_refine(obs[i, ], m, sched, taus[i, ], FALSE,
                                tau_range))
    }

    spread <- (curve_evaluate(out, m) - par_evaluate(out, sched)) * 1e4
    moved <- max(abs(y + spread - pseudo))
    pseudo <- y + spread
    if (moved < tol) break
  }

  tau_cols <- if (sv) c("tau1", "tau2") else "tau1"
  for (i in seq_len(n)) {
    keep_better(i, par_refine(obs[i, ], m, sched, out[i, tau_cols], TRUE,
                              tau_range))
  }
  out
}

#' Drop each interior tenor, refit, and predict it
#'
#' For par yields the prediction is the held-out tenor's par yield, priced off
#' the refitted curve, so the error is in the quantity that was observed.
#'
#' @return A data frame with one row per interior tenor: `maturity` (months),
#'   `n_dates`, and `rmse` and `max_abs` in decimals. `NULL` when there is no
#'   interior tenor or no tenor can be dropped without leaving too few.
#' @keywords internal
#' @noRd
holdout_check <- function(y, maturities, model, tau, tau_range, min_obs,
                          freq = NULL) {
  ord <- sort(maturities)
  inner <- ord[-c(1L, length(ord))]
  if (!length(inner)) return(NULL)

  rows <- lapply(inner, function(h) {
    j <- match(h, maturities)
    p <- curve_params(y[, -j, drop = FALSE] * 1e4, maturities[-j] / 12,
                      model, tau, tau_range, min_obs, freq)
    predicted <- if (is.null(freq)) {
      curve_evaluate(p, h / 12)
    } else {
      par_evaluate(p, par_schedule(h / 12, freq))
    }
    err <- drop(predicted) - y[, j]
    ok <- !is.na(err)
    data.frame(
      maturity = h,
      n_dates = sum(ok),
      rmse = if (any(ok)) sqrt(mean(err[ok]^2)) else NA_real_,
      max_abs = if (any(ok)) max(abs(err[ok])) else NA_real_
    )
  })
  out <- do.call(rbind, rows)
  out[out$n_dates > 0L, , drop = FALSE]
}


# Methods -----------------------------------------------------------------

#' Evaluate fitted curves on a maturity grid
#'
#' @param object A `svensson_fit` from [svensson_fit()].
#' @param maturities Maturities in months. `NULL` gives every whole month from
#'   1 to each curve's longest observed maturity. A numeric vector applies to
#'   every curve; a named list sets each curve's grid separately.
#' @param ... Unused.
#'
#' @return A [yield_panel] of continuously compounded zero-coupon yields --
#'   also when the curve was fitted to par yields -- with the same dates,
#'   curves, instruments and issuers as the panel that was fitted. Maturities
#'   outside the observed
#'   range on each date are flagged as extrapolated, which [atsm()] reads.
#'   Dates that could not be fitted are kept as missing rows rather than
#'   dropped, so that a gap cannot pass unnoticed into the factor VAR.
#'
#' @seealso [svensson_fit()]
#' @export
predict.svensson_fit <- function(object, maturities = NULL, ...) {
  meta <- object$meta
  tt <- length(object$dates)

  long <- lapply(meta$curve, function(nm) {
    cf <- object$curves[[nm]]
    mats <- predict_maturities(maturities, nm, cf$maturities)

    y <- curve_evaluate(cf$params, mats / 12)
    flag <- outer(cf$shortest, mats, ">") | outer(cf$longest, mats, "<")
    flag[is.na(flag)] <- FALSE

    data.frame(
      date = rep(object$dates, times = length(mats)),
      maturity = rep(mats, each = tt),
      yield = as.vector(y),
      curve = nm,
      extrapolated = as.vector(flag),
      stringsAsFactors = FALSE
    )
  })
  # Unfitted dates stay in as missing rows. Dropping them would hand atsm() a
  # panel that silently skips a month, which its factor VAR would read as
  # consecutive; a missing row it refuses outright.
  long <- do.call(rbind, long)

  new_yield_panel_from_long(
    long,
    instrument = stats::setNames(meta$instrument, meta$curve),
    issuer = stats::setNames(meta$issuer, meta$curve)
  )
}

#' @keywords internal
#' @noRd
predict_maturities <- function(maturities, name, observed) {
  if (is.null(maturities)) return(seq_len(floor(max(observed))))

  if (is.list(maturities)) {
    if (!name %in% names(maturities)) {
      stop("`maturities` is a list with no entry for curve '", name, "'.",
           call. = FALSE)
    }
    maturities <- maturities[[name]]
  }
  if (!is.numeric(maturities) || !length(maturities) || anyNA(maturities) ||
      any(!is.finite(maturities)) || any(maturities <= 0)) {
    stop("`maturities` must be positive numbers of months.", call. = FALSE)
  }
  sort(unique(as.numeric(maturities)))
}

#' Parameters of a fitted curve
#'
#' @param object A `svensson_fit`.
#' @param curve Curve name. Defaults to the first.
#' @param ... Unused.
#'
#' @return A data frame with `date`, `beta0`-`beta3` (decimals), `tau1` and
#'   `tau2` (years), the in-sample `rmse` (decimals) and `n_obs`. It can be
#'   passed straight to [svensson_curve()], which returns decimal yields from
#'   it. Nelson-Siegel fits have `beta3 = 0` and `tau2 = NA`.
#'
#' @export
coef.svensson_fit <- function(object, curve = NULL, ...) {
  curve <- curve %||% object$meta$curve[1L]
  cf <- object$curves[[curve]]
  if (is.null(cf)) {
    stop("No curve named '", curve, "'. Available: ",
         paste(names(object$curves), collapse = ", "), call. = FALSE)
  }
  data.frame(date = object$dates, cf$params, n_obs = cf$n_obs,
             row.names = NULL)
}

#' @export
print.svensson_fit <- function(x, ...) {
  cat("<svensson_fit>\n")
  cat("  dates      : ", length(x$dates), " (", x$frequency, ", ",
      format(min(x$dates)), " to ", format(x$dates[length(x$dates)]), ")\n",
      sep = "")

  bp <- function(v) sprintf("%.1fbp", v * 1e4)

  for (nm in names(x$curves)) {
    cf <- x$curves[[nm]]
    p <- cf$params
    ok <- !is.na(p[, "rmse"])
    label <- if (cf$model == "svensson") "Svensson" else "Nelson-Siegel"

    cat("\n  curve '", nm, "'\n", sep = "")
    cat("    model      : ", label, ", decays ",
        if (is.null(cf$tau)) {
          sprintf("estimated per date within [%g, %g] years",
                  x$tau_range[1L], x$tau_range[2L])
        } else {
          paste0("fixed at ", paste(format(cf$tau), collapse = " and "),
                 " years")
        }, "\n", sep = "")
    par <- identical(cf$yield_type, "par")
    cat("    yields     : ", if (par) {
      paste0("par, ", coupon_label(cf$coupon_frequency), " coupons -- ",
             "fitted by pricing par bonds off the curve")
    } else {
      "zero-coupon"
    }, "\n", sep = "")
    cat("    tenors     : ", length(cf$maturities), " (",
        paste(sort(cf$maturities), collapse = ", "), " months)\n",
        sep = "")
    if (any(ok)) {
      worst <- which.max(p[, "rmse"])
      cat("    fit        : RMSE median ", bp(stats::median(p[ok, "rmse"])),
          ", worst ", bp(p[worst, "rmse"]), " on ",
          format(x$dates[worst]), if (par) ", in par yields", "\n", sep = "")
    }
    if (is.null(cf$tau) && any(ok)) {
      cat("    decays     : tau1 median ",
          format(stats::median(p[ok, "tau1"]), digits = 3), "y",
          if (cf$model == "svensson") {
            paste0(", tau2 median ",
                   format(stats::median(p[ok, "tau2"]), digits = 3), "y")
          }, "\n", sep = "")
    }
    if (length(cf$dropped)) {
      cat("    unfitted   : ", length(cf$dropped), " date(s) with fewer than ",
          cf$min_obs, " tenors\n", sep = "")
    }
    cat("    factors    : ", if (is.null(cf$tau)) {
      "not limited by the curve model"
    } else {
      k <- if (cf$model == "svensson") 4L else 3L
      paste0("at most ", k, " -- fixed decays give every date the same ", k,
             " loadings")
    }, "\n", sep = "")

    if (!is.null(cf$check) && nrow(cf$check)) {
      cat("    hold-out   : each interior tenor dropped, every date refitted ",
          "and the tenor predicted\n", sep = "")
      for (i in seq_len(nrow(cf$check))) {
        r <- cf$check[i, ]
        cat(sprintf("      %5gm   RMSE %7s   max %7s\n", r$maturity,
                    bp(r$rmse), bp(r$max_abs)))
      }
    }
  }
  invisible(x)
}

#' @keywords internal
#' @noRd
coupon_label <- function(freq) {
  switch(as.character(freq), "1" = "annual", "2" = "semi-annual",
         "4" = "quarterly", "12" = "monthly", paste(freq, "a year"))
}
