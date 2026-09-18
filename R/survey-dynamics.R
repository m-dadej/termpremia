# Survey-augmented real-world dynamics.
#
# The problem. ACM uses yield information and nothing else. Cohen, Hordahl and
# Xia put the consequence plainly: such models are "prone to overreacting to
# changes in the general level of interest rates... they may tend to interpret
# a change in interest rates as evidence that the steady-state interest rate
# has changed correspondingly," producing "exaggerated movements in
# distant-horizon interest rate projections." This is the shifting-endpoint
# problem, and it is the core model's main known weakness rather than a detail.
#
# The fix. Professional forecasters publish what they expect the short rate to
# be at horizons the model also has an opinion about. Disciplining the factor
# VAR to reproduce those forecasts anchors the distant-horizon projections to
# something outside the yield curve. Kim and Orphanides (2012) and Kim and
# Wright (2005) do this inside a Kalman filter by treating surveys as noisy
# observations of model-implied expectations; here the same idea is applied as
# a penalty on the step-1 regression, which keeps steps 2-3 closed-form.
#
# This is emphatically "Kim-Wright style" and not a replication of Kim-Wright.
# KW use Blue Chip Financial Forecasts and estimate by maximum likelihood; both
# differ here, and the licensing section of the README explains why the first
# one has to.
#
# Units. Monthly rate units throughout, as everywhere else in the estimator.

#' Control parameters for survey-augmented factor dynamics
#'
#' Settings for the survey-disciplined VAR used when `atsm()` is called with
#' `p_dynamics = "survey"`.
#'
#' @section What is being fitted:
#' The factor VAR is chosen to minimise
#'
#' \deqn{\sum_t v_t' \Sigma^{-1} v_t + \sigma_s^{-2} \sum_i (s_i -
#'   \hat{s}_i(\mu, \Phi))^2}
#'
#' where \eqn{v_t} are the usual VAR innovations and \eqn{s_i} is survey
#' observation `i`. The first term is the ordinary least squares criterion, in
#' the metric of the innovation covariance; the second is the cost of
#' disagreeing with the forecasters. `sd` is \eqn{\sigma_s}, so it is not an
#' abstract tuning weight: it is how far the model is allowed to sit from a
#' survey before the disagreement counts as much as a one-standard-deviation
#' VAR innovation. With `sd` very large the criterion reduces to plain OLS,
#' which is a property the tests check.
#'
#' @section What the model says a survey should be:
#' A forecast is described by three things: how far ahead it looks, the tenor
#' of the instrument being forecast, and how many months the forecast is
#' averaged over. The SPF's `TBILL` is the *quarterly average* of the
#' *three-month* bill, so `tenor = 3` and `average_months = 3`; its `BILL10` is
#' the *ten-year average* of the same bill, so `average_months = 120`.
#'
#' Because a useful survey set mixes horizons like those two, `tenor` and
#' `average_months` may also be given as **columns of the survey data frame**,
#' in which case they are used per forecast and these arguments serve only as
#' defaults for rows that lack them. [spf_tbill()] and [spf_bill10()] both
#' return the columns, so their output can simply be `rbind()`-ed.
#'
#' Under the expectations hypothesis the `n`-month rate is the average of the
#' next `n` expected one-month rates, so the model's counterpart to such a
#' forecast is a weighted average of \eqn{E_t[r_{t+p}]} over a window of
#' `average_months + tenor - 1` months, with the weights given by how many ways
#' each month is covered. For the SPF's quarterly forecasts those weights are
#' `(1, 2, 3, 2, 1) / 9`.
#'
#' The approximation this makes is to ignore the term premium on the forecast
#' instrument itself: a forecaster predicting the three-month bill is
#' predicting a rate that contains a three-month term premium, while the model
#' quantity is a pure expectation. At a three-month tenor that premium is a few
#' basis points against a survey disagreement measured in tens, but it is an
#' approximation and not a derivation.
#'
#' @section Short-horizon surveys alone are not enough:
#' This is the finding most likely to catch someone out. Disciplining the model
#' with the SPF's quarterly forecasts -- one to four quarters ahead -- does
#' improve the near-term fit, from 49 to 28 basis points on US data. It also
#' makes the shifting-endpoint problem **worse**: to track the near horizons
#' the VAR becomes more persistent, and the expected short rate ten years out
#' swings *more* than under plain OLS, not less.
#'
#' The reason is simple once seen. Forecasts that reach twelve months say
#' almost nothing about where a model thinks rates settle, and fitting them
#' better is not the same as anchoring the far end. Anchoring the far end needs
#' a far-end forecast. [spf_bill10()] is the only free one -- ten-year average
#' bill expectations, asked once a year -- and 35 observations of it do more
#' for the endpoint than 700 quarterly ones. Used together, near and far, both
#' horizons are fitted and the endpoint is anchored; that is the configuration
#' worth using, and the one Kim and Wright approximate with proprietary Blue
#' Chip long-horizon forecasts.
#'
#' @section Why the intercept is free here and not under BRW:
#' The Bauer-Rudebusch-Wu correction deliberately holds the unconditional mean
#' of the factors fixed, because small-sample bias is about persistence and not
#' about the level things revert to. Surveys are the opposite case: the
#' steady-state level is exactly what they are being asked to discipline, so
#' `mu` is left free to move.
#'
#' @param units Units of the survey `value` column, interpreted as for
#'   [yield_panel()].
#' @param tenor Maturity, in months, of the instrument being forecast. `1`
#'   means the forecast is of the one-month rate itself. Overridden by a
#'   `tenor` column in the survey data frame.
#' @param average_months Number of months the forecast is averaged over. `1` is
#'   a point forecast for a single month; `3` is a quarterly average; `120` is
#'   a ten-year average. Overridden by an `average_months` column in the survey
#'   data frame.
#' @param sd Assumed standard deviation of the gap between a survey forecast
#'   and the model's expectation, as an annualised decimal. Smaller values
#'   trust the surveys more. The default of 30 basis points is in the range Kim
#'   and Orphanides report, and is roughly self-consistent on US data: fitted
#'   against the SPF the model cannot get below about 27 basis points of
#'   disagreement however hard it is pushed, so assuming much more than that
#'   discards information and assuming much less buys nothing. The results are
#'   not knife-edge -- anything from 10 to 50 basis points costs the VAR's own
#'   one-step fit under 3% -- but they are not invariant either, and a serious
#'   application should report the sensitivity.
#' @param max_gap Largest number of days between a survey date and the panel
#'   observation it is attached to. Surveys run on their own calendar and will
#'   rarely land on a month end.
#' @param max_iter Maximum optimiser iterations.
#' @param reltol Relative convergence tolerance passed to [stats::optim()].
#'
#' @return A list of validated control parameters.
#'
#' @references
#' Kim, D. H. and A. Orphanides (2012). "Term structure estimation with survey
#' data on interest rate forecasts." *Journal of Financial and Quantitative
#' Analysis* 47(1), 241-272.
#'
#' Kim, D. H. and J. H. Wright (2005). "An arbitrage-free three-factor term
#' structure model and the recent behavior of long-term yields and
#' distant-horizon forward rates." *Finance and Economics Discussion Series*
#' 2005-33, Federal Reserve Board.
#'
#' @seealso [atsm()], [spf_tbill()]
#'
#' @examples
#' # Defaults describe a forecast of the one-month rate for a single month.
#' survey_control()
#'
#' # Trust the surveys more than the default 50 basis points.
#' survey_control(sd = 0.002)
#'
#' @export
survey_control <- function(units = c("auto", "percent", "decimal"),
                           tenor = 1L,
                           average_months = 1L,
                           sd = 0.003,
                           max_gap = 45L,
                           max_iter = 500L,
                           reltol = 1e-10) {
  units <- match.arg(units)

  check <- function(value, name, min, max, integer = FALSE) {
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value)) {
      stop("`", name, "` must be a single finite number.", call. = FALSE)
    }
    if (value < min || value > max) {
      stop("`", name, "` must be between ", min, " and ", max, ".",
           call. = FALSE)
    }
    if (integer) as.integer(value) else as.numeric(value)
  }

  list(
    units = units,
    tenor = check(tenor, "tenor", 1, 600, TRUE),
    average_months = check(average_months, "average_months", 1, 600, TRUE),
    sd = check(sd, "sd", .Machine$double.eps, Inf),
    max_gap = check(max_gap, "max_gap", 0, Inf, TRUE),
    max_iter = check(max_iter, "max_iter", 1, Inf, TRUE),
    reltol = check(reltol, "reltol", .Machine$double.eps, 1)
  )
}


#' Weights on expected future short rates implied by a survey's definition
#'
#' A forecast of the average, over `average_months` months, of an
#' `tenor`-month rate is a weighted average of one-month expectations over
#' `average_months + tenor - 1` months. Each month of the window is counted
#' once for every (forecast month, tenor offset) pair that reaches it, which
#' gives a trapezoidal weight profile.
#'
#' @param tenor Maturity of the forecast instrument, in months.
#' @param average_months Months the forecast is averaged over.
#'
#' @return A numeric vector of weights summing to one.
#' @keywords internal
#' @noRd
survey_weights <- function(tenor, average_months) {
  counts <- numeric(average_months + tenor - 1L)
  for (u in seq_len(average_months) - 1L) {
    span <- u + seq_len(tenor)
    counts[span] <- counts[span] + 1
  }
  counts / (average_months * tenor)
}


#' Whole months between two dates
#'
#' Calendar months, not days divided by 30. A survey taken in mid-February that
#' forecasts the second quarter is looking two months ahead in the sense the
#' model needs, whatever the day of the month.
#'
#' @param from,to Dates.
#' @return Integer vector of month differences.
#' @keywords internal
#' @noRd
months_between <- function(from, to) {
  f <- as.POSIXlt(from)
  t <- as.POSIXlt(to)
  (t$year - f$year) * 12L + (t$mon - f$mon)
}


#' Attach survey forecasts to panel dates and horizons
#'
#' Turns a survey data frame into the three things the objective needs: which
#' panel observation each forecast was made at, how many months ahead its
#' target window starts, and the forecast value in monthly units.
#'
#' Forecasts whose target window has already begun by the time they are matched
#' to a panel date are dropped. This is not fussiness: the SPF is collected in
#' the middle of the quarter, so the current-quarter forecast is partly a
#' statement about months that have already happened, and the model has no
#' expectation to compare it against. It also carries almost no forecasting
#' content -- against realised data it is accurate to about 16 basis points.
#'
#' @param survey Data frame with `date`, `target_date` and `value`.
#' @param dates The panel's dates.
#' @param control Output of `survey_control()`.
#'
#' @return A list with `index`, `horizon`, `value` and `dropped` counts.
#' @keywords internal
#' @noRd
survey_targets <- function(survey, dates, control) {
  if (!is.data.frame(survey)) {
    stop("`survey` must be a data frame with `date`, `target_date` and ",
         "`value` columns.", call. = FALSE)
  }
  missing_cols <- setdiff(c("date", "target_date", "value"), names(survey))
  if (length(missing_cols)) {
    stop("`survey` is missing column(s): ", paste(missing_cols, collapse = ", "),
         ". See ?atsm for the expected shape.", call. = FALSE)
  }

  s_date <- as_date_strict(survey$date, arg = "survey$date")
  t_date <- as_date_strict(survey$target_date, arg = "survey$target_date")
  s_val <- as_yield_decimal(survey$value, control$units, arg = "survey$value")

  # Window shape may vary by forecast, so that short-horizon and long-horizon
  # surveys can be fitted together. The control supplies the default.
  spec_col <- function(name, default) {
    if (is.null(survey[[name]])) return(rep(default, nrow(survey)))
    v <- survey[[name]]
    if (!is.numeric(v) || anyNA(v) || any(v < 1) || any(v != as.integer(v))) {
      stop("`survey$", name, "` must be whole numbers of months, at least 1.",
           call. = FALSE)
    }
    as.integer(v)
  }
  tenor <- spec_col("tenor", control$tenor)
  avg <- spec_col("average_months", control$average_months)

  keep <- !is.na(s_val) & !is.na(s_date) & !is.na(t_date)
  s_date <- s_date[keep]
  t_date <- t_date[keep]
  s_val <- s_val[keep]
  tenor <- tenor[keep]
  avg <- avg[keep]
  n_missing <- sum(!keep)

  if (!length(s_date)) {
    stop("`survey` has no usable rows.", call. = FALSE)
  }

  # Nearest panel observation to the date the forecast was made.
  idx <- vapply(s_date, function(d) {
    gaps <- abs(as.numeric(dates - d))
    if (min(gaps) > control$max_gap) NA_integer_ else which.min(gaps)
  }, integer(1))

  n_unmatched <- sum(is.na(idx))
  ok <- !is.na(idx)

  # Horizon is measured from the matched panel date, not from the survey date,
  # so that the few weeks between the two do not quietly shift every forecast.
  horizon <- rep(NA_integer_, length(idx))
  horizon[ok] <- months_between(dates[idx[ok]], t_date[ok])

  past <- ok & horizon < 1L
  n_past <- sum(past)
  use <- ok & !past

  if (!any(use)) {
    stop(
      "No survey forecast could be used: ", n_unmatched, " had no panel ",
      "observation within ", control$max_gap, " days and ", n_past,
      " referred to a period that had already begun. Check that `survey` and ",
      "the panel overlap, and that `target_date` is the first month the ",
      "forecast refers to.",
      call. = FALSE
    )
  }

  list(
    index = idx[use],
    horizon = horizon[use],
    value = s_val[use] / 12,          # annualised decimal -> monthly units
    tenor = tenor[use],
    average_months = avg[use],
    survey_date = s_date[use],
    n_used = sum(use),
    n_missing = n_missing,
    n_unmatched = n_unmatched,
    n_past = n_past
  )
}


#' Precompute everything about the survey fit that does not depend on the VAR
#'
#' The optimiser evaluates the objective a few thousand times, so the weights,
#' the horizon each one applies to, and the set of dates that need iterating
#' forward at all are built once here.
#'
#' Only the panel dates a survey was actually matched to are carried forward.
#' Quarterly surveys against a monthly panel touch roughly a third of the
#' sample, and long-horizon forecasts a great deal less, so this is most of the
#' cost of the fit.
#'
#' @param targets Output of `survey_targets()`.
#' @return A list with `rows`, `weights` and `p_max`.
#' @keywords internal
#' @noRd
survey_design <- function(targets) {
  rows <- sort(unique(targets$index))
  row_of <- match(targets$index, rows)

  span <- targets$average_months + targets$tenor - 1L
  p_max <- max(targets$horizon + span - 1L)

  # weights[i, p + 1] is the weight forecast i puts on E_t[r_{t+p}].
  weights <- matrix(0, nrow = length(targets$index), ncol = p_max + 1L)
  for (i in seq_along(targets$index)) {
    w <- survey_weights(targets$tenor[i], targets$average_months[i])
    cols <- targets$horizon[i] + seq_along(w)      # horizon h -> column h + 1
    weights[i, cols] <- w
  }

  list(rows = rows, row_of = row_of, weights = weights, p_max = p_max)
}


#' Model-implied counterparts of a set of survey forecasts
#'
#' Iterates the VAR forward to build \eqn{E_t[r_{t+p}]} for every panel date
#' and every horizon any survey needs, then takes the weighted average over
#' each survey's window.
#'
#' The forward iteration is done for all dates at once -- the state is a
#' `T x k` matrix advanced one horizon per matrix product -- because the
#' objective is evaluated a few hundred times by the optimiser and a loop over
#' dates inside a loop over horizons inside the optimiser is three loops too
#' many.
#'
#' @param mu,phi VAR parameters.
#' @param x0 Factor values at the dates surveys were matched to, `n_rows x k`.
#' @param delta0,delta1 Short-rate loadings.
#' @param design Output of `survey_design()`.
#'
#' @return Vector of model-implied forecasts, in monthly units.
#' @keywords internal
#' @noRd
survey_implied <- function(mu, phi, x0, delta0, delta1, design) {
  n_rows <- nrow(x0)
  f <- matrix(0, nrow = n_rows, ncol = design$p_max + 1L)
  f[, 1L] <- delta0 + drop(x0 %*% delta1)

  ex <- x0
  phi_t <- t(phi)
  mu_rep <- rep(mu, each = n_rows)
  for (p in seq_len(design$p_max)) {
    ex <- ex %*% phi_t + mu_rep
    f[, p + 1L] <- delta0 + drop(ex %*% delta1)
  }

  rowSums(design$weights * f[design$row_of, , drop = FALSE])
}


#' Fit a factor VAR disciplined by survey forecasts
#'
#' @param x `T x k` factor matrix.
#' @param mu0,phi0 OLS starting values.
#' @param sigma Innovation covariance, held fixed.
#' @param delta0,delta1 Short-rate loadings, held fixed.
#' @param targets Output of `survey_targets()`.
#' @param control Output of `survey_control()`.
#'
#' @return A list with the fitted `mu` and `phi` and convergence diagnostics.
#' @keywords internal
#' @noRd
survey_fit_var <- function(x, mu0, phi0, sigma, delta0, delta1, targets,
                           control) {
  tt <- nrow(x)
  k <- ncol(x)
  x_lag <- x[-tt, , drop = FALSE]
  x_led <- x[-1L, , drop = FALSE]
  n_var <- nrow(x_lag)

  sigma_inv <- solve(sigma)
  sd_m <- control$sd / 12                       # annualised decimal -> monthly

  design <- survey_design(targets)
  x0 <- x[design$rows, , drop = FALSE]
  s_obs <- targets$value

  # Parameters are scaled so the optimiser sees comparable magnitudes: in
  # monthly units on principal component factors, mu is of order 1e-5 while
  # the elements of phi are of order 1, and an unscaled search effectively
  # ignores the intercept.
  mu_scale <- max(abs(mu0), .Machine$double.eps^0.5)

  unpack <- function(theta) {
    list(mu = theta[seq_len(k)] * mu_scale,
         phi = matrix(theta[-seq_len(k)], k, k))
  }

  objective <- function(theta) {
    p <- unpack(theta)
    e <- x_led - x_lag %*% t(p$phi) - rep(p$mu, each = n_var)
    var_part <- sum((e %*% sigma_inv) * e)

    s_hat <- survey_implied(p$mu, p$phi, x0, delta0, delta1, design)
    if (!all(is.finite(s_hat))) return(.Machine$double.xmax)
    survey_part <- sum((s_obs - s_hat)^2) / sd_m^2

    (var_part + survey_part) / n_var
  }

  theta0 <- c(mu0 / mu_scale, as.vector(phi0))
  fit <- stats::optim(
    theta0, objective, method = "BFGS",
    control = list(maxit = control$max_iter, reltol = control$reltol)
  )

  out <- unpack(fit$par)

  # Survey fit before and after, in annualised basis points, which is the
  # number that says whether the discipline actually bit.
  rmse <- function(pars) {
    sqrt(mean((s_obs - survey_implied(pars$mu, pars$phi, x0, delta0, delta1,
                                      design))^2)) * 12 * 1e4
  }

  list(
    mu = out$mu,
    phi = out$phi,
    converged = fit$convergence == 0L,
    convergence = fit$convergence,
    iterations = unname(fit$counts[["function"]]),
    objective = fit$value,
    objective_ols = objective(theta0),
    rmse_ols_bp = rmse(list(mu = mu0, phi = phi0)),
    rmse_survey_bp = rmse(out)
  )
}


#' Apply survey discipline to a fitted set of ACM parameters
#'
#' @param pars Output of `acm_three_step()`.
#' @param x `T x k` factor matrix.
#' @param survey Survey data frame.
#' @param dates Panel dates.
#' @param control Output of `survey_control()`.
#'
#' @return A list with the updated `pars` and a `diagnostics` list.
#' @keywords internal
#' @noRd
acm_survey_correct <- function(pars, x, survey, dates, control) {
  targets <- survey_targets(survey, dates, control)

  if (targets$n_unmatched || targets$n_past || targets$n_missing) {
    message(
      "Survey: using ", targets$n_used, " forecast(s); dropped ",
      targets$n_missing, " missing, ", targets$n_unmatched,
      " outside the panel's dates and ", targets$n_past,
      " referring to a period already under way."
    )
  }

  fit <- survey_fit_var(x, pars$mu, pars$phi, pars$sigma,
                        pars$delta0, pars$delta1, targets, control)

  if (!fit$converged) {
    warning(
      "The survey-augmented VAR did not converge (optim code ",
      fit$convergence, "). Raise `max_iter` in survey_control(), or loosen ",
      "`reltol`. The reported dynamics are the best point reached.",
      call. = FALSE
    )
  }

  phi_ols <- pars$phi
  mu_ols <- pars$mu
  pars <- adopt_p_dynamics(pars, fit$mu, fit$phi)

  list(
    pars = pars,
    diagnostics = list(
      n_forecasts = targets$n_used,
      horizons = sort(unique(targets$horizon)),
      first_survey = min(targets$survey_date),
      last_survey = max(targets$survey_date),
      dropped = c(missing = targets$n_missing,
                  unmatched = targets$n_unmatched,
                  already_started = targets$n_past),
      rmse_bp = c(ols = fit$rmse_ols_bp, survey = fit$rmse_survey_bp),
      spectral_radius = c(ols = spectral_radius(phi_ols),
                          survey = spectral_radius(fit$phi)),
      phi_ols = phi_ols,
      mu_ols = mu_ols,
      converged = fit$converged,
      iterations = fit$iterations,
      control = control
    )
  )
}
