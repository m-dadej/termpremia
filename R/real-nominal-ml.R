# The constrained maximum likelihood step of the joint real-nominal model.
#
# Supplementary Appendix Section 1. R/real-nominal.R has the pricing, the
# factors and the closed-form estimator that this starts from; the two are
# kept apart because the closed form is a complete estimator in its own right
# and most of the package's reasoning about it does not involve any of this.
#
# What this step is for. The closed-form estimator fits excess RETURNS.
# Nothing in it ties the LEVEL of either fitted curve to anything, which is
# why the UK replication reproduces the appendix's standard deviations while
# sitting about fourteen basis points off in the mean (HANDOVER 4d). The
# constraints below tie the fitted curves to the observed principal
# components, and that is what pins the levels down.


#' Settings for the maximum likelihood step
#'
#' @param max_outer Maximum augmented-Lagrangian rounds.
#' @param max_inner Maximum `optim()` iterations within a round.
#' @param penalty Initial penalty on constraint violation.
#' @param penalty_factor Multiplier applied to the penalty in a round that
#'   fails to reduce the violation materially.
#' @param tol Constraint violation below which the fit is taken as converged,
#'   in units of a factor's standard deviation.
#' @param reltol Relative convergence tolerance passed to [stats::optim()].
#' @param ndeps Step used for the numerically differenced gradient. The
#'   `optim()` default of 1e-3 is an ABSOLUTE step, which is far too coarse
#'   once the penalty term makes the merit surface steep: the gradient stops
#'   being informative and the inner optimiser reports convergence after a
#'   handful of evaluations without having moved.
#'
#' @return A list of validated settings.
#'
#' @seealso [atsm_real()]
#'
#' @examples
#' # Push harder on the constraints than the default.
#' acmy_control(tol = 1e-10, max_outer = 20)
#'
#' @export
acmy_control <- function(max_outer = 20L,
                         max_inner = 400L,
                         penalty = 1,
                         penalty_factor = 10,
                         tol = 1e-6,
                         reltol = 1e-10,
                         ndeps = 1e-6) {
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
    max_outer = check(max_outer, "max_outer", 1, 1000, TRUE),
    max_inner = check(max_inner, "max_inner", 1, 100000, TRUE),
    penalty = check(penalty, "penalty", .Machine$double.eps, Inf),
    penalty_factor = check(penalty_factor, "penalty_factor", 1, 1e6),
    tol = check(tol, "tol", 0, 1),
    reltol = check(reltol, "reltol", .Machine$double.eps, 1),
    ndeps = check(ndeps, "ndeps", .Machine$double.eps, 0.1)
  )
}


# The free parameters ------------------------------------------------------

#' Pack and unpack the free parameters of the likelihood
#'
#' `Sigma` and the physical dynamics are not among them. The likelihood of
#' Eq. 22 separates into a measurement block and a VAR block, and the VAR
#' block is maximised by ordinary least squares on the factors -- which are
#' observable, so no pricing parameter changes it. That is why the paper can
#' replace both covariances with consistent estimators "at only a minimal loss
#' of asymptotic efficiency", and why `mu`, `Phi` and `Sigma` are carried
#' through this step untouched.
#'
#' @param pars Parameter list.
#' @param k Number of factors.
#' @keywords internal
#' @noRd
acmy_pack <- function(pars, k) {
  c(pars$mu_tilde, as.vector(pars$phi_tilde), pars$delta0, pars$delta1,
    pars$pi0, pars$pi1)
}

#' @param theta Packed vector.
#' @param pi0_fixed Optional value overriding the packed `pi0`.
#' @rdname acmy_pack
#' @keywords internal
#' @noRd
acmy_unpack <- function(theta, k, pi0_fixed = NULL) {
  taken <- 0L
  take <- function(n) {
    out <- theta[taken + seq_len(n)]
    taken <<- taken + n
    out
  }

  mu_tilde <- take(k)
  phi_tilde <- matrix(take(k * k), k, k)
  delta0 <- take(1L)
  delta1 <- take(k)
  pi0 <- take(1L)
  pi1 <- take(k)

  if (!is.null(pi0_fixed)) pi0 <- pi0_fixed

  list(mu_tilde = mu_tilde, phi_tilde = phi_tilde, delta0 = delta0,
       delta1 = delta1, pi0 = pi0, pi1 = pi1)
}

#' Put the parameter blocks on comparable scales
#'
#' In monthly units `delta0` is of order 1e-3 while the entries of `Phi~` are
#' of order 1. A numerically differenced BFGS on the raw vector effectively
#' ignores whichever block is smaller, which here would mean never moving the
#' short-rate and inflation loadings at all. Each block is therefore divided
#' by the magnitude of its own starting values.
#'
#' @keywords internal
#' @noRd
acmy_theta_scale <- function(start, k) {
  block <- rep(seq_len(6L), c(k, k * k, 1L, k, 1L, k))
  scale <- rep(1, length(start))

  for (b in seq_len(6L)) {
    j <- block == b
    s <- max(abs(start[j]), na.rm = TRUE)
    if (is.finite(s) && s > 0) scale[j] <- s
  }
  scale
}


# The estimator ------------------------------------------------------------

#' Constrained maximum likelihood for the joint real-nominal model
#'
#' Maximises the measurement block of Eq. 22 subject to the requirement that
#' factors extracted from the model's own fitted yields equal the observed
#' factors, starting from the closed-form estimates of Section 1.1.
#'
#' @section Why an augmented Lagrangian:
#' The constraints are linear in `A` and `B` and so nonlinear in the
#' parameters, because `A` and `B` come out of the pricing recursion. Base R
#' has no general nonlinear-equality-constrained optimiser and the package has
#' no dependency beyond `stats`, so one is assembled from [stats::optim()]: a
#' sequence of unconstrained problems
#' \deqn{L(\theta) - \lambda'c(\theta) + \tfrac{\rho}{2}\|c(\theta)\|^2}
#' with multipliers updated by \eqn{\lambda \leftarrow \lambda - \rho
#' c(\theta)}. A plain quadratic penalty would need \eqn{\rho\to\infty} and
#' would be hopelessly ill-conditioned long before it got there; the
#' multipliers are what let it converge at a finite penalty.
#'
#' The penalty is raised only in a round that fails to cut the violation by a
#' quarter. Raising it every round conditions the inner problem badly for no
#' gain, which shows up as the inner optimiser stalling while the violation
#' sits still.
#'
#' @param pars Output of `acmy_closed_form()`, the starting value.
#' @param x,rx_nom,rx_real,inflation,r As in `acmy_closed_form()`.
#' @param fac Output of `acmy_factors()`.
#' @param maturities,real_maturities The two yield grids.
#' @param return_maturities,real_return_maturities The two return grids.
#' @param n_max Longest maturity.
#' @param control Output of [acmy_control()].
#' @param fix_pi0 When not `NULL`, `pi0` is held at the value already in
#'   `pars` rather than estimated.
#'
#' @return `pars`, updated, carrying an `ml` element of diagnostics.
#' @keywords internal
#' @noRd
acmy_maximum_likelihood <- function(pars, x, rx_nom, rx_real, inflation, r,
                                    fac, maturities, real_maturities,
                                    return_maturities, real_return_maturities,
                                    n_max, control, fix_pi0 = NULL) {
  tt <- nrow(x)
  k <- ncol(x)
  kk <- fac$k_nominal + fac$k_real

  x_lag <- x[-tt, , drop = FALSE]
  x_led <- x[-1L, , drop = FALSE]
  n_obs <- nrow(x_lag)

  # Real returns enter augmented by realised inflation, exactly as in the
  # closed form: that is what gives the stacked system one common B.
  rr <- cbind(rx_nom, rx_real + inflation[-1L])
  r_obs <- r[-1L]
  pi_obs <- inflation[-1L]
  sigma <- pars$sigma

  idx_nom <- return_maturities - 1L
  idx_real <- real_return_maturities - 1L
  pi0_fixed <- if (is.null(fix_pi0)) NULL else pars$pi0

  # One model evaluation. Both the likelihood and the constraints read from
  # it, so the 120-step recursion runs once per candidate rather than twice.
  evaluate <- function(p) {
    coefs <- acmy_coefficients(n_max, p$mu_tilde, p$phi_tilde, sigma,
                               p$delta0, p$delta1, p$pi0, p$pi1)
    if (!all(is.finite(coefs$b)) || !all(is.finite(coefs$b_real))) {
      return(NULL)
    }

    b_bar <- rbind(coefs$b[idx_nom, , drop = FALSE],
                   sweep(coefs$b_real[idx_real, , drop = FALSE], 2L,
                         p$pi1, "+"))
    quad <- rowSums((b_bar %*% sigma) * b_bar)
    alpha <- -(drop(b_bar %*% p$mu_tilde) + 0.5 * quad)

    w <- x_led - x_lag %*% t(p$phi_tilde)
    eps <- cbind(
      rr - rep(alpha, each = n_obs) - w %*% t(b_bar),
      r_obs - p$delta0 - drop(x_led %*% p$delta1),
      pi_obs - p$pi0 - drop(x_led %*% p$pi1)
    )
    if (!all(is.finite(eps))) return(NULL)

    list(coefs = coefs, eps = eps)
  }

  start <- acmy_pack(pars, k)
  scale <- acmy_theta_scale(start, k)
  at <- function(z) acmy_unpack(z * scale, k, pi0_fixed)

  first <- evaluate(at(start / scale))
  if (is.null(first)) {
    stop("The closed-form estimates do not produce a finite model, so the ",
         "likelihood has nowhere to start. Fit with ",
         "`method = \"closed_form\"` and inspect the result.", call. = FALSE)
  }

  # Sigma_epsilon is held at a consistent estimator, as the paper does, which
  # is what leaves a weighted least-squares problem with no log-determinant.
  sigma_eps <- crossprod(first$eps) / n_obs
  whiten <- sigma_e_whitener(sigma_eps)

  objective <- function(z) {
    ev <- evaluate(at(z))
    if (is.null(ev)) return(NULL)
    list(
      value = sum(whiten(t(ev$eps))^2) / n_obs,
      cons = acmy_constraints(ev$coefs, fac, maturities, real_maturities,
                              k, kk)
    )
  }

  theta <- start / scale
  base <- objective(theta)

  # A note on what NOT to do here. Driving the constraint violation to zero
  # on its own first, as a pure least-squares problem, is tempting and works:
  # it reaches a violation of 1e-7 in one pass. But the feasible set has
  # twelve dimensions, and minimising the violation alone lands on an
  # arbitrary point of it -- one that fits the UK nominal curve at 9bp rather
  # than 6bp and leaves the risk-adjusted dynamics explosive. The subsequent
  # rounds cannot recover, because they start from there. Beginning at the
  # closed-form estimates, which are already a good likelihood point, and
  # letting the penalty pull them onto the feasible set finds a feasible
  # point that is also a good fit.
  lambda <- rep(0, length(base$cons))
  rho <- control$penalty

  best <- list(theta = theta, violation = max(abs(base$cons)),
               value = base$value)
  prev_violation <- best$violation
  trace <- NULL
  converged <- FALSE

  merit <- function(z) {
    o <- objective(z)
    if (is.null(o)) return(.Machine$double.xmax)
    o$value - sum(lambda * o$cons) + 0.5 * rho * sum(o$cons^2)
  }

  for (round in seq_len(control$max_outer)) {
    fit <- stats::optim(theta, merit, method = "BFGS",
                        control = list(maxit = control$max_inner,
                                       reltol = control$reltol,
                                       ndeps = rep(control$ndeps,
                                                   length(theta))))
    theta <- fit$par
    o <- objective(theta)
    if (is.null(o)) break

    violation <- max(abs(o$cons))
    trace <- rbind(trace, data.frame(
      round = round, objective = o$value, violation = violation,
      penalty = rho, code = fit$convergence,
      evaluations = unname(fit$counts[["function"]]),
      stringsAsFactors = FALSE
    ))

    if (violation < best$violation) {
      best <- list(theta = theta, violation = violation, value = o$value)
    }
    if (violation < control$tol) {
      converged <- TRUE
      break
    }

    lambda <- lambda - rho * o$cons
    if (violation > 0.25 * prev_violation) rho <- rho * control$penalty_factor
    prev_violation <- violation
  }

  fitted_pars <- acmy_unpack(best$theta * scale, k, pi0_fixed)

  pars$mu_tilde <- fitted_pars$mu_tilde
  pars$phi_tilde <- fitted_pars$phi_tilde
  pars$delta0 <- fitted_pars$delta0
  pars$delta1 <- fitted_pars$delta1
  pars$pi0 <- fitted_pars$pi0
  pars$pi1 <- fitted_pars$pi1
  pars$lambda0 <- pars$mu - fitted_pars$mu_tilde
  pars$lambda1 <- pars$phi - fitted_pars$phi_tilde
  pars$sigma_eps <- sigma_eps

  pars$ml <- list(
    converged = converged,
    rounds = if (is.null(trace)) 0L else nrow(trace),
    violation = best$violation,
    violation_start = max(abs(base$cons)),
    objective = best$value,
    objective_start = base$value,
    n_constraints = length(base$cons),
    n_parameters = length(start),
    penalty = rho,
    trace = trace
  )
  pars
}
