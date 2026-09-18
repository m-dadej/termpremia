# Bauer, Rudebusch and Wu (2012) small-sample bias correction for the
# real-world factor dynamics.
#
# OLS estimates of an autoregressive coefficient are biased towards zero in
# finite samples (Kendall 1954; Pope 1990). Interest rate factors are highly
# persistent and term structure samples are short, so the bias is large exactly
# where it hurts most: a Phi that mean-reverts too fast makes expected future
# short rates converge to their unconditional mean too quickly, the
# expectations component of a long yield goes flat, and everything the model
# cannot explain is swept into the term premium. BRW show this accounts for a
# substantial part of the variation conventionally attributed to term premia.
#
# What is corrected, and what is not. The bias is a property of the P-dynamics
# -- the time series regression in step 1. The risk-adjusted (Q) dynamics are
# identified by the cross-section of yields, which is large, precisely fitted
# and essentially unaffected by the persistence of the VAR; BRW make this point
# explicitly and it is why their bias-corrected models fit the yield curve
# identically to the uncorrected ones. This implementation preserves that
# property exactly: see `acm_brw_correct()`.
#
# Units follow the rest of the package: everything here is in monthly rate
# units, though nothing in the bias correction depends on the scale.

#' Control parameters for the Bauer-Rudebusch-Wu bias correction
#'
#' Settings for the bootstrap bias correction applied when `atsm()` is called
#' with `p_dynamics = "brw"`.
#'
#' @section Method:
#' Write the factor VAR in deviations from its unconditional mean,
#' \eqn{Z_{t+1} = \Phi Z_t + v_{t+1}}, and let \eqn{\hat\Phi} be the OLS
#' estimate. Simulating from a known \eqn{\Phi} and re-estimating by OLS gives
#' the mean of the OLS estimator at that parameter, \eqn{m(\Phi)}, and hence
#' the bias \eqn{m(\Phi) - \Phi}.
#'
#' `method = "indirect"` is indirect inference: find the \eqn{\tilde\Phi} whose
#' OLS estimates are centred on what was actually observed, that is solve
#' \eqn{m(\tilde\Phi) = \hat\Phi}. The fixed point is found by the iteration
#' \eqn{\tilde\Phi_{j+1} = \tilde\Phi_j + \textrm{step} \times (\hat\Phi -
#' m(\tilde\Phi_j))}, started at \eqn{\hat\Phi}.
#'
#' `method = "bootstrap"` is the cheaper one-step correction
#' \eqn{\tilde\Phi = \hat\Phi - (m(\hat\Phi) - \hat\Phi)}, which evaluates the
#' bias only at the OLS estimate. It is exactly one iteration of the above with
#' a unit step, and it under-corrects when the bias function is steep -- which
#' it is near a unit root, the case the correction exists for.
#'
#' @section Why the random draws are reused:
#' The innovation draws and the starting values are made once and held fixed
#' across iterations. With fresh draws each time, \eqn{m(\cdot)} is a noisy
#' function and the iteration rattles around the fixed point at the amplitude
#' of the Monte Carlo error instead of converging. Holding them fixed makes the
#' simulated map deterministic, so `tol` means what it says. It also means the
#' answer depends on `seed`, which is why one is set by default.
#'
#' @param method `"indirect"` for indirect inference (the default and BRW's
#'   preferred estimator) or `"bootstrap"` for the one-step correction.
#' @param replications Number of bootstrap paths used to evaluate the mean of
#'   the OLS estimator. The default of 1000 gives a Monte Carlo standard error
#'   on each element of \eqn{m(\Phi)} of roughly 3% of its own sampling
#'   standard deviation.
#' @param step Damping factor on the indirect inference update, in `(0, 1]`.
#'   Lower is slower but more stable.
#' @param tol Convergence tolerance on the largest absolute element of
#'   \eqn{\hat\Phi - m(\tilde\Phi_j)}.
#' @param max_iter Maximum number of iterations. Ignored by
#'   `method = "bootstrap"`.
#' @param rho_max Largest permitted spectral radius of the corrected \eqn{\Phi}.
#'   The default sits just inside the unit circle: it allows persistence
#'   indistinguishable from a random walk while keeping the unconditional mean
#'   -- and so the long-horizon expected short rate -- defined.
#' @param block_length Length of the moving blocks used to resample VAR
#'   innovations. `1` is an i.i.d. residual bootstrap. Longer blocks preserve
#'   dependence in the innovations, such as volatility clustering, at the cost
#'   of a coarser resample.
#' @param seed Integer seed for the bootstrap, or `NA` to use the ambient
#'   random number stream. The default fixes the seed so that two calls on the
#'   same data give the same answer; the caller's RNG state is saved and
#'   restored either way, so fitting a model never disturbs a simulation
#'   running around it.
#'
#' @return A list of validated control parameters.
#'
#' @references
#' Bauer, M. D., G. D. Rudebusch and J. C. Wu (2012). "Correcting estimation
#' bias in dynamic term structure models." *Journal of Business & Economic
#' Statistics* 30(3), 454-467.
#'
#' @seealso [atsm()]
#'
#' @examples
#' # Fewer replications for a quick look; the default is 1000.
#' brw_control(replications = 200)
#'
#' @export
brw_control <- function(method = c("indirect", "bootstrap"),
                        replications = 1000L,
                        step = 0.5,
                        tol = 1e-5,
                        max_iter = 50L,
                        rho_max = 0.9999,
                        block_length = 1L,
                        seed = 20120401L) {
  method <- match.arg(method)

  check_scalar <- function(value, name, min, max, integer = FALSE) {
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value)) {
      stop("`", name, "` must be a single finite number.", call. = FALSE)
    }
    if (value < min || value > max) {
      stop("`", name, "` must be between ", min, " and ", max, ".",
           call. = FALSE)
    }
    if (integer) as.integer(value) else as.numeric(value)
  }

  if (length(seed) != 1L || (!is.na(seed) && !is.numeric(seed))) {
    stop("`seed` must be a single number or NA.", call. = FALSE)
  }

  list(
    method = method,
    replications = check_scalar(replications, "replications", 2, Inf, TRUE),
    step = check_scalar(step, "step", .Machine$double.eps, 1),
    tol = check_scalar(tol, "tol", .Machine$double.eps, Inf),
    max_iter = check_scalar(max_iter, "max_iter", 1, Inf, TRUE),
    rho_max = check_scalar(rho_max, "rho_max", .Machine$double.eps, 1),
    block_length = check_scalar(block_length, "block_length", 1, Inf, TRUE),
    seed = if (is.na(seed)) NA_integer_ else as.integer(seed)
  )
}


#' Draw bootstrap indices into the rows of a residual matrix
#'
#' Returns an `n_draw x m` integer matrix, one column per bootstrap path.
#' With `block_length = 1` this is an i.i.d. resample. With longer blocks it is
#' a moving-block bootstrap: block starting points are drawn uniformly and each
#' expands into a run of consecutive residuals, so dependence within a block
#' survives the resample.
#'
#' @param n_v Number of residual rows available.
#' @param n_draw Number of innovations needed per path.
#' @param m Number of paths.
#' @param block_length Block length.
#'
#' @return An `n_draw x m` integer matrix.
#' @keywords internal
#' @noRd
brw_draw_index <- function(n_v, n_draw, m, block_length) {
  if (block_length <= 1L) {
    return(matrix(sample.int(n_v, n_draw * m, replace = TRUE), nrow = n_draw))
  }
  if (block_length > n_v) {
    stop("`block_length` (", block_length, ") exceeds the number of VAR ",
         "residuals (", n_v, ").", call. = FALSE)
  }

  n_block <- ceiling(n_draw / block_length)
  starts <- matrix(sample.int(n_v - block_length + 1L, n_block * m,
                              replace = TRUE),
                   nrow = n_block)

  # Expand each start into block_length consecutive rows. The offset vector
  # recycles down the columns, which is what we want: every path uses the same
  # within-block offsets against its own starts.
  idx <- starts[rep(seq_len(n_block), each = block_length), , drop = FALSE] +
    rep(seq_len(block_length) - 1L, times = n_block)

  idx[seq_len(n_draw), , drop = FALSE]
}


#' Mean of the OLS estimator under a candidate Phi
#'
#' Simulates `m` paths of the mean-zero VAR \eqn{Z_{t+1} = \Phi Z_t + v_{t+1}}
#' and returns the average of the OLS estimates of `Phi` across them.
#'
#' The paths are simulated jointly rather than one at a time: the state is an
#' `m x k` matrix advanced by a single matrix product per period, which turns
#' `m * (tt - 1)` interpreted iterations into `tt - 1` of them. At the default
#' settings this is the difference between seconds and minutes.
#'
#' Two details matter for the answer rather than the speed. Each path starts at
#' a randomly drawn observed state, so the initial condition is representative
#' of the data even when the process is close to a unit root and a burn-in from
#' the unconditional mean would be both slow and unrepresentative. And the OLS
#' regression is run *with* an intercept, because the estimate being corrected
#' was.
#'
#' @param phi Candidate `k x k` autoregressive matrix.
#' @param z0 `m x k` matrix of starting states, in deviations from the mean.
#' @param resid `n_v x k` matrix of VAR innovations to resample.
#' @param idx `(tt - 1) x m` matrix of resampling indices.
#' @param tt Length of each simulated path.
#'
#' @return The `k x k` average OLS estimate.
#' @keywords internal
#' @noRd
brw_ols_mean <- function(phi, z0, resid, idx, tt) {
  m <- nrow(z0)
  k <- ncol(z0)

  # Paths are stored with the replication index last so that each path is a
  # contiguous slice when the regressions are run.
  sim <- array(0, dim = c(tt, k, m))
  s <- z0
  sim[1L, , ] <- t(s)

  phi_t <- t(phi)
  for (t in seq_len(tt - 1L)) {
    s <- s %*% phi_t + resid[idx[t, ], , drop = FALSE]
    sim[t + 1L, , ] <- t(s)
  }

  if (!all(is.finite(sim))) {
    stop(
      "The bias correction simulated a divergent VAR. The candidate ",
      "autoregressive matrix has spectral radius ",
      format(spectral_radius(phi), digits = 4),
      ", which overflows over a path of ", tt, " periods. This usually means ",
      "the factor VAR is too poorly identified for the correction to be ",
      "meaningful; try fewer factors.",
      call. = FALSE
    )
  }

  acc <- matrix(0, nrow = k, ncol = k)
  for (b in seq_len(m)) {
    zb <- matrix(sim[, , b], nrow = tt, ncol = k)
    z1 <- cbind(1, zb[-tt, , drop = FALSE])
    cf <- qr.solve(crossprod(z1), crossprod(z1, zb[-1L, , drop = FALSE]))
    acc <- acc + t(cf[-1L, , drop = FALSE])
  }

  acc / m
}


#' Shrink a corrected Phi back towards OLS until it is stationary
#'
#' The bias correction raises persistence, and on a persistent sample it can
#' raise it past a unit root -- at which point expected short rates never
#' converge and the expectations component of a long yield is meaningless. BRW
#' handle this by shrinking along the segment between the OLS and corrected
#' estimates, \eqn{\Phi(\delta) = \hat\Phi + \delta(\tilde\Phi - \hat\Phi)},
#' taking the largest \eqn{\delta} that keeps the process stationary. It is a
#' line search rather than a solve because the spectral radius is not monotone
#' in `delta` in general.
#'
#' If OLS is *already* explosive there is no `delta` that helps and the
#' correction cannot be applied at all; the caller is told.
#'
#' @param phi_hat The OLS estimate.
#' @param phi_tilde The bias-corrected estimate.
#' @param rho_max Largest permitted spectral radius.
#' @param n_grid Number of points on the `delta` grid.
#'
#' @return A list with `phi`, `delta` and `rho`.
#' @keywords internal
#' @noRd
brw_shrink_stationary <- function(phi_hat, phi_tilde, rho_max, n_grid = 1000L) {
  rho <- spectral_radius(phi_tilde)
  if (rho <= rho_max) {
    return(list(phi = phi_tilde, delta = 1, rho = rho))
  }

  rho_ols <- spectral_radius(phi_hat)
  if (rho_ols > rho_max) {
    return(list(phi = phi_hat, delta = 0, rho = rho_ols))
  }

  d <- phi_tilde - phi_hat
  for (delta in seq(1, 0, length.out = n_grid + 1L)) {
    cand <- phi_hat + delta * d
    r <- spectral_radius(cand)
    if (r <= rho_max) return(list(phi = cand, delta = delta, rho = r))
  }

  list(phi = phi_hat, delta = 0, rho = rho_ols)
}


#' Bias-correct the autoregressive matrix of a factor VAR
#'
#' @param x `T x k` factor matrix.
#' @param mu Length-`k` OLS intercept.
#' @param phi `k x k` OLS autoregressive matrix.
#' @param resid `(T-1) x k` OLS residuals.
#' @param control Output of `brw_control()`.
#'
#' @return A list with the corrected `phi`, the anchor mean `x_bar`, and
#'   convergence diagnostics.
#' @keywords internal
#' @noRd
brw_bias_correct <- function(x, mu, phi, resid, control) {
  tt <- nrow(x)
  k <- ncol(x)
  rho_ols <- spectral_radius(phi)

  # The correction is about persistence, not about the level the factors
  # revert to. Everything is therefore done in deviations from a fixed anchor
  # mean, and the intercept is rebuilt around that anchor afterwards. The
  # natural anchor is the unconditional mean the OLS fit itself implies,
  # (I - Phi)^-1 mu, which exists only if OLS is stationary; the sample mean
  # is the fallback. (With principal component factors the sample mean is zero
  # by construction, so the two nearly coincide anyway.)
  x_bar <- if (rho_ols < 1) {
    drop(solve(diag(k) - phi, mu))
  } else {
    colMeans(x)
  }
  z <- sweep(x, 2L, x_bar, "-")

  cells <- as.numeric(tt) * k * control$replications
  if (cells > 5e7) {
    stop(
      "The bias correction would need ", format(round(cells * 8 / 2^20)),
      " MB to simulate ", control$replications, " paths of ", tt,
      " observations. Reduce `replications` in brw_control().",
      call. = FALSE
    )
  }

  # Common random numbers -- drawn once, reused at every iteration. See the
  # note in ?brw_control.
  idx <- brw_draw_index(nrow(resid), tt - 1L, control$replications,
                        control$block_length)
  z0 <- z[sample.int(tt, control$replications, replace = TRUE), , drop = FALSE]

  if (control$method == "bootstrap") {
    # One step, unit step size: correct by the bias evaluated at OLS.
    mean_phi <- brw_ols_mean(phi, z0, resid, idx, tt)
    gap <- max(abs(phi - mean_phi))
    phi_tilde <- phi + (phi - mean_phi)
    iterations <- 1L
    converged <- TRUE
  } else {
    phi_tilde <- phi
    converged <- FALSE
    iterations <- 0L
    gap <- NA_real_

    for (j in seq_len(control$max_iter)) {
      iterations <- j
      mean_phi <- brw_ols_mean(phi_tilde, z0, resid, idx, tt)
      d <- phi - mean_phi
      gap <- max(abs(d))
      if (gap < control$tol) {
        converged <- TRUE
        break
      }
      phi_tilde <- phi_tilde + control$step * d
    }
  }

  shrunk <- brw_shrink_stationary(phi, phi_tilde, control$rho_max)

  list(
    phi = shrunk$phi,
    phi_unshrunk = phi_tilde,
    x_bar = x_bar,
    delta = shrunk$delta,
    rho_ols = rho_ols,
    rho_brw = shrunk$rho,
    rho_unshrunk = spectral_radius(phi_tilde),
    iterations = iterations,
    converged = converged,
    gap = gap
  )
}


#' Apply the BRW correction to a fitted set of ACM parameters
#'
#' @section What moves and what does not:
<<<<<<< HEAD
#' The corrected \eqn{\tilde\Phi} goes in through `adopt_p_dynamics()`, so the
#' risk-adjusted dynamics and every fitted yield are left exactly where the
#' cross-section put them and only the expectations component moves. That is
#' precisely BRW's point: the correction reallocates variation between
#' expectations and term premium without touching the fit.
=======
#' The corrected \eqn{\tilde\Phi} replaces \eqn{\hat\Phi} in the P-dynamics,
#' and the prices of risk absorb the change so that the risk-adjusted dynamics
#' are left exactly where the cross-section put them:
#'
#' \deqn{\tilde\lambda_1 = \hat\lambda_1 + (\tilde\Phi - \hat\Phi)
#'   \quad\Longrightarrow\quad \tilde\Phi - \tilde\lambda_1 =
#'   \hat\Phi - \hat\lambda_1}
#'
#' and likewise for the intercepts. Fitted yields are therefore *identical*
#' before and after the correction -- as they should be, since the cross-section
#' of yields contains no information about how fast the factors mean-revert and
#' is fitted to well under two basis points either way. What changes is the
#' expectations component, which is the only part that uses \eqn{\Phi} directly,
#' and hence the term premium that is its residual. That is precisely BRW's
#' point: the correction reallocates variation between expectations and term
#' premium without touching the fit.
>>>>>>> 2ec34793201fcb4b52491e543e62932d42cf2034
#'
#' Two quantities are deliberately *not* re-estimated. The innovation
#' covariance \eqn{\Sigma} and the return-regression coefficients keep their
#' OLS values: recomputing them from bias-corrected residuals would perturb the
#' convexity term and so the fitted yields, destroying the property above, in
#' exchange for a second-order adjustment to a quantity that is estimated
#' precisely. BRW make the same argument for holding the cross-sectional
#' parameters fixed.
#'
#' @param pars Output of `acm_three_step()`.
#' @param x `T x k` factor matrix.
#' @param control Output of `brw_control()`.
#'
#' @return A list with the updated `pars` and a `diagnostics` list.
#' @keywords internal
#' @noRd
acm_brw_correct <- function(pars, x, control) {
  if (!is.na(control$seed)) {
    # Save and restore the caller's random stream: fitting a model should not
    # silently advance somebody else's simulation.
    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      stats::runif(1)
    }
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    set.seed(control$seed)
  }

  bc <- brw_bias_correct(x, pars$mu, pars$phi, pars$var_residuals, control)

  phi_ols <- pars$phi
  mu_ols <- pars$mu
  phi_new <- bc$phi
  # Preserves the VAR's conditional forecast at the anchor mean, and hence the
  # unconditional mean of the factors, so only the speed of mean reversion
  # changes. Reduces to mu exactly when phi_new == phi.
  mu_new <- pars$mu + drop((pars$phi - phi_new) %*% bc$x_bar)

<<<<<<< HEAD
  pars <- adopt_p_dynamics(pars, mu_new, phi_new)
=======
  pars$lambda1 <- pars$lambda1 + (phi_new - pars$phi)
  pars$lambda0 <- pars$lambda0 + (mu_new - pars$mu)
  pars$phi <- phi_new
  pars$mu <- mu_new
>>>>>>> 2ec34793201fcb4b52491e543e62932d42cf2034

  if (!bc$converged) {
    warning(
      "The BRW bias correction did not converge in ", control$max_iter,
      " iterations (largest remaining discrepancy ",
      format(bc$gap, digits = 3), " against a tolerance of ", control$tol,
      "). Raise `max_iter` or lower `step` in brw_control().",
      call. = FALSE
    )
  }
  if (bc$delta == 0 && bc$rho_unshrunk > control$rho_max) {
    warning(
      "The BRW bias correction was not applied: the OLS factor VAR is already ",
      "at or above the persistence ceiling (spectral radius ",
      format(bc$rho_ols, digits = 5), " against rho_max = ", control$rho_max,
      "), and shrinking the corrected estimate back towards it cannot bring ",
      "it inside. The fit is the uncorrected one.",
      call. = FALSE
    )
  }

  list(
    pars = pars,
    diagnostics = list(
      method = control$method,
      phi_ols = phi_ols,
      mu_ols = mu_ols,
      phi_unshrunk = bc$phi_unshrunk,
      unconditional_mean = bc$x_bar,
      spectral_radius = c(ols = bc$rho_ols, corrected = bc$rho_brw,
                          unshrunk = bc$rho_unshrunk),
      shrinkage = bc$delta,
      iterations = bc$iterations,
      converged = bc$converged,
      gap = bc$gap,
      control = control
    )
  )
}


#' Half-life implied by a spectral radius, formatted for printing
#'
#' The number of months a shock to the slowest-decaying factor combination
#' takes to halve. It is the interpretable face of the spectral radius: 0.96
#' and 0.99 look similar and mean 17 months against 69.
#'
#' @param rho A spectral radius.
#' @return A short character label.
#' @keywords internal
#' @noRd
half_life_label <- function(rho) {
  if (!is.finite(rho) || rho <= 0) return("n/a")
  if (rho >= 1) return("infinite")

  h <- log(0.5) / log(rho)
  if (h >= 240) {
    paste0(format(round(h / 12)), "y")
  } else if (h >= 24) {
    paste0(format(round(h / 12, 1)), "y")
  } else {
    paste0(format(round(h)), "m")
  }
}
