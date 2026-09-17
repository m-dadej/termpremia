# Bauer-Rudebusch-Wu bias correction.
#
# The bootstrap is the expensive part, so the fits used across several tests
# are built once and reused, at a lower replication count than the default.

brw_panel <- function() {
  yield_panel(gsw_monthly, units = "percent", maturity_unit = "months",
              instrument = "government", issuer = "US")
}

brw_fits <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      panel <- brw_panel()
      cache <<- list(
        ols = suppressWarnings(atsm(panel, n_factors = 5)),
        brw = suppressWarnings(atsm(
          panel, n_factors = 5, p_dynamics = "brw",
          brw_control = brw_control(replications = 200)
        ))
      )
    }
    cache
  }
})

# A small OLS VAR fit, in the form brw_bias_correct() expects.
var_ols <- function(z) {
  tt <- nrow(z)
  z1 <- cbind(1, z[-tt, , drop = FALSE])
  cf <- qr.solve(crossprod(z1), crossprod(z1, z[-1L, , drop = FALSE]))
  list(mu = cf[1L, ], phi = t(cf[-1L, , drop = FALSE]),
       resid = z[-1L, , drop = FALSE] - z1 %*% cf)
}

sim_ar <- function(phi, tt, burn = 300) {
  k <- nrow(phi)
  z <- matrix(0, tt + burn, k)
  for (t in 2:(tt + burn)) z[t, ] <- drop(phi %*% z[t - 1L, ]) + rnorm(k)
  z[-seq_len(burn), , drop = FALSE]
}


# Controls ----------------------------------------------------------------

test_that("brw_control returns validated defaults", {
  ctl <- brw_control()
  expect_equal(ctl$method, "indirect")
  expect_equal(ctl$replications, 1000L)
  expect_true(ctl$rho_max < 1)
  expect_false(is.na(ctl$seed))
})

test_that("brw_control rejects impossible settings", {
  expect_error(brw_control(replications = 1), "between 2")
  expect_error(brw_control(step = 0), "between")
  expect_error(brw_control(step = 2), "between")
  expect_error(brw_control(rho_max = 1.5), "between")
  expect_error(brw_control(max_iter = 0), "between 1")
  expect_error(brw_control(replications = c(10, 20)), "single finite number")
  expect_error(brw_control(replications = NA), "single finite number")
  expect_error(brw_control(seed = c(1, 2)), "single number or NA")
  expect_error(brw_control(method = "magic"), "should be one of")
})


# Bootstrap resampling ----------------------------------------------------

test_that("iid draws have the right shape and stay in range", {
  set.seed(1)
  idx <- brw_draw_index(n_v = 50, n_draw = 30, m = 7, block_length = 1L)
  expect_equal(dim(idx), c(30L, 7L))
  expect_true(all(idx >= 1 & idx <= 50))
  expect_gt(length(unique(as.vector(idx))), 20)   # genuinely resampling
})

test_that("block draws are runs of consecutive residuals", {
  set.seed(2)
  idx <- brw_draw_index(n_v = 60, n_draw = 24, m = 5, block_length = 6L)
  expect_equal(dim(idx), c(24L, 5L))
  expect_true(all(idx >= 1 & idx <= 60))

  # Within each block of 6, indices must increase by exactly one.
  for (b in seq_len(4L)) {
    rows <- (b - 1L) * 6L + seq_len(6L)
    expect_true(all(diff(idx[rows, , drop = FALSE]) == 1))
  }
})

test_that("blocks longer than the residual sample are refused", {
  expect_error(brw_draw_index(10, 20, 3, block_length = 11L),
               "exceeds the number of VAR residuals")
})

test_that("a block length that does not divide the sample is truncated, not padded", {
  set.seed(3)
  idx <- brw_draw_index(n_v = 40, n_draw = 25, m = 2, block_length = 7L)
  expect_equal(nrow(idx), 25L)
  expect_true(all(is.finite(idx)))
})


# The bias itself ---------------------------------------------------------

test_that("simulating from a known Phi reproduces the downward OLS bias", {
  # This is the premise of the whole correction: the mean OLS estimate on a
  # short sample sits below the parameter that generated the data.
  set.seed(4)
  phi <- matrix(0.95, 1, 1)
  z <- sim_ar(phi, 150)
  f <- var_ols(z)

  idx <- brw_draw_index(nrow(f$resid), 149L, 400L, 1L)
  z0 <- z[sample.int(150, 400, replace = TRUE), , drop = FALSE]
  mean_phi <- brw_ols_mean(phi, z0, f$resid, idx, 150L)

  expect_lt(mean_phi[1, 1], 0.95)
  expect_gt(mean_phi[1, 1], 0.88)
})

test_that("the bias shrinks as the sample lengthens", {
  set.seed(5)
  phi <- matrix(0.95, 1, 1)
  bias_at <- function(tt) {
    z <- sim_ar(phi, tt)
    f <- var_ols(z)
    idx <- brw_draw_index(nrow(f$resid), tt - 1L, 300L, 1L)
    z0 <- z[sample.int(tt, 300, replace = TRUE), , drop = FALSE]
    0.95 - brw_ols_mean(phi, z0, f$resid, idx, tt)[1, 1]
  }
  expect_gt(bias_at(120L), bias_at(1200L))
})

test_that("the correction removes most of the bias on simulated data", {
  # The substantive test. Repeated samples from a known persistent VAR: the
  # average OLS estimate is well below the truth, the average corrected
  # estimate is much closer.
  set.seed(6)
  phi <- matrix(0.95, 1, 1)
  ctl <- brw_control(replications = 200, seed = NA, max_iter = 30)

  est <- vapply(seq_len(40L), function(i) {
    z <- sim_ar(phi, 150L)
    f <- var_ols(z)
    bc <- brw_bias_correct(z, f$mu, f$phi, f$resid, ctl)
    c(ols = f$phi[1, 1], brw = bc$phi[1, 1])
  }, numeric(2))

  bias_ols <- mean(est["ols", ]) - 0.95
  bias_brw <- mean(est["brw", ]) - 0.95

  expect_lt(bias_ols, -0.01)                      # OLS is biased down
  expect_lt(abs(bias_brw), abs(bias_ols) / 2)     # and the correction fixes it
})

test_that("the one-step correction under-corrects relative to indirect inference", {
  set.seed(7)
  phi <- matrix(0.97, 1, 1)
  z <- sim_ar(phi, 200L)
  f <- var_ols(z)

  one <- brw_bias_correct(z, f$mu, f$phi, f$resid,
                          brw_control(method = "bootstrap", replications = 300))
  ind <- brw_bias_correct(z, f$mu, f$phi, f$resid,
                          brw_control(method = "indirect", replications = 300))

  expect_gt(one$phi[1, 1], f$phi[1, 1])           # both correct upwards
  expect_gt(ind$phi[1, 1], one$phi[1, 1])         # indirect goes further
  expect_equal(one$iterations, 1L)
})

test_that("indirect inference converges to its fixed point", {
  set.seed(8)
  z <- sim_ar(matrix(0.9, 1, 1), 200L)
  f <- var_ols(z)
  bc <- brw_bias_correct(z, f$mu, f$phi, f$resid,
                         brw_control(replications = 200))
  expect_true(bc$converged)
  expect_lt(bc$gap, brw_control()$tol)
  expect_lt(bc$iterations, 50L)
})

test_that("failure to converge is reported rather than hidden", {
  set.seed(9)
  z <- sim_ar(matrix(0.97, 1, 1), 200L)
  f <- var_ols(z)
  expect_warning(
    acm_brw_correct(list(mu = f$mu, phi = f$phi, var_residuals = f$resid,
                         lambda0 = 0, lambda1 = matrix(0, 1, 1)),
                    z, brw_control(replications = 100, max_iter = 1L,
                                   step = 0.1)),
    "did not converge"
  )
})


# Stationarity adjustment -------------------------------------------------

test_that("a stationary correction is left alone", {
  phi_hat <- matrix(0.9, 1, 1)
  phi_tilde <- matrix(0.95, 1, 1)
  out <- brw_shrink_stationary(phi_hat, phi_tilde, rho_max = 0.99)
  expect_equal(out$phi, phi_tilde)
  expect_equal(out$delta, 1)
})

test_that("an explosive correction is shrunk back towards OLS", {
  phi_hat <- matrix(0.9, 1, 1)
  phi_tilde <- matrix(1.1, 1, 1)
  out <- brw_shrink_stationary(phi_hat, phi_tilde, rho_max = 0.99)

  expect_lte(out$rho, 0.99)
  expect_gt(out$delta, 0)
  expect_lt(out$delta, 1)
  # It lies on the segment between the two estimates.
  expect_equal(drop(out$phi), 0.9 + out$delta * 0.2)
})

test_that("shrinkage cannot rescue an already-explosive OLS estimate", {
  phi_hat <- matrix(1.02, 1, 1)
  phi_tilde <- matrix(1.2, 1, 1)
  out <- brw_shrink_stationary(phi_hat, phi_tilde, rho_max = 0.9999)
  expect_equal(out$delta, 0)
  expect_equal(out$phi, phi_hat)
})

test_that("a binding rho_max caps the corrected persistence", {
  set.seed(10)
  z <- sim_ar(matrix(0.97, 1, 1), 200L)
  f <- var_ols(z)

  # Just above what OLS found, so the cap binds on the correction but leaves
  # the uncorrected estimate admissible.
  cap <- spectral_radius(f$phi) + 0.005
  bc <- brw_bias_correct(z, f$mu, f$phi, f$resid,
                         brw_control(replications = 200, rho_max = cap))

  expect_lte(spectral_radius(bc$phi), cap + 1e-12)
  expect_gt(bc$rho_unshrunk, cap)
  expect_gt(bc$delta, 0)
  expect_lt(bc$delta, 1)
})

test_that("a rho_max below the OLS estimate abandons the correction loudly", {
  set.seed(11)
  z <- sim_ar(matrix(0.97, 1, 1), 200L)
  f <- var_ols(z)
  cap <- spectral_radius(f$phi) - 0.01

  pars <- list(mu = f$mu, phi = f$phi, var_residuals = f$resid,
               lambda0 = 0, lambda1 = matrix(0, 1, 1))
  expect_warning(
    out <- acm_brw_correct(pars, z, brw_control(replications = 100,
                                                rho_max = cap)),
    "already at or above the persistence ceiling"
  )
  # Nothing was changed: the fit is the uncorrected one.
  expect_equal(out$pars$phi, f$phi)
  expect_equal(out$diagnostics$shrinkage, 0)
})


# Integration with atsm() -------------------------------------------------

test_that("the bias correction leaves the yield curve fit untouched", {
  # The defining property: the cross-section says nothing about how fast the
  # factors mean-revert, so correcting the P-dynamics must not move a single
  # fitted yield.
  f <- brw_fits()
  expect_equal(f$brw$fitted, f$ols$fitted, tolerance = 1e-12)
  expect_equal(f$brw$pars$phi - f$brw$pars$lambda1,
               f$ols$pars$phi - f$ols$pars$lambda1, tolerance = 1e-12)
  expect_equal(f$brw$spectral_radius[["risk_adjusted"]],
               f$ols$spectral_radius[["risk_adjusted"]])
})

test_that("the bias correction raises the persistence of the factor VAR", {
  f <- brw_fits()
  expect_gt(f$brw$spectral_radius[["real_world"]],
            f$ols$spectral_radius[["real_world"]])
  expect_lte(f$brw$spectral_radius[["real_world"]], brw_control()$rho_max)
})

test_that("more persistent dynamics move variation out of the term premium", {
  # BRW's substantive claim: OLS understates how long the market expects rate
  # levels to last, so the expectations component is too flat and the term
  # premium too variable.
  f <- brw_fits()
  j <- match(120L, f$ols$maturities)

  expect_lt(sd(f$brw$term_premium[, j]), sd(f$ols$term_premium[, j]))
  expect_gt(sd(f$brw$risk_neutral[, j]), sd(f$ols$risk_neutral[, j]))
})

test_that("the decomposition stays exact under bias correction", {
  f <- brw_fits()
  expect_equal(f$brw$fitted, f$brw$risk_neutral + f$brw$term_premium,
               tolerance = 1e-14)
})

test_that("the sigma2 cancellation survives bias correction", {
  # sigma2 enters both recursions identically and must still drop out.
  f <- brw_fits()
  pars <- f$brw$pars
  pars0 <- pars; pars0$sigma2 <- 0

  tp <- function(p) {
    n <- 120L
    acm_yields(acm_recursion(n, p, FALSE), f$brw$factors, n) -
      acm_yields(acm_recursion(n, p, TRUE), f$brw$factors, n)
  }
  expect_equal(tp(pars), tp(pars0), tolerance = 1e-14)
})

test_that("a fitted model carries its bias-correction diagnostics", {
  f <- brw_fits()
  expect_null(f$ols$brw)
  expect_equal(f$brw$p_dynamics, "brw")

  d <- f$brw$brw
  expect_equal(d$method, "indirect")
  expect_true(d$converged)
  expect_named(d$spectral_radius, c("ols", "corrected", "unshrunk"))
  expect_equal(d$spectral_radius[["ols"]],
               f$ols$spectral_radius[["real_world"]])
  expect_equal(dim(d$phi_ols), c(5L, 5L))
  expect_equal(d$phi_ols, f$ols$pars$phi)
})

test_that("print reports the change in persistence", {
  f <- brw_fits()
  out <- paste(capture.output(print(f$brw)), collapse = "\n")
  expect_match(out, "P-dynamics: brw")
  expect_match(out, "bias corr")
  expect_match(out, "half-life")
  expect_no_match(paste(capture.output(print(f$ols)), collapse = "\n"),
                  "bias corr")
})

test_that("predict carries the corrected dynamics through", {
  f <- brw_fits()
  p <- predict(f$brw, brw_panel())
  expect_equal(p$term_premium, as.vector(f$brw$term_premium),
               tolerance = 1e-12)
  # And differs from what the uncorrected model would have said.
  p_ols <- predict(f$ols, brw_panel())
  expect_false(isTRUE(all.equal(p$term_premium, p_ols$term_premium)))
})

test_that("atsm rejects an unknown p_dynamics", {
  expect_error(atsm(brw_panel(), n_factors = 3, p_dynamics = "kalman"),
               "should be one of")
})


# Reproducibility ---------------------------------------------------------

test_that("the same seed gives the same answer", {
  panel <- brw_panel()
  ctl <- brw_control(replications = 50)
  a <- suppressWarnings(atsm(panel, n_factors = 3, brw_control = ctl,
                             p_dynamics = "brw"))
  b <- suppressWarnings(atsm(panel, n_factors = 3, brw_control = ctl,
                             p_dynamics = "brw"))
  expect_equal(a$pars$phi, b$pars$phi)

  c_ <- suppressWarnings(atsm(panel, n_factors = 3, p_dynamics = "brw",
                              brw_control = brw_control(replications = 50,
                                                        seed = 999)))
  expect_false(isTRUE(all.equal(a$pars$phi, c_$pars$phi)))
  # Different seed, same neighbourhood: this is Monte Carlo noise, not a
  # different estimator.
  expect_lt(max(abs(a$pars$phi - c_$pars$phi)), 0.05)
})

test_that("fitting does not disturb the caller's random stream", {
  panel <- brw_panel()
  set.seed(123)
  before <- .Random.seed
  invisible(suppressWarnings(atsm(panel, n_factors = 3, p_dynamics = "brw",
                                  brw_control = brw_control(replications = 50))))
  expect_identical(before, .Random.seed)

  # And seed = NA genuinely uses the ambient stream.
  set.seed(321)
  d1 <- suppressWarnings(atsm(panel, n_factors = 3, p_dynamics = "brw",
                              brw_control = brw_control(replications = 50,
                                                        seed = NA)))
  d2 <- suppressWarnings(atsm(panel, n_factors = 3, p_dynamics = "brw",
                              brw_control = brw_control(replications = 50,
                                                        seed = NA)))
  expect_false(isTRUE(all.equal(d1$pars$phi, d2$pars$phi)))
})


# Guards ------------------------------------------------------------------

test_that("an absurd replication count is refused before it allocates", {
  panel <- brw_panel()
  expect_error(
    atsm(panel, n_factors = 5, p_dynamics = "brw",
         brw_control = brw_control(replications = 1e6)),
    "Reduce `replications`"
  )
})
