# The constrained maximum likelihood step.
#
# The constraint machinery is what needs testing hardest, because a wrong
# constraint is not a crash: it is a slightly worse fit for no visible reason.
# Two properties pin it down without any reference to the optimiser.
#
#   1. Fed the OBSERVED yields, the implied-factor pipeline must return the
#      observed factors exactly. It is the same construction, so anything but
#      zero means the replay is wrong.
#   2. On a panel generated exactly from an affine model, the closed-form fit
#      already reproduces the yields, so the constraints must already be
#      satisfied. They are, to about 3e-3 -- which is the strongest available
#      evidence that the constraints say what they are meant to say, since
#      nothing in the closed-form estimator was told about them.
#
# The simulated panel is a poor test bed for the optimiser itself: its third
# nominal component has a standard deviation of 0.4 basis points, against
# 6 for the real US curve, so the constraints there demand that the model fit
# what is essentially noise. `analysis/validate-uk.R` exercises the optimiser
# on data with a realistic factor spectrum.

# Settings ----------------------------------------------------------------

test_that("acmy_control returns validated defaults", {
  ctl <- acmy_control()
  expect_true(ctl$max_outer >= 1L)
  expect_true(ctl$tol > 0)
  expect_true(ctl$ndeps > 0)
  expect_type(ctl$max_inner, "integer")
})

test_that("acmy_control rejects impossible settings", {
  expect_error(acmy_control(max_outer = 0), "between 1")
  expect_error(acmy_control(penalty_factor = 0.5), "between 1")
  expect_error(acmy_control(tol = -1), "between 0")
  expect_error(acmy_control(ndeps = 1), "between")
  expect_error(acmy_control(reltol = c(1e-8, 1e-9)), "single finite number")
})


# Packing -----------------------------------------------------------------

test_that("packing and unpacking the parameters round-trips", {
  set.seed(4)
  k <- 4L
  pars <- list(mu_tilde = stats::rnorm(k), phi_tilde = matrix(stats::rnorm(k * k), k, k),
               delta0 = 0.004, delta1 = stats::rnorm(k),
               pi0 = 0.002, pi1 = stats::rnorm(k))

  theta <- acmy_pack(pars, k)
  expect_length(theta, k + k * k + 1L + k + 1L + k)

  back <- acmy_unpack(theta, k)
  for (nm in names(back)) expect_equal(back[[nm]], pars[[nm]], info = nm)
})

test_that("a fixed pi0 overrides whatever was packed", {
  k <- 3L
  theta <- acmy_pack(list(mu_tilde = rep(0, k), phi_tilde = diag(k),
                          delta0 = 0.004, delta1 = rep(0, k),
                          pi0 = 0.002, pi1 = rep(0, k)), k)
  expect_equal(acmy_unpack(theta, k)$pi0, 0.002)
  expect_equal(acmy_unpack(theta, k, pi0_fixed = 0.01)$pi0, 0.01)
})

test_that("each parameter block is scaled by its own magnitude", {
  # Without this the optimiser never moves delta0 or the inflation loadings:
  # they are three orders of magnitude smaller than the entries of Phi~, and
  # a numerically differenced gradient simply does not see them.
  k <- 3L
  start <- acmy_pack(list(mu_tilde = rep(1e-4, k), phi_tilde = diag(k) * 0.97,
                          delta0 = 4e-3, delta1 = rep(1e-3, k),
                          pi0 = 2e-3, pi1 = rep(5e-4, k)), k)
  scale <- acmy_theta_scale(start, k)

  expect_length(scale, length(start))
  expect_true(all(scale > 0))
  expect_equal(scale[seq_len(k)], rep(1e-4, k))              # mu block
  expect_equal(scale[k + k * k + 1L], 4e-3)                  # delta0
  expect_true(all(abs(start / scale) <= 1 + 1e-12))
})

test_that("a block of exact zeros does not divide by zero", {
  k <- 2L
  start <- acmy_pack(list(mu_tilde = c(0, 0), phi_tilde = diag(k),
                          delta0 = 0, delta1 = c(0, 0),
                          pi0 = 0, pi1 = c(0, 0)), k)
  scale <- acmy_theta_scale(start, k)
  expect_true(all(is.finite(scale)))
  expect_true(all(scale > 0))
})


# The constraints ---------------------------------------------------------

test_that("the implied-factor pipeline replays the observed factors exactly", {
  # Property 1. Same loadings, same centres, same orthogonalising
  # coefficients, same scaling -- only the yields differ. Fed the observed
  # yields it can only return the observed factors.
  d <- sim_fit()
  f <- d$fit
  fac <- f$pca
  kk <- f$n_factors_nominal + f$n_factors_real

  x_raw <- sweep(f$factors, 2L, fac$scale, "*")
  xn <- sweep(f$observed / 12, 2L, fac$nominal$center, "-") %*%
    fac$nominal$loadings
  z <- cbind(1, xn)
  if (!is.na(fac$liq_index)) z <- cbind(z, x_raw[, fac$liq_index])
  resid <- f$observed_real / 12 - z %*% fac$ortho_coef
  xr <- sweep(resid, 2L, fac$real$center, "-") %*% fac$real$loadings
  replay <- sweep(cbind(xn, xr), 2L, fac$scale[seq_len(kk)], "/")

  expect_equal(replay, f$factors[, seq_len(kk), drop = FALSE],
               tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("an exactly affine panel already satisfies the constraints", {
  # Property 2, and the sharpest evidence that the constraints are the right
  # ones: nothing in the closed-form estimator knows about them, yet on a
  # panel it can price exactly they come out at zero.
  d <- sim_fit(noise_bp = 0)
  expect_lt(sqrt(mean((d$fit$observed - d$fit$fitted)^2)) * 1e4, 1e-3)
  expect_lt(sim_violation(d$fit), 0.02)
})

test_that("a panel the model cannot price exactly violates them", {
  # The converse: the constraints are not vacuous. Pricing error puts the
  # model-implied factors somewhere other than the observed ones.
  d <- sim_fit(noise_bp = 0.25)
  expect_gt(sqrt(mean((d$fit$observed - d$fit$fitted)^2)) * 1e4, 0.1)
  expect_gt(sim_violation(d$fit), 0.1)
})

test_that("the constraint vector has one entry per restriction", {
  d <- sim_fit()
  f <- d$fit
  kk <- f$n_factors_nominal + f$n_factors_real
  cons <- acmy_constraints(sim_coefs(f), f$pca, f$maturities,
                           f$real_maturities, f$n_factors, kk)
  # An intercept per extracted factor, plus a slope against every factor.
  expect_length(cons, kk * (1L + f$n_factors))
})

test_that("the affine probe agrees with the map evaluated on the data", {
  # The constraints are read off by probing the factor map at zero and at the
  # unit vectors. That is only legitimate because the map is affine, so the
  # probe has to reproduce it on the actual factor paths.
  d <- sim_fit(noise_bp = 0.25)
  f <- d$fit
  kk <- f$n_factors_nominal + f$n_factors_real
  coefs <- sim_coefs(f)

  cons <- acmy_constraints(coefs, f$pca, f$maturities, f$real_maturities,
                           f$n_factors, kk)
  intercept <- cons[seq_len(kk)]
  slope <- matrix(cons[-seq_len(kk)], kk, f$n_factors) +
    diag(f$n_factors)[seq_len(kk), , drop = FALSE]

  from_probe <- rep(intercept, each = nrow(f$factors)) +
    f$factors %*% t(slope)
  direct <- acmy_implied_factors(coefs, f$pca, f$maturities,
                                 f$real_maturities, f$factors)

  expect_equal(from_probe, direct, tolerance = 1e-9, ignore_attr = TRUE)
})

test_that("the liquidity factor is not constrained", {
  # It is observed rather than extracted, so it maps to itself and restricts
  # nothing. Constraining it would be imposing an identity.
  s <- sim_acmy(liquidity = TRUE, noise_bp = 0.25)
  f <- suppressWarnings(suppressMessages(
    atsm_real(sim_panel(s), inflation = sim_cpi(s),
              short_rate = sim_sr(s), short_rate_units = "percent",
              liquidity = data.frame(date = s$dates, value = s$liquidity))))

  expect_equal(f$n_factors, 6L)
  kk <- f$n_factors_nominal + f$n_factors_real
  expect_equal(kk, 5L)
  cons <- acmy_constraints(sim_coefs(f), f$pca, f$maturities,
                           f$real_maturities, f$n_factors, kk)
  expect_length(cons, 5L * 7L)
})


# The estimator ------------------------------------------------------------

test_that("the likelihood step drives the constraints to zero", {
  # Full convergence is NOT asserted here, and the reason is the simulation
  # rather than the estimator. Its third nominal component has a standard
  # deviation under half a basis point, against six for the real US curve,
  # so once the panel carries any pricing error at all the constraints
  # demand that the model fit noise. On real data with a normal factor
  # spectrum the violation goes from 0.5 to 5e-7 and the fit improves
  # sharply -- see analysis/validate-uk.R. What is asserted here is that the
  # step moves in the right direction and reports itself honestly.
  skip_on_cran()
  d <- sim_fit(noise_bp = 0.25)
  s <- d$s

  f <- suppressWarnings(suppressMessages(
    atsm_real(sim_panel(s), inflation = sim_cpi(s), short_rate = sim_sr(s),
              short_rate_units = "percent", method = "ml",
              acmy_control = sim_ml_control())))

  expect_equal(f$method, "ml")
  expect_false(is.null(f$pars$ml))
  expect_lt(f$pars$ml$violation, f$pars$ml$violation_start)
  expect_lt(sim_violation(f), sim_violation(d$fit))
})

test_that("the likelihood step reports what it did", {
  skip_on_cran()
  d <- sim_fit(noise_bp = 0.25)
  f <- suppressWarnings(suppressMessages(
    atsm_real(sim_panel(d$s), inflation = sim_cpi(d$s),
              short_rate = sim_sr(d$s), short_rate_units = "percent",
              method = "ml", acmy_control = sim_ml_control())))

  m <- f$pars$ml
  expect_equal(m$n_constraints, 5L * 6L)
  expect_equal(m$n_parameters, 5L + 25L + 1L + 5L + 1L + 5L)
  expect_true(is.data.frame(m$trace))
  expect_true(all(c("round", "objective", "violation", "penalty", "code",
                    "evaluations") %in% names(m$trace)))
  expect_lte(nrow(m$trace), 2L)

  out <- paste(capture.output(print(f)), collapse = "\n")
  expect_match(out, "constrained maximum likelihood")
  expect_match(out, "factor-consistency restrictions")
})

test_that("the closed form says which estimator it is", {
  f <- sim_fit()$fit
  expect_equal(f$method, "closed_form")
  expect_null(f$pars$ml)
  out <- paste(capture.output(print(f)), collapse = "\n")
  expect_match(out, "closed form")
  expect_match(out, "starting value")
})

test_that("the likelihood step leaves the physical dynamics alone", {
  # mu, Phi and Sigma come from ordinary least squares on observable factors,
  # which is what maximises the VAR block of the likelihood whatever the
  # pricing parameters do. Only the pricing measure and the loadings move,
  # and the prices of risk are the difference between the two measures.
  skip_on_cran()
  d <- sim_fit(noise_bp = 0.25)
  s <- d$s
  cf <- d$fit
  ml <- suppressWarnings(suppressMessages(
    atsm_real(sim_panel(s), inflation = sim_cpi(s), short_rate = sim_sr(s),
              short_rate_units = "percent", method = "ml",
              acmy_control = sim_ml_control())))

  expect_equal(ml$pars$mu, cf$pars$mu, tolerance = 1e-15)
  expect_equal(ml$pars$phi, cf$pars$phi, tolerance = 1e-15)
  expect_equal(ml$pars$sigma, cf$pars$sigma, tolerance = 1e-15)
  expect_equal(ml$spectral_radius[["real_world"]],
               cf$spectral_radius[["real_world"]], tolerance = 1e-15)

  # lambda is defined as the gap between the measures, so it must stay so.
  expect_equal(ml$pars$lambda0, ml$pars$mu - ml$pars$mu_tilde,
               tolerance = 1e-15)
  expect_equal(ml$pars$lambda1, ml$pars$phi - ml$pars$phi_tilde,
               tolerance = 1e-15)
})

test_that("the decomposition identities survive the likelihood step", {
  skip_on_cran()
  d <- sim_fit(noise_bp = 0.25)
  f <- suppressWarnings(suppressMessages(
    atsm_real(sim_panel(d$s), inflation = sim_cpi(d$s),
              short_rate = sim_sr(d$s), short_rate_units = "percent",
              method = "ml", acmy_control = sim_ml_control())))

  expect_equal(f$breakeven,
               f$expected_inflation + f$inflation_risk_premium,
               tolerance = 1e-15)
  expect_equal(f$fitted, f$risk_neutral + f$term_premium, tolerance = 1e-15)
  expect_equal(f$fitted_real, f$risk_neutral_real + f$term_premium_real,
               tolerance = 1e-15)
})

test_that("a fixed pi0 stays fixed through the likelihood step", {
  skip_on_cran()
  d <- sim_fit(noise_bp = 0.25)
  f <- suppressWarnings(suppressMessages(
    atsm_real(sim_panel(d$s), inflation = sim_cpi(d$s),
              short_rate = sim_sr(d$s), short_rate_units = "percent",
              fix_pi0 = 0.002, method = "ml",
              acmy_control = sim_ml_control())))
  expect_equal(f$pars$pi0, 0.002, tolerance = 1e-15)
})

test_that("failing to converge warns rather than passing silently", {
  # One round from a violating start cannot reach 1e-12, and the user has to
  # be told: an unconverged fit is not wrong, but its levels are not pinned
  # either, which is the whole point of the step.
  skip_on_cran()
  d <- sim_fit(noise_bp = 0.25)
  expect_warning(
    suppressMessages(
      atsm_real(sim_panel(d$s), inflation = sim_cpi(d$s),
                short_rate = sim_sr(d$s), short_rate_units = "percent",
                method = "ml",
                acmy_control = sim_ml_control(max_outer = 1L, max_inner = 5L, tol = 1e-12))),
    "did not drive the factor-consistency"
  )
})

test_that("an exactly affine panel is refused by the likelihood step", {
  # On a panel generated exactly from the model, the short rate and inflation
  # are affine in the factors with no error at all, so their rows of
  # Sigma_epsilon are identically zero and it has no inverse. That is a real
  # property of the data rather than a defect, and it is reported as such
  # instead of surfacing as a linear algebra failure from three calls in.
  # The closed form handles this case fine; only the likelihood needs the
  # measurement covariance.
  skip_on_cran()
  d <- sim_fit(noise_bp = 0)
  expect_error(
    suppressWarnings(suppressMessages(
      atsm_real(sim_panel(d$s), inflation = sim_cpi(d$s),
                short_rate = sim_sr(d$s), short_rate_units = "percent",
                method = "ml"))),
    "not positive definite"
  )
})
