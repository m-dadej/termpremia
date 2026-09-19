# The joint real-nominal model (ACMY 2016).
#
# The evidence here is of two kinds and it is worth keeping them apart.
#
# The algebraic tests are exact. The pricing recursions imply their own
# excess-return expressions as identities, with no expectation taken and no
# estimation involved, so any discrepancy beyond machine precision is a bug.
# One such bug was found this way: `phi_adj %*% bt` in place of
# `t(phi_adj) %*% bt` violated Eq. 18 at 1e-5, which is small enough to pass
# for rounding error on a near-diagonal Phi and is not.
#
# The recovery tests simulate a panel FROM the model at known parameters and
# check the estimator gets them back. That is strong evidence about the
# arithmetic and weak evidence about real data, exactly as in
# test-survey-dynamics.R. `analysis/replicate-acmy.R` does the real-data half.


# Simulation helper -------------------------------------------------------

# A stationary five-factor truth with a stable pricing measure. Correlated
# innovations are drawn through a Cholesky factor rather than MASS::mvrnorm,
# so the tests need no package that is not already a dependency.
# `noise_bp` adds i.i.d. pricing error to the simulated curves. It defaults
# to zero, which makes recovery exact and is the sharpest test available --
# but only because the GLS steps whiten by a Cholesky factor of the return
# error covariance instead of inverting it. A panel generated exactly from
# the model has pricing errors that are pure floating-point noise, and
# solving the normal equations against that covariance fails outright.
#
# When it is non-zero, keep it small. The noise is i.i.d. ACROSS maturities
# and a log price multiplies it by n, so two basis points on a ten-year yield
# is 240 basis points of noise in its log price, swamping the one-month
# return the estimator is built on. A fitted curve's own error is smooth in
# maturity and does no such thing.
sim_acmy <- function(tt = 260L, seed = 7L, liquidity = FALSE,
                     noise_bp = 0) {
  set.seed(seed)

  phi <- diag(c(0.990, 0.960, 0.930, 0.970, 0.940))
  phi[1L, 2L] <- 0.010
  phi[2L, 3L] <- -0.008
  phi[4L, 5L] <- 0.006
  mu <- c(4e-5, -1e-5, 5e-6, 2e-6, -1e-6)

  sd_f <- c(3e-3, 2e-3, 1e-3, 8e-4, 5e-4)
  sigma <- diag(sd_f^2)
  sigma[1L, 2L] <- sigma[2L, 1L] <- 0.2 * sd_f[1L] * sd_f[2L]

  lambda1 <- diag(c(0.02, 0.015, 0.01, 0.012, 0.008))
  delta1 <- c(1, -0.6, 0.3, 0.05, -0.02) * 1e-2

  # Inflation loadings large enough that inflation actually varies: at a
  # tenth of this, simulated annual inflation has a standard deviation of
  # under two basis points, and pi1 is then so weakly identified that the
  # recovery test measures rounding rather than arithmetic.
  pi1 <- c(0.15, -0.1, 0.05, 0.4, -0.25) * 1e-1

  if (liquidity) {
    # A sixth factor standing in for an observable liquidity index. Three
    # properties are needed for it to be a fair test rather than a
    # decoration:
    #
    #  - It is O(1) and mostly positive, like a standardised index shifted to
    #    be non-negative, which is what `tips_liquidity_factor()` returns.
    #    This is also the scale mismatch against O(1e-5) yield components
    #    that broke the GLS steps before the factors were normalised.
    #  - `delta1` puts no weight on it, so it leaves nominal yields alone --
    #    the paper's own finding (Section 3.3).
    #  - `pi1` does put weight on it, which is the only channel through which
    #    an extra state variable can reach real yields in this model, so
    #    TIPS returns load on it and its price of risk is identified. A
    #    liquidity series the bonds do not price leaves the estimator with a
    #    rank-deficient B, which it now refuses by name.
    phi <- rbind(cbind(phi, 0), c(rep(0, 5L), 0.950))
    mu <- c(mu, 0.05)
    sigma <- rbind(cbind(sigma, 0), c(rep(0, 5L), 0.12^2))
    lambda1 <- rbind(cbind(lambda1, 0), c(rep(0, 5L), 0.010))
    delta1 <- c(delta1, 0)
    pi1 <- c(pi1, 3e-4)
  }

  k <- length(mu)
  lambda0 <- rep(0, k)
  delta0 <- 0.0035
  pi0 <- 0.0018

  mu_t <- mu - lambda0
  phi_t <- phi - lambda1

  chol_s <- chol(sigma)
  x <- matrix(0, tt, k)
  if (liquidity) x[1L, 6L] <- mu[6L] / (1 - phi[6L, 6L])
  for (t in 2L:tt) {
    x[t, ] <- mu + drop(phi %*% x[t - 1L, ]) +
      drop(crossprod(chol_s, stats::rnorm(k)))
  }

  mat_n <- 1:120
  mat_r <- 23:120
  co <- acmy_coefficients(120L, mu_t, phi_t, sigma, delta0, delta1, pi0, pi1)

  firsts <- seq(as.Date("1999-01-01"), by = "month", length.out = tt + 1L)
  dates <- firsts[-1L] - 1L

  noise <- function(nr, nc) {
    if (noise_bp <= 0) return(matrix(0, nr, nc))
    matrix(stats::rnorm(nr * nc, sd = noise_bp / 1e4), nr, nc)
  }

  list(
    x = x, dates = dates, mat_n = mat_n, mat_r = mat_r, co = co,
    noise_bp = noise_bp,
    y_nom = acmy_yields(co$a, co$b, x, mat_n) * 12 +
      noise(tt, length(mat_n)),
    y_real = acmy_yields(co$a_real, co$b_real, x, mat_r) * 12 +
      noise(tt, length(mat_r)),
    inflation = pi0 + drop(x %*% pi1),
    short_rate = delta0 + drop(x %*% delta1),
    liquidity = if (liquidity) x[, 6L] else NULL,
    truth = list(mu = mu, phi = phi, sigma = sigma, mu_tilde = mu_t,
                 phi_tilde = phi_t, lambda0 = lambda0, lambda1 = lambda1,
                 delta0 = delta0, delta1 = delta1, pi0 = pi0, pi1 = pi1)
  )
}

sim_panel <- function(s, end = NULL) {
  long <- rbind(
    data.frame(date = rep(s$dates, length(s$mat_n)),
               maturity = rep(s$mat_n, each = length(s$dates)),
               yield = as.vector(s$y_nom) * 100, curve = "nominal",
               stringsAsFactors = FALSE),
    data.frame(date = rep(s$dates, length(s$mat_r)),
               maturity = rep(s$mat_r, each = length(s$dates)),
               yield = as.vector(s$y_real) * 100, curve = "tips",
               stringsAsFactors = FALSE)
  )
  if (!is.null(end)) long <- long[long$date <= end, , drop = FALSE]

  yield_panel(long, curve = "curve", units = "percent",
              maturity_unit = "months",
              instrument = c(nominal = "government", tips = "tips"),
              issuer = "SIM")
}

sim_cpi <- function(s) {
  data.frame(date = s$dates, value = 100 * exp(cumsum(s$inflation)),
             stringsAsFactors = FALSE)
}

sim_sr <- function(s) {
  data.frame(date = s$dates, value = s$short_rate * 12 * 100,
             stringsAsFactors = FALSE)
}

sim_fit <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      s <- sim_acmy()
      cache <<- list(
        s = s,
        fit = suppressWarnings(suppressMessages(
          atsm_real(sim_panel(s), inflation = sim_cpi(s),
                    short_rate = sim_sr(s), short_rate_units = "percent")))
      )
    }
    cache
  }
})


# Pricing recursions: exact algebra ---------------------------------------

test_that("the one-month nominal yield is the short rate", {
  set.seed(3)
  k <- 4L
  delta0 <- 0.004
  delta1 <- stats::rnorm(k) * 1e-3

  co <- acmy_coefficients(24L, stats::rnorm(k) * 1e-4, diag(k) * 0.97,
                          diag(k) * 1e-6, delta0, delta1,
                          0.002, stats::rnorm(k) * 5e-4)
  x <- matrix(stats::rnorm(20L * k), 20L, k) * 0.01

  expect_equal(-(co$a[1L] + drop(x %*% co$b[1L, ])),
               delta0 + drop(x %*% delta1), tolerance = 1e-15)
})

test_that("the one-month real yield is the nominal rate less expected inflation", {
  # A sign error anywhere in the pi terms shows up here first: at n = 1 the
  # real yield must be r_t - E^Q[pi_{t+1}] - (1/2) pi1' Sigma pi1, the last
  # term being the Jensen correction from pricing a log-normal payoff.
  set.seed(11)
  k <- 4L
  mu_adj <- rnorm(k) * 1e-4
  phi_adj <- diag(k) * 0.97 + matrix(rnorm(k * k), k, k) * 0.01
  sigma <- crossprod(matrix(rnorm(k * k), k, k)) * 1e-6
  delta0 <- 0.004
  delta1 <- rnorm(k) * 1e-3
  pi0 <- 0.002
  pi1 <- rnorm(k) * 5e-4

  co <- acmy_coefficients(120L, mu_adj, phi_adj, sigma, delta0, delta1,
                          pi0, pi1)
  x <- matrix(rnorm(50L * k), 50L, k) * 0.01

  y1r <- -(co$a_real[1L] + drop(x %*% co$b_real[1L, ]))
  r_t <- delta0 + drop(x %*% delta1)
  e_pi <- pi0 + drop(sweep(x %*% t(phi_adj), 2L, mu_adj, "+") %*% pi1)
  jensen <- 0.5 * drop(crossprod(pi1, sigma %*% pi1))

  expect_equal(y1r, r_t - e_pi - jensen, tolerance = 1e-15)
})

test_that("modelled excess returns are an exact identity, not an approximation", {
  # Eq. 15 and Eq. 18 follow from the recursions by substitution. They must
  # reproduce p(n-1, t+1) - p(n, t) - r_t to machine precision. An asymmetric
  # phi_adj is essential: with a symmetric one, transposing it is invisible.
  set.seed(5)
  k <- 4L
  mu_adj <- rnorm(k) * 1e-4
  phi_adj <- diag(k) * 0.96 + matrix(rnorm(k * k), k, k) * 0.03
  sigma <- crossprod(matrix(rnorm(k * k), k, k)) * 1e-6
  delta0 <- 0.004
  delta1 <- rnorm(k) * 1e-3
  pi0 <- 0.002
  pi1 <- rnorm(k) * 5e-4

  co <- acmy_coefficients(120L, mu_adj, phi_adj, sigma, delta0, delta1,
                          pi0, pi1)
  tt <- 40L
  x <- matrix(rnorm(tt * k), tt, k) * 0.01
  x_lag <- x[-tt, , drop = FALSE]
  x_led <- x[-1L, , drop = FALSE]
  mats <- c(24L, 36L, 60L, 120L)

  modelled <- acmy_real_returns(co, mats, x_lag, x_led, mu_adj, phi_adj,
                                sigma, pi0, pi1)
  direct <- vapply(mats, function(n) {
    (co$a_real[n - 1L] + drop(x_led %*% co$b_real[n - 1L, ])) -
      (co$a_real[n] + drop(x_lag %*% co$b_real[n, ])) -
      (delta0 + drop(x_lag %*% delta1))
  }, numeric(tt - 1L))

  expect_equal(modelled, direct, tolerance = 1e-14, ignore_attr = TRUE)
})

test_that("asymmetric dynamics are required for that test to have teeth", {
  # Guards the test above: if phi_adj were symmetric, Phi'B and PhiB agree and
  # the identity would hold even with the transpose wrong.
  set.seed(5)
  k <- 4L
  phi_adj <- diag(k) * 0.96 + matrix(rnorm(k * k), k, k) * 0.03
  expect_false(isTRUE(all.equal(phi_adj, t(phi_adj))))
})


# Factor construction -----------------------------------------------------

test_that("real factors are orthogonal to the nominal ones by construction", {
  s <- sim_acmy(tt = 120L)
  fac <- acmy_factors(s$y_nom / 12, s$y_real / 12, 3L, 2L)

  x_nom <- fac$x[, 1:3, drop = FALSE]
  x_real <- fac$x[, 4:5, drop = FALSE]
  expect_equal(crossprod(x_nom, x_real),
               matrix(0, 3L, 2L), tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(ncol(fac$x), 5L)
  expect_true(is.na(fac$liq_index))
})

test_that("a liquidity factor is appended last and projected out of the real block", {
  s <- sim_acmy(tt = 120L)
  liq <- abs(sin(seq_len(120L) / 7)) * 2
  fac <- acmy_factors(s$y_nom / 12, s$y_real / 12, 3L, 2L, liquidity = liq)

  expect_equal(ncol(fac$x), 6L)
  expect_equal(fac$liq_index, 6L)
  expect_equal(fac$x_unscaled[, 6L], liq, ignore_attr = TRUE)

  # Orthogonal to the liquidity factor as well, which is the point of
  # including it in the projection rather than only in the state.
  x_real <- fac$x[, 4:5, drop = FALSE]
  expect_equal(drop(crossprod(liq - mean(liq), x_real)),
               c(0, 0), tolerance = 1e-9, ignore_attr = TRUE)
})

test_that("factors are normalised to unit scale without being centred", {
  # Scale, because the GLS steps are otherwise insoluble once an O(1)
  # liquidity index sits beside O(1e-5) yield components. Not location,
  # because centring would destroy the non-negativity that makes the
  # liquidity component's sign interpretable.
  s <- sim_acmy(tt = 120L)
  liq <- abs(sin(seq_len(120L) / 7)) * 2 + 0.5
  fac <- acmy_factors(s$y_nom / 12, s$y_real / 12, 3L, 2L, liquidity = liq)

  expect_equal(apply(fac$x, 2L, stats::sd), rep(1, 6L), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_true(all(fac$x[, 6L] > 0))
  expect_equal(fac$x[, 6L] * fac$scale[[6L]], liq, ignore_attr = TRUE)
})

test_that("rescaling the state leaves the model unchanged", {
  # The invariance that makes normalising the factors safe: yields and the
  # whole decomposition are identical under an affine change of state.
  set.seed(19)
  k <- 4L
  mu <- stats::rnorm(k) * 1e-4
  phi <- diag(k) * 0.96 + matrix(stats::rnorm(k * k), k, k) * 0.02
  sigma <- crossprod(matrix(stats::rnorm(k * k), k, k)) * 1e-6
  delta0 <- 0.004
  delta1 <- stats::rnorm(k) * 1e-3
  pi0 <- 0.002
  pi1 <- stats::rnorm(k) * 5e-4
  x <- matrix(stats::rnorm(60L * k), 60L, k) * 0.01

  d <- c(1, 250, 0.004, 40)          # an arbitrary rescaling
  dd <- diag(d)
  dd_inv <- diag(1 / d)

  co <- acmy_coefficients(120L, mu, phi, sigma, delta0, delta1, pi0, pi1)
  co_s <- acmy_coefficients(
    120L, drop(dd_inv %*% mu), dd_inv %*% phi %*% dd,
    dd_inv %*% sigma %*% dd_inv, delta0, drop(dd %*% delta1),
    pi0, drop(dd %*% pi1)
  )
  x_s <- x %*% dd_inv

  mats <- c(24L, 60L, 120L)
  expect_equal(acmy_yields(co$a, co$b, x, mats),
               acmy_yields(co_s$a, co_s$b, x_s, mats), tolerance = 1e-14)
  expect_equal(acmy_yields(co$a_real, co$b_real, x, mats),
               acmy_yields(co_s$a_real, co_s$b_real, x_s, mats),
               tolerance = 1e-14)
})

test_that("a factor the bonds do not price is refused by name", {
  # A liquidity series that does not move bond returns leaves B column-rank
  # deficient and every GLS step insoluble. Before this check the failure was
  # LAPACK's "singular matrix 'a' in solve", three calls deep, naming
  # nothing. Tested on the loading matrix directly: routed through a full fit,
  # the simulated pricing error is enough to give a junk factor a spurious
  # loading and restore the rank, which hides the case being guarded.
  set.seed(31)
  b <- cbind(matrix(stats::rnorm(20L * 5L), 20L, 5L), 0)
  nms <- c(paste0("N", 1:3), paste0("R", 1:2), "liquidity")

  expect_error(check_return_loadings(b, 6L, nms), "rank 5")
  expect_error(check_return_loadings(b, 6L, nms),
               "weakest exposure is to .liquidity.")
  expect_error(check_return_loadings(b, 6L, nms), "mis-aligned or irrelevant")

  # Names the right culprit when it is not the liquidity factor.
  b2 <- cbind(matrix(stats::rnorm(20L * 4L), 20L, 4L), 0,
              stats::rnorm(20L))
  expect_error(check_return_loadings(b2, 6L, nms),
               "weakest exposure is to .R2.")
  expect_error(check_return_loadings(b2, 6L, nms), "n_factors_real")

  expect_true(check_return_loadings(matrix(stats::rnorm(20L * 6L), 20L, 6L),
                                    6L, nms))
})


# Recovery ----------------------------------------------------------------

test_that("the estimator prices a panel generated by the model", {
  f <- sim_fit()$fit

  # An exactly affine panel must be priced exactly: a thousandth of a basis
  # point, not a basis point.
  expect_lt(sqrt(mean((f$observed - f$fitted)^2)) * 1e4, 1e-3)
  expect_lt(sqrt(mean((f$observed_real - f$fitted_real)^2)) * 1e4, 1e-3)
})

test_that("the pricing measure is recovered from simulated data", {
  d <- sim_fit()
  truth <- max(Mod(eigen(d$s$truth$phi_tilde, only.values = TRUE)$values))
  expect_equal(d$fit$spectral_radius[["risk_adjusted"]], truth,
               tolerance = 1e-5)
})

test_that("the inflation process is recovered from real bond returns", {
  # pi0 and pi1 are identified only through TIPS returns (Eq. S31), so this
  # is the test that the criterion and the recursion agree with each other.
  d <- sim_fit()
  implied <- d$fit$pars$pi0 + drop(d$fit$factors %*% d$fit$pars$pi1)
  rmse_bp <- sqrt(mean((implied - d$fit$inflation)^2)) * 12 * 1e4
  # Judged against how much simulated inflation actually varies, not against
  # zero: recovering a series with 19bp of annual variation to within 3bp is
  # the claim.
  expect_lt(rmse_bp, 0.1)
  expect_gt(stats::sd(d$fit$inflation) * 12 * 1e4, 10)
  expect_true(d$fit$pars$inflation_fit$converged)
})

test_that("the factor VAR is estimated with the usual downward bias", {
  # Not a defect: OLS on 259 observations of a near-unit-root VAR
  # under-states persistence, which is exactly what p_dynamics = "brw"
  # exists to correct on the single-curve side. Recorded so that a future
  # change in this direction is not mistaken for an improvement.
  d <- sim_fit()
  truth <- max(Mod(eigen(d$s$truth$phi, only.values = TRUE)$values))
  expect_lt(d$fit$spectral_radius[["real_world"]], truth)
  expect_gt(d$fit$spectral_radius[["real_world"]], truth - 0.05)
})


# The decomposition identities --------------------------------------------

test_that("breakeven splits exactly into expected inflation and risk premium", {
  f <- sim_fit()$fit
  expect_equal(f$breakeven,
               f$expected_inflation + f$inflation_risk_premium,
               tolerance = 1e-15)
})

test_that("the inflation risk premium is the difference of the term premia", {
  # Follows from expected inflation being the physical-measure breakeven.
  # Exact, not approximate.
  f <- sim_fit()$fit
  jn <- match(f$shared_maturities, f$maturities)
  jr <- match(f$shared_maturities, f$real_maturities)

  expect_equal(
    f$inflation_risk_premium,
    f$term_premium[, jn, drop = FALSE] -
      f$term_premium_real[, jr, drop = FALSE],
    tolerance = 1e-15, ignore_attr = TRUE
  )
})

test_that("each curve splits into its own expectations and premium", {
  f <- sim_fit()$fit
  expect_equal(f$fitted, f$risk_neutral + f$term_premium, tolerance = 1e-15)
  expect_equal(f$fitted_real, f$risk_neutral_real + f$term_premium_real,
               tolerance = 1e-15)
})

test_that("expected inflation is the breakeven computed without prices of risk", {
  d <- sim_fit()
  f <- d$fit
  co <- acmy_coefficients(max(f$maturities), f$pars$mu, f$pars$phi,
                          f$pars$sigma, f$pars$delta0, f$pars$delta1,
                          f$pars$pi0, f$pars$pi1)
  shared <- f$shared_maturities
  yn <- acmy_yields(co$a, co$b, f$factors, shared) * 12
  yr <- acmy_yields(co$a_real, co$b_real, f$factors, shared) * 12

  expect_equal(f$expected_inflation, yn - yr, tolerance = 1e-15,
               ignore_attr = TRUE)
})


# Interface ---------------------------------------------------------------

test_that("atsm_real needs two curves", {
  panel <- yield_panel(gsw_monthly, units = "percent",
                       maturity_unit = "months", issuer = "US")
  expect_error(atsm_real(panel, inflation = data.frame(date = Sys.Date(),
                                                       value = 100)),
               "needs two curves")
})

test_that("the inflation-indexed curve is identified from the metadata", {
  s <- sim_acmy(tt = 120L)
  p <- sim_panel(s)
  expect_equal(resolve_real_nominal_curves(p, NULL, NULL),
               list(nominal = "nominal", real = "tips"))
})

test_that("an ambiguous or absent tips curve is reported, not guessed", {
  s <- sim_acmy(tt = 120L)
  p <- sim_panel(s)
  p$meta$instrument <- c("government", "government")
  expect_error(resolve_real_nominal_curves(p, NULL, NULL),
               "no curve has instrument")

  p$meta$instrument <- c("tips", "tips")
  expect_error(resolve_real_nominal_curves(p, NULL, NULL), "2 curves do")
})

test_that("inflation is required and must be a price level or a rate", {
  s <- sim_acmy(tt = 120L)
  p <- sim_panel(s)

  expect_error(atsm_real(p), "`inflation` is required")
  expect_error(atsm_real(p, inflation = 1:10), "data frame")
  expect_error(
    atsm_real(p, inflation = data.frame(date = s$dates, value = -1)),
    "non-positive"
  )
})

test_that("a gap in the usable sample is an error, not a silent shortening", {
  # The factor VAR reads consecutive rows as consecutive months, so a hole
  # would corrupt every persistence estimate with nothing reported.
  s <- sim_acmy(tt = 120L)
  cpi <- sim_cpi(s)
  cpi <- cpi[-50L, , drop = FALSE]

  expect_error(
    suppressMessages(atsm_real(sim_panel(s), inflation = cpi,
                               short_rate = sim_sr(s),
                               short_rate_units = "percent")),
    "internal gap"
  )
})

test_that("too short a usable overlap is refused", {
  s <- sim_acmy(tt = 120L)
  cpi <- sim_cpi(s)
  cpi <- cpi[1:30, , drop = FALSE]

  expect_error(
    suppressMessages(atsm_real(sim_panel(s), inflation = cpi)),
    "needs a usable sample"
  )
})

test_that("off-grid maturities are refused on either curve", {
  s <- sim_acmy(tt = 120L)
  p <- sim_panel(s)
  cpi <- sim_cpi(s)

  expect_error(atsm_real(p, inflation = cpi, maturities = c(1L, 500L)),
               "not on the 'nominal' grid")
  expect_error(atsm_real(p, inflation = cpi, real_maturities = c(24L, 500L)),
               "not on the 'tips' grid")
})


# Extractors and display --------------------------------------------------

test_that("each component is served on its own maturity grid", {
  f <- sim_fit()$fit

  expect_equal(sort(unique(term_premium(f)$maturity)), f$maturities)
  expect_equal(sort(unique(real_term_premium(f)$maturity)), f$real_maturities)
  expect_equal(sort(unique(breakeven(f)$maturity)), f$shared_maturities)
  expect_equal(sort(unique(expected_inflation(f)$maturity)),
               f$shared_maturities)
})

test_that("asking for a maturity a component does not have is an error", {
  f <- sim_fit()$fit
  # 12 months exists on the nominal curve but not on the real one, so it
  # cannot have a breakeven.
  expect_silent(term_premium(f, maturity = 12L))
  expect_error(breakeven(f, maturity = 12L), "not available")
  expect_error(real_term_premium(f, maturity = 12L), "not available")
})

test_that("the generic extractors dispatch on the joint fit", {
  f <- sim_fit()$fit
  expect_s3_class(f, "atsm_real_fit")
  expect_equal(nrow(term_premium(f, maturity = 120L)), length(f$dates))
  expect_equal(risk_neutral(f, maturity = 120L)$value,
               expected_short_rate(f, maturity = 120L)$value)
})

test_that("liquidity output is absent, and refused, when none was supplied", {
  f <- sim_fit()$fit
  expect_false(f$has_liquidity)
  expect_null(f$liquidity_real)
  expect_error(expected_inflation(f, liquidity_adjusted = TRUE),
               "needs a model fitted with")
  expect_error(liquidity_premium(f), "no 'liquidity_real' component")
})

test_that("a fit with a liquidity factor reports and adjusts for it", {
  s <- sim_acmy(liquidity = TRUE)
  liq <- data.frame(date = s$dates, value = s$liquidity)
  f <- suppressWarnings(suppressMessages(
    atsm_real(sim_panel(s), inflation = sim_cpi(s), short_rate = sim_sr(s),
              short_rate_units = "percent", liquidity = liq)))

  expect_true(f$has_liquidity)
  expect_equal(f$n_factors, 6L)
  expect_false(is.null(f$liquidity_real))

  # The adjustment moves the breakeven and the risk premium by the same
  # amount and leaves expected inflation alone, so the identity survives it.
  be <- breakeven(f, liquidity_adjusted = TRUE)$value
  ei <- expected_inflation(f, liquidity_adjusted = TRUE)$value
  irp <- inflation_risk_premium(f, liquidity_adjusted = TRUE)$value
  expect_equal(be, ei + irp, tolerance = 1e-14)
  expect_equal(ei, expected_inflation(f)$value, tolerance = 1e-15)
  expect_false(isTRUE(all.equal(be, breakeven(f)$value)))
})

test_that("print and summary describe the fit", {
  f <- sim_fit()$fit
  out <- paste(capture.output(print(f)), collapse = "\n")
  expect_match(out, "atsm_real_fit")
  expect_match(out, "ACMY joint real-nominal")
  expect_match(out, "breakeven inflation")
  expect_match(out, "expected inflation")
  expect_match(out, "without a liquidity factor")

  sm <- summary(f)
  expect_s3_class(sm, "summary.atsm_real_fit")
  expect_true(all(c("breakeven_bp", "exp_infl_bp", "irp_bp") %in%
                    names(sm$per_maturity)))
  expect_silent(invisible(capture.output(print(sm))))
})

test_that("coef, fitted and residuals work", {
  f <- sim_fit()$fit
  cf <- coef(f)
  expect_true(all(c("pi0", "pi1", "mu_tilde", "phi_tilde") %in% names(cf)))
  expect_equal(nrow(fitted(f)), length(f$dates) * length(f$maturities))
  expect_lt(max(abs(residuals(f)$value)) * 1e4, 0.01)
})


# Liquidity factor construction -------------------------------------------

test_that("tips_liquidity_factor standardises, averages and shifts positive", {
  d <- seq(as.Date("2000-01-31"), by = "month", length.out = 24L)
  a <- data.frame(date = d, value = seq_len(24L))
  b <- data.frame(date = d, value = rev(seq_len(24L)) * 100)

  out <- tips_liquidity_factor(a, b)
  expect_equal(nrow(out), 24L)
  expect_equal(min(out$value), 0)
  expect_true(all(out$value >= 0))

  # Equal weights on standardised inputs: two exactly opposing ramps cancel.
  expect_lt(stats::sd(out$value), 1e-12)
})

test_that("tips_liquidity_factor accepts one indicator and rejects bad input", {
  d <- seq(as.Date("2000-01-31"), by = "month", length.out = 12L)
  ind <- data.frame(date = d, value = c(rep(1, 6), rep(3, 6)))

  out <- tips_liquidity_factor(ind)
  expect_equal(min(out$value), 0)
  expect_gt(max(out$value), 0)

  expect_error(tips_liquidity_factor(), "at least one")
  expect_error(tips_liquidity_factor(data.frame(date = d, value = rep(1, 12))),
               "zero or undefined variance")
  expect_error(tips_liquidity_factor(1:5), "`dates` is required")
  expect_error(tips_liquidity_factor(data.frame(x = 1)), "`date` and `value`")
})


# Reading the published TIPS curve ----------------------------------------

# A fixture standing in for feds200805.csv. It reproduces the three traps the
# real file contains: a preamble of a different length from the nominal
# file's, Nelson-Siegel rows encoded with NA in BETA3/TAU2 rather than the
# nominal file's -999.99 sentinel, and holiday rows with no curve at all --
# one of which is the last calendar day of its month.
tips_fixture <- function(path) {
  rows <- data.frame(
    Date = c("2003-05-29", "2003-05-30",           # NS rows, month-end on 30th
             "2003-06-26", "2003-06-27", "2003-06-30",  # no curve on the 30th
             "2004-01-02", "2004-01-30"),         # Svensson rows
    BETA0 = c(2.1, 2.2, 2.3, 2.4, NA, 2.5, 2.6),
    BETA1 = c(-1.0, -1.1, -1.2, -1.3, NA, -1.4, -1.5),
    BETA2 = c(0.5, 0.6, 0.7, 0.8, NA, 0.9, 1.0),
    BETA3 = c(NA, NA, NA, NA, NA, 0.4, 0.5),
    TAU1 = c(2.0, 2.1, 2.2, 2.3, NA, 2.4, 2.5),
    TAU2 = c(NA, NA, NA, NA, NA, 9.0, 9.5),
    stringsAsFactors = FALSE
  )
  rows$TIPSY02 <- vapply(seq_len(nrow(rows)), function(i) {
    svensson_yield(2, rows$BETA0[i], rows$BETA1[i], rows$BETA2[i],
                   rows$BETA3[i], rows$TAU1[i], rows$TAU2[i])
  }, numeric(1L))
  rows$TIPSY10 <- vapply(seq_len(nrow(rows)), function(i) {
    svensson_yield(10, rows$BETA0[i], rows$BETA1[i], rows$BETA2[i],
                   rows$BETA3[i], rows$TAU1[i], rows$TAU2[i])
  }, numeric(1L))

  con <- file(path, open = "wt")
  on.exit(close(con))
  writeLines(c("Note: not an official statistical release.", "", "Series,x",
               "TIPS Yields", "Parameters,N/A,BETA0 to TAU2", ""), con)
  utils::write.csv(rows, con, row.names = FALSE, na = "NA")
  invisible(path)
}

test_that("the TIPS file's header is found rather than assumed", {
  path <- withr::local_tempfile(fileext = ".csv")
  tips_fixture(path)
  raw <- read_fed_curve_file(path)
  expect_equal(nrow(raw), 7L)
  expect_true(all(c("BETA0", "TAU2", "TIPSY02") %in% names(raw)))
  expect_s3_class(raw$Date, "Date")
})

test_that("Nelson-Siegel rows survive and no-curve rows are dropped first", {
  # The load-bearing assertion is 2003-06: the month's last row carries no
  # curve, so month-end selection has to fall back to the 27th rather than
  # dropping June. Doing this in the other order is the bug that cost twenty
  # months on the nominal curve.
  path <- withr::local_tempfile(fileext = ".csv")
  tips_fixture(path)

  out <- gsw_tips(maturities = 23:120, cache = path, check = TRUE)
  got <- sort(unique(out$date))

  expect_equal(format(got, "%Y-%m"), c("2003-05", "2003-06", "2004-01"))
  expect_equal(got[2L], as.Date("2003-06-27"))
  expect_equal(got[1L], as.Date("2003-05-30"))
  expect_false(anyNA(out$yield))
})

test_that("the evaluated curve is checked against the file's own tenors", {
  path <- withr::local_tempfile(fileext = ".csv")
  tips_fixture(path)

  # Corrupt the published column on every row -- including the month-end rows
  # that survive filtering -- and the check must refuse the data.
  raw <- readLines(path)
  hdr <- grep("^\"?Date\"?,", raw)
  cols <- gsub('"', "", strsplit(raw[hdr], ",")[[1L]])
  jj <- which(cols == "TIPSY10")
  expect_length(jj, 1L)

  for (i in seq(hdr + 1L, length(raw))) {
    body <- strsplit(raw[i], ",")[[1L]]
    if (length(body) < jj) next
    body[jj] <- "99.9"
    raw[i] <- paste(body, collapse = ",")
  }
  writeLines(raw, path)

  expect_error(gsw_tips(maturities = 23:120, cache = path, check = TRUE),
               "disagree; do not use this data")
  expect_silent(gsw_tips(maturities = 23:120, cache = path, check = FALSE))
})

test_that("a check that could only pass vacuously is refused", {
  path <- withr::local_tempfile(fileext = ".csv")
  tips_fixture(path)

  # 25:30 months contains no 24-month point, so the two-year column cannot be
  # compared -- and 60 and 120 are off the grid entirely. Nothing is
  # checkable, which must be an error rather than a silent pass.
  expect_error(gsw_tips(maturities = 25:30, cache = path, check = TRUE),
               "could only pass vacuously")

  # 23:24 does contain the 24-month point, so it is checkable and passes.
  expect_silent(gsw_tips(maturities = 23:24, cache = path, check = TRUE))
})

test_that("align_to relabels months and drops those it cannot match", {
  path <- withr::local_tempfile(fileext = ".csv")
  tips_fixture(path)

  target <- as.Date(c("2003-05-31", "2003-06-30"))
  expect_message(
    out <- gsw_tips(maturities = 23:120, cache = path, align_to = target,
                    check = FALSE),
    "dropped 1 month"
  )
  expect_equal(sort(unique(out$date)), target)
})


# Degenerate inputs the estimator must refuse rather than mangle -----------

test_that("a panel carrying pricing error is still fitted, less exactly", {
  # The realistic case: real curves always have fitting error, and the ACMY
  # likelihood is built around modelling it. Recovery degrades from exact to
  # roughly the size of the error, and must not fall apart.
  s <- sim_acmy(tt = 220L, noise_bp = 0.25)
  f <- suppressWarnings(suppressMessages(
    atsm_real(sim_panel(s), inflation = sim_cpi(s), short_rate = sim_sr(s),
              short_rate_units = "percent")))

  expect_lt(sqrt(mean((f$observed - f$fitted)^2)) * 1e4, 3)
  expect_equal(f$spectral_radius[["risk_adjusted"]],
               max(Mod(eigen(s$truth$phi_tilde,
                             only.values = TRUE)$values)),
               tolerance = 0.02)
  expect_equal(f$breakeven, f$expected_inflation + f$inflation_risk_premium,
               tolerance = 1e-15)
})

test_that("a non positive definite error covariance is reported as such", {
  bad <- matrix(1, 4L, 4L)          # rank one
  expect_error(sigma_e_whitener(bad), "not positive definite")
  expect_error(sigma_e_whitener(bad), "return series")
})

test_that("the whitener turns GLS into ordinary least squares", {
  set.seed(23)
  n <- 6L
  sigma_e <- crossprod(matrix(stats::rnorm(n * n * 3L), n * 3L, n)) / (n * 3L)
  whiten <- sigma_e_whitener(sigma_e)

  b <- matrix(stats::rnorm(n * 3L), n, 3L)
  y <- stats::rnorm(n)

  # Whitened OLS and the explicit GLS normal equations agree when the latter
  # is well conditioned; the point of whitening is that it keeps working when
  # the latter is not.
  gls <- drop(solve(crossprod(b, solve(sigma_e, b)),
                    crossprod(b, solve(sigma_e, y))))
  expect_equal(drop(qr.solve(whiten(b), whiten(y))), gls, tolerance = 1e-8)
})
