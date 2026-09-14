us_panel <- function() {
  yield_panel(gsw_monthly, units = "percent", maturity_unit = "months",
              instrument = "government", issuer = "US")
}

us_fit <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- atsm(us_panel(), n_factors = 5)
    cache
  }
})


# Factor extraction -------------------------------------------------------

test_that("principal component signs are deterministic", {
  set.seed(1)
  y <- matrix(rnorm(200 * 10), 200, 10)

  a <- acm_factors(y, k = 3)
  b <- acm_factors(y, k = 3)
  expect_equal(a$loadings, b$loadings)

  # Largest-magnitude loading is positive in every component
  for (j in 1:3) {
    expect_gt(a$loadings[which.max(abs(a$loadings[, j])), j], 0)
  }
})

test_that("factors reconstruct the centred panel", {
  y <- as.matrix(matrix(rnorm(100 * 6), 100, 6))
  f <- acm_factors(y, k = 6)
  recon <- f$scores %*% t(f$loadings)
  expect_equal(recon, sweep(y, 2, f$center, "-"), ignore_attr = TRUE)
})

test_that("acm_factors rejects missing data and bad k", {
  y <- matrix(rnorm(50), 10, 5)
  y[3, 2] <- NA
  expect_error(acm_factors(y, 2), "missing values")
  expect_error(acm_factors(matrix(rnorm(50), 10, 5), 99), "between 1 and")
})


# Excess returns ----------------------------------------------------------

test_that("excess returns match their definition", {
  mats <- 1:5
  set.seed(2)
  y <- matrix(runif(20 * 5, 0.001, 0.005), 20, 5)
  p <- -sweep(y, 2, mats, "*")
  colnames(p) <- mats

  rx <- acm_excess_returns(p, mats, 2:5)

  # rx[t+1]^(n-1) = p[t+1]^(n-1) - p[t]^(n) - r[t]
  expect_equal(unname(rx[1, 1]), unname(p[2, 1] - p[1, 2] - (-p[1, 1])))
  expect_equal(unname(rx[5, 3]), unname(p[6, 3] - p[5, 4] - (-p[5, 1])))
  expect_equal(dim(rx), c(19L, 4L))
})

test_that("excess returns need the predecessor maturity on the grid", {
  mats <- c(1, 12, 24)
  p <- matrix(-0.01, 10, 3, dimnames = list(NULL, mats))
  expect_error(acm_excess_returns(p, mats, c(12, 24)), "both n and n-1")
})


# Recursion ---------------------------------------------------------------

fake_pars <- function(k = 2, sigma2 = 1e-8) {
  list(
    mu = rep(0.001, k),
    phi = diag(0.95, k),
    sigma = diag(1e-6, k),
    sigma2 = sigma2,
    lambda0 = rep(0.01, k),
    lambda1 = diag(0.02, k),
    delta0 = 0.004,
    delta1 = rep(0.001, k)
  )
}

test_that("the recursion is seeded so that a one-period bond is riskless", {
  # A1 = -delta0 and B1 = -delta1 exactly: the sigma2 pricing-error term must
  # NOT enter at n = 1, because a one-period bond held one period is riskless.
  p <- fake_pars(sigma2 = 0.5)
  rec <- acm_recursion(5, p)

  expect_equal(rec$a[1], -p$delta0)
  expect_equal(rec$b[1, ], -p$delta1)
})

test_that("with sigma2 = 0 the recursion is the textbook affine one", {
  p <- fake_pars(sigma2 = 0)
  rec <- acm_recursion(8, p)

  a <- numeric(8); b <- matrix(0, 8, 2)
  a[1] <- -p$delta0; b[1, ] <- -p$delta1
  for (n in 2:8) {
    bp <- b[n - 1, ]
    a[n] <- a[n - 1] + sum(bp * (p$mu - p$lambda0)) +
      0.5 * drop(crossprod(bp, p$sigma %*% bp)) - p$delta0
    b[n, ] <- drop(bp %*% (p$phi - p$lambda1)) - p$delta1
  }

  expect_equal(rec$a, a)
  expect_equal(rec$b, b)
})

test_that("sigma2 shifts only the intercepts, never the loadings", {
  # The pricing-error term enters A but not B, so it moves the level of fitted
  # yields without changing their factor sensitivity.
  r0 <- acm_recursion(10, fake_pars(sigma2 = 0))
  r1 <- acm_recursion(10, fake_pars(sigma2 = 1e-4))

  expect_equal(r0$b, r1$b)
  expect_false(isTRUE(all.equal(r0$a, r1$a)))

  # and the shift accumulates linearly in sigma2 from n = 2 onwards
  expect_equal(r1$a[2] - r0$a[2], 0.5e-4)
})

test_that("sigma2 cancels out of the term premium entirely", {
  # sigma2 enters A_n and A^RF_n identically, so it shifts the level of fitted
  # and risk-neutral yields by the same amount and vanishes from their
  # difference. Useful to know when chasing level discrepancies against another
  # implementation: the pricing-error variance cannot be the explanation.
  p0 <- fake_pars(sigma2 = 0)
  p1 <- fake_pars(sigma2 = 1e-3)
  x <- matrix(c(0.01, -0.02, 0.005, 0.001), 2, 2)

  tp <- function(p) {
    acm_yields(acm_recursion(120, p, FALSE), x, 120) -
      acm_yields(acm_recursion(120, p, TRUE), x, 120)
  }
  expect_equal(tp(p0), tp(p1))

  # ... while the fitted level genuinely does move
  lvl <- function(p) acm_yields(acm_recursion(120, p, FALSE), x, 120)
  expect_false(isTRUE(all.equal(lvl(p0), lvl(p1))))
})

test_that("zero prices of risk make the two recursions identical", {
  p <- fake_pars()
  p$lambda0 <- rep(0, 2)
  p$lambda1 <- matrix(0, 2, 2)

  expect_equal(acm_recursion(10, p, risk_neutral = FALSE),
               acm_recursion(10, p, risk_neutral = TRUE))
})

test_that("spectral_radius detects explosive dynamics", {
  expect_lt(spectral_radius(diag(0.9, 3)), 1)
  expect_gt(spectral_radius(diag(1.1, 3)), 1)
  expect_equal(spectral_radius(matrix(c(0, -1, 1, 0), 2)), 1)
})


# End-to-end --------------------------------------------------------------

test_that("the decomposition is exact", {
  # fitted = risk_neutral + term_premium, to machine precision. This is the
  # identity the published NY Fed series satisfies, and the one bug in this
  # literature that most often goes unnoticed.
  f <- us_fit()
  expect_lt(max(abs(f$fitted - (f$risk_neutral + f$term_premium))), 1e-12)
})

test_that("the model fits observed yields closely", {
  f <- us_fit()
  rmse_bp <- sqrt(mean((f$observed - f$fitted)^2)) * 1e4
  # BIS report single-digit basis point residuals for models of this class
  expect_lt(rmse_bp, 10)
})

test_that("estimated dynamics are stationary on US data", {
  f <- us_fit()
  expect_lt(f$spectral_radius[["risk_adjusted"]], 1)
  expect_lt(f$spectral_radius[["real_world"]], 1)
})

test_that("the US term premium reproduces the published NY Fed series", {
  # The credibility check. Levels differ by a roughly constant offset, but the
  # dynamics must track very closely; see the M4 notes on the remaining gap.
  f <- us_fit()
  mine <- term_premium(f, maturity = 120)
  pub <- acm_published[acm_published$maturity == 120, c("date", "term_premium")]

  m <- merge(mine, pub, by = "date")
  expect_equal(nrow(m), 763L)

  mine_bp <- m$value * 1e4
  pub_bp <- m$term_premium * 100

  expect_gt(cor(mine_bp, pub_bp), 0.99)
  expect_gt(cor(diff(mine_bp), diff(pub_bp)), 0.98)
  expect_lt(abs(mean(mine_bp - pub_bp)), 30)
})

test_that("the long-run mean term premium is plausible", {
  # BIS report ~160bp for the US over 1961-2018; our sample runs to 2024 and
  # includes the low-premium late 2010s, so a somewhat lower figure is right.
  f <- us_fit()
  tp10 <- f$term_premium[, match(120, f$maturities)] * 1e4
  expect_gt(mean(tp10), 80)
  expect_lt(mean(tp10), 200)
})

test_that("extractors return tidy frames and validate maturities", {
  f <- us_fit()

  tp <- term_premium(f, maturity = c(24, 120))
  expect_named(tp, c("date", "maturity", "value"))
  expect_equal(nrow(tp), 2L * length(f$dates))
  expect_setequal(unique(tp$maturity), c(24, 120))

  expect_error(term_premium(f, maturity = 999), "not fitted")
  expect_equal(risk_neutral(f, 120), expected_short_rate(f, 120))
})

test_that("residuals are observed minus fitted", {
  f <- us_fit()
  r <- residuals(f)
  expect_equal(r$value,
               fitted(f)$value * -1 + extract_component(f, "observed", NULL)$value)
})

test_that("predict reproduces the fit on its own estimation sample", {
  # Applying the fitted parameters to the same panel must return the same
  # decomposition; this is what makes higher-frequency evaluation valid.
  f <- us_fit()
  p <- predict(f, us_panel())

  expect_equal(nrow(p), length(f$dates) * length(f$maturities))
  expect_equal(p$term_premium, as.vector(f$term_premium), tolerance = 1e-10)
  expect_equal(p$fitted, as.vector(f$fitted), tolerance = 1e-10)
})

test_that("atsm validates its inputs", {
  panel <- us_panel()

  expect_error(atsm("not a panel"), "must be a yield_panel")
  expect_error(atsm(panel, maturities = c(1, 999)), "not on the panel")
  expect_error(atsm(panel, maturities = c(12, 24, 36)), "one-month maturity")
})

test_that("print and summary run", {
  f <- us_fit()
  expect_output(print(f), "atsm_fit")
  expect_output(print(f), "ACM")
  expect_output(print(summary(f)), "By maturity")
})
