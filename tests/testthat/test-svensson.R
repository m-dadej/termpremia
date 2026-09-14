test_that("basis functions have the right limits at zero", {
  # (1 - exp(-x))/x -> 1 as x -> 0; the expression is 0/0 there
  expect_equal(ns_level_factor(0, tau = 2), 1)
  expect_equal(ns_hump_factor(0, tau = 2), 0)

  # and are continuous approaching zero
  expect_equal(ns_level_factor(1e-12, tau = 2), 1, tolerance = 1e-9)
  expect_lt(abs(ns_hump_factor(1e-12, tau = 2)), 1e-9)
})

test_that("yields converge to the level parameter at long maturities", {
  y <- function(n) {
    svensson_yield(n, beta0 = 4, beta1 = -2, beta2 = 1, beta3 = 0.5,
                   tau1 = 1.5, tau2 = 10)
  }

  # The basis functions decay as tau/n, so convergence to beta0 is O(1/n) --
  # slow enough that exact equality is the wrong assertion. Test the limiting
  # behaviour instead: close at large n, and monotonically closer beyond that.
  expect_equal(y(1e5), 4, tolerance = 1e-4)
  expect_lt(abs(y(1e7) - 4), abs(y(1e5) - 4))
})

test_that("zero-maturity limit is beta0 + beta1", {
  # f(0) = 1 and the hump factors vanish, so y(0) = beta0 + beta1
  y <- svensson_yield(0, beta0 = 4, beta1 = -2, beta2 = 1, beta3 = 0.5,
                      tau1 = 1.5, tau2 = 10)
  expect_equal(y, 2)
})

test_that("beta3 = 0 reduces Svensson to Nelson-Siegel", {
  mats <- c(0.25, 1, 5, 10, 30)
  sv <- svensson_yield(mats, beta0 = 4, beta1 = -2, beta2 = 1, beta3 = 0,
                       tau1 = 1.5, tau2 = 10)
  ns <- svensson_yield(mats, beta0 = 4, beta1 = -2, beta2 = 1, tau1 = 1.5)
  expect_equal(sv, ns)
})

test_that("GSW's pre-1980 tau2 sentinel is handled, not propagated", {
  # GSW encode "Nelson-Siegel only" as beta3 = 0 with tau2 = -999.99.
  # The fourth term must be skipped rather than evaluated.
  mats <- c(1, 5, 10)
  sentinel <- svensson_yield(mats, beta0 = 4, beta1 = -2, beta2 = 1,
                             beta3 = 0, tau1 = 1.5, tau2 = -999.99)
  plain <- svensson_yield(mats, beta0 = 4, beta1 = -2, beta2 = 1, tau1 = 1.5)

  expect_equal(sentinel, plain)
  expect_false(anyNA(sentinel))
})

test_that("a non-zero beta3 with an invalid tau2 is an error", {
  # Silently ignoring this would mean quietly dropping a real curve component
  expect_error(
    svensson_yield(1:5, beta0 = 4, beta1 = -2, beta2 = 1, beta3 = 0.5,
                   tau1 = 1.5, tau2 = -999.99),
    "tau2"
  )
  expect_error(
    svensson_yield(1:5, beta0 = 4, beta1 = -2, beta2 = 1, beta3 = 0.5,
                   tau1 = 1.5, tau2 = NA),
    "tau2"
  )
})

test_that("missing parameter rows yield NA rather than erroring", {
  y <- svensson_yield(1:5, beta0 = NA, beta1 = -2, beta2 = 1, tau1 = 1.5)
  expect_true(all(is.na(y)))
  expect_length(y, 5)
})

test_that("invalid maturities and tau1 are rejected", {
  expect_error(svensson_yield(-1, 4, -2, 1, tau1 = 1.5), "non-negative")
  expect_error(svensson_yield(NA, 4, -2, 1, tau1 = 1.5), "NA")
  expect_error(svensson_yield(Inf, 4, -2, 1, tau1 = 1.5), "finite")
  expect_error(svensson_yield(1, 4, -2, 1, tau1 = 0), "positive")
  expect_error(svensson_yield(1, 4, -2, 1, tau1 = -1), "positive")
})

test_that("svensson_curve evaluates a panel consistently with the scalar form", {
  params <- data.frame(
    beta0 = c(4, 4.1), beta1 = c(-2, -1.9), beta2 = c(1, 1.1),
    beta3 = c(0.5, 0), tau1 = c(1.5, 1.5), tau2 = c(10, -999.99)
  )
  mats <- c(1, 5, 10)

  got <- svensson_curve(params, mats)

  expect_equal(dim(got), c(2L, 3L))
  for (i in 1:2) {
    expect_equal(
      unname(got[i, ]),
      svensson_yield(mats, params$beta0[i], params$beta1[i], params$beta2[i],
                     params$beta3[i], params$tau1[i], params$tau2[i])
    )
  }
})

test_that("svensson_curve handles a single maturity and missing columns", {
  params <- data.frame(beta0 = 4, beta1 = -2, beta2 = 1, beta3 = 0,
                       tau1 = 1.5, tau2 = 10)

  expect_equal(dim(svensson_curve(params, 10)), c(1L, 1L))
  expect_error(svensson_curve(params[, 1:3], 10), "missing required column")
})
