# Forward rates, forward curves, expected excess returns, and the joint fit's
# predict() and plot() methods.
#
# The forward tests lean on one exact identity. A forward rate is a difference
# of two log prices, and a fitted yield is a log price divided by its
# maturity, so
#
#     f(m1 -> m2) * (m2 - m1) = m2 * y(m2) - m1 * y(m1)
#
# holds with no expectation taken and no estimation involved. Anything beyond
# machine precision there is a bug, and it ties the forward code to the yield
# code rather than letting the two drift apart.

us_small <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      p <- yield_panel(gsw_monthly, units = "percent",
                       maturity_unit = "months", issuer = "US")
      cache <<- suppressWarnings(atsm(p, n_factors = 3L))
    }
    cache
  }
})


# Forward rates -----------------------------------------------------------

test_that("a forward rate is the difference of the two log prices", {
  f <- us_small()
  fw <- forward_rate(f, start = 60L, end = 120L)
  got <- fw$value[fw$component == "fitted"]

  y60 <- f$fitted[, match(60L, f$maturities)]
  y120 <- f$fitted[, match(120L, f$maturities)]
  expect_equal(got, (120 * y120 - 60 * y60) / 60, tolerance = 1e-14,
               ignore_attr = TRUE)
})

test_that("the risk-neutral forward comes from the risk-neutral curve", {
  f <- us_small()
  fw <- forward_rate(f, start = 24L, end = 60L)
  got <- fw$value[fw$component == "risk_neutral"]

  r24 <- f$risk_neutral[, match(24L, f$maturities)]
  r60 <- f$risk_neutral[, match(60L, f$maturities)]
  expect_equal(got, (60 * r60 - 24 * r24) / 36, tolerance = 1e-14,
               ignore_attr = TRUE)
})

test_that("the forward decomposition adds up exactly", {
  f <- us_small()
  fw <- forward_rate(f, start = 60L, end = 120L)
  w <- split(fw$value, fw$component)

  expect_equal(w$fitted, w$risk_neutral + w$term_premium, tolerance = 1e-15)
  expect_equal(unique(fw$start), 60L)
  expect_equal(unique(fw$end), 120L)
  expect_equal(nrow(fw), 3L * length(f$dates))
})

test_that("a one-month forward window agrees with the forward curve", {
  # The two code paths compute the same object from different directions:
  # forward_rate() differences recursion coefficients, forward_curve()
  # differences fitted yields.
  f <- us_small()
  d <- f$dates[100L]
  fw <- forward_rate(f, start = 59L, end = 60L)
  point <- fw$value[fw$component == "fitted"][100L]

  fc <- forward_curve(f, date = d)
  curve_point <- fc$value[fc$component == "fitted" & fc$maturity == 60L]

  expect_equal(point, curve_point, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("forward windows are validated", {
  f <- us_small()

  expect_error(forward_rate(f, start = 120L, end = 60L), "must be before")
  expect_error(forward_rate(f, start = 60L, end = 60L), "must be before")
  expect_error(forward_rate(f, start = 60L, end = 999L), "not fitted")
  expect_error(forward_rate(f, start = 0L, end = 60L), "at least 1")
  expect_error(forward_rate(f, start = 1.5, end = 60L), "whole number")
  expect_error(forward_rate(f, start = c(12L, 24L), end = 60L), "single whole")
})

test_that("the forward curve recovers the observed curve it came from", {
  f <- us_small()
  fc <- forward_curve(f, date = f$dates[50L])
  obs <- fc[fc$component == "observed", ]

  y <- f$observed[50L, ]
  mats <- f$maturities
  j <- match(obs$maturity, mats)
  expect_equal(obs$value,
               mats[j] * y[j] - (mats[j] - 1L) * y[j - 1L],
               tolerance = 1e-14, ignore_attr = TRUE)

  # The first maturity has no predecessor, so it is absent rather than NA.
  expect_false(min(mats) %in% obs$maturity)
  expect_equal(min(obs$maturity), min(mats) + 1L)
})

test_that("the forward curve defaults to the last date and validates others", {
  f <- us_small()
  expect_equal(nrow(forward_curve(f)), nrow(forward_curve(f, max(f$dates))))
  expect_error(forward_curve(f, "1800-01-01"), "No fitted observation")
})


# Forwards on the joint fit -----------------------------------------------

test_that("the joint forward decomposition keeps every identity", {
  f <- sim_fit()$fit
  fw <- forward_rate(f, start = 60L, end = 120L)
  w <- split(fw$value, fw$component)

  expect_equal(w$fitted, w$risk_neutral + w$term_premium, tolerance = 1e-15)
  expect_equal(w$fitted_real,
               w$risk_neutral_real + w$term_premium_real, tolerance = 1e-15)
  expect_equal(w$breakeven,
               w$expected_inflation + w$inflation_risk_premium,
               tolerance = 1e-15)
  expect_equal(w$inflation_risk_premium,
               w$term_premium - w$term_premium_real, tolerance = 1e-15)
  expect_equal(w$breakeven, w$fitted - w$fitted_real, tolerance = 1e-15)
})

test_that("a joint forward matches the two curves it is built from", {
  f <- sim_fit()$fit
  fw <- forward_rate(f, start = 60L, end = 120L)
  w <- split(fw$value, fw$component)

  jr <- match(c(60L, 120L), f$real_maturities)
  expect_equal(w$fitted_real,
               (120 * f$fitted_real[, jr[2L]] - 60 * f$fitted_real[, jr[1L]]) / 60,
               tolerance = 1e-14, ignore_attr = TRUE)
})

test_that("a window off the real grid is refused by name", {
  f <- sim_fit()$fit
  # 12 months is on the nominal grid but the real curve starts at 23.
  expect_error(forward_rate(f, start = 12L, end = 120L),
               "not fitted on the real grid")
})

test_that("liquidity components appear only when the fit has them", {
  f <- sim_fit()$fit
  expect_false(f$has_liquidity)
  expect_false("liquidity_nominal" %in% forward_rate(f, 60L, 120L)$component)

  s <- sim_acmy(liquidity = TRUE, noise_bp = 0.25)
  g <- suppressWarnings(suppressMessages(
    atsm_real(sim_panel(s), inflation = sim_cpi(s), short_rate = sim_sr(s),
              short_rate_units = "percent",
              liquidity = data.frame(date = s$dates, value = s$liquidity))))
  fw <- forward_rate(g, start = 60L, end = 120L)
  w <- split(fw$value, fw$component)

  expect_true(all(c("liquidity_nominal", "liquidity_real",
                    "breakeven_liquidity_adjusted") %in% fw$component))
  # Section 2.2: the adjusted breakeven is the breakeven net of the liquidity
  # factor's contribution to each curve.
  expect_equal(w$breakeven_liquidity_adjusted,
               w$breakeven - (w$liquidity_nominal - w$liquidity_real),
               tolerance = 1e-15)
})

test_that("the joint forward curve carries both curves", {
  f <- sim_fit()$fit
  fc <- forward_curve(f)
  expect_true(all(c("observed", "fitted", "risk_neutral", "observed_real",
                    "fitted_real", "risk_neutral_real") %in% fc$component))
  expect_equal(min(fc$maturity[fc$component == "fitted_real"]),
               min(f$real_maturities) + 1L)
})


# Expected excess returns -------------------------------------------------

test_that("a one-month bond has exactly zero expected excess return", {
  # It is riskless over the holding period, so B_0 = 0 and the whole
  # expression vanishes. This is what pins the indexing: the loading is
  # B_{n-1}, for the bond one month closer to maturity, not B_n.
  f <- us_small()
  er <- expected_excess_return(f, maturity = 1L)
  expect_equal(er$value, rep(0, length(f$dates)), tolerance = 0,
               ignore_attr = TRUE)
})

test_that("expected excess returns rise with maturity and are affine in X", {
  f <- us_small()
  er <- expected_excess_return(f, maturity = c(12L, 60L, 120L))
  m <- tapply(er$value, er$maturity, mean)
  expect_true(all(diff(m) > 0))

  # Affine in the factors: regressing on them leaves nothing.
  v <- er$value[er$maturity == 120L]
  r <- stats::lm(v ~ f$factors)
  expect_lt(max(abs(stats::residuals(r))), 1e-12)
})

test_that("risk compensation and convexity add to the total", {
  f <- us_small()
  er <- expected_excess_return(f, maturity = c(24L, 120L), convexity = TRUE)
  w <- split(er$value, er$component)

  expect_equal(w$expected, w$risk_compensation + w$convexity,
               tolerance = 1e-15)
  # The convexity adjustment is a subtraction and does not vary over time.
  expect_true(all(w$convexity < 0))
  conv120 <- er$value[er$component == "convexity" & er$maturity == 120L]
  expect_lt(stats::sd(conv120), 1e-15)
})

test_that("the loadings are the risk-price exposure scaled by factor sd", {
  f <- us_small()
  ld <- expected_excess_return_loadings(f, maturity = 120L)

  b <- f$recursion$p$b[119L, ]
  want <- drop(b %*% f$pars$lambda1) * apply(f$factors, 2L, stats::sd)
  expect_equal(ld$value, want, tolerance = 1e-14, ignore_attr = TRUE)
  expect_equal(nrow(ld), f$n_factors)

  expect_equal(expected_excess_return_loadings(f, maturity = 1L)$value,
               rep(0, f$n_factors), tolerance = 0, ignore_attr = TRUE)
})

test_that("expected excess returns work on a joint fit, on either curve", {
  f <- sim_fit()$fit

  n <- expected_excess_return(f, maturity = 120L)
  r <- expected_excess_return(f, maturity = 120L, curve = "real")
  expect_equal(nrow(n), length(f$dates))
  expect_equal(nrow(r), length(f$dates))
  expect_false(isTRUE(all.equal(n$value, r$value)))

  expect_error(expected_excess_return(f, maturity = 12L, curve = "real"),
               "not fitted")
  expect_error(expected_excess_return(f, maturity = 999L), "not fitted")
})


# predict() on the joint fit ----------------------------------------------

test_that("predicting on the estimation sample reproduces the fit", {
  # The strongest available check on project_real_factors(): it must replay
  # the factor construction exactly, or the parameters are being applied to
  # factors on a different basis from the one they were estimated against.
  d <- sim_fit()
  f <- d$fit
  p <- sim_panel(d$s)

  out <- predict(f, p)

  # predict() reaches one date FURTHER BACK than the fit, and should: the
  # decomposition at t needs only the factors at t, while estimation also
  # needs the previous month's price index to form inflation. Compare on the
  # dates they share.
  expect_gt(length(unique(out$date)), length(f$dates))
  expect_true(all(f$dates %in% out$date))

  j <- out$maturity == 120L & out$date %in% f$dates
  expect_equal(out$fitted[j], f$fitted[, match(120L, f$maturities)],
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(out$term_premium[j],
               f$term_premium[, match(120L, f$maturities)],
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(out$breakeven[j],
               f$breakeven[, match(120L, f$shared_maturities)],
               tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("maturities the real curve does not reach come back as NA", {
  d <- sim_fit()
  out <- predict(d$fit, sim_panel(d$s))

  short <- out[out$maturity == 12L, ]
  expect_true(all(is.na(short$fitted_real)))
  expect_true(all(is.na(short$breakeven)))
  expect_false(anyNA(short$fitted))
})

test_that("predict validates what it is given", {
  d <- sim_fit()
  expect_error(predict(d$fit, "not a panel"), "must be a yield_panel")

  single <- yield_panel(gsw_monthly, units = "percent",
                        maturity_unit = "months", issuer = "US")
  expect_error(predict(d$fit, single), "has no curve named")
})

test_that("a liquidity-fitted model demands liquidity for new dates", {
  s <- sim_acmy(liquidity = TRUE, noise_bp = 0.25)
  liq <- data.frame(date = s$dates, value = s$liquidity)
  f <- suppressWarnings(suppressMessages(
    atsm_real(sim_panel(s), inflation = sim_cpi(s), short_rate = sim_sr(s),
              short_rate_units = "percent", liquidity = liq)))

  expect_error(predict(f, sim_panel(s)), "needs one for the new dates")

  out <- predict(f, sim_panel(s), liquidity = liq)
  j <- out$maturity == 120L & out$date %in% f$dates
  expect_equal(out$fitted[j], f$fitted[, match(120L, f$maturities)],
               tolerance = 1e-12, ignore_attr = TRUE)
})


# plot() on the joint fit -------------------------------------------------

test_that("every plot type runs", {
  f <- sim_fit()$fit
  tmp <- withr::local_tempfile(fileext = ".png")

  for (ty in c("inflation", "real", "decomposition", "premium", "residuals")) {
    grDevices::png(tmp)
    expect_silent(plot(f, type = ty))
    grDevices::dev.off()
  }
})

test_that("plot refuses a maturity the chosen type does not have", {
  f <- sim_fit()$fit
  tmp <- withr::local_tempfile(fileext = ".png")
  grDevices::png(tmp)
  on.exit(grDevices::dev.off(), add = TRUE)

  # 12 months exists on the nominal curve but there is no breakeven for it.
  expect_error(plot(f, maturity = 12L, type = "inflation"), "not available")
  expect_silent(plot(f, maturity = 12L, type = "premium"))
})
