# Helpers -----------------------------------------------------------------

standard_tenors <- c(3, 6, 12, 24, 36, 60, 84, 120)

# Curves generated FROM known parameters, observed at a handful of tenors, so
# that a fit has an exact answer to recover.
sim_sparse <- function(tenors = standard_tenors, n = 24L, beta3 = TRUE,
                       seed = 1L) {
  set.seed(seed)
  params <- data.frame(
    beta0 = 0.045 + rnorm(n, sd = 0.003),
    beta1 = -0.020 + rnorm(n, sd = 0.003),
    beta2 = -0.010 + rnorm(n, sd = 0.004),
    beta3 = if (beta3) 0.015 + rnorm(n, sd = 0.004) else 0,
    tau1 = stats::runif(n, 1, 3),
    tau2 = if (beta3) stats::runif(n, 6, 12) else NA_real_
  )
  dates <- seq(as.Date("2010-02-01"), by = "month", length.out = n) - 1
  y <- svensson_curve(params, tenors / 12)
  list(
    panel = yield_panel(y, dates = dates, maturities = tenors,
                        units = "decimal", maturity_unit = "months"),
    params = params
  )
}

# The last ten years of the US curve, observed only at the standard tenors.
gsw_sparse <- function() {
  x <- gsw_monthly[gsw_monthly$maturity %in% standard_tenors &
                   gsw_monthly$date >= as.Date("2015-01-01"), ]
  yield_panel(x, units = "percent", maturity_unit = "months", issuer = "US")
}

warnings_from <- function(expr) {
  w <- character()
  withCallingHandlers(expr, warning = function(e) {
    w <<- c(w, conditionMessage(e))
    invokeRestart("muffleWarning")
  })
  w
}


# Recovery ----------------------------------------------------------------

test_that("a Svensson curve is recovered from eight tenors", {
  s <- sim_sparse()
  fit <- svensson_fit(s$panel, yield_type = "zero", check = FALSE)
  expect_equal(fit$curves$yields$model, "svensson")

  dense <- curve_matrix(predict(fit, maturities = 1:120))
  truth <- svensson_curve(s$params, (1:120) / 12)

  # Not exact, and not because of the optimiser. Eight tenors leave Svensson
  # two residual degrees of freedom, and distinct decay pairs then fit them
  # equally well: on one of these dates decays of (2.75, 0.81) match the
  # tenors to 0.008bp against a truth of (1.22, 7.58). What is recovered is
  # the curve, to a small fraction of a basis point.
  expect_lt(max(abs(dense - truth)), 1e-5)
})

test_that("a Nelson-Siegel curve is recovered with the decay estimated", {
  s <- sim_sparse(tenors = c(3, 12, 24, 60, 120), beta3 = FALSE)
  fit <- svensson_fit(s$panel, yield_type = "zero", check = FALSE)
  expect_equal(fit$curves$yields$model, "nelson_siegel")

  dense <- curve_matrix(predict(fit, maturities = 3:120))
  truth <- svensson_curve(s$params, (3:120) / 12)
  expect_lt(max(abs(dense - truth)), 1e-6)
})

test_that("fixed decays are linear and recover their betas exactly", {
  s <- sim_sparse(beta3 = FALSE)
  s$params$tau1 <- 1.5
  y <- svensson_curve(s$params, standard_tenors / 12)
  panel <- yield_panel(y, dates = s$panel$dates, maturities = standard_tenors,
                       units = "decimal", maturity_unit = "months")

  fit <- svensson_fit(panel, yield_type = "zero", tau = 1.5, check = FALSE)
  p <- coef(fit)
  expect_equal(p$beta0, s$params$beta0, tolerance = 1e-10)
  expect_equal(p$beta1, s$params$beta1, tolerance = 1e-10)
  expect_equal(p$beta2, s$params$beta2, tolerance = 1e-10)
  expect_true(all(p$beta3 == 0))
  expect_true(all(p$tau1 == 1.5))
})

test_that("coef() feeds svensson_curve() and reproduces predict()", {
  s <- sim_sparse()
  fit <- svensson_fit(s$panel, yield_type = "zero", check = FALSE)

  via_coef <- svensson_curve(coef(fit), (1:120) / 12)
  via_predict <- curve_matrix(predict(fit, maturities = 1:120))
  expect_equal(unname(via_coef), unname(via_predict))
})

test_that("a real curve is rebuilt from eight tenors to within a basis point", {
  panel <- gsw_sparse()
  fit <- svensson_fit(panel, yield_type = "zero", check = FALSE)
  dense <- curve_matrix(predict(fit, maturities = 1:120))

  truth <- gsw_monthly[gsw_monthly$date >= as.Date("2015-01-01"), ]
  truth <- curve_matrix(yield_panel(truth, units = "percent",
                                    maturity_unit = "months"))
  between <- setdiff(3:120, standard_tenors)
  expect_lt(sqrt(mean((dense - truth)[, between]^2)), 1e-4)
})


# Model choice ------------------------------------------------------------

test_that("auto chooses the model from the fewest tenors on any date", {
  pick <- function(tenors) {
    s <- sim_sparse(tenors = tenors, n = 6L)
    fit <- svensson_fit(s$panel, yield_type = "zero", check = FALSE)
    fit$curves$yields$model
  }
  expect_equal(pick(c(3, 6, 12, 24, 60, 84, 120)), "svensson")
  expect_equal(pick(c(3, 12, 24, 60, 84, 120)), "nelson_siegel")
  expect_equal(pick(c(3, 12, 24, 60, 120)), "nelson_siegel")
  expect_error(pick(c(3, 24, 60, 120)), "as few as 4")
})

test_that("four tenors are enough once the decay is fixed", {
  s <- sim_sparse(tenors = c(3, 24, 60, 120), beta3 = FALSE)
  fit <- svensson_fit(s$panel, yield_type = "zero", tau = 1.5, check = FALSE)
  expect_equal(fit$curves$yields$model, "nelson_siegel")
  expect_false(anyNA(coef(fit)$beta0))
})

test_that("dates with too few tenors are left unfitted, loudly", {
  s <- sim_sparse(n = 6L)
  s$panel$curves$yields[2, c("6", "36")] <- NA

  expect_warning(
    fit <- svensson_fit(s$panel, yield_type = "zero", model = "svensson",
                        check = FALSE),
    "1 date\\(s\\).*left unfitted"
  )
  p <- coef(fit)
  expect_true(is.na(p$beta0[2]))
  expect_false(anyNA(p$beta0[-2]))
  expect_equal(p$n_obs[2], 6L)

  # The unfitted date stays in as a missing row. Dropping it would leave a
  # panel that skips a month without anything noticing.
  dense <- predict(fit, maturities = 1:120)
  expect_equal(dense$dates, s$panel$dates)
  expect_true(all(is.na(curve_matrix(dense)[2, ])))
  expect_error(atsm(dense, n_factors = 3), "missing values")
})

test_that("dates with no observations at all are not reported as failures", {
  s <- sim_sparse(n = 6L)
  s$panel$curves$yields[1:2, ] <- NA   # a curve that starts later

  expect_no_warning(
    fit <- svensson_fit(s$panel, yield_type = "zero", check = FALSE)
  )
  expect_length(fit$curves$yields$dropped, 0L)
  expect_equal(nrow(curve_matrix(predict(fit))), 6L)
})


# Arguments ---------------------------------------------------------------

test_that("the yield type must be declared, and par yields are refused", {
  s <- sim_sparse(n = 3L)
  expect_error(svensson_fit(s$panel), "yield_type = \"zero\"")
  expect_error(svensson_fit(s$panel, yield_type = "par"), "Par yields")
})

test_that("inconsistent decay arguments are refused", {
  s <- sim_sparse(n = 3L)
  fit <- function(...) {
    svensson_fit(s$panel, yield_type = "zero", check = FALSE, ...)
  }

  expect_error(fit(model = "svensson", tau = 1.5), "needs two")
  expect_error(fit(model = "nelson_siegel", tau = c(1, 5)), "needs one")
  expect_error(fit(tau = c(2, 2)), "must differ")
  expect_error(fit(tau = -1), "positive")
  expect_error(fit(tau_range = c(5, 1)), "increasing")
  expect_error(svensson_fit(gsw_monthly, yield_type = "zero"), "yield_panel")
})


# Extrapolation flags -----------------------------------------------------

test_that("maturities outside each date's observed range are flagged", {
  s <- sim_sparse(n = 4L)
  s$panel$curves$yields[3, "3"] <- NA   # shortest tenor on date 3 is 6m

  fit <- svensson_fit(s$panel, yield_type = "zero", check = FALSE)
  dense <- predict(fit, maturities = c(1, 2, 3, 5, 6, 60, 120, 150))
  flags <- dense$extrapolated$yields

  expect_equal(unname(flags[1, ]), c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
                                     FALSE, TRUE))
  expect_equal(unname(flags[3, ]), c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE,
                                     FALSE, TRUE))
})

test_that("the default grid runs from one month to the longest tenor", {
  s <- sim_sparse(n = 3L)
  fit <- svensson_fit(s$panel, yield_type = "zero", check = FALSE)
  expect_equal(predict(fit)$maturities$yields, 1:120)
})

test_that("yield_panel carries an extrapolation column and validates it", {
  panel <- yield_panel(gsw_monthly, units = "percent",
                       maturity_unit = "months", extrapolated = "extrapolated")
  expect_equal(sum(panel$extrapolated$yields), sum(gsw_monthly$extrapolated))
  expect_output(print(panel), "extrapolated")

  bad <- gsw_monthly
  bad$extrapolated <- as.integer(bad$extrapolated)
  expect_error(
    yield_panel(bad, units = "percent", maturity_unit = "months",
                extrapolated = "extrapolated"),
    "logical"
  )
  expect_null(yield_panel(gsw_monthly, units = "percent",
                          maturity_unit = "months")$extrapolated)
})

test_that("several curves are fitted separately and keep their metadata", {
  as_long <- function(p, curve) {
    y <- curve_matrix(p)
    data.frame(date = rep(p$dates, ncol(y)),
               maturity = rep(p$maturities[[1]], each = nrow(y)),
               yield = as.vector(y), curve = curve)
  }
  a <- sim_sparse(seed = 1L)
  b <- sim_sparse(tenors = c(24, 36, 60, 84, 96, 108, 120), seed = 2L)
  long <- rbind(as_long(a$panel, "nominal"), as_long(b$panel, "real"))
  panel <- yield_panel(long, curve = "curve", units = "decimal",
                       instrument = c(nominal = "government", real = "tips"),
                       issuer = "US")

  fit <- svensson_fit(panel, yield_type = "zero", check = FALSE)
  dense <- predict(fit, maturities = list(nominal = 1:120, real = 23:120))

  expect_equal(dense$meta$instrument, c("government", "tips"))
  expect_equal(dense$maturities$real, 23:120)
  expect_true(all(dense$extrapolated$real[, "23"]))
  expect_false(any(dense$extrapolated$real[, "24"]))
  expect_error(predict(fit, maturities = list(nominal = 1:120)),
               "no entry for curve 'real'")
})


# Hold-out check ----------------------------------------------------------

test_that("the hold-out check refits rather than comparing a fit with itself", {
  s <- sim_sparse(n = 12L)
  clean <- svensson_fit(s$panel, yield_type = "zero")$curves$yields$check
  expect_equal(clean$maturity, c(6, 12, 24, 36, 60, 84))
  expect_lt(max(clean$rmse), 1e-5)

  # Shift one tenor by 50bp. A refit without it predicts the smooth curve, so
  # the check must report the shift; one that compared the full fit with its
  # own inputs would report almost nothing.
  s$panel$curves$yields[, "36"] <- s$panel$curves$yields[, "36"] + 0.005
  bumped <- svensson_fit(s$panel, yield_type = "zero")$curves$yields$check
  expect_gt(bumped$rmse[bumped$maturity == 36], 0.004)
})

test_that("print() reports the model, the fit and the check", {
  s <- sim_sparse(n = 6L)
  expect_output(print(svensson_fit(s$panel, yield_type = "zero")),
                "Svensson.*hold-out")
  expect_output(
    print(svensson_fit(s$panel, yield_type = "zero", tau = 1.5,
                       check = FALSE)),
    "at most 3"
  )
})


# What atsm() does with the result ---------------------------------------

test_that("atsm() refuses more factors than a fixed-decay rebuild contains", {
  fit <- svensson_fit(gsw_sparse(), yield_type = "zero", tau = 1.5,
                      check = FALSE)
  dense <- predict(fit, maturities = 1:120)
  sr <- data.frame(date = dense$dates, value = 0.02)

  expect_error(
    atsm(dense, n_factors = 5, short_rate = sr, short_rate_units = "decimal"),
    "only 3 independent direction"
  )
  expect_s3_class(
    suppressWarnings(atsm(dense, n_factors = 3, short_rate = sr,
                          short_rate_units = "decimal")),
    "atsm_fit"
  )
})

test_that("atsm() warns when its default short rate is extrapolated", {
  fit <- svensson_fit(gsw_sparse(), yield_type = "zero",
                      model = "nelson_siegel", check = FALSE)
  dense <- predict(fit, maturities = 1:120)

  w <- warnings_from(atsm(dense, n_factors = 3))
  expect_true(any(grepl("one-month yield.*extrapolated.*120 of 120", w)))

  sr <- data.frame(date = dense$dates, value = 0.02)
  w <- warnings_from(atsm(dense, n_factors = 3, short_rate = sr,
                          short_rate_units = "decimal"))
  expect_false(any(grepl("extrapolated", w)))

  # GSW flags only its long end, so the bundled data does not trigger it.
  us <- yield_panel(gsw_monthly, units = "percent", maturity_unit = "months",
                    extrapolated = "extrapolated")
  w <- warnings_from(atsm(us, n_factors = 5))
  expect_false(any(grepl("one-month yield", w)))
})
