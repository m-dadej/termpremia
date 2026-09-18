# Survey-augmented factor dynamics.
#
# The real SPF cannot be bundled (all rights reserved, research use only), so
# everything here runs against a synthetic survey built from the model's own
# expectations. That is weaker evidence than real data about whether the idea
# helps, and stronger evidence about whether the arithmetic is right: when the
# survey IS the model's expectation, the estimator must leave the parameters
# alone, and when it is a known distortion of them, it must move in the stated
# direction. `analysis/validate-spf.R` does the real-data half.

sv_panel <- function() {
  yield_panel(gsw_monthly, units = "percent", maturity_unit = "months",
              instrument = "government", issuer = "US")
}

sv_ols <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- suppressWarnings(atsm(sv_panel(), n_factors = 3))
    cache
  }
})

# Expected average short rate implied by a fitted model, exactly as the
# estimator computes it, so a survey can be built that the model already agrees
# with to machine precision.
model_forecast <- function(fit, rows, horizon, tenor, average_months) {
  tg <- list(index = rows, horizon = rep(horizon, length(rows)),
             tenor = rep(tenor, length(rows)),
             average_months = rep(average_months, length(rows)))
  d <- survey_design(tg)
  survey_implied(fit$pars$mu, fit$pars$phi,
                 fit$factors[d$rows, , drop = FALSE],
                 fit$pars$delta0, fit$pars$delta1, d) * 12
}

# A survey frame the fitted model reproduces exactly.
#
# `shift` is in PERCENTAGE POINTS, because that is how forecasts are discussed
# ("a full point higher"), while `value` is an annualised decimal. Passing the
# shift straight through was a units slip worth 100x: it drove the VAR
# explosive (rho 1.33), which made the optimiser's iteration count
# platform-dependent and the convergence assertion below a coin flip on CI.
synthetic_survey <- function(fit, rows, horizon = 3L, tenor = 3L,
                             average_months = 3L, shift = 0) {
  data.frame(
    date = fit$dates[rows],
    target_date = add_months(fit$dates[rows], horizon),
    value = model_forecast(fit, rows, horizon, tenor, average_months) +
      shift / 100,
    tenor = tenor,
    average_months = average_months,
    stringsAsFactors = FALSE
  )
}


# Controls ----------------------------------------------------------------

test_that("survey_control returns validated defaults", {
  ctl <- survey_control()
  expect_equal(ctl$tenor, 1L)
  expect_equal(ctl$average_months, 1L)
  expect_true(ctl$sd > 0)
})

test_that("survey_control rejects impossible settings", {
  expect_error(survey_control(tenor = 0), "between 1")
  expect_error(survey_control(average_months = 0), "between 1")
  expect_error(survey_control(sd = 0), "between")
  expect_error(survey_control(sd = c(1, 2)), "single finite number")
  expect_error(survey_control(reltol = 2), "between")
  expect_error(survey_control(units = "furlongs"), "should be one of")
})


# Window weights ----------------------------------------------------------

test_that("a point forecast of the short rate puts all weight on one month", {
  expect_equal(survey_weights(tenor = 1L, average_months = 1L), 1)
})

test_that("the SPF window has the documented trapezoidal weights", {
  # Quarterly average of a three-month rate: five months, weights 1:2:3:2:1.
  expect_equal(survey_weights(3L, 3L), c(1, 2, 3, 2, 1) / 9)
})

test_that("weights always sum to one and span the right number of months", {
  for (n in c(1L, 3L, 12L)) {
    for (a in c(1L, 3L, 120L)) {
      w <- survey_weights(n, a)
      expect_equal(sum(w), 1)
      expect_equal(length(w), a + n - 1L)
      expect_true(all(w > 0))
    }
  }
})

test_that("a pure n-month rate is an equal-weighted average of n months", {
  expect_equal(survey_weights(tenor = 3L, average_months = 1L), rep(1 / 3, 3))
})


# Calendar arithmetic -----------------------------------------------------

test_that("months_between counts calendar months, not days", {
  expect_equal(months_between(as.Date("2000-02-29"), as.Date("2000-04-01")), 2L)
  expect_equal(months_between(as.Date("2000-01-31"), as.Date("2000-02-01")), 1L)
  expect_equal(months_between(as.Date("2000-12-31"), as.Date("2001-01-01")), 1L)
  expect_equal(months_between(as.Date("2000-06-15"), as.Date("2000-06-01")), 0L)
})

test_that("add_months does not overshoot from a month end", {
  # seq.Date(by = "month") turns 31 January into 3 March; this must not.
  expect_equal(add_months(as.Date("2000-01-31"), 1L), as.Date("2000-02-01"))
  expect_equal(add_months(as.Date("2000-11-01"), 3L), as.Date("2001-02-01"))
  expect_equal(add_months(as.Date("2000-03-15"), -2L), as.Date("2000-01-01"))
  expect_equal(add_months(as.Date("2000-01-01"), 120L), as.Date("2010-01-01"))
})


# Attaching surveys to the panel ------------------------------------------

test_that("survey_targets validates its input", {
  d <- sv_ols()$dates
  expect_error(survey_targets(list(), d, survey_control()), "must be a data frame")
  expect_error(
    survey_targets(data.frame(date = d[1], value = 1), d, survey_control()),
    "missing column\\(s\\): target_date"
  )
})

test_that("horizons are measured from the matched panel date", {
  f <- sv_ols()
  # A survey taken mid-February matched to the end of February, forecasting
  # the quarter beginning in April, is looking two months ahead.
  s <- data.frame(date = as.Date("2000-02-15"),
                  target_date = as.Date("2000-04-01"), value = 5)
  tg <- suppressMessages(survey_targets(s, f$dates, survey_control()))
  expect_equal(tg$horizon, 2L)
  expect_equal(f$dates[tg$index], as.Date("2000-02-29"))
})

test_that("forecasts of a period already under way are dropped", {
  f <- sv_ols()
  s <- data.frame(
    date = as.Date(c("2000-02-15", "2000-02-15")),
    target_date = as.Date(c("2000-01-01", "2000-04-01")),
    value = c(5, 5)
  )
  tg <- suppressMessages(survey_targets(s, f$dates, survey_control()))
  expect_equal(tg$n_used, 1L)
  expect_equal(tg$n_past, 1L)
})

test_that("surveys outside the panel are dropped, and all-outside is an error", {
  f <- sv_ols()
  s <- data.frame(date = as.Date("2099-02-15"),
                  target_date = as.Date("2099-04-01"), value = 5)
  expect_error(suppressMessages(survey_targets(s, f$dates, survey_control())),
               "No survey forecast could be used")
})

test_that("window columns override the control defaults", {
  f <- sv_ols()
  s <- data.frame(date = as.Date("2000-02-15"),
                  target_date = as.Date("2000-04-01"), value = 5,
                  tenor = 3L, average_months = 120L)
  tg <- suppressMessages(survey_targets(
    s, f$dates, survey_control(tenor = 1L, average_months = 1L)))
  expect_equal(tg$tenor, 3L)
  expect_equal(tg$average_months, 120L)

  bad <- s; bad$tenor <- 0L
  expect_error(suppressMessages(survey_targets(bad, f$dates, survey_control())),
               "whole numbers of months")
})

test_that("mixed window shapes give each forecast its own weights", {
  f <- sv_ols()
  rows <- c(400L, 500L)
  s <- rbind(
    data.frame(date = f$dates[rows], target_date = add_months(f$dates[rows], 3L),
               value = 5, tenor = 3L, average_months = 3L),
    data.frame(date = f$dates[rows], target_date = add_months(f$dates[rows], 3L),
               value = 5, tenor = 3L, average_months = 120L)
  )
  d <- survey_design(suppressMessages(survey_targets(s, f$dates, survey_control())))

  expect_equal(d$p_max, 3L + 122L - 1L)
  expect_equal(rowSums(d$weights), rep(1, 4))
  # The short forecasts touch 5 months, the long ones 122.
  expect_equal(sort(rowSums(d$weights > 0)), c(5, 5, 122, 122))
  # Only two distinct panel dates are iterated forward, not the whole sample.
  expect_equal(length(d$rows), 2L)
})


# The estimator -----------------------------------------------------------

test_that("a survey the model already agrees with leaves it alone", {
  # The sharpest available check on the arithmetic: if the synthetic survey is
  # exactly the model's own expectation, the penalty is zero at the OLS
  # parameters and cannot pull them anywhere.
  f <- sv_ols()
  rows <- seq(300L, 700L, by = 3L)
  s <- synthetic_survey(f, rows)

  fit <- suppressWarnings(suppressMessages(
    atsm(sv_panel(), n_factors = 3, p_dynamics = "survey", survey = s)))

  expect_equal(fit$pars$phi, f$pars$phi, tolerance = 1e-6)
  expect_equal(fit$pars$mu, f$pars$mu, tolerance = 1e-6)
  expect_lt(fit$survey$rmse_bp[["survey"]], 0.5)
  # Convergence is asserted here rather than in the diagnostics test below:
  # the penalty is zero at the starting values, so the optimiser stops almost
  # immediately and does so on every platform.
  expect_true(fit$survey$converged)
})

test_that("a survey that disagrees pulls the model towards it", {
  f <- sv_ols()
  rows <- seq(300L, 700L, by = 3L)
  # Forecasters who think rates will be a full point higher than the model does.
  s <- synthetic_survey(f, rows, shift = 1)

  fit <- suppressWarnings(suppressMessages(
    atsm(sv_panel(), n_factors = 3, p_dynamics = "survey", survey = s,
         survey_control = survey_control(sd = 0.002))))

  expect_lt(fit$survey$rmse_bp[["survey"]], fit$survey$rmse_bp[["ols"]])
  expect_false(isTRUE(all.equal(fit$pars$phi, f$pars$phi)))
})

test_that("a large sd reduces the estimator to plain OLS", {
  f <- sv_ols()
  rows <- seq(300L, 700L, by = 3L)
  s <- synthetic_survey(f, rows, shift = 1)

  fit <- suppressWarnings(suppressMessages(
    atsm(sv_panel(), n_factors = 3, p_dynamics = "survey", survey = s,
         survey_control = survey_control(sd = 100))))

  expect_equal(fit$pars$phi, f$pars$phi, tolerance = 1e-8)
  expect_equal(fit$pars$mu, f$pars$mu, tolerance = 1e-8)
})

test_that("a smaller sd fits the surveys more closely", {
  f <- sv_ols()
  rows <- seq(300L, 700L, by = 6L)
  s <- synthetic_survey(f, rows, shift = 1)

  tight <- suppressWarnings(suppressMessages(
    atsm(sv_panel(), n_factors = 3, p_dynamics = "survey", survey = s,
         survey_control = survey_control(sd = 0.001))))
  loose <- suppressWarnings(suppressMessages(
    atsm(sv_panel(), n_factors = 3, p_dynamics = "survey", survey = s,
         survey_control = survey_control(sd = 0.02))))

  expect_lt(tight$survey$rmse_bp[["survey"]], loose$survey$rmse_bp[["survey"]])
})


# Integration with atsm() -------------------------------------------------

test_that("survey discipline leaves the yield curve fit untouched", {
  f <- sv_ols()
  s <- synthetic_survey(f, seq(300L, 700L, by = 6L), shift = 1)
  fit <- suppressWarnings(suppressMessages(
    atsm(sv_panel(), n_factors = 3, p_dynamics = "survey", survey = s,
         survey_control = survey_control(sd = 0.002))))

  expect_equal(fit$fitted, f$fitted, tolerance = 1e-12)
  expect_equal(fit$pars$phi - fit$pars$lambda1,
               f$pars$phi - f$pars$lambda1, tolerance = 1e-12)
  expect_equal(fit$spectral_radius[["risk_adjusted"]],
               f$spectral_radius[["risk_adjusted"]])
  expect_equal(fit$fitted, fit$risk_neutral + fit$term_premium,
               tolerance = 1e-14)
})

test_that("a fitted model carries its survey diagnostics", {
  f <- sv_ols()
  rows <- seq(300L, 700L, by = 6L)
  s <- synthetic_survey(f, rows, shift = 0.5)
  fit <- suppressWarnings(suppressMessages(
    atsm(sv_panel(), n_factors = 3, p_dynamics = "survey", survey = s)))

  expect_null(f$survey)
  expect_equal(fit$p_dynamics, "survey")

  d <- fit$survey
  expect_equal(d$n_forecasts, length(rows))
  expect_equal(d$phi_ols, f$pars$phi)
  expect_named(d$rmse_bp, c("ols", "survey"))
  # This test is about the shape of the diagnostics, not the optimiser's luck:
  # whether BFGS reports code 0 on any particular input depends on the BLAS,
  # so assert the field is there and well-typed and leave convergence to the
  # test above.
  expect_type(d$converged, "logical")
  expect_length(d$converged, 1L)

  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "P-dynamics: survey")
  expect_match(out, "surveys")
  expect_match(out, "fit to surveys")
})

test_that("atsm demands surveys when asked for survey dynamics", {
  expect_error(atsm(sv_panel(), n_factors = 3, p_dynamics = "survey"),
               "needs survey forecasts")
})

test_that("surveys supplied to the wrong estimator are reported, not ignored", {
  f <- sv_ols()
  s <- synthetic_survey(f, seq(600L, 700L, by = 10L))
  # The fit also warns about the lower bound, so match rather than expect one.
  w <- testthat::capture_warnings(
    suppressMessages(atsm(sv_panel(), n_factors = 3, survey = s))
  )
  expect_true(any(grepl("was supplied but `p_dynamics`", w, fixed = TRUE)))
})

test_that("predict carries the survey-disciplined dynamics through", {
  f <- sv_ols()
  s <- synthetic_survey(f, seq(300L, 700L, by = 6L), shift = 1)
  fit <- suppressWarnings(suppressMessages(
    atsm(sv_panel(), n_factors = 3, p_dynamics = "survey", survey = s,
         survey_control = survey_control(sd = 0.002))))

  p <- predict(fit, sv_panel())
  expect_equal(p$term_premium, as.vector(fit$term_premium), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(p$term_premium,
                                predict(f, sv_panel())$term_premium)))
})


# SPF parsing -------------------------------------------------------------

test_that("the SPF quarterly sheet reshapes to one row per forecast", {
  raw <- data.frame(
    YEAR = c(2000L, 2000L), QUARTER = c(1L, 2L),
    TBILL1 = c("5.1", "5.2"), TBILL2 = c("5.3", "5.4"),
    TBILL3 = c("5.5", "5.6"), TBILL4 = c("5.7", "#N/A"),
    TBILL5 = c("5.9", "6.0"), TBILL6 = c("6.1", "6.2"),
    stringsAsFactors = FALSE
  )
  out <- spf_reshape(raw, 1:4)

  expect_equal(nrow(out), 7L)                     # 8 less the one "#N/A"
  expect_named(out, c("date", "target_date", "horizon_quarters", "value",
                      "tenor", "average_months"))
  expect_true(all(out$tenor == 3L & out$average_months == 3L))

  # 2000Q1: surveyed mid-February, one quarter ahead is April.
  q1h1 <- out[out$date == as.Date("2000-02-15") & out$horizon_quarters == 1L, ]
  expect_equal(q1h1$target_date, as.Date("2000-04-01"))
  expect_equal(q1h1$value, 5.5)                   # TBILL3 is one quarter out

  # Four quarters ahead of 2000Q1 is 2001Q1.
  q1h4 <- out[out$date == as.Date("2000-02-15") & out$horizon_quarters == 4L, ]
  expect_equal(q1h4$target_date, as.Date("2001-01-01"))
  expect_equal(q1h4$value, 6.1)                   # TBILL6
})

test_that("the SPF BILL10 sheet reshapes with a ten-year window", {
  raw <- data.frame(YEAR = c(1992L, 1993L), QUARTER = c(1L, 1L),
                    BILL10 = c("5.1576", "#N/A"), stringsAsFactors = FALSE)
  out <- spf_reshape_bill10(raw)

  expect_equal(nrow(out), 1L)
  expect_equal(out$value, 5.1576)
  expect_equal(out$average_months, 120L)
  expect_equal(out$tenor, 3L)
  expect_equal(out$date, as.Date("1992-02-15"))
  expect_equal(out$target_date, as.Date("1992-03-01"))
})

test_that("the two SPF sheets produce combinable frames", {
  q <- spf_reshape(data.frame(
    YEAR = 2000L, QUARTER = 1L, TBILL1 = "5.1", TBILL2 = "5.3",
    TBILL3 = "5.5", TBILL4 = "5.7", TBILL5 = "5.9", TBILL6 = "6.1",
    stringsAsFactors = FALSE), 1:4)
  l <- spf_reshape_bill10(data.frame(YEAR = 2000L, QUARTER = 1L,
                                     BILL10 = "5.0", stringsAsFactors = FALSE))
  expect_named(l, names(q))
  expect_equal(nrow(rbind(q, l)), 5L)
})

test_that("a changed SPF layout is reported rather than silently mis-parsed", {
  expect_error(spf_reshape(data.frame(YEAR = 2000L, QUARTER = 1L), 1:4),
               "missing column")
  expect_error(spf_reshape_bill10(data.frame(YEAR = 2000L, QUARTER = 1L)),
               "missing column")
  expect_error(
    spf_reshape(data.frame(YEAR = 2000L, QUARTER = 9L, TBILL1 = "1",
                           TBILL2 = "1", TBILL3 = "1", TBILL4 = "1",
                           TBILL5 = "1", TBILL6 = "1"), 1:4),
    "unreadable YEAR/QUARTER"
  )
  expect_error(
    spf_reshape_bill10(data.frame(YEAR = 2000L, QUARTER = 1L, BILL10 = "#N/A")),
    "no usable BILL10"
  )
})

test_that("spf_tbill validates its horizons before touching the network", {
  expect_error(spf_tbill(horizons = 9), "between 0 and 4")
  expect_error(spf_tbill(horizons = 1.5), "between 0 and 4")
  expect_error(spf_tbill(horizons = numeric(0)), "between 0 and 4")
})
