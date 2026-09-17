us_panel3 <- function() {
  yield_panel(gsw_monthly, units = "percent", maturity_unit = "months",
              instrument = "government", issuer = "US")
}

fit3 <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- suppressWarnings(atsm(us_panel3(), n_factors = 5))
    cache
  }
})


# Plotting ----------------------------------------------------------------

test_that("every plot type renders without error", {
  f <- fit3()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  expect_silent(plot(f))
  expect_silent(plot(f, type = "premium"))
  expect_silent(plot(f, type = "residuals"))
  expect_silent(plot(f, maturity = 24))
})

test_that("plot returns the fit invisibly and validates maturity", {
  f <- fit3()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  expect_identical(plot(f), f)
  expect_error(plot(f, maturity = 999), "was not fitted")
})

test_that("plot restores graphical parameters it changed", {
  # The residuals panel splits the device; leaving mfrow altered would corrupt
  # whatever the user plots next.
  f <- fit3()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  before <- par("mfrow")
  plot(f, type = "residuals")
  expect_equal(par("mfrow"), before)
})


# Lower-bound diagnostic --------------------------------------------------

test_that("zlb_check reports without warning for plausible rates", {
  rn <- matrix(seq(0.01, 0.05, length.out = 40), 10, 4)
  expect_silent(z <- zlb_check(rn, maturities = c(1, 12, 60, 120)))

  expect_equal(z$frac_negative, 0)
  expect_equal(z$min, min(rn))
})

test_that("zlb_check tolerates mildly negative rates", {
  # Polish, euro area and Swiss policy rates were genuinely below zero, so a
  # small negative value is data, not a modelling failure.
  rn <- matrix(c(rep(0.02, 30), rep(-0.004, 10)), 10, 4)
  expect_silent(z <- zlb_check(rn, maturities = c(1, 12, 60, 120)))
  expect_gt(z$frac_negative, 0)
  expect_equal(z$frac_below_threshold, 0)
})

test_that("zlb_check warns only on implausibly negative rates", {
  rn <- matrix(c(rep(0.02, 30), rep(-0.05, 10)), 10, 4)
  expect_warning(z <- zlb_check(rn, maturities = c(1, 12, 60, 120)),
                 "no lower bound")
  expect_gt(z$frac_below_threshold, 0)
  expect_equal(z$worst_maturity, 120)
})

test_that("a fitted model carries its lower-bound diagnostic", {
  f <- fit3()
  expect_true(all(c("min", "frac_negative", "worst_maturity") %in% names(f$zlb)))
  expect_equal(f$zlb$min, min(f$risk_neutral))
})

test_that("the US fit breaches the lower bound, and says so", {
  # Not a defect in the implementation but a documented limitation of Gaussian
  # affine models: fitted on 1961-2024 data spanning the zero lower bound, the
  # model projects expected short rates to about -3%, a policy rate the US has
  # never delivered. The BIS 2018 review flags exactly this. Asserted here so
  # the behaviour is a known, visible property rather than a surprise.
  expect_warning(atsm(us_panel3(), n_factors = 5), "no lower bound")

  f <- fit3()
  expect_lt(f$zlb$min, -0.01)
  expect_gt(f$zlb$frac_negative, 0)
})

test_that("print surfaces the lower-bound breach", {
  expect_output(print(fit3()), "negative")
})


# Survey-based term premium -----------------------------------------------

make_survey <- function(n = 20, start = "2000-02-15", value = 3.5) {
  data.frame(
    date = seq(as.Date(start), by = "year", length.out = n),
    value = rep(value, n)
  )
}

test_that("the survey premium is yield minus expected short rate", {
  p <- us_panel3()
  s <- make_survey()

  tp <- term_premium_survey(p, s, maturity = 120, survey_units = "percent")

  expect_named(tp, c("date", "survey_date", "maturity", "yield",
                     "expected_short_rate", "term_premium"))
  expect_equal(tp$term_premium, tp$yield - tp$expected_short_rate)
  expect_true(all(tp$expected_short_rate == 0.035))
  expect_equal(unique(tp$maturity), 120)
})

test_that("one row is returned per matched survey date, not per panel date", {
  # An annual survey must not silently become a monthly series.
  p <- us_panel3()
  s <- make_survey(n = 20)

  tp <- term_premium_survey(p, s, maturity = 120, survey_units = "percent")

  expect_lte(nrow(tp), 20L)
  expect_lt(nrow(tp), length(p$dates))
})

test_that("survey dates are matched to nearby panel observations", {
  p <- us_panel3()
  s <- make_survey(n = 5)

  tp <- term_premium_survey(p, s, maturity = 120, survey_units = "percent",
                            max_gap = 45)
  expect_true(all(abs(as.numeric(tp$date - tp$survey_date)) <= 45))
})

test_that("survey dates outside the panel are dropped and reported", {
  p <- us_panel3()
  s <- data.frame(date = as.Date(c("2010-02-15", "2100-02-15")), value = 3.5)

  expect_message(
    tp <- term_premium_survey(p, s, maturity = 120, survey_units = "percent"),
    "no panel observation within"
  )
  expect_equal(nrow(tp), 1L)
})

test_that("a survey with no overlap at all is an error", {
  p <- us_panel3()
  s <- data.frame(date = as.Date(c("2100-01-01", "2101-01-01")), value = 3.5)
  expect_error(
    suppressMessages(term_premium_survey(p, s, maturity = 120,
                                         survey_units = "percent")),
    "No survey date could be matched"
  )
})

test_that("interpolation is off by default and loud when requested", {
  p <- us_panel3()
  s <- make_survey(n = 20)

  plain <- term_premium_survey(p, s, maturity = 120, survey_units = "percent")

  expect_message(
    interp <- term_premium_survey(p, s, maturity = 120,
                                  survey_units = "percent",
                                  interpolate = TRUE),
    "constructed, not observed"
  )
  expect_gt(nrow(interp), nrow(plain))
  expect_true(all(is.na(interp$survey_date)))
})

test_that("term_premium_survey validates its inputs", {
  p <- us_panel3()
  s <- make_survey()

  expect_error(term_premium_survey("nope", s), "must be a yield_panel")
  expect_error(term_premium_survey(p, list(a = 1)), "must be a data frame")
  expect_error(term_premium_survey(p, data.frame(date = Sys.Date())),
               "missing column")
  expect_error(term_premium_survey(p, s, maturity = 999), "not on the panel")
})

test_that("survey and model premia are comparable in magnitude", {
  # Not a precise agreement test -- the point of the model-free measure is that
  # it is independent -- but they should be in the same ballpark.
  f <- fit3()
  p <- us_panel3()
  s <- make_survey(n = 25, start = "1995-02-15", value = 4.0)

  surv <- term_premium_survey(p, s, maturity = 120, survey_units = "percent")
  model <- term_premium(f, 120)
  m <- merge(surv[, c("date", "term_premium")], model, by = "date")

  expect_gt(nrow(m), 10L)
  expect_lt(abs(mean(m$term_premium) - mean(m$value)) * 1e4, 300)
})
