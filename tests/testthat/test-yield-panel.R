make_long <- function(n_dates = 12, mats = c(12, 60, 120), curve = NULL,
                      pct = TRUE) {
  d <- expand.grid(
    date = seq(as.Date("2020-01-31"), by = "month", length.out = n_dates),
    maturity = mats,
    KEEP.OUT.ATTRS = FALSE
  )
  d$yield <- 2 + d$maturity / 120
  if (!pct) d$yield <- d$yield / 100
  if (!is.null(curve)) d$curve <- curve
  d
}

test_that("a long data frame round-trips into canonical form", {
  p <- yield_panel(make_long(), units = "percent", maturity_unit = "months")

  expect_s3_class(p, "yield_panel")
  expect_length(p$dates, 12)
  expect_equal(p$maturities[[1]], c(12, 60, 120))

  m <- curve_matrix(p)
  expect_equal(dim(m), c(12L, 3L))

  # 2% + maturity/120 in percent becomes decimals
  expect_equal(unname(m[1, 1]), (2 + 12 / 120) / 100)
})

test_that("maturities in years are converted to months", {
  d <- make_long(mats = c(1, 5, 10))
  p <- yield_panel(d, units = "percent", maturity_unit = "years")
  expect_equal(p$maturities[[1]], c(12, 60, 120))
})

test_that("percent and decimal inputs give identical panels", {
  a <- yield_panel(make_long(pct = TRUE), units = "percent",
                   maturity_unit = "months")
  b <- yield_panel(make_long(pct = FALSE), units = "decimal",
                   maturity_unit = "months")
  expect_equal(curve_matrix(a), curve_matrix(b))
})

test_that("unit inference reports what it inferred", {
  expect_message(
    yield_panel(make_long(pct = TRUE), maturity_unit = "months"),
    "percent"
  )
  expect_message(
    yield_panel(make_long(pct = FALSE), maturity_unit = "months"),
    "decimal"
  )
})

test_that("implausible yields after conversion raise a warning", {
  # Decimal-looking data declared as decimals but actually in percent would be
  # silently wrong; declaring percent data as decimals is catchable.
  d <- make_long(pct = TRUE)
  d$yield <- d$yield * 100 # now ~200-300, decimal would mean 20000%
  expect_warning(
    yield_panel(d, units = "decimal", maturity_unit = "months"),
    "implausible"
  )
})

test_that("a wide matrix is accepted via dimnames", {
  m <- matrix(c(2.1, 2.2, 2.5, 2.6, 3.0, 3.1), nrow = 2,
              dimnames = list(c("2020-01-31", "2020-02-29"), c("12", "60", "120")))
  p <- yield_panel(m, units = "percent", maturity_unit = "months")

  expect_equal(p$maturities[[1]], c(12, 60, 120))
  expect_equal(unname(curve_matrix(p)[1, 1]), 0.021)
})

test_that("a matrix without dimnames or explicit labels is rejected", {
  m <- matrix(1:6, nrow = 2)
  expect_error(yield_panel(m, units = "percent", maturity_unit = "months"),
               "rownames")
  expect_error(
    yield_panel(m, dates = as.Date(c("2020-01-31", "2020-02-29")),
                units = "percent", maturity_unit = "months"),
    "colnames"
  )
})

test_that("mismatched dates or maturities are rejected", {
  m <- matrix(1:6, nrow = 2)
  expect_error(
    yield_panel(m, dates = as.Date("2020-01-31"), maturities = c(12, 60, 120),
                units = "percent", maturity_unit = "months"),
    "length 1 but"
  )
})

test_that("multiple curves share a date index but not a maturity grid", {
  a <- make_long(mats = c(12, 60, 120), curve = "nominal")
  b <- make_long(mats = c(60, 120), curve = "tips")
  p <- yield_panel(rbind(a, b), curve = "curve", units = "percent",
                   maturity_unit = "months",
                   instrument = c(nominal = "government", tips = "tips"),
                   issuer = "US")

  expect_equal(nrow(p$meta), 2L)
  expect_equal(p$maturities$nominal, c(12, 60, 120))
  expect_equal(p$maturities$tips, c(60, 120))
  expect_equal(p$meta$instrument[p$meta$curve == "tips"], "tips")
  expect_equal(dim(curve_matrix(p, "tips")), c(12L, 2L))
})

test_that("ragged date coverage across curves becomes explicit NA", {
  a <- make_long(n_dates = 12, curve = "nominal")
  b <- make_long(n_dates = 6, curve = "tips")
  p <- yield_panel(rbind(a, b), curve = "curve", units = "percent",
                   maturity_unit = "months")

  expect_length(p$dates, 12)
  expect_true(all(is.na(curve_matrix(p, "tips")[7:12, ])))
  expect_false(anyNA(curve_matrix(p, "nominal")))
})

test_that("duplicated observations are rejected", {
  d <- make_long()
  expect_error(
    yield_panel(rbind(d, d[1, ]), units = "percent", maturity_unit = "months"),
    "duplicated"
  )
})

test_that("invalid instruments and unmatched metadata are rejected", {
  d <- make_long(curve = "nominal")

  expect_error(
    yield_panel(d, curve = "curve", units = "percent",
                maturity_unit = "months", instrument = "equity"),
    "Unknown `instrument`"
  )
  expect_error(
    yield_panel(d, curve = "curve", units = "percent",
                maturity_unit = "months",
                instrument = c(other_curve = "government")),
    "no entry for curve"
  )
})

test_that("bad dates and maturities are rejected with useful messages", {
  d <- make_long()

  d_num <- d
  d_num$date <- as.numeric(d_num$date)
  expect_error(yield_panel(d_num, units = "percent", maturity_unit = "months"),
               "ambiguous")

  d_zero <- d
  d_zero$maturity[1] <- 0
  expect_error(yield_panel(d_zero, units = "percent", maturity_unit = "months"),
               "strictly positive")

  expect_error(
    yield_panel(d, date = "nope", units = "percent", maturity_unit = "months"),
    "missing column"
  )
})

test_that("frequency is detected from the date spacing", {
  monthly <- yield_panel(make_long(n_dates = 24), units = "percent",
                         maturity_unit = "months")
  expect_equal(monthly$frequency, "monthly")

  d <- expand.grid(date = seq(as.Date("2020-01-01"), by = "day", length.out = 30),
                   maturity = 12)
  d$yield <- 2
  expect_equal(
    yield_panel(d, units = "percent", maturity_unit = "months")$frequency,
    "daily"
  )
})

test_that("print and summary run and report the canonical unit", {
  p <- yield_panel(make_long(), units = "percent", maturity_unit = "months")

  expect_output(print(p), "yield_panel")
  expect_output(print(p), "decimal")
  expect_s3_class(summary(p), "summary.yield_panel")
  expect_output(print(summary(p)), "dates")
})

test_that("curve_matrix rejects an unknown curve name", {
  p <- yield_panel(make_long(), units = "percent", maturity_unit = "months")
  expect_error(curve_matrix(p, "nonexistent"), "No curve named")
})
