test_that("gsw_monthly has the documented shape", {
  expect_s3_class(gsw_monthly, "data.frame")
  expect_named(gsw_monthly, c("date", "maturity", "yield", "extrapolated"))
  expect_s3_class(gsw_monthly$date, "Date")
  expect_type(gsw_monthly$extrapolated, "logical")

  expect_equal(sort(unique(gsw_monthly$maturity)), 1:120)
  expect_true(all(gsw_monthly$yield > -5 & gsw_monthly$yield < 25))
})

test_that("gsw_monthly is a complete date-by-maturity grid", {
  n_dates <- length(unique(gsw_monthly$date))
  expect_equal(nrow(gsw_monthly), n_dates * 120L)
  expect_false(any(duplicated(gsw_monthly[, c("date", "maturity")])))
  expect_false(anyNA(gsw_monthly$yield))
})

test_that("the extrapolation flag behaves monotonically in maturity", {
  # If a maturity is extrapolated on a given date, every longer one must be too
  by_date <- split(gsw_monthly, gsw_monthly$date)
  sampled <- by_date[seq(1, length(by_date), length.out = 25)]

  for (d in sampled) {
    d <- d[order(d$maturity), ]
    expect_false(is.unsorted(d$extrapolated))
  }
})

test_that("short maturities are never extrapolated", {
  short <- gsw_monthly[gsw_monthly$maturity <= 12, ]
  expect_false(any(short$extrapolated))
})

test_that("acm_published has the documented shape", {
  expect_s3_class(acm_published, "data.frame")
  expect_named(acm_published,
               c("date", "maturity", "fitted", "term_premium", "risk_neutral"))
  expect_s3_class(acm_published$date, "Date")
  expect_equal(sort(unique(acm_published$maturity)), seq(12L, 120L, by = 12L))
})

test_that("the published ACM decomposition is exact", {
  # fitted = risk_neutral + term_premium, to machine precision. This is the
  # identity our own implementation must also satisfy.
  gap <- with(acm_published, fitted - (risk_neutral + term_premium))
  expect_lt(max(abs(gap)), 1e-12)
})

test_that("gsw_monthly and acm_published share an identical date index", {
  # M4 compares our ACM fit on gsw_monthly against acm_published, so any
  # mismatch here would silently misalign the benchmark.
  g <- sort(unique(gsw_monthly$date))
  a <- sort(unique(acm_published$date))
  expect_identical(g, a)
})

test_that("bundled data feeds straight into yield_panel", {
  p <- yield_panel(gsw_monthly, units = "percent", maturity_unit = "months",
                   instrument = "government", issuer = "US")

  expect_s3_class(p, "yield_panel")
  expect_equal(p$frequency, "monthly")
  expect_equal(p$maturities[[1]], 1:120)
  expect_equal(p$meta$instrument, "government")

  m <- curve_matrix(p)
  expect_equal(dim(m), c(length(unique(gsw_monthly$date)), 120L))
  expect_false(anyNA(m))

  # converted to decimals on the way in
  expect_lt(max(m), 0.25)
})
