# The UK data layer: Bank of England curves and the ONS retail price index.
#
# Nothing here touches the network. The Bank publishes .xlsx inside .zip, and
# writing an .xlsx back out would need a package that is not a dependency, so
# the workbook reader is exercised through the internal functions that take
# the parsed structures -- which is where the logic worth testing lives
# anyway. `analysis/validate-uk.R` does the real-data half.


# Merging the two published sheets ---------------------------------------

test_that("the short-end sheet wins on maturities both sheets publish", {
  # The Bank publishes a monthly grid to five years on one sheet and a
  # half-yearly grid to twenty-five years on another. They overlap, and the
  # finer one has to take precedence or the monthly points are thrown away.
  dates <- as.Date(c("2000-01-31", "2000-02-29"))
  short <- list(dates = dates, maturity = c(1L, 2L, 3L),
                yield = rbind(c(5.0, 5.1, 5.2), c(4.0, 4.1, 4.2)))
  long <- list(dates = dates, maturity = c(3L, 6L),
               yield = rbind(c(9.9, 5.5), c(9.9, 4.5)))

  out <- boe_merge_sheets(list(short, long))

  expect_equal(out$maturity, c(1L, 2L, 3L, 6L))
  expect_equal(out$yield[1L, ], c(5.0, 5.1, 5.2, 5.5))
  expect_equal(out$yield[2L, ], c(4.0, 4.1, 4.2, 4.5))
})

test_that("a long-end value fills in where the short end has none", {
  dates <- as.Date("2000-01-31")
  short <- list(dates = dates, maturity = c(1L, 2L),
                yield = matrix(c(5.0, NA), 1L, 2L))
  long <- list(dates = dates, maturity = c(2L, 6L),
               yield = matrix(c(4.8, 5.5), 1L, 2L))

  out <- boe_merge_sheets(list(short, long))
  expect_equal(out$yield[1L, ], c(5.0, 4.8, 5.5))
})

test_that("sheets that disagree on their dates are refused", {
  a <- list(dates = as.Date("2000-01-31"), maturity = 1L,
            yield = matrix(5, 1L, 1L))
  b <- list(dates = as.Date("2000-02-29"), maturity = 1L,
            yield = matrix(5, 1L, 1L))
  expect_error(boe_merge_sheets(list(a, b)), "disagree on their dates")
})

test_that("workbooks covering different eras stack in date order", {
  # The archive splits the history across three files, and they arrive in
  # whatever order the zip lists them.
  later <- list(dates = as.Date(c("2016-01-31", "2016-02-29")),
                maturity = c(1L, 2L),
                yield = rbind(c(1.0, 1.1), c(1.2, 1.3)))
  earlier <- list(dates = as.Date(c("2015-11-30", "2015-12-31")),
                  maturity = c(1L, 2L, 3L),
                  yield = rbind(c(2.0, 2.1, 2.2), c(2.3, 2.4, 2.5)))

  out <- boe_bind(list(later, earlier))
  expect_equal(out$dates, as.Date(c("2015-11-30", "2015-12-31",
                                    "2016-01-31", "2016-02-29")))
  expect_equal(out$maturity, c(1L, 2L, 3L))
  expect_equal(out$yield[, 1L], c(2.0, 2.3, 1.0, 1.2))
  expect_true(is.na(out$yield[3L, 3L]))    # later file has no 3-month point
})

test_that("a date appearing in two workbooks is kept once", {
  a <- list(dates = as.Date(c("2015-12-31", "2016-01-31")), maturity = 1L,
            yield = matrix(c(1, 2), 2L, 1L))
  b <- list(dates = as.Date("2016-01-31"), maturity = 1L,
            yield = matrix(9, 1L, 1L))
  out <- boe_bind(list(a, b))
  expect_equal(out$dates, as.Date(c("2015-12-31", "2016-01-31")))
  expect_equal(nrow(out$yield), 2L)
})


# Interpolation ----------------------------------------------------------

test_that("interpolation reproduces the points it passes through", {
  curve <- list(dates = as.Date("2000-01-31"),
                maturity = c(12L, 24L, 36L, 48L, 60L),
                yield = matrix(c(5.0, 5.2, 5.35, 5.45, 5.5), 1L, 5L))

  out <- boe_interpolate(curve, c(12L, 36L, 60L))
  expect_equal(drop(out$yield), c(5.0, 5.35, 5.5), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_false(any(out$extrapolated))
})

test_that("maturities outside the published range are flagged, not hidden", {
  curve <- list(dates = as.Date("2000-01-31"),
                maturity = c(12L, 24L, 36L, 48L),
                yield = matrix(c(5.0, 5.2, 5.35, 5.45), 1L, 4L))

  out <- boe_interpolate(curve, c(1L, 12L, 48L, 120L))
  expect_equal(drop(out$extrapolated), c(TRUE, FALSE, FALSE, TRUE),
               ignore_attr = TRUE)
  expect_true(all(is.finite(out$yield)))
})

test_that("extrapolation is linear so it cannot run away", {
  # Natural rather than the default spline. A third of the UK nominal sample
  # needs the short end extended, and a cubic tail can leave the plausible
  # range entirely over a few months.
  curve <- list(dates = as.Date("2000-01-31"),
                maturity = c(12L, 24L, 36L, 48L, 60L),
                yield = matrix(c(5.0, 5.2, 5.35, 5.45, 5.5), 1L, 5L))

  # Three equally spaced points, all BELOW the shortest published maturity,
  # so the whole span is extrapolated: equal steps in maturity must give
  # equal steps in yield.
  out <- drop(boe_interpolate(curve, c(2L, 6L, 10L))$yield)
  expect_equal(unname(out[2L] - out[1L]), unname(out[3L] - out[2L]),
               tolerance = 1e-8)
  expect_true(all(out > 4 & out < 6))
})

test_that("a date with too few published points yields NA, not nonsense", {
  curve <- list(dates = as.Date(c("2000-01-31", "2000-02-29")),
                maturity = c(12L, 24L, 36L, 48L),
                yield = rbind(c(5.0, 5.2, 5.35, 5.45),
                              c(5.0, NA, NA, NA)))

  out <- boe_interpolate(curve, c(12L, 24L))
  expect_true(all(is.finite(out$yield[1L, ])))
  expect_true(all(is.na(out$yield[2L, ])))
})

test_that("the interpolation is checked against a held-out maturity", {
  # Comparing interpolated values with the points they were interpolated
  # through would compare a number with itself. Dropping a column first is
  # the only version of this check with any content.
  mats <- c(12L, 24L, 36L, 48L, 60L, 72L, 84L, 96L, 120L)
  smooth <- function(n) 5 + 0.5 * log(n / 12)
  curve <- list(dates = as.Date(c("2000-01-31", "2000-02-29")),
                maturity = mats,
                yield = rbind(smooth(mats), smooth(mats) + 0.1))

  expect_error(boe_check_interpolation(curve), "Only 2 date")

  many <- list(dates = seq(as.Date("2000-01-31"), by = "month",
                           length.out = 24L),
               maturity = mats,
               yield = matrix(smooth(mats), 24L, length(mats), byrow = TRUE))
  expect_lt(boe_check_interpolation(many), 0.01)

  # Corrupt the held-out column and the check must refuse the data.
  bad <- many
  bad$yield[, match(84L, mats)] <- 99
  expect_error(boe_check_interpolation(bad), "do not use this data")
})

test_that("a hold-out check with no column to hold out is refused", {
  curve <- list(dates = seq(as.Date("2000-01-31"), by = "month",
                            length.out = 24L),
                maturity = c(12L, 24L, 36L, 120L),
                yield = matrix(5, 24L, 4L))
  expect_error(boe_check_interpolation(curve), "would go unverified")
})


# Chaining the retail price index ----------------------------------------

test_that("the index is extended backwards through the percentage change", {
  # Build a known index, derive its exact twelve-month change, throw the early
  # years away, and check they come back.
  dates <- seq(as.Date("1984-01-01"), by = "month", length.out = 60L)
  level <- 100 * cumprod(c(1, rep(1.003, 59L)))
  full <- data.frame(date = dates, value = level, stringsAsFactors = FALSE)

  pct <- data.frame(
    date = dates[13:60],
    value = 100 * (level[13:60] / level[1:48] - 1),
    stringsAsFactors = FALSE
  )
  published <- full[full$date >= as.Date("1987-01-01"), , drop = FALSE]

  out <- ons_chain_back(published, pct, "1984-01-01")

  expect_equal(min(out$date), as.Date("1984-01-01"))
  expect_equal(nrow(out), 60L)
  expect_equal(out$value, level, tolerance = 1e-10)
  expect_equal(sum(out$chained), 36L)
  expect_false(any(out$chained[out$date >= as.Date("1987-01-01")]))
})

test_that("chaining stops and warns when the change series runs out", {
  dates <- seq(as.Date("1990-01-01"), by = "month", length.out = 36L)
  published <- data.frame(date = dates, value = 100 + seq_along(dates),
                          stringsAsFactors = FALSE)
  pct <- data.frame(date = dates[13:36], value = rep(3, 24L),
                    stringsAsFactors = FALSE)

  expect_warning(out <- ons_chain_back(published, pct, "1960-01-01"),
                 "could only be extended back")
  expect_gt(min(out$date), as.Date("1960-01-01"))
})

test_that("the index and its percentage change are checked for consistency", {
  dates <- seq(as.Date("1984-01-01"), by = "month", length.out = 120L)
  level <- 100 * cumprod(c(1, rep(1.004, 119L)))
  index <- data.frame(date = dates, value = level, stringsAsFactors = FALSE)
  pct <- data.frame(date = dates[13:120],
                    value = 100 * (level[13:120] / level[1:108] - 1),
                    stringsAsFactors = FALSE)

  expect_lt(ons_check_chaining(index, pct), 1e-10)

  # A factor-of-100 slip is exactly what this guards against.
  expect_error(ons_check_chaining(index, transform(pct, value = value / 100)),
               "cannot be chained backwards")
})

test_that("the consistency check is relative, so rounding does not fail it", {
  # The percentage change is published to one decimal place. On an index above
  # 400 that rounding is worth 0.2 index points, which an absolute tolerance
  # tight enough to be useful would reject.
  dates <- seq(as.Date("1990-01-01"), by = "month", length.out = 120L)
  level <- 400 * cumprod(c(1, rep(1.003, 119L)))
  index <- data.frame(date = dates, value = level, stringsAsFactors = FALSE)
  exact <- 100 * (level[13:120] / level[1:108] - 1)
  pct <- data.frame(date = dates[13:120], value = round(exact, 1L),
                    stringsAsFactors = FALSE)

  worst_abs <- max(abs(level[1:108] * (1 + pct$value / 100) - level[13:120]))
  expect_gt(worst_abs, 0.15)            # would fail an absolute 0.15 tolerance
  expect_silent(ons_check_chaining(index, pct))
})

test_that("a check with too little overlap is refused", {
  dates <- seq(as.Date("1990-01-01"), by = "month", length.out = 20L)
  index <- data.frame(date = dates, value = 100 + seq_along(dates),
                      stringsAsFactors = FALSE)
  pct <- data.frame(date = dates[13:20], value = rep(3, 8L),
                    stringsAsFactors = FALSE)
  expect_error(ons_check_chaining(index, pct), "refusing to pass vacuously")
})
