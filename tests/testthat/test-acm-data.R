# acm_published() downloads the New York Fed's ACM series rather than
# bundling it, so its parser is tested on a synthetic sheet laid out like the
# published one. The real file is exercised in test-data.R, off CRAN.

fake_acm_sheet <- function(dates = c("31-Jan-2000", "29-Feb-2000",
                                     "31-Mar-2000")) {
  n <- length(dates)
  sheet <- data.frame(DATE = dates, stringsAsFactors = FALSE)
  for (k in 1:10) {
    tp <- 0.1 * k + seq_len(n) / 100
    rn <- 2 + 0.05 * k + seq_len(n) / 50
    sheet[[sprintf("ACMY%02d", k)]] <- tp + rn
    sheet[[sprintf("ACMTP%02d", k)]] <- tp
    sheet[[sprintf("ACMRNY%02d", k)]] <- rn
  }
  sheet
}

test_that("the sheet is reshaped to one row per date and maturity", {
  out <- acm_reshape(fake_acm_sheet())

  expect_named(out,
               c("date", "maturity", "fitted", "term_premium", "risk_neutral"))
  expect_s3_class(out$date, "Date")
  expect_equal(nrow(out), 30L)
  expect_equal(unique(out$maturity), seq(12L, 120L, by = 12L))
  expect_false(is.unsorted(out$date))

  # The 5-year term premium on the second date: 0.1 * 5 + 2 / 100
  row <- out[out$date == as.Date("2000-02-29") & out$maturity == 60L, ]
  expect_equal(row$term_premium, 0.52)
  expect_equal(row$risk_neutral, 2 + 0.25 + 2 / 50)
  expect_equal(row$fitted, row$term_premium + row$risk_neutral)
})

test_that("rows come back in date order whatever the sheet's order", {
  sheet <- fake_acm_sheet()
  out <- acm_reshape(sheet[3:1, ])
  expect_equal(out, acm_reshape(sheet))
})

test_that("dates are read without depending on the locale", {
  # %b would follow LC_TIME and fail on a non-English system.
  expect_equal(parse_acm_date(c("30-Jun-1961", "1-JAN-2001", "31-dec-2024")),
               as.Date(c("1961-06-30", "2001-01-01", "2024-12-31")))
  expect_equal(parse_acm_date("2024-12-31"), as.Date("2024-12-31"))
  expect_equal(parse_acm_date(as.POSIXct("2024-12-31", tz = "UTC")),
               as.Date("2024-12-31"))
  expect_equal(parse_acm_date(c("31-Foo-2024", "yesterday")),
               as.Date(c(NA, NA)))
})

test_that("`end` truncates the sample, inclusively", {
  out <- acm_reshape(fake_acm_sheet(), end = "2000-02-29")
  expect_equal(unique(out$date), as.Date(c("2000-01-31", "2000-02-29")))

  expect_error(acm_reshape(fake_acm_sheet(), end = "1999-12-31"),
               "no observations")
  expect_error(acm_reshape(fake_acm_sheet(), end = c("2000-01-31",
                                                     "2000-02-29")),
               "single date")
})

test_that("a changed layout is an error, not a silent misread", {
  sheet <- fake_acm_sheet()

  expect_error(acm_reshape(sheet[setdiff(names(sheet), "ACMTP07")]),
               "missing column.*ACMTP07")

  bad_date <- sheet
  bad_date$DATE[2] <- "February 2000"
  expect_error(acm_reshape(bad_date), "Could not read 1 date")

  # Columns that no longer add up mean they are not what they are assumed to be
  swapped <- sheet
  swapped$ACMTP03 <- sheet$ACMRNY03
  expect_error(acm_reshape(swapped), "do not add up")
})
