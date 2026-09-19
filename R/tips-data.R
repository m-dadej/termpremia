# Inputs for the joint real-nominal model that the package fetches rather
# than bundles.
#
# None of these are withheld for licensing reasons -- the Federal Reserve
# Board's TIPS curve and FRED's own series carry the same terms as the
# bundled nominal curve. They are fetched because the bundled snapshots are
# deliberately kept to the minimum needed to demonstrate and test the
# package, and because a real-nominal fit wants a current sample rather than
# a frozen one.

#' Download the Gurkaynak-Sack-Wright TIPS zero-coupon curve
#'
#' Fetches the Federal Reserve Board's fitted real yield curve (series
#' `feds200805`) and evaluates it onto a monthly maturity grid, which is the
#' form [atsm_real()] needs.
#'
#' @section Why the published columns cannot be used directly:
#' The same three traps as the nominal file, plus two that are specific to
#' this one.
#'
#' The published `TIPSY` columns run `TIPSY02` to `TIPSY20` in **annual**
#' steps from two years. A one-month holding return needs maturities `n` and
#' `n - 1` adjacent on the grid, so the curve has to be re-evaluated from the
#' Svensson parameters rather than read off.
#'
#' Rows switch between functional forms, and the encoding is **not** the
#' nominal file's `-999.99` sentinel: here a Nelson-Siegel row carries a
#' literal `NA` in `BETA3` and `TAU2`. All of 1999-2003 is Nelson-Siegel,
#' because too few TIPS were outstanding to identify the second hump, and
#' scattered later rows are too. A missing `TAU2` therefore does not mean a
#' missing curve, and code that treats it that way discards five years of
#' data.
#'
#' Separately, rows with **no** curve at all -- every parameter `NA` -- fall
#' on US federal holidays. Nine of them are the last calendar day of their
#' month, concentrated in March and May (Good Friday and Memorial Day), so
#' selecting month-ends before dropping them silently loses those months.
#' This is the same order-of-operations bug that cost twenty months on the
#' nominal curve, and the filtering here is deliberately in the safe order.
#'
#' @section Joining to a nominal curve:
#' Each file's month-end is whichever business day that month last carried a
#' fitted curve, and the two files do not always agree: over 1999-2024 the
#' nominal and TIPS curves share 311 of 312 calendar months exactly, but in
#' June 2003 the nominal file has a curve on the 30th and the TIPS file's last
#' one is the 27th. A [yield_panel] keys on exact dates, so that single
#' mismatch is enough to put a one-month hole in the joint panel -- which
#' [atsm_real()] then refuses to fit, correctly, because its factor VAR reads
#' consecutive rows as consecutive months.
#'
#' Pass `align_to` with the nominal curve's dates to avoid this. Each TIPS
#' observation is relabelled to the nominal date for the same calendar month,
#' which is what a monthly panel means in any case; months absent from
#' `align_to` are dropped, and the count is reported.
#'
#' @param maturities Maturity grid in months. The default starts at 23 rather
#'   than 24 so that a one-month return on a two-year TIPS can be formed;
#'   see [atsm_real()].
#' @param align_to Optional vector of dates, normally the nominal curve's
#'   month-ends, to relabel observations onto. See below.
#' @param end Optional last date to keep.
#' @param url Source location, exposed so a local copy can be used.
#' @param cache Optional path to a downloaded copy. When it exists it is read
#'   instead of downloading; when it does not, the download is saved there.
#' @param quiet Suppress the download progress bar.
#' @param check Verify the evaluated curve against the file's own published
#'   `TIPSY` columns and stop if they disagree by more than a tenth of a basis
#'   point. Verifying against tenors that are off the grid would compare
#'   nothing and pass, so off-grid check tenors are an error.
#'
#' @return A data frame with `date`, `maturity` (months), `yield` (percent)
#'   and `extrapolated`, flagging maturities beyond the longest tenor the file
#'   itself publishes on that date.
#'
#' @references
#' Gurkaynak, R. S., B. Sack and J. H. Wright (2010). "The TIPS yield curve
#' and inflation compensation." *American Economic Journal: Macroeconomics*
#' 2(1), 70-92.
#'
#' @seealso [atsm_real()], [fred_series()]
#'
#' @examples
#' \dontrun{
#' tips <- gsw_tips(cache = "feds200805.csv")
#' str(tips)
#' }
#'
#' @export
gsw_tips <- function(maturities = 23:120,
                     align_to = NULL,
                     end = NULL,
                     url = paste0("https://www.federalreserve.gov/data/",
                                  "yield-curve-tables/feds200805.csv"),
                     cache = NULL,
                     quiet = TRUE,
                     check = TRUE) {
  if (!is.numeric(maturities) || !length(maturities) ||
      any(maturities < 1 | maturities != as.integer(maturities))) {
    stop("`maturities` must be whole numbers of months, at least 1.",
         call. = FALSE)
  }
  maturities <- as.integer(maturities)

  path <- fetch_to_file(url, cache, quiet,
                        what = "the GSW TIPS curve (feds200805)")
  raw <- read_fed_curve_file(path)

  needed <- c("BETA0", "BETA1", "BETA2", "BETA3", "TAU1", "TAU2")
  missing_cols <- setdiff(needed, names(raw))
  if (length(missing_cols)) {
    stop("The TIPS file has no column(s) ",
         paste(missing_cols, collapse = ", "),
         ". The Board have changed the layout before; check the file.",
         call. = FALSE)
  }

  raw <- raw[order(raw$Date), , drop = FALSE]
  if (!is.null(end)) {
    end <- as.Date(end)
    raw <- raw[raw$Date <= end, , drop = FALSE]
  }

  # ORDER MATTERS. Drop the no-curve rows first; only then pick month-ends.
  # BETA3 and TAU2 are deliberately not required: their absence marks a
  # Nelson-Siegel fit, which is a curve, not a gap.
  has_curve <- stats::complete.cases(raw[, c("BETA0", "BETA1", "BETA2",
                                             "TAU1")])
  if (!any(has_curve)) {
    stop("No row of the TIPS file has a fitted curve.", call. = FALSE)
  }
  raw <- raw[has_curve, , drop = FALSE]

  ym <- format(raw$Date, "%Y-%m")
  eom <- raw[sort(unname(tapply(seq_len(nrow(raw)), ym, max))), , drop = FALSE]

  if (!is.null(align_to)) {
    align_to <- sort(unique(as.Date(align_to)))
    target <- align_to[match(format(eom$Date, "%Y-%m"),
                             format(align_to, "%Y-%m"))]
    dropped <- sum(is.na(target))
    if (dropped) {
      message("gsw_tips(): dropped ", dropped, " month(s) absent from ",
              "`align_to`.")
    }
    eom <- eom[!is.na(target), , drop = FALSE]
    if (!nrow(eom)) {
      stop("`align_to` shares no calendar month with the TIPS file.",
           call. = FALSE)
    }
    eom$Date <- target[!is.na(target)]
  }

  params <- data.frame(
    beta0 = eom$BETA0, beta1 = eom$BETA1, beta2 = eom$BETA2,
    beta3 = eom$BETA3, tau1 = eom$TAU1, tau2 = eom$TAU2
  )
  yields_pct <- svensson_curve(params, maturity = maturities / 12)

  tipsy <- sprintf("TIPSY%02d", 2:20)
  tipsy <- tipsy[tipsy %in% names(eom)]
  published <- as.matrix(eom[, tipsy, drop = FALSE])
  longest_m <- 12 * apply(published, 1L, function(r) {
    ok <- which(is.finite(r))
    if (!length(ok)) NA_integer_ else as.integer(sub("TIPSY", "", tipsy[max(ok)]))
  })

  if (check) {
    check_tips_curve(yields_pct, eom, maturities)
  }

  out <- data.frame(
    date = rep(eom$Date, times = length(maturities)),
    maturity = rep(maturities, each = nrow(eom)),
    yield = as.vector(yields_pct),
    extrapolated = rep(maturities, each = nrow(eom)) >
      rep(longest_m, times = length(maturities)),
    stringsAsFactors = FALSE
  )
  out <- out[order(out$date, out$maturity), , drop = FALSE]
  rownames(out) <- NULL

  attr(out, "source_url") <- url
  attr(out, "retrieved") <- Sys.Date()
  out
}

#' Verify the evaluated TIPS curve against the file's own published tenors
#'
#' Copied in spirit from `data-raw/gsw.R`: a check tenor that is not on the
#' grid, or that has no overlapping finite observations, is an error rather
#' than a silent pass. An earlier version of the nominal builder "validated"
#' 20- and 30-year tenors against a 120-month grid, comparing nothing.
#'
#' @keywords internal
#' @noRd
check_tips_curve <- function(yields_pct, eom, maturities) {
  # A tenor is only checkable if it is BOTH on the requested grid and
  # published in the file. Testing `yr * 12 <= max(maturities)` is not enough
  # -- a grid of 25:30 months passes that for two years while having no
  # 24-month column to compare -- and skipping such tenors one by one lets
  # the whole check pass having compared nothing, which is the failure mode
  # this function exists to prevent.
  check_years <- c(2, 5, 10)
  checkable <- check_years[
    (check_years * 12L) %in% maturities &
      sprintf("TIPSY%02d", check_years) %in% names(eom)
  ]
  if (!length(checkable)) {
    stop("`check = TRUE` but none of the 2, 5 or 10 year tenors is both on ",
         "the requested `maturities` grid and published in the file, so the ",
         "check could only pass vacuously. Include 24, 60 or 120 months in ",
         "`maturities`, or pass `check = FALSE`.", call. = FALSE)
  }

  compared <- 0L
  for (yr in checkable) {
    col <- match(yr * 12L, maturities)
    nm <- sprintf("TIPSY%02d", yr)

    mine <- yields_pct[, col]
    theirs <- eom[[nm]]
    ok <- is.finite(mine) & is.finite(theirs)
    if (!sum(ok)) next

    worst <- max(abs(mine[ok] - theirs[ok]))
    if (worst > 1e-3) {
      stop("Evaluated ", yr, "y TIPS yields differ from the published ", nm,
           " by up to ", format(worst, digits = 3), " percentage points ",
           "across ", sum(ok), " observations. The evaluator and the file ",
           "disagree; do not use this data.", call. = FALSE)
    }
    compared <- compared + 1L
  }

  if (!compared) {
    stop("No checkable tenor had a single overlapping finite observation, so ",
         "the curve was never actually verified. Refusing to pass vacuously.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Read a Federal Reserve Board fitted-curve CSV
#'
#' The preamble length differs between files -- nine lines before the header
#' in the nominal file, eighteen in the TIPS one -- and has changed over time,
#' so the header is located by looking for it rather than assumed.
#'
#' @keywords internal
#' @noRd
read_fed_curve_file <- function(path) {
  head_lines <- readLines(path, n = 200L, warn = FALSE)
  skip <- grep("^\"?Date\"?,", head_lines)
  if (!length(skip)) {
    stop("Could not find the header row (a line beginning \"Date,\") in the ",
         "first 200 lines of\n  ", path,
         "\nThe file layout may have changed.", call. = FALSE)
  }

  raw <- utils::read.csv(path, skip = skip[1L] - 1L, stringsAsFactors = FALSE,
                         na.strings = c("NA", "", "#N/A"))
  raw$Date <- as.Date(raw$Date)
  if (anyNA(raw$Date)) {
    raw <- raw[!is.na(raw$Date), , drop = FALSE]
  }
  raw
}

#' @keywords internal
#' @noRd
fetch_to_file <- function(url, cache, quiet, what) {
  if (!is.null(cache) && file.exists(cache)) return(cache)

  dest <- cache %||% tempfile(fileext = ".csv")
  tryCatch(
    utils::download.file(url, dest, mode = "wb", quiet = quiet),
    error = function(e) {
      stop("Could not download ", what, " from\n  ", url, "\n",
           conditionMessage(e),
           "\nPass a local copy through `cache` if you already have the file.",
           call. = FALSE)
    }
  )
  dest
}

#' Download a series from FRED
#'
#' A small convenience for the two macro inputs [atsm_real()] needs that are
#' not yield curves: a price index and a short rate.
#'
#' @section Which series:
#' The paper uses seasonally **unadjusted** CPI-U as the price index the TIPS
#' payouts are tied to (`"CPIAUCNS"`), and the effective federal funds rate as
#' the nominal short rate (`"DFF"` daily, or `"FEDFUNDS"` monthly).
#'
#' Pass `at_dates` for a short rate. `r_t` in the model is the riskless return
#' earned from `t` to `t + 1`, so what is wanted is the rate *prevailing on*
#' each panel date, not a backward-looking average of the month that just
#' ended. With `at_dates` the most recent observation at or before each date
#' is returned, which also sidesteps the fact that a month-end panel date can
#' be a weekend with no FRED observation on it.
#'
#' A price index needs no such alignment: [atsm_real()] matches it by calendar
#' month, because an index dated to the first of the month and a curve dated
#' to the last business day would never join exactly.
#'
#' @param series_id A FRED series identifier, e.g. `"CPIAUCNS"`.
#' @param at_dates Optional vector of dates. When supplied, one row is
#'   returned per date, carrying the last observation at or before it forward.
#' @param url Template for the download, with `%s` for the series id.
#' @param quiet Suppress the download progress bar.
#'
#' @return A data frame with `date` and `value`.
#'
#' @section Requirements:
#' A working internet connection. FRED is not a dependency of the package and
#' nothing here is redistributed.
#'
#' @seealso [atsm_real()], [gsw_tips()]
#'
#' @examples
#' \dontrun{
#' cpi <- fred_series("CPIAUCNS")
#' ff <- fred_series("DFF", at_dates = unique(gsw_monthly$date))
#' }
#'
#' @export
fred_series <- function(series_id,
                        at_dates = NULL,
                        url = "https://fred.stlouisfed.org/graph/fredgraph.csv?id=%s",
                        quiet = TRUE) {
  if (!is.character(series_id) || length(series_id) != 1L ||
      !nzchar(series_id)) {
    stop("`series_id` must be a single non-empty string.", call. = FALSE)
  }

  path <- fetch_to_file(sprintf(url, series_id), NULL, quiet,
                        what = paste0("FRED series ", series_id))
  raw <- utils::read.csv(path, stringsAsFactors = FALSE,
                         na.strings = c("NA", ".", ""))

  if (ncol(raw) < 2L) {
    stop("FRED returned an unexpected layout for ", series_id, " (",
         ncol(raw), " column(s)). Check that the series id is correct.",
         call. = FALSE)
  }

  out <- data.frame(
    date = as.Date(raw[[1L]]),
    value = suppressWarnings(as.numeric(raw[[2L]])),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$date), , drop = FALSE]
  out <- out[order(out$date), , drop = FALSE]

  if (all(is.na(out$value))) {
    stop("Every value of FRED series ", series_id, " parsed as NA.",
         call. = FALSE)
  }
  if (is.null(at_dates)) {
    rownames(out) <- NULL
    return(out)
  }

  at_dates <- as.Date(at_dates)
  have <- out[!is.na(out$value), , drop = FALSE]
  idx <- findInterval(at_dates, have$date)
  if (any(idx == 0L)) {
    stop(series_id, " begins ", format(min(have$date)), ", after ",
         sum(idx == 0L), " of the requested dates (first: ",
         format(min(at_dates)), "). Restrict `at_dates` to the period the ",
         "series covers.", call. = FALSE)
  }

  data.frame(date = at_dates, value = have$value[idx],
             stringsAsFactors = FALSE)
}

#' Build the TIPS liquidity factor
#'
#' Combines observable indicators of how illiquid inflation-indexed bonds are
#' relative to nominal ones into the single pricing factor that
#' [atsm_real()] takes as `liquidity`, following Section 3.1 of the paper:
#' standardise each indicator, average them with equal weights, then shift the
#' result so its minimum is zero.
#'
#' @section Why positivity:
#' The shift is not cosmetic. With the factor constrained to be
#' non-negative, the liquidity term can only ever raise a yield, never lower
#' it, which is the only economically sensible sign for a discount demanded
#' for illiquidity.
#'
#' @section The paper's two indicators:
#' The first is the average absolute fitting error of the
#' Nelson-Siegel-Svensson curve fitted to TIPS: large errors mean the market
#' is not arbitraging small mispricings, which is what illiquidity looks like
#' from the outside. **This series is not part of the published `feds200805`
#' file** and has to be obtained from the Board directly, which is why
#' `atsm_real()` treats `liquidity` as optional.
#'
#' The second is a 13-week moving average of the ratio of primary dealers'
#' nominal Treasury to TIPS transaction volumes, from the New York Fed's
#' FR2004 release. That one is public.
#'
#' Supplying only the volume ratio is a defensible fallback and the function
#' allows it, but it is a different factor from the paper's and should be
#' described as such.
#'
#' @param ... One or more data frames with `date` and `value` columns, or
#'   numeric vectors of equal length. Each is one indicator.
#' @param dates Dates corresponding to the rows, required when the indicators
#'   are given as bare numeric vectors.
#'
#' @return A data frame with `date` and `value`, ready to pass as
#'   `liquidity`.
#'
#' @seealso [atsm_real()]
#'
#' @examples
#' # A single synthetic indicator, standardised and shifted to be positive.
#' d <- seq(as.Date("2000-01-31"), by = "month", length.out = 24)
#' ind <- data.frame(date = d, value = c(rep(1, 12), rep(3, 12)))
#' liq <- tips_liquidity_factor(ind)
#' range(liq$value)
#'
#' @export
tips_liquidity_factor <- function(..., dates = NULL) {
  parts <- list(...)
  if (!length(parts)) {
    stop("Supply at least one liquidity indicator.", call. = FALSE)
  }

  frames <- lapply(seq_along(parts), function(i) {
    p <- parts[[i]]
    if (is.data.frame(p)) {
      if (!all(c("date", "value") %in% names(p))) {
        stop("Indicator ", i, " must have `date` and `value` columns.",
             call. = FALSE)
      }
      data.frame(date = as_date_strict(p$date, arg = "date"),
                 value = as.numeric(p$value), stringsAsFactors = FALSE)
    } else if (is.numeric(p)) {
      if (is.null(dates)) {
        stop("Indicator ", i, " is a bare numeric vector, so `dates` is ",
             "required.", call. = FALSE)
      }
      if (length(p) != length(dates)) {
        stop("Indicator ", i, " has ", length(p), " values but `dates` has ",
             length(dates), ".", call. = FALSE)
      }
      data.frame(date = as.Date(dates), value = as.numeric(p),
                 stringsAsFactors = FALSE)
    } else {
      stop("Indicator ", i, " must be a data frame or a numeric vector.",
           call. = FALSE)
    }
  })

  all_dates <- sort(unique(do.call(c, lapply(frames, `[[`, "date"))))

  # Standardise each indicator on its own before averaging: they are measured
  # in unrelated units (basis points of fitting error against a volume ratio),
  # so an unstandardised average would be whichever one happens to be larger.
  z <- vapply(frames, function(f) {
    v <- f$value[match(all_dates, f$date)]
    s <- stats::sd(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) {
      stop("An indicator has zero or undefined variance, so it cannot be ",
           "standardised.", call. = FALSE)
    }
    (v - mean(v, na.rm = TRUE)) / s
  }, numeric(length(all_dates)))

  if (length(frames) == 1L) z <- matrix(z, ncol = 1L)
  composite <- rowMeans(z, na.rm = TRUE)

  usable <- is.finite(composite)
  if (!any(usable)) {
    stop("No date has any indicator observed.", call. = FALSE)
  }
  composite[!usable] <- NA_real_

  data.frame(
    date = all_dates,
    value = composite - min(composite, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
