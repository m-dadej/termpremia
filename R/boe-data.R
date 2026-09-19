# UK inputs for the joint real-nominal model: the Bank of England's fitted
# government liability curves and the ONS retail price index.
#
# Fetched, never bundled, for the same reason as the TIPS curve: bundled data
# is kept to the minimum needed to test and demonstrate the package.

#' Download a Bank of England fitted government yield curve
#'
#' Fetches the Bank of England's month-end government liability curve --
#' nominal or inflation-indexed -- and puts it on a monthly maturity grid,
#' which is the form [atsm_real()] needs. This is the curve the UK
#' specification in the Supplementary Appendix of Abrahams et al. (2016) is
#' estimated on.
#'
#' @section What the published files contain:
#' Two sheets matter and both are needed. `"3. spot, short end"` carries a
#' genuinely **monthly** grid -- 1 to 60 months for the nominal curve, 25 to
#' 60 for the real one -- and `"4. spot curve"` carries **half-yearly** steps
#' out to 25 years. Merging them gives monthly coverage to five years and
#' six-monthly beyond, which is why the paper's UK real return maturities run
#' `60, 66, ..., 120`: those are exactly the published points.
#'
#' Everything between the half-yearly points is interpolated here. That is not
#' a liberty being taken with the data: the Bank's curve *is* a cubic spline
#' fitted to bond prices with a variable roughness penalty, so interpolating
#' its output with a spline recovers the same function rather than inventing
#' structure. A one-month holding return needs maturities `n` and `n - 1`
#' adjacent, so a monthly grid is unavoidable in any case.
#'
#' @section The ragged short end, which is a real limitation:
#' The nominal curve is not published all the way down on every date. Over
#' 1985-2012 the one-month point exists on only 190 of 336 month-ends,
#' starting in March 1997, and the shortest published maturity before then
#' ranges from three to twelve months. The Supplementary Appendix nonetheless
#' names the one-month nominal yield as its short rate over the full sample,
#' so reproducing it requires extending the curve below its shortest fitted
#' point.
#'
#' Maturities below the shortest published point on a given date are therefore
#' **extrapolated**, linearly rather than cubically so that the result cannot
#' fly off, and every such cell is flagged in the `extrapolated` column. Check
#' it before trusting anything at the short end. This is the same hazard as
#' the fitted one-month yield on the US nominal curve, which runs about 28
#' basis points above actual bill rates because bills are excluded from the
#' fit.
#'
#' @param type `"nominal"` or `"real"`.
#' @param maturities Maturity grid in months. The nominal default runs from 1
#'   month; the real default starts at 59 so that a one-month return on a
#'   five-year bond can be formed, which is where the paper's UK real returns
#'   begin.
#' @param align_to Optional vector of dates to relabel observations onto, as
#'   in [gsw_tips()]. Both Bank of England curves are published on the same
#'   month-ends, so this is only needed when joining to something else.
#' @param start,end Optional date bounds.
#' @param cache Optional path to the downloaded `.zip`. Read when it exists,
#'   written when it does not.
#' @param quiet Suppress the download progress bar.
#' @param check Verify the interpolation by holding out a published maturity
#'   and comparing against it. See below.
#'
#' @section The interpolation check:
#' Comparing interpolated values against the published points they were
#' interpolated *through* would compare a number with itself and pass
#' regardless. The check therefore drops the 84-month column, re-interpolates
#' without it, and compares -- a genuine out-of-sample test of whether the
#' spline reproduces the Bank's curve. It fails loudly above a basis point.
#'
#' @return A data frame with `date`, `maturity` (months), `yield` (percent)
#'   and `extrapolated`.
#'
#' @section Requirements:
#' Needs the `readxl` package and a working internet connection. Neither is a
#' hard dependency of termpremia, and nothing here is redistributed.
#'
#' @references
#' Anderson, N. and J. Sleath (2001). "New estimates of the UK real and
#' nominal yield curves." *Bank of England Quarterly Bulletin*, Spring.
#'
#' @seealso [atsm_real()], [ons_rpi()], [gsw_tips()]
#'
#' @examples
#' \dontrun{
#' nom <- boe_yield_curve("nominal", end = "2012-12-31")
#' rea <- boe_yield_curve("real", end = "2012-12-31")
#' mean(nom$extrapolated)
#' }
#'
#' @export
boe_yield_curve <- function(type = c("nominal", "real"),
                            maturities = NULL,
                            align_to = NULL,
                            start = NULL,
                            end = NULL,
                            cache = NULL,
                            quiet = TRUE,
                            check = TRUE) {
  type <- match.arg(type)
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("`boe_yield_curve()` needs the readxl package: ",
         "install.packages(\"readxl\").", call. = FALSE)
  }

  maturities <- maturities %||% if (type == "nominal") 1:120 else 59:120
  if (!is.numeric(maturities) || !length(maturities) ||
      any(maturities < 1 | maturities != as.integer(maturities))) {
    stop("`maturities` must be whole numbers of months, at least 1.",
         call. = FALSE)
  }
  maturities <- as.integer(maturities)

  zip_name <- if (type == "nominal") "glcnominalmonthedata.zip" else
    "glcrealmonthedata.zip"
  url <- paste0("https://www.bankofengland.co.uk/-/media/boe/files/",
                "statistics/yield-curves/", zip_name)

  path <- fetch_to_file(url, cache, quiet,
                        what = paste0("the Bank of England ", type,
                                      " yield curve"))

  ex <- file.path(tempdir(), paste0("boe-", type, "-", as.integer(Sys.time())))
  dir.create(ex, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(ex, recursive = TRUE), add = TRUE)

  files <- tryCatch(utils::unzip(path, exdir = ex),
                    error = function(e) {
                      stop("Could not unpack the Bank of England archive at\n  ",
                           path, "\n", conditionMessage(e), call. = FALSE)
                    })
  files <- files[grepl("\\.xlsx?$", files) & !grepl("^~\\$", basename(files))]
  if (!length(files)) {
    stop("The Bank of England archive contains no workbook.", call. = FALSE)
  }

  parts <- lapply(sort(files), boe_read_workbook)
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) {
    stop("No usable spot-curve sheet in the Bank of England archive. The ",
         "workbook layout may have changed.", call. = FALSE)
  }
  curve <- boe_bind(parts)

  keep <- rep(TRUE, length(curve$dates))
  if (!is.null(start)) keep <- keep & curve$dates >= as.Date(start)
  if (!is.null(end)) keep <- keep & curve$dates <= as.Date(end)
  curve$dates <- curve$dates[keep]
  curve$yield <- curve$yield[keep, , drop = FALSE]

  if (!length(curve$dates)) {
    stop("No observations left after applying `start`/`end`.", call. = FALSE)
  }

  if (check) boe_check_interpolation(curve)

  if (!is.null(align_to)) {
    align_to <- sort(unique(as.Date(align_to)))
    target <- align_to[match(format(curve$dates, "%Y-%m"),
                             format(align_to, "%Y-%m"))]
    dropped <- sum(is.na(target))
    if (dropped) {
      message("boe_yield_curve(): dropped ", dropped,
              " month(s) absent from `align_to`.")
    }
    curve$yield <- curve$yield[!is.na(target), , drop = FALSE]
    curve$dates <- target[!is.na(target)]
    if (!length(curve$dates)) {
      stop("`align_to` shares no calendar month with the Bank of England ",
           "curve.", call. = FALSE)
    }
  }

  grid <- boe_interpolate(curve, maturities)

  out <- data.frame(
    date = rep(curve$dates, times = length(maturities)),
    maturity = rep(maturities, each = length(curve$dates)),
    yield = as.vector(grid$yield),
    extrapolated = as.vector(grid$extrapolated),
    stringsAsFactors = FALSE
  )
  out <- out[order(out$date, out$maturity), , drop = FALSE]
  rownames(out) <- NULL

  attr(out, "source_url") <- url
  attr(out, "retrieved") <- Sys.Date()
  out
}

#' Read the two spot-curve sheets of one Bank of England workbook
#'
#' @keywords internal
#' @noRd
boe_read_workbook <- function(path) {
  sheets <- tryCatch(readxl::excel_sheets(path), error = function(e) character())
  want <- c("3. spot, short end", "4. spot curve")
  have <- intersect(want, sheets)
  if (!length(have)) return(NULL)

  parts <- lapply(have, function(s) boe_read_sheet(path, s))
  parts <- Filter(function(p) !is.null(p) && length(p$maturity), parts)
  if (!length(parts)) return(NULL)
  boe_merge_sheets(parts)
}

#' Read one sheet
#'
#' The maturity header is found by its label rather than by row number,
#' because the two sheets put it in different rows and the Bank has moved it
#' before. `"months:"` is preferred over `"years:"` where both exist: the
#' months row holds clean integers, while the years row on the short-end
#' sheet carries at least one value off by a factor of a hundred.
#'
#' @keywords internal
#' @noRd
boe_read_sheet <- function(path, sheet) {
  raw <- as.data.frame(readxl::read_excel(
    path, sheet = sheet, col_names = FALSE, .name_repair = "minimal"
  ))
  if (!ncol(raw) || !nrow(raw)) return(NULL)

  label <- trimws(as.character(raw[[1L]]))
  row_m <- which(label == "months:")[1L]
  row_y <- which(label == "years:")[1L]
  header <- if (!is.na(row_m)) row_m else row_y
  if (is.na(header)) return(NULL)

  mats <- suppressWarnings(as.numeric(unlist(raw[header, -1L])))
  if (is.na(row_m)) mats <- mats * 12

  serial <- suppressWarnings(as.numeric(label))
  rows <- which(!is.na(serial) & serial > 10000)
  if (!length(rows)) return(NULL)

  y <- as.matrix(raw[rows, -1L, drop = FALSE])
  storage.mode(y) <- "double"

  ok <- !is.na(mats)
  list(dates = as.Date(serial[rows], origin = "1899-12-30"),
       maturity = as.integer(round(mats[ok])),
       yield = y[, ok, drop = FALSE])
}

#' Merge the short-end and long-end sheets of one workbook
#'
#' They overlap; the short-end sheet wins on shared maturities because it is
#' the finer grid.
#'
#' @keywords internal
#' @noRd
boe_merge_sheets <- function(parts) {
  dates <- parts[[1L]]$dates
  for (p in parts) {
    if (!identical(p$dates, dates)) {
      stop("The spot-curve sheets of one Bank of England workbook disagree ",
           "on their dates, so they cannot be merged.", call. = FALSE)
    }
  }

  mats <- sort(unique(unlist(lapply(parts, `[[`, "maturity"))))
  out <- matrix(NA_real_, length(dates), length(mats))

  for (p in rev(parts)) {              # first part applied last, so it wins
    j <- match(p$maturity, mats)
    block <- out[, j, drop = FALSE]
    fill <- !is.na(p$yield)
    block[fill] <- p$yield[fill]
    out[, j] <- block
  }
  list(dates = dates, maturity = mats, yield = out)
}

#' Stack the workbooks covering different eras
#'
#' @keywords internal
#' @noRd
boe_bind <- function(parts) {
  mats <- sort(unique(unlist(lapply(parts, `[[`, "maturity"))))
  dates <- do.call(c, lapply(parts, `[[`, "dates"))

  out <- matrix(NA_real_, length(dates), length(mats))
  offset <- 0L
  for (p in parts) {
    rows <- offset + seq_along(p$dates)
    out[rows, match(p$maturity, mats)] <- p$yield
    offset <- offset + length(p$dates)
  }

  o <- order(dates)
  dates <- dates[o]
  out <- out[o, , drop = FALSE]

  dup <- duplicated(dates)
  if (any(dup)) {
    dates <- dates[!dup]
    out <- out[!dup, , drop = FALSE]
  }
  list(dates = dates, maturity = mats, yield = out)
}

#' Put a curve on a monthly maturity grid
#'
#' One natural cubic spline per date through that date's published points.
#' Natural rather than the default: outside the published range it continues
#' linearly instead of letting a cubic run away, which matters because the
#' nominal short end genuinely has to be extended on more than a third of the
#' sample.
#'
#' @keywords internal
#' @noRd
boe_interpolate <- function(curve, maturities) {
  n_dates <- length(curve$dates)
  out <- matrix(NA_real_, n_dates, length(maturities))
  extrap <- matrix(FALSE, n_dates, length(maturities))

  for (i in seq_len(n_dates)) {
    have <- which(is.finite(curve$yield[i, ]))
    if (length(have) < 4L) next        # too few points to spline through

    x <- curve$maturity[have]
    y <- curve$yield[i, have]
    out[i, ] <- stats::spline(x, y, xout = maturities, method = "natural")$y
    extrap[i, ] <- maturities < min(x) | maturities > max(x)
  }

  dimnames(out) <- dimnames(extrap) <-
    list(format(curve$dates), as.character(maturities))
  list(yield = out, extrapolated = extrap)
}

#' Verify the interpolation against a held-out published maturity
#'
#' Interpolated values agree with the points they were interpolated through by
#' construction, so comparing them proves nothing. Dropping a column first and
#' predicting it does.
#'
#' @keywords internal
#' @noRd
boe_check_interpolation <- function(curve, holdout = 84L, tol = 0.01) {
  j <- match(holdout, curve$maturity)
  if (is.na(j)) {
    stop("`check = TRUE` but the ", holdout, "-month column needed for the ",
         "hold-out test is not published, so the interpolation would go ",
         "unverified. Pass `check = FALSE` if that is acceptable.",
         call. = FALSE)
  }

  reduced <- list(dates = curve$dates,
                  maturity = curve$maturity[-j],
                  yield = curve$yield[, -j, drop = FALSE])
  predicted <- boe_interpolate(reduced, holdout)$yield[, 1L]
  actual <- curve$yield[, j]

  ok <- is.finite(predicted) & is.finite(actual)
  if (sum(ok) < 12L) {
    stop("Only ", sum(ok), " date(s) could be used for the interpolation ",
         "hold-out test; refusing to pass vacuously.", call. = FALSE)
  }

  worst <- max(abs(predicted[ok] - actual[ok]))
  if (worst > tol) {
    stop("Interpolating the Bank of England curve reproduces its own ",
         holdout, "-month point only to ", format(worst, digits = 3),
         " percentage points across ", sum(ok), " dates. The grid is too ",
         "coarse for a monthly model; do not use this data.", call. = FALSE)
  }
  invisible(worst)
}


# ONS ---------------------------------------------------------------------

#' Download the UK retail price index
#'
#' Fetches the Office for National Statistics' RPI All Items index, which is
#' the price series UK index-linked gilts are tied to and the one the
#' Supplementary Appendix of Abrahams et al. (2016) uses as \eqn{Q_t}.
#'
#' @section Why two series are downloaded:
#' The published index (`CHAW`) begins in **January 1987**, because that is
#' when it was rebased to 100, but the UK specification starts in January
#' 1985 and its first monthly inflation reading needs December 1984. The
#' twelve-month percentage change (`CZBH`) runs from 1948, so the index is
#' chained backwards through it:
#' \deqn{RPI_t = RPI_{t+12} / (1 + CZBH_{t+12}/100)}
#'
#' The chaining is verified, not assumed. The same relation applied *forwards*
#' from 1987 predicts the published index for later years, and that
#' reconstruction is checked against the real thing before any extension is
#' returned.
#'
#' @section How accurate the extension is:
#' `CZBH` is published to one decimal place, so each chaining step carries
#' about half a basis point of a percent of rounding error in the level. That
#' is harmless for the *mean* of inflation but not negligible in a single
#' month's change: two adjacent chained months round independently, which puts
#' roughly 80 basis points of annualised noise on the monthly inflation
#' reading. Reconstructed months are flagged in the `chained` column so the
#' sensitivity can be checked -- 25 of the 336 months of the UK sample are
#' affected, and `analysis/validate-uk.R` reports the fit with and without
#' them.
#'
#' @param extend_to Earliest date the index is required to cover. Backward
#'   chaining runs until it is reached. `NULL` returns the published index
#'   only.
#' @param quiet Suppress download progress bars.
#' @param check Verify the chaining relation against the published index.
#'
#' @return A data frame with `date` (first of the month), `value` (index) and
#'   `chained`, flagging months reconstructed rather than published.
#'
#' @seealso [atsm_real()], [boe_yield_curve()]
#'
#' @examples
#' \dontrun{
#' rpi <- ons_rpi(extend_to = "1984-12-01")
#' sum(rpi$chained)
#' }
#'
#' @export
ons_rpi <- function(extend_to = NULL, quiet = TRUE, check = TRUE) {
  index <- ons_timeseries("chaw", quiet)
  pct <- ons_timeseries("czbh", quiet)

  if (check) ons_check_chaining(index, pct)

  index$chained <- FALSE
  if (is.null(extend_to)) return(index)

  ons_chain_back(index, pct, extend_to)
}

#' Extend a price index backwards through its twelve-month percentage change
#'
#' \eqn{RPI_t = RPI_{t+12} / (1 + \Delta_{t+12}/100)}, applied a year at a
#' time until `extend_to` is reached. Kept separate from the download so the
#' arithmetic can be tested without a network.
#'
#' @param index Data frame with `date` and `value`, ordered.
#' @param pct Data frame with `date` and `value`, the twelve-month change in
#'   percent.
#' @param extend_to Earliest date required.
#' @return `index` with a `chained` column and earlier rows prepended.
#' @keywords internal
#' @noRd
ons_chain_back <- function(index, pct, extend_to) {
  extend_to <- as.Date(extend_to)
  if (is.null(index$chained)) index$chained <- FALSE

  while (min(index$date) > extend_to) {
    block <- seq(add_months(min(index$date), -12L), by = "month",
                 length.out = 12L)
    block <- block[block < min(index$date)]
    if (!length(block)) break

    ahead <- add_months(block, 12L)
    i_idx <- match(ahead, index$date)
    i_pct <- match(ahead, pct$date)
    usable <- !is.na(i_idx) & !is.na(i_pct) & !is.na(pct$value[i_pct])
    if (!any(usable)) break

    index <- rbind(
      data.frame(
        date = block[usable],
        value = index$value[i_idx[usable]] /
          (1 + pct$value[i_pct[usable]] / 100),
        chained = TRUE,
        stringsAsFactors = FALSE
      ),
      index
    )
    index <- index[order(index$date), , drop = FALSE]
  }

  if (min(index$date) > extend_to) {
    warning("The retail price index could only be extended back to ",
            format(min(index$date)), ", short of ", format(extend_to),
            ". The percentage-change series does not reach further.",
            call. = FALSE)
  }
  rownames(index) <- NULL
  index
}

#' Fetch one ONS monthly time series by its CDID
#'
#' @keywords internal
#' @noRd
ons_timeseries <- function(cdid, quiet = TRUE) {
  url <- sprintf(paste0("https://www.ons.gov.uk/generator?format=csv&uri=",
                        "/economy/inflationandpriceindices/timeseries/%s/mm23"),
                 tolower(cdid))

  dest <- tempfile(fileext = ".csv")
  on.exit(unlink(dest), add = TRUE)

  # The ONS rate-limits with a 429 rather than refusing outright, so a couple
  # of unhurried retries is the difference between working and not.
  ok <- FALSE
  for (attempt in seq_len(3L)) {
    ok <- tryCatch({
      utils::download.file(url, dest, mode = "wb", quiet = quiet)
      TRUE
    }, error = function(e) FALSE)
    if (ok) break
    Sys.sleep(5 * attempt)
  }
  if (!ok) {
    stop("Could not download ONS series '", cdid, "' from\n  ", url,
         "\nThe ONS rate-limits repeated requests; try again shortly.",
         call. = FALSE)
  }

  raw <- utils::read.csv(dest, header = FALSE, stringsAsFactors = FALSE,
                         colClasses = "character")
  if (ncol(raw) < 2L) {
    stop("ONS series '", cdid, "' came back in an unexpected layout.",
         call. = FALSE)
  }

  label <- trimws(raw[[1L]])
  monthly <- grepl("^[12][0-9]{3} [A-Za-z]{3}$", label)
  if (!any(monthly)) {
    stop("ONS series '", cdid, "' contains no monthly observations.",
         call. = FALSE)
  }

  month_abb <- toupper(base::month.abb)
  parts <- strsplit(label[monthly], " ", fixed = TRUE)
  year <- as.integer(vapply(parts, `[`, character(1L), 1L))
  mon <- match(toupper(vapply(parts, `[`, character(1L), 2L)), month_abb)

  out <- data.frame(
    date = as.Date(sprintf("%d-%02d-01", year, mon)),
    value = suppressWarnings(as.numeric(raw[[2L]][monthly])),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$date) & !is.na(out$value), , drop = FALSE]
  out <- out[order(out$date), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Check that the index and the percentage change describe the same series
#'
#' @keywords internal
#' @noRd
ons_check_chaining <- function(index, pct, tol = 0.002) {
  ahead <- add_months(index$date, -12L)
  i_prev <- match(ahead, index$date)
  i_pct <- match(index$date, pct$date)

  ok <- !is.na(i_prev) & !is.na(i_pct)
  if (sum(ok) < 24L) {
    stop("Only ", sum(ok), " month(s) let the RPI index and its percentage ",
         "change be compared; refusing to pass vacuously.", call. = FALSE)
  }

  implied <- index$value[i_prev[ok]] * (1 + pct$value[i_pct[ok]] / 100)

  # Relative, not absolute. The percentage change is published to one decimal
  # place, so a rounding of half a basis point of a percent lands as 0.2 index
  # points once the index itself is above 400 -- which is arithmetic, not
  # disagreement. A relative tolerance is scale-free and still catches the
  # failure this guards against: the wrong series, the wrong base year, or a
  # factor-of-100 slip.
  worst <- max(abs(implied / index$value[ok] - 1))

  if (worst > tol) {
    stop("The published RPI index and its twelve-month percentage change ",
         "disagree by up to ", format(100 * worst, digits = 3), "% across ",
         sum(ok), " months, so the index cannot be chained backwards through ",
         "the percentage change. Do not use this data.", call. = FALSE)
  }
  invisible(worst)
}
