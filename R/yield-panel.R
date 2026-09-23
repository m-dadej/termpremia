#' Construct a zero-coupon yield panel
#'
#' Builds the object every estimator in this package consumes: one or more
#' zero-coupon yield curves observed on a common set of dates. The constructor
#' is deliberately strict, because the mistakes it catches -- yields supplied in
#' percent when the code expects decimals, maturities in years when it expects
#' months -- do not produce errors downstream. They produce plausible-looking
#' term premia that are wrong by orders of magnitude.
#'
#' @section Canonical form:
#' Whatever you supply, the panel stores yields as **decimals** (0.045, not 4.5)
#' and maturities in **months**. Accessors return the canonical form.
#'
#' @section Multiple curves:
#' A panel may hold several curves -- nominal, TIPS, OIS, a foreign sovereign.
#' They must share a common date index. Curves are *not* required to share a
#' maturity grid: TIPS are typically quoted over a shorter range than nominals.
#'
#' Each curve carries an `instrument` and an `issuer`. This is not decoration.
#' The choice of benchmark materially moves estimated term premia: Hordahl and
#' Tristani use French government yields for the euro area, Adrian, Crump and
#' Moench use German government yields, and the Bank of France uses EONIA OIS,
#' differing by roughly 25 basis points on average before any modelling is done.
#' Recording the instrument lets later comparisons refuse to mix them silently.
#'
#' @param x A data frame in long format (one row per date-maturity), or a
#'   numeric matrix in wide format (rows = dates, columns = maturities).
#' @param ... Passed to methods.
#'
#' @return An object of class `yield_panel`.
#'
#' @examples
#' long <- expand.grid(date = seq(as.Date("2020-01-31"), by = "month",
#'                                length.out = 24),
#'                     maturity = c(12, 24, 60, 120))
#' long$yield <- 2 + long$maturity / 120 + rnorm(nrow(long), sd = 0.05)
#' yield_panel(long, maturity_unit = "months", units = "percent")
#'
#' @export
yield_panel <- function(x, ...) {
  UseMethod("yield_panel")
}

#' @rdname yield_panel
#'
#' @param date,maturity,yield,curve Column names in `x`. `curve` may be `NULL`
#'   for a single-curve panel.
#' @param units One of `"auto"`, `"percent"` or `"decimal"`. `"auto"` infers
#'   from the magnitude of the data and always reports what it inferred --
#'   set it explicitly if the inference is wrong.
#' @param maturity_unit Either `"months"` or `"years"`. No default: years and
#'   months are not distinguishable from the data when maturities are small.
#' @param instrument Instrument type, recycled across curves or a named
#'   character vector keyed by curve name. One of `"government"`, `"ois"`,
#'   `"swap"`, `"tips"` or `"other"`.
#' @param issuer Free-text issuer or market label, for example `"US"` or
#'   `"DE"`. Recycled or named like `instrument`.
#' @param extrapolated Optional name of a logical column marking yields that
#'   lie outside the maturities actually observed on that date, such as the
#'   `extrapolated` column of [gsw_monthly] or [boe_yield_curve()]. [atsm()]
#'   reads it to warn when its default short rate is an extrapolation.
#'   [predict.svensson_fit()] sets it for you.
#'
#' @export
yield_panel.data.frame <- function(x,
                                   date = "date",
                                   maturity = "maturity",
                                   yield = "yield",
                                   curve = NULL,
                                   units = c("auto", "percent", "decimal"),
                                   maturity_unit = c("months", "years"),
                                   instrument = "government",
                                   issuer = NA_character_,
                                   extrapolated = NULL,
                                   ...) {
  units <- match.arg(units)
  maturity_unit <- match.arg(maturity_unit)

  needed <- c(date, maturity, yield, curve, extrapolated)
  missing_cols <- setdiff(needed, names(x))
  if (length(missing_cols)) {
    stop("`x` is missing column(s): ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  dates <- as_date_strict(x[[date]], arg = date)
  mats <- as_maturity_months(x[[maturity]], maturity_unit, arg = maturity)
  yld <- as_yield_decimal(x[[yield]], units, arg = yield)

  curve_id <- if (is.null(curve)) {
    rep("yields", nrow(x))
  } else {
    as.character(x[[curve]])
  }
  if (anyNA(curve_id)) {
    stop("`", curve, "` must not contain NA.", call. = FALSE)
  }

  long <- data.frame(
    date = dates, maturity = mats, yield = yld, curve = curve_id,
    stringsAsFactors = FALSE
  )

  if (!is.null(extrapolated)) {
    flag <- x[[extrapolated]]
    if (!is.logical(flag) || anyNA(flag)) {
      stop("`", extrapolated, "` must be a logical column with no NA.",
           call. = FALSE)
    }
    long$extrapolated <- flag
  }

  dup <- duplicated(long[, c("date", "maturity", "curve")])
  if (any(dup)) {
    stop(
      sum(dup), " duplicated date/maturity/curve combination(s) in `x`; ",
      "each must appear at most once.",
      call. = FALSE
    )
  }

  new_yield_panel_from_long(long, instrument = instrument, issuer = issuer)
}

#' @rdname yield_panel
#'
#' @param dates A `Date` vector of length `nrow(x)`. If `NULL`, taken from
#'   `rownames(x)`.
#' @param maturities A numeric vector of length `ncol(x)`. If `NULL`, taken
#'   from `colnames(x)`.
#' @param curve_name Name for the single curve a matrix represents.
#'
#' @export
yield_panel.matrix <- function(x,
                               dates = NULL,
                               maturities = NULL,
                               curve_name = "yields",
                               units = c("auto", "percent", "decimal"),
                               maturity_unit = c("months", "years"),
                               instrument = "government",
                               issuer = NA_character_,
                               ...) {
  units <- match.arg(units)
  maturity_unit <- match.arg(maturity_unit)

  if (!is.numeric(x)) {
    stop("`x` must be a numeric matrix.", call. = FALSE)
  }

  if (is.null(dates)) {
    if (is.null(rownames(x))) {
      stop("`dates` is NULL and `x` has no rownames to take dates from.",
           call. = FALSE)
    }
    dates <- rownames(x)
  }
  if (is.null(maturities)) {
    if (is.null(colnames(x))) {
      stop("`maturities` is NULL and `x` has no colnames to take maturities ",
           "from.", call. = FALSE)
    }
    maturities <- colnames(x)
  }

  dates <- as_date_strict(dates, arg = "dates")
  maturities <- as_maturity_months(maturities, maturity_unit, arg = "maturities")

  if (length(dates) != nrow(x)) {
    stop("`dates` has length ", length(dates), " but `x` has ", nrow(x),
         " rows.", call. = FALSE)
  }
  if (length(maturities) != ncol(x)) {
    stop("`maturities` has length ", length(maturities), " but `x` has ",
         ncol(x), " columns.", call. = FALSE)
  }

  long <- data.frame(
    date = rep(dates, times = ncol(x)),
    maturity = rep(maturities, each = nrow(x)),
    yield = as_yield_decimal(as.vector(x), units, arg = "x"),
    curve = curve_name,
    stringsAsFactors = FALSE
  )

  new_yield_panel_from_long(long, instrument = instrument, issuer = issuer)
}

# data.table and tibble both inherit from data.frame, so they dispatch to the
# data.frame method without needing their own.


# Construction ------------------------------------------------------------

#' Assemble a yield_panel from validated long data
#' @noRd
new_yield_panel_from_long <- function(long, instrument, issuer) {
  curve_names <- unique(long$curve)
  dates <- sort(unique(long$date))

  meta <- build_curve_meta(curve_names, instrument, issuer)

  curves <- list()
  maturities <- list()
  # Extrapolation flags, when supplied, are kept as a parallel logical matrix
  # per curve. A cell with no row in `long` is missing, not extrapolated.
  flags <- if ("extrapolated" %in% names(long)) list()

  for (nm in curve_names) {
    sub <- long[long$curve == nm, , drop = FALSE]
    mats <- sort(unique(sub$maturity))

    m <- matrix(
      NA_real_,
      nrow = length(dates), ncol = length(mats),
      dimnames = list(as.character(dates), as.character(mats))
    )
    cell <- cbind(match(sub$date, dates), match(sub$maturity, mats))
    m[cell] <- sub$yield

    curves[[nm]] <- m
    maturities[[nm]] <- mats

    if (!is.null(flags)) {
      f <- matrix(FALSE, nrow(m), ncol(m), dimnames = dimnames(m))
      f[cell] <- sub$extrapolated
      flags[[nm]] <- f
    }
  }

  structure(
    list(
      curves = curves,
      dates = dates,
      maturities = maturities,
      meta = meta,
      frequency = detect_frequency(dates),
      extrapolated = flags
    ),
    class = "yield_panel"
  )
}

#' Build and validate per-curve metadata
#' @noRd
build_curve_meta <- function(curve_names, instrument, issuer) {
  valid <- c("government", "ois", "swap", "tips", "other")

  expand <- function(v, nm) {
    if (!is.null(names(v))) {
      idx <- match(curve_names, names(v))
      if (anyNA(idx)) {
        stop("`", nm, "` is named but has no entry for curve(s): ",
             paste(curve_names[is.na(idx)], collapse = ", "), call. = FALSE)
      }
      return(unname(v[idx]))
    }
    if (length(v) == 1L) return(rep(v, length(curve_names)))
    if (length(v) == length(curve_names)) return(v)
    stop("`", nm, "` must be length 1, length ", length(curve_names),
         ", or a named vector keyed by curve name.", call. = FALSE)
  }

  instrument <- as.character(expand(instrument, "instrument"))
  issuer <- as.character(expand(issuer, "issuer"))

  bad <- setdiff(instrument, valid)
  if (length(bad)) {
    stop("Unknown `instrument` value(s): ", paste(unique(bad), collapse = ", "),
         ". Must be one of: ", paste(valid, collapse = ", "), ".", call. = FALSE)
  }

  data.frame(
    curve = curve_names, instrument = instrument, issuer = issuer,
    stringsAsFactors = FALSE
  )
}


# Coercion and validation helpers -----------------------------------------

#' @noRd
as_date_strict <- function(x, arg) {
  if (inherits(x, "Date")) {
    d <- x
  } else if (inherits(x, "POSIXt")) {
    d <- as.Date(x)
  } else if (is.character(x) || is.factor(x)) {
    d <- as.Date(as.character(x))
  } else {
    stop("`", arg, "` must be a Date, POSIXt, or character vector of dates; ",
         "got ", class(x)[1], ". Numeric input is refused because its origin ",
         "is ambiguous.", call. = FALSE)
  }

  if (anyNA(d)) {
    stop("`", arg, "` contains ", sum(is.na(d)), " value(s) that could not be ",
         "parsed as dates.", call. = FALSE)
  }
  d
}

#' @noRd
as_maturity_months <- function(x, maturity_unit, arg) {
  m <- suppressWarnings(as.numeric(as.character(x)))

  if (anyNA(m)) {
    stop("`", arg, "` contains value(s) that are not numeric maturities.",
         call. = FALSE)
  }
  if (any(m <= 0)) {
    stop("`", arg, "` must be strictly positive; maturity 0 is not a bond.",
         call. = FALSE)
  }
  if (any(!is.finite(m))) {
    stop("`", arg, "` must be finite.", call. = FALSE)
  }

  if (maturity_unit == "years") m <- m * 12

  m
}

#' @noRd
as_yield_decimal <- function(x, units, arg) {
  y <- as.numeric(x)

  if (all(is.na(y))) {
    stop("`", arg, "` is entirely missing.", call. = FALSE)
  }

  rng <- range(y, na.rm = TRUE)
  if (any(is.infinite(y), na.rm = TRUE)) {
    stop("`", arg, "` contains infinite values.", call. = FALSE)
  }

  if (units == "auto") {
    units <- if (max(abs(rng)) > 1) "percent" else "decimal"
    # Named after the argument rather than after yield_panel(), because survey
    # forecasts and external short rates come through here too.
    message(
      "Inferred that `", arg, "` is in ", units,
      " (observed range ", format(rng[1], digits = 3), " to ",
      format(rng[2], digits = 3), "). ",
      "Set the units explicitly if that is wrong."
    )
  }

  if (units == "percent") y <- y / 100

  # After conversion, plausible nominal yields sit well inside +/-100%.
  post <- range(y, na.rm = TRUE)
  if (max(abs(post)) > 1) {
    warning(
      "Yields range from ", format(post[1], digits = 3), " to ",
      format(post[2], digits = 3),
      " after conversion to decimals, which is implausible for a bond yield. ",
      "Check the `units` argument.",
      call. = FALSE
    )
  }

  y
}

#' @noRd
detect_frequency <- function(dates) {
  if (length(dates) < 3L) return("unknown")

  gap <- stats::median(as.numeric(diff(sort(dates))))

  if (gap <= 4) {
    "daily"
  } else if (gap <= 10) {
    "weekly"
  } else if (gap <= 45) {
    "monthly"
  } else if (gap <= 135) {
    "quarterly"
  } else if (gap <= 420) {
    "annual"
  } else {
    "irregular"
  }
}


# Methods -----------------------------------------------------------------

#' @export
print.yield_panel <- function(x, ...) {
  cat("<yield_panel>\n")
  cat("  dates      : ", length(x$dates), " (", x$frequency, ", ",
      format(min(x$dates)), " to ", format(max(x$dates)), ")\n", sep = "")
  cat("  curves     : ", nrow(x$meta), "\n", sep = "")

  for (i in seq_len(nrow(x$meta))) {
    nm <- x$meta$curve[i]
    mats <- x$maturities[[nm]]
    m <- x$curves[[nm]]
    pct_na <- 100 * mean(is.na(m))

    flags <- x$extrapolated[[nm]]
    cat(sprintf(
      "    %-12s %s%s  maturities %g-%gm (n=%d)  missing %.1f%%%s\n",
      nm,
      x$meta$instrument[i],
      if (is.na(x$meta$issuer[i])) "" else paste0("/", x$meta$issuer[i]),
      min(mats), max(mats), length(mats), pct_na,
      if (is.null(flags)) "" else
        sprintf("  extrapolated %.1f%%", 100 * mean(flags))
    ))
  }

  cat("  yields     : decimal (e.g. 0.045 = 4.5%)\n")
  invisible(x)
}

#' @export
summary.yield_panel <- function(object, ...) {
  rows <- lapply(seq_len(nrow(object$meta)), function(i) {
    nm <- object$meta$curve[i]
    m <- object$curves[[nm]]
    data.frame(
      curve = nm,
      instrument = object$meta$instrument[i],
      issuer = object$meta$issuer[i],
      n_maturities = ncol(m),
      min_maturity_m = min(object$maturities[[nm]]),
      max_maturity_m = max(object$maturities[[nm]]),
      pct_missing = round(100 * mean(is.na(m)), 2),
      min_yield = round(min(m, na.rm = TRUE), 5),
      max_yield = round(max(m, na.rm = TRUE), 5),
      stringsAsFactors = FALSE
    )
  })

  structure(
    list(
      n_dates = length(object$dates),
      date_range = range(object$dates),
      frequency = object$frequency,
      curves = do.call(rbind, rows)
    ),
    class = "summary.yield_panel"
  )
}

#' @export
print.summary.yield_panel <- function(x, ...) {
  cat("<yield_panel> summary\n")
  cat("  ", x$n_dates, " dates (", x$frequency, "), ",
      format(x$date_range[1]), " to ", format(x$date_range[2]), "\n\n", sep = "")
  print(x$curves, row.names = FALSE)
  invisible(x)
}

#' Extract a curve as a matrix
#'
#' @param x A [yield_panel].
#' @param curve Name of the curve. Defaults to the first (and usually only) one.
#' @return A numeric matrix of decimal yields, dates in rows and maturities
#'   (in months) in columns.
#' @export
curve_matrix <- function(x, curve = NULL) {
  stopifnot(inherits(x, "yield_panel"))

  curve <- curve %||% x$meta$curve[1]
  if (!curve %in% names(x$curves)) {
    stop("No curve named '", curve, "'. Available: ",
         paste(names(x$curves), collapse = ", "), call. = FALSE)
  }
  x$curves[[curve]]
}

#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
