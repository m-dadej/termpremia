# predict() and plot() for the joint real-nominal fit.
#
# Kept apart from R/atsm-real.R, which is already long, and because both of
# these are about applying a finished fit rather than producing one.

#' Apply a fitted joint model to new observations
#'
#' Recomputes the decomposition at dates outside the estimation sample,
#' holding every estimated parameter fixed. This is how a model estimated on
#' month-end data is evaluated at a higher frequency: the parameters need
#' monthly excess returns to estimate, but the decomposition itself needs only
#' the factors on the date in question.
#'
#' @section What is held fixed:
#' Everything the estimation produced, including the whole factor
#' construction: the principal component loadings and centres, the
#' orthogonalising regression's coefficients, and the scaling. New yields are
#' pushed through that fixed pipeline rather than re-extracted, because
#' re-extracting would put the factors on a different basis from the one the
#' parameters were estimated against.
#'
#' @param object An `atsm_real_fit`.
#' @param newdata A [yield_panel] carrying both curves on the estimation
#'   grids.
#' @param liquidity Required when the model was fitted with a liquidity
#'   factor: a data frame with `date` and `value`, or a numeric vector
#'   matching `newdata`'s dates.
#' @param ... Unused.
#'
#' @return A data frame with `date`, `maturity` and one column per component.
#'   Maturities the real curve does not reach carry `NA` in the real and
#'   inflation columns rather than being dropped, so the frame stays
#'   rectangular.
#'
#' @seealso [atsm_real()], [forward_rate()]
#'
#' @examples
#' \dontrun{
#' fit <- atsm_real(monthly_panel, inflation = cpi)
#' predict(fit, weekly_panel)
#' }
#'
#' @export
predict.atsm_real_fit <- function(object, newdata, liquidity = NULL, ...) {
  if (!inherits(newdata, "yield_panel")) {
    stop("`newdata` must be a yield_panel.", call. = FALSE)
  }

  nom <- object$curves[["nominal"]]
  rea <- object$curves[["real"]]
  for (nm in c(nom, rea)) {
    if (!nm %in% newdata$meta$curve) {
      stop("`newdata` has no curve named '", nm, "'. The fit used '", nom,
           "' and '", rea, "'.", call. = FALSE)
    }
  }

  grid_n <- newdata$maturities[[nom]]
  grid_r <- newdata$maturities[[rea]]
  if (!all(object$maturities %in% grid_n) ||
      !all(object$real_maturities %in% grid_r)) {
    stop("`newdata` is missing maturities the model was fitted on.",
         call. = FALSE)
  }

  y_n <- curve_matrix(newdata, nom)[, match(object$maturities, grid_n),
                                    drop = FALSE] / 12
  y_r <- curve_matrix(newdata, rea)[, match(object$real_maturities, grid_r),
                                    drop = FALSE] / 12

  liq <- resolve_liquidity(liquidity, newdata$dates)
  if (object$has_liquidity && is.null(liq)) {
    stop("This model was fitted with a liquidity factor, so `predict()` needs ",
         "one for the new dates too. Pass it through `liquidity`.",
         call. = FALSE)
  }

  keep <- stats::complete.cases(y_n) & stats::complete.cases(y_r)
  if (!is.null(liq)) keep <- keep & !is.na(liq)
  if (!any(keep)) {
    stop("No date in `newdata` has both curves (and the liquidity factor) ",
         "observed.", call. = FALSE)
  }

  dates <- newdata$dates[keep]
  x <- project_real_factors(object, y_n[keep, , drop = FALSE],
                            y_r[keep, , drop = FALSE],
                            if (is.null(liq)) NULL else liq[keep])

  co_q <- object$recursion$q
  co_p <- object$recursion$p
  mats <- object$maturities
  mats_r <- object$real_maturities

  fit_n <- acmy_yields(co_q$a, co_q$b, x, mats) * 12
  rn_n <- acmy_yields(co_p$a, co_p$b, x, mats) * 12
  fit_r <- acmy_yields(co_q$a_real, co_q$b_real, x, mats_r) * 12
  rn_r <- acmy_yields(co_p$a_real, co_p$b_real, x, mats_r) * 12

  n <- length(dates)
  widen <- function(m) {
    j <- match(mats, mats_r)
    out <- matrix(NA_real_, n, length(mats))
    ok <- !is.na(j)
    out[, ok] <- m[, j[ok], drop = FALSE]
    out
  }

  fit_rw <- widen(fit_r)
  rn_rw <- widen(rn_r)

  data.frame(
    date = rep(dates, times = length(mats)),
    maturity = rep(mats, each = n),
    fitted = as.vector(fit_n),
    risk_neutral = as.vector(rn_n),
    term_premium = as.vector(fit_n - rn_n),
    fitted_real = as.vector(fit_rw),
    risk_neutral_real = as.vector(rn_rw),
    term_premium_real = as.vector(fit_rw - rn_rw),
    breakeven = as.vector(fit_n - fit_rw),
    expected_inflation = as.vector(rn_n - rn_rw),
    inflation_risk_premium = as.vector((fit_n - rn_n) - (fit_rw - rn_rw)),
    stringsAsFactors = FALSE
  )
}

#' Push new yields through the fitted factor construction
#'
#' The same pipeline `acmy_factors()` ran, with every learned quantity held at
#' its estimated value. Shares its arithmetic with `acmy_implied_factors()`,
#' which does this for model-implied yields rather than new observed ones.
#'
#' @keywords internal
#' @noRd
project_real_factors <- function(object, y_nom, y_real, liquidity) {
  fac <- object$pca

  x_nom <- sweep(y_nom, 2L, fac$nominal$center, "-") %*% fac$nominal$loadings

  z <- cbind(1, x_nom)
  if (!is.na(fac$liq_index)) z <- cbind(z, liquidity)

  resid <- y_real - z %*% fac$ortho_coef
  x_real <- sweep(resid, 2L, fac$real$center, "-") %*% fac$real$loadings

  x <- cbind(x_nom, x_real)
  if (!is.na(fac$liq_index)) x <- cbind(x, liquidity)

  sweep(x, 2L, fac$scale, "/")
}


# Plotting ----------------------------------------------------------------

#' Plot a fitted joint real-nominal decomposition
#'
#' @param x An `atsm_real_fit`.
#' @param maturity Maturity in months. Defaults to the longest maturity the
#'   chosen `type` supports.
#' @param type One of:
#'   \describe{
#'     \item{`"inflation"`}{Breakeven inflation split into expected inflation
#'       and the inflation risk premium. The chart the joint model exists to
#'       draw, and the default.}
#'     \item{`"real"`}{Both term premia together, nominal and real, with the
#'       inflation risk premium as the gap between them.}
#'     \item{`"decomposition"`}{The nominal curve alone, as [plot.atsm_fit()]
#'       draws it.}
#'     \item{`"premium"`}{The nominal term premium alone.}
#'     \item{`"residuals"`}{Nominal pricing errors.}
#'   }
#' @param col Length-3 vector of colours.
#' @param ... Passed to the underlying plotting calls.
#'
#' @return `x`, invisibly. Called for the side effect.
#'
#' @examples
#' \dontrun{
#' fit <- atsm_real(panel, inflation = cpi)
#' plot(fit)                     # breakeven decomposition
#' plot(fit, type = "real")      # the two term premia
#' }
#'
#' @export
plot.atsm_real_fit <- function(x,
                               maturity = NULL,
                               type = c("inflation", "real", "decomposition",
                                        "premium", "residuals"),
                               col = c("grey25", "steelblue", "#E8A0A0"),
                               ...) {
  type <- match.arg(type)
  shared_types <- c("inflation", "real")

  grid <- if (type %in% shared_types) x$shared_maturities else x$maturities
  maturity <- maturity %||% max(grid)

  idx <- match(maturity, grid)
  if (is.na(idx)) {
    stop("Maturity ", maturity, " is not available for type = \"", type,
         "\". Available: ", min(grid), "-", max(grid), " months.",
         call. = FALSE)
  }

  lab <- if (maturity %% 12 == 0) {
    paste0(maturity %/% 12, "y")
  } else {
    paste0(maturity, "m")
  }

  if (type == "inflation") {
    plot_inflation(x, idx, lab, col, ...)
  } else if (type == "real") {
    plot_both_premia(x, maturity, lab, col, ...)
  } else {
    # The nominal curve of a joint fit plots exactly as a single-curve fit
    # does, so the existing panels are reused rather than duplicated.
    shim <- structure(
      list(dates = x$dates, maturities = x$maturities, observed = x$observed,
           fitted = x$fitted, risk_neutral = x$risk_neutral,
           term_premium = x$term_premium, curve = x$curves[["nominal"]],
           # curve_tag() reads these for the panel title; a joint fit stores
           # them per curve, so the nominal curve's are the ones to pass on.
           issuer = x$issuer %||% NA_character_,
           instrument = unname(x$instrument["nominal"]) %||% NA_character_),
      class = "atsm_fit"
    )
    switch(type,
           decomposition = plot_decomposition(shim, idx, lab, col, ...),
           premium       = plot_premium(shim, idx, lab, col, ...),
           residuals     = plot_residuals(shim, idx, lab, col, ...))
  }

  invisible(x)
}

#' Breakeven inflation as expectations plus risk premium
#'
#' @keywords internal
#' @noRd
plot_inflation <- function(x, idx, lab, col, ...) {
  be <- x$breakeven[, idx] * 100
  ei <- x$expected_inflation[, idx] * 100
  irp <- x$inflation_risk_premium[, idx] * 100

  graphics::plot(x$dates, be, type = "n", xlab = "", ylab = "percent",
                 main = paste0(lab, " breakeven inflation"),
                 ylim = range(c(be, ei, irp, 0), na.rm = TRUE), ...)
  graphics::abline(h = 0, col = "grey80")
  graphics::polygon(c(x$dates, rev(x$dates)), c(ei, rev(be)),
                    col = col[3L], border = NA)
  graphics::lines(x$dates, be, col = col[1L], lwd = 2)
  graphics::lines(x$dates, ei, col = col[2L], lwd = 2)
  graphics::lines(x$dates, irp, col = col[2L], lty = 3)
  graphics::legend("topright", bty = "n", lwd = c(2, 2, 1), lty = c(1, 1, 3),
                   col = c(col[1L], col[2L], col[2L]),
                   legend = c("breakeven", "expected inflation",
                              "inflation risk premium"))
  graphics::box()
}

#' The two term premia, with the inflation risk premium as the gap
#'
#' @keywords internal
#' @noRd
plot_both_premia <- function(x, maturity, lab, col, ...) {
  tp_n <- x$term_premium[, match(maturity, x$maturities)] * 100
  tp_r <- x$term_premium_real[, match(maturity, x$real_maturities)] * 100

  graphics::plot(x$dates, tp_n, type = "n", xlab = "", ylab = "percent",
                 main = paste0(lab, " term premia"),
                 ylim = range(c(tp_n, tp_r, 0), na.rm = TRUE), ...)
  graphics::abline(h = 0, col = "grey80")
  graphics::polygon(c(x$dates, rev(x$dates)), c(tp_r, rev(tp_n)),
                    col = col[3L], border = NA)
  graphics::lines(x$dates, tp_n, col = col[1L], lwd = 2)
  graphics::lines(x$dates, tp_r, col = col[2L], lwd = 2)
  graphics::legend("topright", bty = "n", lwd = 2,
                   col = c(col[1L], col[2L]),
                   legend = c("nominal", "real (gap = inflation risk premium)"))
  graphics::box()
}
