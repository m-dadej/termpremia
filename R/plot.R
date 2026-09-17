#' Plot a fitted term structure decomposition
#'
#' Base graphics, so the package carries no plotting dependency.
#'
#' @param x An `atsm_fit`.
#' @param maturity Maturity in months. Defaults to the longest fitted.
#' @param type One of:
#'   \describe{
#'     \item{`"decomposition"`}{The observed yield and the risk-neutral
#'       (expected average short rate) component, with the term premium shaded
#'       between them. The canonical chart: the premium *is* the gap.}
#'     \item{`"premium"`}{The term premium alone, with a zero line. Negative
#'       premia are common post-2010 and worth being able to see.}
#'     \item{`"residuals"`}{Pricing errors over time and fit quality by
#'       maturity, for checking that the model actually fits the curve.}
#'   }
#' @param col Length-3 vector of colours for the yield, risk-neutral component
#'   and premium shading.
#' @param ... Passed to the underlying plotting calls.
#'
#' @return `x`, invisibly. Called for the side effect.
#'
#' @examples
#' panel <- yield_panel(gsw_monthly, units = "percent",
#'                      maturity_unit = "months", issuer = "US")
#' fit <- atsm(panel, n_factors = 5)
#'
#' plot(fit)                       # 10y decomposition
#' plot(fit, type = "premium")
#' plot(fit, maturity = 24, type = "decomposition")
#'
#' @export
plot.atsm_fit <- function(x,
                          maturity = NULL,
                          type = c("decomposition", "premium", "residuals"),
                          col = c("grey25", "steelblue", "#E8A0A0"),
                          ...) {
  type <- match.arg(type)
  maturity <- maturity %||% max(x$maturities)

  idx <- match(maturity, x$maturities)
  if (is.na(idx)) {
    stop("Maturity ", maturity, " was not fitted. Available: ",
         min(x$maturities), "-", max(x$maturities), " months.", call. = FALSE)
  }

  lab <- if (maturity %% 12 == 0) paste0(maturity %/% 12, "y") else paste0(maturity, "m")

  switch(
    type,
    decomposition = plot_decomposition(x, idx, lab, col, ...),
    premium       = plot_premium(x, idx, lab, col, ...),
    residuals     = plot_residuals(x, idx, lab, col, ...)
  )

  invisible(x)
}

#' @noRd
plot_decomposition <- function(x, idx, lab, col, ...) {
  obs <- x$observed[, idx] * 100
  rn <- x$risk_neutral[, idx] * 100
  d <- x$dates

  plot(d, obs, type = "n", ylim = range(c(obs, rn, 0), na.rm = TRUE),
       xlab = "", ylab = "percent",
       main = paste0(lab, " yield decomposition", curve_tag(x)), ...)

  # The premium is the area between the yield and its risk-neutral component.
  polygon(c(d, rev(d)), c(obs, rev(rn)), col = col[3], border = NA)
  abline(h = 0, col = "grey85", lty = 2)

  lines(d, obs, lwd = 2, col = col[1])
  lines(d, rn, lwd = 1.5, col = col[2])

  legend("topright", bty = "n", cex = 0.8,
         legend = c("observed yield", "expected average short rate",
                    "term premium"),
         col = c(col[1], col[2], col[3]),
         lwd = c(2, 1.5, NA), pch = c(NA, NA, 15), pt.cex = 1.4)
}

#' @noRd
plot_premium <- function(x, idx, lab, col, ...) {
  tp <- x$term_premium[, idx] * 1e4
  d <- x$dates

  plot(d, tp, type = "n", xlab = "", ylab = "basis points",
       main = paste0(lab, " term premium", curve_tag(x)), ...)

  # Negative premia are a real and much-discussed feature of the post-2010
  # period, so make the sign visible rather than incidental.
  usr <- par("usr")
  if (usr[3] < 0) {
    rect(usr[1], usr[3], usr[2], 0, col = "grey95", border = NA)
  }
  abline(h = 0, col = "grey60", lty = 2)

  lines(d, tp, lwd = 1.8, col = col[2])
  box()
}

#' @noRd
plot_residuals <- function(x, idx, lab, col, ...) {
  op <- par(mfrow = c(2, 1), mar = c(3, 4, 2.5, 1))
  on.exit(par(op), add = TRUE)

  resid_bp <- (x$observed - x$fitted) * 1e4

  plot(x$dates, resid_bp[, idx], type = "l", col = col[2],
       xlab = "", ylab = "basis points",
       main = paste0(lab, " pricing error", curve_tag(x)), ...)
  abline(h = 0, col = "grey60", lty = 2)

  rmse <- sqrt(colMeans(resid_bp^2))
  plot(x$maturities, rmse, type = "l", lwd = 2, col = col[1],
       xlab = "maturity (months)", ylab = "bp",
       main = "RMSE by maturity", ylim = c(0, max(rmse) * 1.05))
  points(x$maturities[idx], rmse[idx], pch = 19, col = col[2])
}

#' @noRd
curve_tag <- function(x) {
  if (is.na(x$issuer) || !nzchar(x$issuer)) "" else paste0(" (", x$issuer, ")")
}
