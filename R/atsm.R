#' Fit an affine term structure model
#'
#' Estimates a Gaussian affine term structure model on a zero-coupon yield
#' panel and decomposes yields into an expected-average-short-rate component
#' and a term premium.
#'
#' @section The decomposition:
#' \deqn{y^{(n)}_t = \underbrace{-(A^{RF}_n + B^{RF\prime}_n X_t)/n}_{\textrm{risk
#' neutral}} + \textrm{term premium}}
#'
#' The risk-neutral component is what the yield would be if investors demanded
#' no compensation for interest rate risk; it is the expected average short rate
#' over the life of the bond. The term premium is the residual. By construction
#' `fitted = risk_neutral + term_premium` exactly, and this is asserted in the
#' package's tests.
#'
#' @section A caution about levels:
#' Different term structure models disagree substantially about the *level* of
#' the term premium while agreeing closely on its *direction*. Cohen, Hordahl
#' and Xia (2018) find gaps of up to 200 basis points between published
#' estimates for the same market and date, alongside monthly-change correlations
#' of 0.77 to 0.92. Treat a single model's level with scepticism.
#'
#' ACM is additionally known to over-react to changes in the general level of
#' rates, because it uses yield information only: it can read a persistent fall
#' in rates as evidence that the steady state itself has moved, exaggerating
#' distant-horizon expectations. Survey-augmented dynamics are the usual remedy.
#'
#' @param panel A [yield_panel].
#' @param pricing Pricing model. Currently only `"acm"`, the three-step
#'   regression estimator of Adrian, Crump and Moench (2013).
#' @param n_factors Number of principal components used as pricing factors.
#'   ACM's baseline, and the New York Fed's published model, use 5.
#' @param curve Which curve in `panel` to fit. Defaults to the first.
#' @param maturities Maturities, in months, entering the factor extraction.
#'   Defaults to every maturity on the panel's grid.
#' @param return_maturities Maturities, in months, for which excess returns are
#'   formed. Each requires its predecessor on the grid. Defaults to every
#'   maturity from 2 upwards.
#' @param p_dynamics Estimator for the real-world factor dynamics. Currently
#'   only `"ols"`.
#'
#' @return An object of class `atsm_fit`.
#'
#' @references
#' Adrian, T., R. K. Crump and E. Moench (2013). "Pricing the term structure
#' with linear regressions." *Journal of Financial Economics* 110(1), 110-138.
#'
#' Cohen, B., P. Hordahl and D. Xia (2018). "Term premia: models and some
#' stylised facts." *BIS Quarterly Review*, September.
#'
#' @examples
#' panel <- yield_panel(gsw_monthly, units = "percent",
#'                      maturity_unit = "months", issuer = "US")
#' fit <- atsm(panel, n_factors = 5)
#' fit
#'
#' tp <- term_premium(fit, maturity = 120)
#' tail(tp)
#'
#' @export
atsm <- function(panel,
                 pricing = c("acm"),
                 n_factors = 5L,
                 curve = NULL,
                 maturities = NULL,
                 return_maturities = NULL,
                 p_dynamics = c("ols")) {
  pricing <- match.arg(pricing)
  p_dynamics <- match.arg(p_dynamics)

  if (!inherits(panel, "yield_panel")) {
    stop("`panel` must be a yield_panel; see ?yield_panel.", call. = FALSE)
  }

  curve <- curve %||% panel$meta$curve[1L]
  y_ann <- curve_matrix(panel, curve)
  grid <- panel$maturities[[curve]]

  maturities <- maturities %||% grid
  bad <- setdiff(maturities, grid)
  if (length(bad)) {
    stop("Maturities not on the panel's grid: ", paste(bad, collapse = ", "),
         call. = FALSE)
  }

  keep <- match(maturities, grid)
  y_ann <- y_ann[, keep, drop = FALSE]

  if (anyNA(y_ann)) {
    stop("The selected yields contain missing values. ACM requires a complete ",
         "panel; drop incomplete maturities or dates first.", call. = FALSE)
  }
  if (!1 %in% maturities) {
    stop("The one-month maturity is required as the short rate. Include it in ",
         "`maturities`.", call. = FALSE)
  }

  return_maturities <- return_maturities %||%
    maturities[maturities >= 2 & (maturities - 1) %in% maturities]
  if (!length(return_maturities)) {
    stop("No maturity has its predecessor on the grid, so no excess return ",
         "can be formed. ACM needs a dense (ideally monthly) grid.",
         call. = FALSE)
  }

  # Monthly rate units throughout; see R/acm.R.
  y <- y_ann / 12
  p <- -sweep(y, 2L, maturities, "*")
  colnames(p) <- as.character(maturities)

  fac <- acm_factors(y, k = n_factors)
  x <- fac$scores
  rownames(x) <- rownames(y)

  rx <- acm_excess_returns(p, maturities, return_maturities)
  r <- y[, match(1L, maturities)]

  pars <- acm_three_step(x, rx, r)

  # The pricing recursion iterates up to n_max times. Unstable risk-adjusted
  # dynamics make it diverge geometrically and silently -- a sparse set of
  # return maturities can produce yields of order 1e267 with no error raised.
  rho_q <- spectral_radius(pars$phi - pars$lambda1)
  rho_p <- spectral_radius(pars$phi)

  if (rho_q >= 1) {
    warning(
      "Risk-adjusted dynamics are explosive (spectral radius ",
      format(rho_q, digits = 4), " >= 1): fitted yields will diverge as ",
      "maturity grows. This usually means the cross-section of excess returns ",
      "is too sparse or ill-conditioned to identify the prices of risk. Try a ",
      "denser `return_maturities`.",
      call. = FALSE
    )
  }
  if (rho_p >= 1) {
    warning(
      "The factor VAR is non-stationary (spectral radius ",
      format(rho_p, digits = 4), " >= 1): expected short rates will not ",
      "converge and risk-neutral yields at long maturities are unreliable.",
      call. = FALSE
    )
  }

  n_max <- max(maturities)
  rec_p <- acm_recursion(n_max, pars, risk_neutral = FALSE)
  rec_q <- acm_recursion(n_max, pars, risk_neutral = TRUE)

  fitted_ann <- acm_yields(rec_p, x, maturities) * 12
  rn_ann <- acm_yields(rec_q, x, maturities) * 12

  zlb <- zlb_check(rn_ann, maturities)

  structure(
    list(
      pricing = pricing,
      p_dynamics = p_dynamics,
      n_factors = n_factors,
      curve = curve,
      instrument = panel$meta$instrument[panel$meta$curve == curve],
      issuer = panel$meta$issuer[panel$meta$curve == curve],
      dates = panel$dates,
      maturities = maturities,
      return_maturities = return_maturities,
      factors = x,
      pca = fac,
      pars = pars,
      recursion = list(p = rec_p, q = rec_q),
      spectral_radius = c(risk_adjusted = rho_q, real_world = rho_p),
      zlb = zlb,
      observed = y_ann,
      fitted = fitted_ann,
      risk_neutral = rn_ann,
      term_premium = fitted_ann - rn_ann,
      frequency = panel$frequency
    ),
    class = "atsm_fit"
  )
}


#' Check model-implied expected short rates against the lower bound
#'
#' A Gaussian affine model has no lower bound: nothing stops it projecting
#' expected short rates arbitrarily far below zero. Near the effective lower
#' bound that is a real problem rather than a curiosity, because a model that
#' cannot represent the bound will place probability on rates that policy will
#' not deliver, and bias the resulting term premium.
#'
#' Mildly negative values are not automatically wrong -- Polish, euro area and
#' Swiss policy rates were genuinely negative, and the Polish curve in this
#' package contains negative fitted yields. The diagnostic therefore reports
#' rather than forbids, and only warns when the breach is large enough that it
#' cannot be an actual policy rate.
#'
#' @param rn Matrix of annualised risk-neutral yields (decimals).
#' @param maturities Maturities in months.
#' @param threshold Annualised decimal below which a value is treated as
#'   implausible for a realised policy rate. Default -1%.
#'
#' @return A list with `min`, `frac_negative`, `frac_below_threshold` and
#'   `worst_maturity`.
#' @keywords internal
#' @noRd
zlb_check <- function(rn, maturities, threshold = -0.01) {
  worst <- min(rn)
  frac_neg <- mean(rn < 0)
  frac_bad <- mean(rn < threshold)

  if (worst < threshold) {
    j <- which(apply(rn, 2L, min) < threshold)
    warning(
      "Model-implied expected short rates fall to ",
      format(round(worst * 100, 2), nsmall = 2), "% at ",
      maturities[j[1]], "-month maturity, below the ",
      format(threshold * 100), "% plausibility threshold. Gaussian affine ",
      "models impose no lower bound, so estimates near the effective lower ",
      "bound may be unreliable. See Wu and Xia (2016) on shadow-rate models.",
      call. = FALSE
    )
  }

  list(
    min = worst,
    frac_negative = frac_neg,
    frac_below_threshold = frac_bad,
    threshold = threshold,
    worst_maturity = maturities[which.min(apply(rn, 2L, min))]
  )
}


# Extractors --------------------------------------------------------------

#' Extract components of a fitted decomposition
#'
#' @param object An `atsm_fit`.
#' @param maturity Maturities in months. Defaults to every fitted maturity.
#' @param ... Unused.
#'
#' @return A data frame with columns `date`, `maturity` and `value`. Yields are
#'   annualised decimals, so multiply by 100 for percent or 10000 for basis
#'   points.
#'
#' @name atsm-extractors
NULL

#' @rdname atsm-extractors
#' @export
term_premium <- function(object, maturity = NULL, ...) {
  extract_component(object, "term_premium", maturity)
}

#' @rdname atsm-extractors
#' @export
risk_neutral <- function(object, maturity = NULL, ...) {
  extract_component(object, "risk_neutral", maturity)
}

#' @rdname atsm-extractors
#'
#' @details `expected_short_rate()` is the average expected short rate over the
#'   life of an `n`-month bond, which is the risk-neutral yield -- not an
#'   instantaneous forecast of the policy rate.
#' @export
expected_short_rate <- function(object, maturity = NULL, ...) {
  extract_component(object, "risk_neutral", maturity)
}

#' @noRd
extract_component <- function(object, what, maturity = NULL) {
  stopifnot(inherits(object, "atsm_fit"))

  maturity <- maturity %||% object$maturities
  bad <- setdiff(maturity, object$maturities)
  if (length(bad)) {
    stop("Maturity(ies) not fitted: ", paste(bad, collapse = ", "),
         call. = FALSE)
  }

  m <- object[[what]][, match(maturity, object$maturities), drop = FALSE]

  data.frame(
    date = rep(object$dates, times = length(maturity)),
    maturity = rep(maturity, each = length(object$dates)),
    value = as.vector(m),
    stringsAsFactors = FALSE
  )
}

#' @export
fitted.atsm_fit <- function(object, ...) {
  extract_component(object, "fitted", NULL)
}

#' @export
residuals.atsm_fit <- function(object, ...) {
  out <- extract_component(object, "observed", NULL)
  out$value <- out$value - extract_component(object, "fitted", NULL)$value
  out
}

#' @export
coef.atsm_fit <- function(object, ...) {
  object$pars[c("mu", "phi", "sigma", "sigma2",
                "lambda0", "lambda1", "delta0", "delta1")]
}

#' Apply a fitted model to new factor observations
#'
#' Recomputes the decomposition at dates outside the estimation sample, holding
#' the estimated parameters fixed. This is how a model estimated on month-end
#' data is evaluated at a higher frequency: the parameters need monthly excess
#' returns to estimate, but the decomposition itself needs only the factors on
#' the date in question.
#'
#' @param object An `atsm_fit`.
#' @param newdata A [yield_panel] whose curve shares the estimation grid.
#' @param ... Unused.
#'
#' @return A data frame with `date`, `maturity`, `fitted`, `risk_neutral` and
#'   `term_premium`.
#'
#' @export
predict.atsm_fit <- function(object, newdata, ...) {
  if (!inherits(newdata, "yield_panel")) {
    stop("`newdata` must be a yield_panel.", call. = FALSE)
  }

  grid <- newdata$maturities[[newdata$meta$curve[1L]]]
  if (!all(object$maturities %in% grid)) {
    stop("`newdata` is missing maturities the model was fitted on.",
         call. = FALSE)
  }

  y_ann <- curve_matrix(newdata)[, match(object$maturities, grid), drop = FALSE]
  y <- y_ann / 12

  # Project onto the estimation-sample principal components: centre with the
  # training means and rotate with the training loadings, so the factors are
  # on the same scale the parameters were estimated against.
  xc <- sweep(y, 2L, object$pca$center, "-")
  x <- xc %*% object$pca$loadings

  fit <- acm_yields(object$recursion$p, x, object$maturities) * 12
  rn <- acm_yields(object$recursion$q, x, object$maturities) * 12

  data.frame(
    date = rep(newdata$dates, times = length(object$maturities)),
    maturity = rep(object$maturities, each = length(newdata$dates)),
    fitted = as.vector(fit),
    risk_neutral = as.vector(rn),
    term_premium = as.vector(fit - rn),
    stringsAsFactors = FALSE
  )
}


# Display -----------------------------------------------------------------

#' @export
print.atsm_fit <- function(x, ...) {
  resid <- x$observed - x$fitted
  rmse_bp <- sqrt(mean(resid^2)) * 1e4

  cat("<atsm_fit>\n")
  cat("  model      : ", toupper(x$pricing), " (", x$n_factors,
      " factors, P-dynamics: ", x$p_dynamics, ")\n", sep = "")
  cat("  curve      : ", x$curve, " [", x$instrument,
      if (is.na(x$issuer)) "" else paste0("/", x$issuer), "]\n", sep = "")
  cat("  sample     : ", length(x$dates), " ", x$frequency, " obs, ",
      format(min(x$dates)), " to ", format(max(x$dates)), "\n", sep = "")
  cat("  maturities : ", min(x$maturities), "-", max(x$maturities),
      " months (", length(x$maturities), ")\n", sep = "")
  cat("  fit        : ", sprintf("%.2f bp RMSE", rmse_bp), "\n", sep = "")

  long <- max(x$maturities)
  tp <- x$term_premium[, match(long, x$maturities)] * 1e4
  cat("  ", long %/% 12, "y term premium: mean ", sprintf("%.0f bp", mean(tp)),
      ", last ", sprintf("%.0f bp", tp[length(tp)]), "\n", sep = "")

  if (x$zlb$frac_negative > 0) {
    cat("  note       : expected short rates are negative in ",
        sprintf("%.1f%%", 100 * x$zlb$frac_negative),
        " of cells (min ", sprintf("%.2f%%", 100 * x$zlb$min),
        "); no lower bound is imposed\n", sep = "")
  }
  if (x$spectral_radius[["risk_adjusted"]] >= 1) {
    cat("  warning    : risk-adjusted dynamics explosive (rho = ",
        sprintf("%.4f", x$spectral_radius[["risk_adjusted"]]), ")\n", sep = "")
  }

  invisible(x)
}

#' @export
summary.atsm_fit <- function(object, ...) {
  resid <- object$observed - object$fitted

  per_mat <- data.frame(
    maturity = object$maturities,
    rmse_bp = round(sqrt(colMeans(resid^2)) * 1e4, 2),
    mean_tp_bp = round(colMeans(object$term_premium) * 1e4, 1),
    sd_tp_bp = round(apply(object$term_premium, 2L, stats::sd) * 1e4, 1),
    stringsAsFactors = FALSE
  )

  structure(
    list(fit = object, per_maturity = per_mat),
    class = "summary.atsm_fit"
  )
}

#' @export
print.summary.atsm_fit <- function(x, ...) {
  print(x$fit)
  cat("\nBy maturity:\n")

  show <- x$per_maturity
  if (nrow(show) > 12L) {
    keep <- unique(c(1L, which(show$maturity %% 12 == 0)))
    show <- show[keep, , drop = FALSE]
  }
  print(show, row.names = FALSE)
  invisible(x)
}
