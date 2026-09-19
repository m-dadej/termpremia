# User-facing joint real-nominal model. The estimator itself is in
# R/real-nominal.R; this file is the interface, the input plumbing and the
# decomposition arithmetic.

#' Joint real and nominal term structure model
#'
#' Fits the Abrahams, Adrian, Crump, Moench and Yu (2016) affine model to a
#' nominal and an inflation-indexed yield curve at once, and decomposes both
#' into expectations and risk premium components. The breakeven inflation rate
#' implied by the two curves is split into expected inflation, an inflation
#' risk premium and -- when a liquidity factor is supplied -- a liquidity
#' component.
#'
#' @section What you get that `atsm()` cannot give you:
#' [atsm()] fits one curve and splits its yields into expected average future
#' short rates and a term premium. Running it separately on a nominal and a
#' real curve and differencing the results does **not** give a
#' decomposition of breakeven inflation, because the two fits would have
#' unrelated factor spaces and unrelated prices of risk. The joint model
#' prices both curves off one state vector and one stochastic discount factor,
#' which is what makes the inflation risk premium a well-defined object.
#'
#' @section The factor structure:
#' Six factors in the paper's US specification, and the ordering is load
#' bearing: `n_factors_nominal` principal components of the nominal curve,
#' then `n_factors_real` principal components of inflation-indexed yields
#' *after* projecting out the nominal components and the liquidity factor,
#' then the liquidity factor itself. The orthogonalisation is what keeps the
#' joint VAR identified; see Section 3.1 of the paper.
#'
#' @section Liquidity:
#' TIPS trade less liquidly than nominal Treasuries, especially before 2004
#' and after the Lehman failure, and that illiquidity is priced. The paper
#' handles it by putting an observable liquidity index into the state vector.
#' Supply it through `liquidity`, optionally built with
#' [tips_liquidity_factor()].
#'
#' `liquidity` is optional because the paper's own UK specification omits it
#' -- no comparable measure of inflation-indexed gilt liquidity was available
#' to the authors -- and because one of its two US inputs, the average
#' absolute TIPS curve fitting error, is not part of the published
#' Gurkaynak-Sack-Wright file. Without it the model still identifies expected
#' inflation and the inflation risk premium; what it can no longer do is
#' separate a liquidity premium from the inflation risk premium, so the latter
#' absorbs it. That is a real limitation and it is reported in `print()`.
#'
#' @section Sample alignment:
#' A panel holding a 1961-onwards nominal curve and a 1999-onwards TIPS curve
#' has no complete rows before 1999. The fit is therefore restricted to dates
#' where both curves, the short rate and inflation are all observed, and the
#' retained range is reported. An internal gap in that range is an error
#' rather than a warning: the factor VAR treats consecutive rows as
#' consecutive months, so a missing month would silently corrupt every
#' persistence estimate downstream.
#'
#' @param panel A [yield_panel] containing both curves.
#' @param nominal,real Curve names within `panel`. Default to the curve whose
#'   `instrument` metadata is `"tips"` for `real`, and the first other curve
#'   for `nominal`.
#' @param inflation A data frame with `date` and `value` columns giving the
#'   consumer price index used to index the real bonds' payouts. The paper
#'   uses seasonally **unadjusted** CPI-U; see [fred_series()].
#' @param inflation_units `"index"` (the default) treats `value` as a price
#'   level and differences its logarithm; `"log_change"` takes it as
#'   already-differenced monthly log inflation.
#' @param n_factors_nominal,n_factors_real Components from each block. The
#'   paper's US specification uses 3 and 2.
#' @param liquidity Optional data frame with `date` and `value` columns, or a
#'   numeric vector matching the retained dates.
#' @param maturities,real_maturities Maturity grids to use, in months.
#'   `real_maturities` defaults to the real curve's grid from 23 months up:
#'   23 rather than 24 because a one-month holding return on a two-year bond
#'   needs the 23-month point, and the paper's own 24-month floor comes from
#'   the CPI indexation lag distorting shorter tenors.
#' @param return_maturities,real_return_maturities Maturities whose one-month
#'   excess returns enter the estimation. Default to the paper's
#'   `c(6, 12, seq(24, 120, 12))` and `seq(24, 120, 12)`, intersected with
#'   what the grids can support.
#' @param short_rate Optional data frame with `date` and `value` columns. The
#'   paper uses the effective federal funds rate; without it the one-month
#'   fitted nominal yield is used, which for a curve fitted without bills is
#'   an extrapolation. See the `short_rate` discussion in [atsm()].
#' @param short_rate_units Units of `short_rate$value`.
#'
#' @return An object of class `atsm_real_fit`.
#'
#' @references
#' Abrahams, M., T. Adrian, R. K. Crump, E. Moench and R. Yu (2016).
#' "Decomposing real and nominal yield curves." *Journal of Monetary
#' Economics* 84, 182-200.
#'
#' @seealso [atsm()] for a single curve, [expected_inflation()],
#'   [inflation_risk_premium()], [breakeven()], [gsw_tips()]
#'
#' @examples
#' \dontrun{
#' nom <- gsw_monthly
#' nom$curve <- "nominal"
#' tips <- gsw_tips()
#' tips$curve <- "tips"
#'
#' panel <- yield_panel(
#'   rbind(nom[, c("date", "maturity", "yield", "curve")],
#'         tips[, c("date", "maturity", "yield", "curve")]),
#'   curve = "curve", units = "percent", maturity_unit = "months",
#'   instrument = c(nominal = "government", tips = "tips"), issuer = "US"
#' )
#'
#' cpi <- fred_series("CPIAUCNS")
#' fit <- atsm_real(panel, inflation = cpi,
#'                  short_rate = fred_series("DFF", frequency = "monthly"))
#' fit
#' }
#'
#' @export
atsm_real <- function(panel,
                      nominal = NULL,
                      real = NULL,
                      inflation,
                      inflation_units = c("index", "log_change"),
                      n_factors_nominal = 3L,
                      n_factors_real = 2L,
                      liquidity = NULL,
                      maturities = NULL,
                      real_maturities = NULL,
                      return_maturities = NULL,
                      real_return_maturities = NULL,
                      short_rate = NULL,
                      short_rate_units = c("auto", "percent", "decimal")) {
  inflation_units <- match.arg(inflation_units)
  short_rate_units <- match.arg(short_rate_units)

  if (!inherits(panel, "yield_panel")) {
    stop("`panel` must be a yield_panel; see ?yield_panel.", call. = FALSE)
  }
  if (nrow(panel$meta) < 2L) {
    stop("`atsm_real()` needs two curves in `panel`, a nominal one and an ",
         "inflation-indexed one; it found only '", panel$meta$curve[1L],
         "'. Use atsm() for a single curve.", call. = FALSE)
  }

  curves <- resolve_real_nominal_curves(panel, nominal, real)
  nominal <- curves$nominal
  real <- curves$real

  grid_nom <- panel$maturities[[nominal]]
  grid_real <- panel$maturities[[real]]

  maturities <- maturities %||% grid_nom
  real_maturities <- real_maturities %||% grid_real[grid_real >= 23]

  check_on_grid(maturities, grid_nom, nominal)
  check_on_grid(real_maturities, grid_real, real)

  if (!1 %in% maturities) {
    stop("The one-month nominal maturity is required as the short rate ",
         "fallback. Include it in `maturities`.", call. = FALSE)
  }

  y_nom_all <- curve_matrix(panel, nominal)[, match(maturities, grid_nom),
                                            drop = FALSE]
  y_real_all <- curve_matrix(panel, real)[, match(real_maturities, grid_real),
                                          drop = FALSE]

  # --- align the sample ----------------------------------------------------
  infl_all <- resolve_inflation(inflation, inflation_units, panel$dates)
  liq_all <- resolve_liquidity(liquidity, panel$dates)

  ok <- stats::complete.cases(y_nom_all) & stats::complete.cases(y_real_all) &
    !is.na(infl_all)
  if (!is.null(liq_all)) ok <- ok & !is.na(liq_all)

  if (sum(ok) < 60L) {
    stop("Only ", sum(ok), " dates have both curves, inflation and (if ",
         "supplied) liquidity observed. The joint model needs a usable ",
         "sample; check that the curves overlap and that `inflation` covers ",
         "them.", call. = FALSE)
  }

  keep <- which(ok)
  gaps <- which(diff(keep) != 1L)
  if (length(gaps)) {
    stop(
      "The usable sample has ", length(gaps), " internal gap(s), the first ",
      "after ", format(panel$dates[keep[gaps[1L]]]), ". The factor VAR reads ",
      "consecutive rows as consecutive periods, so a gap would corrupt every ",
      "persistence estimate without any error being raised. Drop the ",
      "incomplete maturities or restrict the panel to an unbroken run.",
      call. = FALSE
    )
  }

  dates <- panel$dates[keep]
  y_nom_ann <- y_nom_all[keep, , drop = FALSE]
  y_real_ann <- y_real_all[keep, , drop = FALSE]
  infl <- infl_all[keep]
  liq <- if (is.null(liq_all)) NULL else liq_all[keep]

  # --- monthly units, log prices ------------------------------------------
  y_nom <- y_nom_ann / 12
  y_real <- y_real_ann / 12
  p_nom <- -sweep(y_nom, 2L, maturities, "*")
  p_real <- -sweep(y_real, 2L, real_maturities, "*")

  # --- factors -------------------------------------------------------------
  fac <- acmy_factors(y_nom, y_real, n_factors_nominal, n_factors_real, liq)
  x <- fac$x
  rownames(x) <- format(dates)

  # --- returns -------------------------------------------------------------
  r <- resolve_short_rate(short_rate, short_rate_units, dates,
                          y_nom[, match(1L, maturities)])

  return_maturities <- return_maturities %||%
    intersect(c(6L, 12L, seq(24L, 120L, by = 12L)), maturities)
  real_return_maturities <- real_return_maturities %||%
    intersect(seq(24L, 120L, by = 12L), real_maturities)

  rx_nom <- acm_excess_returns(p_nom, maturities, return_maturities, r)
  rx_real <- acm_excess_returns(p_real, real_maturities,
                                real_return_maturities, r)

  # --- estimate ------------------------------------------------------------
  n_max <- max(maturities, real_maturities)
  pars <- acmy_closed_form(x, rx_nom, rx_real, infl, r,
                           real_return_maturities, n_max)

  rho_q <- spectral_radius(pars$phi_tilde)
  rho_p <- spectral_radius(pars$phi)

  if (rho_q >= 1) {
    warning(
      "Risk-adjusted dynamics are explosive (spectral radius ",
      format(rho_q, digits = 4), " >= 1): fitted yields on both curves will ",
      "diverge as maturity grows. With six factors this usually means the ",
      "joint cross-section of excess returns is too sparse to identify the ",
      "prices of risk.", call. = FALSE
    )
  }
  if (rho_p >= 1) {
    warning(
      "The factor VAR is non-stationary (spectral radius ",
      format(rho_p, digits = 4), " >= 1): expected short rates and expected ",
      "inflation will not converge, so long-horizon expectations components ",
      "are unreliable.", call. = FALSE
    )
  }

  # --- price under both measures ------------------------------------------
  # Q: the pricing measure, giving fitted yields.
  # P: the physical measure, i.e. both prices of risk set to zero. Section 2.1
  #    -- replacing mu-tilde and Phi-tilde with mu and Phi is exactly what
  #    turns the model's breakeven into expected inflation.
  co_q <- acmy_coefficients(n_max, pars$mu_tilde, pars$phi_tilde, pars$sigma,
                            pars$delta0, pars$delta1, pars$pi0, pars$pi1)
  co_p <- acmy_coefficients(n_max, pars$mu, pars$phi, pars$sigma,
                            pars$delta0, pars$delta1, pars$pi0, pars$pi1)

  fit_nom <- acmy_yields(co_q$a, co_q$b, x, maturities) * 12
  fit_real <- acmy_yields(co_q$a_real, co_q$b_real, x, real_maturities) * 12
  rn_nom <- acmy_yields(co_p$a, co_p$b, x, maturities) * 12
  rn_real <- acmy_yields(co_p$a_real, co_p$b_real, x, real_maturities) * 12

  # --- inflation decomposition on the shared maturities -------------------
  shared <- intersect(maturities, real_maturities)
  i_nom <- match(shared, maturities)
  i_real <- match(shared, real_maturities)

  be <- fit_nom[, i_nom, drop = FALSE] - fit_real[, i_real, drop = FALSE]
  ei <- rn_nom[, i_nom, drop = FALSE] - rn_real[, i_real, drop = FALSE]
  irp <- be - ei
  colnames(be) <- colnames(ei) <- colnames(irp) <- as.character(shared)

  liq_nom <- liquidity_component(co_q$b, x, maturities, fac$liq_index)
  liq_real <- liquidity_component(co_q$b_real, x, real_maturities,
                                  fac$liq_index)

  structure(
    list(
      pricing = "acmy",
      n_factors = ncol(x),
      n_factors_nominal = n_factors_nominal,
      n_factors_real = n_factors_real,
      has_liquidity = !is.null(liq),
      curves = c(nominal = nominal, real = real),
      dates = dates,
      maturities = maturities,
      real_maturities = real_maturities,
      shared_maturities = shared,
      return_maturities = return_maturities,
      real_return_maturities = real_return_maturities,
      factors = x,
      pca = fac,
      pars = pars,
      recursion = list(q = co_q, p = co_p),
      spectral_radius = c(risk_adjusted = rho_q, real_world = rho_p),
      inflation = infl,
      short_rate = r,
      liquidity = liq,
      observed = set_dimnames(y_nom_ann, dates, maturities),
      observed_real = set_dimnames(y_real_ann, dates, real_maturities),
      fitted = fit_nom,
      fitted_real = fit_real,
      risk_neutral = rn_nom,
      risk_neutral_real = rn_real,
      term_premium = fit_nom - rn_nom,
      term_premium_real = fit_real - rn_real,
      breakeven = be,
      expected_inflation = ei,
      inflation_risk_premium = irp,
      liquidity_nominal = liq_nom,
      liquidity_real = liq_real,
      frequency = panel$frequency
    ),
    class = "atsm_real_fit"
  )
}


# Input plumbing ----------------------------------------------------------

#' Work out which curve is nominal and which is inflation-indexed
#'
#' The panel already records an `instrument` per curve, so `"tips"` identifies
#' the real curve without the user having to name it. Guessing is confined to
#' this function and it refuses to guess when the metadata is ambiguous.
#'
#' @keywords internal
#' @noRd
resolve_real_nominal_curves <- function(panel, nominal, real) {
  meta <- panel$meta

  if (is.null(real)) {
    tips <- meta$curve[meta$instrument == "tips"]
    if (length(tips) != 1L) {
      stop(
        "Could not tell which curve is inflation-indexed: ",
        if (!length(tips)) {
          "no curve has instrument = \"tips\""
        } else {
          paste0(length(tips), " curves do (",
                 paste(tips, collapse = ", "), ")")
        },
        ". Name it with `real = `, or set the instrument metadata in ",
        "yield_panel().", call. = FALSE
      )
    }
    real <- tips
  }
  if (is.null(nominal)) {
    others <- setdiff(meta$curve, real)
    if (!length(others)) {
      stop("`panel` has no curve other than '", real, "'.", call. = FALSE)
    }
    nominal <- others[1L]
  }

  for (nm in c(nominal, real)) {
    if (!nm %in% meta$curve) {
      stop("No curve named '", nm, "' in `panel`. Available: ",
           paste(meta$curve, collapse = ", "), call. = FALSE)
    }
  }
  if (identical(nominal, real)) {
    stop("`nominal` and `real` must be different curves.", call. = FALSE)
  }

  list(nominal = nominal, real = real)
}

#' @keywords internal
#' @noRd
set_dimnames <- function(m, dates, maturities) {
  dimnames(m) <- list(format(dates), as.character(maturities))
  m
}

#' @keywords internal
#' @noRd
check_on_grid <- function(maturities, grid, what) {
  bad <- setdiff(maturities, grid)
  if (length(bad)) {
    stop("Maturities not on the '", what, "' grid: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

#' One-period log inflation aligned to the panel dates
#'
#' Matching is by calendar month, not by exact date. A price index is dated to
#' the first of the month while a month-end yield curve is dated to the last
#' business day, so an exact join would find nothing at all.
#'
#' @keywords internal
#' @noRd
resolve_inflation <- function(inflation, units, dates) {
  if (missing(inflation) || is.null(inflation)) {
    stop("`inflation` is required: the real bonds' payouts are indexed to a ",
         "price level, so the model cannot be identified without one. ",
         "Seasonally unadjusted CPI-U is what the paper uses; see ",
         "?fred_series.", call. = FALSE)
  }
  if (!is.data.frame(inflation) ||
      !all(c("date", "value") %in% names(inflation))) {
    stop("`inflation` must be a data frame with `date` and `value` columns.",
         call. = FALSE)
  }

  idate <- as_date_strict(inflation$date, arg = "inflation$date")
  value <- as.numeric(inflation$value)

  o <- order(idate)
  idate <- idate[o]
  value <- value[o]

  ym <- format(idate, "%Y-%m")
  target <- format(dates, "%Y-%m")

  if (units == "log_change") {
    out <- value[match(target, ym)]
    return(out)
  }

  # Price level: the monthly log change into the panel date's own month.
  prev <- format(add_months(dates, -1L), "%Y-%m")
  now_idx <- match(target, ym)
  prev_idx <- match(prev, ym)

  out <- rep(NA_real_, length(dates))
  have <- !is.na(now_idx) & !is.na(prev_idx)
  if (any(value[stats::na.omit(c(now_idx, prev_idx))] <= 0, na.rm = TRUE)) {
    stop("`inflation` contains non-positive values, so it cannot be a price ",
         "index. Pass `inflation_units = \"log_change\"` if it is already a ",
         "rate of change.", call. = FALSE)
  }
  out[have] <- log(value[now_idx[have]] / value[prev_idx[have]])
  out
}

#' @keywords internal
#' @noRd
resolve_liquidity <- function(liquidity, dates) {
  if (is.null(liquidity)) return(NULL)

  if (is.numeric(liquidity)) {
    if (length(liquidity) != length(dates)) {
      stop("`liquidity` has ", length(liquidity), " values but the panel has ",
           length(dates), " dates. Supply a data frame with `date` and ",
           "`value` columns to have it matched instead.", call. = FALSE)
    }
    return(as.numeric(liquidity))
  }
  if (!is.data.frame(liquidity) ||
      !all(c("date", "value") %in% names(liquidity))) {
    stop("`liquidity` must be a numeric vector or a data frame with `date` ",
         "and `value` columns.", call. = FALSE)
  }

  ldate <- as_date_strict(liquidity$date, arg = "liquidity$date")
  ym <- format(ldate, "%Y-%m")
  as.numeric(liquidity$value)[match(format(dates, "%Y-%m"), ym)]
}

#' The part of each fitted yield attributable to the liquidity factor
#'
#' Section 2.2 reports liquidity-adjusted inflation measures by stripping the
#' liquidity factor's contribution out of both curves before differencing
#' them. That contribution is just the factor's own term in the affine
#' expression for the yield.
#'
#' @keywords internal
#' @noRd
liquidity_component <- function(b, x, maturities, liq_index) {
  if (is.na(liq_index)) return(NULL)

  out <- matrix(NA_real_, nrow = nrow(x), ncol = length(maturities),
                dimnames = list(rownames(x), as.character(maturities)))
  for (j in seq_along(maturities)) {
    n <- maturities[j]
    out[, j] <- -(b[n, liq_index] * x[, liq_index]) / n * 12
  }
  out
}


# Extractors --------------------------------------------------------------

#' Extract components of a joint real-nominal decomposition
#'
#' @section The decomposition:
#' For a maturity `n` the model splits the two curves and their difference as
#'
#' \deqn{y^{(n)} = \underbrace{E[\bar r]}_{\mathrm{risk\ neutral}} +
#'   \mathrm{TP}^{(n)}, \qquad
#'   y^{(n)}_R = E[\bar r_R] + \mathrm{TP}^{(n)}_R}
#' \deqn{\underbrace{y^{(n)} - y^{(n)}_R}_{\mathrm{breakeven}} =
#'   \underbrace{E[\bar\pi]}_{\mathrm{expected\ inflation}} +
#'   \underbrace{\mathrm{TP}^{(n)} - \mathrm{TP}^{(n)}_R}_{
#'     \mathrm{inflation\ risk\ premium}}}
#'
#' Both identities hold exactly, by construction rather than approximately,
#' and the tests assert them. The inflation risk premium being the *difference*
#' of the two term premia is not an extra assumption: it follows from
#' expected inflation being the breakeven computed under the physical measure.
#'
#' @section Liquidity adjustment:
#' With `liquidity_adjusted = TRUE` the liquidity factor's own contribution is
#' stripped from both curves before they are differenced, which is how
#' Section 2.2 of the paper reports its inflation measures. It is an error to
#' ask for it when the model was fitted without a liquidity factor, rather
#' than silently returning the unadjusted number.
#'
#' @param object An `atsm_real_fit` from [atsm_real()].
#' @param maturity Maturities in months. Defaults to every maturity the
#'   component is defined on -- the nominal grid for nominal quantities, the
#'   real grid for real ones, and their intersection for inflation ones.
#' @param liquidity_adjusted Strip the liquidity factor's contribution first.
#' @param ... Unused.
#'
#' @return A data frame with columns `date`, `maturity` and `value`, in
#'   annualised decimals.
#'
#' @name atsm-real-extractors
NULL

#' @rdname atsm-real-extractors
#' @export
breakeven <- function(object, maturity = NULL, liquidity_adjusted = FALSE,
                      ...) {
  extract_real_component(object, "breakeven", maturity, liquidity_adjusted)
}

#' @rdname atsm-real-extractors
#' @export
expected_inflation <- function(object, maturity = NULL,
                               liquidity_adjusted = FALSE, ...) {
  extract_real_component(object, "expected_inflation", maturity,
                         liquidity_adjusted)
}

#' @rdname atsm-real-extractors
#' @export
inflation_risk_premium <- function(object, maturity = NULL,
                                   liquidity_adjusted = FALSE, ...) {
  extract_real_component(object, "inflation_risk_premium", maturity,
                         liquidity_adjusted)
}

#' @rdname atsm-real-extractors
#' @export
real_term_premium <- function(object, maturity = NULL, ...) {
  extract_real_component(object, "term_premium_real", maturity)
}

#' @rdname atsm-real-extractors
#' @export
real_risk_neutral <- function(object, maturity = NULL, ...) {
  extract_real_component(object, "risk_neutral_real", maturity)
}

#' @rdname atsm-real-extractors
#' @export
liquidity_premium <- function(object, maturity = NULL, ...) {
  extract_real_component(object, "liquidity_real", maturity)
}

#' @export
term_premium.atsm_real_fit <- function(object, maturity = NULL, ...) {
  extract_real_component(object, "term_premium", maturity)
}

#' @export
risk_neutral.atsm_real_fit <- function(object, maturity = NULL, ...) {
  extract_real_component(object, "risk_neutral", maturity)
}

#' @export
expected_short_rate.atsm_real_fit <- function(object, maturity = NULL, ...) {
  extract_real_component(object, "risk_neutral", maturity)
}

#' Which maturity grid a component of the joint fit lives on
#'
#' Three grids coexist in one fit and mixing them up would silently
#' misalign a series by a couple of years: nominal quantities run over the
#' whole nominal grid, real ones start where inflation-indexed data does, and
#' the inflation decomposition exists only where the two overlap.
#'
#' @keywords internal
#' @noRd
real_component_grid <- function(object, what) {
  switch(
    what,
    observed = ,
    fitted = ,
    risk_neutral = ,
    term_premium = ,
    liquidity_nominal = object$maturities,
    observed_real = ,
    fitted_real = ,
    risk_neutral_real = ,
    term_premium_real = ,
    liquidity_real = object$real_maturities,
    breakeven = ,
    expected_inflation = ,
    inflation_risk_premium = object$shared_maturities,
    stop("Unknown component '", what, "'.", call. = FALSE)
  )
}

#' @keywords internal
#' @noRd
extract_real_component <- function(object, what, maturity = NULL,
                                   liquidity_adjusted = FALSE) {
  if (!inherits(object, "atsm_real_fit")) {
    stop("`object` must be an atsm_real_fit; see ?atsm_real.", call. = FALSE)
  }

  value <- object[[what]]
  if (is.null(value)) {
    stop("This fit has no '", what, "' component. ",
         if (grepl("liquidity", what)) {
           "It was fitted without a liquidity factor; see ?atsm_real."
         } else "",
         call. = FALSE)
  }

  if (liquidity_adjusted) {
    if (!object$has_liquidity) {
      stop("`liquidity_adjusted = TRUE` needs a model fitted with a ",
           "`liquidity` factor. This one was not, so there is no liquidity ",
           "component to remove and the inflation risk premium already ",
           "absorbs it.", call. = FALSE)
    }
    value <- value - liquidity_adjustment(object, what)
  }

  grid <- real_component_grid(object, what)
  maturity <- maturity %||% grid
  bad <- setdiff(maturity, grid)
  if (length(bad)) {
    stop("Maturity(ies) not available for '", what, "': ",
         paste(bad, collapse = ", "), ". Available: ", min(grid), "-",
         max(grid), " months.", call. = FALSE)
  }

  m <- value[, match(maturity, grid), drop = FALSE]

  data.frame(
    date = rep(object$dates, times = length(maturity)),
    maturity = rep(maturity, each = length(object$dates)),
    value = as.vector(m),
    stringsAsFactors = FALSE
  )
}

#' The liquidity contribution to be removed from an inflation component
#'
#' Breakeven inflation carries the liquidity loading of both curves, so the
#' adjustment is the difference of the two. Expected inflation is a
#' physical-measure quantity but the liquidity factor is in the state vector
#' either way, so the same difference applies; what is left over goes to the
#' risk premium.
#'
#' @keywords internal
#' @noRd
liquidity_adjustment <- function(object, what) {
  shared <- object$shared_maturities
  ln <- object$liquidity_nominal[, match(shared, object$maturities),
                                 drop = FALSE]
  lr <- object$liquidity_real[, match(shared, object$real_maturities),
                              drop = FALSE]
  diff_l <- ln - lr

  switch(
    what,
    breakeven = diff_l,
    expected_inflation = matrix(0, nrow(diff_l), ncol(diff_l)),
    inflation_risk_premium = diff_l,
    stop("'", what, "' has no liquidity adjustment defined.", call. = FALSE)
  )
}


# Display -----------------------------------------------------------------

#' @export
print.atsm_real_fit <- function(x, ...) {
  rmse_n <- sqrt(mean((x$observed - x$fitted)^2)) * 1e4
  rmse_r <- sqrt(mean((x$observed_real - x$fitted_real)^2)) * 1e4

  cat("<atsm_real_fit>\n")
  cat("  model      : ACMY joint real-nominal (", x$n_factors, " factors: ",
      x$n_factors_nominal, " nominal + ", x$n_factors_real, " real",
      if (x$has_liquidity) " + liquidity" else "", ")\n", sep = "")
  cat("  curves     : ", x$curves[["nominal"]], " (nominal) + ",
      x$curves[["real"]], " (real)\n", sep = "")
  cat("  sample     : ", length(x$dates), " ", x$frequency, " obs, ",
      format(min(x$dates)), " to ", format(max(x$dates)), "\n", sep = "")
  cat("  maturities : nominal ", min(x$maturities), "-", max(x$maturities),
      "m, real ", min(x$real_maturities), "-", max(x$real_maturities),
      "m\n", sep = "")
  cat("  fit        : ", sprintf("%.2f bp", rmse_n), " nominal, ",
      sprintf("%.2f bp", rmse_r), " real (RMSE)\n", sep = "")

  long <- max(x$shared_maturities)
  j <- match(long, x$shared_maturities)
  jn <- match(long, x$maturities)
  jr <- match(long, x$real_maturities)
  lab <- paste0(long %/% 12, "y")

  cat("\n  ", lab, " decomposition, mean (sd) in bp:\n", sep = "")
  row <- function(name, v) {
    cat(sprintf("    %-26s %8.0f (%5.0f)\n", name,
                mean(v) * 1e4, stats::sd(v) * 1e4))
  }
  row("nominal term premium", x$term_premium[, jn])
  row("real term premium", x$term_premium_real[, jr])
  row("breakeven inflation", x$breakeven[, j])
  row("  expected inflation", x$expected_inflation[, j])
  row("  inflation risk premium", x$inflation_risk_premium[, j])
  if (x$has_liquidity) {
    row("  TIPS liquidity premium", x$liquidity_real[, jr])
  }

  if (!x$has_liquidity) {
    cat("\n  note       : fitted without a liquidity factor, so any TIPS ",
        "liquidity\n               premium is absorbed into the inflation ",
        "risk premium\n", sep = "")
  }
  if (!x$pars$inflation_fit$converged) {
    cat("  warning    : the inflation loadings did not converge (optim code ",
        x$pars$inflation_fit$convergence, ")\n", sep = "")
  }
  if (x$spectral_radius[["risk_adjusted"]] >= 1) {
    cat("  warning    : risk-adjusted dynamics explosive (rho = ",
        sprintf("%.4f", x$spectral_radius[["risk_adjusted"]]), ")\n", sep = "")
  }
  if (x$spectral_radius[["real_world"]] >= 1) {
    cat("  warning    : factor VAR non-stationary (rho = ",
        sprintf("%.4f", x$spectral_radius[["real_world"]]), ")\n", sep = "")
  }

  invisible(x)
}

#' @export
summary.atsm_real_fit <- function(object, ...) {
  shared <- object$shared_maturities
  jn <- match(shared, object$maturities)
  jr <- match(shared, object$real_maturities)

  per_mat <- data.frame(
    maturity = shared,
    rmse_nom_bp = round(sqrt(colMeans(
      (object$observed[, jn, drop = FALSE] -
         object$fitted[, jn, drop = FALSE])^2)) * 1e4, 2),
    rmse_real_bp = round(sqrt(colMeans(
      (object$observed_real[, jr, drop = FALSE] -
         object$fitted_real[, jr, drop = FALSE])^2)) * 1e4, 2),
    tp_nom_bp = round(colMeans(object$term_premium[, jn, drop = FALSE]) * 1e4, 0),
    tp_real_bp = round(colMeans(
      object$term_premium_real[, jr, drop = FALSE]) * 1e4, 0),
    breakeven_bp = round(colMeans(object$breakeven) * 1e4, 0),
    exp_infl_bp = round(colMeans(object$expected_inflation) * 1e4, 0),
    irp_bp = round(colMeans(object$inflation_risk_premium) * 1e4, 0),
    stringsAsFactors = FALSE
  )

  structure(list(fit = object, per_maturity = per_mat),
            class = "summary.atsm_real_fit")
}

#' @export
print.summary.atsm_real_fit <- function(x, ...) {
  print(x$fit)
  cat("\nBy maturity (means, bp):\n")

  show <- x$per_maturity
  if (nrow(show) > 12L) {
    keep <- which(show$maturity %% 12 == 0)
    if (length(keep)) show <- show[keep, , drop = FALSE]
  }
  print(show, row.names = FALSE)
  invisible(x)
}

#' @export
coef.atsm_real_fit <- function(object, ...) {
  object$pars[c("mu", "phi", "sigma", "mu_tilde", "phi_tilde",
                "lambda0", "lambda1", "delta0", "delta1", "pi0", "pi1")]
}

#' @export
fitted.atsm_real_fit <- function(object, ...) {
  extract_real_component(object, "fitted", NULL)
}

#' @export
residuals.atsm_real_fit <- function(object, ...) {
  out <- extract_real_component(object, "observed", NULL)
  out$value <- out$value - extract_real_component(object, "fitted", NULL)$value
  out
}
