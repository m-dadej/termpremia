# The Abrahams, Adrian, Crump, Moench and Yu (2016) joint real-nominal model.
#
# Equation numbers refer to "Decomposing real and nominal yield curves",
# Journal of Monetary Economics 84, 182-200, and to its Supplementary
# Appendix where marked (S).
#
# Units. As in R/acm.R, everything here is in MONTHLY rate units: an
# annualised decimal yield of 0.048 enters as 0.004, log prices are
# p(n) = -n * y(n) with n in months, and one-period log inflation is a monthly
# log change. Results are re-annualised on the way out.
#
# Relationship to R/acm.R. This is a different estimator, not a wrapper. Two
# differences are easy to miss and both change the numbers:
#
#   1. There is no `+ sigma2 / 2` term. ACM (2013) carry one because they
#      allow a maturity-specific return pricing error that is folded into the
#      no-arbitrage recursion (see R/acm.R). ACMY instead put additive
#      Gaussian errors on returns in the likelihood, and their published
#      recursions (Eq. 4-5, 10-11) have no such term. Borrowing
#      `acm_recursion()` here would give plausible but wrong premia.
#
#   2. The recursion is seeded at A_0 = B_0 = 0 rather than ACM's A_1, B_1.
#      With B_0 = 0 the two agree at n = 1 anyway, because the quadratic term
#      vanishes; they diverge from n = 2 through the sigma2 term above.


# Pricing -----------------------------------------------------------------

#' Joint pricing recursions for nominal and inflation-indexed bonds
#'
#' Nominal bonds (Eq. 4-5):
#' \deqn{A_n = A_{n-1} + B_{n-1}'\tilde\mu +
#'   \tfrac12 B_{n-1}'\Sigma B_{n-1} - \delta_0, \quad
#'   B_n' = B_{n-1}'\tilde\Phi - \delta_1'}
#'
#' Inflation-indexed bonds (Eq. 10-11) are the same recursion with two
#' changes, both of which come from the fact that the bond pays
#' \eqn{Q_{t+n}/Q_t} rather than 1:
#' \deqn{A_{n,R} = A_{n-1,R} + \tilde B_{n-1,R}'\tilde\mu +
#'   \tfrac12 \tilde B_{n-1,R}'\Sigma\tilde B_{n-1,R} - (\delta_0 - \pi_0),
#'   \quad B_{n,R}' = \tilde B_{n-1,R}'\tilde\Phi - \delta_1'}
#'
#' where \eqn{\tilde B_{n,R} = B_{n,R} + \pi_1} throughout. One period of
#' inflation \eqn{\pi_{t+1} = \pi_0 + \pi_1'X_{t+1}} enters the payoff, so it
#' shifts the loading that gets carried into the next step and adds its own
#' intercept; it does **not** appear as a separate term in `B_{n,R}`, because
#' it reaches the price only through \eqn{X_{t+1}}.
#'
#' A useful check on the signs: at `n = 1` the real yield comes out as
#' \eqn{r_t - E^Q_t[\pi_{t+1}] - \tfrac12\pi_1'\Sigma\pi_1}, the one-period
#' nominal rate less risk-neutral expected inflation less a Jensen term. This
#' is asserted in the tests.
#'
#' Both curves are returned from one call because they share `mu_adj`,
#' `phi_adj` and `sigma`, and because keeping them side by side is the only
#' way the `pi` terms stay auditable.
#'
#' @param n_max Longest maturity, in months.
#' @param mu_adj,phi_adj Intercept and autoregressive matrix of the measure
#'   being priced under: \eqn{\tilde\mu = \mu - \lambda_0} and
#'   \eqn{\tilde\Phi = \Phi - \lambda_1} for the pricing measure, or plain
#'   `mu` and `phi` for the physical measure, which is what turns the
#'   breakeven into expected inflation (Section 2.1).
#' @param sigma Innovation covariance, `k x k`.
#' @param delta0,delta1 Nominal short rate loadings,
#'   \eqn{r_t = \delta_0 + \delta_1'X_t}.
#' @param pi0,pi1 One-period log inflation loadings,
#'   \eqn{\pi_t = \pi_0 + \pi_1'X_t}.
#'
#' @return A list with `a`, `b` (nominal) and `a_real`, `b_real`, each
#'   indexed by maturity in months from 1 to `n_max`.
#' @keywords internal
#' @noRd
acmy_coefficients <- function(n_max, mu_adj, phi_adj, sigma,
                              delta0, delta1, pi0, pi1) {
  k <- length(mu_adj)

  a <- numeric(n_max)
  b <- matrix(0, nrow = n_max, ncol = k)
  a_real <- numeric(n_max)
  b_real <- matrix(0, nrow = n_max, ncol = k)

  # A_0 = B_0 = 0; the loop below starts from those rather than from n = 1, so
  # the two curves are written once.
  a_prev <- 0
  b_prev <- numeric(k)
  ar_prev <- 0
  br_prev <- numeric(k)

  for (n in seq_len(n_max)) {
    a[n] <- a_prev + sum(b_prev * mu_adj) +
      0.5 * drop(crossprod(b_prev, sigma %*% b_prev)) - delta0
    b[n, ] <- drop(b_prev %*% phi_adj) - delta1

    # The only two differences from the nominal block: the loading carried
    # forward is B + pi1, and the intercept is delta0 - pi0.
    bt <- br_prev + pi1
    a_real[n] <- ar_prev + sum(bt * mu_adj) +
      0.5 * drop(crossprod(bt, sigma %*% bt)) - delta0 + pi0
    b_real[n, ] <- drop(bt %*% phi_adj) - delta1

    a_prev <- a[n]
    b_prev <- b[n, ]
    ar_prev <- a_real[n]
    br_prev <- b_real[n, ]
  }

  list(a = a, b = b, a_real = a_real, b_real = b_real)
}

#' Model-implied excess returns on inflation-indexed bonds
#'
#' Eq. 18: \eqn{rx^{(n-1,R)}_{t+1} = \alpha_{n-1,R} -
#'   \tilde B_{n-1,R}'\tilde\Phi X_t + B_{n-1,R}'X_{t+1}} with
#' \eqn{\alpha_{n-1,R} = -[\pi_0 + \tilde B_{n-1,R}'\tilde\mu +
#'   \tfrac12\tilde B_{n-1,R}'\Sigma\tilde B_{n-1,R}]}.
#'
#' This is the `g()` of Eq. S31, the criterion used to pin down `pi0` and
#' `pi1`: real bond returns are linear in `pi0` and quadratic in `pi1`, so the
#' cross-section and time series of TIPS returns identify both.
#'
#' @param coefs Output of `acmy_coefficients()` under the pricing measure.
#' @param maturities Maturities `n` whose one-month returns are wanted.
#' @param x_lag,x_led `(T-1) x k` matrices of \eqn{X_t} and \eqn{X_{t+1}}.
#' @param mu_adj,phi_adj,sigma,pi0,pi1 As in `acmy_coefficients()`.
#'
#' @return A `(T-1) x length(maturities)` matrix.
#' @keywords internal
#' @noRd
acmy_real_returns <- function(coefs, maturities, x_lag, x_led,
                              mu_adj, phi_adj, sigma, pi0, pi1) {
  out <- matrix(NA_real_, nrow = nrow(x_lag), ncol = length(maturities),
                dimnames = list(rownames(x_lag), as.character(maturities)))

  for (j in seq_along(maturities)) {
    n <- maturities[j]
    b_prev <- coefs$b_real[n - 1L, ]
    bt <- b_prev + pi1

    alpha <- -(pi0 + sum(bt * mu_adj) +
                 0.5 * drop(crossprod(bt, sigma %*% bt)))

    # The coefficient on X_t is Phi~' B~, which as a row-vector product is
    # `bt %*% phi_adj` and NOT `phi_adj %*% bt`. The two differ whenever
    # Phi~ is asymmetric, which it always is, and for a near-diagonal Phi~
    # the difference is small enough to pass for rounding error: it shows up
    # as a violation of Eq. 18 at 1e-5 rather than 1e-16.
    coef_lag <- drop(bt %*% phi_adj)

    out[, j] <- alpha -
      drop(x_lag %*% coef_lag) +
      drop(x_led %*% b_prev)
  }
  out
}

#' Yields from recursion coefficients
#'
#' @param a,b Intercepts and loadings indexed by maturity.
#' @param x A `T x k` factor matrix.
#' @param maturities Maturities in months.
#' @return A `T x length(maturities)` matrix of yields in monthly units.
#' @keywords internal
#' @noRd
acmy_yields <- function(a, b, x, maturities) {
  out <- matrix(NA_real_, nrow = nrow(x), ncol = length(maturities),
                dimnames = list(rownames(x), as.character(maturities)))
  for (j in seq_along(maturities)) {
    n <- maturities[j]
    out[, j] <- -(a[n] + drop(x %*% b[n, ])) / n
  }
  out
}


# Factors -----------------------------------------------------------------

#' Pricing factors for the joint model
#'
#' Section 3.1. Three blocks, in this order:
#'
#' 1. `k_nominal` principal components of nominal yields.
#' 2. `k_real` principal components of the **residuals** from regressing each
#'    inflation-indexed yield on a constant, the nominal components and (when
#'    supplied) the liquidity factor.
#' 3. The liquidity factor itself, if supplied.
#'
#' The orthogonalisation in step 2 is not cosmetic. Principal components of
#' TIPS yields taken directly are strongly collinear with the nominal
#' components -- both curves are dominated by a level factor -- and a VAR on
#' the union is close to unidentified. Regressing first leaves the real
#' components spanning only what the nominal curve and liquidity do not
#' explain, which is what makes six factors estimable at all.
#'
#' The regression coefficients and the loadings are returned because
#' reproducing the factors out of sample needs them, and because the
#' model-implied-factor constraints of Eq. S(page 2) are functions of exactly
#' these quantities.
#'
#' @param y_nom,y_real `T x N` matrices of nominal and real yields in monthly
#'   units, dates in rows.
#' @param k_nominal,k_real Number of components from each block.
#' @param liquidity Optional length-`T` observable liquidity factor. The US
#'   specification of the paper includes one; the UK specification in the
#'   Supplementary Appendix does not, because no comparable measure of
#'   inflation-indexed gilt liquidity was available.
#'
#' @return A list with `x` (`T x k`), the per-block pieces, and the index of
#'   the liquidity factor within `x` (or `NA`).
#' @keywords internal
#' @noRd
acmy_factors <- function(y_nom, y_real, k_nominal, k_real, liquidity = NULL) {
  if (anyNA(y_nom) || anyNA(y_real)) {
    stop("Yield matrices contain missing values; the joint model needs a ",
         "complete panel on both curves.", call. = FALSE)
  }

  nominal <- acm_factors(y_nom, k = k_nominal)
  x_nom <- nominal$scores

  z <- cbind(1, x_nom)
  if (!is.null(liquidity)) z <- cbind(z, liquidity)

  ortho_coef <- qr.solve(z, y_real)
  resid_real <- y_real - z %*% ortho_coef

  real <- acm_factors(resid_real, k = k_real)

  x <- cbind(x_nom, real$scores)
  colnames(x) <- c(sprintf("N%d", seq_len(k_nominal)),
                   sprintf("R%d", seq_len(k_real)))

  liq_index <- NA_integer_
  if (!is.null(liquidity)) {
    x <- cbind(x, liquidity = liquidity)
    liq_index <- ncol(x)
  }

  # Normalise each factor to unit standard deviation.
  #
  # This is not cosmetic and it is not optional. In monthly rate units the
  # principal components of yields are of order 1e-5, while an observable
  # liquidity index is of order 1. The GLS steps of the estimator solve
  # `B' Sigma_e^-1 B`, a k x k system whose entries then span ten orders of
  # magnitude, and `solve()` declares it singular -- not as a warning, as an
  # error, and only once a liquidity factor is present.
  #
  # Scaling is safe because an affine change of state variables leaves a
  # Gaussian affine model unchanged: rescaling X by D rescales delta1, pi1
  # and every B by D in the opposite direction, and every yield, premium and
  # decomposition comes out identical. The tests assert that invariance.
  #
  # Only the scale is touched, not the location. Centring would destroy the
  # non-negativity that `tips_liquidity_factor()` deliberately imposes so
  # that illiquidity can only ever raise a yield.
  scale <- apply(x, 2L, stats::sd)
  scale[!is.finite(scale) | scale <= 0] <- 1
  x_scaled <- sweep(x, 2L, scale, "/")
  dimnames(x_scaled) <- dimnames(x)

  list(
    x = x_scaled,
    x_unscaled = x,
    scale = scale,
    nominal = nominal,
    real = real,
    ortho_coef = ortho_coef,
    residual_real = resid_real,
    liq_index = liq_index,
    k_nominal = k_nominal,
    k_real = k_real
  )
}


# Estimation --------------------------------------------------------------

#' Closed-form estimator for the joint real-nominal model
#'
#' Supplementary Appendix Section 1.1. Stacking nominal excess returns (Eq. 15)
#' with inflation-augmented real excess returns (Eq. 23) gives one system with
#' a common coefficient matrix (Eq. S25-S28),
#' \deqn{R_{t+1} = \alpha + \bar B(X_{t+1} - \tilde\Phi X_t) +
#'   \varepsilon_{t+1}}
#'
#' which is estimated in four passes: an unconstrained OLS to get
#' \eqn{\hat\Sigma_e}, a GLS solve for \eqn{\tilde\Phi} (Eq. S29), a
#' constrained regression for \eqn{\alpha} and \eqn{\bar B}, and a GLS solve
#' for \eqn{\tilde\mu} (Eq. S30). The inflation loadings then come from
#' minimising real return pricing errors (Eq. S31).
#'
#' Real returns enter **augmented by realised inflation**. That is what makes
#' the stacked system have one common `B`: \eqn{rx^{R} + \pi_{t+1}} has
#' exactly the functional form of a nominal return with loading
#' \eqn{\tilde B_{n,R}} (Eq. 23-24), whereas \eqn{rx^{R}} alone does not.
#' Feeding unaugmented real returns in here is silent and wrong.
#'
#' The prices of risk come out as residuals of the two measures,
#' \eqn{\lambda_0 = \mu - \tilde\mu} and \eqn{\lambda_1 = \Phi - \tilde\Phi},
#' which is the same decomposition `adopt_p_dynamics()` relies on in
#' R/acm.R, read in the other direction.
#'
#' @param x A `T x k` factor matrix.
#' @param rx_nom `(T-1) x N_N` nominal excess returns.
#' @param rx_real `(T-1) x N_R` real excess returns, **not** augmented; the
#'   inflation is added here so that the caller cannot forget.
#' @param inflation Length-`T` one-period log inflation, monthly units.
#' @param r Length-`T` one-period nominal short rate, monthly units.
#' @param real_maturities Maturities of the real returns, in months.
#' @param n_max Longest maturity needed by the recursion.
#'
#' @return A list of parameters and diagnostics.
#' @keywords internal
#' @noRd
acmy_closed_form <- function(x, rx_nom, rx_real, inflation, r,
                             real_maturities, n_max, fix_pi0 = NULL) {
  tt <- nrow(x)
  k <- ncol(x)

  x_lag <- x[-tt, , drop = FALSE]
  x_led <- x[-1L, , drop = FALSE]
  n_obs <- nrow(x_lag)

  # --- P-dynamics: the factor VAR ------------------------------------------
  z1 <- cbind(1, x_lag)
  var_coef <- qr.solve(z1, x_led)
  mu <- var_coef[1L, ]
  phi <- t(var_coef[-1L, , drop = FALSE])
  v <- x_led - z1 %*% var_coef
  sigma <- crossprod(v) / n_obs

  # --- the stacked return system -------------------------------------------
  # Real returns are augmented by next period's realised inflation; see above.
  rx_real_aug <- rx_real + inflation[-1L]
  rr <- cbind(rx_nom, rx_real_aug)
  n_nom <- ncol(rx_nom)
  n_real <- ncol(rx_real)

  # Pass 1: unconstrained, purely to get Sigma_e.
  z2 <- cbind(1, x_lag, x_led)
  coef2 <- qr.solve(z2, rr)
  e_ols <- rr - z2 %*% coef2
  sigma_e <- crossprod(e_ols) / n_obs

  b_ols <- t(coef2[(k + 2L):(2L * k + 1L), , drop = FALSE])   # N x k, = Bbar
  c_ols <- t(coef2[2L:(k + 1L), , drop = FALSE])              # N x k, = -Bbar Phi~

  # Every GLS step below solves `B' Sigma_e^-1 B`, which needs B to have full
  # column rank: each factor must move bond returns somehow, or its price of
  # risk is not identified. LAPACK's report of this is "singular matrix 'a' in
  # solve" from three calls deeper, which says nothing about which factor is
  # at fault, so it is caught here instead.
  #
  # The usual cause is a supplied `liquidity` series that the bonds do not
  # actually price -- an irrelevant or mis-aligned index.
  check_return_loadings(b_ols, k, colnames(x))

  # Pass 2: GLS for Phi-tilde (Eq. S29), as a whitened least-squares problem
  # rather than by forming `B' Sigma_e^-1 B`.
  #
  # The Gram matrix squares the condition number, and Sigma_e is the one
  # matrix here that is genuinely close to singular: twenty return series are
  # driven by six factors, so their pricing errors are nearly collinear
  # whenever the model fits well. Solving the normal equations then fails, or
  # worse, succeeds inaccurately.
  whiten <- sigma_e_whitener(sigma_e)
  phi_tilde <- qr.solve(whiten(b_ols), whiten(-c_ols))

  # Pass 3: impose the constraint that the same Bbar multiplies X[t+1] and
  # -Phi~ X[t], by regressing on the single composite regressor.
  w <- x_led - x_lag %*% t(phi_tilde)
  z3 <- cbind(1, w)
  coef3 <- qr.solve(z3, rr)
  alpha <- coef3[1L, ]
  b_bar <- t(coef3[-1L, , drop = FALSE])                      # N x k

  e_gls <- rr - z3 %*% coef3
  sigma_e <- crossprod(e_gls) / n_obs

  # Pass 4: GLS for mu-tilde (Eq. S30), from alpha = -(Bbar mu~ + 1/2 B).
  quad <- rowSums((b_bar %*% sigma) * b_bar)
  whiten <- sigma_e_whitener(sigma_e)
  mu_tilde <- drop(qr.solve(whiten(b_bar), whiten(-(alpha + 0.5 * quad))))

  # --- short rate ----------------------------------------------------------
  zr <- cbind(1, x)
  dcoef <- qr.solve(zr, r)
  delta0 <- dcoef[1L]
  delta1 <- dcoef[-1L]

  # --- inflation loadings (Eq. S31) ----------------------------------------
  pi_fit <- acmy_fit_inflation(x, rx_real, inflation, real_maturities, n_max,
                               mu_tilde, phi_tilde, sigma, delta0, delta1,
                               fix_pi0)

  list(
    mu = mu, phi = phi, sigma = sigma,
    mu_tilde = mu_tilde, phi_tilde = phi_tilde,
    lambda0 = mu - mu_tilde,
    lambda1 = phi - phi_tilde,
    delta0 = delta0, delta1 = delta1,
    pi0 = pi_fit$pi0, pi1 = pi_fit$pi1,
    alpha = alpha, b_bar = b_bar,
    sigma_e = sigma_e,
    n_nominal_returns = n_nom,
    n_real_returns = n_real,
    var_residuals = v,
    return_residuals = e_gls,
    inflation_fit = pi_fit
  )
}

#' Pin down the inflation loadings from real bond returns
#'
#' Eq. S31. Everything else is held at its GLS value and `(pi0, pi1)` is
#' chosen to minimise the sum of squared real return pricing errors. Only
#' `1 + k` parameters, and the objective is smooth, so a plain BFGS from a
#' sensible start is enough.
#'
#' The starting value matters more than the optimiser. `pi1` is initialised by
#' regressing realised inflation on the factors, which is the projection the
#' model would use if inflation were priced by ordinary least squares, and
#' `pi0` at the sample mean residual. Started from zero the search can settle
#' in a flat region where real bonds carry no inflation exposure at all.
#'
#' @inheritParams acmy_closed_form
#' @param mu_tilde,phi_tilde,sigma,delta0,delta1 Held fixed.
#' @return A list with `pi0`, `pi1` and convergence diagnostics.
#' @keywords internal
#' @noRd
acmy_fit_inflation <- function(x, rx_real, inflation, real_maturities, n_max,
                               mu_tilde, phi_tilde, sigma, delta0, delta1,
                               fix_pi0 = NULL) {
  tt <- nrow(x)
  k <- ncol(x)
  x_lag <- x[-tt, , drop = FALSE]
  x_led <- x[-1L, , drop = FALSE]

  zi <- cbind(1, x)
  icoef <- qr.solve(zi, inflation)

  # `fix_pi0` holds the inflation intercept at a supplied value instead of
  # estimating it, which is what the UK specification of the Supplementary
  # Appendix does ("Average RPI inflation during this sample period is 2.48%
  # ... we fix pi0 in the estimation"). With a shorter real curve and fewer
  # return maturities, the level and the slope of inflation are weakly
  # separated, and pinning the level at its sample mean is the sane response.
  # `TRUE` means the sample mean of realised inflation.
  pi0_fixed <- NULL
  if (!is.null(fix_pi0) && !identical(fix_pi0, FALSE)) {
    pi0_fixed <- if (isTRUE(fix_pi0)) mean(inflation) else as.numeric(fix_pi0)
    if (length(pi0_fixed) != 1L || !is.finite(pi0_fixed)) {
      stop("`fix_pi0` must be TRUE, FALSE, or a single finite number in ",
           "monthly rate units.", call. = FALSE)
    }
  }

  start <- if (is.null(pi0_fixed)) c(icoef[1L], icoef[-1L]) else icoef[-1L]

  unpack <- function(theta) {
    if (is.null(pi0_fixed)) list(pi0 = theta[1L], pi1 = theta[-1L])
    else list(pi0 = pi0_fixed, pi1 = theta)
  }

  objective <- function(theta) {
    p <- unpack(theta)
    coefs <- acmy_coefficients(n_max, mu_tilde, phi_tilde, sigma,
                               delta0, delta1, p$pi0, p$pi1)
    if (!all(is.finite(coefs$b_real))) return(.Machine$double.xmax)

    model <- acmy_real_returns(coefs, real_maturities, x_lag, x_led,
                               mu_tilde, phi_tilde, sigma, p$pi0, p$pi1)
    if (!all(is.finite(model))) return(.Machine$double.xmax)
    sum((rx_real - model)^2)
  }

  fit <- stats::optim(start, objective, method = "BFGS",
                      control = list(maxit = 2000L, reltol = 1e-12))
  out <- unpack(fit$par)

  list(
    pi0 = out$pi0,
    pi1 = out$pi1,
    pi0_fixed = !is.null(pi0_fixed),
    pi0_ols = icoef[1L],
    pi1_ols = icoef[-1L],
    converged = fit$convergence == 0L,
    convergence = fit$convergence,
    iterations = unname(fit$counts[["gradient"]]),
    objective = fit$value
  )
}

#' Refuse to proceed when a factor does not move bond returns
#'
#' The prices of risk are identified from the exposure of excess returns to
#' the factors. A factor with no such exposure leaves `B` column-rank
#' deficient and every GLS step of the estimator insoluble. Reporting which
#' factor is responsible is the whole point: the error LAPACK raises names
#' none of them.
#'
#' @param b `N x k` matrix of return loadings.
#' @param k Number of factors.
#' @param nms Factor names, for the message.
#' @keywords internal
#' @noRd
check_return_loadings <- function(b, k, nms = NULL) {
  nms <- nms %||% paste0("factor ", seq_len(k))
  rank <- qr(b)$rank
  if (rank >= k) return(invisible(TRUE))

  # Rank deficiency is the condition that matters, but on its own it does not
  # say which factor to remove. The weakest columns, measured against the
  # strongest, almost always are the culprits, so they are named.
  norms <- sqrt(colSums(b^2))
  rel <- norms / max(norms)
  suspects <- order(rel)[seq_len(min(k - rank, k))]

  stop(
    "Excess returns do not identify all ", k, " factors: their ", nrow(b),
    " x ", k, " loading matrix has rank ", rank, ", so the prices of risk ",
    "cannot be separated.\n",
    "The weakest exposure is to ",
    paste(sprintf("'%s' (%.1e of the largest loading)",
                  nms[suspects], rel[suspects]), collapse = ", "), ".\n",
    if ("liquidity" %in% nms[suspects]) {
      paste0("A liquidity factor the bonds do not price is usually a ",
             "mis-aligned or irrelevant series: check that it lines up with ",
             "the panel dates and actually moves with TIPS valuations. ",
             "Refitting without `liquidity` will proceed regardless.")
    } else {
      paste0("Lower `n_factors_nominal` or `n_factors_real`, or widen the ",
             "return maturities so the cross-section can tell the factors ",
             "apart.")
    },
    call. = FALSE
  )
}

#' A whitening transform for the return pricing-error covariance
#'
#' Returns a function that left-multiplies by `L^-1`, where `Sigma_e = L L'`.
#' Applying it to both sides of a GLS problem turns it into ordinary least
#' squares, which can then be solved from the design matrix by QR without ever
#' forming a Gram matrix.
#'
#' `Sigma_e` is the matrix in this estimator most likely to be numerically
#' awkward. It is `N x N` for `N` return series -- twenty in the paper's US
#' specification -- all driven by the same handful of factors, so the better
#' the model fits, the more nearly collinear the pricing errors are and the
#' closer it comes to singular. A simulated panel generated exactly from the
#' model is the extreme case: the errors are pure floating-point noise and
#' `Sigma_e` has no usable inverse at all.
#'
#' A Cholesky factorisation is therefore attempted, and a failure is reported
#' as what it is rather than allowed to surface as a linear algebra error from
#' somewhere further in.
#'
#' @param sigma_e An `N x N` covariance matrix.
#' @return A function of a matrix or vector.
#' @keywords internal
#' @noRd
sigma_e_whitener <- function(sigma_e) {
  ch <- tryCatch(chol(sigma_e), error = function(e) NULL)

  if (is.null(ch)) {
    stop(
      "The covariance of return pricing errors is not positive definite, so ",
      "the GLS steps of the estimator have no solution. With ", nrow(sigma_e),
      " return series driven by a few factors this happens when the errors ",
      "are collinear -- most often because the panel was generated exactly ",
      "from an affine model and carries no pricing error at all, or because ",
      "`return_maturities` asks for more series than the cross-section can ",
      "support. Use fewer return maturities.",
      call. = FALSE
    )
  }

  lower <- t(ch)
  function(m) {
    if (is.matrix(m)) forwardsolve(lower, m) else drop(forwardsolve(lower, m))
  }
}
