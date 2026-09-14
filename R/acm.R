# The Adrian, Crump and Moench (2013) three-step regression estimator.
#
# Equation numbers refer to Federal Reserve Bank of New York Staff Report 340.
#
# Units. Everything in this file works in MONTHLY rate units: an annualised
# decimal yield of 0.048 enters as 0.004. Log prices are then p(n) = -n * y(n)
# with n counted in months, and the one-period short rate is the one-month
# yield. Results are annualised again on the way out.

#' Extract principal component factors from a yield panel
#'
#' Principal components of yields are the observable pricing factors in ACM's
#' baseline specification.
#'
#' The sign and ordering of principal components are not identified by the
#' decomposition itself: flipping an eigenvector gives an equally valid answer
#' with the opposite-signed factor. Left alone this makes results irreproducible
#' across platforms and even across LAPACK versions. Each component is therefore
#' signed so that its largest-magnitude loading is positive, which is
#' deterministic and independent of the numerical path taken.
#'
#' @param y A `T x N` numeric matrix of yields, dates in rows.
#' @param k Number of components to retain.
#'
#' @return A list with `scores` (`T x k`), `loadings` (`N x k`), `center`
#'   (length `N`) and `sdev` (length `k`).
#'
#' @keywords internal
#' @noRd
acm_factors <- function(y, k) {
  if (anyNA(y)) {
    stop("Yield matrix contains missing values; ACM needs a complete panel.",
         call. = FALSE)
  }
  if (k < 1 || k > min(dim(y))) {
    stop("`k` must be between 1 and ", min(dim(y)), ".", call. = FALSE)
  }

  center <- colMeans(y)
  yc <- sweep(y, 2L, center, "-")

  sv <- svd(yc, nu = 0, nv = k)
  loadings <- sv$v
  sdev <- sv$d[seq_len(k)] / sqrt(max(nrow(y) - 1L, 1L))

  # Deterministic sign convention: largest-magnitude loading is positive.
  for (j in seq_len(k)) {
    lead <- which.max(abs(loadings[, j]))
    if (loadings[lead, j] < 0) loadings[, j] <- -loadings[, j]
  }

  list(
    scores = yc %*% loadings,
    loadings = loadings,
    center = center,
    sdev = sdev
  )
}

#' One-period excess log holding returns
#'
#' \eqn{rx^{(n-1)}_{t+1} = p^{(n-1)}_{t+1} - p^{(n)}_t - r_t} (Eq. 6): buy an
#' `n`-month bond, hold it one month, sell it as an `(n-1)`-month bond, and
#' subtract the riskless one-month return.
#'
#' Forming this requires both `n` and `n - 1` to be present in the maturity
#' grid, which is why ACM is run on a dense monthly grid rather than the sparse
#' set of tenors that curves are usually published at.
#'
#' @param p A `T x N` matrix of log prices in monthly units, dates in rows,
#'   with `colnames` giving maturities in months.
#' @param maturities Integer maturities, in months, matching `colnames(p)`.
#' @param return_maturities Maturities `n` for which to form returns.
#'
#' @return A `(T-1) x length(return_maturities)` matrix of excess returns.
#'
#' @keywords internal
#' @noRd
acm_excess_returns <- function(p, maturities, return_maturities) {
  idx_n <- match(return_maturities, maturities)
  idx_nm1 <- match(return_maturities - 1L, maturities)

  if (anyNA(idx_n) || anyNA(idx_nm1)) {
    bad <- return_maturities[is.na(idx_n) | is.na(idx_nm1)]
    stop(
      "Excess returns need both n and n-1 on the maturity grid. Missing for ",
      "maturity(ies): ", paste(bad, collapse = ", "), ".",
      call. = FALSE
    )
  }

  tt <- nrow(p)
  short_idx <- match(1L, maturities)
  if (is.na(short_idx)) {
    stop("The one-month maturity is required as the short rate but is absent ",
         "from the grid.", call. = FALSE)
  }
  r <- -p[, short_idx]  # one-month yield in monthly units

  rx <- p[-1L, idx_nm1, drop = FALSE] -
    p[-tt, idx_n, drop = FALSE] -
    r[-tt]

  colnames(rx) <- as.character(return_maturities)
  rx
}

#' ACM three-step estimator
#'
#' @param x A `T x k` matrix of pricing factors, dates in rows.
#' @param rx A `(T-1) x N` matrix of excess returns.
#' @param r Length-`T` vector of one-period short rates, monthly units.
#'
#' @return A list of estimated parameters.
#' @keywords internal
#' @noRd
acm_three_step <- function(x, rx, r) {
  tt <- nrow(x)
  k <- ncol(x)

  # --- Step 1: VAR(1) on the factors (Eq. 1) --------------------------------
  x_lag <- x[-tt, , drop = FALSE]
  x_led <- x[-1L, , drop = FALSE]

  z1 <- cbind(1, x_lag)
  var_coef <- qr.solve(crossprod(z1), crossprod(z1, x_led))

  mu <- var_coef[1L, ]
  phi <- t(var_coef[-1L, , drop = FALSE])   # X[t+1] = mu + phi %*% X[t]
  v <- x_led - z1 %*% var_coef              # (T-1) x k innovations

  # ACM normalise by the number of VAR observations, not T - 1 degrees of
  # freedom: Sigma-hat = V-hat V-hat' / T.
  n_var <- nrow(v)
  sigma <- crossprod(v) / n_var

  # --- Step 2: excess returns on constant, innovations, lagged factors ------
  # (Eq. 14-15)
  z <- cbind(1, v, x_lag)
  coef <- qr.solve(crossprod(z), crossprod(z, rx))   # (1 + 2k) x N

  a <- coef[1L, ]                                    # N
  beta <- coef[2L:(k + 1L), , drop = FALSE]          # k x N
  cc <- t(coef[(k + 2L):(2L * k + 1L), , drop = FALSE])  # N x k

  e <- rx - z %*% coef
  n_mat <- ncol(rx)
  sigma2 <- sum(e^2) / (n_mat * n_var)

  # --- Step 3: prices of risk (Eq. 16-17) ----------------------------------
  # B* has rows vec(beta_n beta_n'), so B* vec(Sigma) is the convexity term
  # bond by bond.
  b_star <- t(apply(beta, 2L, function(b) as.vector(tcrossprod(b))))
  if (k == 1L) b_star <- matrix(b_star, ncol = 1L)

  bb <- tcrossprod(beta)                             # k x k
  adj <- a + 0.5 * (b_star %*% as.vector(sigma) + sigma2)

  lambda0 <- qr.solve(bb, beta %*% adj)              # k x 1
  lambda1 <- qr.solve(bb, beta %*% cc)               # k x k

  # --- Short rate loadings: r_t = delta0 + delta1' X_t ----------------------
  zr <- cbind(1, x)
  dcoef <- qr.solve(crossprod(zr), crossprod(zr, r))
  delta0 <- dcoef[1L]
  delta1 <- dcoef[-1L]

  list(
    mu = mu, phi = phi, sigma = sigma, sigma2 = sigma2,
    a = a, beta = beta, cc = cc,
    lambda0 = drop(lambda0), lambda1 = lambda1,
    delta0 = delta0, delta1 = delta1,
    var_residuals = v, return_residuals = e
  )
}

#' Largest eigenvalue modulus of a matrix
#'
#' The pricing recursion iterates `B[n]' = B[n-1]' (Phi - lambda1) - delta1'`
#' up to 120 times. If `Phi - lambda1` has any eigenvalue outside the unit
#' circle, `B` grows geometrically and yields diverge -- silently, to values of
#' order 1e267, because nothing in the arithmetic errors. A sparse or
#' ill-conditioned cross-section of excess returns is enough to cause it, so
#' this is checked rather than assumed.
#'
#' @param m A square matrix.
#' @return The largest absolute eigenvalue.
#' @keywords internal
#' @noRd
spectral_radius <- function(m) {
  max(Mod(eigen(m, only.values = TRUE)$values))
}

#' Bond pricing recursion
#'
#' Equations 25-27. Yields are affine in the factors,
#' \eqn{y^{(n)}_t = -(A_n + B_n' X_t)/n}.
#'
#' Two details that are easy to get wrong:
#'
#' The `+ sigma2 / 2` term does **not** appear in the textbook affine
#' recursion. It is specific to ACM, arising because they allow a
#' maturity-specific return pricing error that is conditionally orthogonal to
#' the state innovations, and so fold pricing errors into the no-arbitrage
#' recursion. Substituting a standard recursion gives plausible but wrong
#' term premia.
#'
#' The recursion is seeded at `n = 1` with \eqn{A_1 = -\delta_0},
#' \eqn{B_1 = -\delta_1} rather than iterated from \eqn{A_0 = B_0 = 0}. A
#' one-period bond held for one period is riskless and carries no return
#' pricing error, so the `sigma2` term must not enter at the first step.
#'
#' @param n_max Longest maturity, in months.
#' @param pars Output of `acm_three_step()`.
#' @param risk_neutral If `TRUE`, set both prices of risk to zero, giving the
#'   expected-average-short-rate component.
#'
#' @return A list with `a` (length `n_max`) and `b` (`n_max x k`).
#' @keywords internal
#' @noRd
acm_recursion <- function(n_max, pars, risk_neutral = FALSE) {
  k <- length(pars$mu)

  lambda0 <- if (risk_neutral) rep(0, k) else pars$lambda0
  lambda1 <- if (risk_neutral) matrix(0, k, k) else pars$lambda1

  a <- numeric(n_max)
  b <- matrix(0, nrow = n_max, ncol = k)

  a[1L] <- -pars$delta0
  b[1L, ] <- -pars$delta1

  mu_adj <- pars$mu - lambda0
  phi_adj <- pars$phi - lambda1

  if (n_max >= 2L) {
    for (n in 2L:n_max) {
      bp <- b[n - 1L, ]
      a[n] <- a[n - 1L] +
        sum(bp * mu_adj) +
        0.5 * (drop(crossprod(bp, pars$sigma %*% bp)) + pars$sigma2) -
        pars$delta0
      b[n, ] <- drop(bp %*% phi_adj) - pars$delta1
    }
  }

  list(a = a, b = b)
}

#' Model-implied yields from recursion coefficients
#'
#' @param rec Output of `acm_recursion()`.
#' @param x A `T x k` factor matrix.
#' @param maturities Maturities in months to return.
#' @return A `T x length(maturities)` matrix of yields in monthly units.
#' @keywords internal
#' @noRd
acm_yields <- function(rec, x, maturities) {
  out <- matrix(NA_real_, nrow = nrow(x), ncol = length(maturities),
                dimnames = list(rownames(x), as.character(maturities)))

  for (j in seq_along(maturities)) {
    n <- maturities[j]
    out[, j] <- -(rec$a[n] + drop(x %*% rec$b[n, ])) / n
  }
  out
}
