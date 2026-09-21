# Shared simulation helpers for the joint real-nominal model.
#
# A panel generated FROM the model at known parameters, so the estimator can
# be asked to get them back. Used by test-real-nominal.R and by
# test-real-nominal-ml.R, which is why these live in a helper rather than in
# either file.


# Simulation helper -------------------------------------------------------

# A stationary five-factor truth with a stable pricing measure. Correlated
# innovations are drawn through a Cholesky factor rather than MASS::mvrnorm,
# so the tests need no package that is not already a dependency.
# `noise_bp` adds i.i.d. pricing error to the simulated curves. It defaults
# to zero, which makes recovery exact and is the sharpest test available --
# but only because the GLS steps whiten by a Cholesky factor of the return
# error covariance instead of inverting it. A panel generated exactly from
# the model has pricing errors that are pure floating-point noise, and
# solving the normal equations against that covariance fails outright.
#
# When it is non-zero, keep it small. The noise is i.i.d. ACROSS maturities
# and a log price multiplies it by n, so two basis points on a ten-year yield
# is 240 basis points of noise in its log price, swamping the one-month
# return the estimator is built on. A fitted curve's own error is smooth in
# maturity and does no such thing.
sim_acmy <- function(tt = 260L, seed = 7L, liquidity = FALSE,
                     noise_bp = 0) {
  set.seed(seed)

  phi <- diag(c(0.990, 0.960, 0.930, 0.970, 0.940))
  phi[1L, 2L] <- 0.010
  phi[2L, 3L] <- -0.008
  phi[4L, 5L] <- 0.006
  mu <- c(4e-5, -1e-5, 5e-6, 2e-6, -1e-6)

  sd_f <- c(3e-3, 2e-3, 1e-3, 8e-4, 5e-4)
  sigma <- diag(sd_f^2)
  sigma[1L, 2L] <- sigma[2L, 1L] <- 0.2 * sd_f[1L] * sd_f[2L]

  lambda1 <- diag(c(0.02, 0.015, 0.01, 0.012, 0.008))
  delta1 <- c(1, -0.6, 0.3, 0.05, -0.02) * 1e-2

  # Inflation loadings large enough that inflation actually varies: at a
  # tenth of this, simulated annual inflation has a standard deviation of
  # under two basis points, and pi1 is then so weakly identified that the
  # recovery test measures rounding rather than arithmetic.
  pi1 <- c(0.15, -0.1, 0.05, 0.4, -0.25) * 1e-1

  if (liquidity) {
    # A sixth factor standing in for an observable liquidity index. Three
    # properties are needed for it to be a fair test rather than a
    # decoration:
    #
    #  - It is O(1) and mostly positive, like a standardised index shifted to
    #    be non-negative, which is what `tips_liquidity_factor()` returns.
    #    This is also the scale mismatch against O(1e-5) yield components
    #    that broke the GLS steps before the factors were normalised.
    #  - `delta1` puts no weight on it, so it leaves nominal yields alone --
    #    the paper's own finding (Section 3.3).
    #  - `pi1` does put weight on it, which is the only channel through which
    #    an extra state variable can reach real yields in this model, so
    #    TIPS returns load on it and its price of risk is identified. A
    #    liquidity series the bonds do not price leaves the estimator with a
    #    rank-deficient B, which it now refuses by name.
    phi <- rbind(cbind(phi, 0), c(rep(0, 5L), 0.950))
    mu <- c(mu, 0.05)
    sigma <- rbind(cbind(sigma, 0), c(rep(0, 5L), 0.12^2))
    lambda1 <- rbind(cbind(lambda1, 0), c(rep(0, 5L), 0.010))
    delta1 <- c(delta1, 0)
    pi1 <- c(pi1, 3e-4)
  }

  k <- length(mu)
  lambda0 <- rep(0, k)
  delta0 <- 0.0035
  pi0 <- 0.0018

  mu_t <- mu - lambda0
  phi_t <- phi - lambda1

  chol_s <- chol(sigma)
  x <- matrix(0, tt, k)
  if (liquidity) x[1L, 6L] <- mu[6L] / (1 - phi[6L, 6L])
  for (t in 2L:tt) {
    x[t, ] <- mu + drop(phi %*% x[t - 1L, ]) +
      drop(crossprod(chol_s, stats::rnorm(k)))
  }

  mat_n <- 1:120
  mat_r <- 23:120
  co <- acmy_coefficients(120L, mu_t, phi_t, sigma, delta0, delta1, pi0, pi1)

  firsts <- seq(as.Date("1999-01-01"), by = "month", length.out = tt + 1L)
  dates <- firsts[-1L] - 1L

  noise <- function(nr, nc) {
    if (noise_bp <= 0) return(matrix(0, nr, nc))
    matrix(stats::rnorm(nr * nc, sd = noise_bp / 1e4), nr, nc)
  }

  list(
    x = x, dates = dates, mat_n = mat_n, mat_r = mat_r, co = co,
    noise_bp = noise_bp,
    y_nom = acmy_yields(co$a, co$b, x, mat_n) * 12 +
      noise(tt, length(mat_n)),
    y_real = acmy_yields(co$a_real, co$b_real, x, mat_r) * 12 +
      noise(tt, length(mat_r)),
    inflation = pi0 + drop(x %*% pi1),
    short_rate = delta0 + drop(x %*% delta1),
    liquidity = if (liquidity) x[, 6L] else NULL,
    truth = list(mu = mu, phi = phi, sigma = sigma, mu_tilde = mu_t,
                 phi_tilde = phi_t, lambda0 = lambda0, lambda1 = lambda1,
                 delta0 = delta0, delta1 = delta1, pi0 = pi0, pi1 = pi1)
  )
}

sim_panel <- function(s, end = NULL) {
  long <- rbind(
    data.frame(date = rep(s$dates, length(s$mat_n)),
               maturity = rep(s$mat_n, each = length(s$dates)),
               yield = as.vector(s$y_nom) * 100, curve = "nominal",
               stringsAsFactors = FALSE),
    data.frame(date = rep(s$dates, length(s$mat_r)),
               maturity = rep(s$mat_r, each = length(s$dates)),
               yield = as.vector(s$y_real) * 100, curve = "tips",
               stringsAsFactors = FALSE)
  )
  if (!is.null(end)) long <- long[long$date <= end, , drop = FALSE]

  yield_panel(long, curve = "curve", units = "percent",
              maturity_unit = "months",
              instrument = c(nominal = "government", tips = "tips"),
              issuer = "SIM")
}

sim_cpi <- function(s) {
  data.frame(date = s$dates, value = 100 * exp(cumsum(s$inflation)),
             stringsAsFactors = FALSE)
}

sim_sr <- function(s) {
  data.frame(date = s$dates, value = s$short_rate * 12 * 100,
             stringsAsFactors = FALSE)
}

# Fits are cached per noise level: the estimator is the slow part of the
# suite and several tests want the same one.
sim_fit <- local({
  cache <- list()
  function(noise_bp = 0, ...) {
    key <- paste(noise_bp, ...)
    if (is.null(cache[[key]])) {
      s <- sim_acmy(noise_bp = noise_bp, ...)
      cache[[key]] <<- list(
        s = s,
        fit = suppressWarnings(suppressMessages(
          atsm_real(sim_panel(s), inflation = sim_cpi(s),
                    short_rate = sim_sr(s), short_rate_units = "percent")))
      )
    }
    cache[[key]]
  }
})

# The pricing-measure recursion coefficients of a fitted model.
sim_coefs <- function(f) {
  p <- f$pars
  acmy_coefficients(max(f$maturities, f$real_maturities), p$mu_tilde,
                    p$phi_tilde, p$sigma, p$delta0, p$delta1, p$pi0, p$pi1)
}

# The largest factor-consistency violation of a fitted model.
sim_violation <- function(f) {
  kk <- f$n_factors_nominal + f$n_factors_real
  max(abs(acmy_constraints(sim_coefs(f), f$pca, f$maturities,
                           f$real_maturities, f$n_factors, kk)))
}



# A deliberately tiny optimiser budget for tests.
#
# The maximum likelihood tests check plumbing -- that the step runs, records
# what it did, leaves the physical dynamics alone, honours a fixed pi0 and
# preserves the decomposition identities. None of that needs convergence, and
# the default budget of 20 rounds x 400 inner iterations over 42 parameters
# costs seven minutes of the suite on its own. Convergence quality is
# demonstrated on real data in analysis/validate-uk.R, where it belongs.
sim_ml_control <- function(max_outer = 2L, max_inner = 25L, ...) {
  acmy_control(max_outer = max_outer, max_inner = max_inner, ...)
}
