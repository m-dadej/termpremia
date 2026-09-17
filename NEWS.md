# termpremia 0.0.0.9000

First public release, for feedback. **The API is not yet stable** — it may
change before the first CRAN submission. Pin a commit if you build on it.

## Estimation

* `atsm()` fits the Adrian, Crump & Moench (2013) three-step regression
  estimator to a zero-coupon yield panel, returning a decomposition into a
  risk-neutral (expected average short rate) component and a term premium.
  Validated against the New York Fed's published series: correlation 0.9997 in
  levels and 0.9935 in monthly changes over 763 month-ends, fitting the curve
  to 1.48bp RMSE. See `vignette("validation")`.
* `term_premium()`, `risk_neutral()`, `expected_short_rate()` extract
  components; `fitted()`, `residuals()`, `coef()`, `predict()`, `print()`,
  `summary()` and `plot()` methods are provided.
* `predict()` applies fitted parameters to new factor observations, which is
  how a model estimated monthly is evaluated at a higher frequency.
* `short_rate` argument accepts an external one-period rate. This matters more
  than it sounds: substituting an actual bill rate for the GSW fitted one-month
  yield moves the estimated 10-year US term premium by about 27bp.
* `p_dynamics = "brw"` applies the Bauer, Rudebusch & Wu (2012) small-sample
  bias correction to the factor VAR, configured with `brw_control()`. Least
  squares understates the persistence of near-unit-root factors, which flattens
  the expectations component and pushes the variation it should have carried
  into the term premium instead. The correction holds the risk-adjusted
  dynamics fixed, so fitted yields are unchanged to machine precision and only
  the decomposition moves: on US data it takes the factor half-life from 6.5
  years to effectively permanent, and cuts the standard deviation of the
  10-year term premium from 139bp to 120bp.

## Data handling

* `yield_panel()` accepts long data frames or wide matrices and normalises to
  decimal yields and maturities in months. Validation is deliberately loud,
  because supplying percent for decimal, or years for months, does not error
  downstream — it produces plausible but wrong term premia.
* Curves carry `instrument` and `issuer` metadata, since benchmark choice moves
  estimates materially.
* `svensson_curve()` and `svensson_yield()` evaluate Nelson-Siegel and Svensson
  parameters onto any maturity grid, verified against the Fed's own published
  tenors to 5e-05.

## Diagnostics

* `atsm()` warns when the risk-adjusted dynamics are explosive. This is not
  hypothetical: a sparse cross-section was found producing fitted yields of
  order 1e267 with no error raised.
* `atsm()` reports when model-implied expected short rates breach the lower
  bound. Fitted on US data spanning the ZLB, the model reaches −3.15%.
* `term_premium_survey()` gives a model-free benchmark from survey
  expectations. It returns one row per survey date and never silently
  interpolates.

## Data

* `gsw_monthly` — US Treasury zero-coupon yields, month-end, 1961–2024,
  maturities 1–120 months, with extrapolated observations flagged.
* `acm_published` — the New York Fed's published ACM decomposition, for
  validation. Governed by the New York Fed's terms, not this package's MIT
  licence; see `LICENSE.note`.
* `tbill_3m` — 3-month Treasury bill rate, as an alternative short rate.
