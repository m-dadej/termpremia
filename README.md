# termpremia

Decompose zero-coupon government bond yields into expected average future short
rates and a term premium, using the Adrian, Crump & Moench (2013) three-step
regression estimator.

```
yield(n) = expected average future short rate(n) + term premium(n)
```

The package is working and validated but the API is not yet stable.

## Installation

```r
# install.packages("remotes")
remotes::install_github("m-dadej/termpremia", build_vignettes = TRUE)
```

No compiled code and no hard dependencies beyond base R, so it installs
anywhere without a toolchain.

## Quick start

```r
library(termpremia)

panel <- yield_panel(gsw_monthly, units = "percent", maturity_unit = "months",
                     instrument = "government", issuer = "US")

fit <- atsm(panel, n_factors = 5)

term_premium(fit, maturity = 120)   # 10-year term premium
plot(fit)                           # yield, expected short rate, premium
```

Bring your own curve — the package is curve-agnostic and takes any zero-coupon
panel, in long or wide form:

```r
my_panel <- yield_panel(my_data, date = "date", maturity = "tenor",
                        yield = "zero_rate", units = "percent",
                        maturity_unit = "years", issuer = "DE")
```

If you hold Nelson-Siegel or Svensson parameters rather than yields (as
published by the Fed, Bundesbank, and others), `svensson_curve()` evaluates
them onto whatever maturity grid you need.

## Reproducibility

Against the New York Fed's published ACM series, 763 month-ends, 1961–2024:

| | |
|---|---|
| Correlation, levels | **0.9997** |
| Correlation, monthly changes | **0.9935** |
| Fit to the observed curve | **1.48 bp** RMSE |
| `fitted = risk_neutral + term_premium` | exact to 7e-18 |

For scale, the BIS report change correlations of 0.77–0.92 between published
implementations of *different* term structure models.

`vignette("validation")` reproduces this and then investigates the part that
does **not** match: a residual ~15bp level difference, which the vignette
localises entirely to the risk-neutral component and traces to the choice of
short rate. A validation that reports only the agreement is marketing.

## Stability and scope

**The API may change before the first CRAN release.** The point of publishing
here is to get feedback while changes are still cheap. If you build something
on it, pin a commit. Breaking changes will be noted in `NEWS.md`.

In scope: ACM estimation; term premia by tenor, risk-neutral yields, expected
short-rate paths, expected excess returns; multiple curves; a model-free survey
benchmark.

Not in scope: fitting curves from bond prices, shadow-rate models, joint
multi-country (GVAR) estimation, credit, derivatives.

On the roadmap: Bauer-Rudebusch-Wu bias correction, survey-augmented dynamics,
and the joint real-nominal decomposition.

## Licence and bundled data

The package is MIT licensed. **The bundled datasets are not.** Each is governed
by its own source's terms, recorded in
[`LICENSE.note`](LICENSE.note) — read it before redistributing.

In particular, `acm_published` is redistributed under the New York Fed's Terms
of Use, which require attribution, clear labelling of modifications, and that
**the same permissions pass to anyone you pass the data to**:

> © Federal Reserve Bank of New York. Content from the New York Fed subject to
> the Terms of Use at newyorkfed.org.

`gsw_monthly` and `tbill_3m` derive from Federal Reserve Board releases, which
assert no copyright.

## References

- Adrian, T., R. K. Crump & E. Moench (2013). "Pricing the term structure with
  linear regressions." *Journal of Financial Economics* 110(1), 110–138.
- Cohen, B., P. Hördahl & D. Xia (2018). "Term premia: models and some stylised
  facts." *BIS Quarterly Review*, September.
- Gürkaynak, R. S., B. Sack & J. H. Wright (2007). "The U.S. Treasury yield
  curve: 1961 to the present." *Journal of Monetary Economics* 54(8), 2291–2304.
