# termpremia

Decompose zero-coupon government bond yields into expected average future short
rates and a term premium, using the Adrian, Crump & Moench (2013) three-step
regression estimator. With a second, inflation-indexed curve it goes further
and splits breakeven inflation into expected inflation and an inflation risk
premium, following Abrahams, Adrian, Crump, Moench & Yu (2016).

```
yield(n) = expected average future short rate(n) + term premium(n)

breakeven(n) = expected inflation(n) + inflation risk premium(n)
```

The package is working and validated, but the API is not yet stable.

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

Least squares understates how persistent yield factors are, which flattens the
expectations component and inflates the term premium — badly so on short
samples. `p_dynamics = "brw"` applies the Bauer, Rudebusch & Wu (2012)
bootstrap correction, which leaves the fit to the yield curve untouched and
changes only the split:

```r
atsm(panel, n_factors = 5, p_dynamics = "brw")
```

A model that sees only yields also has no way to know where rates settle, so it
reads every persistent move as a shift in the steady state. `p_dynamics =
"survey"` anchors the dynamics to published forecasts instead —
Kim–Wright style, not a replication of Kim–Wright:

```r
surveys <- rbind(spf_tbill(), spf_bill10())   # downloaded, not bundled
atsm(panel, n_factors = 5, p_dynamics = "survey", survey = surveys)
```

Use both horizons. Short-horizon forecasts alone make the far end *worse*; see
[`analysis/validate-spf.R`](analysis/validate-spf.R).

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

## Forwards, and why you probably want them

A ten-year spot rate mixes the next five years with the five after that. The
five-to-ten year forward removes the near end, where policy is known and
expectations are sharp, and isolates the part of the curve where risk premia
dominate. The source papers state their conclusions in forward terms, so
comparisons against them should be made there too.

```r
forward_rate(fit, start = 60, end = 120)  # fitted, risk-neutral, term premium
forward_curve(fit)                        # the whole one-month forward curve
```

Every component of the spot decomposition is affine in the factors, so its
forward counterpart is that same difference of coefficients — the identities
that hold for spot rates hold exactly here as well.

Expected one-month excess returns, the quantity the estimator is actually
fitted to, come out of the same object:

```r
expected_excess_return(fit, maturity = 120, convexity = TRUE)
```

## Real and nominal together

Running `atsm()` separately on a nominal and a real curve and differencing the
results does **not** decompose breakeven inflation: the two fits would have
unrelated factor spaces and unrelated prices of risk. `atsm_real()` prices both
curves off one state vector and one stochastic discount factor, which is what
makes the inflation risk premium a well-defined object.

```r
nom  <- gsw_monthly; nom$curve  <- "nominal"
tips <- gsw_tips();  tips$curve <- "tips"        # downloaded, not bundled

panel <- yield_panel(
  rbind(nom[, c("date", "maturity", "yield", "curve")],
        tips[, c("date", "maturity", "yield", "curve")]),
  curve = "curve", units = "percent", maturity_unit = "months",
  instrument = c(nominal = "government", tips = "tips"), issuer = "US"
)

fit <- atsm_real(panel,
                 inflation  = fred_series("CPIAUCNS"),
                 short_rate = fred_series("DFF", frequency = "monthly"),
                 method     = "ml")

breakeven(fit, maturity = 120)
expected_inflation(fit, maturity = 120)
inflation_risk_premium(fit, maturity = 120)
real_term_premium(fit, maturity = 120)
```

Both identities — nominal yield into risk-neutral plus term premium, breakeven
into expected inflation plus inflation risk premium — hold by construction
rather than approximately, and the tests assert it.

`method = "closed_form"` is the regression estimator the paper uses as its
*starting value*; it fits excess returns and leaves yield levels unconstrained.
`method = "ml"` goes on to the paper's actual constrained maximum likelihood,
which is what pins the levels down. It costs a few minutes and it is worth
them: on US data it takes the nominal RMSE from 11.05bp to 4.10bp, the TIPS
RMSE from 12.06bp to 4.16bp, and the risk-adjusted spectral radius from
1.0044 — explosive — to 0.9915.

TIPS trade less liquidly than nominal Treasuries and that illiquidity is
priced. Supply a liquidity index through `liquidity`, optionally built with
`tips_liquidity_factor()`. Without one the model still identifies expected
inflation and the inflation risk premium, but the latter absorbs the liquidity
premium — a real limitation, and one `print()` reports rather than leaves
implicit.

None of this is US-only. `boe_yield_curve()` and `ons_rpi()` fetch the UK
inputs, and the paper's own UK specification runs without a liquidity factor;
see [`analysis/validate-uk.R`](analysis/validate-uk.R).

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
short rate.

The joint real-nominal model has no published benchmark series to correlate
against, so it is checked against every number its source paper reports. On the
paper's own US sample — 1999:01–2014:11, T = 191, matched exactly —
[`analysis/replicate-acmy.R`](analysis/replicate-acmy.R) keeps nominal and TIPS
pricing errors inside the stated tolerances and reproduces both headline
findings: expected inflation stable at 2.06% on the five-to-ten year forward
with a 17bp standard deviation, against their "2.1 to 2.5%" and "quite stable",
and a variance decomposition led by the real term premium with the inflation
risk premium a minor contributor.
[`analysis/validate-uk.R`](analysis/validate-uk.R) does the same against the
Supplementary Appendix's UK specification, landing on its T = 336 sample
without tuning anything — the start date turns out to be a property of the Bank
of England's real curve — and fitting the real curve to 1.17bp RMSE.

## Stability and scope

**The API may change before the first CRAN release.** The point of publishing
here is to get feedback while changes are still cheap. If you build something
on it, pin a commit. Breaking changes will be noted in `NEWS.md`.

In scope: ACM estimation; Bauer-Rudebusch-Wu bias correction; survey-augmented
dynamics; the joint real-nominal decomposition; term premia by tenor,
risk-neutral yields, expected short-rate paths, expected excess returns,
forward rates and their decomposition; multiple curves; a model-free survey
benchmark.

Not in scope: fitting curves from bond prices, shadow-rate models, joint
multi-country (GVAR) estimation, credit, derivatives.

On the roadmap: model comparison diagnostics and expanding-window
re-estimation.

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

Everything else is fetched rather than bundled, either because its licence does
not permit redistribution or because the file is too large to carry:
`gsw_tips()`, `fred_series()`, `spf_tbill()` and `spf_bill10()` for the US,
`boe_yield_curve()` and `ons_rpi()` for the UK. They all take a `cache` path,
so a script need only download once.

## References

- Abrahams, M., T. Adrian, R. K. Crump, E. Moench & R. Yu (2016). "Decomposing
  real and nominal yield curves." *Journal of Monetary Economics* 84, 182–200.
- Adrian, T., R. K. Crump & E. Moench (2013). "Pricing the term structure with
  linear regressions." *Journal of Financial Economics* 110(1), 110–138.
- Bauer, M. D., G. D. Rudebusch & J. C. Wu (2012). "Correcting estimation bias
  in dynamic term structure models." *Journal of Business & Economic
  Statistics* 30(3), 454–467.
- Kim, D. H. & A. Orphanides (2012). "Term structure estimation with survey
  data on interest rate forecasts." *Journal of Financial and Quantitative
  Analysis* 47(1), 241–272.
- Cohen, B., P. Hördahl & D. Xia (2018). "Term premia: models and some stylised
  facts." *BIS Quarterly Review*, September.
- Gürkaynak, R. S., B. Sack & J. H. Wright (2007). "The U.S. Treasury yield
  curve: 1961 to the present." *Journal of Monetary Economics* 54(8), 2291–2304.
- Gürkaynak, R. S., B. Sack & J. H. Wright (2010). "The TIPS yield curve and
  inflation compensation." *American Economic Journal: Macroeconomics* 2(1),
  70–92.
