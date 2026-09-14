# termpremia

Decompose zero-coupon government bond yields into expected average future short rates
and a term premium:

```
yield(n) = expected average future short rate(n) + term premium(n)
```

**Status: early development. Not yet usable.**

## Why this package

R can fit Nelson–Siegel and Svensson curves (`YieldCurve`, `yieldcurves`), but curve
fitting is interpolation, not no-arbitrage pricing — it never produces a term premium.
`MultiATSM` estimates affine models by maximum likelihood and does report term premia,
risk-neutral yields and expected short rates, but it targets multi-country
unspanned-macro-risk research and its interface is built around matrices and Excel input.

`termpremia` fills a narrower gap:

- **ACM three-step regressions** (Adrian, Crump & Moench 2013) — closed-form OLS, so
  estimation is fast, deterministic, and cannot fail to converge. Not currently
  available in R.
- **A desk-usable API** — hand it a curve, get a term premium.
- **Validated against the New York Fed's published ACM series.**

## Scope

One question: given a zero-coupon curve, what is the no-arbitrage decomposition of yields?

**In scope.** ACM three-step regressions; term premia by tenor, risk-neutral yields,
expected short-rate paths, expected excess returns; a model-free survey benchmark;
pluggable real-world dynamics (OLS, survey-augmented); joint real–nominal decomposition
yielding inflation risk premia and model-implied expected inflation.

**Out of scope.** Fitting curves from bond prices, shadow-rate models, multi-country
GVAR estimation, credit, derivatives.

The package is curve-agnostic: it accepts any zero-coupon curve, and never requires a
particular data vendor.

## A caution on levels

Term premium models disagree, and the disagreement is not small. Cohen, Hördahl & Xia
(*BIS Quarterly Review*, September 2018) find that ACM and Hördahl–Tristani estimates
can differ by as much as 200bp in level while agreeing closely on direction — monthly
change correlations of 0.77–0.92 for the US.

Treat any single model's *level* with scepticism. This package reports level agreement
and change agreement separately for that reason, and provides a model-free survey
benchmark as an independent anchor.

## References

- Adrian, T., R. Crump and E. Moench (2013). "Pricing the term structure with linear
  regressions." *Journal of Financial Economics* 110(1), 110–138.
- Abrahams, M., T. Adrian, R. Crump, E. Moench and R. Yu (2016). "Decomposing real and
  nominal yield curves." *Journal of Monetary Economics*.
- Cohen, B., P. Hördahl and D. Xia (2018). "Term premia: models and some stylised facts."
  *BIS Quarterly Review*, September.

## Licence

MIT. Any third-party data distributed with the package is governed by its own source's
terms rather than by this licence; those terms are recorded alongside the data.
