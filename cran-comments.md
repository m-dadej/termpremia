## Submission

This is a new submission.

## Test environments

* Local: Windows 11, R 4.5.2 and R 4.3.2
* GitHub Actions: macOS (R 4.6.2), Windows (R 4.6.1), Ubuntu (R 4.6.2,
  R-devel 4.7.0, R 4.5.3 and R 4.3.3, the declared minimum)

<!-- Before submitting: run devtools::check_win_devel(), add its result
     above ("win-builder: R-devel"), then delete this comment. -->

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

Any words flagged as possibly misspelled in DESCRIPTION are author names
(Adrian, Crump, Moench, Rudebusch, Orphanides, Abrahams, Yu, Siegel,
Svensson).

## Notes for reviewers

* Internet resources. Several functions download public data on the user's
  behalf: `acm_published()`, `gsw_tips()`, `fred_series()`,
  `boe_yield_curve()`, `ons_rpi()`, `spf_tbill()` and `spf_bill10()`. They
  fail with an informative error when the resource is unavailable. Their
  examples, and those of `atsm_real()` and its methods, which need the
  downloaded inflation-indexed curve, are wrapped in
  `@examplesIf interactive()`. No test run on CRAN uses the network.

* Data licensing. The bundled datasets derive from Federal Reserve Board
  releases that assert no copyright; `LICENSE.note` records their sources.
  Data whose terms do not fit the package's MIT licence (the New York Fed's
  published term premia, the Philadelphia Fed's Survey of Professional
  Forecasters) is downloaded on request and never redistributed.

* The vignette is precomputed, because it compares the package's estimates
  with the New York Fed's published series, which is downloaded rather than
  bundled. Building it needs neither the data nor a network connection.
