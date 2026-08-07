# darkniche 0.1.0

First release.

## Environmental space and metrics

* `dn_space()`, `dn_density()`, `dn_centres()`, `dn_moments()` — one shared
  lattice for every object; fixed-bandwidth kernel densities reproducing
  `scipy.ndimage.gaussian_filter` to within 2e-7.
* `dn_envelope()`, `dn_entryq()`, `dn_weights()` — q-envelopes, cumulative
  density rank, and the weight vectors behind each coverage statistic.
* `dn_metrics()` — the directional family (DNC_q, DNC_auc, NSI_q, NSI_min,
  DNC_min, CFI, CFI_w, NOV) with Schoener's D, Warren's I and the
  Bhattacharyya coefficient computed alongside for comparison.

## Null models

* `dn_null_breadth()` (N1), `dn_null_shift()` (N2, reproducible via
  injectable shift targets), `dn_null_relocate()` (N3).
* `dn_ses()`, `dn_null_summary()` — standardised effect sizes with a floored
  null standard deviation.

## Calibration, accounting and inference

* `dn_match()`, `dn_auc()`, `dn_delong()`, `dn_calibrate()` — realm- and
  area-matched negatives, fast DeLong AUC and paired tests, Youden
  calibration.
* `dn_dark()`, `dn_cddr()`, `dn_report()` — conjunction rule, area-adjusted
  dark-diversity ratio, one-call region report.
* `dn_aggregate()` — Poisson-binomial region-level test with `expected >= 5`
  reliability suppression.
* `dn_sweep()`, `dn_direction_census()` — the threshold-sweep diagnostic with
  reversal detection in `print` and `plot`.

## Reporting and data

* Seven plotting functions and two table exporters (CSV, Markdown, LaTeX).
* `dn_data()` — published results for all 619 WWF terrestrial ecoregions.

## Quality

* 169 unit tests: 56 golden assertions against an independent Python
  reference implementation (`tools/make_golden.py` regenerates the fixtures),
  113 behavioural.
* `R CMD check` clean; CI on Linux, Windows and macOS.
* `analysis/` regenerates every figure and table of the software paper.

## Installation notes

R does not install a package's tests by default. To run the suite from an
installed copy:

```r
remotes::install_github("ehsanrahimi666/darkniche",
                        INSTALL_opts = "--install-tests")
testthat::test_package("darkniche")
```

## Distribution

Installed from GitHub with `remotes::install_github("ehsanrahimi666/darkniche")`.
Not submitted to CRAN: the existing CRAN package `DarkDiv` estimates dark
diversity from co-occurrence data, and a near-name collision would confuse
users of both. Documentation is published to GitHub Pages by CI.

## Roadmap (0.2.0)

* Compiled backends for `dn_density()` and `dn_metrics()`, validated by the
  existing golden tests.
* Three-dimensional environmental space behind the same interface.
* Tail-refined permutation p-values.
* A GLM establishment-model helper closing the calibration loop in R.
* Further pkgdown articles (a full worked reanalysis).
