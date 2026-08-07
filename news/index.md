# Changelog

## darkniche 0.1.0

First release.

### Environmental space and metrics

- [`dn_space()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_space.md),
  [`dn_density()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_density.md),
  [`dn_centres()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_centres.md),
  [`dn_moments()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_moments.md)
  — one shared lattice for every object; fixed-bandwidth kernel
  densities reproducing `scipy.ndimage.gaussian_filter` to within 2e-7.
- [`dn_envelope()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_envelope.md),
  [`dn_entryq()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_entryq.md),
  [`dn_weights()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_weights.md)
  — q-envelopes, cumulative density rank, and the weight vectors behind
  each coverage statistic.
- [`dn_metrics()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_metrics.md)
  — the directional family (DNC_q, DNC_auc, NSI_q, NSI_min, DNC_min,
  CFI, CFI_w, NOV) with Schoener’s D, Warren’s I and the Bhattacharyya
  coefficient computed alongside for comparison.

### Null models

- [`dn_null_breadth()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_null_breadth.md)
  (N1),
  [`dn_null_shift()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_null_shift.md)
  (N2, reproducible via injectable shift targets),
  [`dn_null_relocate()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_null_relocate.md)
  (N3).
- [`dn_ses()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_ses.md),
  [`dn_null_summary()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_null_summary.md)
  — standardised effect sizes with a floored null standard deviation.

### Calibration, accounting and inference

- [`dn_match()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_match.md),
  [`dn_auc()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_auc.md),
  [`dn_delong()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_delong.md),
  [`dn_calibrate()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_calibrate.md)
  — realm- and area-matched negatives, fast DeLong AUC and paired tests,
  Youden calibration.
- [`dn_dark()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_dark.md),
  [`dn_cddr()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_cddr.md),
  [`dn_report()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_report.md)
  — conjunction rule, area-adjusted dark-diversity ratio, one-call
  region report.
- [`dn_aggregate()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_aggregate.md)
  — Poisson-binomial region-level test with `expected >= 5` reliability
  suppression.
- [`dn_sweep()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_sweep.md),
  [`dn_direction_census()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_direction_census.md)
  — the threshold-sweep diagnostic with reversal detection in `print`
  and `plot`.

### Reporting and data

- Seven plotting functions and two table exporters (CSV, Markdown,
  LaTeX).
- [`dn_data()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_data.md)
  — published results for all 619 WWF terrestrial ecoregions.

### Quality

- 162 unit tests: 56 golden assertions against an independent Python
  reference implementation (`tools/make_golden.py` regenerates the
  fixtures), 106 behavioural.
- `R CMD check` clean; CI on Linux, Windows and macOS.
- `analysis/` regenerates every figure and table of the software paper.

### Installation notes

R does not install a package’s tests by default. To run the suite from
an installed copy:

``` r

remotes::install_github("ehsanrahimi666/darkniche",
                        INSTALL_opts = "--install-tests")
testthat::test_package("darkniche")
```

### Distribution

Installed from GitHub with
`remotes::install_github("ehsanrahimi666/darkniche")`. Not submitted to
CRAN: the existing CRAN package `DarkDiv` estimates dark diversity from
co-occurrence data, and a near-name collision would confuse users of
both. Documentation is published to GitHub Pages by CI.

### Roadmap (0.2.0)

- Compiled backends for
  [`dn_density()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_density.md)
  and
  [`dn_metrics()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_metrics.md),
  validated by the existing golden tests.
- Three-dimensional environmental space behind the same interface.
- Tail-refined permutation p-values.
- A GLM establishment-model helper closing the calibration loop in R.
- Further pkgdown articles (a full worked reanalysis).
