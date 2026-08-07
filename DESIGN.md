# darkniche — package design

Directional niche coverage, climatic dark diversity, and threshold-sweep
null-model diagnostics. R port of the `ddniche` Python pipeline (v2.1)
underlying Rahimi & Jung (Iran companion study, in review) and the
global ecoregion study.

## Design principles

1.  **Data-agnostic core.** The package ships *no* range polygons, no
    ecoregion shapefile, no rasters. Users bring point sets or cell
    scores in a shared environmental space; the package never downloads
    anything. IUCN polygons and the WWF shapefile are inputs to the
    *user’s* workflow, not package data (both are too large and
    licence-encumbered to bundle). The vignette demonstrates the full
    workflow, and the published *results* for 619 ecoregions are bundled
    via
    [`dn_data()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_data.md)
    so every example runs on real output.
2.  **Numerical fidelity.** Every core function reproduces the Python
    reference to float32 precision, enforced by golden tests generated
    from the actual `ddniche` code (`tests/testthat/test-golden.R`,
    `tools/make_golden.py`). Passing golden tests are the package’s
    correctness claim; any future refactor must keep them green.
3.  **Guardrails are API, not prose.** The two field-tested failure
    modes are encoded in the objects themselves:
    [`dn_aggregate()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_aggregate.md)
    returns `reliable = FALSE` (and suppresses SES/enrichment) when the
    expected count is \< 5, and
    [`dn_sweep()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_sweep.md)
    prints a direction-reversal warning when the verdict flips across
    thresholds.
    [`dn_cddr()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_cddr.md)
    computes the area-adjusted residual alongside the raw ratio and
    documents that raw CDDR must never be used to rank regions.
4.  **Base R first.** v0.1 has zero hard dependencies beyond base +
    stats. `terra` is Suggests-only (for users projecting rasters). An
    `Rcpp` backend for
    [`dn_density()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_density.md)/[`dn_metrics()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_metrics.md)
    is a v0.2 optimisation once the API is frozen — profiling shows the
    pure-R KDE handles the package’s intended per-study scale
    (10^(2–10)4 objects) comfortably; the 21.7M-pair global run remains
    the Python pipeline’s job.

## Module map

| File | Exports | Source |
|----|----|----|
| `R/espace.R` | `dn_space`, `dn_centres`, `dn_density`, `dn_moments` | `niche.build_pcspace`, `niche.kde_grid`, `niche.density_moments` |
| `R/metrics.R` | `dn_envelope`, `dn_entryq`, `dn_metrics`, `dn_weights` | `metrics.envelopes_flat`, `metrics.entry_quantiles`, `metrics.pair_metrics`, `nulls._coverage_weights` |
| `R/nulls.R` | `dn_ses`, `dn_null_breadth`, `dn_null_shift`, `dn_null_relocate`, `dn_null_summary` | `utils.ses`, `nulls.breadth_matched_null`, `nulls.shift_null_dnc`, `nulls.relocation_null_dnc`, `nulls.summarise_null` |
| `R/aggregate.R` | `dn_aggregate`, `dn_sweep` (+print/plot), `dn_dark`, `dn_cddr`, `dn_auc`, `dn_delong`, `dn_match`, `dn_calibrate` | `nulls.aggregate_admissibility_test`, stage-07 sweep protocol, `darkniche`, `validate` |
| `R/data.R` | `dn_data`, `dn_direction_census`, `dn_report` (+print) | new in the R package; no Python counterpart |
| `R/plots.R` | `dn_plot_*` (7), `dn_table_extremes`, `dn_table_export` | figure designs of the source papers, reimplemented in base graphics |
| `R/darkniche-package.R` | package documentation | — |

## Key API decisions

- **`dn_null_shift(targets=)`** accepts explicit target cells so a draw
  can be reproduced exactly (used by the golden tests; also lets users
  share seeds across languages). Default behaviour draws
  availability-weighted targets.
- **`dn_dark(p_min = 0.177)`**: the default is the externally calibrated
  Youden threshold of the source studies’ species-blocked, climate-only
  establishment model. It is documented as *their* calibration; users
  applying different metrics/models must recalibrate via
  [`dn_calibrate()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_calibrate.md)
  on their own validation positives.
- **[`dn_sweep()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_sweep.md)
  returns a classed object** so that printing it in a console or an Rmd
  cannot silently drop the reversal warning.
- **Array layout**: densities are `R x R` matrices indexed `[x, y]`
  (Python’s first axis = x). Flat vectors inside the package are R
  column-major; the golden tests handle the C-order conversion
  explicitly and document it.

## Deliberate exclusions (v0.1)

- The GBM establishment model (`sklearn` in the pipeline): out of scope;
  the calibrated threshold is provided as a documented constant, and
  [`dn_calibrate()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_calibrate.md)
  covers recalibration. A `glm`-based equivalent is a candidate for
  v0.2.
- Tail-refined permutation p-values (`utils.permutation_p` refinement):
  [`dn_null_breadth()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_null_breadth.md)/[`dn_null_summary()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_null_summary.md)
  return the floored empirical p only. SES, null mean and SD — the
  quantities the dark-diversity rule uses — are golden-tested; the
  refined p is not yet ported. **TODO(v0.2).**
- The 3-D lattice (`niche.PCSpaceND`): sensitivity-analysis machinery,
  not core. TODO(v0.2) behind the same API (`dn_space(ndim = 3)`).
- Figures: the package plots only the sweep; publication figures stay in
  user scripts.

## Test strategy

1.  **Golden tests** (56 assertions): lattice construction, KDE vs scipy
    `gaussian_filter` (\< 2e-7), envelopes and entry quantiles (exact),
    all twelve pair metrics (\< 5e-6; float32 reference), shift null
    with injected targets (exact for `q`; \< 2e-4 for `auc`, limited by
    the reference’s own float32 weight storage), SES flooring (exact),
    aggregate test incl. reliability flags (exact), breadth-matched null
    (\< 1e-6), and the containment phenomenon of the papers’ Figure 1
    (DNC_q ≈ 1, DNC_auc ≈ 0.9, D ≈ 0.08 on the engulfing-niche pair).
2.  **Real-data reproduction** (pre-release gate): run
    `tools/export_golden_real.py` against the actual pipeline store to
    dump ~20 species x 10 ecoregions of real densities and metrics, and
    confirm the R port reproduces published table values. This is the
    step that needs files from the authors’ machine (see below).

## Repository layout

    darkniche/
      R/            package source (7 files, 33 exported functions)
      man/          generated Rd (roxygen2; do not edit by hand)
      inst/extdata/
        golden/     fixtures from the Python reference implementation
        ecoregions/ published results for 619 WWF terrestrial ecoregions
      tests/        162 assertions (test-golden.R numerical, test-workflow.R behavioural)
      vignettes/    quickstart tutorial
      analysis/     scripts regenerating every figure and table of the paper
        figures/    generated PNGs
        output/     generated CSVs
      manual/       darkniche-manual.pdf (R CMD Rd2pdf)
      manuscript/   the software paper (.docx and .pdf)
      tools/        fixture generators (Python)
      .github/      R-CMD-check workflow

## Development workflow

``` bash
Rscript -e 'roxygen2::roxygenise()'          # regenerate man/ and NAMESPACE
R CMD INSTALL --no-docs .
Rscript -e 'testthat::test_dir("tests/testthat", package = "darkniche")'
R CMD build . && R CMD check darkniche_0.1.0.tar.gz
R CMD Rd2pdf --no-preview --force \
  --title="darkniche: Climatic Dark Diversity and Directional Niche Coverage in R -- Reference Manual" \
  -o manual/darkniche-manual.pdf .
Rscript analysis/01_simulation.R && Rscript analysis/02_ecoregions.R \
  && Rscript analysis/03_schematics.R
```

## Defect log

Kept deliberately, because each was found by a test rather than by
reading:

1.  **[`dn_sweep()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_sweep.md)
    used unreliable thresholds in its verdict** (fixed in 0.1.0 before
    release). The reversal detection, the printed ratio range and the
    plotted curve were computed over all finite ratios, including rows
    where
    [`dn_aggregate()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_aggregate.md)
    had suppressed the effect size for an expected count below five — so
    the diagnostic could report exactly the numbers the test refuses to
    report. All three now use reliable rows only, and
    `attr(x, "n_reliable")` exposes the count. Regression test: *“the
    sweep verdict never uses cells the test calls unreliable”*.
2.  **[`dn_auc()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_auc.md)
    returned `NaN` for single-class labels** (fixed in 0.1.0). A
    validation set with no negatives produced `NaN` silently, which in a
    scripted workflow propagates into a results table without any signal
    that something went wrong.
    [`dn_auc()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_auc.md),
    [`dn_delong()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_delong.md)
    and
    [`dn_calibrate()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_calibrate.md)
    now stop with a message naming the problem, and
    [`dn_auc()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_auc.md)
    also checks that `score` and `label` have equal length. Regression
    test: *“degenerate inputs fail loudly rather than returning
    nonsense”*.
3.  **[`dn_cddr()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_cddr.md)
    warned on a constant ratio.**
    [`cor()`](https://rdrr.io/r/stats/cor.html) emits a
    zero-standard-deviation warning when every region has the same
    ratio; `confound_r` is now `NA` in that case and the report and plot
    say “not estimable”.
4.  **README verification snippet did not run as written** (fixed
    post-release, in documentation only). Two problems, both found by a
    user pasting it into a fresh session rather than by any test: the
    metric example printed values from an unseeded draw, so the printed
    output never matched the comment; and it told the reader to run
    `testthat::test_local(system.file(package = "darkniche"))`, which
    fails because R does not install a package’s tests unless
    `INSTALL_opts = "--install-tests"` is passed. Documentation examples
    that quote expected output must either be seeded or must not quote
    output.
5.  **`CFI_w` cell-ordering.** The Mahalanobis weight grid must be
    enumerated in the same order as the flattened density. R densities
    flatten column-major (x fastest, matching `expand.grid`); the Python
    fixtures are C-order (x slowest). The golden test builds its own
    grid for the fixture ordering and documents why.

## Roadmap

- **v0.1** (this release): core verified against the reference
  implementation, bundled ecoregion data, two case studies, PDF manual,
  `R CMD check` clean.
- **Distribution**: GitHub only, installed with
  `remotes::install_github()`. A pkgdown site is built and deployed to
  GitHub Pages by CI on every push to `main`. CRAN submission is
  deliberately not a goal: an existing CRAN package, `DarkDiv` (Carmona
  & Partel), estimates dark diversity from co-occurrence, and a
  near-name collision would confuse users of both.
- **v0.2**: compiled backends for
  [`dn_density()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_density.md)
  and
  [`dn_metrics()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_metrics.md);
  permutation-p refinement; three-dimensional space; GLM establishment
  helper; pkgdown content expansion (further articles).
- **Software paper**: Methods in Ecology and Evolution (Applications).
  The R-versus-Python agreement doubles as an
  independent-reimplementation robustness statement citable in both
  empirical papers’ revisions.
