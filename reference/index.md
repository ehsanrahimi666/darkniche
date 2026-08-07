# Package index

## Environmental space

One shared lattice, one density per object.

- [`dn_space()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_space.md)
  : Define a shared environmental-space lattice
- [`dn_centres()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_centres.md)
  : Lattice cell centres
- [`dn_density()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_density.md)
  : Kernel density of a point set on the lattice
- [`dn_moments()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_moments.md)
  : First and second moments of a lattice density

## Coverage geometry

Envelopes, cumulative density rank, and the coverage weight vectors.

- [`dn_envelope()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_envelope.md)
  : q-quantile niche envelope
- [`dn_entryq()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_entryq.md)
  : Entry quantile of every lattice cell
- [`dn_weights()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_weights.md)
  : Coverage weight vector for one species

## Metrics

The directional family, with the symmetric statistics computed alongside
for comparison.

- [`dn_metrics()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_metrics.md)
  : All pairwise niche metrics for a block of species and regions

## Null models

Three ways of asking whether a coverage value is more than chance.

- [`dn_null_breadth()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_null_breadth.md)
  : Breadth-matched null (N1)
- [`dn_null_shift()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_null_shift.md)
  : Environmental-shift null (N2)
- [`dn_null_relocate()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_null_relocate.md)
  : Geographic relocation null (N3)
- [`dn_null_summary()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_null_summary.md)
  : Summarise a null matrix
- [`dn_ses()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_ses.md)
  : Standardised effect size with a floored null SD

## Calibration

Fixing what “admissible” means, against documented introductions.

- [`dn_match()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_match.md)
  : Realm- and area-matched negatives for introduced-range validation
- [`dn_auc()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_auc.md)
  : AUC with DeLong variance (fast algorithm of Sun & Xu 2014)
- [`dn_delong()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_delong.md)
  : Paired DeLong test between two correlated AUCs
- [`dn_calibrate()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_calibrate.md)
  : Calibrate the admissibility threshold (Youden)

## Dark-diversity accounting

The conjunction rule and the area-adjusted ratio.

- [`dn_dark()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_dark.md)
  : Dark-diversity classification of species-region pairs
- [`dn_cddr()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_cddr.md)
  : Climatic dark-diversity ratio and its area adjustment

## Region-level inference and its diagnostic

The Poisson-binomial count test with reliability flags, and the
threshold sweep that shows when its direction is an artefact.

- [`dn_aggregate()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_aggregate.md)
  : Region-level aggregate admissibility test (single threshold)
- [`dn_sweep()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_sweep.md)
  : Threshold-sweep diagnostic for the aggregate null test
- [`dn_direction_census()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_direction_census.md)
  : Direction census of an aggregate null test across thresholds
- [`dn_report()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_report.md)
  : One-call region report

## Data

- [`dn_data()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_data.md)
  : Load bundled ecoregion example data

## Figures and tables

- [`dn_plot_cddr_area()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_plot_cddr_area.md)
  : CDDR against area with the confound displayed (papers' Figure 7
  style)
- [`dn_plot_containment()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_plot_containment.md)
  : Conceptual containment figure (papers' Figure 1)
- [`dn_plot_density()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_plot_density.md)
  : Plot a density in environmental space
- [`dn_plot_matrix()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_plot_matrix.md)
  : Region-by-taxon matrix heatmap (papers' Figure 9 style)
- [`dn_plot_roc()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_plot_roc.md)
  : ROC curves for metric comparison (papers' Figure 3 style)
- [`dn_plot_sources()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_plot_sources.md)
  : Stacked source-decomposition bars (papers' Figure 11 style)
- [`dn_plot_sweeps()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_plot_sweeps.md)
  : Multi-region threshold-sweep figure (papers' Figure 10 style)
- [`dn_table_extremes()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_table_extremes.md)
  : Extremes table (papers' Table 6 style)
- [`dn_table_export()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_table_export.md)
  : Export a table as CSV, Markdown or LaTeX

## Package

- [`darkniche`](https://ehsanrahimi666.github.io/darkniche/reference/darkniche-package.md)
  [`darkniche-package`](https://ehsanrahimi666.github.io/darkniche/reference/darkniche-package.md)
  : darkniche: Climatic Dark Diversity and Directional Niche Coverage
