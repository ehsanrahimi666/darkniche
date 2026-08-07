# Conceptual containment figure (papers' Figure 1)

Generates the three-panel demonstration that symmetric overlap cannot
answer a containment question: a region nested in a niche core, nested
at the margin, and outside the niche, with D, DNC_q and DNC_auc
annotated.

## Usage

``` r
dn_plot_containment(R = 100, broad = 14, narrow = 4, offsets = c(0, 18, 34))
```

## Arguments

- R:

  lattice resolution (default 100).

- broad, narrow:

  Gaussian widths (lattice cells) of niche and region.

- offsets:

  horizontal offsets of the region for the three panels.

## Value

invisibly, a data.frame of the annotated metric values.

## Examples

``` r
dn_plot_containment()
```
