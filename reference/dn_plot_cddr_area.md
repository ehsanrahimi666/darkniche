# CDDR against area with the confound displayed (papers' Figure 7 style)

Raw ratio against log area with the fitted confound line, the top
positive and negative residuals labelled – the display that shows the
confound rather than hiding it.

## Usage

``` r
dn_plot_cddr_area(cddr_tab, area, labels = NULL, n_label = 5)
```

## Arguments

- cddr_tab:

  output of
  [`dn_cddr()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_cddr.md)
  (with `area` supplied).

- area:

  the area vector used.

- labels:

  optional region names for labelling extremes.

- n_label:

  how many extremes to label on each side (default 5).

## Value

invisibly, the labelled subset.

## Examples

``` r
set.seed(1)
area <- 10^runif(80, 4, 6.7)
nd <- rpois(80, 30 + 200 / log10(area)); np <- rpois(80, 60)
x <- dn_cddr(nd, np, area)
dn_plot_cddr_area(x, area, labels = paste0("eco", 1:80))
```
