# Multi-region threshold-sweep figure (papers' Figure 10 style)

Draws every region's observed/expected curve on a log axis with the
cross-region median, colouring curves by their direction at the lowest
threshold – the display in which the bimodal-vs-unimodal crossing is
visible at a glance.

## Usage

``` r
dn_plot_sweeps(sweeps, main = "aggregate null test across thresholds")
```

## Arguments

- sweeps:

  a list of
  [`dn_sweep()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_sweep.md)
  results (one per region).

- main:

  title.

## Value

invisibly, the share of regions whose direction reverses.

## Examples

``` r
set.seed(1)
sweeps <- lapply(1:15, function(i) {
  obs <- c(rbeta(250, 0.3, 3), rbeta(80, 4, 1.5))
  dn_sweep(obs, rbeta(330, 1, 2.2), 0.08 + 0.1 * runif(330))
})
dn_plot_sweeps(sweeps)
```
