# Extremes table (papers' Table 6 style)

The `n` highest and lowest regions by the area-adjusted residual, with
the raw ratio shown alongside – never the ranking key.

## Usage

``` r
dn_table_extremes(region, cddr_tab, n_obs = NULL, n_dark = NULL, n = 8)
```

## Arguments

- region:

  region names.

- cddr_tab:

  output of
  [`dn_cddr()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_cddr.md)
  (with `area`).

- n_obs, n_dark:

  optional richness columns to carry.

- n:

  rows per side (default 8).

## Value

data.frame, top block then bottom block.

## Examples

``` r
set.seed(1)
area <- 10^runif(40, 4, 6.5)
x <- dn_cddr(rpois(40, 40), rpois(40, 60), area)
dn_table_extremes(paste0("eco", 1:40), x)
#>    region      cddr    residual
#> 31  eco31 0.5048544  0.11052285
#> 30  eco30 0.4851485  0.08939036
#> 40  eco40 0.4651163  0.07007204
#> 20  eco20 0.4536082  0.06224981
#> 39  eco39 0.4495413  0.05764197
#> 14  eco14 0.4509804  0.05566266
#> 17  eco17 0.4356436  0.04368292
#> 25  eco25 0.4382022  0.04170800
#> 28  eco28 0.3529412 -0.04239383
#> 19  eco19 0.3518519 -0.04350683
#> 13  eco13 0.3483146 -0.04395401
#> 22  eco22 0.3500000 -0.04704866
#> 32  eco32 0.3363636 -0.05678530
#> 38  eco38 0.3333333 -0.06476417
#> 4    eco4 0.3139535 -0.07608872
#> 23  eco23 0.2941176 -0.09850678
```
