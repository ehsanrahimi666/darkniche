# Climatic dark-diversity ratio and its area adjustment

`CDDR = n_dark / (n_dark + n_present)`. The raw ratio is confounded with
region area (a small region occupies little environmental space, so more
niches enclose it, while its smaller area supports fewer natives; in the
global analysis r = -0.44), so cross-region comparison must use the
residual from the CDDR ~ log(area) regression. **Never rank regions by
raw CDDR.**

## Usage

``` r
dn_cddr(n_dark, n_present, area = NULL)
```

## Arguments

- n_dark, n_present:

  integer vectors per region.

- area:

  region areas (same order); required for `adjusted`.

## Value

data.frame with `dn_cddr`, and when `area` is given `cddr_residual` plus
the fitted confound as attributes `fit`.

## Examples

``` r
set.seed(1)
area <- 10^runif(50, 4, 6.7)
nd <- rpois(50, 40); np <- rpois(50, 60)
x <- dn_cddr(nd, np, area)
attr(x, "confound_r")     # never rank regions by raw cddr
#> [1] -0.02210871
head(x)
#>        cddr cddr_residual
#> 1 0.3861386   -0.01699436
#> 2 0.3644860   -0.03825963
#> 3 0.3448276   -0.05718871
#> 4 0.3809524   -0.01984547
#> 5 0.3518519   -0.05151302
#> 6 0.4044944    0.00366086
```
