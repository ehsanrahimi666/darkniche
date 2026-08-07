# AUC with DeLong variance (fast algorithm of Sun & Xu 2014)

AUC with DeLong variance (fast algorithm of Sun & Xu 2014)

## Usage

``` r
dn_auc(score, label)
```

## Arguments

- score:

  numeric vector; higher = more positive.

- label:

  logical or 0/1 vector.

## Value

list with `auc`, `var`, `ci` (95%).

## Examples

``` r
set.seed(1)
lab <- rep(c(TRUE, FALSE), each = 200)
score <- ifelse(lab, rnorm(200, 1), rnorm(200, 0))
dn_auc(score, lab)
#> $auc
#> [1] 0.760675
#> 
#> $var
#> [1] 0.0005545369
#> 
#> $ci
#> [1] 0.7145206 0.8068294
#> 
```
