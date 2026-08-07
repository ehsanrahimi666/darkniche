# Entry quantile of every lattice cell

For each cell, the smallest envelope level q at which it enters the
q-envelope – its cumulative density rank within the niche. Cells at the
mode enter near 0; peripheral cells enter near 1; unoccupied cells are
assigned 1. `1 - dn_entryq(Z)` is the weight vector whose inner product
with a region density gives `DNC_auc`.

## Usage

``` r
dn_entryq(Z)
```

## Arguments

- Z:

  density vector (or matrix, flattened internally).

## Value

numeric vector of entry quantiles in \[0, 1\].

## Examples

``` r
sp <- dn_space(cbind(rnorm(3000), rnorm(3000)), R = 50)
Z <- dn_density(cbind(rnorm(200), rnorm(200)), sp)
q <- dn_entryq(Z)
range(q)                                 # mode enters near 0, edges near 1
#> [1] 0.004590167 1.000000000
```
