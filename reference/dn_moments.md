# First and second moments of a lattice density

Mean vector and covariance matrix of a density in environmental space,
used by the core-fit index (CFI).

## Usage

``` r
dn_moments(Z, sp)
```

## Arguments

- Z:

  `R x R` density matrix.

- sp:

  a `dn_space`.

## Value

list with `mu` (length 2) and `cov` (2 x 2).

## Examples

``` r
sp <- dn_space(cbind(rnorm(3000), rnorm(3000)), R = 50)
Z <- dn_density(cbind(rnorm(200, 1, 0.3), rnorm(200, 0, 0.4)), sp)
dn_moments(Z, sp)
#> $mu
#> [1] 0.9787128 0.0456484
#> 
#> $cov
#>              [,1]         [,2]
#> [1,]  0.121159283 -0.004367732
#> [2,] -0.004367732  0.156977085
#> 
```
