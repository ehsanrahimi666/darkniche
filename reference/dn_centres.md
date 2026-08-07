# Lattice cell centres

Lattice cell centres

## Usage

``` r
dn_centres(sp)
```

## Arguments

- sp:

  a `dn_space`

## Value

list with `x`, `y` centre vectors

## Examples

``` r
sp <- dn_space(cbind(rnorm(1000), rnorm(1000)), R = 40)
ctr <- dn_centres(sp)
str(ctr)
#> List of 2
#>  $ x: num [1:40] -3.3 -3.14 -2.97 -2.8 -2.64 ...
#>  $ y: num [1:40] -2.98 -2.82 -2.66 -2.5 -2.34 ...
```
