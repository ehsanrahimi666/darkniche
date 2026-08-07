# q-quantile niche envelope

The smallest set of lattice cells containing a fraction `q` of the
density mass. This replaces the naive `Z > 0` support, which after
kernel smoothing covers most of the lattice and makes coverage metrics
bandwidth-driven rather than ecology-driven.

## Usage

``` r
dn_envelope(Z, q = 0.95, flat = FALSE)
```

## Arguments

- Z:

  density as a matrix or vector (one object), or a matrix with one *row
  per object* if `flat = TRUE`.

- q:

  envelope level in (0, 1\]; `q >= 1` returns the positive support.

- flat:

  if `TRUE`, treat `Z` as (m x P) flattened densities and return an (m
  x P) logical matrix.

## Value

logical vector/matrix marking envelope membership.

## Examples

``` r
sp <- dn_space(cbind(rnorm(3000), rnorm(3000)), R = 50)
Z <- dn_density(cbind(rnorm(200), rnorm(200)), sp)
E <- dn_envelope(Z, q = 0.95)
mean(E)                                  # share of lattice in the envelope
#> [1] 0.302
```
