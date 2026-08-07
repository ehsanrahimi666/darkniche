# Kernel density of a point set on the lattice

Bins the points on the lattice (half-open bins, right edge of the last
bin closed, exactly as `numpy.histogram2d`), smooths with a Gaussian
kernel of `sp$bw_cells` cells, and normalises to unit mass.

## Usage

``` r
dn_density(points, sp, weights = NULL)
```

## Arguments

- points:

  numeric matrix (n x 2) of scores (e.g. a species' occupied cells
  projected into the space, or a region's cells).

- sp:

  a `dn_space`.

- weights:

  optional per-point weights (e.g. cell areas).

## Value

numeric `R x R` matrix summing to 1 (all-zero if no point falls inside
the lattice).

## Examples

``` r
set.seed(1)
sp <- dn_space(cbind(rnorm(5000), rnorm(5000)), R = 60)
pts <- cbind(rnorm(300, 0.5, 0.4), rnorm(300, -0.2, 0.5))
Z <- dn_density(pts, sp)
sum(Z)                                   # unit mass
#> [1] 1
image(Z, main = "species niche density")
```
