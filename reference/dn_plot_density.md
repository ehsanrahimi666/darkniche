# Plot a density in environmental space

Filled image of a lattice density with optional envelope contour and a
second density overlaid as contour lines – the building block of
species-vs-region panels.

## Usage

``` r
dn_plot_density(Z, sp, q = 0.95, overlay = NULL, main = "", ...)
```

## Arguments

- Z:

  `R x R` density (species or region).

- sp:

  a `dn_space`.

- q:

  optional envelope level to outline (default 0.95; `NA` for none).

- overlay:

  optional second density drawn as contours (e.g. the region over a
  species niche).

- main:

  title.

- ...:

  passed to
  [`graphics::image()`](https://rdrr.io/r/graphics/image.html).

## Value

invisibly, `Z`.

## Examples

``` r
set.seed(1)
sp <- dn_space(cbind(rnorm(4000), rnorm(4000)), R = 60)
Zs <- dn_density(cbind(rnorm(800, 0, 1), rnorm(800, 0, 1)), sp)
Ze <- dn_density(cbind(rnorm(300, 0.4, 0.25), rnorm(300, 0, 0.25)), sp)
dn_plot_density(Zs, sp, overlay = Ze, main = "region inside niche")
```
