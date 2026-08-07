# Define a shared environmental-space lattice

Builds the common `R x R` lattice on which all species niches and all
region climate availabilities are estimated, from the global
distribution of environmental scores (typically the first two principal
components of climate over all terrestrial cells).

## Usage

``` r
dn_space(scores, R = 100L, clip_q = c(0.001, 0.999), bw_frac = 0.02)
```

## Arguments

- scores:

  numeric matrix (n x 2) of environmental scores for the *global*
  background (all terrestrial grid cells), not for any one species or
  region. The lattice must be shared or no two densities are comparable.

- R:

  lattice resolution per axis (default 100, as calibrated in the source
  studies; sensitivity there: Spearman \>= 0.996 across 60-150).

- clip_q:

  lower/upper quantiles at which the axes are clipped before 5% padding
  (default `c(0.001, 0.999)`).

- bw_frac:

  fixed kernel bandwidth as a fraction of the axis range (default 0.02;
  one bandwidth for every object so that estimated niche breadth
  reflects ecology, not sample size).

## Value

an object of class `dn_space`: list with `x_edges`, `y_edges`, `R`,
`bw_cells`.

## Examples

``` r
set.seed(1)
bg <- cbind(rnorm(5000), rnorm(5000))   # global background scores
sp <- dn_space(bg, R = 60)
sp
#> <dn_space> 60 x 60 lattice, bandwidth 1.20 cells
#>   PC1 [-3.508, 3.225]  PC2 [-3.425, 3.302]
```
