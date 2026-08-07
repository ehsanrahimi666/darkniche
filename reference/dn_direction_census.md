# Direction census of an aggregate null test across thresholds

Summarises, over many regions, how many are called enriched and how many
depleted at each threshold, and what share of regions change direction
somewhere in the sweep. This is the quantity a study must report before
any single-threshold enrichment claim can be trusted.

## Usage

``` r
dn_direction_census(sweep, reliable_only = TRUE)
```

## Arguments

- sweep:

  a long data.frame with columns `region`, `quantile` and `enrichment`
  (observed/expected), and optionally `reliable`.

- reliable_only:

  if `TRUE` (default) only rows flagged reliable are used; unreliable
  cells (expected count \< 5) are not interpretable.

## Value

list with `census` (per threshold: n enriched, n depleted), `n_regions`,
`n_changed`, `frac_changed`, `n_monotone`.

## Examples

``` r
cen <- dn_direction_census(dn_data("sweep"))
cen$census
#>   quantile n_enriched n_depleted
#> 1     0.80         71        548
#> 2     0.90        344        275
#> 3     0.95        520         81
#> 4     0.99        449          3
round(cen$frac_changed, 3)     # share of ecoregions that flip direction
#> [1] 0.971
```
