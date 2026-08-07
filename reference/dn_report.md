# One-call region report

Assembles the quantities a region-level dark-diversity study must report
together: raw ratio, area-adjusted residual, rank on each, and the
direction census of the region's own threshold sweep. Printing the
result surfaces the two failure modes – an area-driven ranking and a
threshold-driven direction – instead of leaving them to a footnote.

## Usage

``` r
dn_report(eco, sweep = NULL)
```

## Arguments

- eco:

  data.frame with `region`, `n_obs`, `n_dark`, `area_km2`.

- sweep:

  optional long sweep table (as
  [`dn_data()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_data.md)
  `"sweep"`).

## Value

object of class `dn_report`.

## Examples

``` r
rep <- dn_report(dn_data("ecoregion"), dn_data("sweep"))
rep
#> Climatic dark-diversity report
#>   regions: 619
#>   CDDR: 0.022 - 0.924 (median 0.579)
#>   area confound: r = -0.441 -- rank on the residual, never on raw CDDR
#>   286 regions move more than 10% of the ranking when area is removed
#>   threshold sweep: 439 of 452 regions (97%) change direction
#>   direction census by threshold:
#>  quantile n_enriched n_depleted
#>      0.80         71        548
#>      0.90        344        275
#>      0.95        520         81
#>      0.99        449          3
#>   ** a single-threshold enrichment claim is not supportable here **
head(rep$table)
#>                                         region n_obs n_dark area_km2      cddr
#> 1                         Sinaloan dry forests   864   1651    93200 0.6564612
#> 2                                   Everglades   267   1684    22000 0.8631471
#> 3                       Veracruz moist forests  1051   1350    92800 0.5622657
#> 4                    Hispaniolan moist forests   277   1267    49600 0.8205959
#> 5                      Hispaniolan dry forests   302   1382    23600 0.8206651
#> 6 Trans-Mexican Volcanic Belt pine-oak forests  1177   2020    87200 0.6318424
#>       residual rank_raw rank_adjusted rank_shift
#> 1  0.092831186      219           194        -25
#> 2  0.199141643       11            73         62
#> 3 -0.001663362      334           346         12
#> 4  0.213111647       34            61         27
#> 5  0.161540661       33           106         73
#> 6  0.063585793      248           231        -17
```
