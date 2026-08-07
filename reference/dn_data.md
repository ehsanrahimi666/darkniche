# Load bundled ecoregion example data

Results of the global analysis of Rahimi & Jung (in review) for all 619
WWF terrestrial ecoregions (35,087 species; 21,718,853 species-ecoregion
pairs), shipped so that every example, test and vignette in `darkniche`
runs on real output rather than on toy numbers.

## Usage

``` r
dn_data(
  which = c("ecoregion", "sweep", "sources", "taxon", "region_taxon", "auc",
    "containment")
)
```

## Arguments

- which:

  name of the dataset (see above).

## Value

a data.frame.

## Details

Available datasets:

- `ecoregion`:

  619 rows: `region`, `n_obs`, `n_dark`, `cddr`, `area_km2`.

- `sweep`:

  2,476 rows: `region`, `quantile`, `threshold`, `enrichment`,
  `reliable` – the aggregate null test at four thresholds per ecoregion.

- `sources`:

  619 rows: `near`, `regional`, `distant` shares of each ecoregion's
  dark diversity by native-range distance.

- `taxon`:

  dark-diversity rate per taxonomic group with 95% CI.

- `region_taxon`:

  3,095 rows: CDDR for each ecoregion x taxon.

- `auc`:

  AUC of each metric against documented introduced ranges.

- `containment`:

  AUC by niche-to-region size stratum.

## Examples

``` r
eco <- dn_data("ecoregion")
str(eco)
#> 'data.frame':    619 obs. of  5 variables:
#>  $ region  : chr  "Sinaloan dry forests" "Everglades" "Veracruz moist forests" "Hispaniolan moist forests" ...
#>  $ n_obs   : int  864 267 1051 277 302 1177 874 553 1532 1802 ...
#>  $ n_dark  : int  1651 1684 1350 1267 1382 2020 1363 1652 793 958 ...
#>  $ cddr    : num  0.656 0.863 0.562 0.821 0.821 ...
#>  $ area_km2: num  93200 22000 92800 49600 23600 ...
summary(eco$cddr)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.02183 0.43932 0.57880 0.55674 0.69851 0.92447 

# metric performance against documented introductions
dn_data("auc")
#>     metric     n n_pos       auc          se        lo        hi
#> 1    DNC_q 18899  3480 0.7049489 0.004504626 0.6961200 0.7137778
#> 2  DNC_auc 18899  3480 0.7022837 0.004324280 0.6938083 0.7107591
#> 3  DNC_min 18899  3480 0.6883640 0.004339415 0.6798589 0.6968691
#> 4    NSI_q 18899  3480 0.6833160 0.004325208 0.6748388 0.6917933
#> 5      CFI 18899  3480 0.6974737 0.004457697 0.6887367 0.7062106
#> 6    CFI_w 18899  3480 0.7006914 0.004426512 0.6920156 0.7093672
#> 7      NOV 18899  3480 0.2855750 0.004264521 0.2772167 0.2939333
#> 8        D 18899  3480 0.6883598 0.004339569 0.6798544 0.6968652
#> 9        I 18899  3480 0.6946552 0.004308862 0.6862100 0.7031005
#> 10      BC 18899  3480 0.6946552 0.004308862 0.6862100 0.7031005
```
