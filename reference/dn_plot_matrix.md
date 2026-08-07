# Region-by-taxon matrix heatmap (papers' Figure 9 style)

Region-by-taxon matrix heatmap (papers' Figure 9 style)

## Usage

``` r
dn_plot_matrix(M, main = "", values = TRUE, key = TRUE)
```

## Arguments

- M:

  numeric matrix, regions in rows, taxa in columns (e.g. CDDR).

- main:

  title.

- values:

  logical, print cell values (default TRUE).

- key:

  logical, draw a colour key in the right margin (default TRUE).

## Value

invisibly, `M`.

## Examples

``` r
M <- matrix(runif(6 * 5), 6, 5,
            dimnames = list(paste("region", 1:6),
                            c("amph", "birds", "mamm", "plants", "rept")))
dn_plot_matrix(M, main = "CDDR by region and taxon")
```
