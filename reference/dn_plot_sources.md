# Stacked source-decomposition bars (papers' Figure 11 style)

Stacked source-decomposition bars (papers' Figure 11 style)

## Usage

``` r
dn_plot_sources(
  shares,
  main = "source of dark diversity",
  legend_pos = "topright",
  cex_names = 0.7
)
```

## Arguments

- shares:

  matrix or data.frame with columns `near`, `regional`, `distant`
  (proportions summing to 1), regions in rows.

- main:

  title.

- legend_pos:

  legend position, or `NA` to suppress it.

- cex_names:

  label size for the region names.

## Value

invisibly, `shares`.

## Examples

``` r
s <- t(sapply(1:8, function(i) { x <- runif(3); x / sum(x) }))
colnames(s) <- c("near", "regional", "distant")
rownames(s) <- paste("region", 1:8)
dn_plot_sources(s)
```
