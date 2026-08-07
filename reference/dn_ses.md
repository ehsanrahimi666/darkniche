# Standardised effect size with a floored null SD

Reference distributions for coverage metrics are heavily zero-inflated;
when a stratum piles up at zero the null SD collapses and SES values in
the hundreds are arithmetic, not evidence. The null SD is floored at
`rel_floor * scale` (2% of the between-species spread by default), which
preserves ranking and compresses only the absurd tail.

## Usage

``` r
dn_ses(obs, null_mean, null_sd, scale = NULL, rel_floor = 0.02, min_sd = 1e-09)
```

## Arguments

- obs, null_mean, null_sd:

  numeric vectors.

- scale:

  spread of the metric across all species (per region), used for the
  floor; `NULL` disables flooring.

- rel_floor:

  floor as a fraction of `scale` (default 0.02).

- min_sd:

  SDs below this become `NA` (default 1e-9).

## Value

numeric vector of SES values.

## Examples

``` r
obs <- c(0.4, 0.1, 0.8)
dn_ses(obs, null_mean = c(0.1, 0.1, 0.1), null_sd = c(0.05, 1e-6, 0.2),
       scale = 0.3)   # the 1e-6 SD is floored, not divided by
#> [1] 6.0 0.0 3.5
```
