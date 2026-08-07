# Summarise a null matrix

Summarise a null matrix

## Usage

``` r
dn_null_summary(obs, null, scale = NULL, rel_floor = 0.02)
```

## Arguments

- obs:

  length-k observed values.

- null:

  (k x B) null matrix.

- scale, rel_floor:

  SES flooring, see
  [`dn_ses()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_ses.md).

## Value

data.frame with null mean, sd, SES, empirical p, n valid.

## Examples

``` r
obs <- c(0.5, 0.2)
nul <- matrix(runif(2 * 99, 0, 0.4), 2, 99)
dn_null_summary(obs, nul)
#>   null_mean   null_sd        ses p_emp n_valid
#> 1 0.1943104 0.1197118 2.55354506  0.01      99
#> 2 0.1926844 0.1180602 0.06196529  0.41      99
```
