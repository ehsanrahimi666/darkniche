# Paired DeLong test between two correlated AUCs

Paired DeLong test between two correlated AUCs

## Usage

``` r
dn_delong(s1, s2, label)
```

## Arguments

- s1, s2:

  two score vectors on the same cases.

- label:

  logical or 0/1 vector.

## Value

list with `auc1`, `auc2`, `z`, `p` (two-sided).

## Examples

``` r
set.seed(1)
lab <- rep(c(TRUE, FALSE), each = 300)
s1 <- ifelse(lab, rnorm(300, 1.0), rnorm(300, 0))
s2 <- s1 + rnorm(600, 0, 0.8)      # a noisier version of the same signal
dn_delong(s1, s2, lab)
#> $auc1
#> [1] 0.7697778
#> 
#> $auc2
#> [1] 0.7110111
#> 
#> $z
#> [1] 4.094606
#> 
#> $p
#> [1] 4.228874e-05
#> 
```
