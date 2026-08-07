# Calibrate the admissibility threshold (Youden)

Calibrate the admissibility threshold (Youden)

## Usage

``` r
dn_calibrate(score, label)
```

## Arguments

- score:

  establishment-model probabilities (or any coverage score).

- label:

  validation labels.

## Value

the Youden-optimal threshold.

## Examples

``` r
set.seed(1)
lab <- rep(c(TRUE, FALSE), each = 200)
p <- ifelse(lab, rbeta(200, 4, 2), rbeta(200, 2, 4))
dn_calibrate(p, lab)     # Youden-optimal threshold
#> [1] 0.5210077
```
