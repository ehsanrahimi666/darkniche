# Dark-diversity classification of species-region pairs

Applies the full conjunction rule of the source studies: a pair is
*climatically admissible* if the calibrated establishment probability
exceeds `p_min` and every supplied null screen has SES \>= `ses_min`; it
is *dark* if additionally the species is absent from the region.

## Usage

``` r
dn_dark(pairs, p_min = 0.177, ses_cols = character(), ses_min = 2)
```

## Arguments

- pairs:

  data.frame with at least `species`, `region`, `present` (logical), a
  probability column `p_establish`, and one or more SES columns given in
  `ses_cols`.

- p_min:

  admissibility threshold on the establishment probability. The source
  studies' externally calibrated value (Youden threshold of a
  species-blocked, climate-only model validated on documented introduced
  ranges) is 0.177; recalibrate with
  [`dn_calibrate()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_calibrate.md)
  when using a different model or metric.

- ses_cols:

  character vector of SES column names (e.g.
  `c("n1_ses", "n2_ses", "n3_ses")`).

- ses_min:

  SES criterion (default 2).

## Value

`pairs` with logical columns `admissible` and `dark`.

## Examples

``` r
pairs <- data.frame(species = c("a", "a", "b"), region = c("r1", "r2", "r1"),
                    present = c(TRUE, FALSE, FALSE),
                    p_establish = c(0.6, 0.4, 0.05),
                    n2_ses = c(3.1, 2.5, 4.0))
dn_dark(pairs, ses_cols = "n2_ses")
#>   species region present p_establish n2_ses admissible  dark
#> 1       a     r1    TRUE        0.60    3.1       TRUE FALSE
#> 2       a     r2   FALSE        0.40    2.5       TRUE  TRUE
#> 3       b     r1   FALSE        0.05    4.0      FALSE FALSE
```
