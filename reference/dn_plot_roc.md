# ROC curves for metric comparison (papers' Figure 3 style)

ROC curves for metric comparison (papers' Figure 3 style)

## Usage

``` r
dn_plot_roc(scores, label, main = "introduced-range validation")
```

## Arguments

- scores:

  named list of score vectors on the same cases.

- label:

  logical/0-1 vector.

- main:

  title.

## Value

invisibly, data.frame of AUCs with 95% CIs.

## Examples

``` r
set.seed(1)
lab <- rep(c(TRUE, FALSE), each = 300)
s1 <- ifelse(lab, rnorm(300, 1), rnorm(300, 0))
s2 <- s1 + rnorm(600, 0, 1)
dn_plot_roc(list(DNC_q = s1, D = s2), lab)
```
