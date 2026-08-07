# Coverage weight vector for one species

Per-cell weights whose inner product with a region density yields the
coverage statistic: the boolean q-envelope for `statistic = "q"`
(product = `DNC_q`), or `1 - dn_entryq` for `statistic = "auc"` (product
= `DNC_auc`). Both depend only on the density *ranking*, which is why
the environmental-shift null can translate a precomputed weight array
instead of re-deriving it per replicate.

## Usage

``` r
dn_weights(Z, q = 0.95, statistic = c("auc", "q"))
```

## Arguments

- Z:

  species density (matrix or vector).

- q:

  envelope level.

- statistic:

  `"q"` or `"auc"`.

## Value

numeric vector of per-cell weights.

## Examples

``` r
sp <- dn_space(cbind(rnorm(3000), rnorm(3000)), R = 50)
Z <- dn_density(cbind(rnorm(200), rnorm(200)), sp)
w_q   <- dn_weights(Z, statistic = "q")     # boolean envelope
w_auc <- dn_weights(Z, statistic = "auc")   # graded, threshold-free
cor(w_q, w_auc)
#> [1] 0.7434002
```
