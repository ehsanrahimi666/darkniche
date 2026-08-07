# Region-level aggregate admissibility test (single threshold)

Tests whether a region holds more climatically admissible species than
the null expects: the observed count of species whose coverage exceeds
`threshold` against the Poisson-binomial expectation of the same count
under each species' own null. One rule, both sides – applying a stricter
criterion to the observed side than the null side guarantees a negative
result regardless of the data.

## Usage

``` r
dn_aggregate(obs, null_mean, null_sd, threshold)
```

## Arguments

- obs:

  length-n observed coverage values (one per species).

- null_mean, null_sd:

  length-n per-species null moments.

- threshold:

  coverage threshold.

## Value

one-row data.frame: `n_species`, `n_admissible_obs`, `expected`, `sd`,
`dn_ses`, `p`, `enrichment`, `reliable`.

## Details

Only nulls whose marginals are genuinely counterfactual support this
test (the environmental-shift null N2 or relocation null N3); the
breadth-matched null reproduces the observed marginal by construction
and has no power here.

Results with expected count \< 5 are flagged `reliable = FALSE`: the
normal approximation to the Poisson-binomial count is invalid there and
the SES/enrichment are arithmetic, not evidence. **Never report
unreliable cells.**

## Examples

``` r
set.seed(1)
obs <- c(rbeta(300, 0.3, 3), rbeta(100, 4, 1.5))   # bimodal coverage
mu <- rbeta(400, 1, 3); sdv <- 0.05 + 0.1 * runif(400)
dn_aggregate(obs, mu, sdv, threshold = quantile(obs, 0.95))
#>   n_species n_admissible_obs expected        sd ses  p enrichment reliable
#> 1       400               20 1.162139 0.9684571  NA NA         NA    FALSE
```
