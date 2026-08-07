# Breadth-matched null (N1)

Standardises an (n_species x n_region) metric matrix against reference
species of similar niche breadth *within taxon*: for each species, the
`K` nearest neighbours in log10 breadth form the reference set. This
null asks whether a given species is exceptional relative to species of
its own niche breadth. It cannot support the region-level aggregate test
(its marginal reproduces the observed distribution by construction); use
the shift or relocation null there.

## Usage

``` r
dn_null_breadth(values, breadth, taxon, K = 200L, rel_floor = 0.02)
```

## Arguments

- values:

  (n x k) metric matrix.

- breadth:

  length-n niche breadths (any monotone measure).

- taxon:

  length-n taxon labels.

- K:

  reference window size (default 200).

- rel_floor:

  SES floor, see
  [`dn_ses()`](https://ehsanrahimi666.github.io/darkniche/reference/dn_ses.md).

## Value

list with `dn_ses`, `null_mean`, `null_sd`, `p_emp` (empirical
permutation p with floor 1/(n_ref+1)), each (n x k).

## Examples

``` r
set.seed(1)
V <- matrix(runif(30 * 3), 30, 3)          # 30 species x 3 regions
br <- rlnorm(30)
nb <- dn_null_breadth(V, br, taxon = rep("bird", 30), K = 10)
head(nb$ses)
#>            [,1]        [,2]       [,3]
#> [1,] -1.2358852 -0.20010347  1.2844100
#> [2,] -0.7254392  0.82968937 -0.5872127
#> [3,]  0.4213323  0.01948876  0.1060114
#> [4,]  1.1796359 -0.96390259 -1.5572793
#> [5,] -0.9180253  1.47554973  0.6396285
#> [6,]  1.3906353  1.13066835 -0.7073664
```
