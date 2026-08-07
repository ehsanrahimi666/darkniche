# All pairwise niche metrics for a block of species and regions

Computes, for every (species, region) pair, the symmetric baselines
(Schoener's D, Warren's I, the Bhattacharyya coefficient), the
directional coverage family (`DNC_q`, its threshold-free integral
`DNC_auc`, the converse share `NSI_q`), the min-overlap identities
(`DNC_min`, `NSI_min`; with unit-mass densities `DNC_min == D`, the
algebraic result that no symmetric functional of the min family can
express containment), the novelty share `NOV`, and – when moments are
supplied – the Mahalanobis core-fit index `CFI` and its
availability-weighted form `CFI_w`.

## Usage

``` r
dn_metrics(
  Zs,
  Ze,
  q = 0.95,
  sp_mu = NULL,
  sp_cov = NULL,
  eco_mu = NULL,
  centres = NULL
)
```

## Arguments

- Zs:

  (m x P) matrix: one flattened species density per row.

- Ze:

  (k x P) matrix: one flattened region density per row.

- q:

  envelope level for the threshold form (default 0.95).

- sp_mu, sp_cov:

  optional species moments: (m x 2) means and a list of m 2x2 covariance
  matrices (or an (m x 4) matrix, row-major).

- eco_mu:

  optional (k x 2) region centroids.

- centres:

  optional `dn_centres(sp)` list, required for `CFI_w`.

## Value

named list of (m x k) matrices.

## Examples

``` r
set.seed(1)
sp <- dn_space(cbind(rnorm(5000, 0, 2), rnorm(5000, 0, 2)), R = 60)
# broad species niche engulfing a small region: the containment case
broad  <- dn_density(cbind(rnorm(1500, 0, 1.5), rnorm(1500, 0, 1.5)), sp)
region <- dn_density(cbind(rnorm(400, 0, 0.25), rnorm(400, 0, 0.25)), sp)
m <- dn_metrics(matrix(as.numeric(broad), 1), matrix(as.numeric(region), 1))
# D is near zero although the region is trivially habitable:
round(c(D = m$D, DNC_q = m$DNC_q, DNC_auc = m$DNC_auc), 3)
#>       D   DNC_q DNC_auc 
#>   0.205   1.000   0.860 
```
