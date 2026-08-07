# Envelopes, entry quantiles, and the pairwise metric block.
# Faithful port of ddarkniche.metrics; verified against golden fixtures.

# Stable descending order of a numeric vector (ties broken by index, matching
# numpy's stable argsort of -Z).
.ord_desc <- function(z) order(-z, seq_along(z), method = "radix")

#' q-quantile niche envelope
#'
#' The smallest set of lattice cells containing a fraction `q` of the density
#' mass. This replaces the naive `Z > 0` support, which after kernel smoothing
#' covers most of the lattice and makes coverage metrics bandwidth-driven
#' rather than ecology-driven.
#'
#' @param Z density as a matrix or vector (one object), or a matrix with one
#'   *row per object* if `flat = TRUE`.
#' @param q envelope level in (0, 1]; `q >= 1` returns the positive support.
#' @param flat if `TRUE`, treat `Z` as (m x P) flattened densities and return
#'   an (m x P) logical matrix.
#' @return logical vector/matrix marking envelope membership.
#' @examples
#' sp <- dn_space(cbind(rnorm(3000), rnorm(3000)), R = 50)
#' Z <- dn_density(cbind(rnorm(200), rnorm(200)), sp)
#' E <- dn_envelope(Z, q = 0.95)
#' mean(E)                                  # share of lattice in the envelope
#' @export
dn_envelope <- function(Z, q = 0.95, flat = FALSE) {
  one <- function(z) {
    tot <- sum(z)
    if (tot <= 0) return(rep(FALSE, length(z)))
    if (q >= 1) return(z > 0)
    o <- .ord_desc(z)
    cs <- cumsum(z[o]) / tot
    keep <- cs < q
    keep[1L] <- TRUE
    keep[which.max(cs >= q)] <- TRUE     # include the first crossing cell
    out <- logical(length(z)); out[o] <- keep
    out & (z > 0)
  }
  if (flat) {
    Zm <- as.matrix(Z)
    t(apply(Zm, 1L, one))
  } else if (is.matrix(Z) && !flat) {
    matrix(one(as.vector(Z)), nrow(Z), ncol(Z))
  } else one(as.numeric(Z))
}

#' Entry quantile of every lattice cell
#'
#' For each cell, the smallest envelope level q at which it enters the
#' q-envelope -- its cumulative density rank within the niche. Cells at the
#' mode enter near 0; peripheral cells enter near 1; unoccupied cells are
#' assigned 1. `1 - dn_entryq(Z)` is the weight vector whose inner
#' product with a region density gives `DNC_auc`.
#'
#' @param Z density vector (or matrix, flattened internally).
#' @return numeric vector of entry quantiles in \[0, 1].
#' @examples
#' sp <- dn_space(cbind(rnorm(3000), rnorm(3000)), R = 50)
#' Z <- dn_density(cbind(rnorm(200), rnorm(200)), sp)
#' q <- dn_entryq(Z)
#' range(q)                                 # mode enters near 0, edges near 1
#' @export
dn_entryq <- function(Z) {
  z <- as.numeric(Z)
  tot <- sum(z)
  Qv <- rep(1.0, length(z))
  if (tot > 0) {
    o <- .ord_desc(z)
    Qv[o] <- cumsum(z[o]) / tot
  }
  Qv[z <= 0] <- 1.0
  Qv
}

#' All pairwise niche metrics for a block of species and regions
#'
#' Computes, for every (species, region) pair, the symmetric baselines
#' (Schoener's D, Warren's I, the Bhattacharyya coefficient), the directional
#' coverage family (`DNC_q`, its threshold-free integral `DNC_auc`, the
#' converse share `NSI_q`), the min-overlap identities (`DNC_min`, `NSI_min`;
#' with unit-mass densities `DNC_min == D`, the algebraic result that no
#' symmetric functional of the min family can express containment), the
#' novelty share `NOV`, and -- when moments are supplied -- the Mahalanobis
#' core-fit index `CFI` and its availability-weighted form `CFI_w`.
#'
#' @param Zs (m x P) matrix: one flattened species density per row.
#' @param Ze (k x P) matrix: one flattened region density per row.
#' @param q envelope level for the threshold form (default 0.95).
#' @param sp_mu,sp_cov optional species moments: (m x 2) means and a list of
#'   m 2x2 covariance matrices (or an (m x 4) matrix, row-major).
#' @param eco_mu optional (k x 2) region centroids.
#' @param centres optional `dn_centres(sp)` list, required for `CFI_w`.
#' @return named list of (m x k) matrices.
#' @examples
#' set.seed(1)
#' sp <- dn_space(cbind(rnorm(5000, 0, 2), rnorm(5000, 0, 2)), R = 60)
#' # broad species niche engulfing a small region: the containment case
#' broad  <- dn_density(cbind(rnorm(1500, 0, 1.5), rnorm(1500, 0, 1.5)), sp)
#' region <- dn_density(cbind(rnorm(400, 0, 0.25), rnorm(400, 0, 0.25)), sp)
#' m <- dn_metrics(matrix(as.numeric(broad), 1), matrix(as.numeric(region), 1))
#' # D is near zero although the region is trivially habitable:
#' round(c(D = m$D, DNC_q = m$DNC_q, DNC_auc = m$DNC_auc), 3)
#' @export
dn_metrics <- function(Zs, Ze, q = 0.95, sp_mu = NULL, sp_cov = NULL,
                         eco_mu = NULL, centres = NULL) {
  Zs <- as.matrix(Zs); Ze <- as.matrix(Ze)
  m <- nrow(Zs); k <- nrow(Ze); P <- ncol(Zs)
  stopifnot(ncol(Ze) == P)
  out <- list()

  # symmetric baselines
  MIN <- matrix(0, m, k); L1 <- matrix(0, m, k)
  ss <- rowSums(Zs); se <- rowSums(Ze)
  for (j in seq_len(k)) {
    d <- abs(sweep(Zs, 2L, Ze[j, ], "-"))
    L1[, j] <- rowSums(d)
    MIN[, j] <- 0.5 * (ss + se[j] - L1[, j])
  }
  out$D <- 1 - 0.5 * L1
  BC <- sqrt(Zs) %*% t(sqrt(Ze))
  out$BC <- BC
  out$I <- 1 - 0.5 * pmax(2 - 2 * BC, 0)

  # directional coverage
  Es <- dn_envelope(Zs, q, flat = TRUE)
  Ee <- dn_envelope(Ze, q, flat = TRUE)
  out$DNC_q <- (Es + 0) %*% t(Ze)
  out$NSI_q <- Zs %*% t(Ee + 0)
  W <- 1 - t(apply(Zs, 1L, dn_entryq))
  out$DNC_auc <- W %*% t(Ze)

  out$DNC_min <- sweep(MIN, 2L, ifelse(se > 0, se, NA_real_), "/")
  out$NSI_min <- MIN / ifelse(ss > 0, ss, NA_real_)

  Es_full <- dn_envelope(Zs, 1.0, flat = TRUE)
  out$NOV <- 1 - (Es_full + 0) %*% t(Ze)

  # core fit
  na_mk <- matrix(NA_real_, m, k)
  out$mahal_d2 <- na_mk; out$CFI <- na_mk; out$CFI_w <- na_mk
  if (!is.null(sp_mu) && !is.null(sp_cov) && !is.null(eco_mu)) {
    covs <- if (is.list(sp_cov)) sp_cov else
      lapply(seq_len(m), function(i) matrix(sp_cov[i, ], 2, 2, byrow = TRUE))
    d2 <- matrix(NA_real_, m, k)
    for (i in seq_len(m)) {
      Ci <- tryCatch(solve(covs[[i]] + 1e-12 * diag(2)),
                     error = function(e) NULL)
      if (is.null(Ci)) next
      dv <- sweep(as.matrix(eco_mu), 2L, sp_mu[i, ], "-")
      d2[i, ] <- rowSums((dv %*% Ci) * dv)
    }
    out$mahal_d2 <- d2
    out$CFI <- stats::pchisq(d2, df = 2, lower.tail = FALSE)
    if (!is.null(centres)) {
      G <- as.matrix(expand.grid(x = centres$x, y = centres$y))
      # expand.grid varies x fastest; Python meshgrid(indexing="ij") varies the
      # second axis fastest with row-major ravel. Both orderings enumerate the
      # same P cells; what matters is agreement with the density flattening.
      # Densities here are flattened R-style (column-major on the R x R matrix,
      # x index fastest), matching expand.grid.
      CFIw <- matrix(NA_real_, m, k)
      for (i in seq_len(m)) {
        Ci <- tryCatch(solve(covs[[i]] + 1e-12 * diag(2)),
                       error = function(e) NULL)
        if (is.null(Ci)) next
        dv <- sweep(G, 2L, sp_mu[i, ], "-")
        wcell <- stats::pchisq(rowSums((dv %*% Ci) * dv), df = 2,
                               lower.tail = FALSE)
        CFIw[i, ] <- Ze %*% wcell
      }
      out$CFI_w <- CFIw
    }
  }
  out
}

#' Coverage weight vector for one species
#'
#' Per-cell weights whose inner product with a region density yields the
#' coverage statistic: the boolean q-envelope for `statistic = "q"`
#' (product = `DNC_q`), or `1 - dn_entryq` for `statistic = "auc"`
#' (product = `DNC_auc`). Both depend only on the density *ranking*, which is
#' why the environmental-shift null can translate a precomputed weight array
#' instead of re-deriving it per replicate.
#'
#' @param Z species density (matrix or vector).
#' @param q envelope level.
#' @param statistic `"q"` or `"auc"`.
#' @return numeric vector of per-cell weights.
#' @examples
#' sp <- dn_space(cbind(rnorm(3000), rnorm(3000)), R = 50)
#' Z <- dn_density(cbind(rnorm(200), rnorm(200)), sp)
#' w_q   <- dn_weights(Z, statistic = "q")     # boolean envelope
#' w_auc <- dn_weights(Z, statistic = "auc")   # graded, threshold-free
#' cor(w_q, w_auc)
#' @export
dn_weights <- function(Z, q = 0.95, statistic = c("auc", "q")) {
  statistic <- match.arg(statistic)
  if (statistic == "q") as.numeric(dn_envelope(as.numeric(Z), q))
  else 1 - dn_entryq(Z)
}
