# Environmental space: lattice definition and kernel density estimation.
#
# Numerical contract: every function in this file reproduces the Python
# reference implementation (ddarkniche.niche) to float32 precision, verified by
# the golden tests in tests/testthat/test-golden.R.

#' Define a shared environmental-space lattice
#'
#' Builds the common `R x R` lattice on which all species niches and all
#' region climate availabilities are estimated, from the global distribution
#' of environmental scores (typically the first two principal components of
#' climate over all terrestrial cells).
#'
#' @param scores numeric matrix (n x 2) of environmental scores for the
#'   *global* background (all terrestrial grid cells), not for any one
#'   species or region. The lattice must be shared or no two densities are
#'   comparable.
#' @param R lattice resolution per axis (default 100, as calibrated in the
#'   source studies; sensitivity there: Spearman >= 0.996 across 60-150).
#' @param clip_q lower/upper quantiles at which the axes are clipped before
#'   5% padding (default `c(0.001, 0.999)`).
#' @param bw_frac fixed kernel bandwidth as a fraction of the axis range
#'   (default 0.02; one bandwidth for every object so that estimated niche
#'   breadth reflects ecology, not sample size).
#' @return an object of class `dn_space`: list with `x_edges`, `y_edges`,
#'   `R`, `bw_cells`.
#' @examples
#' set.seed(1)
#' bg <- cbind(rnorm(5000), rnorm(5000))   # global background scores
#' sp <- dn_space(bg, R = 60)
#' sp
#' @export
dn_space <- function(scores, R = 100L, clip_q = c(0.001, 0.999),
                     bw_frac = 0.02) {
  s <- as.matrix(scores)[, 1:2, drop = FALSE]
  storage.mode(s) <- "double"
  xq <- stats::quantile(s[, 1], clip_q, names = FALSE, type = 7)
  yq <- stats::quantile(s[, 2], clip_q, names = FALSE, type = 7)
  padx <- 0.05 * (xq[2] - xq[1]); pady <- 0.05 * (yq[2] - yq[1])
  sp <- list(
    x_edges = seq(xq[1] - padx, xq[2] + padx, length.out = R + 1L),
    y_edges = seq(yq[1] - pady, yq[2] + pady, length.out = R + 1L),
    R = as.integer(R),
    bw_cells = bw_frac * R
  )
  class(sp) <- "dn_space"
  sp
}

#' @export
print.dn_space <- function(x, ...) {
  cat(sprintf("<dn_space> %d x %d lattice, bandwidth %.2f cells\n",
              x$R, x$R, x$bw_cells))
  cat(sprintf("  PC1 [%.3f, %.3f]  PC2 [%.3f, %.3f]\n",
              min(x$x_edges), max(x$x_edges),
              min(x$y_edges), max(x$y_edges)))
  invisible(x)
}

#' Lattice cell centres
#' @param sp a `dn_space`
#' @return list with `x`, `y` centre vectors
#' @examples
#' sp <- dn_space(cbind(rnorm(1000), rnorm(1000)), R = 40)
#' ctr <- dn_centres(sp)
#' str(ctr)
#' @export
dn_centres <- function(sp) {
  list(x = (sp$x_edges[-1] + sp$x_edges[-length(sp$x_edges)]) / 2,
       y = (sp$y_edges[-1] + sp$y_edges[-length(sp$y_edges)]) / 2)
}

# 1-D Gaussian kernel identical to scipy.ndimage.gaussian_filter1d:
# radius = floor(truncate * sigma + 0.5), weights exp(-i^2 / 2 sigma^2)
# normalised to sum 1.
.gauss_kernel <- function(sigma, truncate = 4.0) {
  r <- as.integer(truncate * sigma + 0.5)
  x <- (-r):r
  w <- exp(-(x^2) / (2 * sigma^2))
  w / sum(w)
}

# Separable Gaussian blur with zero (constant) boundary, matching
# scipy.ndimage.gaussian_filter(mode = "constant", cval = 0, truncate = 4).
.gauss_blur <- function(H, sigma) {
  w <- .gauss_kernel(sigma)
  r <- (length(w) - 1L) / 2L
  n <- nrow(H); m <- ncol(H)
  # correlate along rows (axis 1) then columns (axis 2); kernel symmetric so
  # correlation == convolution. Zero-pad, convolve, crop.
  conv1 <- function(M) {                       # blur down each column
    P <- rbind(matrix(0, r, ncol(M)), M, matrix(0, r, ncol(M)))
    out <- matrix(0, nrow(M), ncol(M))
    for (k in seq_along(w))
      out <- out + w[k] * P[(k - 1L) + seq_len(nrow(M)), , drop = FALSE]
    out
  }
  t(conv1(t(conv1(H))))
}

#' Kernel density of a point set on the lattice
#'
#' Bins the points on the lattice (half-open bins, right edge of the last bin
#' closed, exactly as `numpy.histogram2d`), smooths with a Gaussian kernel of
#' `sp$bw_cells` cells, and normalises to unit mass.
#'
#' @param points numeric matrix (n x 2) of scores (e.g. a species' occupied
#'   cells projected into the space, or a region's cells).
#' @param sp a `dn_space`.
#' @param weights optional per-point weights (e.g. cell areas).
#' @return numeric `R x R` matrix summing to 1 (all-zero if no point falls
#'   inside the lattice).
#' @examples
#' set.seed(1)
#' sp <- dn_space(cbind(rnorm(5000), rnorm(5000)), R = 60)
#' pts <- cbind(rnorm(300, 0.5, 0.4), rnorm(300, -0.2, 0.5))
#' Z <- dn_density(pts, sp)
#' sum(Z)                                   # unit mass
#' image(Z, main = "species niche density")
#' @export
dn_density <- function(points, sp, weights = NULL) {
  Z0 <- matrix(0, sp$R, sp$R)
  if (is.null(points) || NROW(points) == 0L) return(Z0)
  p <- as.matrix(points); storage.mode(p) <- "double"
  xe <- sp$x_edges; ye <- sp$y_edges
  ix <- findInterval(p[, 1], xe, rightmost.closed = TRUE)
  iy <- findInterval(p[, 2], ye, rightmost.closed = TRUE)
  keep <- ix >= 1L & ix <= sp$R & iy >= 1L & iy <= sp$R
  if (!any(keep)) return(Z0)
  w <- if (is.null(weights)) rep(1, sum(keep)) else as.numeric(weights)[keep]
  H <- Z0
  idx <- cbind(ix[keep], iy[keep])
  for (i in seq_len(nrow(idx))) H[idx[i, 1], idx[i, 2]] <- H[idx[i, 1], idx[i, 2]] + w[i]
  if (sum(H) == 0) return(Z0)
  Z <- .gauss_blur(H, sp$bw_cells)
  s <- sum(Z)
  if (s <= 0) Z0 else Z / s
}

#' First and second moments of a lattice density
#'
#' Mean vector and covariance matrix of a density in environmental space,
#' used by the core-fit index (CFI).
#'
#' @param Z `R x R` density matrix.
#' @param sp a `dn_space`.
#' @return list with `mu` (length 2) and `cov` (2 x 2).
#' @examples
#' sp <- dn_space(cbind(rnorm(3000), rnorm(3000)), R = 50)
#' Z <- dn_density(cbind(rnorm(200, 1, 0.3), rnorm(200, 0, 0.4)), sp)
#' dn_moments(Z, sp)
#' @export
dn_moments <- function(Z, sp) {
  ctr <- dn_centres(sp)
  w <- Z / sum(Z)
  wx <- rowSums(w); wy <- colSums(w)
  mu <- c(sum(ctr$x * wx), sum(ctr$y * wy))
  X <- outer(ctr$x - mu[1], rep(1, sp$R))
  Y <- outer(rep(1, sp$R), ctr$y - mu[2])
  cxx <- sum(w * X * X); cyy <- sum(w * Y * Y); cxy <- sum(w * X * Y)
  list(mu = mu, cov = matrix(c(cxx, cxy, cxy, cyy), 2, 2))
}
