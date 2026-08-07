# Null models: breadth-matched (N1), environmental-shift (N2), geographic
# relocation (N3), and the SES machinery. Ports of ddarkniche.nulls / utils.

#' Standardised effect size with a floored null SD
#'
#' Reference distributions for coverage metrics are heavily zero-inflated;
#' when a stratum piles up at zero the null SD collapses and SES values in
#' the hundreds are arithmetic, not evidence. The null SD is floored at
#' `rel_floor * scale` (2% of the between-species spread by default), which
#' preserves ranking and compresses only the absurd tail.
#'
#' @param obs,null_mean,null_sd numeric vectors.
#' @param scale spread of the metric across all species (per region), used
#'   for the floor; `NULL` disables flooring.
#' @param rel_floor floor as a fraction of `scale` (default 0.02).
#' @param min_sd SDs below this become `NA` (default 1e-9).
#' @return numeric vector of SES values.
#' @examples
#' obs <- c(0.4, 0.1, 0.8)
#' dn_ses(obs, null_mean = c(0.1, 0.1, 0.1), null_sd = c(0.05, 1e-6, 0.2),
#'        scale = 0.3)   # the 1e-6 SD is floored, not divided by
#' @export
dn_ses <- function(obs, null_mean, null_sd, scale = NULL, rel_floor = 0.02,
                min_sd = 1e-9) {
  sd <- as.numeric(null_sd)
  if (!is.null(scale)) sd <- pmax(sd, rel_floor * as.numeric(scale))
  sd[sd < min_sd] <- NA_real_
  (as.numeric(obs) - as.numeric(null_mean)) / sd
}

#' Breadth-matched null (N1)
#'
#' Standardises an (n_species x n_region) metric matrix against reference
#' species of similar niche breadth *within taxon*: for each species, the
#' `K` nearest neighbours in log10 breadth form the reference set. This null
#' asks whether a given species is exceptional relative to species of its own
#' niche breadth. It cannot support the region-level aggregate test (its
#' marginal reproduces the observed distribution by construction); use the
#' shift or relocation null there.
#'
#' @param values (n x k) metric matrix.
#' @param breadth length-n niche breadths (any monotone measure).
#' @param taxon length-n taxon labels.
#' @param K reference window size (default 200).
#' @param rel_floor SES floor, see [dn_ses()].
#' @return list with `dn_ses`, `null_mean`, `null_sd`, `p_emp` (empirical
#'   permutation p with floor 1/(n_ref+1)), each (n x k).
#' @examples
#' set.seed(1)
#' V <- matrix(runif(30 * 3), 30, 3)          # 30 species x 3 regions
#' br <- rlnorm(30)
#' nb <- dn_null_breadth(V, br, taxon = rep("bird", 30), K = 10)
#' head(nb$ses)
#' @export
dn_null_breadth <- function(values, breadth, taxon, K = 200L, rel_floor = 0.02) {
  V <- as.matrix(values); n <- nrow(V); k <- ncol(V)
  lb <- log10(as.numeric(breadth) + 1e-12)
  col_scale <- apply(V, 2L, stats::sd, na.rm = TRUE)
  SES <- MU <- SD <- PE <- matrix(NA_real_, n, k)
  for (tx in unique(taxon)) {
    idx <- which(taxon == tx)
    if (length(idx) < 10L) next
    ord <- idx[order(lb[idx])]
    m <- length(ord); half <- max(5L, K %/% 2L)
    for (pos in seq_len(m)) {
      lo <- max(1L, pos - half); hi <- min(m, pos + half)
      sel <- setdiff(lo:hi, pos)
      if (length(sel) < 5L) next
      ref <- V[ord[sel], , drop = FALSE]
      mu <- colMeans(ref, na.rm = TRUE)
      sd <- apply(ref, 2L, stats::sd, na.rm = TRUE)
      i <- ord[pos]
      MU[i, ] <- mu; SD[i, ] <- sd
      SES[i, ] <- dn_ses(V[i, ], mu, sd, scale = col_scale,
                      rel_floor = rel_floor)
      ge <- colSums(sweep(ref, 2L, V[i, ], ">="), na.rm = TRUE)
      PE[i, ] <- pmax(ge, 1) / (nrow(ref) + 1)
    }
  }
  list(ses = SES, null_mean = MU, null_sd = SD, p_emp = PE)
}

#' Environmental-shift null (N2)
#'
#' Null distribution of directional coverage obtained by translating the
#' species' coverage-weight array to random, availability-weighted positions
#' in environmental space. Translation permutes lattice cells without
#' altering density values, so the envelope of a shifted niche is exactly the
#' shifted envelope; the weights are computed once and rolled.
#'
#' @param Z_sp `R x R` species density.
#' @param Ze (k x P) region densities.
#' @param Z_global `R x R` global availability density; shift targets are
#'   drawn proportional to it.
#' @param q envelope level.
#' @param B number of replicates (default 499).
#' @param statistic `"q"` or `"auc"`.
#' @param targets optional integer vector of length B giving flat target cell
#'   indices (1-based, column-major on the R x R lattice); supply to
#'   reproduce a specific draw exactly (e.g. from the Python pipeline).
#' @param seed RNG seed used when `targets` is NULL.
#' @return (k x B) matrix of null coverage values.
#' @examples
#' set.seed(1)
#' bg <- cbind(rnorm(4000), rnorm(4000))
#' sp <- dn_space(bg, R = 50)
#' Zg <- dn_density(bg, sp)
#' Zs <- dn_density(cbind(rnorm(300, 1, 0.4), rnorm(300, 0, 0.4)), sp)
#' Ze <- matrix(as.numeric(dn_density(cbind(rnorm(300, 1, 0.3),
#'                                          rnorm(300, 0, 0.3)), sp)), 1)
#' obs <- as.numeric(Ze %*% dn_weights(Zs))
#' nul <- dn_null_shift(Zs, Ze, Zg, B = 99, seed = 1)
#' dn_null_summary(obs, nul)
#' @export
dn_null_shift <- function(Z_sp, Ze, Z_global, q = 0.95, B = 499L,
                       statistic = c("auc", "q"), targets = NULL,
                       seed = NULL) {
  statistic <- match.arg(statistic)
  Rn <- nrow(Z_sp)
  w <- as.numeric(Z_global); w <- if (sum(w) > 0) w / sum(w) else
    rep(1 / length(w), length(w))
  if (is.null(targets)) {
    if (!is.null(seed)) set.seed(seed)
    targets <- sample.int(length(w), B, replace = TRUE, prob = w)
  }
  ti <- ((targets - 1L) %% Rn) + 1L          # row index, column-major
  tj <- ((targets - 1L) %/% Rn) + 1L
  cc <- .centroid_cell(Z_sp)
  E0 <- matrix(dn_weights(Z_sp, q, statistic), Rn, Rn)
  Wm <- matrix(0, Rn * Rn, length(targets))
  for (b in seq_along(targets)) {
    Wm[, b] <- as.numeric(.roll2(E0, ti[b] - cc[1], tj[b] - cc[2]))
  }
  as.matrix(Ze) %*% Wm
}

.centroid_cell <- function(Z) {
  w <- Z / sum(Z)
  Rn <- nrow(Z)
  i <- sum(seq_len(Rn) * rowSums(w))
  j <- sum(seq_len(Rn) * colSums(w))
  c(round(i), round(j))
}

.roll2 <- function(M, di, dj) {
  n <- nrow(M); m <- ncol(M)
  di <- ((di %% n) + n) %% n; dj <- ((dj %% m) + m) %% m
  if (di > 0) M <- rbind(M[(n - di + 1):n, , drop = FALSE],
                         M[1:(n - di), , drop = FALSE])
  if (dj > 0) M <- cbind(M[, (m - dj + 1):m, drop = FALSE],
                         M[, 1:(m - dj), drop = FALSE])
  M
}

#' Geographic relocation null (N3)
#'
#' Null coverage obtained by stamping the species' occupied-cell footprint at
#' random positions on the geographic grid and recomputing coverage from the
#' relocated range. Requires the geographic cell adjacency of the analysis
#' grid; supply `project` to map relocated geographic cells to environmental
#' scores.
#'
#' @param cells integer matrix (n x 2) of occupied geographic cell indices.
#' @param grid_dim geographic grid dimensions `c(nrow, ncol)`.
#' @param valid logical matrix marking terrestrial cells.
#' @param project function taking an (n x 2) cell-index matrix and returning
#'   (n x 2) environmental scores (NA rows for cells without climate).
#' @param sp a `dn_space`.
#' @param Ze (k x P) region densities.
#' @param q,statistic,B,seed as in [dn_null_shift()].
#' @return (k x B) matrix of null coverage values (NA where a relocation
#'   could not be placed).
#' @examples
#' # Requires a geographic grid and a projection function; sketch:
#' \dontrun{
#' nul <- dn_null_relocate(cells, grid_dim = dim(valid), valid = valid,
#'                         project = my_lookup, sp = sp, Ze = Ze, B = 99)
#' }
#' @export
dn_null_relocate <- function(cells, grid_dim, valid, project, sp, Ze,
                          q = 0.95, statistic = c("auc", "q"),
                          B = 199L, seed = 1L) {
  statistic <- match.arg(statistic)
  set.seed(seed)
  cells <- as.matrix(cells)
  ctr <- round(colMeans(cells))
  offs <- sweep(cells, 2L, ctr, "-")
  cand <- which(valid, arr.ind = TRUE)
  out <- matrix(NA_real_, nrow(Ze), B)
  for (b in seq_len(B)) {
    tgt <- cand[sample.int(nrow(cand), 1L), ]
    newc <- sweep(offs, 2L, tgt, "+")
    ok <- newc[, 1] >= 1 & newc[, 1] <= grid_dim[1] &
      newc[, 2] >= 1 & newc[, 2] <= grid_dim[2]
    if (mean(ok) < 0.5) next
    newc <- newc[ok, , drop = FALSE]
    inl <- valid[newc]
    if (mean(inl) < 0.5) next
    scr <- project(newc[inl, , drop = FALSE])
    scr <- scr[stats::complete.cases(scr), , drop = FALSE]
    if (nrow(scr) < 5L) next
    Zb <- dn_density(scr, sp)
    wv <- dn_weights(Zb, q, statistic)
    out[, b] <- as.matrix(Ze) %*% wv
  }
  out
}

#' Summarise a null matrix
#'
#' @param obs length-k observed values.
#' @param null (k x B) null matrix.
#' @param scale,rel_floor SES flooring, see [dn_ses()].
#' @return data.frame with null mean, sd, SES, empirical p, n valid.
#' @examples
#' obs <- c(0.5, 0.2)
#' nul <- matrix(runif(2 * 99, 0, 0.4), 2, 99)
#' dn_null_summary(obs, nul)
#' @export
dn_null_summary <- function(obs, null, scale = NULL, rel_floor = 0.02) {
  mu <- rowMeans(null, na.rm = TRUE)
  sd <- apply(null, 1L, stats::sd, na.rm = TRUE)
  nv <- rowSums(is.finite(null))
  ge <- rowSums(sweep(null, 1L, obs, ">="), na.rm = TRUE)
  data.frame(null_mean = mu, null_sd = sd,
             ses = dn_ses(obs, mu, sd, scale = scale, rel_floor = rel_floor),
             p_emp = pmax(ge, 1) / (nv + 1), n_valid = nv)
}
