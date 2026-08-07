# Region-level aggregate test, the threshold-sweep diagnostic, dark-diversity
# accounting, and introduced-range validation helpers.

#' Region-level aggregate admissibility test (single threshold)
#'
#' Tests whether a region holds more climatically admissible species than the
#' null expects: the observed count of species whose coverage exceeds
#' `threshold` against the Poisson-binomial expectation of the same count
#' under each species' own null. One rule, both sides -- applying a stricter
#' criterion to the observed side than the null side guarantees a negative
#' result regardless of the data.
#'
#' Only nulls whose marginals are genuinely counterfactual support this test
#' (the environmental-shift null N2 or relocation null N3); the
#' breadth-matched null reproduces the observed marginal by construction and
#' has no power here.
#'
#' Results with expected count < 5 are flagged `reliable = FALSE`: the normal
#' approximation to the Poisson-binomial count is invalid there and the
#' SES/enrichment are arithmetic, not evidence. **Never report unreliable
#' cells.**
#'
#' @param obs length-n observed coverage values (one per species).
#' @param null_mean,null_sd length-n per-species null moments.
#' @param threshold coverage threshold.
#' @return one-row data.frame: `n_species`, `n_admissible_obs`, `expected`,
#'   `sd`, `dn_ses`, `p`, `enrichment`, `reliable`.
#' @examples
#' set.seed(1)
#' obs <- c(rbeta(300, 0.3, 3), rbeta(100, 4, 1.5))   # bimodal coverage
#' mu <- rbeta(400, 1, 3); sdv <- 0.05 + 0.1 * runif(400)
#' dn_aggregate(obs, mu, sdv, threshold = quantile(obs, 0.95))
#' @export
dn_aggregate <- function(obs, null_mean, null_sd, threshold) {
  obs <- as.numeric(obs); mu <- as.numeric(null_mean); sd <- as.numeric(null_sd)
  ok <- is.finite(obs) & is.finite(mu) & is.finite(sd) & sd > 0
  if (sum(ok) < 30L)
    return(data.frame(n_species = sum(ok), n_admissible_obs = NA_integer_,
                      expected = NA_real_, sd = NA_real_, ses = NA_real_,
                      p = NA_real_, enrichment = NA_real_, reliable = FALSE))
  n_obs <- sum(obs[ok] >= threshold)
  p_i <- stats::pnorm((threshold - mu[ok]) / sd[ok], lower.tail = FALSE)
  expd <- sum(p_i); v <- sum(p_i * (1 - p_i))
  z <- if (v > 0) (n_obs - expd) / sqrt(v) else NA_real_
  reliable <- expd >= 5
  data.frame(
    n_species = sum(ok), n_admissible_obs = n_obs, expected = expd,
    sd = sqrt(v),
    ses = if (is.finite(z) && reliable) z else NA_real_,
    p = if (is.finite(z) && reliable)
      stats::pnorm(z, lower.tail = FALSE) else NA_real_,
    enrichment = if (expd > 0 && reliable) n_obs / expd else NA_real_,
    reliable = reliable)
}

#' Threshold-sweep diagnostic for the aggregate null test
#'
#' The central methodological caution of the source studies: the direction of
#' the single-threshold aggregate test is, in general, a property of the
#' threshold, not of the region -- in the global analysis its verdict
#' reversed with the threshold in 89% of the world's ecoregions, because
#' observed coverage is bimodal while null coverage is unimodal, so the two
#' distributions must cross. Any aggregate null test on a thresholded,
#' bimodally distributed similarity statistic should either demonstrate that
#' the sign of its result is stable across a stated threshold range, or
#' report the full sweep.
#'
#' Thresholds are placed at quantiles of the pooled observed coverage
#' distribution, so the sweep is comparable across regions.
#'
#' @param obs length-n observed coverage values for one region.
#' @param null_mean,null_sd per-species null moments (as [dn_aggregate()]).
#' @param quantiles threshold quantiles of the observed distribution
#'   (default `seq(0.80, 0.99, by = 0.01)`).
#' @param pooled_obs optional vector from which threshold quantiles are
#'   computed (e.g. all regions pooled); defaults to `obs`.
#' @return object of class `dn_sweep_result`: data.frame with one row per
#'   threshold plus attributes `direction_change` (logical) and `monotone`.
#' @examples
#' set.seed(1)
#' obs <- c(rbeta(300, 0.3, 3), rbeta(100, 4, 1.5))   # bimodal, like real data
#' mu <- rbeta(400, 1, 2.2); sdv <- 0.08 + 0.1 * runif(400)
#' sw <- dn_sweep(obs, mu, sdv)
#' sw          # prints a warning if the verdict reverses across thresholds
#' plot(sw)
#' @export
dn_sweep <- function(obs, null_mean, null_sd,
                            quantiles = seq(0.80, 0.99, by = 0.01),
                            pooled_obs = NULL) {
  src <- if (is.null(pooled_obs)) obs else pooled_obs
  thr <- stats::quantile(src, quantiles, names = FALSE, na.rm = TRUE)
  rows <- lapply(seq_along(thr), function(i) {
    r <- dn_aggregate(obs, null_mean, null_sd, thr[i])
    r$threshold_q <- quantiles[i]; r$threshold <- thr[i]; r
  })
  out <- do.call(rbind, rows)
  ratio <- out$n_admissible_obs / out$expected
  # Only reliable thresholds may contribute to the verdict. An expected count
  # below 5 makes the ratio arithmetic rather than evidence, so including such
  # rows here would let the diagnostic report exactly the numbers the test
  # itself refuses to report.
  ok <- out$reliable & is.finite(ratio)
  attr(out, "ratio") <- ratio
  attr(out, "n_reliable") <- sum(ok)
  attr(out, "direction_change") <-
    any(ok) && any(ratio[ok] < 1) && any(ratio[ok] > 1)
  attr(out, "monotone") <-
    if (sum(ok) >= 2) all(diff(ratio[ok]) >= 0) || all(diff(ratio[ok]) <= 0)
    else NA
  class(out) <- c("dn_sweep_result", class(out))
  out
}

#' @export
print.dn_sweep_result <- function(x, ...) {
  ratio <- x$n_admissible_obs / x$expected
  cat("Threshold sweep of the aggregate null test\n")
  cat(sprintf("  thresholds: q%.2f - q%.2f (%d points)\n",
              min(x$threshold_q), max(x$threshold_q), nrow(x)))
  fin <- x$reliable & is.finite(ratio)
  if (any(fin)) {
    rr <- range(ratio[fin])
    cat(sprintf("  obs/expected ratio: %.3f - %.3f\n", rr[1], rr[2]))
  } else {
    cat("  no reliable thresholds (all expected counts < 5 or n < 30);\n",
        "   nothing can be concluded for this region at these thresholds.\n")
  }
  if (isTRUE(attr(x, "direction_change")))
    cat("  ** DIRECTION REVERSES across the sweep: a single-threshold\n",
        "    verdict for this region is an artefact of the threshold. **\n")
  else if (any(fin))
    cat("  direction is stable across the sweep.\n")
  if (any(!x$reliable))
    cat(sprintf("  %d/%d thresholds unreliable (expected < 5); not shown.\n",
                sum(!x$reliable), nrow(x)))
  invisible(x)
}

#' @export
plot.dn_sweep_result <- function(x, log = "y", ...) {
  ratio <- x$n_admissible_obs / x$expected
  ok <- is.finite(ratio) & x$reliable
  graphics::plot(x$threshold_q[ok], ratio[ok], type = "b", pch = 16,
                 log = log, xlab = "threshold (quantile of coverage)",
                 ylab = "observed / expected admissible species", ...)
  graphics::abline(h = 1, lty = 2)
  invisible(x)
}

#' Dark-diversity classification of species-region pairs
#'
#' Applies the full conjunction rule of the source studies: a pair is
#' *climatically admissible* if the calibrated establishment probability
#' exceeds `p_min` and every supplied null screen has SES >= `ses_min`; it is
#' *dark* if additionally the species is absent from the region.
#'
#' @param pairs data.frame with at least `species`, `region`, `present`
#'   (logical), a probability column `p_establish`, and one or more SES
#'   columns given in `ses_cols`.
#' @param p_min admissibility threshold on the establishment probability.
#'   The source studies' externally calibrated value (Youden threshold of a
#'   species-blocked, climate-only model validated on documented introduced
#'   ranges) is 0.177; recalibrate with [dn_calibrate()] when using a
#'   different model or metric.
#' @param ses_cols character vector of SES column names (e.g.
#'   `c("n1_ses", "n2_ses", "n3_ses")`).
#' @param ses_min SES criterion (default 2).
#' @return `pairs` with logical columns `admissible` and `dark`.
#' @examples
#' pairs <- data.frame(species = c("a", "a", "b"), region = c("r1", "r2", "r1"),
#'                     present = c(TRUE, FALSE, FALSE),
#'                     p_establish = c(0.6, 0.4, 0.05),
#'                     n2_ses = c(3.1, 2.5, 4.0))
#' dn_dark(pairs, ses_cols = "n2_ses")
#' @export
dn_dark <- function(pairs, p_min = 0.177,
                           ses_cols = character(), ses_min = 2) {
  adm <- pairs$p_establish >= p_min
  for (cn in ses_cols) adm <- adm & (pairs[[cn]] >= ses_min)
  pairs$admissible <- adm & !is.na(adm)
  pairs$dark <- pairs$admissible & !pairs$present
  pairs
}

#' Climatic dark-diversity ratio and its area adjustment
#'
#' `CDDR = n_dark / (n_dark + n_present)`. The raw ratio is confounded with
#' region area (a small region occupies little environmental space, so more
#' niches enclose it, while its smaller area supports fewer natives; in the
#' global analysis r = -0.44), so cross-region comparison must use the
#' residual from the CDDR ~ log(area) regression. **Never rank regions by raw
#' CDDR.**
#'
#' @param n_dark,n_present integer vectors per region.
#' @param area region areas (same order); required for `adjusted`.
#' @return data.frame with `dn_cddr`, and when `area` is given `cddr_residual`
#'   plus the fitted confound as attributes `fit`.
#' @examples
#' set.seed(1)
#' area <- 10^runif(50, 4, 6.7)
#' nd <- rpois(50, 40); np <- rpois(50, 60)
#' x <- dn_cddr(nd, np, area)
#' attr(x, "confound_r")     # never rank regions by raw cddr
#' head(x)
#' @export
dn_cddr <- function(n_dark, n_present, area = NULL) {
  ratio <- n_dark / (n_dark + n_present)
  out <- data.frame(cddr = ratio)
  if (!is.null(area)) {
    la <- log10(area)
    fit <- stats::lm(ratio ~ la)
    out$cddr_residual <- stats::resid(fit)
    attr(out, "fit") <- fit
    # a constant ratio (or a constant area) has no confound to report, and
    # cor() would only warn about a zero standard deviation
    attr(out, "confound_r") <-
      if (stats::sd(ratio) == 0 || stats::sd(la) == 0) NA_real_
      else stats::cor(ratio, la)
  }
  out
}

# ---------------------------------------------------------------------------
# Validation against documented introduced ranges
# ---------------------------------------------------------------------------

#' AUC with DeLong variance (fast algorithm of Sun & Xu 2014)
#'
#' @param score numeric vector; higher = more positive.
#' @param label logical or 0/1 vector.
#' @return list with `auc`, `var`, `ci` (95%).
#' @examples
#' set.seed(1)
#' lab <- rep(c(TRUE, FALSE), each = 200)
#' score <- ifelse(lab, rnorm(200, 1), rnorm(200, 0))
#' dn_auc(score, lab)
#' @export
dn_auc <- function(score, label) {
  y <- as.logical(label)
  if (length(score) != length(y))
    stop("`score` and `label` must have the same length")
  ok <- is.finite(score) & !is.na(y)
  score <- score[ok]; y <- y[ok]
  if (!any(y) || all(y))
    stop("AUC is undefined: `label` must contain both classes (got ",
         sum(y), " positives and ", sum(!y), " negatives)")
  x <- score[y]; z <- score[!y]
  m <- length(x); n <- length(z)
  r <- rank(c(x, z), ties.method = "average")
  rx <- rank(x, ties.method = "average")
  rz <- rank(z, ties.method = "average")
  V10 <- (r[seq_len(m)] - rx) / n
  V01 <- 1 - (r[m + seq_len(n)] - rz) / m
  auc <- mean(V10)
  v <- stats::var(V10) / m + stats::var(V01) / n
  list(auc = auc, var = v,
       ci = auc + c(-1, 1) * stats::qnorm(0.975) * sqrt(v))
}

#' Paired DeLong test between two correlated AUCs
#'
#' @param s1,s2 two score vectors on the same cases.
#' @param label logical or 0/1 vector.
#' @return list with `auc1`, `auc2`, `z`, `p` (two-sided).
#' @examples
#' set.seed(1)
#' lab <- rep(c(TRUE, FALSE), each = 300)
#' s1 <- ifelse(lab, rnorm(300, 1.0), rnorm(300, 0))
#' s2 <- s1 + rnorm(600, 0, 0.8)      # a noisier version of the same signal
#' dn_delong(s1, s2, lab)
#' @export
dn_delong <- function(s1, s2, label) {
  y <- as.logical(label)
  if (length(s1) != length(s2) || length(s1) != length(y))
    stop("`s1`, `s2` and `label` must have the same length")
  ok <- is.finite(s1) & is.finite(s2) & !is.na(y)
  s1 <- s1[ok]; s2 <- s2[ok]; y <- y[ok]
  if (!any(y) || all(y))
    stop("a paired AUC test is undefined: `label` must contain both classes")
  comp <- function(s) {
    x <- s[y]; z <- s[!y]; m <- length(x); n <- length(z)
    r <- rank(c(x, z), ties.method = "average")
    rx <- rank(x, ties.method = "average")
    rz <- rank(z, ties.method = "average")
    list(V10 = (r[seq_len(m)] - rx) / n,
         V01 = 1 - (r[m + seq_len(n)] - rz) / m)
  }
  a <- comp(s1); b <- comp(s2)
  auc1 <- mean(a$V10); auc2 <- mean(b$V10)
  m <- length(a$V10); n <- length(a$V01)
  S10 <- stats::cov(cbind(a$V10, b$V10)); S01 <- stats::cov(cbind(a$V01, b$V01))
  S <- S10 / m + S01 / n
  L <- c(1, -1)
  v <- drop(t(L) %*% S %*% L)
  z <- (auc1 - auc2) / sqrt(v)
  list(auc1 = auc1, auc2 = auc2, z = z,
       p = 2 * stats::pnorm(-abs(z)))
}

#' Realm- and area-matched negatives for introduced-range validation
#'
#' For each positive (region a species actually reached), samples `ratio`
#' negatives from unreached regions in the same biogeographic realm with the
#' nearest log-areas -- matching without which negatives are dominated by
#' climatically implausible regions and every metric rejects them equally.
#'
#' @param positives data.frame with `species`, `region`.
#' @param regions data.frame with `region`, `realm`, `area`.
#' @param reached list mapping species to reached regions.
#' @param ratio negatives per positive (default 5).
#' @param seed RNG seed.
#' @return data.frame `species`, `region`, `label`.
#' @examples
#' regs <- data.frame(region = paste0("r", 1:8), realm = rep(c("PA", "NT"), 4),
#'                    area = 10^runif(8, 4, 6))
#' pos <- data.frame(species = "sp1", region = "r1")
#' dn_match(pos, regs, reached = list(sp1 = "r1"), ratio = 3)
#' @export
dn_match <- function(positives, regions, reached, ratio = 5L,
                            seed = 1L) {
  set.seed(seed)
  la <- log10(regions$area); names(la) <- regions$region
  realm <- regions$realm; names(realm) <- regions$region
  out <- vector("list", nrow(positives))
  for (i in seq_len(nrow(positives))) {
    spx <- positives$species[i]; rg <- positives$region[i]
    pool <- regions$region[realm == realm[rg]]
    pool <- setdiff(pool, reached[[spx]])
    if (!length(pool)) next
    d <- abs(la[pool] - la[rg])
    neg <- pool[order(d)][seq_len(min(ratio, length(pool)))]
    out[[i]] <- data.frame(
      species = spx, region = c(rg, neg),
      label = c(TRUE, rep(FALSE, length(neg))))
  }
  do.call(rbind, out)
}

#' Calibrate the admissibility threshold (Youden)
#'
#' @param score establishment-model probabilities (or any coverage score).
#' @param label validation labels.
#' @return the Youden-optimal threshold.
#' @examples
#' set.seed(1)
#' lab <- rep(c(TRUE, FALSE), each = 200)
#' p <- ifelse(lab, rbeta(200, 4, 2), rbeta(200, 2, 4))
#' dn_calibrate(p, lab)     # Youden-optimal threshold
#' @export
dn_calibrate <- function(score, label) {
  y <- as.logical(label)
  if (!any(y) || all(y))
    stop("calibration needs both classes in `label`")
  th <- sort(unique(score))
  sens <- vapply(th, function(t) mean(score[y] >= t), 0)
  spec <- vapply(th, function(t) mean(score[!y] < t), 0)
  th[which.max(sens + spec - 1)]
}
