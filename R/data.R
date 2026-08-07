# Bundled example data and high-level workflow wrappers.

#' Load bundled ecoregion example data
#'
#' Results of the global analysis of Rahimi & Jung (in review) for all 619
#' WWF terrestrial ecoregions (35,087 species; 21,718,853 species-ecoregion
#' pairs), shipped so that every example, test and vignette in `darkniche`
#' runs on real output rather than on toy numbers.
#'
#' Available datasets:
#'
#' \describe{
#'   \item{`ecoregion`}{619 rows: `region`, `n_obs`, `n_dark`, `cddr`,
#'     `area_km2`.}
#'   \item{`sweep`}{2,476 rows: `region`, `quantile`, `threshold`,
#'     `enrichment`, `reliable` -- the aggregate null test at four
#'     thresholds per ecoregion.}
#'   \item{`sources`}{619 rows: `near`, `regional`, `distant` shares of each
#'     ecoregion's dark diversity by native-range distance.}
#'   \item{`taxon`}{dark-diversity rate per taxonomic group with 95% CI.}
#'   \item{`region_taxon`}{3,095 rows: CDDR for each ecoregion x taxon.}
#'   \item{`auc`}{AUC of each metric against documented introduced ranges.}
#'   \item{`containment`}{AUC by niche-to-region size stratum.}
#' }
#'
#' @param which name of the dataset (see above).
#' @return a data.frame.
#' @examples
#' eco <- dn_data("ecoregion")
#' str(eco)
#' summary(eco$cddr)
#'
#' # metric performance against documented introductions
#' dn_data("auc")
#' @export
dn_data <- function(which = c("ecoregion", "sweep", "sources", "taxon",
                              "region_taxon", "auc", "containment")) {
  which <- match.arg(which)
  f <- c(ecoregion = "ecoregion_darkniche.csv", sweep = "threshold_sweep.csv",
         sources = "source_shares.csv", taxon = "taxon_summary.csv",
         region_taxon = "ecoregion_taxon.csv", auc = "metric_auc.csv",
         containment = "auc_by_containment.csv")[[which]]
  p <- system.file("extdata", "ecoregions", f, package = "darkniche")
  if (!nzchar(p)) stop("bundled dataset not found: ", which)
  utils::read.csv(p, stringsAsFactors = FALSE)
}

#' Direction census of an aggregate null test across thresholds
#'
#' Summarises, over many regions, how many are called enriched and how many
#' depleted at each threshold, and what share of regions change direction
#' somewhere in the sweep. This is the quantity a study must report before
#' any single-threshold enrichment claim can be trusted.
#'
#' @param sweep a long data.frame with columns `region`, `quantile` and
#'   `enrichment` (observed/expected), and optionally `reliable`.
#' @param reliable_only if `TRUE` (default) only rows flagged reliable are
#'   used; unreliable cells (expected count < 5) are not interpretable.
#' @return list with `census` (per threshold: n enriched, n depleted),
#'   `n_regions`, `n_changed`, `frac_changed`, `n_monotone`.
#' @examples
#' cen <- dn_direction_census(dn_data("sweep"))
#' cen$census
#' round(cen$frac_changed, 3)     # share of ecoregions that flip direction
#' @export
dn_direction_census <- function(sweep, reliable_only = TRUE) {
  s <- sweep
  if (reliable_only && "reliable" %in% names(s))
    s <- s[s$reliable %in% c(TRUE, "True", "TRUE"), , drop = FALSE]
  qs <- sort(unique(s$quantile))
  census <- data.frame(
    quantile = qs,
    n_enriched = vapply(qs, function(q)
      sum(s$enrichment[s$quantile == q] > 1, na.rm = TRUE), 0L),
    n_depleted = vapply(qs, function(q)
      sum(s$enrichment[s$quantile == q] < 1, na.rm = TRUE), 0L))
  sp <- split(s, s$region)
  full <- vapply(sp, function(d) length(unique(d$quantile)) == length(qs),
                 TRUE)
  sp <- sp[full]
  chg <- vapply(sp, function(d) {
    e <- d$enrichment[order(d$quantile)]
    any(e < 1, na.rm = TRUE) && any(e > 1, na.rm = TRUE)
  }, TRUE)
  mono <- vapply(sp, function(d) {
    e <- d$enrichment[order(d$quantile)]
    e <- e[is.finite(e)]
    length(e) < 2 || all(diff(e) >= 0) || all(diff(e) <= 0)
  }, TRUE)
  list(census = census, n_regions = length(sp), n_changed = sum(chg),
       frac_changed = mean(chg), n_monotone = sum(mono))
}

#' One-call region report
#'
#' Assembles the quantities a region-level dark-diversity study must report
#' together: raw ratio, area-adjusted residual, rank on each, and the
#' direction census of the region's own threshold sweep. Printing the result
#' surfaces the two failure modes -- an area-driven ranking and a
#' threshold-driven direction -- instead of leaving them to a footnote.
#'
#' @param eco data.frame with `region`, `n_obs`, `n_dark`, `area_km2`.
#' @param sweep optional long sweep table (as [dn_data()] `"sweep"`).
#' @return object of class `dn_report`.
#' @examples
#' rep <- dn_report(dn_data("ecoregion"), dn_data("sweep"))
#' rep
#' head(rep$table)
#' @export
dn_report <- function(eco, sweep = NULL) {
  x <- dn_cddr(eco$n_dark, eco$n_obs, eco$area_km2)
  tab <- data.frame(region = eco$region, n_obs = eco$n_obs,
                    n_dark = eco$n_dark, area_km2 = eco$area_km2,
                    cddr = x$cddr, residual = x$cddr_residual,
                    rank_raw = rank(-x$cddr),
                    rank_adjusted = rank(-x$cddr_residual))
  tab$rank_shift <- tab$rank_adjusted - tab$rank_raw
  out <- list(table = tab, cddr = x, confound_r = attr(x, "confound_r"),
              fit = attr(x, "fit"),
              census = if (!is.null(sweep)) dn_direction_census(sweep) else NULL)
  class(out) <- "dn_report"
  out
}

#' @export
print.dn_report <- function(x, ...) {
  cat("Climatic dark-diversity report\n")
  cat(sprintf("  regions: %d\n", nrow(x$table)))
  cat(sprintf("  CDDR: %.3f - %.3f (median %.3f)\n",
              min(x$table$cddr), max(x$table$cddr),
              stats::median(x$table$cddr)))
  if (is.na(x$confound_r))
    cat("  area confound: not estimable (constant ratio or area)\n")
  else
    cat(sprintf("  area confound: r = %.3f -- rank on the residual, never on raw CDDR\n",
                x$confound_r))
  big <- sum(abs(x$table$rank_shift) > nrow(x$table) / 10)
  cat(sprintf("  %d regions move more than 10%% of the ranking when area is removed\n",
              big))
  if (!is.null(x$census)) {
    cn <- x$census
    cat(sprintf("  threshold sweep: %d of %d regions (%.0f%%) change direction\n",
                cn$n_changed, cn$n_regions, 100 * cn$frac_changed))
    cat("  direction census by threshold:\n")
    print(cn$census, row.names = FALSE)
    if (cn$frac_changed > 0.5)
      cat("  ** a single-threshold enrichment claim is not supportable here **\n")
  }
  invisible(x)
}
