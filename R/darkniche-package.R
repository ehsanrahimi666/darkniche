#' darkniche: Climatic Dark Diversity and Directional Niche Coverage
#'
#' Which species are absent from a region although its climate is inside
#' their tolerance? `darkniche` answers that question end to end: it estimates
#' niches and regional climate availabilities on one shared environmental
#' lattice, scores climatic admissibility with a *directional* metric family
#' that expresses containment (a question symmetric overlap statistics such
#' as Schoener's *D* provably cannot express), screens each pair against
#' three analytical null models, calibrates the admissibility threshold
#' against documented introduced ranges, and reports the resulting dark
#' diversity with the two guardrails the underlying studies found to be
#' indispensable: an area adjustment on the dark-diversity ratio, and a
#' threshold-sweep diagnostic that reveals when the direction of an
#' aggregate null-model test is an artefact of the threshold rather than a
#' property of the region.
#'
#' @section Workflow:
#' \enumerate{
#'   \item [dn_space()], [dn_density()] -- one lattice, one density per
#'     species and per region.
#'   \item [dn_metrics()] -- the directional family (`DNC_q`, `DNC_auc`,
#'     `NSI_q`, `CFI`, `CFI_w`, `NOV`) alongside `D`, `I` and `BC`.
#'   \item [dn_null_breadth()], [dn_null_shift()], [dn_null_relocate()] --
#'     null screens; [dn_ses()] standardises them.
#'   \item [dn_match()], [dn_auc()], [dn_delong()], [dn_calibrate()] --
#'     validation against documented introduced ranges.
#'   \item [dn_dark()], [dn_cddr()], [dn_report()] -- accounting.
#'   \item [dn_aggregate()], [dn_sweep()], [dn_direction_census()] --
#'     region-level inference and its diagnostic.
#' }
#'
#' @section Guardrails:
#' [dn_aggregate()] suppresses the effect size, p-value and enrichment ratio
#' whenever the expected count falls below five; [dn_cddr()] returns the
#' area-adjusted residual beside the raw ratio and documents that raw CDDR
#' must not rank regions; [dn_sweep()] prints a reversal warning from the
#' object itself, computed only from thresholds the test considers reliable.
#'
#' @section Example data:
#' [dn_data()] gives the published results for all 619 WWF terrestrial
#' ecoregions, so every example and test runs on real output.
#'
#' @references
#' Rahimi, E. & Jung, C. (in review) Directional niche coverage and the
#' climatic dark diversity of Iranian ecoregions.
#'
#' Rahimi, E. & Jung, C. (in review) Climatic dark diversity of the world's
#' terrestrial ecoregions and the threshold dependence of aggregate
#' null-model tests.
#'
#' @keywords internal
"_PACKAGE"
