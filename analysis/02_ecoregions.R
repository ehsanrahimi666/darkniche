# ---------------------------------------------------------------------------
# darkniche -- Case study 2: the world's 619 terrestrial ecoregions
#
# Runs the package's reporting and diagnostic layer over the bundled results
# of the global analysis (35,087 species; 21,718,853 species-ecoregion pairs)
# and reproduces the published quantities from the package alone.
#
# Run:  Rscript analysis/02_ecoregions.R
# Out:  analysis/output/eco_*.csv, analysis/figures/fig4_*.png, fig5_*.png
# ---------------------------------------------------------------------------

library(darkniche)
dir.create("analysis/output",  showWarnings = FALSE, recursive = TRUE)
dir.create("analysis/figures", showWarnings = FALSE, recursive = TRUE)

eco   <- dn_data("ecoregion")
sweep <- dn_data("sweep")
src   <- dn_data("sources")
rt    <- dn_data("region_taxon")
auc   <- dn_data("auc")
cont  <- dn_data("containment")

# --- 1. the report object ---------------------------------------------------
rep <- dn_report(eco, sweep)
print(rep)
dn_table_export(rep$table, "analysis/output/eco_report.csv", digits = 4)

# --- 2. what the area adjustment changes ------------------------------------
tab <- rep$table
top_raw <- head(tab[order(-tab$cddr), c("region", "cddr", "area_km2")], 8)
top_adj <- head(tab[order(-tab$residual), c("region", "cddr", "residual",
                                            "area_km2")], 8)
cat("\nTop 8 by RAW cddr (area-confounded, must not be used):\n")
print(top_raw, row.names = FALSE)
cat("\nTop 8 by AREA-ADJUSTED residual (the comparable quantity):\n")
print(top_adj, row.names = FALSE)
ext <- dn_table_extremes(tab$region, rep$cddr, n_obs = tab$n_obs,
                         n_dark = tab$n_dark, n = 8)
dn_table_export(ext, "analysis/output/eco_extremes.csv", digits = 4)
cat(sprintf("\nSpearman(raw rank, adjusted rank) = %.3f; median |rank shift| = %.0f\n",
            stats::cor(tab$rank_raw, tab$rank_adjusted, method = "spearman"),
            stats::median(abs(tab$rank_shift))))

# --- 3. the direction census ------------------------------------------------
cen_all <- dn_direction_census(sweep, reliable_only = FALSE)
cen_rel <- dn_direction_census(sweep, reliable_only = TRUE)
census <- data.frame(
  quantile = cen_all$census$quantile,
  enriched_all = cen_all$census$n_enriched,
  depleted_all = cen_all$census$n_depleted,
  enriched_reliable = cen_rel$census$n_enriched,
  depleted_reliable = cen_rel$census$n_depleted)
print(census, row.names = FALSE)
dn_table_export(census, "analysis/output/eco_direction_census.csv", digits = 2)
cat(sprintf("direction reverses in %d/%d regions (all) and %d/%d (reliable only)\n",
            cen_all$n_changed, cen_all$n_regions,
            cen_rel$n_changed, cen_rel$n_regions))

# --- 4. source composition and taxonomic structure --------------------------
share <- colMeans(src[, c("near", "regional", "distant")])
cat(sprintf("\nmean share of dark diversity by origin: near %.3f, regional %.3f, distant %.3f\n",
            share["near"], share["regional"], share["distant"]))
tx <- stats::aggregate(cddr ~ taxon, data = rt, FUN = stats::median)
tx <- tx[order(-tx$cddr), ]
print(tx, row.names = FALSE)
dn_table_export(tx, "analysis/output/eco_taxon_cddr.csv", digits = 3)

# --- 5. figure 4: global structure ------------------------------------------
png("analysis/figures/fig4_global.png", 2400, 1750, res = 210)
op <- par(mfrow = c(2, 2), mar = c(4.2, 4.4, 2.8, 1))

dn_plot_cddr_area(rep$cddr, tab$area_km2)
title(main = "(a) the area confound, shown not hidden")

plot(tab$rank_raw, tab$rank_adjusted, pch = 16, cex = 0.5, col = "#4F81BD66",
     xlab = "rank on raw CDDR", ylab = "rank on area-adjusted residual",
     main = "(b) what the adjustment changes")
abline(0, 1, lty = 2, col = "grey40")

ord <- order(src$distant)
sel <- src[ord[round(seq(1, nrow(src), length.out = 14))], ]
rownames(sel) <- substr(sel$region, 1, 28)
dn_plot_sources(sel, main = "(c) dark diversity is a distant pool")

M <- stats::xtabs(cddr ~ region + taxon, data = rt)
sel_r <- names(sort(tapply(rt$cddr, rt$region, mean)))
sel_r <- sel_r[round(seq(1, length(sel_r), length.out = 14))]
dn_plot_matrix(as.matrix(M[sel_r, ]), main = "(d) dispersal ability orders taxa",
               values = FALSE)
par(op); dev.off()

# --- 6. figure 5: the diagnostic -------------------------------------------
png("analysis/figures/fig5_sweep.png", 2400, 900, res = 215)
op <- par(mfrow = c(1, 3), mar = c(4.4, 4.4, 2.8, 1))

qs <- sort(unique(sweep$quantile))
rel <- sweep[sweep$reliable %in% c(TRUE, "True"), ]
wide <- stats::reshape(rel[, c("region", "quantile", "enrichment")],
                       idvar = "region", timevar = "quantile",
                       direction = "wide")
E <- as.matrix(wide[, -1]); rownames(E) <- wide$region
full <- stats::complete.cases(E)
E <- E[full, , drop = FALSE]
plot(range(qs), range(E[E > 0], na.rm = TRUE), type = "n", log = "y",
     xlab = "threshold (quantile of coverage)",
     ylab = "observed / expected admissible species",
     main = "(a) every ecoregion's sweep")
flip <- apply(E, 1, function(e) any(e < 1) && any(e > 1))
for (i in seq_len(nrow(E)))
  lines(qs, E[i, ], col = if (flip[i]) "#B2182B22" else "#4F81BD44")
lines(qs, apply(E, 2, stats::median), lwd = 3)
abline(h = 1, lty = 2)
legend("topleft", bty = "n", cex = 0.85,
       legend = sprintf("%d of %d reliable ecoregions (%.0f%%) cross the line",
                        sum(flip), nrow(E), 100 * mean(flip)))

bp <- t(as.matrix(census[, c("enriched_reliable", "depleted_reliable")]))
barplot(bp, beside = FALSE, names.arg = paste0("q", census$quantile),
        col = c("#2166AC", "#B2182B"), border = NA,
        ylab = "ecoregions", main = "(b) the verdict is the threshold")
legend("topright", bty = "n", fill = c("#2166AC", "#B2182B"), border = NA,
       legend = c("enriched", "depleted"), cex = 0.85)

with(auc[auc$metric %in% c("DNC_q", "DNC_auc", "CFI_w", "CFI", "I", "D",
                           "NSI_q"), ], {
  o <- order(auc)
  plot(auc[o], seq_along(o), pch = 16, xlim = c(0.66, 0.73), yaxt = "n",
       xlab = "AUC (documented introductions)", ylab = "",
       main = "(c) validation on real introductions")
  segments(lo[o], seq_along(o), hi[o], seq_along(o))
  axis(2, seq_along(o), metric[o], las = 1, cex.axis = 0.8)
  abline(v = auc[metric == "D"], lty = 3, col = "#B2182B")
})
par(op); dev.off()

saveRDS(list(report = rep, census = census, taxon = tx, share = share),
        "analysis/output/eco_results.rds")
cat("\ncase study 2 complete\n")
