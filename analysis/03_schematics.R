# ---------------------------------------------------------------------------
# darkniche -- Figures 1 and 2 of the software paper: architecture and engine
#
# Run:  Rscript analysis/03_schematics.R
# Out:  analysis/figures/fig1_architecture.png, fig2_engine.png
# ---------------------------------------------------------------------------

library(darkniche)
dir.create("analysis/figures", showWarnings = FALSE, recursive = TRUE)

box <- function(x, y, w, h, lab, fill, txt = "black", cex = 0.72, lwd = 1) {
  rect(x - w/2, y - h/2, x + w/2, y + h/2, col = fill, border = "grey35",
       lwd = lwd)
  text(x, y, lab, cex = cex, col = txt)
}
arrow <- function(x0, y0, x1, y1, col = "grey30", lwd = 1.6) {
  arrows(x0, y0, x1, y1, length = 0.06, col = col, lwd = lwd)
}

# ---------------------------------------------------------------------------
# Figure 1: modular architecture and data flow
# ---------------------------------------------------------------------------
png("analysis/figures/fig1_architecture.png", 2400, 1500, res = 220)
par(mar = c(0.4, 0.4, 1.6, 0.4))
plot(NA, xlim = c(0, 100), ylim = c(1, 61), axes = FALSE, xlab = "", ylab = "",
     main = "")

g_in  <- "#DCE9F5"; g_eng <- "#F4E3D3"; g_nul <- "#E8E2F0"
g_out <- "#DDEEDD"; g_val <- "#FBE3E4"; g_rep <- "#EDEDED"

# inputs
text(14, 59, "Inputs", font = 2, cex = 0.8)
box(14, 53, 30, 6, "environmental scores (climate PCs)\nfor the global background", g_in)
box(14, 45, 30, 6, "occupied cells per species\noccupied cells per region", g_in)

# engine
text(14, 38.2, "Shared engine", font = 2, cex = 0.8)
box(14, 32, 30, 6.5, "dn_space()  \u2014 one lattice for every object\ndn_density() \u2014 fixed-bandwidth KDE", g_eng)
box(14, 23, 30, 6.5, "dn_envelope() / dn_entryq()\nq-envelopes and cumulative rank", g_eng)

# metrics
text(50, 59, "Directional metrics", font = 2, cex = 0.8)
box(50, 51, 30, 9,
    "dn_metrics()\nDNC_q  DNC_auc  NSI_q  CFI  CFI_w  NOV\nSchoener's D, Warren's I, BC for comparison",
    g_eng)

# nulls
text(50, 42.5, "Null screens", font = 2, cex = 0.8)
box(50, 35, 30, 8,
    "dn_null_breadth()  N1 breadth-matched\ndn_null_shift()    N2 environmental shift\ndn_null_relocate() N3 geographic relocation",
    g_nul)

# validation
text(86, 59, "Calibration", font = 2, cex = 0.8)
box(86, 51, 26, 9,
    "dn_match()  matched negatives\ndn_auc() / dn_delong()\ndn_calibrate()  Youden threshold",
    g_val)

# accounting
text(50, 26, "Dark-diversity accounting", font = 2, cex = 0.8)
box(50, 19, 30, 7, "dn_dark()  conjunction rule\ndn_cddr()  area-adjusted ratio", g_out)

# diagnostics
text(86, 34, "Diagnostics", font = 2, cex = 0.8)
box(86, 26, 26, 9,
    "dn_aggregate()  region-level test\n(reliability flags)\ndn_sweep() / dn_direction_census()",
    g_val)

# reporting
text(50, 10.5, "Reporting", font = 2, cex = 0.8)
box(50, 4.5, 62, 6.5,
    "dn_report()  |  dn_plot_containment, dn_plot_density, dn_plot_cddr_area, dn_plot_matrix, dn_plot_sources, dn_plot_roc, dn_plot_sweeps  |  dn_table_extremes, dn_table_export",
    g_rep, cex = 0.6)

arrow(14, 50, 14, 48)
arrow(6, 42, 6, 35.5)
arrow(6, 28.7, 6, 26.5)
arrow(29, 23, 35, 45)
arrow(29, 32, 35, 48)
arrow(40, 46.5, 40, 39.2)
arrow(65, 51, 73, 51)
arrow(86, 46.4, 86, 31)
arrow(74, 47, 63, 21.5)        # calibrated threshold into accounting
arrow(40, 31, 40, 22.6)
arrow(65, 19, 73, 24)
arrow(40, 15.4, 40, 8)
arrow(86, 21.4, 64, 6.5)
text(71.5, 32, "calibrated\nthreshold", cex = 0.56, col = "grey25")
dev.off()

# ---------------------------------------------------------------------------
# Figure 2: the shared engine, and why direction matters
# ---------------------------------------------------------------------------
set.seed(7)
png("analysis/figures/fig2_engine.png", 2400, 1650, res = 215)
layout(matrix(c(1, 2, 3, 4, 4, 4), 2, 3, byrow = TRUE), heights = c(1.05, 1))
par(mar = c(4, 4, 2.8, 1))

bg <- rbind(cbind(rnorm(9000, -1.4, 0.9), rnorm(9000, 0.4, 0.7)),
            cbind(rnorm(7000, 1.1, 1.1), rnorm(7000, -0.9, 0.9)),
            cbind(rnorm(4000, 0.2, 1.7), rnorm(4000, 1.7, 0.7)))
sp <- dn_space(bg, R = 90)
Zg <- dn_density(bg, sp)
Zs <- dn_density(cbind(rnorm(1200, -0.3, 1.25), rnorm(1200, 0.2, 1.05)), sp)
Ze <- dn_density(cbind(rnorm(500, 0.55, 0.28), rnorm(500, -0.15, 0.24)), sp)

dn_plot_density(Zg, sp, q = NA, main = "(a) one shared environmental space")
dn_plot_density(Zs, sp, q = 0.95,
                main = "(b) a niche and its 95% envelope")
dn_plot_density(Zs, sp, q = 0.95, overlay = Ze,
                main = "(c) a region inside that envelope")

m <- dn_metrics(rbind(as.numeric(Zs)), rbind(as.numeric(Ze)), q = 0.95)
legend("bottomleft", bty = "n", cex = 0.72, text.col = "grey15",
       legend = c(sprintf("D = %.2f", m$D),
                  sprintf("DNC_q = %.2f", m$DNC_q),
                  sprintf("DNC_auc = %.2f", m$DNC_auc)))

par(mar = c(4.4, 4.6, 3.4, 1))
ctr <- dn_centres(sp)
row0 <- which.min(abs(ctr$y - 0.2))                 # transect through the niche
Wq   <- matrix(dn_weights(Zs, q = 0.95, statistic = "q"),   sp$R, sp$R)
Wauc <- matrix(dn_weights(Zs, statistic = "auc"),           sp$R, sp$R)
zs <- Zs[, row0] / max(Zs[, row0])
ze <- Ze[, row0] / max(Ze[, row0])

plot(ctr$x, zs, type = "n", ylim = c(0, 1.18), xlab = "PC1",
     ylab = "relative density  /  coverage weight",
     main = "(d) what the two coverage statistics integrate, along a transect through the niche")
polygon(c(ctr$x[1], ctr$x, ctr$x[sp$R]), c(0, zs, 0), col = "#E6E6E6",
        border = NA)
polygon(c(ctr$x[1], ctr$x, ctr$x[sp$R]), c(0, ze, 0), col = "#B2182B33",
        border = "#B2182B")
lines(ctr$x, Wq[, row0], col = "#2166AC", lwd = 2.4, type = "s")
lines(ctr$x, Wauc[, row0] / max(Wauc), col = "#EF8A62", lwd = 2.4)
legend("topright", bty = "n", cex = 0.78, lwd = c(NA, NA, 2.4, 2.4),
       pch = c(15, 15, NA, NA), col = c("#E6E6E6", "#B2182B33", "#2166AC",
                                        "#EF8A62"),
       legend = c("species niche density", "region climate density",
                  "DNC_q weight (0/1 at the 95% envelope)",
                  "DNC_auc weight (graded, threshold-free)"))
mtext("DNC is the region-density-weighted integral of these weights: a step function fixes one envelope, the graded form integrates over all of them",
      side = 3, line = 0.15, cex = 0.6, col = "grey25")
dev.off()

cat("schematics complete\n")
