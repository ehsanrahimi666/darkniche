# ---------------------------------------------------------------------------
# darkniche -- Case study 1: a controlled simulation with known ground truth
#
# Symmetric overlap statistics cannot express containment. Here we build a
# landscape in which admissibility is known by construction, so the claim can
# be tested rather than argued.
#
# The ground truth is defined entirely from the *generating parameters* --
# never from an estimated density and never from any metric the analysis then
# evaluates. Species i is generated with tolerance centre c_i and scale w_i;
# its physiological tolerance is the disc of radius 2 w_i around c_i. Region
# j is generated with climate centre r_j and scale v_j. The pair is
# admissible when at least 90% of the region's climate mass falls inside the
# species' tolerance disc, estimated by Monte Carlo on fresh draws that are
# not the points used for density estimation.
#
# Run:  Rscript analysis/01_simulation.R
# Out:  analysis/output/sim_*.csv, analysis/figures/fig3_simulation.png
# ---------------------------------------------------------------------------

library(darkniche)
set.seed(20260807)

dir.create("analysis/output",  showWarnings = FALSE, recursive = TRUE)
dir.create("analysis/figures", showWarnings = FALSE, recursive = TRUE)

# --- 1. environmental space -------------------------------------------------
n_bg <- 30000
comp <- sample(1:3, n_bg, TRUE, prob = c(0.45, 0.35, 0.20))
mu_bg <- rbind(c(-1.4, 0.4), c(1.1, -0.9), c(0.2, 1.7))
sd_bg <- rbind(c(0.9, 0.7), c(1.1, 0.9), c(0.6, 0.7))
bg <- mu_bg[comp, ] + matrix(rnorm(2 * n_bg), n_bg, 2) * sd_bg[comp, ]
sp <- dn_space(bg, R = 80)
Zg <- dn_density(bg, sp)

# --- 2. species with controlled niche breadth -------------------------------
n_sp <- 120
sp_ctr <- cbind(runif(n_sp, -2.5, 2.5), runif(n_sp, -2.0, 2.5))
sp_wid <- exp(runif(n_sp, log(0.20), log(1.60)))       # narrow .. very broad
Zs <- matrix(0, n_sp, sp$R^2)
sp_mu <- matrix(0, n_sp, 2); sp_cov <- vector("list", n_sp)
for (i in seq_len(n_sp)) {
  pts <- cbind(rnorm(600, sp_ctr[i, 1], sp_wid[i]),
               rnorm(600, sp_ctr[i, 2], sp_wid[i]))
  Zi <- dn_density(pts, sp)
  Zs[i, ] <- as.numeric(Zi)
  mm <- dn_moments(Zi, sp); sp_mu[i, ] <- mm$mu; sp_cov[[i]] <- mm$cov
}

# --- 3. regions with controlled size ----------------------------------------
n_rg <- 60
rg_ctr <- cbind(runif(n_rg, -2.2, 2.2), runif(n_rg, -1.8, 2.2))
rg_wid <- exp(runif(n_rg, log(0.12), log(0.90)))       # small .. large
Ze <- matrix(0, n_rg, sp$R^2); eco_mu <- matrix(0, n_rg, 2)
for (j in seq_len(n_rg)) {
  pts <- cbind(rnorm(600, rg_ctr[j, 1], rg_wid[j]),
               rnorm(600, rg_ctr[j, 2], rg_wid[j]))
  Zj <- dn_density(pts, sp)
  Ze[j, ] <- as.numeric(Zj)
  eco_mu[j, ] <- dn_moments(Zj, sp)$mu
}

# --- 4. ground truth from the generating parameters only --------------------
n_mc <- 4000
truth <- matrix(FALSE, n_sp, n_rg)
for (j in seq_len(n_rg)) {
  draw <- cbind(rnorm(n_mc, rg_ctr[j, 1], rg_wid[j]),
                rnorm(n_mc, rg_ctr[j, 2], rg_wid[j]))
  for (i in seq_len(n_sp)) {
    d <- sqrt((draw[, 1] - sp_ctr[i, 1])^2 + (draw[, 2] - sp_ctr[i, 2])^2)
    truth[i, j] <- mean(d <= 2 * sp_wid[i]) >= 0.90
  }
}
lab <- as.vector(truth)
cat(sprintf("ground truth: %d of %d pairs admissible (%.1f%%)\n",
            sum(lab), length(lab), 100 * mean(lab)))

# --- 5. every metric scored against that truth ------------------------------
m <- dn_metrics(Zs, Ze, q = 0.95, sp_mu = sp_mu, sp_cov = sp_cov,
                eco_mu = eco_mu, centres = dn_centres(sp))
cands <- c("DNC_q", "DNC_auc", "NSI_q", "CFI", "CFI_w", "D", "I")
auc_tab <- do.call(rbind, lapply(cands, function(k) {
  v <- as.vector(m[[k]]); a <- dn_auc(v, lab)
  dl <- dn_delong(v, as.vector(m$D), lab)
  data.frame(metric = k, auc = a$auc, lo = a$ci[1], hi = a$ci[2],
             z_vs_D = if (k == "D") NA else dl$z,
             p_vs_D = if (k == "D") NA else dl$p)
}))
auc_tab <- auc_tab[order(-auc_tab$auc), ]
print(auc_tab, row.names = FALSE)
dn_table_export(auc_tab, "analysis/output/sim_metric_auc.csv")

# --- 6. where does the symmetric statistic fail? ----------------------------
ratio <- outer(sp_wid, rg_wid, "/")
strat <- cut(as.vector(ratio),
             stats::quantile(as.vector(ratio), c(0, 1/3, 2/3, 1)),
             labels = c("niche ~ region", "intermediate", "niche >> region"),
             include.lowest = TRUE)
safe_auc <- function(v, y) if (length(unique(y)) < 2) NA_real_ else dn_auc(v, y)$auc
safe_p <- function(a, b, y) if (length(unique(y)) < 2) NA_real_ else
  dn_delong(a, b, y)$p
strat_tab <- do.call(rbind, lapply(levels(strat), function(s) {
  k <- strat == s
  data.frame(stratum = s, n = sum(k), n_admissible = sum(lab[k]),
             DNC_q   = safe_auc(as.vector(m$DNC_q)[k],   lab[k]),
             DNC_auc = safe_auc(as.vector(m$DNC_auc)[k], lab[k]),
             D       = safe_auc(as.vector(m$D)[k],       lab[k]),
             I       = safe_auc(as.vector(m$I)[k],       lab[k]),
             p_DNCauc_vs_D = safe_p(as.vector(m$DNC_auc)[k],
                                    as.vector(m$D)[k], lab[k]))
}))
print(strat_tab, row.names = FALSE)
dn_table_export(strat_tab, "analysis/output/sim_by_stratum.csv")

# --- 7. cost of getting it wrong: the recovered admissible pool -------------
recov <- do.call(rbind, lapply(c("DNC_q", "DNC_auc", "D", "I"), function(k) {
  v <- as.vector(m[[k]]); th <- dn_calibrate(v, lab)
  data.frame(metric = k, threshold = th,
             recall = mean(v[lab] >= th),
             precision = mean(lab[v >= th]),
             pool_size = sum(v >= th), true_pool = sum(lab),
             pool_error = sum(v >= th) / sum(lab) - 1)
}))
print(recov, row.names = FALSE)
dn_table_export(recov, "analysis/output/sim_recovery.csv")

# --- 8. figure --------------------------------------------------------------
png("analysis/figures/fig3_simulation.png", 2400, 900, res = 220)
op <- par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2.6, 1))

plot(as.vector(m$D), as.vector(m$DNC_auc), pch = 16, cex = 0.35,
     col = ifelse(lab, "#B2182B66", "#4F81BD33"),
     xlab = "Schoener's D", ylab = expression(DNC[auc]),
     main = "(a) admissible pairs that D scores low")
legend("bottomright", bty = "n", pch = 16, col = c("#B2182B", "#4F81BD"),
       legend = c("truly admissible", "not admissible"), cex = 0.8)

dn_plot_roc(list("DNC_q" = as.vector(m$DNC_q),
                 "DNC_auc" = as.vector(m$DNC_auc),
                 "D" = as.vector(m$D)), lab,
            main = "(b) recovery of known admissibility")

keep <- !is.na(strat_tab$D)
xs <- which(keep)
yr <- range(c(strat_tab$DNC_q, strat_tab$DNC_auc, strat_tab$D, strat_tab$I),
            na.rm = TRUE)
plot(range(seq_len(nrow(strat_tab))) + c(-0.35, 0.35),
     c(yr[1] - 0.02, 1.005), type = "n",
     xaxt = "n", xlab = "", ylab = "AUC against known admissibility",
     main = "(c) the gap is a containment effect")
axis(1, seq_len(nrow(strat_tab)), strat_tab$stratum, cex.axis = 0.72)
cols <- c(DNC_q = "#2166AC", DNC_auc = "#EF8A62", D = "#B2182B")
for (k in names(cols))
  lines(xs, strat_tab[[k]][keep], type = "b", pch = 16, col = cols[k], lwd = 2)
legend("bottomleft", bty = "n", lwd = 2, pch = 16, col = cols,
       legend = c(expression(DNC[q]), expression(DNC[auc]), "D"), cex = 0.8)
par(op); dev.off()

saveRDS(list(auc = auc_tab, strat = strat_tab, recov = recov),
        "analysis/output/sim_results.rds")
cat("case study 1 complete\n")
