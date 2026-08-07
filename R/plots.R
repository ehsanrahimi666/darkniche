# Publication-quality figures and tables for directional-niche studies.
# Every plot works from package objects alone (base graphics, no extra
# dependencies) and follows the figure designs of the source papers.

.dn_pal <- function(n) grDevices::hcl.colors(n, "Zissou 1")
.dn_heat <- function(n) grDevices::hcl.colors(n, "Inferno", rev = TRUE)
.dn_div <- function(n) grDevices::hcl.colors(n, "Blue-Red 3")

#' Plot a density in environmental space
#'
#' Filled image of a lattice density with optional envelope contour and a
#' second density overlaid as contour lines -- the building block of
#' species-vs-region panels.
#'
#' @param Z `R x R` density (species or region).
#' @param sp a `dn_space`.
#' @param q optional envelope level to outline (default 0.95; `NA` for none).
#' @param overlay optional second density drawn as contours (e.g. the region
#'   over a species niche).
#' @param main title.
#' @param ... passed to [graphics::image()].
#' @return invisibly, `Z`.
#' @examples
#' set.seed(1)
#' sp <- dn_space(cbind(rnorm(4000), rnorm(4000)), R = 60)
#' Zs <- dn_density(cbind(rnorm(800, 0, 1), rnorm(800, 0, 1)), sp)
#' Ze <- dn_density(cbind(rnorm(300, 0.4, 0.25), rnorm(300, 0, 0.25)), sp)
#' dn_plot_density(Zs, sp, overlay = Ze, main = "region inside niche")
#' @export
dn_plot_density <- function(Z, sp, q = 0.95, overlay = NULL,
                            main = "", ...) {
  ctr <- dn_centres(sp)
  graphics::image(ctr$x, ctr$y, Z, col = .dn_heat(64), useRaster = TRUE,
                  xlab = "PC1", ylab = "PC2", main = main, ...)
  if (!is.na(q)) {
    E <- matrix(as.numeric(dn_envelope(Z, q)), sp$R, sp$R)
    graphics::contour(ctr$x, ctr$y, E, levels = 0.5, add = TRUE,
                      drawlabels = FALSE, lwd = 1.4, col = "grey25")
  }
  if (!is.null(overlay))
    graphics::contour(ctr$x, ctr$y, overlay, nlevels = 5, add = TRUE,
                      drawlabels = FALSE, col = "#B2182B", lwd = 1.4)
  invisible(Z)
}

#' Conceptual containment figure (papers' Figure 1)
#'
#' Generates the three-panel demonstration that symmetric overlap cannot
#' answer a containment question: a region nested in a niche core, nested at
#' the margin, and outside the niche, with D, DNC_q and DNC_auc annotated.
#'
#' @param R lattice resolution (default 100).
#' @param broad,narrow Gaussian widths (lattice cells) of niche and region.
#' @param offsets horizontal offsets of the region for the three panels.
#' @return invisibly, a data.frame of the annotated metric values.
#' @examples
#' dn_plot_containment()
#' @export
dn_plot_containment <- function(R = 100, broad = 14, narrow = 4,
                                offsets = c(0, 18, 34)) {
  x <- seq_len(R) - R / 2
  G <- function(cx, s) {
    M <- exp(-(outer(x - cx, x * 0, "+")^2 + outer(x * 0, x, "+")^2) /
               (2 * s^2))
    M / sum(M)
  }
  op <- graphics::par(mfrow = c(1, 3), mar = c(2, 2, 3, 1))
  on.exit(graphics::par(op))
  titles <- c("region in niche core", "region at niche margin",
              "region outside niche")
  out <- vector("list", 3)
  Zs <- G(0, broad)
  for (i in seq_along(offsets)) {
    Ze <- G(offsets[i], narrow)
    m <- dn_metrics(matrix(as.numeric(Zs), 1), matrix(as.numeric(Ze), 1))
    graphics::image(Zs, col = grDevices::grey.colors(64, 0.97, 0.25),
                    axes = FALSE, main = titles[i], useRaster = TRUE)
    graphics::contour(matrix(as.numeric(Ze), R, R), nlevels = 4, add = TRUE,
                      col = "#B2182B", drawlabels = FALSE, lwd = 1.5)
    graphics::mtext(sprintf("D = %.2f   DNCq = %.2f   DNCauc = %.2f",
                            m$D, m$DNC_q, m$DNC_auc),
                    side = 1, line = 0.5, cex = 0.75)
    out[[i]] <- data.frame(panel = titles[i], D = m$D[1, 1],
                           DNC_q = m$DNC_q[1, 1], DNC_auc = m$DNC_auc[1, 1])
  }
  invisible(do.call(rbind, out))
}

#' CDDR against area with the confound displayed (papers' Figure 7 style)
#'
#' Raw ratio against log area with the fitted confound line, the top
#' positive and negative residuals labelled -- the display that shows the
#' confound rather than hiding it.
#'
#' @param cddr_tab output of [dn_cddr()] (with `area` supplied).
#' @param area the area vector used.
#' @param labels optional region names for labelling extremes.
#' @param n_label how many extremes to label on each side (default 5).
#' @return invisibly, the labelled subset.
#' @examples
#' set.seed(1)
#' area <- 10^runif(80, 4, 6.7)
#' nd <- rpois(80, 30 + 200 / log10(area)); np <- rpois(80, 60)
#' x <- dn_cddr(nd, np, area)
#' dn_plot_cddr_area(x, area, labels = paste0("eco", 1:80))
#' @export
dn_plot_cddr_area <- function(cddr_tab, area, labels = NULL, n_label = 5) {
  la <- log10(area)
  res <- if (!is.null(cddr_tab$cddr_residual)) cddr_tab$cddr_residual else
    cddr_tab$residual
  fit <- attr(cddr_tab, "fit")
  if (is.null(fit)) fit <- cddr_tab$fit
  if (is.null(fit)) fit <- stats::lm(cddr_tab$cddr ~ la)
  r <- attr(cddr_tab, "confound_r")
  if (is.null(r)) r <- suppressWarnings(stats::cor(cddr_tab$cddr, la))
  graphics::plot(la, cddr_tab$cddr, pch = 16, col = "#00000066",
                 xlab = expression(log[10] ~ "area (km"^2 * ")"),
                 ylab = "climatic dark-diversity ratio")
  graphics::abline(fit, lwd = 2, col = "#B2182B")
  graphics::legend("topright", bty = "n",
                   legend = if (is.na(r)) "confound not estimable" else
                     sprintf("r = %.2f (confound)", r))
  lab_out <- NULL
  if (!is.null(labels)) {
    o <- order(res)
    ex <- c(utils::head(o, n_label), utils::tail(o, n_label))
    graphics::text(la[ex], cddr_tab$cddr[ex], labels[ex], cex = 0.6, pos = 3,
                   col = "#2166AC")
    lab_out <- data.frame(region = labels[ex], residual = res[ex])
  }
  invisible(lab_out)
}

#' Region-by-taxon matrix heatmap (papers' Figure 9 style)
#'
#' @param M numeric matrix, regions in rows, taxa in columns (e.g. CDDR).
#' @param main title.
#' @param values logical, print cell values (default TRUE).
#' @param key logical, draw a colour key in the right margin (default TRUE).
#' @return invisibly, `M`.
#' @examples
#' M <- matrix(runif(6 * 5), 6, 5,
#'             dimnames = list(paste("region", 1:6),
#'                             c("amph", "birds", "mamm", "plants", "rept")))
#' dn_plot_matrix(M, main = "CDDR by region and taxon")
#' @export
dn_plot_matrix <- function(M, main = "", values = TRUE, key = TRUE) {
  n <- nrow(M); k <- ncol(M)
  op <- graphics::par(mar = c(5, 12, 3, 5))
  on.exit(graphics::par(op))
  graphics::image(seq_len(k), seq_len(n), t(M[n:1, , drop = FALSE]),
                  col = .dn_heat(64), axes = FALSE, xlab = "", ylab = "",
                  main = main, zlim = c(0, max(M, na.rm = TRUE)))
  graphics::axis(1, at = seq_len(k), labels = colnames(M), las = 2,
                 cex.axis = 0.8)
  graphics::axis(2, at = seq_len(n), labels = rev(rownames(M)), las = 1,
                 cex.axis = 0.7)
  if (values)
    for (i in seq_len(n)) for (j in seq_len(k))
      graphics::text(j, n - i + 1, sprintf("%.2f", M[i, j]), cex = 0.6,
                     col = if (M[i, j] > 0.6 * max(M, na.rm = TRUE))
                       "white" else "black")
  if (key) {
    zmax <- max(M, na.rm = TRUE)
    usr <- graphics::par("usr")
    xw <- (usr[2] - usr[1]) * 0.035
    x0 <- usr[2] + (usr[2] - usr[1]) * 0.03
    ys <- seq(usr[3] + (usr[4] - usr[3]) * 0.15,
              usr[3] + (usr[4] - usr[3]) * 0.85, length.out = 65)
    cols <- .dn_heat(64)
    graphics::rect(x0, ys[-65], x0 + xw, ys[-1], col = cols, border = NA,
                   xpd = NA)
    graphics::rect(x0, ys[1], x0 + xw, ys[65], border = "grey40", xpd = NA)
    graphics::text(x0 + xw * 1.4, c(ys[1], ys[65]),
                   sprintf("%.2f", c(0, zmax)), cex = 0.65, adj = 0, xpd = NA)
  }
  invisible(M)
}

#' Stacked source-decomposition bars (papers' Figure 11 style)
#'
#' @param shares matrix or data.frame with columns `near`, `regional`,
#'   `distant` (proportions summing to 1), regions in rows.
#' @param main title.
#' @param legend_pos legend position, or `NA` to suppress it.
#' @param cex_names label size for the region names.
#' @return invisibly, `shares`.
#' @examples
#' s <- t(sapply(1:8, function(i) { x <- runif(3); x / sum(x) }))
#' colnames(s) <- c("near", "regional", "distant")
#' rownames(s) <- paste("region", 1:8)
#' dn_plot_sources(s)
#' @export
dn_plot_sources <- function(shares, main = "source of dark diversity",
                            legend_pos = "topright", cex_names = 0.7) {
  S <- as.matrix(shares[, c("near", "regional", "distant")])
  op <- graphics::par(mar = c(4.2, 12, 3, 1))
  on.exit(graphics::par(op))
  graphics::barplot(t(S[nrow(S):1, , drop = FALSE]), horiz = TRUE,
                    col = c("#B2182B", "#EF8A62", "#2166AC"), border = NA,
                    las = 1, cex.names = cex_names, xlim = c(0, 1.28),
                    xaxt = "n",
                    xlab = "proportion of climatic dark diversity",
                    main = main)
  graphics::axis(1, at = seq(0, 1, 0.25))
  if (!is.na(legend_pos))
    graphics::legend(legend_pos, fill = c("#B2182B", "#EF8A62", "#2166AC"),
                     border = NA, bty = "n", cex = 0.8,
                     legend = c("near", "regional", "distant"))
  invisible(shares)
}

#' ROC curves for metric comparison (papers' Figure 3 style)
#'
#' @param scores named list of score vectors on the same cases.
#' @param label logical/0-1 vector.
#' @param main title.
#' @return invisibly, data.frame of AUCs with 95% CIs.
#' @examples
#' set.seed(1)
#' lab <- rep(c(TRUE, FALSE), each = 300)
#' s1 <- ifelse(lab, rnorm(300, 1), rnorm(300, 0))
#' s2 <- s1 + rnorm(600, 0, 1)
#' dn_plot_roc(list(DNC_q = s1, D = s2), lab)
#' @export
dn_plot_roc <- function(scores, label, main = "introduced-range validation") {
  y <- as.logical(label)
  cols <- .dn_pal(length(scores))
  graphics::plot(0:1, 0:1, type = "l", lty = 2, col = "grey60",
                 xlab = "false positive rate", ylab = "true positive rate",
                 main = main)
  out <- vector("list", length(scores))
  for (i in seq_along(scores)) {
    s <- scores[[i]]
    th <- sort(unique(s), decreasing = TRUE)
    tpr <- vapply(th, function(t) mean(s[y] >= t), 0)
    fpr <- vapply(th, function(t) mean(s[!y] >= t), 0)
    graphics::lines(c(0, fpr, 1), c(0, tpr, 1), col = cols[i], lwd = 2)
    a <- dn_auc(s, y)
    out[[i]] <- data.frame(metric = names(scores)[i], auc = a$auc,
                           lo = a$ci[1], hi = a$ci[2])
  }
  res <- do.call(rbind, out)
  graphics::legend("bottomright", bty = "n", lwd = 2, col = cols,
                   legend = sprintf("%s  AUC %.3f", res$metric, res$auc))
  invisible(res)
}

#' Multi-region threshold-sweep figure (papers' Figure 10 style)
#'
#' Draws every region's observed/expected curve on a log axis with the
#' cross-region median, colouring curves by their direction at the lowest
#' threshold -- the display in which the bimodal-vs-unimodal crossing is
#' visible at a glance.
#'
#' @param sweeps a list of [dn_sweep()] results (one per region).
#' @param main title.
#' @return invisibly, the share of regions whose direction reverses.
#' @examples
#' set.seed(1)
#' sweeps <- lapply(1:15, function(i) {
#'   obs <- c(rbeta(250, 0.3, 3), rbeta(80, 4, 1.5))
#'   dn_sweep(obs, rbeta(330, 1, 2.2), 0.08 + 0.1 * runif(330))
#' })
#' dn_plot_sweeps(sweeps)
#' @export
dn_plot_sweeps <- function(sweeps, main = "aggregate null test across thresholds") {
  get_ratio <- function(s) s$n_admissible_obs / s$expected
  qx <- sweeps[[1]]$threshold_q
  Rm <- sapply(sweeps, get_ratio)
  ok <- is.finite(Rm)
  graphics::plot(range(qx), range(Rm[ok]), type = "n", log = "y",
                 xlab = "threshold (quantile of coverage distribution)",
                 ylab = "observed / expected admissible species",
                 main = main)
  flips <- 0L
  for (j in seq_len(ncol(Rm))) {
    r <- Rm[, j]; f <- is.finite(r)
    if (!any(f)) next
    col <- if (r[which(f)[1]] < 1) "#B2182B55" else "#2166AC55"
    if (any(r[f] < 1) && any(r[f] > 1)) flips <- flips + 1L
    graphics::lines(qx[f], r[f], col = col)
  }
  med <- apply(Rm, 1L, stats::median, na.rm = TRUE)
  graphics::lines(qx, med, lwd = 3)
  graphics::abline(h = 1, lty = 2)
  graphics::legend("topleft", bty = "n",
                   legend = sprintf("%d / %d regions change direction",
                                    flips, ncol(Rm)))
  invisible(flips / ncol(Rm))
}

# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------

#' Extremes table (papers' Table 6 style)
#'
#' The `n` highest and lowest regions by the area-adjusted residual, with the
#' raw ratio shown alongside -- never the ranking key.
#'
#' @param region region names.
#' @param cddr_tab output of [dn_cddr()] (with `area`).
#' @param n_obs,n_dark optional richness columns to carry.
#' @param n rows per side (default 8).
#' @return data.frame, top block then bottom block.
#' @examples
#' set.seed(1)
#' area <- 10^runif(40, 4, 6.5)
#' x <- dn_cddr(rpois(40, 40), rpois(40, 60), area)
#' dn_table_extremes(paste0("eco", 1:40), x)
#' @export
dn_table_extremes <- function(region, cddr_tab, n_obs = NULL, n_dark = NULL,
                              n = 8) {
  res <- if (!is.null(cddr_tab$cddr_residual)) cddr_tab$cddr_residual else
    cddr_tab$residual
  if (is.null(res))
    stop("cddr_tab needs a 'cddr_residual' (or 'residual') column: ",
         "call dn_cddr() with an 'area' argument")
  d <- data.frame(region = region, cddr = cddr_tab$cddr, residual = res)
  if (!is.null(n_obs)) d$n_obs <- n_obs
  if (!is.null(n_dark)) d$n_dark <- n_dark
  o <- order(-d$residual)
  rbind(utils::head(d[o, ], n), utils::tail(d[o, ], n))
}

#' Export a table as CSV, Markdown or LaTeX
#'
#' @param x data.frame.
#' @param file output path; format from extension (`.csv`, `.md`, `.tex`).
#' @param digits rounding for numeric columns (default 3).
#' @param caption optional caption (md/tex).
#' @return invisibly, the path.
#' @examples
#' tf <- tempfile(fileext = ".md")
#' dn_table_export(data.frame(a = 1:2, b = c(0.123456, 2)), tf,
#'                 caption = "demo")
#' cat(readLines(tf), sep = "\n")
#' @export
dn_table_export <- function(x, file, digits = 3, caption = NULL) {
  num <- vapply(x, is.numeric, TRUE)
  x[num] <- lapply(x[num], round, digits)
  ext <- tolower(tools::file_ext(file))
  if (ext == "csv") {
    utils::write.csv(x, file, row.names = FALSE)
  } else if (ext == "md") {
    hdr <- paste("|", paste(names(x), collapse = " | "), "|")
    sep <- paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|")
    rows <- apply(x, 1L, function(r) paste("|", paste(r, collapse = " | "),
                                           "|"))
    txt <- c(if (!is.null(caption)) paste0("**", caption, "**\n"),
             hdr, sep, rows)
    writeLines(txt, file)
  } else if (ext == "tex") {
    al <- paste(ifelse(num, "r", "l"), collapse = "")
    rows <- apply(x, 1L, function(r) paste(paste(r, collapse = " & "),
                                           "\\\\"))
    txt <- c("\\begin{table}[ht]\\centering",
             if (!is.null(caption)) sprintf("\\caption{%s}", caption),
             sprintf("\\begin{tabular}{%s}\\hline", al),
             paste(paste(names(x), collapse = " & "), "\\\\ \\hline"),
             rows, "\\hline\\end{tabular}\\end{table}")
    writeLines(txt, file)
  } else stop("unsupported extension: ", ext)
  invisible(file)
}
