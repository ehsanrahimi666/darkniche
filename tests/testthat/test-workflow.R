# Behavioural tests for the layers built on top of the verified core:
# bundled data, the report object, the direction census, plots and exporters.

test_that("bundled ecoregion data loads and matches the published summary", {
  eco <- dn_data("ecoregion")
  expect_equal(nrow(eco), 619)
  expect_true(all(c("region", "n_obs", "n_dark", "cddr", "area_km2") %in%
                    names(eco)))
  expect_equal(round(min(eco$cddr), 3), 0.022)
  expect_equal(round(max(eco$cddr), 3), 0.924)
  expect_equal(round(stats::median(eco$cddr), 3), 0.579)
  expect_equal(sum(eco$cddr > 0.5), 411)
  expect_equal(round(stats::cor(eco$cddr, log10(eco$area_km2)), 3), -0.441)
})

test_that("every bundled dataset is readable and non-empty", {
  for (nm in c("ecoregion", "sweep", "sources", "taxon", "region_taxon",
               "auc", "containment")) {
    d <- dn_data(nm)
    expect_s3_class(d, "data.frame")
    expect_gt(nrow(d), 0)
  }
  expect_error(dn_data("nonexistent"))
})

test_that("published validation AUCs are reproduced from bundled data", {
  a <- dn_data("auc")
  get <- function(m) round(a$auc[a$metric == m], 3)
  expect_equal(get("DNC_q"), 0.705)
  expect_equal(get("DNC_auc"), 0.702)
  expect_equal(get("D"), 0.688)
  # the algebraic identity DNC_min == D must hold in the published output
  expect_equal(get("DNC_min"), get("D"))
  expect_gt(get("DNC_q"), get("D"))
})

test_that("direction census reproduces the published threshold reversal", {
  cen_all <- dn_direction_census(dn_data("sweep"), reliable_only = FALSE)
  expect_equal(cen_all$census$n_enriched, c(71, 344, 538, 602))
  expect_equal(cen_all$census$n_depleted, c(548, 275, 81, 17))
  expect_equal(cen_all$n_changed, 553)
  expect_equal(cen_all$n_regions, 619)
  cen_rel <- dn_direction_census(dn_data("sweep"), reliable_only = TRUE)
  expect_equal(cen_rel$n_changed, 439)
  expect_equal(cen_rel$n_regions, 452)
  expect_gt(cen_rel$frac_changed, 0.95)
})

test_that("dn_report assembles and prints without error", {
  rep <- dn_report(dn_data("ecoregion"), dn_data("sweep"))
  expect_s3_class(rep, "dn_report")
  expect_equal(nrow(rep$table), 619)
  expect_equal(round(rep$confound_r, 3), -0.441)
  expect_true(all(c("rank_raw", "rank_adjusted", "rank_shift") %in%
                    names(rep$table)))
  # the adjustment must actually move regions
  expect_gt(stats::median(abs(rep$table$rank_shift)), 10)
  expect_output(print(rep), "area confound")
  expect_output(print(rep), "change direction")
  # report works without a sweep
  expect_silent(r2 <- dn_report(dn_data("ecoregion")))
})

test_that("dn_cddr residual is orthogonal to log area by construction", {
  eco <- dn_data("ecoregion")
  x <- dn_cddr(eco$n_dark, eco$n_obs, eco$area_km2)
  expect_lt(abs(stats::cor(x$cddr_residual, log10(eco$area_km2))), 1e-10)
  expect_equal(x$cddr, eco$cddr, tolerance = 1e-9)
})

test_that("dn_dark applies the conjunction rule and respects absence", {
  p <- data.frame(species = c("a", "a", "b", "b"),
                  region = c("r1", "r2", "r1", "r2"),
                  present = c(TRUE, FALSE, FALSE, FALSE),
                  p_establish = c(0.9, 0.9, 0.9, 0.05),
                  n2_ses = c(5, 5, 1, 5))
  out <- dn_dark(p, p_min = 0.177, ses_cols = "n2_ses")
  expect_equal(out$admissible, c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(out$dark, c(FALSE, TRUE, FALSE, FALSE))
})

test_that("aggregate test suppresses unreliable cells", {
  set.seed(1)
  obs <- runif(400, 0, 0.2)
  r <- dn_aggregate(obs, rep(0.02, 400), rep(0.01, 400), threshold = 0.95)
  expect_false(r$reliable)
  expect_true(is.na(r$ses))
  expect_true(is.na(r$p))
  expect_true(is.na(r$enrichment))
  # too few usable species -> nothing is reported
  r2 <- dn_aggregate(runif(10), runif(10), runif(10, .1, .2), 0.5)
  expect_false(r2$reliable)
  expect_true(is.na(r2$expected))
})

test_that("threshold sweep flags reversal and survives degenerate input", {
  set.seed(2)
  # bimodal observed coverage against a unimodal null: the configuration in
  # which the two distributions must cross somewhere in the sweep
  obs <- c(stats::rbeta(300, 0.3, 5), stats::rbeta(100, 9, 1.1))
  mu <- stats::rbeta(400, 2.5, 2.5)
  sdv <- 0.18 + 0.06 * stats::runif(400)
  sw <- dn_sweep(obs, mu, sdv)
  expect_s3_class(sw, "dn_sweep_result")
  expect_equal(attr(sw, "n_reliable"), nrow(sw))
  expect_true(attr(sw, "direction_change"))
  expect_output(print(sw), "DIRECTION REVERSES")
  # a threshold-stable region must not be flagged: observed coverage is the
  # null shifted up by a constant, so it is enriched at every threshold
  set.seed(2)
  mu2 <- stats::rbeta(1000, 3, 3)
  sw3 <- dn_sweep(pmin(mu2 + 0.05, 1), mu2, rep(0.05, 1000))
  expect_equal(attr(sw3, "n_reliable"), nrow(sw3))
  expect_false(attr(sw3, "direction_change"))
  expect_true(all(sw3$n_admissible_obs / sw3$expected > 1))
  expect_output(print(sw3), "direction is stable")
  # fewer than 30 usable species: nothing is reliable, printing must not crash
  sw2 <- dn_sweep(stats::runif(20), stats::runif(20),
                  stats::runif(20, .05, .1))
  expect_equal(attr(sw2, "n_reliable"), 0)
  expect_output(print(sw2), "no reliable thresholds")
})

test_that("the sweep verdict never uses cells the test calls unreliable", {
  set.seed(6)
  obs <- stats::runif(400, 0, 0.25)
  # nulls so tight that high thresholds have expected counts far below 5
  sw <- dn_sweep(obs, rep(0.02, 400), rep(0.008, 400))
  ratio <- attr(sw, "ratio")
  unreliable <- !sw$reliable
  expect_true(any(unreliable))
  # ratios exist for unreliable rows, but must not drive the verdict
  if (all(unreliable)) {
    expect_false(attr(sw, "direction_change"))
    expect_equal(attr(sw, "n_reliable"), 0)
  }
  expect_true(all(is.na(sw$ses[unreliable])))
  expect_true(all(is.na(sw$enrichment[unreliable])))
})

test_that("DeLong test agrees with the AUC difference and is symmetric", {
  set.seed(3)
  lab <- rep(c(TRUE, FALSE), each = 400)
  s1 <- ifelse(lab, stats::rnorm(400, 1), stats::rnorm(400, 0))
  s2 <- s1 + stats::rnorm(800, 0, 1.2)
  d <- dn_delong(s1, s2, lab)
  expect_equal(d$auc1, dn_auc(s1, lab)$auc, tolerance = 1e-12)
  expect_equal(d$auc2, dn_auc(s2, lab)$auc, tolerance = 1e-12)
  d2 <- dn_delong(s2, s1, lab)
  expect_equal(d$z, -d2$z, tolerance = 1e-10)
  expect_equal(d$p, d2$p, tolerance = 1e-12)
  # AUC of a perfect and a useless score
  expect_equal(dn_auc(as.numeric(lab), lab)$auc, 1)
  expect_equal(dn_auc(rep(1, 800), lab)$auc, 0.5)
})

test_that("matched negatives respect realm and exclude reached regions", {
  regs <- data.frame(region = paste0("r", 1:10),
                     realm = rep(c("PA", "NT"), 5),
                     area = 10^seq(4, 6, length.out = 10))
  pos <- data.frame(species = "sp1", region = "r1")
  out <- dn_match(pos, regs, reached = list(sp1 = c("r1", "r3")), ratio = 3)
  expect_true(all(out$region[!out$label] %in%
                    regs$region[regs$realm == "PA"]))
  expect_false("r3" %in% out$region[!out$label])
  expect_equal(sum(out$label), 1)
})

test_that("calibration returns a usable threshold", {
  set.seed(4)
  lab <- rep(c(TRUE, FALSE), each = 300)
  p <- ifelse(lab, stats::rbeta(300, 5, 2), stats::rbeta(300, 2, 5))
  th <- dn_calibrate(p, lab)
  expect_true(th > 0 && th < 1)
  expect_gt(mean(p[lab] >= th), 0.6)
  expect_gt(mean(p[!lab] < th), 0.6)
})

test_that("plot functions run on bundled data without error", {
  f <- tempfile(fileext = ".png")
  grDevices::png(f, 900, 700)
  on.exit({grDevices::dev.off(); unlink(f)}, add = TRUE)
  eco <- dn_data("ecoregion")
  x <- dn_cddr(eco$n_dark, eco$n_obs, eco$area_km2)
  expect_silent(dn_plot_cddr_area(x, eco$area_km2))
  expect_silent(dn_plot_cddr_area(x, eco$area_km2, labels = eco$region))
  src <- dn_data("sources")[1:8, ]
  expect_silent(dn_plot_sources(src))
  rt <- dn_data("region_taxon")
  M <- stats::xtabs(cddr ~ region + taxon, data = rt)[1:6, ]
  expect_silent(dn_plot_matrix(as.matrix(M)))
  expect_silent(dn_plot_matrix(as.matrix(M), values = FALSE, key = FALSE))
  sw <- lapply(1:5, function(i)
    dn_sweep(stats::runif(200), stats::runif(200, 0, .3),
             stats::runif(200, .05, .2)))
  expect_silent(dn_plot_sweeps(sw))
  expect_invisible(dn_plot_containment(R = 40))
})

test_that("table exporters write all three formats", {
  d <- data.frame(a = c("x", "y"), b = c(1.23456, 2))
  for (ext in c(".csv", ".md", ".tex")) {
    f <- tempfile(fileext = ext)
    dn_table_export(d, f, digits = 2, caption = "demo")
    expect_true(file.exists(f))
    expect_gt(length(readLines(f, warn = FALSE)), 1)
    unlink(f)
  }
  expect_error(dn_table_export(d, tempfile(fileext = ".xyz")))
})

test_that("dn_table_extremes needs an area-adjusted residual", {
  eco <- dn_data("ecoregion")
  bad <- dn_cddr(eco$n_dark, eco$n_obs)          # no area supplied
  expect_error(dn_table_extremes(eco$region, bad), "residual")
  good <- dn_cddr(eco$n_dark, eco$n_obs, eco$area_km2)
  ex <- dn_table_extremes(eco$region, good, n = 5)
  expect_equal(nrow(ex), 10)
  expect_true(ex$residual[1] > ex$residual[nrow(ex)])
})

test_that("degenerate inputs fail loudly rather than returning nonsense", {
  # AUC and its relatives are undefined without both classes
  expect_error(dn_auc(stats::runif(10), rep(TRUE, 10)), "both classes")
  expect_error(dn_auc(stats::runif(10), rep(FALSE, 10)), "both classes")
  expect_error(dn_auc(stats::runif(10), rep(c(TRUE, FALSE), 4)),
               "same length")
  expect_error(dn_delong(stats::runif(10), stats::runif(10), rep(TRUE, 10)),
               "both classes")
  expect_error(dn_calibrate(stats::runif(10), rep(TRUE, 10)), "both classes")
  # a perfect and a reversed score
  expect_equal(dn_auc(c(2, 3, 1, 0), c(TRUE, TRUE, FALSE, FALSE))$auc, 1)
  expect_equal(dn_auc(c(0, 1, 3, 2), c(TRUE, TRUE, FALSE, FALSE))$auc, 0)
})

test_that("constant ratios and empty inputs are handled without warnings", {
  x <- dn_cddr(rep(0, 10), stats::rpois(10, 50), 10^stats::runif(10, 4, 6))
  expect_true(is.na(attr(x, "confound_r")))
  expect_silent(dn_cddr(rep(0, 10), stats::rpois(10, 50),
                        10^stats::runif(10, 4, 6)))
  rep0 <- dn_report(data.frame(region = letters[1:10],
                               n_obs = stats::rpois(10, 50),
                               n_dark = rep(0, 10),
                               area_km2 = 10^stats::runif(10, 4, 6)))
  expect_output(print(rep0), "not estimable")
  cen <- dn_direction_census(data.frame(region = character(),
                                        quantile = numeric(),
                                        enrichment = numeric()))
  expect_equal(cen$n_regions, 0)
})

test_that("densities survive empty, singular and out-of-range point sets", {
  sp <- dn_space(cbind(stats::rnorm(500), stats::rnorm(500)), R = 20)
  expect_equal(sum(dn_density(matrix(numeric(0), 0, 2), sp)), 0)
  expect_equal(sum(dn_density(cbind(c(99, -99), c(99, -99)), sp)), 0)
  z1 <- dn_density(cbind(0, 0), sp)
  expect_equal(sum(z1), 1, tolerance = 1e-12)
  expect_false(any(dn_envelope(rep(0, 400), 0.95)))
  expect_true(all(dn_entryq(rep(0, 400)) == 1))
  expect_true(all(dn_envelope(stats::runif(400), 1)))
})
