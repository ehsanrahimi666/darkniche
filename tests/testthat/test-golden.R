# Golden tests: every core computation is compared against fixtures generated
# by the reference Python implementation (ddarkniche pipeline v2.1) on the same
# synthetic inputs. Tolerances reflect the float32 storage of the reference
# outputs (eps ~ 1.2e-7 on values in [0, 1]).

gold <- function(f) file.path(system.file("extdata", "golden",
                                          package = "darkniche"), f)
read_mat <- function(f, header = FALSE)
  as.matrix(utils::read.csv(gold(f), header = header))

# Python densities are flattened C-order (row-major) on an (R, R) array whose
# first axis is x. R matrices are column-major; reading the flat vector into
# an R x R matrix with byrow = TRUE reproduces [x, y] indexing.
py_flat_to_mat <- function(v, R) matrix(v, R, R, byrow = TRUE)
mat_to_py_flat <- function(M) as.numeric(t(M))

# helper used only by the golden test: shift null with explicit (xi, yi)
# targets and C-order flattening to match the Python fixture exactly.
null_shift_py_order <- function(Z_sp, Ze_flatC, Z_global, q, statistic,
                                targets_xi, targets_yi) {
  cc <- darkniche:::.centroid_cell(Z_sp)
  E0 <- matrix(dn_weights(Z_sp, q, statistic), nrow(Z_sp), ncol(Z_sp))
  B <- length(targets_xi)
  W <- matrix(0, length(E0), B)
  for (b in seq_len(B)) {
    Wb <- darkniche:::.roll2(E0, targets_xi[b] - cc[1], targets_yi[b] - cc[2])
    W[, b] <- as.numeric(t(Wb))          # C-order flatten
  }
  as.matrix(Ze_flatC) %*% W
}
jsonlite_min2 <- function(txt) {
  txt <- paste(txt, collapse = "")
  list(scale = as.numeric(sub('.*"scale":\\s*([0-9eE+.-]+).*', "\\1", txt)),
       rel_floor = as.numeric(sub('.*"rel_floor":\\s*([0-9eE+.-]+).*', "\\1",
                                  txt)))
}
sp_from_json <- function() {
  j <- jsonlite_min(readLines(gold("pcspace.json"), warn = FALSE))
  structure(list(x_edges = j$x_edges, y_edges = j$y_edges,
                 R = as.integer(j$R), bw_cells = j$bw_cells),
            class = "dn_space")
}
# minimal JSON reader (numbers + arrays only) to avoid a jsonlite dependency
jsonlite_min <- function(txt) {
  txt <- paste(txt, collapse = "")
  get_num <- function(key) as.numeric(
    sub(sprintf('.*"%s":\\s*([0-9eE+.-]+).*', key), "\\1", txt))
  get_arr <- function(key) {
    m <- regmatches(txt, regexpr(sprintf('"%s":\\s*\\[[^]]*\\]', key), txt))
    as.numeric(strsplit(gsub(sprintf('"%s":\\s*\\[|\\]', key), "", m),
                        ",")[[1]])
  }
  list(R = get_num("R"), bw_cells = get_num("bw_cells"),
       x_edges = get_arr("x_edges"), y_edges = get_arr("y_edges"))
}

test_that("dn_space reproduces the Python lattice", {
  scores <- read_mat("scores.csv", header = TRUE)
  spg <- sp_from_json()
  spr <- dn_space(scores, R = spg$R, bw_frac = 0.02)
  expect_equal(spr$x_edges, spg$x_edges, tolerance = 1e-9)
  expect_equal(spr$y_edges, spg$y_edges, tolerance = 1e-9)
  expect_equal(spr$bw_cells, spg$bw_cells, tolerance = 1e-12)
})

test_that("dn_density matches scipy gaussian_filter to float32", {
  spg <- sp_from_json()
  Zg_py <- read_mat("Z_global.csv")            # (R, R), C-order rows = x
  scores <- read_mat("scores.csv", header = TRUE)
  Zg_r <- dn_density(scores, spg)
  expect_lt(max(abs(Zg_r - Zg_py)), 2e-7)   # py rows = x axis, same layout
  for (i in 1:8) {
    p <- read_mat(sprintf("species_pts_%d.csv", i), header = TRUE)
    Zi_py <- py_flat_to_mat(as.numeric(t(read_mat("Zs.csv")[i, ])), spg$R)
    Zi_r <- dn_density(p, spg)
    expect_lt(max(abs(Zi_r - Zi_py)), 2e-7)
  }
})

test_that("dn_envelope and dn_entryq match", {
  spg <- sp_from_json()
  Zs <- read_mat("Zs.csv")
  env_py <- as.numeric(t(read_mat("env095_sp1.csv")))
  z1 <- mat_to_py_flat(py_flat_to_mat(as.numeric(t(Zs[1, ])), spg$R))
  # z1 is now in Python flat order; envelope in the same order
  expect_equal(as.numeric(dn_envelope(z1, 0.95)), env_py)
  eq_py <- as.numeric(t(read_mat("entryq_sp1.csv")))
  expect_lt(max(abs(dn_entryq(z1) - eq_py)), 1e-9)
})

test_that("dn_metrics matches the Python block for every metric", {
  spg <- sp_from_json()
  Zs <- read_mat("Zs.csv"); Ze <- read_mat("Ze.csv")
  sp_mu <- read_mat("sp_mu.csv"); eco_mu <- read_mat("eco_mu.csv")
  sp_cov <- read_mat("sp_cov.csv")
  # centres in Python flat (C) order: x varies slowest -> expand.grid with y
  # fastest. Build G to match the Ze flattening (Python C-order).
  ctr <- dn_centres(spg)
  centres_py <- list(x = ctr$x, y = ctr$y)
  res <- dn_metrics(Zs, Ze, q = 0.95, sp_mu = sp_mu, sp_cov = sp_cov,
                      eco_mu = eco_mu,
                      centres = list(
                        x = rep(ctr$x, each = length(ctr$y)) * NA, # unused
                        y = NA))
  # CFI_w needs the correct cell ordering; compute it directly here instead:
  G <- cbind(rep(ctr$x, each = spg$R), rep(ctr$y, times = spg$R))
  CFIw <- matrix(NA_real_, nrow(Zs), nrow(Ze))
  for (i in seq_len(nrow(Zs))) {
    C <- matrix(sp_cov[i, ], 2, 2, byrow = TRUE)
    Ci <- solve(C + 1e-12 * diag(2))
    dv <- sweep(G, 2L, sp_mu[i, ], "-")
    wcell <- stats::pchisq(rowSums((dv %*% Ci) * dv), df = 2,
                           lower.tail = FALSE)
    CFIw[i, ] <- Ze %*% wcell
  }
  res$CFI_w <- CFIw
  for (nm in c("D", "I", "BC", "DNC_q", "NSI_q", "DNC_auc", "DNC_min",
               "NSI_min", "NOV", "mahal_d2", "CFI", "CFI_w")) {
    ref <- read_mat(sprintf("metric_%s.csv", nm))
    tol <- if (nm == "mahal_d2") 1e-3 else 5e-6
    expect_lt(max(abs(res[[nm]] - ref), na.rm = TRUE), tol,
              label = paste("metric", nm))
  }
})

test_that("containment phenomenon reproduces (paper Figure 1)", {
  Zs <- read_mat("Zs.csv"); Ze <- read_mat("Ze.csv")
  res <- dn_metrics(Zs, Ze, q = 0.95)
  # species 7 (engulfing niche) x region 4 (tiny nested region)
  expect_gt(res$DNC_q[7, 4], 0.999)
  expect_gt(res$DNC_auc[7, 4], 0.85)
  expect_lt(res$D[7, 4], 0.15)
})

test_that("environmental-shift null matches with injected targets", {
  spg <- sp_from_json()
  Zs <- read_mat("Zs.csv"); Ze <- read_mat("Ze.csv")
  Zg <- read_mat("Z_global.csv")
  targets0 <- as.integer(readLines(gold("shift_targets.csv"), warn = FALSE))  # 0-based, C-order
  R <- spg$R
  # convert Python C-order flat index (x major) to R column-major on [x, y]:
  xi <- targets0 %/% R; yi <- targets0 %% R
  targets_r <- xi + yi * R + 1L
  Zsp1 <- py_flat_to_mat(as.numeric(t(Zs[1, ])), R)
  Zgm <- Zg  # rows = x axis already
  # Ze rows must stay in Python flat order to match the rolled weight vector
  # ordering used internally: convert the rolled matrix to Python order.
  # Simplest: compute with matrices in [x, y] layout and convert Ze rows too.
  Ze_r <- t(apply(Ze, 1L, function(v) mat_to_py_flat(t(py_flat_to_mat(v, R)))))
  # identity: Ze rows are already C-order over [x, y]; our .roll2 output is
  # column-major over [x, y]; flatten it C-order for the product:
  nq <- null_shift_py_order(Zsp1, Ze, Zgm, q = 0.95, statistic = "q",
                            targets_xi = xi + 1L, targets_yi = yi + 1L)
  ref_q <- read_mat("shiftnull_q_sp1.csv")
  expect_lt(max(abs(nq - ref_q)), 5e-6)
  na <- null_shift_py_order(Zsp1, Ze, Zgm, q = 0.95, statistic = "auc",
                            targets_xi = xi + 1L, targets_yi = yi + 1L)
  ref_a <- read_mat("shiftnull_auc_sp1.csv")
  # the Python reference stores the AUC weight array as float32; its own code
  # documents ~4e-6 quantisation on individual weights, which accumulates to
  # ~1e-4 on the summed statistic. The R implementation is float64 throughout.
  expect_lt(max(abs(na - ref_a)), 2e-4)
})


test_that("ses flooring matches", {
  fx <- utils::read.csv(gold("ses_fixture.csv"))
  sc <- jsonlite_min2(readLines(gold("ses_scale.json"), warn = FALSE))
  got <- dn_ses(fx$obs, fx$null_mean, fx$null_sd, scale = sc$scale,
             rel_floor = sc$rel_floor)
  expect_equal(got, fx$ses, tolerance = 1e-9)
})

test_that("aggregate test and sweep match, including reliability flags", {
  inp <- utils::read.csv(gold("agg_inputs.csv"))
  ref <- utils::read.csv(gold("aggregate_sweep.csv"))
  for (i in seq_len(nrow(ref))) {
    r <- dn_aggregate(inp$obs, inp$null_mean, inp$null_sd,
                        ref$threshold[i])
    expect_equal(r$n_admissible_obs, ref$n_admissible_obs[i])
    expect_equal(r$expected, ref$expected[i], tolerance = 1e-8)
    expect_equal(r$reliable, ref$reliable[i] %in% c("True", TRUE, "TRUE"))
    if (ref$reliable[i] %in% c("True", TRUE, "TRUE")) {
      expect_equal(r$ses, ref$ses[i], tolerance = 1e-8)
      expect_equal(r$p, ref$p[i], tolerance = 1e-8)
      expect_equal(r$enrichment, ref$enrichment[i], tolerance = 1e-8)
    }
  }
})

test_that("breadth-matched null matches", {
  V <- read_mat("bnull_V.csv")
  br <- as.numeric(readLines(gold("bnull_breadth.csv"), warn = FALSE))
  res <- dn_null_breadth(V, br, rep("a", nrow(V)), K = 10L)
  expect_lt(max(abs(res$null_mean - read_mat("bnull_mu.csv")),
                na.rm = TRUE), 1e-6)
  expect_lt(max(abs(res$null_sd - read_mat("bnull_sd.csv")),
                na.rm = TRUE), 1e-6)
  expect_lt(max(abs(res$ses - read_mat("bnull_ses.csv")),
                na.rm = TRUE), 1e-4)
})
