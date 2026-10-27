"""Generate golden test fixtures for the R port from the actual ddarkniche code.

Run: python tools/make_golden.py. If the ddarkniche Python package is not installed,
set DDARKNICHE_PATH to the folder that contains it.
"""
import sys, os, json
from pathlib import Path
import numpy as np
if os.environ.get("DDARKNICHE_PATH"):
    sys.path.insert(0, os.environ["DDARKNICHE_PATH"])
from ddarkniche.niche import build_pcspace, kde_grid, envelope_mask, density_moments
from ddarkniche.metrics import pair_metrics, envelopes_flat, entry_quantiles
from ddarkniche.nulls import shift_null_dnc, aggregate_admissibility_test, _availability_weights, breadth_matched_null
from ddarkniche.utils import ses

OUT = str(Path(__file__).resolve().parents[1] / "inst" / "extdata" / "golden")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(42)

# ---- synthetic global scores and lattice ---------------------------------- #
n = 20000
comp = rng.integers(0, 3, n)
mu = np.array([[-1.5, 0.5], [1.0, -1.0], [0.2, 1.8]])
sd = np.array([[0.8, 0.6], [1.1, 0.9], [0.5, 0.7]])
scores = mu[comp] + rng.standard_normal((n, 2)) * sd[comp]
R = 60
sp = build_pcspace(scores, R=R, bw_frac=0.02)
np.savetxt(f"{OUT}/scores.csv", scores, delimiter=",", fmt="%.10e", header="pc1,pc2", comments="")
with open(f"{OUT}/pcspace.json", "w") as f:
    json.dump({"R": R, "bw_cells": sp.bw_cells,
               "x_edges": sp.x_edges.tolist(), "y_edges": sp.y_edges.tolist()}, f)

# global availability density
Zg = kde_grid(scores, sp)
np.savetxt(f"{OUT}/Z_global.csv", Zg.astype(np.float64), delimiter=",", fmt="%.10e")

# ---- species and region point sets ---------------------------------------- #
def cloud(cx, cy, sx, sy, m, rho=0.0):
    z = rng.standard_normal((m, 2))
    x = cx + sx * z[:, 0]
    y = cy + sy * (rho * z[:, 0] + np.sqrt(1 - rho**2) * z[:, 1])
    return np.column_stack([x, y])

species_pts = [
    cloud(-1.5, 0.5, 0.9, 0.7, 400),          # broad, near comp 1
    cloud(1.0, -1.0, 0.3, 0.3, 120),           # narrow
    cloud(0.2, 1.8, 0.6, 0.5, 250, rho=0.4),
    cloud(-0.5, -0.5, 1.4, 1.2, 800),          # very broad generalist
    cloud(2.2, 0.8, 0.25, 0.4, 60),            # marginal narrow
    cloud(-2.5, 1.5, 0.5, 0.5, 150),
    cloud(0.0, 0.0, 2.0, 1.8, 1500),           # engulfing niche (containment)
    cloud(1.5, 1.5, 0.4, 0.3, 90),
]
region_pts = [
    cloud(-1.4, 0.4, 0.35, 0.30, 300),
    cloud(1.1, -0.9, 0.30, 0.25, 300),
    cloud(0.3, 1.7, 0.25, 0.25, 300),
    cloud(-0.2, -0.1, 0.15, 0.15, 300),        # tiny region inside sp7 niche
    cloud(2.4, 1.0, 0.30, 0.30, 300),
]
Zs = np.stack([kde_grid(p, sp).ravel() for p in species_pts])
Ze = np.stack([kde_grid(p, sp).ravel() for p in region_pts])
np.savetxt(f"{OUT}/Zs.csv", Zs.astype(np.float64), delimiter=",", fmt="%.10e")
np.savetxt(f"{OUT}/Ze.csv", Ze.astype(np.float64), delimiter=",", fmt="%.10e")
for i, p in enumerate(species_pts):
    np.savetxt(f"{OUT}/species_pts_{i+1}.csv", p, delimiter=",", fmt="%.10e", header="pc1,pc2", comments="")
for i, p in enumerate(region_pts):
    np.savetxt(f"{OUT}/region_pts_{i+1}.csv", p, delimiter=",", fmt="%.10e", header="pc1,pc2", comments="")

# ---- moments for CFI ------------------------------------------------------- #
mus, covs, eco_mu = [], [], []
for i in range(Zs.shape[0]):
    m0, C = density_moments(Zs[i].reshape(R, R), sp)
    mus.append(m0); covs.append(C)
for j in range(Ze.shape[0]):
    m0, _ = density_moments(Ze[j].reshape(R, R), sp)
    eco_mu.append(m0)
sp_mu = np.array(mus); sp_cov = np.array(covs); eco_mu = np.array(eco_mu)
np.savetxt(f"{OUT}/sp_mu.csv", sp_mu, delimiter=",", fmt="%.10e")
np.savetxt(f"{OUT}/sp_cov.csv", sp_cov.reshape(len(covs), 4), delimiter=",", fmt="%.10e")
np.savetxt(f"{OUT}/eco_mu.csv", eco_mu, delimiter=",", fmt="%.10e")

# ---- pair metrics ----------------------------------------------------------- #
res = pair_metrics(Zs, Ze, q=0.95, sp_mu=sp_mu, sp_cov=sp_cov, eco_mu=eco_mu,
                   pc_centres=(sp.x_centres, sp.y_centres))
for k, v in res.items():
    np.savetxt(f"{OUT}/metric_{k}.csv", v.astype(np.float64), delimiter=",", fmt="%.8e")

# envelope + entry quantiles for species 1 and the engulfing species 7
np.savetxt(f"{OUT}/env095_sp1.csv", envelopes_flat(Zs[0], 0.95).astype(int), delimiter=",", fmt="%d")
np.savetxt(f"{OUT}/entryq_sp1.csv", entry_quantiles(Zs[0]), delimiter=",", fmt="%.10e")

# ---- shift null (targets exported for exact replication) -------------------- #
B = 199
w = _availability_weights(Zg)
r0 = np.random.default_rng(0)
targets = r0.choice(w.size, size=B, p=w)      # replicates the call inside shift_null_dnc
np.savetxt(f"{OUT}/shift_targets.csv", targets, fmt="%d")
null_q = shift_null_dnc(Zs[0].reshape(R, R), Ze, Zg, q=0.95, B=B,
                        rng=np.random.default_rng(0), statistic="q")
null_auc = shift_null_dnc(Zs[0].reshape(R, R), Ze, Zg, q=0.95, B=B,
                          rng=np.random.default_rng(0), statistic="auc")
np.savetxt(f"{OUT}/shiftnull_q_sp1.csv", null_q.astype(np.float64), delimiter=",", fmt="%.8e")
np.savetxt(f"{OUT}/shiftnull_auc_sp1.csv", null_auc.astype(np.float64), delimiter=",", fmt="%.8e")

# ---- ses with floor --------------------------------------------------------- #
obs = res["DNC_auc"][:, 0].astype(float)
mu_n = null_auc.mean(axis=1); sd_n = null_auc.std(axis=1, ddof=1)
# broadcast trick: per-species obs vs region-1 null of species 1 is not meaningful
# ecologically; this fixture tests only the arithmetic of ses().
scale = float(np.nanstd(obs))
s_out = ses(obs[:5], mu_n, sd_n, scale=scale, rel_floor=0.02)
np.savetxt(f"{OUT}/ses_fixture.csv",
           np.column_stack([obs[:5], mu_n, sd_n, s_out]), delimiter=",",
           header="obs,null_mean,null_sd,ses", comments="", fmt="%.10e")
with open(f"{OUT}/ses_scale.json", "w") as f:
    json.dump({"scale": scale, "rel_floor": 0.02}, f)

# ---- aggregate test + sweep ------------------------------------------------- #
rng2 = np.random.default_rng(7)
n_sp = 400
obs_v = np.clip(rng2.beta(0.4, 2.5, n_sp), 0, 1)
mu_v = np.clip(rng2.beta(0.5, 3.0, n_sp), 0, 1) * 0.8
sd_v = 0.05 + 0.1 * rng2.random(n_sp)
rows = []
for thr_q in [0.80, 0.90, 0.95, 0.99]:
    thr = float(np.quantile(obs_v, thr_q))
    d = aggregate_admissibility_test(obs_v, mu_v, sd_v, thr)
    d["threshold_q"] = thr_q; d["threshold"] = thr
    rows.append(d)
import pandas as pd
pd.DataFrame(rows).to_csv(f"{OUT}/aggregate_sweep.csv", index=False)
np.savetxt(f"{OUT}/agg_inputs.csv", np.column_stack([obs_v, mu_v, sd_v]),
           delimiter=",", header="obs,null_mean,null_sd", comments="", fmt="%.10e")

# ---- breadth-matched null ---------------------------------------------------- #
V = res["DNC_auc"].astype(np.float64)          # (8 species, 5 regions)
# need >=10 species per taxon; replicate block to 30 species with jitter
V30 = np.vstack([V + rng2.normal(0, 1e-3, V.shape) for _ in range(4)])[:30]
V30 = np.clip(V30, 0, 1)
breadth = np.abs(rng2.normal(1, 0.5, 30)) + 0.1
taxon = np.array(["a"] * 30)
SES_b, P_b, MU_b, SD_b, PE_b = breadth_matched_null(V30, breadth, taxon, K=10)
np.savetxt(f"{OUT}/bnull_V.csv", V30, delimiter=",", fmt="%.10e")
np.savetxt(f"{OUT}/bnull_breadth.csv", breadth, delimiter=",", fmt="%.10e")
np.savetxt(f"{OUT}/bnull_ses.csv", SES_b, delimiter=",", fmt="%.8e")
np.savetxt(f"{OUT}/bnull_mu.csv", MU_b, delimiter=",", fmt="%.8e")
np.savetxt(f"{OUT}/bnull_sd.csv", SD_b, delimiter=",", fmt="%.8e")

print("golden fixtures written:", len(os.listdir(OUT)), "files")
print("containment check sp7 x region4: DNC_q=%.4f DNC_auc=%.4f D=%.4f"
      % (res["DNC_q"][6, 3], res["DNC_auc"][6, 3], res["D"][6, 3]))
