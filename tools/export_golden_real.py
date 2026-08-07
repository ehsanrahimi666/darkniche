"""Export a real-data golden set from the authors' pipeline store.

Run this next to your completed pipeline outputs (the directory holding the
niche store .npy/.parquet files written by stage 04 and the metric tables of
stage 05). It dumps, for a small sample of real species and ecoregions, the
flattened densities and the pipeline's own metric values, so the R package
can be verified against *published* numbers, not only synthetic ones.

Usage:
    python export_golden_real.py --store /path/to/04_niches \
        --metrics /path/to/05_metrics/pair_metrics.parquet \
        --out golden_real --n-species 20 --n-regions 10 --seed 1

Then copy the golden_real/ folder into darkniche/inst/extdata/ and run
Rscript tools/verify_real.R (added in v0.2).
"""
import argparse, os, sys
import numpy as np

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", required=True)
    ap.add_argument("--metrics", required=True)
    ap.add_argument("--out", default="golden_real")
    ap.add_argument("--n-species", type=int, default=20)
    ap.add_argument("--n-regions", type=int, default=10)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    # Adjust the import below to wherever the ddarkniche package lives locally:
    from ddarkniche.niche import NicheStore  # noqa
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    store = NicheStore(a.store, R=100)
    meta = store.meta()
    pick = rng.choice(len(meta), size=a.n_species, replace=False)
    Z = store.Z()[pick]
    np.savetxt(f"{a.out}/Zs_real.csv", Z.reshape(len(pick), -1),
               delimiter=",", fmt="%.8e")
    meta.iloc[pick].to_csv(f"{a.out}/species_meta.csv", index=False)
    print("wrote", a.out, "- now add ecoregion densities and the matching "
          "rows of the stage-05 metric table for the same species/regions.")

if __name__ == "__main__":
    main()
