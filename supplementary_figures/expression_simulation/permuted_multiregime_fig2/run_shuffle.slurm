#!/bin/bash
#SBATCH --job-name=fig2_perm_shuffle
#SBATCH --partition=cpuq
#SBATCH --time=0:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=log/shuffle_%j.out
#SBATCH --error=log/shuffle_%j.err
#SBATCH --export=NONE

module load EBModules
module load Anaconda3/2024.02-1
source ~/.bashrc
conda activate singlecellstochastics

set -euo pipefail

if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    cd "$SLURM_SUBMIT_DIR"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    cd "$SCRIPT_DIR"
fi

OUTDIR="simulation"

OUTDIR="$OUTDIR" python <<'PY'
import os
import pandas as pd

outdir = os.environ["OUTDIR"]
simulations = ["Base", "Thetas", "Theta1", "Sigma", "BM", "AlphaS", "AlphaL", "rS", "rL", "r0"]
for name in simulations:
    in_path = f"{outdir}/readcounts_{name}.tsv"
    out_path = f"{outdir}/readcounts_{name}_shuffled.tsv"
    if os.path.exists(out_path):
        print(f"skipping {name}: {out_path} already exists")
        continue
    df = pd.read_csv(in_path, sep="\t", index_col=0)
    shuffled = df.apply(lambda col: col.sample(frac=1).values)
    shuffled.to_csv(out_path, sep="\t", header=True, index=True)
    print(f"shuffled {name}")
PY

echo "=== Shuffle complete ==="

