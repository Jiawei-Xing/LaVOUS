#!/bin/bash
#SBATCH --job-name=sim_shuffle
#SBATCH --partition=cpuq
#SBATCH --time=0:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=sim_shuffle_%j.out
#SBATCH --error=sim_shuffle_%j.err
#SBATCH --export=NONE

# Submit after run_simulation.sh:
#   jid=$(sbatch --parsable expression_simulation/run_simulation.sh)
#   sbatch --dependency=afterok:${jid} expression_simulation/run_shuffle.sh

module load EBModules
module load Anaconda3/2024.02-1
source ~/.bashrc
conda activate singlecellstochastics

set -euo pipefail

OUTDIR="expression_simulation/simulation"

OUTDIR="$OUTDIR" python <<'PY'
import os
import pandas as pd

outdir = os.environ["OUTDIR"]
simulations = ["Base", "Thetas", "Theta1", "Sigma", "BM", "AlphaS", "AlphaL", "rS", "rL", "r0"]
for n in simulations:
    df = pd.read_csv(f"{outdir}/readcounts_{n}.tsv", sep="\t", index_col=0)
    shuffled = df.apply(lambda col: col.sample(frac=1).values)
    shuffled.to_csv(f"{outdir}/readcounts_{n}_shuffled.tsv", sep="\t", header=True, index=True)
    print(f"shuffled {n}")
PY

echo "=== Shuffle complete ==="
