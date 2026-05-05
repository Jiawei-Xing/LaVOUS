#!/bin/bash
#SBATCH --job-name=sim_result_pipeline
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:h100:1
#SBATCH --time=10:00:00
#SBATCH --mem=40G
#SBATCH --output=sim_diff_%j.out
#SBATCH --error=sim_diff_%j.err
#SBATCH --export=NONE

module load EBModules
module load Anaconda3/2024.02-1
source ~/.bashrc
conda activate singlecellstochastics

set -euo pipefail

TREE="lineage_simulation/tree.nwk"
REGIME="lineage_simulation/regime.csv"
SIMDIR="expression_simulation/simulation"
OUTDIR="expression_simulation/diff"

mkdir -p "$OUTDIR"

# theta1 in {1,3,5,7} x r in {5,50,0}
LABELS="Base rL r0 Theta1 t3r50 t3r0 t5r5 t5r50 t5r0 t7r5 t7r50 t7r0"

for label in $LABELS; do
    EXPR="${SIMDIR}/readcounts_${label}.tsv"
    OUTFILE="${OUTDIR}/diff_${label}_chi-squared.tsv"
    if [ ! -f "$EXPR" ]; then
        echo "Skipping ${label}: ${EXPR} not found"
        continue
    fi
    if [ -f "$OUTFILE" ]; then
        echo "Skipping ${label}: ${OUTFILE} already exists"
        continue
    fi

    echo "=== Diff ${label} ==="
    lavous-diff \
        --tree "$TREE" \
        --expression "$EXPR" \
        --regime "$REGIME" \
        --null 0 \
        --outdir "$OUTDIR" \
        --prefix "diff_${label}" \
        --batch 50 \
        --wandb sim_diff \
        --resume
done

echo "=== All diff runs complete ==="
