#!/bin/bash
#SBATCH --job-name=fig2_perm_diff
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:h100:1
#SBATCH --time=10:00:00
#SBATCH --mem=40G
#SBATCH --output=log/diff_%j.out
#SBATCH --error=log/diff_%j.err
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

TREE="lineage_simulation/tree.nwk"
REGIME="lineage_simulation/regime.csv"
SIMDIR="simulation"
OUTDIR="diff"

mkdir -p "$OUTDIR" log

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
        --wandb fig2_permuted_multiregime_diff \
        --resume
done

echo "=== All diff runs complete ==="

