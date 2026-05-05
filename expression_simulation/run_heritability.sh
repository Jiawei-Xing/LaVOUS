#!/bin/bash
#SBATCH --job-name=sim_result_pipeline
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:h100:1
#SBATCH --time=5:00:00
#SBATCH --mem=40G
#SBATCH --output=sim_heritability_%j.out
#SBATCH --error=sim_heritability_%j.err
#SBATCH --export=NONE

module load EBModules
module load Anaconda3/2024.02-1
source ~/.bashrc
conda activate singlecellstochastics

set -euo pipefail

TREE="lineage_simulation/tree.nwk"
SIMDIR="expression_simulation/simulation"
OUTDIR="expression_simulation/heritability"

mkdir -p "$OUTDIR"

LABELS="Base Thetas Theta1 Sigma BM AlphaS AlphaL rS rL r0"

for label in $LABELS; do
    EXPR="${SIMDIR}/readcounts_${label}.tsv"
    OUTFILE="${OUTDIR}/heritability_${label}.tsv"
    if [ ! -f "$EXPR" ]; then
        echo "Skipping ${label}: ${EXPR} not found"
        continue
    fi
    if [ -f "$OUTFILE" ]; then
        echo "Skipping ${label}: ${OUTFILE} already exists"
        continue
    fi

    echo "=== Heritability ${label} ==="
    lavous-heritability \
        --tree "$TREE" \
        --expression "$EXPR" \
        --outfile "$OUTFILE"
done

for label in $LABELS; do
    EXPR="${SIMDIR}/readcounts_${label}_shuffled.tsv"
    OUTFILE="${OUTDIR}/heritability_${label}_shuffled.tsv"
    if [ ! -f "$EXPR" ]; then
        echo "Skipping ${label}_shuffled: ${EXPR} not found"
        continue
    fi
    if [ -f "$OUTFILE" ]; then
        echo "Skipping ${label}_shuffled: ${OUTFILE} already exists"
        continue
    fi

    echo "=== Heritability ${label}_shuffled ==="
    lavous-heritability \
        --tree "$TREE" \
        --expression "$EXPR" \
        --outfile "$OUTFILE"
done

echo "=== All heritability runs complete ==="
