#!/bin/bash
#SBATCH --job-name=fig2_perm_reconst
#SBATCH --partition=cpuq
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=log/reconstruct_%j.out
#SBATCH --error=log/reconstruct_%j.err
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
OUTDIR="reconst"
DIFFDIR="diff"

mkdir -p "$OUTDIR" log

LABELS="Base rL r0 Theta1 t3r50 t3r0 t5r5 t5r50 t5r0 t7r5 t7r50 t7r0"

for label in $LABELS; do
    OUT_TSV="${OUTDIR}/history_${label}_gene500.tsv"
    OUT_FIG="${OUTDIR}/history_${label}_gene500.png"
    Q_PARAMS="${DIFFDIR}/diff_${label}_h1_q-mean-std_0.tsv"
    MODEL_PARAMS="${DIFFDIR}/diff_${label}_model-params.tsv"
    READ_COUNTS="${SIMDIR}/readcounts_${label}.tsv"

    if [ ! -f "$Q_PARAMS" ] || [ ! -f "$MODEL_PARAMS" ] || [ ! -f "$READ_COUNTS" ]; then
        echo "Skipping ${label}: missing diff or simulation inputs"
        continue
    fi
    if [ -f "$OUT_TSV" ] && [ -f "$OUT_FIG" ]; then
        echo "Skipping ${label}: reconstruction already exists"
        continue
    fi

    echo "=== Reconstruct ${label} ==="
    lavous-reconstruct \
        --tree "$TREE" \
        --q_params "$Q_PARAMS" \
        --read_counts "$READ_COUNTS" \
        --gene Gene_500 \
        --model ou \
        --regime "$REGIME" \
        --ou "$MODEL_PARAMS" \
        --out_tsv "$OUT_TSV" \
        --out_fig "$OUT_FIG"
done

echo "=== All reconstruction runs complete ==="

