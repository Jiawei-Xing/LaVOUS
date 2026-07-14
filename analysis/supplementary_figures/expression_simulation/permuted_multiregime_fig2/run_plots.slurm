#!/bin/bash
#SBATCH --job-name=fig2_perm_plots
#SBATCH --partition=cpuq
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=log/plots_%j.out
#SBATCH --error=log/plots_%j.err
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
SCRIPT_DIR="$(pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

mkdir -p log

echo "=== Plot heritability summaries ==="
(
    cd heritability
    python "$REPO_ROOT/expression_simulation/heritability/plot_plasticity_results.py"
)

echo "=== Plot reconstruction comparison ==="
python "$REPO_ROOT/expression_simulation/plot_sim-vs-reconst.py" \
    --simdir "$SCRIPT_DIR/simulation" \
    --reconstdir "$SCRIPT_DIR/reconst" \
    --tree "$SCRIPT_DIR/lineage_simulation/tree.nwk"

echo "=== Plot differential ROC ==="
python plot_diff.py

echo "=== Plotting complete ==="

