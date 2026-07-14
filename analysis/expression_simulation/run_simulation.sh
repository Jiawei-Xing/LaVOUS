#!/bin/bash
#SBATCH --job-name=sim_array
#SBATCH --partition=cpuq
#SBATCH --time=5:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --array=0-17
#SBATCH --output=sim_%A_%a.out
#SBATCH --error=sim_%A_%a.err
#SBATCH --export=NONE

# Submit:
#   jid=$(sbatch --parsable expression_simulation/run_simulation.sh)
#   sbatch --dependency=afterok:${jid} expression_simulation/run_shuffle.sh

module load EBModules
module load Anaconda3/2024.02-1
source ~/.bashrc
conda activate singlecellstochastics

set -euo pipefail

TREE="lineage_simulation/tree.nwk"
REGIME="lineage_simulation/regime.csv"
OUTDIR="expression_simulation/simulation"
N_GENES=500
TEST_REGIME="1"

mkdir -p "$OUTDIR"

#         label    alpha   sigma   theta0  theta1  r       model
LABELS=(  Base     Thetas  Theta1  Sigma   BM      AlphaS  AlphaL  rS      rL      r0      t3r50   t3r0    t5r5    t5r50   t5r0    t7r5    t7r50   t7r0  )
ALPHAS=(  1        1       1       1       none    0.3     3       1       1       1       1       1       1       1       1       1       1       1     )
SIGMAS=(  3        3       3       5       3       3       3       3       3       3       3       3       3       3       3       3       3       3     )
THETA0S=( 1        3       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1     )
THETA1S=( 1        3       3       1       none    1       1       1       1       1       3       3       5       5       5       7       7       7     )
RS=(      5        5       5       5       5       5       5       0.5     50      0       50      0       5       50      0       5       50      0     )
MODELS=(  OU       OU      OU      OU      BM      OU      OU      OU      OU      OU      OU      OU      OU      OU      OU      OU      OU      OU    )

i=${SLURM_ARRAY_TASK_ID}
label=${LABELS[$i]}
alpha=${ALPHAS[$i]}
sigma=${SIGMAS[$i]}
theta0=${THETA0S[$i]}
theta1=${THETA1S[$i]}
r=${RS[$i]}
model=${MODELS[$i]}

OUTFILE="${OUTDIR}/readcounts_${label}.tsv"
if [ -f "$OUTFILE" ]; then
    echo "Skipping ${label}: ${OUTFILE} already exists"
    exit 0
fi

echo "=== Scenario ${label} (task ${i}) ==="

cmd="lavous-simulate --tree ${TREE} --regime ${REGIME} --background ${theta0} --n_genes ${N_GENES} --sigma ${sigma} --out ${OUTDIR} --label ${label} --tree_plot"
[ "$alpha"  != "none" ] && cmd+=" --alpha ${alpha}"
[ "$theta1" != "none" ] && cmd+=" --optim ${theta1}"
[ "$r"      != "none" ] && cmd+=" --dispersion ${r}"
[ "$model"  != "OU"   ] && cmd+=" --bg ${model}"
# pure-BM has no OU test regime; otherwise apply test regime
[ "$model"  == "OU"   ] && cmd+=" --test ${TEST_REGIME}"

echo "$cmd"
eval "$cmd"
