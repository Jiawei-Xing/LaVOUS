#!/bin/bash
#SBATCH --job-name=egx_array
#SBATCH --time=20:00:00
#SBATCH --output=egx_%A_%a.out  # %A is the array job ID, %a is the specific task ID
#SBATCH --error=egx_%A_%a.err
#SBATCH --export=NONE
#SBATCH --array=0-11             # Creates 12 tasks numbered 0 through 11
#SBATCH --cpus-per-task=4        # Request 1 core per task (adjust if your python/bash scripts use more)
#SBATCH --mem=16G                 # Request 8GB of memory per task (adjust as needed)

module load EBModules
module load Anaconda3/2024.02-1
source ~/.bashrc
conda activate EGX

set -euo pipefail

TREE="lineage_simulation/tree_egx.nwk"
REGIME="lineage_simulation/regime_mrca.csv"
SIMDIR="expression_simulation/simulation"
OUTDIR="expression_simulation/evogenex"

# Ensure output directory exists
mkdir -p "$OUTDIR"

# 1. Define all your labels in a bash array
LABELS=(Base rL r0 Theta1 t3r50 t3r0 t5r5 t5r50 t5r0 t7r5 t7r50 t7r0)

# 2. Extract the specific label for this Slurm task using the array task ID
label=${LABELS[$SLURM_ARRAY_TASK_ID]}

counts_file="$SIMDIR/readcounts_${label}.tsv"
long="$OUTDIR/readcounts_${label}_long.csv"
REGIME0="$OUTDIR/regime0_${label}.csv"
output_csv="$OUTDIR/egx_${label}.csv"

# python - "$counts_file" "$long" <<'PY'
# import sys
# import pandas as pd

# in_file = sys.argv[1]
# out_file = sys.argv[2]

# # count matrix: rows = cells, columns = genes
# df = pd.read_csv(in_file, sep="\t", index_col=0)

# with open(out_file, "w") as f:
#     f.write("gene,species,replicate,exprval\n")
#     for gene_idx, gene in enumerate(df.columns, start=1):
#         for cell in df.index:
#             f.write(f"{gene_idx},{cell},R1,{df.loc[cell, gene]}\n")
# PY

# python - "$REGIME" "$REGIME0" <<'PY'
# import sys
# import pandas as pd

# in_file = sys.argv[1]
# out_file = sys.argv[2]

# df = pd.read_csv(in_file, header=0)

# # Set the last column to all 0
# df.iloc[:, -1] = 0

# df.to_csv(out_file, index=False)
# PY

# 3. Run the pipeline just for this specific label
echo "Task ID $SLURM_ARRAY_TASK_ID is processing label: $label"

Rscript "$OUTDIR/adaptive_evogenex.R" \
    "$TREE" \
    "$REGIME0" \
    "$REGIME" \
    "$long" \
    "$output_csv"
