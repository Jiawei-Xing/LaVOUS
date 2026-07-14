#!/bin/bash
#SBATCH --job-name=fig2_perm_egx
#SBATCH --time=20:00:00
#SBATCH --output=log/egx_%A_%a.out
#SBATCH --error=log/egx_%A_%a.err
#SBATCH --export=NONE
#SBATCH --array=0-11
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

module load EBModules
module load Anaconda3/2024.02-1
source ~/.bashrc
conda activate EGX

set -euo pipefail

if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    cd "$SLURM_SUBMIT_DIR"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    cd "$SCRIPT_DIR"
fi
SCRIPT_DIR="$(pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

TREE="lineage_simulation/tree.nwk"
REGIME="lineage_simulation/regime.csv"
TREE_EGX="lineage_simulation/tree_egx.nwk"
REGIME_EGX="lineage_simulation/regime_mrca.csv"
SIMDIR="simulation"
OUTDIR="evogenex"
PYTHON_PREP="${PYTHON_PREP:-/grid/siepel/home/xing/.conda/envs/singlecellstochastics/bin/python}"

mkdir -p "$OUTDIR" log

"$PYTHON_PREP" - "$TREE" "$REGIME" "$TREE_EGX" "$REGIME_EGX" <<'PY'
import csv
import sys
import pandas as pd
from Bio import Phylo

tree_path, regime_path, tree_out, regime_out = sys.argv[1:]
tree = Phylo.read(tree_path, "newick")
regimes = pd.read_csv(regime_path).astype(str)
regime_map = dict(zip(regimes["node_name"], regimes["regime"]))

# EvoGeneX inputs used by the original workflow omit the artificial unary c0 root.
if len(tree.root.clades) == 1:
    tree.root = tree.root.clades[0]
Phylo.write(tree, tree_out, "newick")

with open(regime_out, "w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["node", "node2", "regime"])
    for clade in tree.find_clades(order="preorder"):
        regime = regime_map[str(clade.name)]
        if clade.is_terminal():
            writer.writerow([clade.name, "", regime])
        else:
            left = clade.clades[0].get_terminals()[0].name
            right = clade.clades[1].get_terminals()[0].name
            writer.writerow([left, right, regime])
PY

LABELS=(Base rL r0 Theta1 t3r50 t3r0 t5r5 t5r50 t5r0 t7r5 t7r50 t7r0)
label=${LABELS[$SLURM_ARRAY_TASK_ID]}

counts_file="${SIMDIR}/readcounts_${label}.tsv"
long="${OUTDIR}/readcounts_${label}_long.csv"
regime0="${OUTDIR}/regime0_${label}.csv"
output_csv="${OUTDIR}/egx_${label}.csv"

if [ ! -f "$long" ]; then
    "$PYTHON_PREP" - "$counts_file" "$long" <<'PY'
import sys
import pandas as pd

in_file, out_file = sys.argv[1:]
df = pd.read_csv(in_file, sep="\t", index_col=0)
with open(out_file, "w") as handle:
    handle.write("gene,species,replicate,exprval\n")
    for gene_idx, gene in enumerate(df.columns, start=1):
        for cell in df.index:
            handle.write(f"{gene_idx},{cell},R1,{df.loc[cell, gene]}\n")
PY
fi

if [ ! -f "$regime0" ]; then
    "$PYTHON_PREP" - "$REGIME_EGX" "$regime0" <<'PY'
import sys
import pandas as pd

in_file, out_file = sys.argv[1:]
df = pd.read_csv(in_file)
df.iloc[:, -1] = 0
df.to_csv(out_file, index=False)
PY
fi

if [ -f "$output_csv" ]; then
    echo "Skipping ${label}: ${output_csv} already exists"
    exit 0
fi

echo "=== EvoGeneX ${label} ==="
Rscript "$REPO_ROOT/expression_simulation/evogenex/adaptive_evogenex.R" \
    "$TREE_EGX" \
    "$regime0" \
    "$REGIME_EGX" \
    "$long" \
    "$output_csv"

