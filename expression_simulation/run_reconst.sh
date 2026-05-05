set -euo pipefail

TREE="lineage_simulation/tree.nwk"
REGIME="lineage_simulation/regime.csv"
SIMDIR="expression_simulation/simulation"
OUTDIR="expression_simulation/reconst"
DIFFDIR="expression_simulation/diff"

mkdir -p "$OUTDIR"

# theta1 in {1,3,5,7} x r in {5,50,0}
LABELS="Base rL r0 Theta1 t3r50 t3r0 t5r5 t5r50 t5r0 t7r5 t7r50 t7r0"

for label in $LABELS; do
    echo "=== Diff ${label} ==="
    lavous-reconstruct \
        --tree "$TREE" \
        --q_params "${DIFFDIR}/diff_${label}_h1_q-mean-std_0.tsv" \
        --read_counts "${SIMDIR}/readcounts_${label}.tsv" \
        --gene Gene_500 \
        --model ou \
        --regime "$REGIME" \
        --ou "${DIFFDIR}/diff_${label}_model-params.tsv" \
        --out_tsv "${OUTDIR}/history_${label}_gene500.tsv" \
        --out_fig "${OUTDIR}/history_${label}_gene500.png"
done

echo "=== All diff runs complete ==="
