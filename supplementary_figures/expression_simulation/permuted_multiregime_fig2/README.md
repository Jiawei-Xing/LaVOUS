# Fig. 2 rerun on topology-perturbed multi-regime tree

This directory reruns the Fig. 2 simulation workflow using the topology-perturbed
multi-origin regime tree generated in `supplementary_figures/lineage_simulation/permuted_multiregime_tree`.
All outputs are written under this directory and do not overwrite the original
`expression_simulation/{simulation,heritability,diff,reconst}` results.

Submit the LaVOUS pipeline:

```bash
cd supplementary_figures/expression_simulation/permuted_multiregime_fig2
bash submit_fig2_pipeline.sh
```

Optional EvoGeneX benchmark for the ROC panel:

```bash
cd supplementary_figures/expression_simulation/permuted_multiregime_fig2
RUN_EVOGENEX=1 bash submit_fig2_pipeline.sh
```

Main outputs:

- `simulation/readcounts_*.tsv`
- `simulation/sim_history_*.tsv`
- `heritability/heritability_*.tsv`
- `diff/diff_*`
- `reconst/history_*_gene500.*`
- `heritability/heritability_LR_comparison.{png,pdf}`
- `reconst/comparison.png`
- `diff_ROC.{png,pdf}`

