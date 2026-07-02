#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p log simulation heritability diff reconst

sim_jid=$(sbatch --parsable run_simulation.slurm)
shuffle_jid=$(sbatch --parsable --dependency=afterok:${sim_jid} run_shuffle.slurm)
herit_jid=$(sbatch --parsable --dependency=afterok:${shuffle_jid} run_heritability.slurm)
diff_jid=$(sbatch --parsable --dependency=afterok:${sim_jid} run_diff.slurm)
reconst_jid=$(sbatch --parsable --dependency=afterok:${diff_jid} run_reconstruct.slurm)

plot_deps="${herit_jid}:${diff_jid}:${reconst_jid}"
egx_jid=""
if [ "${RUN_EVOGENEX:-0}" = "1" ]; then
    egx_jid=$(sbatch --parsable --dependency=afterok:${sim_jid} run_evogenex.slurm)
    plot_deps="${plot_deps}:${egx_jid}"
fi

plot_jid=$(sbatch --parsable --dependency=afterok:${plot_deps} run_plots.slurm)

echo "Submitted Fig. 2 rerun jobs:"
echo "  simulation:     ${sim_jid}"
echo "  shuffle:        ${shuffle_jid}"
echo "  heritability:   ${herit_jid}"
echo "  diff:           ${diff_jid}"
echo "  reconstruction: ${reconst_jid}"
if [ -n "$egx_jid" ]; then
    echo "  EvoGeneX:       ${egx_jid}"
fi
echo "  plots:          ${plot_jid}"

