#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "${SCRIPT_DIR}/../.." && pwd)
PYTHON_BIN=${PYTHON_BIN:-/grid/siepel/home/xing/.conda/envs/singlecellstochastics/bin/python}

LABEL=${LABEL:-t5r5}
N_GENES=${N_GENES:-40}
N_SAMPLES=${N_SAMPLES:-512}
N_REPLICATES=${N_REPLICATES:-3}
BATCH_SIZE=${BATCH_SIZE:-4}
PROPOSAL=${PROPOSAL:-overdispersed}
MIX=${MIX:-0.8}
Q_SCALE=${Q_SCALE:-2.0}
SEED=${SEED:-20240625}

cd "${REPO_ROOT}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/plot_elbo_importance_sampling_quality.py" \
  --label "${LABEL}" \
  --n-genes "${N_GENES}" \
  --n-samples "${N_SAMPLES}" \
  --n-replicates "${N_REPLICATES}" \
  --batch-size "${BATCH_SIZE}" \
  --proposal "${PROPOSAL}" \
  --mix "${MIX}" \
  --q-scale "${Q_SCALE}" \
  --seed "${SEED}" \
  --force
