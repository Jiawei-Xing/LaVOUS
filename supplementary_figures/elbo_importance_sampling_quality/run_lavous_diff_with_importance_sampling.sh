#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "${SCRIPT_DIR}/../.." && pwd)
LAVOUS_DIFF=${LAVOUS_DIFF:-/grid/siepel/home/xing/.conda/envs/singlecellstochastics/bin/lavous-diff}

LABEL=${LABEL:-t5r5}
OUTDIR=${OUTDIR:-${SCRIPT_DIR}/lavous_diff_${LABEL}_defensive_is}
PREFIX=${PREFIX:-diff_${LABEL}_defensive_is}
IMPORTANCE=${IMPORTANCE:-512}
MIX=${MIX:-1.0}
BATCH=${BATCH:-50}
ITER=${ITER:-10000}
LR=${LR:-0.1}
WINDOW=${WINDOW:-200}
TOL=${TOL:-0.0001}
DTYPE=${DTYPE:-float32}
NULL_REGIME=${NULL_REGIME:-0}

TREE=${TREE:-${REPO_ROOT}/lineage_simulation/tree.nwk}
EXPRESSION=${EXPRESSION:-${REPO_ROOT}/expression_simulation/simulation/readcounts_${LABEL}.tsv}
REGIME=${REGIME:-${REPO_ROOT}/lineage_simulation/regime.csv}

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

"${LAVOUS_DIFF}" \
  --tree "${TREE}" \
  --expression "${EXPRESSION}" \
  --regime "${REGIME}" \
  --null "${NULL_REGIME}" \
  --outdir "${OUTDIR}" \
  --prefix "${PREFIX}" \
  --batch "${BATCH}" \
  --iter "${ITER}" \
  --lr "${LR}" \
  --window "${WINDOW}" \
  --tol "${TOL}" \
  --dtype "${DTYPE}" \
  --importance "${IMPORTANCE}" \
  --mix "${MIX}"
