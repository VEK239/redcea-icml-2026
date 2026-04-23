#!/bin/bash

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-/projects/immunestatus/pogorelyy}"
AIRR_DIR="${AIRR_DIR:-$ROOT_DIR/airr_format}"
OUT_DIR="${OUT_DIR:-$ROOT_DIR/redcea}"
RUNS_DIR="${RUNS_DIR:-$OUT_DIR/runs}"
LOG_DIR="${LOG_DIR:-$OUT_DIR/logs/redcea}"

CPUS_PER_TASK="${CPUS_PER_TASK:-16}"
MEMORY="${MEMORY:-64gb}"
TIME_LIMIT="${TIME_LIMIT:-08:00:00}"
PARTITION="${PARTITION:-medium}"
CONSTRAINT="${CONSTRAINT:-hpc}"
NPROC="${NPROC:-16}"
CHAIN="${CHAIN:-TRB}"
DONORS=(${DONORS:-P1 P2 Q1 Q2 S1 S2})
SAMPLE_REPLICA="${SAMPLE_REPLICA:-F1}"
BACKGROUND_REPLICA="${BACKGROUND_REPLICA:-F1}"
K_NEIGHBORS="${K_NEIGHBORS:-3}"
CLUSTER_MIN_SAMPLES="${CLUSTER_MIN_SAMPLES:-3}"
LEIDEN_RESOLUTION="${LEIDEN_RESOLUTION:-2.0}"
RESOLUTION_TAG="${RESOLUTION_TAG:-${LEIDEN_RESOLUTION//./p}}"
PARAM_TAG="${PARAM_TAG:-k${K_NEIGHBORS}_res${RESOLUTION_TAG}_min${CLUSTER_MIN_SAMPLES}}"

mkdir -p "$OUT_DIR" "$RUNS_DIR" "$LOG_DIR"

for donor in "${DONORS[@]}"; do
  base_prefix="yfv_${donor}_${SAMPLE_REPLICA}"
  prefix="${base_prefix}_leiden_${PARAM_TAG}"
  sample_path="$AIRR_DIR/${donor}_15_${SAMPLE_REPLICA}.tsv"
  background_path="$AIRR_DIR/${donor}_0_${BACKGROUND_REPLICA}.tsv"
  embedding_run_dir="$RUNS_DIR/$base_prefix"
  run_dir="$RUNS_DIR/${base_prefix}_leiden_${PARAM_TAG}"
  sample_embedding_path="$embedding_run_dir/${base_prefix}_sample_embeddings.parquet"
  background_embedding_path="$embedding_run_dir/${base_prefix}_background_embeddings.parquet"

  if [[ ! -f "$sample_path" ]]; then
    sample_path="$AIRR_DIR/${donor}_15_${SAMPLE_REPLICA}.txt"
  fi
  if [[ ! -f "$background_path" ]]; then
    background_path="$AIRR_DIR/${donor}_0_${BACKGROUND_REPLICA}.txt"
  fi

  if [[ ! -f "$sample_path" ]]; then
    echo "Skipping ${donor}: sample AIRR file not found" >&2
    continue
  fi
  if [[ ! -f "$background_path" ]]; then
    echo "Skipping ${donor}: background AIRR file not found" >&2
    continue
  fi
  if [[ ! -f "$sample_embedding_path" ]]; then
    echo "Skipping ${donor}: sample embeddings not found at ${sample_embedding_path}" >&2
    continue
  fi
  if [[ ! -f "$background_embedding_path" ]]; then
    echo "Skipping ${donor}: background embeddings not found at ${background_embedding_path}" >&2
    continue
  fi

  sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${prefix}
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --mem=${MEMORY}
#SBATCH --time=${TIME_LIMIT}
#SBATCH --output=${LOG_DIR}/${prefix}.%j.log
#SBATCH --constraint=${CONSTRAINT}
#SBATCH --partition=${PARTITION}

set -euo pipefail
mkdir -p $(printf '%q' "$run_dir") $(printf '%q' "$OUT_DIR")

redcea \\
  -is $(printf '%q' "$sample_path") \\
  -ib $(printf '%q' "$background_path") \\
  -o $(printf '%q' "$run_dir") \\
  -e $(printf '%q' "$prefix") \\
  -c $(printf '%q' "$CHAIN") \\
  -np $(printf '%q' "$NPROC") \\
  -se $(printf '%q' "$sample_embedding_path") \\
  -be $(printf '%q' "$background_embedding_path") \\
  --cluster-algo leiden \\
  -kn $(printf '%q' "$K_NEIGHBORS") \\
  -ms $(printf '%q' "$CLUSTER_MIN_SAMPLES") \\
  --leiden-resolution $(printf '%q' "$LEIDEN_RESOLUTION")

cp -f $(printf '%q' "$run_dir/${prefix}_tcremp_clusters.tsv") $(printf '%q' "$OUT_DIR/")
cp -f $(printf '%q' "$run_dir/${prefix}_summary_tcrempnet.tsv") $(printf '%q' "$OUT_DIR/")
cp -f $(printf '%q' "$run_dir/${prefix}_enriched_clonotypes_tcremp.tsv") $(printf '%q' "$OUT_DIR/")
EOF
done
