#!/bin/bash

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-/projects/immunestatus/pogorelyy}"
AIRR_DIR="${AIRR_DIR:-$ROOT_DIR/airr_format}"
OUT_DIR="${OUT_DIR:-$ROOT_DIR/vdjmatch2_yfv}"
MATCHES_DIR="${MATCHES_DIR:-$OUT_DIR/day15_f1_matches}"
LOG_DIR="${LOG_DIR:-$OUT_DIR/logs}"

ENV_NAME="${ENV_NAME:-redcea-benchmark}"
CPUS_PER_TASK="${CPUS_PER_TASK:-16}"
MEMORY="${MEMORY:-16G}"
TIME_LIMIT="${TIME_LIMIT:-04:00:00}"
PARTITION="${PARTITION:-medium}"
CONSTRAINT="${CONSTRAINT:-hpc}"
THREADS="${THREADS:-16}"
CHAIN="${CHAIN:-beta}"
PATTERN="${PATTERN:-*_15_F1.txt}"

mkdir -p "$MATCHES_DIR" "$LOG_DIR"

shopt -s nullglob
samples=( "$AIRR_DIR"/$PATTERN )
shopt -u nullglob

if [[ ${#samples[@]} -eq 0 ]]; then
  echo "No files matched '$PATTERN' in $AIRR_DIR" >&2
  exit 1
fi

for query in "${samples[@]}"; do
  query_base="$(basename "$query")"
  query_id="${query_base%.txt}"

  for target in "${samples[@]}"; do
    target_base="$(basename "$target")"
    target_id="${target_base%.txt}"

    if [[ "$query" == "$target" ]]; then
      continue
    fi

    job_name="vdjm2_${query_id}_vs_${target_id}"
    out_file="$MATCHES_DIR/${query_id}__vs__${target_id}.tsv"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --mem=${MEMORY}
#SBATCH --time=${TIME_LIMIT}
#SBATCH --output=${LOG_DIR}/${job_name}.%j.log
#SBATCH --constraint=${CONSTRAINT}
#SBATCH --partition=${PARTITION}

set -euo pipefail

# shellcheck disable=SC1091
source "\$(conda info --base)/etc/profile.d/conda.sh"
conda activate $(printf '%q' "$ENV_NAME")

mkdir -p $(printf '%q' "$MATCHES_DIR")
vdjmatch2 $(printf '%q' "$query") $(printf '%q' "$target") \\
  --out $(printf '%q' "$out_file") \\
  --max-sub 1 \\
  --max-ins 0 \\
  --max-del 0 \\
  --max-edits 1 \\
  --threads $(printf '%q' "$THREADS") \\
  --junction-col junction_aa \\
  --chain-col locus \\
  --gene $(printf '%q' "$CHAIN")
EOF

    echo "Submitted ${job_name}"
  done
done
