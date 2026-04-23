#!/usr/bin/env bash
set -eo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

ENV_NAME="${1:-redcea-benchmark}"
METHODS="${METHODS:-tcrdist3 giana gliph2 tcrnet}"

VDJDB_SLIM="${VDJDB_SLIM:-$PWD/data/vdjdb/vdjdb.slim.txt}"
TRUTH_CSV="${TRUTH_CSV:-$PWD/data/private/01_05_2025_TCRvdb.csv}"
TOOLS_DIR="${TOOLS_DIR:-$PWD/tools}"
WORK_DIR="${WORK_DIR:-$PWD/work/vdjdb_method_benchmark}"
RESULTS_DIR="${RESULTS_DIR:-$PWD/results/vdjdb_method_benchmark}"

GIANA_HOME="${GIANA_HOME:-$TOOLS_DIR/GIANA}"
VDJTOOLS_JAR="${VDJTOOLS_JAR:-$TOOLS_DIR/vdjtools-1.2.1/vdjtools-1.2.1.jar}"

TCRDIST3_RADIUS="${TCRDIST3_RADIUS:-24}"
TCRDIST3_MIN_CLUSTER_SIZE="${TCRDIST3_MIN_CLUSTER_SIZE:-3}"
GIANA_MIN_CLUSTER_SIZE="${GIANA_MIN_CLUSTER_SIZE:-3}"
GLIPH2_MIN_CLUSTER_SIZE="${GLIPH2_MIN_CLUSTER_SIZE:-3}"
TCRNET_MIN_CLUSTER_SIZE="${TCRNET_MIN_CLUSTER_SIZE:-5}"

SHORTS=("GLC" "YLQ")

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "Error: '$1' not found in PATH" >&2
    exit 1
  }
}

require_path() {
  [[ -e "$1" ]] || {
    echo "Error: required path not found: $1" >&2
    exit 1
  }
}

require_cmd conda
require_cmd python
require_cmd Rscript
require_cmd java
require_cmd git

require_path "$VDJDB_SLIM"
require_path "$TRUTH_CSV"
require_path "$GIANA_HOME"
require_path "$VDJTOOLS_JAR"

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

mkdir -p "$WORK_DIR" "$RESULTS_DIR"

python scripts/benchmark_vdjdb_methods.py prepare \
  --vdjdb "$VDJDB_SLIM" \
  --truth "$TRUTH_CSV" \
  --work-dir "$WORK_DIR"

if [[ " $METHODS " == *" tcrdist3 "* ]]; then
  mkdir -p "$RESULTS_DIR/tcrdist3"
  for short_name in "${SHORTS[@]}"; do
    python scripts/benchmark_vdjdb_methods.py run-tcrdist3 \
      --input "$WORK_DIR/inputs/generic/${short_name}.tsv" \
      --output "$RESULTS_DIR/tcrdist3/${short_name}.txt" \
      --short-name "$short_name" \
      --radius "$TCRDIST3_RADIUS" \
      --min-cluster-size "$TCRDIST3_MIN_CLUSTER_SIZE"
  done
  python - <<PY
from pathlib import Path
import pandas as pd
root = Path(r"$RESULTS_DIR") / "tcrdist3"
parts = []
columns = ["gene", "cdr3aa", "v.segm", "j.segm", "cid", "antigen.epitope", "method", "epitope_short"]
for short in ["GLC", "YLQ"]:
    path = root / f"{short}.txt"
    if path.exists():
        try:
            parts.append(pd.read_csv(path, sep="\t"))
        except pd.errors.EmptyDataError:
            pass
out = root / "cluster_members_TRB.txt"
if parts:
    pd.concat(parts, ignore_index=True).to_csv(out, sep="\t", index=False)
else:
    pd.DataFrame(columns=columns).to_csv(out, sep="\t", index=False)
print(f"Saved {out}")
PY
fi

if [[ " $METHODS " == *" giana "* ]]; then
  mkdir -p "$RESULTS_DIR/giana" "$WORK_DIR/giana_raw"
  for short_name in "${SHORTS[@]}"; do
    input_path="$WORK_DIR/inputs/giana/${short_name}.tsv"
    raw_output="$WORK_DIR/giana_raw/${short_name}.txt"
    python "$GIANA_HOME/GIANA4.1.py" -f "$input_path" -O "$raw_output"
    python scripts/benchmark_vdjdb_methods.py convert-giana \
      --input "$WORK_DIR/inputs/generic/${short_name}.tsv" \
      --giana-output "$raw_output" \
      --output "$RESULTS_DIR/giana/${short_name}.txt" \
      --short-name "$short_name" \
      --min-cluster-size "$GIANA_MIN_CLUSTER_SIZE"
  done
  python - <<PY
from pathlib import Path
import pandas as pd
root = Path(r"$RESULTS_DIR") / "giana"
parts = []
columns = ["gene", "cdr3aa", "v.segm", "j.segm", "cid", "antigen.epitope", "method", "epitope_short"]
for short in ["GLC", "YLQ"]:
    path = root / f"{short}.txt"
    if path.exists():
        try:
            parts.append(pd.read_csv(path, sep="\t"))
        except pd.errors.EmptyDataError:
            pass
out = root / "cluster_members_TRB.txt"
if parts:
    pd.concat(parts, ignore_index=True).to_csv(out, sep="\t", index=False)
else:
    pd.DataFrame(columns=columns).to_csv(out, sep="\t", index=False)
print(f"Saved {out}")
PY
fi

if [[ " $METHODS " == *" gliph2 "* ]]; then
  mkdir -p "$RESULTS_DIR/gliph2" "$WORK_DIR/gliph2"
  for short_name in "${SHORTS[@]}"; do
    Rscript scripts/run_gliph2_benchmark.R \
      "$WORK_DIR/inputs/gliph2/${short_name}.tsv" \
      "$RESULTS_DIR/gliph2/${short_name}.txt" \
      "$short_name" \
      "$GLIPH2_MIN_CLUSTER_SIZE" \
      "$WORK_DIR/gliph2/${short_name}"
  done
  python - <<PY
from pathlib import Path
import pandas as pd
root = Path(r"$RESULTS_DIR") / "gliph2"
parts = []
columns = ["gene", "cdr3aa", "v.segm", "j.segm", "cid", "antigen.epitope", "method", "epitope_short"]
for short in ["GLC", "YLQ"]:
    path = root / f"{short}.txt"
    if path.exists():
        try:
            parts.append(pd.read_csv(path, sep="\t"))
        except pd.errors.EmptyDataError:
            pass
out = root / "cluster_members_TRB.txt"
if parts:
    pd.concat(parts, ignore_index=True).to_csv(out, sep="\t", index=False)
else:
    pd.DataFrame(columns=columns).to_csv(out, sep="\t", index=False)
print(f"Saved {out}")
PY
fi

if [[ " $METHODS " == *" tcrnet "* ]]; then
  mkdir -p "$RESULTS_DIR/tcrnet" "$WORK_DIR/tcrnet"
  for short_name in "${SHORTS[@]}"; do
    target="$WORK_DIR/inputs/tcrnet/targets/${short_name}.txt"
    background="$WORK_DIR/inputs/tcrnet/backgrounds/${short_name}.txt"
    out_prefix="$WORK_DIR/tcrnet/${short_name}"
    mkdir -p "$out_prefix"
    java -jar "$VDJTOOLS_JAR" CalcDegreeStats -o 1,0,1 -g2 dummy -b "$background" "$target" "$out_prefix"

    annotated=""
    for candidate in \
      "$out_prefix/$(basename "$target")" \
      "$out_prefix/${short_name}.txt" \
      "$out_prefix/$(basename "$target" .txt).txt"
    do
      if [[ -f "$candidate" ]]; then
        annotated="$candidate"
        break
      fi
    done
    if [[ -z "$annotated" ]]; then
      echo "Error: could not find annotated TCRNet output in $out_prefix" >&2
      exit 1
    fi

    python scripts/benchmark_vdjdb_methods.py convert-tcrnet \
      --input "$WORK_DIR/inputs/generic/${short_name}.tsv" \
      --annotated "$annotated" \
      --output "$RESULTS_DIR/tcrnet/${short_name}.txt" \
      --short-name "$short_name" \
      --min-cluster-size "$TCRNET_MIN_CLUSTER_SIZE"
  done
  python - <<PY
from pathlib import Path
import pandas as pd
root = Path(r"$RESULTS_DIR") / "tcrnet"
parts = []
columns = ["gene", "cdr3aa", "v.segm", "j.segm", "cid", "antigen.epitope", "method", "epitope_short"]
for short in ["GLC", "YLQ"]:
    path = root / f"{short}.txt"
    if path.exists():
        try:
            parts.append(pd.read_csv(path, sep="\t"))
        except pd.errors.EmptyDataError:
            pass
out = root / "cluster_members_TRB.txt"
if parts:
    pd.concat(parts, ignore_index=True).to_csv(out, sep="\t", index=False)
else:
    pd.DataFrame(columns=columns).to_csv(out, sep="\t", index=False)
print(f"Saved {out}")
PY
fi

python scripts/benchmark_vdjdb_methods.py evaluate \
  --truth "$WORK_DIR/truth_glc_ylq.csv" \
  --results-root "$RESULTS_DIR" \
  --output "$RESULTS_DIR/metrics.csv"

echo
echo "Done."
echo "Results: $RESULTS_DIR"
