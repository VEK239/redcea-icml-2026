#!/usr/bin/env bash
set -eo pipefail

ENV_NAME="${1:-redcea-benchmark}"
PYTHON_VERSION="${PYTHON_VERSION:-3.11}"
TOOLS_DIR="${TOOLS_DIR:-$PWD/tools}"

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "Error: '$1' not found in PATH" >&2
    exit 1
  }
}

require_cmd conda
require_cmd git
require_cmd curl
require_cmd unzip

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"

mkdir -p "$TOOLS_DIR"

if ! conda env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
  conda create -y -n "$ENV_NAME" "python=${PYTHON_VERSION}" pip openjdk
else
  conda install -y -n "$ENV_NAME" "python=${PYTHON_VERSION}" pip openjdk
fi

conda activate "$ENV_NAME"

python -m pip install --upgrade pip setuptools wheel
python -m pip install \
  numpy \
  pandas \
  scipy \
  scikit-learn \
  biopython \
  numba \
  cython \
  faiss-cpu \
  tcrdist3 \
  "git+https://github.com/antigenomics/redcea"

if [[ ! -d "$TOOLS_DIR/GIANA" ]]; then
  git clone https://github.com/s175573/GIANA.git "$TOOLS_DIR/GIANA"
fi

if [[ ! -f "$TOOLS_DIR/vdjtools-1.2.1/vdjtools-1.2.1.jar" ]]; then
  curl -L https://github.com/mikessh/vdjtools/releases/download/1.2.1/vdjtools-1.2.1.zip -o "$TOOLS_DIR/vdjtools-1.2.1.zip"
  unzip -o "$TOOLS_DIR/vdjtools-1.2.1.zip" -d "$TOOLS_DIR" >/dev/null
fi

cat <<EOF

Done.

Activate environment:
  conda activate $ENV_NAME

Installed:
  redcea                  -> Python package
  tcrdist3                -> Python package
  GIANA                   -> $TOOLS_DIR/GIANA
  TCRNet / VDJtools       -> $TOOLS_DIR/vdjtools-1.2.1/vdjtools-1.2.1.jar

Examples:
  python -c "import redcea; print('redcea ok')"
  python -c "import tcrdist; print('tcrdist3 ok')"
  python "$TOOLS_DIR/GIANA/GIANA4.1.py" -h
  java -jar "$TOOLS_DIR/vdjtools-1.2.1/vdjtools-1.2.1.jar" CalcDegreeStats -h
EOF
