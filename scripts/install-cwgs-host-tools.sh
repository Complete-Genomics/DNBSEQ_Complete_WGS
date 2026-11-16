#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 /absolute/install/prefix" >&2
    exit 2
fi

PREFIX=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
CWGS_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
RUNTIME_ENV="$PREFIX/runtime"
SPECIMMUNE_ENV="$PREFIX/specimmune-env"
SPECIMMUNE_HOME="$PREFIX/SpecImmune"
SPECIMMUNE_REF="v0.0.3"
SPECIMMUNE_COMMIT="5262824954bce5bb8578b14440f2a4bd3d172a04"

if command -v mamba >/dev/null 2>&1; then
    CONDA_CMD=mamba
elif command -v conda >/dev/null 2>&1; then
    CONDA_CMD=conda
else
    echo "mamba or conda is required" >&2
    exit 2
fi

mkdir -p "$PREFIX"

"$CONDA_CMD" env create --prefix "$RUNTIME_ENV" \
    --file "$CWGS_ROOT/envs/cwgs-runtime.yaml"

git clone --depth 1 --branch "$SPECIMMUNE_REF" \
    https://github.com/deepomicslab/SpecImmune.git "$SPECIMMUNE_HOME"

if [ "$(git -C "$SPECIMMUNE_HOME" rev-parse HEAD)" != "$SPECIMMUNE_COMMIT" ]; then
    echo "Unexpected SpecImmune revision" >&2
    exit 1
fi

SPECIMMUNE_YAML=$(mktemp)
trap 'rm -f "$SPECIMMUNE_YAML"' EXIT
sed '/^prefix: /d' "$SPECIMMUNE_HOME/environment.yml" > "$SPECIMMUNE_YAML"
"$CONDA_CMD" env create --prefix "$SPECIMMUNE_ENV" \
    --file "$SPECIMMUNE_YAML"
"$SPECIMMUNE_ENV/bin/pip" install --no-deps dysgu==1.6.2
chmod +x "$SPECIMMUNE_HOME"/bin/*
"$SPECIMMUNE_ENV/bin/python" "$SPECIMMUNE_HOME/scripts/main.py" -h >/dev/null

mkdir -p "$PREFIX/nxf-home"
cat <<EOF
Installation complete.
export NXF_HOME="$PREFIX/nxf-home"
"$RUNTIME_ENV/bin/nextflow" run modules/main.nf -profile singularity \\
  --specimmune_env "$SPECIMMUNE_ENV" \\
  --specimmune_home "$SPECIMMUNE_HOME" \\
  --specimmune_db /path/to/specimmune-hla-db
EOF
