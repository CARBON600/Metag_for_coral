#!/usr/bin/env bash
#
# Run the pipeline. Before starting, check the conda environments declared by
# the rules (referenced by name) and create any that are missing from the
# matching specification under envs/. Existing environments are reused as-is.
#
# Usage:  bash run.sh [-n] [--cores N] [--keep-going] [...]
set -euo pipefail

cd "$(dirname "$0")"

ENVS=(metag megahit binning comebin_env metadecoder SemiBin gtdbtk-2.3.2)

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: 'conda' not found in PATH; activate your conda base/snakemake env first." >&2
    exit 1
fi

have="$(conda env list | awk '{print $1}')"
for e in "${ENVS[@]}"; do
    if echo "$have" | grep -Fxq -- "$e"; then
        echo "[run.sh] found  ${e}"
    else
        if [ ! -f "envs/${e}.yml" ]; then
            echo "ERROR: env '${e}' is missing and envs/${e}.yml does not exist." >&2
            exit 1
        fi
        echo "[run.sh] create ${e}  <- envs/${e}.yml"
        conda env create -n "${e}" -f "envs/${e}.yml"
    fi
done

exec snakemake --use-conda --conda-frontend conda "$@"
