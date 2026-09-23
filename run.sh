#!/usr/bin/env bash
#
# Run the pipeline. Before starting, check the conda environments declared by
# the rules (referenced by name) and create any that are missing from the
# matching specification under envs/. Existing environments are reused as-is.
#
# Usage:  bash run.sh [-n] [--cores N] [--keep-going] [...]
set -euo pipefail

# run.sh itself relies on these before any rule runs: `dirname` resolves the
# script location and `awk` parses config.yaml / conda output. Fail with a clear
# message instead of a bare "command not found".
for _t in awk dirname; do
    if ! command -v "${_t}" >/dev/null 2>&1; then
        echo "ERROR: required tool '${_t}' not found in PATH; run.sh needs it to parse config.yaml/paths." >&2
        exit 1
    fi
done

cd "$(dirname "$0")"

# Export the configured scratch directory so every tool that honours $TMPDIR
# (samtools, CheckM, GTDB-Tk, ...) uses it instead of node-local /tmp. MEGAHIT
# and SPAdes are additionally given --tmp-dir explicitly in their rules.
# Accept both `tmpdir: "/path"` and `tmpdir: /path` (plus a trailing comment) so a
# forgotten pair of quotes cannot silently disable $TMPDIR. A `#` inside the path
# would still be treated as a comment.
tmpdir="$(awk '
    /^tmpdir:/ {
        v = $0
        sub(/^[^:]*:[ \t]*/, "", v)   # drop the "tmpdir:" key
        sub(/[ \t]*#.*$/, "", v)      # drop a trailing comment
        gsub(/"/, "", v)              # drop surrounding quotes
        sub(/^[ \t]+/, "", v); sub(/[ \t]+$/, "", v)
        print v
        exit
    }' config.yaml || true)"
if [ -n "${tmpdir}" ]; then
    export TMPDIR="${tmpdir}"
else
    echo "WARN: could not parse 'tmpdir' from config.yaml (expected: tmpdir: \"/path\" or tmpdir: /path); TMPDIR not set." >&2
fi

ENVS=(metag megahit binning comebin_env metadecoder SemiBin gtdbtk-2.3.2 drep metawrap-env)

# Pin which conda installation is used: the same named env may exist under more
# than one Anaconda install, and a tool (e.g. bowtie2) resolved inside an env can
# otherwise come from a different install. Prefer the conda already on PATH, then
# an explicit $CONDA_BASE. There is deliberately no hardcoded fallback: a
# machine-specific path here would silently select the wrong install (or none).
if [ -z "${CONDA_BASE:-}" ] && command -v conda >/dev/null 2>&1; then
    CONDA_BASE="$(conda info --base 2>/dev/null || true)"
fi
if [ -n "${CONDA_BASE:-}" ] && [ -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]; then
    # shellcheck disable=SC1090
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
fi

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: 'conda' not found in PATH; activate your conda base/snakemake env first, or export CONDA_BASE=/path/to/miniconda3." >&2
    exit 1
fi

# The pipeline runs *under* Snakemake, so the driver environment has to be active
# already. Its pinned specification is envs/snakemake.yml; the expected version is
# read from that file rather than repeated here, so the pin has one home.
DRIVER_SPEC="${PWD}/envs/snakemake.yml"
if ! command -v snakemake >/dev/null 2>&1; then
    echo "ERROR: 'snakemake' not found in PATH. Create and activate the driver environment first:" >&2
    echo "         mamba env create -f ${DRIVER_SPEC}" >&2
    echo "         conda activate snakemake" >&2
    exit 1
fi
snak_pin="$(awk -F= '/^  - snakemake=/{print $2; exit}' "${DRIVER_SPEC}" 2>/dev/null || true)"
snak_run="$(snakemake --version 2>/dev/null || true)"
if [ -n "${snak_pin}" ] && [ -n "${snak_run}" ] && [ "${snak_pin}" != "${snak_run}" ]; then
    echo "WARN: active Snakemake ${snak_run} != ${snak_pin} pinned in envs/snakemake.yml." >&2
elif [ -n "${snak_pin}" ] && [ -n "${snak_run}" ]; then
    echo "[run.sh] driver: snakemake ${snak_run} (= envs/snakemake.yml)"
else
    echo "[run.sh] note: could not read the snakemake pin from ${DRIVER_SPEC}." >&2
fi

cv="$(conda --version 2>/dev/null | awk '{print $2}' || true)"
if [ -n "${cv}" ] && [ "$(printf '%s\n24.7.1\n' "${cv}" | sort -V | head -n1)" != "24.7.1" ]; then
    echo "WARN: conda ${cv} < 24.7.1; Snakemake requires >=24.7.1 when it creates environments." >&2
fi

have="$(conda env list | awk '{print $1}' || true)"
if [ -z "${have}" ]; then
    echo "ERROR: 'conda env list' returned nothing (conda broken, or no environment activated). Refusing to continue: an empty list would make this script try to create all ${#ENVS[@]} environments." >&2
    exit 1
fi
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

# Default to `--rerun-triggers=mtime` so that editing rule code, params, or conda
# environments does NOT silently trigger a full recompute (file mtimes and input
# checksums are still compared as usual). Pass your own --rerun-triggers to
# override.  NOTE: after applying a patch that edits rule code or envs, run once
# with `-F`/`--forceall` (or delete the affected .done markers), otherwise the fix
# will not take effect.  (`--rerun-triggers all` is NOT valid in Snakemake 7.x.)
#
# The value MUST be passed as a single token (`--rerun-triggers=mtime`). Snakemake
# declares this option with nargs='+', so the space-separated form
# `--rerun-triggers mtime <target>` makes argparse greedily swallow every
# following token and abort with "invalid choice: '<target>'".
rerun_args=""
case " $* " in
    *"--rerun-triggers"*) ;;
    *) rerun_args="--rerun-triggers=mtime" ;;
esac

if [ -n "$rerun_args" ]; then
    echo "[run.sh] rerun-triggers: mtime (default; pass --rerun-triggers=mtime to override, or --rerun-triggers=mtime,params,input,software-env,code for the full set)" >&2
else
    echo "[run.sh] rerun-triggers: user-provided (see the snakemake invocation below)" >&2
fi

# By default keep going: without this a single failed job cancels every pending
# job (one run lost 202 jobs to 3 failures). Pass -k/--keep-going yourself to
# keep the message quiet; the flag is 0-arg so it cannot swallow targets.
keep_args="--keep-going"
for _a in "$@"; do
    case "$_a" in
        -k|--keep-going) keep_args="" ;;
    esac
done
if [ -n "$keep_args" ]; then
    echo "[run.sh] --keep-going: enabled by default (a failed job no longer cancels pending jobs)" >&2
fi

# A second launch into the same directory otherwise dies with only a bare
# "Directory cannot be locked". Snakemake keeps lock files in .snakemake/locks.
if [ -d .snakemake/locks ] && [ -n "$(ls -A .snakemake/locks 2>/dev/null || true)" ]; then
    echo "[run.sh] WARNING: .snakemake/locks is not empty." >&2
    echo "[run.sh]          Another run may be active; if not, clear it with: snakemake --unlock" >&2
fi

exec snakemake --use-conda --conda-frontend conda $rerun_args $keep_args "$@"
