#!/usr/bin/env bash
# One-command regression run for Metag_for_coral-main on a tiny fixture.
#
#   bash tools/mini/validate_run.sh --repo ~/pipeline/Metag_for_coral-main [--cores 16]
#
# It syncs the repo into an isolated clone (default <repo>.minitest) so the real
# output/ tree is never touched, dry-runs to fingerprint the DAG, runs the whole
# pipeline at forced recompute, then asserts the outputs.
#
# Options:
#   --repo DIR       the checkout under test (required)
#   --clone DIR      isolated working copy (default <repo>.minitest)
#   --cores N        default 16
#   --incremental    do not force recompute (fast; only re-runs what changed)
#   --dry-only       stop after the DAG fingerprint
#   --bless          overwrite the expected DAG fingerprint
#   --no-sync        use the clone as-is (skip rsync from --repo)
set -uo pipefail

TOOLS="$(cd "$(dirname "$0")" && pwd)"
REPO=""; CLONE=""; CORES=16; FULL=1; DRY_ONLY=0; BLESS=0; SYNC=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo)        REPO="$2"; shift 2 ;;
    --clone)       CLONE="$2"; shift 2 ;;
    --cores)       CORES="$2"; shift 2 ;;
    --incremental) FULL=0; shift ;;
    --dry-only)    DRY_ONLY=1; shift ;;
    --bless)       BLESS=1; shift ;;
    --no-sync)     SYNC=0; shift ;;
    -h|--help)     sed -n '2,16p' "$0"; exit 0 ;;
    *)             echo "unknown option: $1" >&2; exit 2 ;;
  esac
done
[[ -n "${REPO}" ]] || { echo "ERROR: --repo required (see --help)" >&2; exit 2; }
[[ -d "${REPO}" ]] || { echo "ERROR: no such repo: ${REPO}" >&2; exit 2; }
REPO="$(readlink -f "${REPO}")"
CLONE="${CLONE:-${REPO%/}.minitest}"
CFG_REL="tools/mini/fixture/config.mini.yaml"

for t in python3 tar; do
  command -v "$t" >/dev/null 2>&1 || { echo "ERROR: '$t' not in PATH" >&2; exit 2; }
done
HAVE_RSYNC=0; command -v rsync >/dev/null 2>&1 && HAVE_RSYNC=1

echo "== [1/6] clone =="
echo "   repo  : ${REPO}"
echo "   clone : ${CLONE}"
mkdir -p "${CLONE}"
if (( SYNC == 1 )); then
  # The clone's results/logs/scheduler state are never synced, and the rsync below
  # does not clean them either -- a PASS computed from a previous run is worse than
  # a missing one.  KEEP_OUTPUT=1 preserves them on purpose.
  if [[ "${KEEP_OUTPUT:-0}" != "1" ]]; then
    rm -rf "${CLONE}/output" "${CLONE}/.snakemake" "${CLONE}/logs" "${CLONE}/harvest"
  fi
  # code flows in; the fixture itself stays in the clone
  if (( HAVE_RSYNC == 1 )); then
    rsync -a --delete \
      --exclude 'output/' --exclude '.snakemake/' --exclude 'logs/' \
      --exclude 'harvest/' --exclude 'tools/mini/' --exclude '__pycache__/' \
      "${REPO}/" "${CLONE}/" || { echo "ERROR: rsync failed" >&2; exit 2; }
  else
    echo "   note: rsync not found -> tar copy (files deleted from the repo stay in the clone)"
    ( cd "${REPO}" && tar -cf - --exclude=output --exclude=.snakemake --exclude=logs \
        --exclude=harvest --exclude=tools/mini --exclude=__pycache__ . ) \
      | ( cd "${CLONE}" && tar -xf - ) || { echo "ERROR: tar copy failed" >&2; exit 2; }
  fi
fi
mkdir -p "${CLONE}/tools/mini" "${CLONE}/logs"
for f in check_outputs.py make_mini_fixture.py branch_coverage.py test_all_branches.sh; do
  if [ -f "${TOOLS}/${f}" ]; then cp -f "${TOOLS}/${f}" "${CLONE}/tools/mini/"; fi
done

# Non-fatal L3 preflight: the guard rewrites metawrap-modules/reassemble_bins.sh
# via `patch-l3` so its temp lives under a short path. `verify-l3` needs
# `metawrap` on PATH, hence the metawrap env; a miss is only a WARN because a
# short local_base can make the patch unnecessary.
echo "== guard L3 patch preflight (non-fatal) =="
if command -v conda >/dev/null 2>&1; then
  if conda run -n metawrap-env bash "${REPO}/tools/mw_guard.sh" verify-l3 >/dev/null 2>&1; then
    echo "   verify-l3: STATUS patched"
  else
    echo "   WARN: guard verify-l3 did not report 'patched' (short local_base may not need it); run: bash tools/mw_guard.sh patch-l3" >&2
  fi
else
  echo "   note: conda not on PATH; skipping verify-l3 preflight" >&2
fi

echo "== [2/6] fixture config =="
if [[ ! -f "${CLONE}/${CFG_REL}" ]]; then
  cat >&2 <<EOM
ERROR: no fixture yet at ${CLONE}/${CFG_REL}

Create it once (fixture lives in the clone, so it survives every sync):

  cd ${CLONE}
  python3 tools/mini/make_mini_fixture.py --out tools/mini/fixture --repo ${CLONE} \\
      --synthesize \\
      --guard <repo>/tools/mw_guard.sh --local-base <short node-local scratch>
  # or, for a run that can actually pass the MIMAG gate, use real genomes:
  #   --genome <MAG_A>.fna --genome <MAG_B>.fna

Then re-run this script.
EOM
  exit 2
fi
if grep -q '/path/to' "${CLONE}/${CFG_REL}"; then
  echo "ERROR: ${CFG_REL} still contains '/path/to' placeholders (metawrap.guard /" >&2
  echo "       metawrap.local_base). Fix them, then re-run." >&2
  exit 2
fi
python3 - "${CLONE}/${CFG_REL}" <<'PY' || exit 2
import glob, os, sys, yaml
cfg = yaml.safe_load(open(sys.argv[1]))
dd = cfg["data_dir"]
if not os.path.isdir(dd):
    raise SystemExit("ERROR: data_dir does not exist: {0}".format(dd))
r1 = sorted(glob.glob(os.path.join(dd, "*_1.fq.gz")))
print("   data_dir : {0}".format(dd))
print("   samples  : {0}".format([os.path.basename(p)[:-len('_1.fq.gz')] for p in r1]))
print("   binners  : {0}".format(", ".join(cfg["binners"])))
print("   treats   : {0}".format(", ".join(cfg["assembly_treats"])))
print("   methods  : {0}".format(", ".join(cfg["filter_methods"])))
print("   groups   : {0}".format(", ".join(cfg["groups"])))
if not r1:
    raise SystemExit("ERROR: no <sample>_1.fq.gz under data_dir")
PY

# -F matters: without forced execution, jobs whose outputs already exist drop out of
# the plan, so the job-count fingerprint becomes state-dependent (28 -> 26 once
# check_inputs.done and check_qc_inputs.done are present) and every later run warns
# about a "changed" DAG that did not change.
echo "== [3/6] dry run (DAG fingerprint) =="
( cd "${CLONE}" && bash run.sh -n -F --configfile "${CFG_REL}" --cores "${CORES}" ) \
  > "${CLONE}/logs/dry.log" 2>&1
rc=$?
if (( rc != 0 )); then
  echo "ERROR: dry run failed (rc=${rc}); tail of ${CLONE}/logs/dry.log:" >&2
  tail -30 "${CLONE}/logs/dry.log" >&2
  exit 1
fi
EXPECT="${CLONE}/tools/mini/fixture/expected_jobs.json"
if (( BLESS == 1 )) || [[ ! -f "${EXPECT}" ]]; then
  python3 "${CLONE}/tools/mini/check_outputs.py" fingerprint --log "${CLONE}/logs/dry.log" \
    --expect "${EXPECT}" --bless
else
  python3 "${CLONE}/tools/mini/check_outputs.py" fingerprint --log "${CLONE}/logs/dry.log" \
    --expect "${EXPECT}" --check || echo "WARN: DAG fingerprint changed (see DIFF lines above)"
fi
if (( DRY_ONLY == 1 )); then
  echo "== done (--dry-only) =="
  exit 0
fi

echo "== [4/6] full run (this is the slow part: SPAdes, metaWRAP reassembly, GTDB-Tk) =="
force=""; (( FULL == 1 )) && force="-F"
( cd "${CLONE}" && bash run.sh --configfile "${CFG_REL}" --cores "${CORES}" ${force} ) \
  2>&1 | tee "${CLONE}/logs/validate_run.log"
rc=${PIPESTATUS[0]}

echo "== [5/6] assertions =="
python3 "${CLONE}/tools/mini/check_outputs.py" outputs --output-dir "${CLONE}/output" \
  --config "${CLONE}/${CFG_REL}"
checks=$?

echo "== [6/6] summary =="
echo "   snakemake rc    : ${rc}"
echo "   assertion rc    : ${checks}  (0 = no FAIL)"
echo "   run log         : ${CLONE}/logs/validate_run.log"
echo "   outputs         : ${CLONE}/output"
echo "   fixtures/config : ${CLONE}/tools/mini/fixture/"
if (( rc != 0 )); then
  echo "   -> inspect failed jobs: grep -n 'Error in rule' ${CLONE}/logs/validate_run.log"
fi
(( rc == 0 && checks == 0 )) || exit 1
echo "VALIDATE: OK"
