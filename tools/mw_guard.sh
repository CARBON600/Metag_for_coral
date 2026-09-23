#!/usr/bin/env bash
# =============================================================================
# mw_guard.sh -- one-file guard for a metaWRAP 1.3.2 / CheckM binning pipeline
#   Merges runner + probe + preflight + L3 patch into ONE script: the Python
#   probe is embedded and written to a temp file at run time.
#
#   Modes: selftest | check | batch-dry --plan F | batch-run --plan F | run |
#          postrun --samples ID | verify-l3 | patch-l3 | unpatch-l3
#
# THE FAILURE BEING GUARDED (traced to upstream source; locators are file:line of
# the CPython 2.7 / metaWRAP / CheckM sources, not recalled):
#   reassemble_bins.sh:337,381  mkdir ${out}/tmp
#   reassemble_bins.sh:338,382  checkm lineage_wf ... --tmpdir ${out}/tmp ...
#    -> CheckM bin/checkm:45-52  ChangeTempAction: tempfile.tempdir = values
#       (verbatim; $TMPDIR/TMP/TEMP are BYPASSED inside CheckM afterwards)
#    -> checkm/markerGeneFinder.py:68  mp.Manager().dict()
#       -> py2.7 managers.py:528       self._address = reader.recv()  (parent waits)
#          managers.py:162/541-550     Listener(address=None) -> _run_server
#          connection.py:124-127       family = ... or default_family
#          connection.py:64-65         default_family = 'AF_UNIX' on Linux
#          connection.py:90            tempfile.mktemp(prefix='listener-', dir=get_temp_dir())
#          util.py:135-143             tempfile.mkdtemp(prefix='pymp-') + Finalize(rmtree)
#          connection.py:256           self._socket.bind(address)   <-- kernel check
#          include/uapi/linux/un.h     #define UNIX_PATH_MAX 108 (107 usable + NUL)
#       -> child dies BEFORE writer.send(address); parent recv() raises EOFError
#          => CheckM prints "Unexpected error: EOFError" => metaWRAP sees a missing
#          bin_stats_ext.tsv => "Something went wrong with running CheckM."
#   socket path = <--tmpdir value>/pymp-XXXXXX/listener-XXXXXX
#   metaWRAP passes a RELATIVE tmpdir for binning/refinement (${bin_set}.tmp after
#   `cd $out`, ~39 bytes, never fails) and an ABSOLUTE one for reassembly
#   (${out}/tmp) -- which is exactly why only the final reassembly CheckM died.
#   Separately: on NFS, unlink of a still-open file is silly-renamed to
#   ".nfsXXXXXXXX" (fs/nfs/unlink.c:413) and removing that name returns EBUSY;
#   py2.7's cleanup finalizer prints a traceback but util.py:275-280 catches it and
#   only PRINTS -- COSMETIC. The 0920 log proves it: 16x EBUSY + SUCCESS.
#
# NOTHING ASSUMES A FIXED SUFFIX: selftest MEASURES the suffix through the very
# functions multiprocessing uses, and measures the kernel limit with real bind()
# calls, because py2.7 emits 6-char random names (suffix 28) and py3 emits 8 (32).
# =============================================================================
set -Eeuo pipefail
shopt -s nullglob
VERSION="mw_guard 1.1 (2026-09-21)"

# ------------------------------- Tunables ------------------------------------
METAWRAP_ENV="${METAWRAP_ENV:-metawrap-env}"
THREADS="${THREADS:-24}"
MIN_COMPLETENESS="${MIN_COMPLETENESS:-50}"
MAX_CONTAMINATION="${MAX_CONTAMINATION:-10}"
RUN_CONCOCT="${RUN_CONCOCT:-1}"
REFINE_MEM_GB="${REFINE_MEM_GB:-0}"
REASSEMBLE_MEM_GB="${REASSEMBLE_MEM_GB:-0}"
MEM_HEADROOM_GB="${MEM_HEADROOM_GB:-2}"
BASE="${BASE:-/mnt/storage5/yfdai/software/METAWRAP_OFFICIAL_TEST}"
READS_DIR="${READS_DIR:-${BASE}/RAW_READS}"
ASSEMBLY_DIR="${ASSEMBLY_DIR:-${BASE}/RAW_READS/ERR011347_assembly}"
BASE_FINAL="${BASE_FINAL:-${BASE}/METAWRAP_TEST}"
if [[ -n "${MW_LOCAL_BASE:-}" ]]; then LOCAL_BASE="${MW_LOCAL_BASE}"
elif [[ -n "${SLURM_TMPDIR:-}" && -d "${SLURM_TMPDIR}" ]]; then LOCAL_BASE="${SLURM_TMPDIR}"
else LOCAL_BASE="/tmp/${USER:-$(id -un)}/mw"; fi
ALLOW_NETWORK_WORK="${ALLOW_NETWORK_WORK:-0}"
KEEP_SHORT_WORK="${KEEP_SHORT_WORK:-1}"
MAX_STEP_ATTEMPTS="${MAX_STEP_ATTEMPTS:-2}"
RESUME="${RESUME:-0}"
ARCHIVE_EXISTING="${ARCHIVE_EXISTING:-1}"
INSTALL_RMTREE_PATCH="${INSTALL_RMTREE_PATCH:-auto}"
INPUT_STRICT="${INPUT_STRICT:-1}"
OUT_REPORT="${OUT_REPORT:-${PWD}/mw_batch_dryrun.report}"
# Footprint heuristic: DERIVATION with conservative factors, NOT a literature
# value; override freely.
EST_READS_FACTOR="${EST_READS_FACTOR:-3}"
EST_ASM_FACTOR="${EST_ASM_FACTOR:-6}"
EST_FLAT_MB="${EST_FLAT_MB:-1024}"
# Deterministic: a retry cannot change the outcome (cause is layout/limits).
DETERMINISTIC_RE='AF_UNIX path too long|Unexpected error:.*EOFError|Something went wrong with running CheckM'
# Genuinely transient: a retry may succeed after clearing stale temp.
TRANSIENT_RE='No space left on device|Cannot allocate memory|MemoryError|Killed|Connection reset by peer|Temporary failure in name resolution'
# Cosmetic on NFS: reported, never fatal.
ADVISORY_RE='Device or resource busy|\.nfs'

# ------------------------------- Runtime state -------------------------------
SUN_LIMIT=0; SOCK_SUFFIX=0
PYBIN="${PYBIN:-python}"; PROBE_PY=""
SAMPLE_ID=""; SAMPLE_IDS="${SAMPLE_IDS:-}"
WORK_FS=""; NODE_MEM_GB=0
SHORT_WORK=""; SHORT_TMP=""; INITIAL_BINNING=""; BIN_REFINEMENT=""; BIN_REASSEMBLY=""
FINAL_BINS_SHORT=""; LOGDIR_SHORT=""; MW_TMPDIR=""; FINAL_ROOT=""
REFINED_BINS=""; REFINED_STATS=""; REASSEMBLED_BINS=""; REASSEMBLED_STATS=""; REASSEMBLED_CHECKM_TSV=""
INPUT_ERR=""; IN_READS_BYTES=0; IN_ASM_BYTES=0; IN_RECORDS=0
PLAN_NF=0; RUN_R1=""; RUN_R2=""; RUN_ASM=""

# ------------------------------- Helpers ------------------------------------
ts()   { date '+%Y-%m-%d %H:%M:%S'; }
log()  { printf '[%s] %s\n' "$(ts)" "$*"; }
warn() { printf '[%s] WARNING: %s\n' "$(ts)" "$*" >&2; }
die()  { printf '[%s] ERROR: %s\n' "$(ts)" "$*" >&2; exit 1; }
row()  { printf '%s\n' "$*"; printf '%s\n' "$*" >> "${OUT_REPORT}" 2>/dev/null || true; }
on_error() {
    local rc=$? line=${1:-unknown}
    printf '[%s] ERROR: aborted at %s line %s (exit %s)\n' "$(ts)" "${BASH_SOURCE[0]}" "${line}" "${rc}" >&2
    exit "${rc}"
}
trap 'on_error ${LINENO}' ERR

# Every command substitution below is written so a failing inner command cannot
# abort the run under `set -e -o pipefail`.
fs_type_of() {
    local p=$1 t=""
    [[ -e "${p}" ]] || p=$(dirname -- "${p}")
    t=$(findmnt -no FSTYPE -T "${p}" 2>/dev/null | head -n1 || true)
    [[ -n "${t}" ]] || t=$(stat -f -c %T "${p}" 2>/dev/null || true)
    printf '%s\n' "${t}"
}
is_network_fs() {
    case "$1" in
        nfs|nfs4|nfs*) return 0 ;;
        cifs|smb3|smbfs|smb*) return 0 ;;
        ceph|glusterfs|9p|afs) return 0 ;;
        fuse.sshfs|fuse.glusterfs|fuse.s3fs|fuse.gcsfs|fuse.s3ql) return 0 ;;
        *) return 1 ;;
    esac
}
is_tmpfs() { case "$1" in tmpfs|ramfs) return 0 ;; *) return 1 ;; esac; }
detect_mem_limit_gb() {
    local total=0 cg=0 v=0 m=0
    if [[ -r /proc/meminfo ]]; then
        total=$(awk '/^MemTotal:/{printf "%d", $2/1024/1024}' /proc/meminfo 2>/dev/null || true)
    fi
    if [[ -r /sys/fs/cgroup/memory.max ]]; then
        v=$(cat /sys/fs/cgroup/memory.max 2>/dev/null || true)
        [[ "${v}" =~ ^[0-9]+$ ]] && cg=$(( v / 1024 / 1024 / 1024 ))
    elif [[ -r /sys/fs/cgroup/memory/memory.limit_in_bytes ]]; then
        v=$(cat /sys/fs/cgroup/memory/memory.limit_in_bytes 2>/dev/null || true)
        [[ "${v}" =~ ^[0-9]+$ ]] && cg=$(( v / 1024 / 1024 / 1024 ))
    fi
    m="${total}"
    (( cg > 0 && (total == 0 || cg < total) )) && m="${cg}"
    (( m > 0 )) || m=0
    printf '%s\n' "${m}"
}
df_avail_gb()     { df -Pk -- "$1" 2>/dev/null | awk 'NR==2{printf "%d", $4/1024/1024}' || true; }
df_avail_inodes() { df -Pi -- "$1" 2>/dev/null | awk 'NR==2{printf "%d", $4}' || true; }
human_bytes() {
    awk -v b="${1:-0}" 'BEGIN{s="B";split("K M G T P",u," ");for(i=1;i<=5;i++){if(b<1024)break;b/=1024;s=u[i]}printf "%.2f%s",b,s}' || true
}

# gzip-aware input helpers. metaWRAP/bwa read gzipped FASTQ transparently, but
# the binners require a PLAIN assembly, so assembly gz is normalized in
# run_sample. Detection is by magic bytes (1f 8b), never by file name.
is_gzip() { [[ "$(head -c2 -- "$1" 2>/dev/null | od -An -tx1 2>/dev/null | tr -d ' \n')" == "1f8b" ]]; }
gz_cat()  { if is_gzip "$1"; then gzip -dc -- "$1"; else cat -- "$1"; fi; }
decompressed_bytes() {
    local n=0
    if is_gzip "$1"; then
        n=$(gz_cat "$1" 2>/dev/null | wc -c | tr -d '[:space:]') || n=0
    else
        n=$(stat -c %s -- "$1" 2>/dev/null) || n=0
    fi
    printf '%s\n' "${n:-0}"
    return 0
}

# --------------------------- Embedded Python probe ---------------------------
write_probe() {
    [[ -n "${PROBE_PY}" ]] && return 0
    PROBE_PY="${TMPDIR:-/tmp}/mw_guard_probe.$$.py"
    cat > "${PROBE_PY}" <<'PROBE_EOF'
#!/usr/bin/env python
"""Fidelity probe for CheckM's temp / AF_UNIX behaviour.

  --calibrate              measure this interpreter + kernel; print KEY=VALUE
  <tmpdir>                 expect multiprocessing.Manager() to succeed
  <tmpdir> --expect-fail   expect it to fail BECAUSE OF THE PATH LENGTH
"""
from __future__ import print_function

import os
import shutil
import socket
import sys
import tempfile

RC_OK, RC_BAD = 0, 1


def mkpath_of_len(base, n):
    """A directory path of length exactly n under base (caller creates it)."""
    rem = n - len(base)
    if rem < 2:
        return None
    path = base
    while rem > 0:
        take = min(100, rem - 1)
        if take <= 0:
            return None
        path = path + os.sep + "D" * take
        rem -= 1 + take
    return path


def measure_suffix_len():
    """Exactly the upstream path: util.get_temp_dir() + arbitrary_address()."""
    td = tempfile.mkdtemp(prefix="cal-")
    old = tempfile.tempdir
    try:
        tempfile.tempdir = td
        pymp = tempfile.mkdtemp(prefix="pymp-")
        listener = tempfile.mktemp(prefix="listener-", dir=pymp)
        return len(listener) - len(td)
    finally:
        tempfile.tempdir = old
        shutil.rmtree(td, ignore_errors=True)


def measure_sun_path_limit():
    """Largest socket path length the kernel accepts, by real bind() calls."""
    base = tempfile.mkdtemp(prefix="lim-")
    limit = 0
    reason = ""
    try:
        total = len(base) + 4
        while total < len(base) + 160:
            d = mkpath_of_len(base, total - 2)          # + "/s" == total
            if d is None:
                break
            try:
                os.makedirs(d)
            except OSError as exc:
                reason = "mkdir failed at %d: %s" % (total, exc)
                break
            p = os.path.join(d, "s")
            try:
                sk = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
                sk.bind(p)
                sk.close()
                os.unlink(p)
            except (OSError, socket.error) as exc:
                reason = "bind failed at %d: %s" % (total, exc)
                break
            limit = max(limit, len(p))
            total += 1
        return limit, reason
    finally:
        shutil.rmtree(base, ignore_errors=True)


def calibrate():
    suffix = measure_suffix_len()
    limit, reason = measure_sun_path_limit()
    have_unix = 1 if hasattr(socket, "AF_UNIX") else 0
    print("interpreter=%s" % sys.version.split()[0])
    print("py_major=%d" % sys.version_info[0])
    print("af_unix=%d" % have_unix)
    print("suffix_len=%d" % suffix)
    print("sun_path_limit=%d" % limit)
    print("limit_probe_note=%s" % (reason or "clean"))
    return RC_OK if (limit > 0 and have_unix) else RC_BAD


def one(tmpdir, expect_fail):
    if not os.path.isdir(tmpdir):
        print("ERROR: tmpdir must already exist (CheckM's ChangeTempAction runs "
              "os.path.isdir first): %s" % tmpdir)
        return 2

    suffix = measure_suffix_len()
    limit, _ = measure_sun_path_limit()
    pred = tmpdir + os.sep + "pymp-XXXXXX" + os.sep + "listener-XXXXXX"
    over = len(pred) > limit

    print("interpreter      = %s" % sys.version.split()[0])
    print("tempdir          = %s" % tmpdir)
    print("predicted socket = %s" % pred)
    print("predicted length = %d   (measured kernel limit %d, measured suffix %d)"
          % (len(pred), limit, suffix))
    print("over_limit       = %s" % ("yes" if over else "no"))

    tempfile.tempdir = tmpdir                    # == CheckM ChangeTempAction
    pymp = None
    passed = False
    detail = "no exception"
    try:
        pymp = tempfile.mkdtemp(prefix="pymp-")
        import multiprocessing as mp
        manager = mp.Manager()
        try:
            d = manager.dict()
            d["status"] = "PASS"
            passed = (d["status"] == "PASS")
        finally:
            manager.shutdown()
    except Exception as exc:                     # report verbatim
        detail = "%s: %s" % (type(exc).__name__, exc)
    finally:
        if pymp and os.path.isdir(pymp):
            shutil.rmtree(pymp, ignore_errors=True)

    if expect_fail:
        if passed:
            print("RESULT: UNEXPECTED PASS (a failure was required) -> %s" % detail)
            return RC_BAD
        if over:
            print("RESULT: FAIL as expected, ATTRIBUTED TO LENGTH (predicted %d > "
                  "limit %d); parent saw: %s" % (len(pred), limit, detail))
            return RC_OK
        print("RESULT: FAILED FOR THE WRONG REASON -- predicted %d <= limit %d, so "
              "this is NOT the AF_UNIX length limit; parent saw: %s"
              % (len(pred), limit, detail))
        return RC_BAD

    if passed:
        print("RESULT: PASS")
        return RC_OK
    print("RESULT: FAIL -> %s" % detail)
    return RC_BAD


def main(argv):
    if len(argv) < 2 or argv[1] in ("-h", "--help"):
        print(__doc__)
        return 2
    if argv[1] == "--calibrate":
        return calibrate()
    return one(argv[1], "--expect-fail" in argv[2:])


if __name__ == "__main__":
    sys.exit(main(sys.argv))
PROBE_EOF
    chmod 700 "${PROBE_PY}"
    trap 'rm -f "${PROBE_PY}"' EXIT
}

calibrate() {
    write_probe
    local out af
    out=$("${PYBIN}" "${PROBE_PY}" --calibrate 2>&1) || { printf '%s\n' "${out}"; die "probe calibration failed under ${PYBIN}"; }
    SUN_LIMIT=$(printf '%s\n' "${out}" | awk -F= '/^sun_path_limit=/{print $2}')
    SOCK_SUFFIX=$(printf '%s\n' "${out}" | awk -F= '/^suffix_len=/{print $2}')
    af=$(printf '%s\n' "${out}" | awk -F= '/^af_unix=/{print $2}')
    [[ -n "${SUN_LIMIT}" && -n "${SOCK_SUFFIX}" ]] || { printf '%s\n' "${out}"; die "cannot parse calibration output"; }
    (( SUN_LIMIT > 0 )) || die "no AF_UNIX socket could be bound here (measured limit 0)"
    [[ "${af}" == "1" ]] || die "AF_UNIX unavailable: CheckM cannot run here at all"
    log "calibrated: $(printf '%s\n' "${out}" | awk -F= '/^interpreter=/{print $2}') | measured suffix=${SOCK_SUFFIX} | measured sun_path_limit=${SUN_LIMIT}"
    (( SOCK_SUFFIX == 28 )) || { warn "measured suffix is ${SOCK_SUFFIX}, not 28 -- this is NOT the CPython 2.7 CheckM runs under."; warn "All length verdicts use the MEASURED value so they stay correct; run under ${METAWRAP_ENV}'s python2.7 to certify the production path."; }
}

probe_pair() {
    local tmpdir=$1 tag=${2:-probe} pad neg
    mkdir -p "${tmpdir}"
    log "-- positive probe (a directory the pipeline really uses): ${tmpdir}"
    "${PYBIN}" "${PROBE_PY}" "${tmpdir}" || die "positive probe FAILED at ${tmpdir}"
    pad=$(printf 'x%.0s' $(seq 1 $(( SUN_LIMIT - SOCK_SUFFIX + 10 ))))
    neg="${SHORT_WORK:-${LOCAL_BASE}}/neg_${pad}"
    mkdir -p "${neg}"
    log "-- negative control: an over-long tmpdir must fail FOR THE LENGTH REASON"
    if ! "${PYBIN}" "${PROBE_PY}" "${neg}" --expect-fail; then
        rm -rf "${neg}"
        die "negative control failed for a reason other than length -- the AF_UNIX root cause is NOT confirmed on this box; refusing to claim the guard is meaningful"
    fi
    rm -rf "${neg}"
    log "probe pair OK (${tag})"
}

# ------------------------------ Path planning -------------------------------
# Per-sample and collision-free by construction. FINAL_ROOT is per sample: a
# shared FINAL_ROOT made sample N+1 archive (mv) sample N's results away.
plan_sample() {
    SAMPLE_ID=$1
    # MW_WORK_DIR pins the active tree to an EXISTING run: without it, SHORT_WORK
    # embeds $$, so a later `postrun` in a new shell would plan a different path
    # than the run used and never find its logs.
    if [[ -n "${MW_WORK_DIR:-}" ]]; then
        SHORT_WORK="${MW_WORK_DIR}"
    else
        SHORT_WORK="${LOCAL_BASE}/w_${SAMPLE_ID}_$$"
    fi
    SHORT_TMP="${LOCAL_BASE}/t_${SAMPLE_ID}_$$"
    INITIAL_BINNING="${SHORT_WORK}/01_BIN"
    BIN_REFINEMENT="${SHORT_WORK}/02_REF"
    BIN_REASSEMBLY="${SHORT_WORK}/03_REA"
    FINAL_BINS_SHORT="${SHORT_WORK}/04_MAG"
    LOGDIR_SHORT="${SHORT_WORK}/log"
    MW_TMPDIR="${MW_TMPDIR_OVERRIDE:-${SHORT_WORK}/mt}"
    FINAL_ROOT="${BASE_FINAL}/${SAMPLE_ID}"
    REFINED_BINS="${BIN_REFINEMENT}/metawrap_${MIN_COMPLETENESS}_${MAX_CONTAMINATION}_bins"
    REFINED_STATS="${REFINED_BINS}.stats"
    REASSEMBLED_BINS="${BIN_REASSEMBLY}/reassembled_bins"
    REASSEMBLED_STATS="${REASSEMBLED_BINS}.stats"
    REASSEMBLED_CHECKM_TSV="${REASSEMBLED_BINS}.checkm/storage/bin_stats_ext.tsv"
}
reassemble_socket_len() { printf '%d\n' $(( ${#BIN_REASSEMBLY} + 4 + SOCK_SUFFIX )); }  # 4 = "/tmp"
derived_work_cap()      { printf '%d\n' $(( SUN_LIMIT - 4 - SOCK_SUFFIX - 7 )); }       # 7 = "/03_REA"

# ------------------------------ Env / tools ---------------------------------
activate_conda_env() {
    local conda_base
    if [[ -n "${CONDA_EXE:-}" ]]; then conda_base="$("${CONDA_EXE}" info --base)"
    elif command -v conda >/dev/null 2>&1; then conda_base="$(conda info --base)"
    else warn "conda not in PATH"; return 1; fi
    set +u
    source "${conda_base}/etc/profile.d/conda.sh"
    conda activate "${METAWRAP_ENV}" || { set -u; warn "conda activate ${METAWRAP_ENV} failed"; return 1; }
    set -u
    return 0
}
# Tools the modules actually invoke (binning.sh:215-329 and reassemble_bins.sh:177-234
# use bwa / jgi_summarize_bam_contig_depths / samtools / spades.py; bowtie2 is NOT used).
required_tools() {
    local req=(metawrap metabat2 run_MaxBin.pl checkm bwa jgi_summarize_bam_contig_depths samtools spades.py md5sum awk find gzip od)
    (( RUN_CONCOCT == 1 )) && req+=(concoct cut_up_fasta.py concoct_coverage_table.py merge_cutup_clustering.py)
    printf '%s\n' "${req[@]}"
}
# One candidate list for BOTH detection and patching, so they cannot disagree.
find_mw_module() {
    local name=$1 mwbin cand
    mwbin=$(dirname -- "$(command -v metawrap 2>/dev/null || true)")
    [[ -n "${mwbin}" && "${mwbin}" != "." ]] || return 1
    for cand in \
        "${mwbin}/metawrap-modules/${name}" \
        "${mwbin}/../bin/metawrap-modules/${name}" \
        "${mwbin}/../metawrap-modules/${name}" \
        "${mwbin}/${name}"; do
        [[ -f "${cand}" ]] && { printf '%s\n' "${cand}"; return 0; }
    done
    cand=$(find "$(dirname -- "${mwbin}")" -maxdepth 5 -name "${name}" 2>/dev/null | head -n1 || true)
    [[ -n "${cand}" ]] && { printf '%s\n' "${cand}"; return 0; }
    return 1
}
audit_module_tmpdirs() {
    local name path rel line
    for name in binning.sh bin_refinement.sh reassemble_bins.sh; do
        path=$(find_mw_module "${name}" || true)
        [[ -n "${path}" ]] || { log "  ${name}: not found"; continue; }
        while IFS= read -r line; do
            case "${line}" in
                *'--tmpdir ${out}'*|*'--tmpdir ${1}'*|*'--tmpdir /'*) rel="ABSOLUTE-constrained-by-O-length" ;;
                *'--tmpdir ${'*|*'.tmp'*) rel="relative-never-exceeds" ;;
                *) rel="?" ;;
            esac
            log "  ${name} [${rel}] $(printf '%s' "${line}" | sed 's/^[[:space:]]*//' | cut -c1-100)"
        done < <(grep -n -- '--tmpdir' "${path}" 2>/dev/null || true)
    done
}
mw_patched() {
    local p; p=$(find_mw_module reassemble_bins.sh || true)
    [[ -n "${p}" ]] && grep -q 'MW_TMPDIR' "${p}" 2>/dev/null
}

# ------------------------------- Input checks -------------------------------
# Returns 1 and sets INPUT_ERR instead of dying, so a preflight can report on
# every sample before aborting (a `die` inside a condition would kill the report).
check_inputs() {
    local r1=$1 r2=$2 asm=$3 l1 l2 b1 b2 hdr
    INPUT_ERR=""; IN_READS_BYTES=0; IN_ASM_BYTES=0; IN_RECORDS=0
    [[ -f "${r1}" && -s "${r1}" ]] || { INPUT_ERR="missing/empty READ1: ${r1}"; return 1; }
    [[ -f "${r2}" && -s "${r2}" ]] || { INPUT_ERR="missing/empty READ2: ${r2}"; return 1; }
    [[ -f "${asm}" && -s "${asm}" ]] || { INPUT_ERR="missing/empty assembly: ${asm}"; return 1; }
    b1=$(gz_cat "${r1}" 2>/dev/null | head -c 1 || true)
    b2=$(gz_cat "${r2}" 2>/dev/null | head -c 1 || true)
    hdr=$(gz_cat "${asm}" 2>/dev/null | grep -m1 '^>' || true)
    [[ "${b1}" == '@' ]] || { INPUT_ERR="READ1 is not FASTQ/.fq(.gz): ${r1}"; return 1; }
    [[ "${b2}" == '@' ]] || { INPUT_ERR="READ2 is not FASTQ/.fq(.gz): ${r2}"; return 1; }
    [[ -n "${hdr}" ]]   || { INPUT_ERR="assembly has no FASTA header: ${asm}"; return 1; }
    # Footprint must use the DECOMPRESSED size, otherwise gz inputs are
    # systematically under-estimated (the pipeline expands them on disk).
    IN_READS_BYTES=$(( $(decompressed_bytes "${r1}") + $(decompressed_bytes "${r2}") ))
    IN_ASM_BYTES=$(decompressed_bytes "${asm}")
    if (( INPUT_STRICT == 1 )); then
        l1=$(gz_cat "${r1}" | wc -l) || l1=0
        l2=$(gz_cat "${r2}" | wc -l) || l2=0
        (( l1 % 4 == 0 )) || { INPUT_ERR="READ1 line count not divisible by 4"; return 1; }
        (( l2 % 4 == 0 )) || { INPUT_ERR="READ2 line count not divisible by 4"; return 1; }
        (( l1 == l2 ))   || { INPUT_ERR="R1/R2 line counts differ (${l1} vs ${l2})"; return 1; }
        IN_RECORDS=$(( l1 / 4 ))
    fi
    return 0
}
est_footprint_bytes() {
    printf '%d\n' $(( EST_READS_FACTOR * IN_READS_BYTES + EST_ASM_FACTOR * IN_ASM_BYTES + EST_FLAT_MB * 1024 * 1024 ))
}

# ---------------------------- Pipeline execution ----------------------------
run_logged() {
    local label=$1; shift
    mkdir -p "${LOGDIR_SHORT}"
    local logfile="${LOGDIR_SHORT}/${label}.log" rc=0
    {
        printf '[%s] CWD=%s\n' "$(ts)" "${SHORT_WORK}"
        printf '[%s] TMPDIR=%s MW_TMPDIR=%s\n' "$(ts)" "${SHORT_TMP}" "${MW_TMPDIR}"
        printf '[%s] Command:' "$(ts)"; printf ' %q' "$@"; printf '\n'
    } >> "${logfile}"
    set +e
    (
        cd "${SHORT_WORK}"
        env TMPDIR="${SHORT_TMP}" TMP="${SHORT_TMP}" TEMP="${SHORT_TMP}" TEMPDIR="${SHORT_TMP}" \
            MW_TMPDIR="${MW_TMPDIR}" "$@"
    ) >> "${logfile}" 2>&1
    rc=$?
    set -e
    tail -n 4 "${logfile}" 2>/dev/null | sed 's/^/    | /' || true
    return "${rc}"
}
sig_in_log() { local f="${LOGDIR_SHORT}/$1.log"; [[ -s "${f}" ]] || return 1; grep -Eq "$2" "${f}"; }

clean_transient_state() {
    rm -rf "${SHORT_TMP}"; mkdir -p "${SHORT_TMP}"; chmod 700 "${SHORT_TMP}"
    find "${SHORT_WORK}" -depth -name 'pymp-*' -exec rm -rf {} + 2>/dev/null || true
    find "${SHORT_WORK}" -depth -name '.nfs*' -exec rm -rf {} + 2>/dev/null || true
}
fatal_step() { die "$1 failed (exit $2): $3 -- see ${LOGDIR_SHORT}/$1.log"; }

# Retry ONLY genuinely transient failures; deterministic ones fail fast so a
# multi-hour reassembly is not repeated pointlessly. Each attempt's log is
# rotated aside, otherwise a successful retry still trips the acceptance grep
# (the previous script's `tee -a` made the retry self-defeating).
run_step() {
    local label=$1; shift
    local attempt=1 rc=0
    while :; do
        rc=0; run_logged "${label}" "$@" || rc=$?
        (( rc == 0 )) && { log "completed ${label}"; return 0; }
        if sig_in_log "${label}" "${DETERMINISTIC_RE}"; then
            fatal_step "${label}" "${rc}" "deterministic signature (AF_UNIX/EOF/CheckM) -- a retry cannot help"
        fi
        if (( attempt >= MAX_STEP_ATTEMPTS )) || ! sig_in_log "${label}" "${TRANSIENT_RE}"; then
            fatal_step "${label}" "${rc}" "no retryable transient signature"
        fi
        attempt=$(( attempt + 1 ))
        mv -f "${LOGDIR_SHORT}/${label}.log" "${LOGDIR_SHORT}/${label}.attempt$(( attempt - 1 )).log" 2>/dev/null || true
        warn "transient failure in ${label}; rotated log, clearing temp, retry ${attempt}/${MAX_STEP_ATTEMPTS}"
        clean_transient_state
    done
}
report_advisory() {
    if sig_in_log "$1" "${ADVISORY_RE}"; then
        warn "advisory only: NFS silly-rename (.nfs/EBUSY) noise in $1.log -- non-fatal (util.py:275-280 catches it and only prints; the 0920 log shows 16x EBUSY + SUCCESS). Node-local storage removes it."
    fi
}
install_rmtree_patch() {
    local site="${CONDA_PREFIX:-}/lib/python2.7/site-packages" target
    target="${site}/sitecustomize.py"
    if [[ -z "${CONDA_PREFIX:-}" || ! -d "${site}" ]]; then
        warn "no python2.7 site-packages under '${CONDA_PREFIX:-<unset>}'; skipping rmtree patch"; return 0
    fi
    grep -q 'MW_NFS_RMTREE_PATCH' "${target}" 2>/dev/null && { log "rmtree patch already present"; return 0; }
    [[ -f "${target}" ]] && { cp -a "${target}" "${target}.bak_$(date +%Y%m%d_%H%M%S)"; printf '\n' >> "${target}"; }
    cat >> "${target}" <<'PY'
# MW_NFS_RMTREE_PATCH -- ignore ONLY NFS silly-rename residue (.nfs* still open
# elsewhere) while tearing temp dirs down. multiprocessing/util.py:138 imports
# shutil INSIDE get_temp_dir() and reads the attribute at call time, so replacing
# shutil.rmtree at interpreter startup is sufficient; there is no
# multiprocessing.util.rmtree to patch.
import errno as _errno
import os as _os
import shutil as _shutil

_MW_ORIG_RMTREE = _shutil.rmtree


def _mw_rmtree(path, ignore_errors=False, onerror=None):
    def _onerror(func, p, exc_info):
        exc = exc_info[1]
        if getattr(exc, "errno", None) in (_errno.EBUSY, _errno.ENOENT) \
                and ".nfs" in _os.path.basename(p):
            return
        if onerror is not None:
            onerror(func, p, exc_info)
        elif not ignore_errors:
            raise
    return _MW_ORIG_RMTREE(path, ignore_errors=ignore_errors, onerror=_onerror)


_shutil.rmtree = _mw_rmtree
PY
    log "rmtree patch installed at ${target}"
}
rmtree_patch_wanted() {
    case "${INSTALL_RMTREE_PATCH}" in
        1|true|yes) return 0 ;;
        0|false|no) return 1 ;;
        auto|AUTO|"") is_network_fs "$(fs_type_of "${LOCAL_BASE}")" ;;
        *) return 1 ;;
    esac
}

# ------------------------------ One sample -----------------------------------
run_sample() {
    local r1=$1 r2=$2 asm=$3 n1 n2 n3 f base n_mag broken copied sock
    rm -rf "${SHORT_WORK}" "${SHORT_TMP}" "${MW_TMPDIR}"
    mkdir -p "${SHORT_WORK}" "${SHORT_TMP}" "${LOGDIR_SHORT}" "${MW_TMPDIR}"
    chmod 700 "${SHORT_WORK}" "${SHORT_TMP}" "${MW_TMPDIR}"
    export TMPDIR="${SHORT_TMP}" TMP="${SHORT_TMP}" TEMP="${SHORT_TMP}" TEMPDIR="${SHORT_TMP}" MW_TMPDIR
    [[ -n "${WORK_FS}" ]] || WORK_FS=$(fs_type_of "${LOCAL_BASE}")
    log "sample=${SAMPLE_ID} work=${SHORT_WORK} final=${FINAL_ROOT}"
    log "memory: refine -m=${REFINE_MEM_GB} reassemble -m=${REASSEMBLE_MEM_GB} (usable ${NODE_MEM_GB} GB)"
    # bwa mem reads gzipped FASTQ transparently, so the read symlinks can point at
    # .fq/.fq.gz alike (the _1.fastq/_2.fastq names satisfy metaWRAP's own check).
    ln -sfn "$(readlink -f "${r1}")" "${SHORT_WORK}/x_1.fastq"
    ln -sfn "$(readlink -f "${r2}")" "${SHORT_WORK}/x_2.fastq"
    # The binners (metabat2/maxbin2/concoct) require a PLAIN assembly, so a gz
    # assembly is decompressed once into the work tree instead of symlinked.
    if is_gzip "${asm}"; then
        gz_cat "${asm}" > "${SHORT_WORK}/a.fa"
    else
        ln -sfn "$(readlink -f "${asm}")" "${SHORT_WORK}/a.fa"
    fi
    rmtree_patch_wanted && install_rmtree_patch

    if (( RESUME == 1 )) && compgen -G "${INITIAL_BINNING}/metabat2_bins/*.fa" >/dev/null; then
        log "resume: 01_initial_binning already complete"
    else
        local bcmd=(metawrap binning -o "${INITIAL_BINNING}" -t "${THREADS}" -a "${SHORT_WORK}/a.fa" --metabat2 --maxbin2)
        (( RUN_CONCOCT == 1 )) && bcmd+=(--concoct)
        bcmd+=("${SHORT_WORK}/x_1.fastq" "${SHORT_WORK}/x_2.fastq")
        run_step 01_initial_binning "${bcmd[@]}"
    fi
    n1=$(find "${INITIAL_BINNING}/metabat2_bins" -maxdepth 1 -name '*.fa' 2>/dev/null | wc -l || true)
    n2=$(find "${INITIAL_BINNING}/maxbin2_bins" -maxdepth 1 -name '*.fa' 2>/dev/null | wc -l || true)
    n3=0; (( RUN_CONCOCT == 1 )) && n3=$(find "${INITIAL_BINNING}/concoct_bins" -maxdepth 1 -name '*.fa' 2>/dev/null | wc -l || true)
    (( n1 + n2 + n3 > 0 )) || die "no initial bins produced"
    log "initial bins: metabat2=${n1} maxbin2=${n2} concoct=${n3}"

    if (( RESUME == 1 )) && compgen -G "${BIN_REFINEMENT}/metawrap_*_bins.stats" >/dev/null; then
        log "resume: 02_bin_refinement already complete"
    else
        local rcmd=(metawrap bin_refinement -o "${BIN_REFINEMENT}" -t "${THREADS}" -m "${REFINE_MEM_GB}"
                    -A "${INITIAL_BINNING}/metabat2_bins" -B "${INITIAL_BINNING}/maxbin2_bins"
                    -c "${MIN_COMPLETENESS}" -x "${MAX_CONTAMINATION}")
        (( RUN_CONCOCT == 1 )) && rcmd+=(-C "${INITIAL_BINNING}/concoct_bins")
        run_step 02_bin_refinement "${rcmd[@]}"
    fi
    report_advisory 02_bin_refinement
    # Discover the real stats name instead of trusting one literal (the 0920 run
    # died solely because the wrapper expected metawrap_50_10.stats).
    if [[ ! -s "${REFINED_STATS}" ]]; then
        local found=""
        found=$(find "${BIN_REFINEMENT}" -maxdepth 1 -type f -name 'metawrap_*_bins.stats' 2>/dev/null | sort | head -n1 || true)
        [[ -n "${found}" ]] || die "no metawrap_*_bins.stats in ${BIN_REFINEMENT}"
        REFINED_STATS="${found}"; REFINED_BINS="${found%.stats}"
        log "discovered refined stats: ${REFINED_STATS}"
    fi
    [[ -s "${REFINED_STATS}" ]] || die "refined stats empty: ${REFINED_STATS}"
    [[ -d "${REFINED_BINS}" ]]   || die "refined bins dir missing: ${REFINED_BINS}"

    rm -rf "${SHORT_TMP}"; mkdir -p "${SHORT_TMP}"; chmod 700 "${SHORT_TMP}"
    probe_pair "${MW_TMPDIR}" before-reassembly
    sock=$(reassemble_socket_len)
    if (( sock > SUN_LIMIT )) && ! mw_patched; then
        die "reassemble -o socket path would be ${sock} > ${SUN_LIMIT} (${BIN_REASSEMBLY}/tmp). Shorten the path or run patch-l3."
    elif (( sock > SUN_LIMIT )); then
        log "reassemble -o socket would be ${sock} > ${SUN_LIMIT}, but the L3 patch routes CheckM temp to ${MW_TMPDIR} -- continuing"
    fi
    if (( RESUME == 1 )) && [[ -s "${REASSEMBLED_STATS}" ]] && [[ -s "${REASSEMBLED_CHECKM_TSV}" ]]; then
        log "resume: 03_reassemble_bins already complete"
    else
        run_step 03_reassemble_bins \
            metawrap reassemble_bins -o "${BIN_REASSEMBLY}" \
            -1 "${SHORT_WORK}/x_1.fastq" -2 "${SHORT_WORK}/x_2.fastq" \
            -t "${THREADS}" -m "${REASSEMBLE_MEM_GB}" \
            -c "${MIN_COMPLETENESS}" -x "${MAX_CONTAMINATION}" -b "${REFINED_BINS}"
    fi
    report_advisory 03_reassemble_bins
    if sig_in_log 03_reassemble_bins 'AF_UNIX path too long|Unexpected error:.*EOFError|Something went wrong with running CheckM'; then
        die "fatal CheckM signature in 03_reassemble_bins.log"
    fi
    grep -q 'REASSEMBLY PIPELINE SUCCESSFULLY FINISHED' "${LOGDIR_SHORT}/03_reassemble_bins.log" \
        || die "reassembly SUCCESS banner not found"
    [[ -s "${REASSEMBLED_STATS}" ]]      || die "missing/empty ${REASSEMBLED_STATS}"
    [[ -s "${REASSEMBLED_CHECKM_TSV}" ]] || die "missing/empty ${REASSEMBLED_CHECKM_TSV} (metaWRAP's own CheckM success test, reassemble_bins.sh:339,383)"

    mkdir -p "${FINAL_BINS_SHORT}"
    for f in "${REASSEMBLED_BINS}"/*.fa "${REASSEMBLED_BINS}"/*.fasta "${REASSEMBLED_BINS}"/*.fna; do
        [[ -f "${f}" ]] || continue
        base=$(basename "${f}"); base=${base%.fa}; base=${base%.fasta}; base=${base%.fna}
        ln -sfn "$(readlink -f "${f}")" "${FINAL_BINS_SHORT}/${base}.fa"
    done
    n_mag=$(find "${FINAL_BINS_SHORT}" -maxdepth 1 -type l -name '*.fa' 2>/dev/null | wc -l || true)
    (( n_mag > 0 )) || die "no final MAG links created"
    broken=$(find "${FINAL_BINS_SHORT}" -maxdepth 1 -type l ! -exec test -e {} \; -print 2>/dev/null || true)
    [[ -z "${broken}" ]] || die "broken MAG links: ${broken}"

    local manifest="${SHORT_WORK}/final_bins_manifest.tsv" chk="${SHORT_WORK}/final_bins.md5"
    { printf 'bin_id\tfinal_path\tsource_path\n'
      for f in "${FINAL_BINS_SHORT}"/*.fa; do
          printf '%s\t%s\t%s\n' "$(basename "${f}" .fa)" "${f}" "$(readlink -f "${f}")"
      done; } > "${manifest}"
    ( cd "${FINAL_BINS_SHORT}" && md5sum ./*.fa ) > "${chk}"

    mkdir -p "${FINAL_ROOT}/04_FINAL_BINS_FOR_GTDB" "${FINAL_ROOT}/logs"
    cp -a "${INITIAL_BINNING}" "${FINAL_ROOT}/01_INITIAL_BINNING"
    cp -a "${BIN_REFINEMENT}"  "${FINAL_ROOT}/02_BIN_REFINEMENT"
    cp -a "${BIN_REASSEMBLY}"  "${FINAL_ROOT}/03_BIN_REASSEMBLY"
    cp -L "${FINAL_BINS_SHORT}"/*.fa "${FINAL_ROOT}/04_FINAL_BINS_FOR_GTDB/"
    cp -a "${LOGDIR_SHORT}"/. "${FINAL_ROOT}/logs/"
    cp "${manifest}" "${FINAL_ROOT}/final_bins_manifest.tsv"
    cp "${chk}" "${FINAL_ROOT}/final_bins.md5"
    copied=$(find "${FINAL_ROOT}/04_FINAL_BINS_FOR_GTDB" -maxdepth 1 -type f -name '*.fa' 2>/dev/null | wc -l || true)
    (( copied == n_mag )) || die "copy-back count mismatch: ${copied} != ${n_mag}"
    ( cd "${FINAL_ROOT}/04_FINAL_BINS_FOR_GTDB" && md5sum ./*.fa | sed 's#  \./#  #' | sort -k2 ) > "${FINAL_ROOT}/final_bins.permanent.md5"
    sed 's#  \./#  #' "${chk}" | sort -k2 > "${SHORT_WORK}/expected.sorted.md5"
    cmp -s "${SHORT_WORK}/expected.sorted.md5" "${FINAL_ROOT}/final_bins.permanent.md5" \
        || die "permanent-copy MD5 verification failed"
    { echo "metaWRAP sample summary"; echo "sample: ${SAMPLE_ID}"; echo "completed: $(ts)"
      echo "status: SUCCESS"; echo "gtdb_tk: NOT RUN"
      echo "active work fs: ${WORK_FS:-unknown} (${SHORT_WORK})"
      echo "reassemble -o socket bytes: ${sock}/${SUN_LIMIT} (measured suffix ${SOCK_SUFFIX})"
      echo "l3 patch: $(mw_patched && echo yes || echo no)"
      echo "usable memory GB: ${NODE_MEM_GB}"; echo "reassemble -m GB: ${REASSEMBLE_MEM_GB}"
      echo "initial bins: metabat2=${n1} maxbin2=${n2} concoct=${n3}"
      echo "final MAGs: ${copied}"; echo "MAG dir: ${FINAL_ROOT}/04_FINAL_BINS_FOR_GTDB"; } > "${FINAL_ROOT}/VALIDATION_SUMMARY.txt"
    log "sample ${SAMPLE_ID} SUCCESS -> ${FINAL_ROOT}/04_FINAL_BINS_FOR_GTDB (${copied} MAGs, md5 verified)"
    rm -rf "${SHORT_TMP}"
    (( KEEP_SHORT_WORK == 0 )) && rm -rf "${SHORT_WORK}"
    return 0
}

# ================================ MODES =====================================
mode_selftest() {
    log "${VERSION}"
    command -v conda >/dev/null 2>&1 && { activate_conda_env || true; }
    if [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/python2.7" ]]; then
        PYBIN="${CONDA_PREFIX}/bin/python2.7"; log "using the env's python2.7 as the production interpreter"
    else
        warn "no env python2.7 found; probing with '${PYBIN}' (verdicts stay correct because they are measured)"
    fi
    calibrate
    plan_sample "${SAMPLE_IDS:-SELFTEST}"
    probe_pair "${MW_TMPDIR}" selftest
    log "LOCAL_BASE ${LOCAL_BASE} -> fs '$(fs_type_of "${LOCAL_BASE}")'"
    log "SELFTEST: OK"
}

mode_check() {
    plan_sample "${SAMPLE_IDS:-SAMPLE}"
    local sid="${SAMPLE_ID}"
    local r1="${READS_DIR}/${sid}_1.fastq" r2="${READS_DIR}/${sid}_2.fastq"
    local asm="${ASSEMBLY_DIR}/final.contigs.fa" st=0 k v sock cap need avail wfs
    log "== path plan (derived cap for SHORT_WORK = $(derived_work_cap)) =="
    for k in SHORT_WORK SHORT_TMP INITIAL_BINNING BIN_REFINEMENT BIN_REASSEMBLY MW_TMPDIR FINAL_ROOT; do
        eval "v=\${$k}"
        log "  $(printf '%-17s' "${k}") len=$(printf '%3d' "${#v}")  ${v}"
    done
    sock=$(reassemble_socket_len); cap=$(derived_work_cap)
    if (( sock <= SUN_LIMIT )); then log "  $(printf '%-17s' 'socket(reassembly)') len=$(printf '%3d' "${sock}") limit=${SUN_LIMIT}  OK"
    else log "  socket(reassembly) ${sock} > ${SUN_LIMIT}: patch-l3 or shorten LOCAL_BASE"; st=1; fi
    (( ${#SHORT_WORK} <= cap )) || { log "  ERROR: SHORT_WORK exceeds the derived cap by $(( ${#SHORT_WORK} - cap ))"; st=1; }

    wfs=$(fs_type_of "${LOCAL_BASE}")
    log "== filesystem =="
    log "  LOCAL_BASE ${LOCAL_BASE} -> '${wfs:-unknown}'"
    if is_network_fs "${wfs}"; then
        if (( ALLOW_NETWORK_WORK == 1 )); then log "  WARNING: network FS allowed; .nfs/EBUSY will appear (advisory only)"
        else log "  ERROR: network FS for ACTIVE work; CheckM writes temp inside -o. Move LOCAL_BASE."; st=1; fi
    fi
    is_tmpfs "${wfs}" && log "  WARNING: tmpfs (memory-backed): large .bam / SPAdes output consumes RAM."
    log "  FINAL_ROOT ${FINAL_ROOT} -> '$(fs_type_of "${BASE_FINAL}")'"
    df -h  -- "${LOCAL_BASE}" 2>/dev/null | sed 's/^/  /' || true
    df -Pi -- "${LOCAL_BASE}" 2>/dev/null | sed 's/^/  /' || true

    log "== memory =="
    NODE_MEM_GB=$(detect_mem_limit_gb)
    (( NODE_MEM_GB > 0 )) || { NODE_MEM_GB=40; warn "memory undetectable; assuming ${NODE_MEM_GB} GB"; }
    (( REFINE_MEM_GB <= 0 )) && REFINE_MEM_GB="${NODE_MEM_GB}"
    if (( REASSEMBLE_MEM_GB <= 0 )); then
        if (( NODE_MEM_GB > MEM_HEADROOM_GB )); then REASSEMBLE_MEM_GB=$(( NODE_MEM_GB - MEM_HEADROOM_GB )); else REASSEMBLE_MEM_GB="${NODE_MEM_GB}"; fi
    fi
    log "  usable ${NODE_MEM_GB} GB -> refine -m=${REFINE_MEM_GB}, reassemble -m=${REASSEMBLE_MEM_GB} (spades.py -m)"
    (( REASSEMBLE_MEM_GB <= NODE_MEM_GB )) || { log "  ERROR: reassemble -m exceeds usable RAM -> OOM kill (looks like EOF)"; st=1; }

    log "== inputs =="
    if check_inputs "${r1}" "${r2}" "${asm}"; then
        log "  reads=$(human_bytes "${IN_READS_BYTES}") assembly=$(human_bytes "${IN_ASM_BYTES}") records=${IN_RECORDS}"
        need=$(est_footprint_bytes); avail=$(( $(df_avail_gb "${LOCAL_BASE}") * 1024 * 1024 * 1024 ))
        log "  estimated footprint=$(human_bytes "${need}") (heuristic: ${EST_READS_FACTOR}x reads + ${EST_ASM_FACTOR}x assembly + ${EST_FLAT_MB}MB) avail=$(human_bytes "${avail}")"
        (( need < avail )) || { log "  ERROR: estimated footprint does not fit"; st=1; }
    else
        log "  ERROR: ${INPUT_ERR}"; st=1
    fi

    log "== installed modules =="
    audit_module_tmpdirs
    if mw_patched; then log "  L3 patch (MW_TMPDIR): applied"; else log "  L3 patch (MW_TMPDIR): not applied"; fi
    (( st == 0 )) && log "CHECK: OK" || log "CHECK: FAILED"
    return "${st}"
}

# ------- batch feasibility: "will a big batch work?" without running metaWRAP --
plan_fields() {   # plan_fields <line> -> R1 R2 ASM ; sets PLAN_NF to the column count
    local line=$1 nf=0
    line=${line%$'\r'}                       # CRLF-tolerant
    [[ -n "${line}" ]] && nf=$(printf '%s' "${line}" | awk -F'\t' '{print NF}')
    PLAN_NF=${nf:-0}
    # awk (not cut) so a line WITHOUT a TAB yields EMPTY fields -> defaults apply;
    # `cut -f2` would return the whole line and silently use the sample id as a path.
    R1=$(printf '%s' "${line}" | awk -F'\t' '{print $2}')
    R2=$(printf '%s' "${line}" | awk -F'\t' '{print $3}')
    ASM=$(printf '%s' "${line}" | awk -F'\t' '{print $4}')
    R1=${R1%$'\r'}; R2=${R2%$'\r'}; ASM=${ASM%$'\r'}
    # 1 column = sample id only -> all defaults ; 4 columns = explicit paths ;
    # anything else (2/3/>4) is ambiguous and must be reported by the caller,
    # never silently back-filled with a default assembly.
    if (( PLAN_NF == 1 )); then R1=""; R2=""; ASM=""; fi
    [[ -n "${R1}" ]]  || R1="${READS_DIR}/${SAMPLE_ID}_1.fastq"
    [[ -n "${R2}" ]]  || R2="${READS_DIR}/${SAMPLE_ID}_2.fastq"
    [[ -n "${ASM}" ]] || ASM="${ASSEMBLY_DIR}/final.contigs.fa"
}
mode_batch_dry() {
    local plan=$1 line sid sock cap need avail_gb memok coll ok rc_inputs wfs reasons rh
    local seen_work=" " seen_final=" " n=0 bad=0 tdir insitu="OK" probe_fail=0
    local -a sids=()
    [[ -f "${plan}" ]] || die "plan file not found: ${plan}"
    while IFS= read -r line; do
        line=${line%$'\r'}
        [[ -z "${line}" || "${line}" == \#* ]] && continue
        sids+=("$(printf '%s' "${line}" | cut -f1)")
    done < "${plan}"
    (( ${#sids[@]} > 0 )) || die "no samples in ${plan}"

    NODE_MEM_GB=$(detect_mem_limit_gb); (( NODE_MEM_GB > 0 )) || NODE_MEM_GB=40
    (( REFINE_MEM_GB <= 0 )) && REFINE_MEM_GB="${NODE_MEM_GB}"
    if (( REASSEMBLE_MEM_GB <= 0 )); then
        if (( NODE_MEM_GB > MEM_HEADROOM_GB )); then REASSEMBLE_MEM_GB=$(( NODE_MEM_GB - MEM_HEADROOM_GB )); else REASSEMBLE_MEM_GB="${NODE_MEM_GB}"; fi
    fi
    wfs=$(fs_type_of "${LOCAL_BASE}"); cap=$(derived_work_cap)
    : > "${OUT_REPORT}" 2>/dev/null || warn "cannot write report to ${OUT_REPORT}"
    log "batch-dry: ${#sids[@]} sample(s) | LOCAL_BASE=${LOCAL_BASE} (${wfs:-unknown}) | usable memory ${NODE_MEM_GB} GB | reassemble -m ${REASSEMBLE_MEM_GB} GB"
    log "derived caps: len(SHORT_WORK) <= ${cap}; reassemble socket <= ${SUN_LIMIT} (measured suffix ${SOCK_SUFFIX})"
    row "$(printf '%-14s %5s %5s %5s %8s %7s %5s %5s %10s %s' SAMPLE wlen sock fs availGB needGB memOK inOK collision VERDICT)"

    for sid in "${sids[@]}"; do
        plan_sample "${sid}"
        sock=$(reassemble_socket_len)
        plan_fields "$(awk -F'\t' -v s="${sid}" '$1==s{print;exit}' "${plan}" 2>/dev/null || true)"
        need=0; rc_inputs=OK; INPUT_ERR=""
        if (( PLAN_NF != 1 && PLAN_NF != 4 )); then
            rc_inputs=BAD; INPUT_ERR="plan line for ${sid} has ${PLAN_NF} columns (need 1 or 4)"
        elif ! check_inputs "${R1}" "${R2}" "${ASM}"; then
            rc_inputs=BAD
        fi
        [[ "${rc_inputs}" == OK ]] && need=$(est_footprint_bytes)
        avail_gb=$(df_avail_gb "${LOCAL_BASE}")
        memok=OK; (( REASSEMBLE_MEM_GB > NODE_MEM_GB )) && memok=BAD
        coll=OK
        case "${seen_work} "  in *" ${SHORT_WORK} "*)  coll=DUP_WORK ;; esac
        case "${seen_final} " in *" ${FINAL_ROOT} "*) coll="${coll}+DUP_FINAL" ;; esac
        is_network_fs "${wfs}" && coll="${coll}+NFS"
        seen_work="${seen_work}${SHORT_WORK} "; seen_final="${seen_final}${FINAL_ROOT} "

        # Accumulate ALL failing reasons (never overwrite): a later check must not
        # hide an earlier, more specific one.
        reasons=""
        (( sock <= SUN_LIMIT ))     || reasons="${reasons}+SOCK_TOO_LONG"
        (( ${#SHORT_WORK} <= cap )) || reasons="${reasons}+PATH_TOO_LONG"
        [[ "${rc_inputs}" == OK ]]  || reasons="${reasons}+INPUTS(${INPUT_ERR})"
        [[ "${coll}" == OK ]]       || reasons="${reasons}+COLLISION(${coll})"
        [[ "${memok}" == OK ]]      || reasons="${reasons}+MEM"
        (( need < avail_gb * 1024 * 1024 * 1024 )) || reasons="${reasons}+SPACE"

        # In-situ rehearsal in a SEPARATE tree of identical path length (r_ vs w_),
        # so it never risks a pinned real work dir and still exercises the real FS,
        # permissions and path length. BIN_REFINEMENT is probed as an ABSOLUTE path
        # -- a conservative upper bound on metaWRAP's relative `binsA.tmp`.
        insitu=OK
        rh="${LOCAL_BASE}/r_${SAMPLE_ID}_$$"
        if mkdir -p "${rh}/01_BIN" "${rh}/02_REF" "${rh}/03_REA/tmp" "${rh}/mt" "${rh}/log" 2>/dev/null; then
            for tdir in "${rh}/mt" "${rh}/03_REA/tmp" "${rh}/02_REF"; do
                if ! "${PYBIN}" "${PROBE_PY}" "${tdir}" >/dev/null 2>&1; then
                    insitu="FAIL"; probe_fail=$(( probe_fail + 1 ))
                    log "  ${sid}: in-situ AF_UNIX bind FAILED at ${tdir} -- the real run WOULD die at this step"
                fi
            done
            printf 'rehearsed %s (insitu=%s)\n' "${rh}" "${insitu}" >> "${OUT_REPORT}" 2>/dev/null || true
            rm -rf "${rh}"
        else
            insitu="FAIL"; probe_fail=$(( probe_fail + 1 ))
            log "  ${sid}: cannot create the rehearsal tree under ${LOCAL_BASE}"
        fi
        [[ "${insitu}" == OK ]] || reasons="${reasons}+INSITU_FAIL"

        ok="${reasons#+}"; [[ -n "${ok}" ]] || ok=OK
        [[ "${ok}" == OK ]] || bad=$(( bad + 1 ))
        row "$(printf '%-14s %5d %5d %5s %8d %7d %5s %5s %10s %s' \
            "${sid}" "${#SHORT_WORK}" "${sock}" "${wfs:-?}" "${avail_gb}" \
            "$(( need / 1024 / 1024 / 1024 ))" "${memok}" "${rc_inputs}" "${coll}" "${ok}")"
        n=$(( n + 1 ))
    done

    log "batch-dry complete: ${n} sample(s); ${bad} with a non-OK verdict; ${probe_fail} in-situ socket failure(s)"
    (( bad == 0 )) && log "BATCH-DRY: OK" || log "BATCH-DRY: FAILED"
    log "report: ${OUT_REPORT}"
    if (( bad > 0 )); then return 1; fi
    return 0
}

mode_batch_run() {
    local plan=$1 line sid r1 r2 asm
    local -a oklist=() badlist=()
    local i=0
    log "batch-run: running the read-only preflight first"
    mode_batch_dry "${plan}" || die "batch preflight failed; refusing to start the batch"
    while IFS= read -r line; do
        [[ -z "${line}" || "${line}" == \#* ]] && continue
        line=${line%$'\r'}
        sid=$(printf '%s' "${line}" | cut -f1); i=$(( i + 1 ))
        plan_sample "${sid}"; plan_fields "${line}"
        r1="${R1}"; r2="${R2}"; asm="${ASM}"
        log "=== [${i}] sample ${sid} ==="
        if (( PLAN_NF != 1 && PLAN_NF != 4 )); then
            badlist+=("${sid}"); warn "sample ${sid}: plan line has ${PLAN_NF} columns (need 1 or 4); skipping"
            continue
        fi
        if ! check_inputs "${r1}" "${r2}" "${asm}"; then
            badlist+=("${sid}"); warn "sample ${sid}: ${INPUT_ERR}"
            continue
        fi
        if (( ARCHIVE_EXISTING == 1 )) && [[ -e "${FINAL_ROOT}" ]]; then
            mv "${FINAL_ROOT}" "${FINAL_ROOT}_backup_$(date +%Y%m%d_%H%M%S)"
            log "archived existing ${FINAL_ROOT}"
        fi
        # Subshell: a `die`/exit inside run_sample must fail ONLY this sample so the
        # batch keeps going (input errors above already continue).
        if ( run_sample "${r1}" "${r2}" "${asm}" ); then oklist+=("${sid}")
        else badlist+=("${sid}"); warn "sample ${sid} FAILED (continuing with the next sample)"; fi
    done < "${plan}"
    log "=== batch summary ==="
    log "  succeeded (${#oklist[@]}): ${oklist[*]:-none}"
    log "  failed    (${#badlist[@]}): ${badlist[*]:-none}"
    if (( ${#badlist[@]} == 0 )); then log "BATCH-RUN: OK"; return 0; fi
    log "BATCH-RUN: FAILED"
    return 1
}

mode_postrun() {
    [[ -n "${SAMPLE_IDS}" ]] || die "postrun needs --samples <sample>"
    local sid; sid=$(printf '%s' "${SAMPLE_IDS}" | cut -d, -f1)
    if [[ -z "${MW_WORK_DIR:-}" ]]; then
        local newest
        newest=$(ls -1dt "${LOCAL_BASE}/w_${sid}_"* 2>/dev/null | head -n1 || true)
        if [[ -n "${newest}" ]]; then
            MW_WORK_DIR="${newest}"; log "postrun: using the existing work dir ${MW_WORK_DIR}"
        else
            warn "no ${LOCAL_BASE}/w_${sid}_* work dir found; falling back to a fresh plan (logs will be missing)"
        fi
    fi
    plan_sample "${sid}"
    local st=0 lg="${LOGDIR_SHORT}/03_reassemble_bins.log" nfa
    log "== criterion 1: fatal CheckM signatures (the cosmetic ones are advisory) =="
    if [[ ! -s "${lg}" ]]; then log "  ERROR: log missing/empty: ${lg}"; st=1
    elif grep -nE 'AF_UNIX path too long|Unexpected error:.*EOFError|Something went wrong with running CheckM' "${lg}"; then
        log "  ERROR: fatal signature(s) above"; st=1
    else log "  OK (none)"; fi
    if [[ -s "${lg}" ]] && grep -qE "${ADVISORY_RE}" "${lg}"; then
        log "  ADVISORY: .nfs/EBUSY present -- non-fatal by design (util.py:275-280 catches and prints only)"
    fi
    log "== criterion 2: SUCCESS banner =="
    if grep -q 'REASSEMBLY PIPELINE SUCCESSFULLY FINISHED' "${lg}" 2>/dev/null; then log "  OK"
    else log "  ERROR: banner missing"; st=1; fi
    log "== criterion 3: outputs =="
    nfa=$(find "${REASSEMBLED_BINS}" -maxdepth 1 -type f -name '*.fa' 2>/dev/null | wc -l || true)
    if [[ -s "${REASSEMBLED_STATS}" ]]; then log "  OK ${REASSEMBLED_STATS}"
    else log "  ERROR: missing/empty ${REASSEMBLED_STATS}"; st=1; fi
    if [[ -s "${REASSEMBLED_CHECKM_TSV}" ]]; then log "  OK ${REASSEMBLED_CHECKM_TSV}"
    else log "  ERROR: missing/empty ${REASSEMBLED_CHECKM_TSV}"; st=1; fi
    log "  reassembled .fa: ${nfa}"
    (( nfa > 0 )) || { log "  ERROR: no reassembled .fa"; st=1; }
    (( st == 0 )) && log "POSTRUN: OK" || log "POSTRUN: FAILED"
    return "${st}"
}

mode_verify_l3() {
    local p n
    p=$(find_mw_module reassemble_bins.sh || true)
    [[ -n "${p}" ]] || die "reassemble_bins.sh not found (is 'metawrap' on PATH?)"
    n=$(grep -c 'MW_TMPDIR' "${p}" 2>/dev/null || true)
    log "module        : ${p}"
    log "MW_TMPDIR hits: ${n} (expect 4 = 2x mkdir + 2x --tmpdir)"
    grep -n -- '--tmpdir\|MW_TMPDIR\|mkdir.*tmp' "${p}" || true
    if (( n >= 4 )); then log "STATUS: patched"; return 0; fi
    log "STATUS: unpatched"; return 1
}
mode_patch_l3() {
    local p; p=$(find_mw_module reassemble_bins.sh || true)
    [[ -n "${p}" ]] || die "reassemble_bins.sh not found"
    if grep -q 'MW_TMPDIR' "${p}"; then log "already patched"; mode_verify_l3; return 0; fi
    cp -a "${p}" "${p}.mw_tmpdir.bak"
    sed -i 's|--tmpdir ${out}/tmp|--tmpdir ${MW_TMPDIR:-${out}/tmp}|g' "${p}"
    sed -i 's|mkdir ${out}/tmp$|mkdir -p ${MW_TMPDIR:-${out}/tmp}|g' "${p}"
    if ! bash -n "${p}"; then cp -a "${p}.mw_tmpdir.bak" "${p}"; die "patched module fails bash -n; restored"; fi
    mode_verify_l3
    log "NOTE: scope is reassemble_bins.sh only. binning.sh:62 also passes an ABSOLUTE --tmpdir when invoked with --run-checkm, which this guard never uses."
}
mode_unpatch_l3() {
    local p; p=$(find_mw_module reassemble_bins.sh || true)
    [[ -n "${p}" ]] || die "reassemble_bins.sh not found"
    [[ -f "${p}.mw_tmpdir.bak" ]] || die "no backup at ${p}.mw_tmpdir.bak"
    cp -a "${p}.mw_tmpdir.bak" "${p}"
    log "reverted ${p}"
}

mode_run() {
    plan_sample "${SAMPLE_IDS:-ERR011347}"
    # Explicit inputs (--r1/--r2/--assembly) let a driver such as Snakemake pass
    # arbitrary read/assembly paths; without them the historical defaults apply.
    # Reads may be .fq or .fq.gz; the assembly may be .fa or .fa.gz.
    local r1="${RUN_R1:-${READS_DIR}/${SAMPLE_ID}_1.fastq}"
    local r2="${RUN_R2:-${READS_DIR}/${SAMPLE_ID}_2.fastq}"
    local asm="${RUN_ASM:-${ASSEMBLY_DIR}/final.contigs.fa}" c
    WORK_FS=$(fs_type_of "${LOCAL_BASE}")
    if is_network_fs "${WORK_FS}" && (( ALLOW_NETWORK_WORK == 0 )); then
        die "LOCAL_BASE ${LOCAL_BASE} is on '${WORK_FS}' (network FS). CheckM writes temp inside -o. Move it or set ALLOW_NETWORK_WORK=1."
    fi
    check_inputs "${r1}" "${r2}" "${asm}" || die "${INPUT_ERR}"
    activate_conda_env || die "could not activate ${METAWRAP_ENV}"
    PYBIN="${CONDA_PREFIX}/bin/python2.7"
    if [[ ! -x "${PYBIN}" ]]; then
        warn "no ${CONDA_PREFIX}/bin/python2.7; falling back to 'python' (the suffix is measured, not assumed)"
        PYBIN="$(command -v python || command -v python3)"
    fi
    calibrate
    log "PYBIN=${PYBIN}"
    log "MW_TMPDIR=${MW_TMPDIR} (takes effect only when the L3 patch is applied)"
    while IFS= read -r c; do
        command -v "${c}" >/dev/null 2>&1 || die "missing executable: ${c}"
    done < <(required_tools)
    log "module temp-dir audit:"; audit_module_tmpdirs
    if (( ARCHIVE_EXISTING == 1 )) && [[ -e "${FINAL_ROOT}" ]]; then
        mv "${FINAL_ROOT}" "${FINAL_ROOT}_backup_$(date +%Y%m%d_%H%M%S)"; log "archived existing ${FINAL_ROOT}"
    fi
    NODE_MEM_GB=$(detect_mem_limit_gb); (( NODE_MEM_GB > 0 )) || NODE_MEM_GB=40
    (( REFINE_MEM_GB <= 0 )) && REFINE_MEM_GB="${NODE_MEM_GB}"
    if (( REASSEMBLE_MEM_GB <= 0 )); then
        if (( NODE_MEM_GB > MEM_HEADROOM_GB )); then REASSEMBLE_MEM_GB=$(( NODE_MEM_GB - MEM_HEADROOM_GB )); else REASSEMBLE_MEM_GB="${NODE_MEM_GB}"; fi
    fi
    run_sample "${r1}" "${r2}" "${asm}"
}

usage() {
    cat <<'EOF'
mw_guard.sh <mode> [options]

  selftest                     calibrate this box + reason-checked probe pair
  check                        single-sample preflight
  batch-dry --plan FILE        batch feasibility, runs no metaWRAP
  batch-run --plan FILE        preflight, then run every sample in FILE
  run                          run one sample
  postrun --samples ID         acceptance for a finished sample
  verify-l3 | patch-l3 | unpatch-l3

Options
  --plan FILE      TSV: sample_id [TAB r1] [TAB r2] [TAB assembly]; '#' comments ok
                   (columns 2-4 optional; 1 column = defaults; 2/3 columns = error)
  --samples LIST   comma-separated sample ids (check/run/postrun)
  --work DIR       pin the active tree to an existing run (postrun/check)
  --r1 FILE        explicit forward reads for 'run' (overrides READS_DIR default)
  --r2 FILE        explicit reverse reads for 'run' (overrides READS_DIR default)
  --assembly FILE  explicit assembly for 'run' (overrides ASSEMBLY_DIR default)
                   reads may be .fq or .fq.gz; assembly .fa or .fa.gz

Env: METAWRAP_ENV THREADS READS_DIR ASSEMBLY_DIR BASE_FINAL MW_LOCAL_BASE
     REFINE_MEM_GB REASSEMBLE_MEM_GB INPUT_STRICT RESUME ALLOW_NETWORK_WORK
     INSTALL_RMTREE_PATCH EST_READS_FACTOR EST_ASM_FACTOR EST_FLAT_MB OUT_REPORT

NOTE on batching: batch-run is SEQUENTIAL on purpose (predictable memory, no
shared temp). For N-way parallelism use the scheduler -- sbatch --array, or
xargs -P -- with per task MW_LOCAL_BASE=$SLURM_TMPDIR (or another node-local
dir) so no two samples share a temp dir, then re-run with RESUME=1.
EOF
}

main() {
    local mode="${1:-}"; shift || true
    local plan=""
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --plan)    plan=$2; shift 2 ;;
            --samples) SAMPLE_IDS=$2; shift 2 ;;
            --work)    MW_WORK_DIR=$2; shift 2 ;;
            --r1)      RUN_R1=$2; shift 2 ;;
            --r2)      RUN_R2=$2; shift 2 ;;
            --assembly) RUN_ASM=$2; shift 2 ;;
            -h|--help) usage; exit 0 ;;
            *) die "unknown option: $1" ;;
        esac
    done
    # NOTE: every mode that can legitimately report failure (check / batch-dry /
    # batch-run / postrun / verify-l3) is called as `... || rc=$?` and the code is
    # returned via `exit "${rc}"`. Without this, `set -e` plus the ERR trap turned a
    # normal "PREFLIGHT FAILED" verdict into a bogus
    # "ERROR: aborted at <file> line N" crash message.
    local rc=0
    case "${mode}" in
        selftest)   mode_selftest || rc=$? ;;
        check)      calibrate; mode_check || rc=$? ;;
        batch-dry)  [[ -n "${plan}" ]] || die "batch-dry needs --plan FILE"
                    calibrate; mode_batch_dry "${plan}" || rc=$? ;;
        batch-run)  [[ -n "${plan}" ]] || die "batch-run needs --plan FILE"
                    activate_conda_env || die "could not activate ${METAWRAP_ENV}"
                    PYBIN="${CONDA_PREFIX}/bin/python2.7"
                    [[ -x "${PYBIN}" ]] || PYBIN="$(command -v python || command -v python3)"
                    calibrate; mode_batch_run "${plan}" || rc=$? ;;
        run)        mode_run ;;
        postrun)    mode_postrun || rc=$? ;;
        verify-l3)  mode_verify_l3 || rc=$? ;;
        patch-l3)   mode_patch_l3 ;;
        unpatch-l3) mode_unpatch_l3 ;;
        -h|--help|"") usage; if [[ -z "${mode}" ]]; then exit 2; fi; exit 0 ;;
        *) die "unknown mode: ${mode}" ;;
    esac
    exit "${rc}"
}
main "$@"
