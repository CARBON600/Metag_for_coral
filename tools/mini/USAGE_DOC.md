# mini fixture: end-to-end regression check for Metag_for_coral-main

**In one sentence**: after changing a line of code, prove with a single command on
a real machine that the whole pipeline still works, without reading the code by
eye. The fixture is simulated paired reads plus a partial config override, run in
an **isolated clone**, and never touches your real `output/`. The scripts ship in
the repository under `tools/mini/`, so there is nothing to copy.

---

## 0. Full command sequence (copy-paste)

```bash
# -- Prerequisites ------------------------------------------------
conda activate snakemake
REPO=/mnt/storage5/yfdai/for_nan/SCGC/pipeline/Metag_for_coral-main
GUARD="$REPO/tools/mw_guard.sh"
LOCAL_BASE=/tmp/yfdai/mw          # short, node-local; capacity >= a single sample's needGB
ls -l "$GUARD" && md5sum "$GUARD" # expect 1.1 / 50e6127e...

# -- (1) the scripts ship with the repo (tools/mini/); nothing to copy ----
#   tools/mini/{make_mini_fixture.py,validate_run.sh,check_outputs.py,README.md}

# -- (2) first run: build the fixture (once; it lives in the clone and is skipped on every sync) --
#   recommended: use two of your own high-quality MAGs so CheckM/GTDB are "real"
bash tools/mini/validate_run.sh --repo "$REPO" --dry-only   # first time: reports "no fixture yet"
CLONE="${REPO%/}.minitest"
cd "$CLONE"
python3 tools/mini/make_mini_fixture.py --out tools/mini/fixture --repo "$PWD" \
  --genome "$REPO"/output/drep_cross_workflow/dereplicated_genomes/<MAG_A>.fna \
  --genome "$REPO"/output/drep_cross_workflow/dereplicated_genomes/<MAG_B>.fna \
  --guard "$GUARD" --local-base "$LOCAL_BASE"
#   with no MAGs yet, use dependency-free synthetic genomes (plumbing only; 0 MAGs pass the gate):
#   python3 tools/mini/make_mini_fixture.py --out tools/mini/fixture --repo "$PWD" --synthesize \
#       --guard "$GUARD" --local-base "$LOCAL_BASE"

# -- (3) see the DAG and pin the baseline fingerprint (seconds) ----
cd "$REPO"                         # go back to the repo root before running validate_run
bash tools/mini/validate_run.sh --repo "$REPO" --dry-only --bless

# -- (4) first full run (30-60 min; SPAdes, metaWRAP per-bin reassembly, GTDB-Tk) --
bash tools/mini/validate_run.sh --repo "$REPO" --cores 16

# -- (5) read the result ------------------------------------------
#   the terminal should end with: [check_outputs] NN PASS / 0 WARN / 0 FAIL then VALIDATE: OK
#   on failure:
grep -nE "Error in rule|FAIL" "${CLONE}/logs/validate_run.log" "${CLONE}/logs/dry.log"
python3 "$CLONE/tools/mini/check_outputs.py" outputs --output-dir "$CLONE/output" \
        --config "$CLONE/tools/mini/fixture/config.mini.yaml"

# -- (6) after every code change, this one command -----------------
bash tools/mini/validate_run.sh --repo "$REPO" --cores 16
#   to fill gaps only (fast, but rule-code changes do not trigger a recompute; run.sh defaults to mtime):
bash tools/mini/validate_run.sh --repo "$REPO" --incremental
#   to confirm only that the DAG is intact (seconds):
bash tools/mini/validate_run.sh --repo "$REPO" --dry-only

# -- (7) if you changed the matrix/rule expectations the fingerprint will DIFF; after review, re-pin it --
bash tools/mini/validate_run.sh --repo "$REPO" --dry-only --bless
```

---

## 1. Files in this directory

| File | Purpose |
|---|---|
| `make_mini_fixture.py` | Build the fixture: simulate paired 2×150 bp reads from 1 to 2 genomes (seeded, reproducible) and write `data/`, `genomes/`, `fixture_manifest.tsv`, `config.mini.yaml` (a partial override that changes only data_dir and the matrix/threads; everything else comes from the repository `config.yaml`) |
| `validate_run.sh` | One command: sync to an isolated clone, dry-run for the DAG fingerprint, forced full run, call the asserter, report PASS/FAIL |
| `check_outputs.py` | Asserter (`outputs` mode, 31 checks) plus DAG fingerprint (`fingerprint` mode, `--bless`/`--check`) |
| `README.md` | This document (coverage boundaries and troubleshooting table) |

## 2. How isolation works

`validate_run.sh` syncs the repository into **`<repo>.minitest`** (a sibling
directory by default; change it with `--clone`) and excludes
`output/`, `.snakemake/`, `logs/`, `harvest/`, `tools/mini/`:

* **code** flows in on every sync, so a rule change takes effect immediately;
* **results** stay in the clone, so incremental runs are fast;
* the **fixture** lives in `tools/mini/fixture/` and is never deleted by a sync;
* the real repository's `output/` (especially the **global** tables
  `output/summary/mimag_*.tsv`, `drep_cross_workflow/`) is **completely
  untouched**; this is why the isolated clone is required.

## 3. How the matrix is narrowed

`1 sample (MINI1) × 1 group (exp) × 1 method (bt2) × 1 treat (control_assemble) × 2 binners`
= **2 combinations, 28 jobs** (measured). `--methods bt2,fastqs` gives 48, additionally
covering `fastp_clean`/`fastq_screen`.

Every rule runs at least once:

* preprocessing and assembly: `bowtie2_map → bowtie2_unmapped → bam_to_fastq_bt2 → spades → remap_reads_to_final_contigs → filter_contigs_r2000`
* native arm: `prepare_bins → refinem_bins → checkm_lineage_wf → checkm_qa → gtdbtk_classify → drep_per_workflow`
* metaWRAP arm: `metawrap_bins(guard) → prepare_bins → metawrap_ingest → checkm_lineage_wf → checkm_qa → gtdbtk → drep`
* summary layer: `mimag_gate → drep_cross_workflow → drep_taxonomy`

## 4. What the asserter checks

Each stage marker; bin basenames equal the CheckM keys; metawrap did not run
RefineM (`scaffold_stats.tsv` must not exist); the reported CheckM table is
**not** a copy of the guard's 1.0.12 table (different md5); `results.tsv` carries
the `Strain heterogeneity`/`Marker lineage` headers; the `strain_heterogeneity`
column is **really non-empty**; `mimag_summary` covers every expected combination
with a valid status; `per_workflow.tsv` has a metawrap row; dRep's
`Bdb/Cdb/Wdb.csv` and `final_mags.tsv` all exist.

## 5. What is not covered (do not mistake it for full green)

1. A real community (unless `--genome` passes real MAGs): under `--synthesize`, 0 MAGs pass the MIMAG gate.
2. The `coverm` branch, the `PCR_assemble` branch, the second group's database paths.
3. Scale: real-sample footprint, peak memory, runtime, concurrent resource contention.
4. `mw_guard.sh` `batch-*` modes, the L3 patch, `KEEP_SHORT_WORK` cleanup.
5. `checkm rescored` is an **md5 heuristic** (it proves "not a byte-for-byte copy", not the version).

## 6. Troubleshooting table (these 6 defects are known to be caught)

| FAIL | Meaning | Where to look |
|---|---|---|
| `key==filename …` | bin basenames in refinem are decoupled from CheckM's first column | rename logic, `BINNER_LAYOUT` vs `common.smk:binner_extension` |
| `refinem skipped metawrap …` | RefineM was scheduled for metawrap | `wildcard_constraints` of `rules/binning.smk::refinem_bins` |
| `checkm rescored metawrap …` | the guard's 1.0.12 table was used as the result | whether `checkm_lineage_wf` still covers metawrap; whether `metawrap_ingest` wrote an extra `checkm.done` |
| `strain_heterogeneity populated` | the `results.tsv` path is broken (column/key mismatch) | `rules/binning.smk::checkm_qa`, `harvest_mags._load_checkm_qa` |
| `checkm_qa … missing results.tsv` | the new rule did not run, or its output path changed | the `checkm_qa` output declaration |
| `prepare_bins … missing prepare.done` | the arm was skipped, so `combo_status` reports `qc_empty` | `common.smk::binner_done` |

## 7. Appendix: how to check the two values (`config.yaml`'s `metawrap.guard` and `local_base`)

```bash
# guard: relative path shipped with the repo, plus its fingerprint
cd "$REPO"
md5sum tools/mw_guard.sh                                   # expect 1.1 / 50e6127e...

# local_base: candidate check (filesystem type, path length, free space)
for d in "/tmp/$USER/mw" /mnt/storage3/yfdai/mw /mnt/storage5/yfdai/mw; do
  echo "$d  len=${#d}  fs=$(stat -f -c %T "$(dirname "$d")" 2>/dev/null)"
done
df -hT /tmp /mnt/storage3 /mnt/storage5
MW_LOCAL_BASE="$LOCAL_BASE" bash tools/mw_guard.sh check --samples ERR011347
```

Verdict: a **local** filesystem (`xfs`/`ext4`/`tmpfs`; `nfs`/`lustre`/`beegfs` are
rejected by the guard), a path as short as possible (the AF_UNIX budget is 107
bytes), and free space at least a single sample's `needGB` (`/tmp` on `/` has only
about 47 GB; measure one real sample before batching).

## 8. Verified and not verified in the sandbox

**Verified**: generator self-checks (FASTQ 4-line structure, seq/qual equal
length, pure ACGT, gz, `/1`/`/2` pairing, reads locate in their genome, pairing
geometry); the mini config really narrows the DAG to 28/48 jobs (Snakemake 7.32.4
dry-run); DAG fingerprint bless/check and tamper detection; the asserter reports
31 PASS/0 FAIL on a correct tree and catches the 6 injected defects above.

**Not verified** (must run on your machine): the real pipeline end to end. The
sandbox has none of the nine conda envs and no coral/GTDB/CheckM databases. Allow
30 to 60 minutes for the first run and check `logs/validate_run.log` first.
