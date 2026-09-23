# mini fixture: end-to-end regression check after every code change

A tiny fixture plus an assertion suite for `Metag_for_coral-main`. The goal: after
changing one line of code, prove with a single command on a real machine that the
whole pipeline still works, instead of reading the code by eye.

This directory ships with the repository (`tools/mini/`), so there is nothing to
copy elsewhere.

## Files

| File | Purpose |
|---|---|
| `make_mini_fixture.py` | Build the fixture: simulate paired reads from 1 to 2 genomes and write `data/`, `genomes/`, `fixture_manifest.tsv`, `config.mini.yaml` (a partial override that changes only data_dir and the matrix/threads; everything else comes from the repository `config.yaml`) |
| `validate_run.sh` | One command: sync to an isolated clone, dry-run for the DAG fingerprint, run the full pipeline at forced recompute, then call the asserter |
| `check_outputs.py` | Asserter (`outputs` mode) plus DAG fingerprint (`fingerprint` mode, `--bless`/`--check`) |
| `README.md` | This file |

## One-time setup (on the cluster, after `conda activate snakemake`)

```bash
# From the repository root (or anywhere, using --repo to point at the repo):
bash tools/mini/validate_run.sh --repo <your repos>/Metag_for_coral-main --dry-only
# The first time it reports "no fixture yet" and prints the command to create one:

cd <your repos>/Metag_for_coral-main.minitest
python3 tools/mini/make_mini_fixture.py --out tools/mini/fixture --repo "$PWD" \
    --synthesize \
    --guard "$PWD/tools/mw_guard.sh" \
    --local-base /tmp/$USER/mw
```

The fixture lives in the **clone** (`tools/mini/fixture/`); `validate_run.sh`
skips `tools/mini/` on every sync, so it is created only once.

For a fixture that can actually pass the MIMAG gate (recommended; a much stronger
test), pass two of your own high-quality MAGs:

```bash
python3 tools/mini/make_mini_fixture.py --out tools/mini/fixture --repo "$PWD" \
    --genome <MAG_A>.fna --genome <MAG_B>.fna \
    --guard "$PWD/tools/mw_guard.sh" --local-base /tmp/$USER/mw
```
MAG source: `<repo>/output/drep_cross_workflow/dereplicated_genomes/*.fna` (pick
two different species). `--synthesize` builds random genomes: perfectly valid
format, assembles and bins, but CheckM scores about 0% completeness, so **no bin
can pass the MIMAG gate**. That fixture tests only that the plumbing runs, not
the QC-passing path.

## After every code change

```bash
bash tools/mini/validate_run.sh --repo <your repos>/Metag_for_coral-main   # default --cores 16, forced full run
# Common variants:
#   --incremental   fill gaps only (fast; but rule-code changes do not trigger a recompute, run.sh defaults to mtime)
#   --dry-only      dry-run and DAG fingerprint only (seconds)
#   --bless         rewrite the DAG fingerprint baseline after changing the matrix/rule expectations
```

Normal output looks like:

```
== [3/6] dry run (DAG fingerprint) ==
[fingerprint] MATCH (expected total=28, got 28)
== [5/6] assertions ==
PASS  key==filename metawrap/control_assemble/bt2/exp/MINI1    2 id(s)
...
[check_outputs] 31 PASS / 0 WARN / 0 FAIL
VALIDATE: OK
```

## What this fixture covers

The matrix is narrowed to `1 sample × 1 group × 1 method(bt2) × 1 treat(control_assemble) × 2 binners
(semibin2_single + metawrap)` = **2 combinations, 28 jobs** (`--methods bt2,fastqs` gives 48,
also covering `fastp_clean`/`fastq_screen`). Every rule runs at least once, including:

* preprocessing and assembly: `bowtie2_map → bowtie2_unmapped → bam_to_fastq_bt2 → spades → remap → filter_contigs_r2000`
* native arm: `prepare_bins → refinem_bins → checkm_lineage_wf → checkm_qa → gtdbtk_classify → drep_per_workflow`
* metaWRAP arm: `metawrap_bins(guard) → prepare_bins → metawrap_ingest → checkm_lineage_wf → checkm_qa → gtdbtk → drep`
* summary layer: `mimag_gate → drep_cross_workflow → drep_taxonomy`

The asserter checks each item (a FAIL prints the concrete path and value): every
stage marker; bin basenames equal the CheckM keys; metawrap did not run RefineM
(`scaffold_stats.tsv` must not exist); the reported CheckM table is not a copy of
the guard's 1.0.12 table (different md5); `results.tsv` carries the
`Strain heterogeneity`/`Marker lineage` headers; the `strain_heterogeneity`
column is really non-empty; `mimag_summary` covers every expected combination
with a valid status; `per_workflow.tsv` has a metawrap row; both dRep products
exist.

## What this fixture does NOT cover (do not mistake it for full green)

1. A real community (unless you pass real MAGs with `--genome`): under `--synthesize`, 0 MAGs pass the gate.
2. The `coverm` branch, the `PCR_assemble` branch, and the second `group`'s database paths.
3. Anything scale-related: real-sample footprint, peak memory, runtime, concurrent resource contention.
4. `mw_guard.sh` `batch-*` modes, the L3 patch, `KEEP_SHORT_WORK` cleanup.
5. `checkm rescored` is an **md5 heuristic**: it shows "not a byte-for-byte copy", not the version.

## Troubleshooting table (these 6 defects are known to be caught)

| FAIL | Meaning | Where to look |
|---|---|---|
| `key==filename …` | bin basenames in refinem are decoupled from CheckM's first column | rename logic, `BINNER_LAYOUT` vs `common.smk:binner_extension` |
| `refinem skipped metawrap …` | RefineM was scheduled for metawrap | `wildcard_constraints` of `rules/binning.smk::refinem_bins` |
| `checkm rescored metawrap …` | the guard's 1.0.12 table was used as the result | whether `checkm_lineage_wf` still covers metawrap; whether `metawrap_ingest` wrote an extra `checkm.done` |
| `strain_heterogeneity populated` | the `results.tsv` path is broken (header/column/key mismatch) | `rules/binning.smk::checkm_qa`, `harvest_mags._load_checkm_qa` |
| `checkm_qa … missing results.tsv` | the new rule did not run, or its output path changed | the `checkm_qa` output declaration |
| `prepare_bins … missing prepare.done` | the arm was skipped, so `combo_status` reports `qc_empty` | `common.smk::binner_done` |

## Verified and not verified in the sandbox

**Verified**: generator self-checks (FASTQ 4-line structure, equal lengths, ACGT,
gz, pairing geometry, reads locate in their genome); the mini config really
narrows the DAG to 28/48 jobs (Snakemake 7.32.4 dry-run); DAG fingerprint
bless/check and tamper detection; the asserter reports 31 PASS/0 FAIL on a
correct tree and catches the 6 injected defects above.

**Not verified** (must run on your machine): the real pipeline end to end. The
sandbox has none of the nine conda envs and no coral/GTDB/CheckM databases.

## Measured cost and rerun strategy (2026-09-23 on cl007, 28-job fixture, real MAGs: bin.1.strict + bin.10.permissive)

Total **1 h 39 min** (00:27:58 to 02:07:17): 28/28 jobs, `VALIDATE: OK`,
31 PASS / 0 WARN / 0 FAIL. The cost is dominated by **fixed overhead**, not data
volume:

| Stage | Time |
|---|---|
| `metawrap_bins` (metaWRAP binning, bin refinement, per-bin reassembly) | 55m48s |
| `gtdbtk_classify` ×2 arms | 14m13s + 14m58s |
| `checkm_lineage_wf` ×2 arms | 5m15s + 5m29s |
| the pipeline's own `spades_control_assemble` | 85s |
| everything else (bowtie2, bam_to_fastq, refinem, semibin2, drep, summary, qa) | about 2 min total |

Therefore:

- **Changing code in `tools/` or `rules/` requires the full pass (`-F`), about 1h40m.** Do not
  use `--incremental`: the repository `run.sh` defaults `rerun-triggers` to `mtime`
  only, so changed script content does not trigger a recompute, and you get a
  false green ("it ran, but on the old result").
- **Changing only the config list dimensions (binners/methods/treats/groups) takes minutes with `--incremental`**,
  because the affected jobs are "missing outputs", not "need recompute".
- The first run, or a change of fixture genomes, needs the full pass.
