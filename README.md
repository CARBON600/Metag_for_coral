# Standardizing a Bioinformatic Analysis Pipeline for Coral Metagenomics

This repository is a Snakemake workflow for genome-resolved metagenomics of
coral holobiont samples. It removes coral host and Symbiodiniaceae reads,
assembles the remaining reads, recovers metagenome-assembled genomes (MAGs)
with four binners, and assesses their quality and taxonomy. 

Maintainer: `DAI Yifan` (`yidai204@gmail.com`).

## Experimental design

Each sample corresponds to one wet-lab protocol. The workflow enumerates every
combination of the following factors:

- host-read removal: `bt2`, `coverm`, `fastqs`
- library treatment: `control_assemble` (no PCR deduplication), `PCR_assemble` (PCR deduplication)
- reference database group: `control` (control database), `exp` (custom *Acropora* + Symbiodiniaceae database)
- binner: `unitem`, `comebin`, `metadecoder`, `semibin2_single`

That is `3 × 2 × 2 × 4 = 48` pipeline runs per sample. The four component
lists are defined in `config.yaml` and read by the `Snakefile`, so any subset
can be selected.

Dereplication (dRep) is part of this workflow, see
[dRep layer](#drep-layer-mag-dereplication) below. Functional annotation
(Prokka/KofamScan) and the statistical comparison of MAG sets are still done
downstream, in a separate notebook.

## Pipeline

```
Raw Reads (paired-end)
|
├── FastQ-Screen + fastp ──────────┐
├── Bowtie2 → unmapped → FASTQ ────┤
└── CoverM → inverse filter ───────┤
                                    ▼
                           Filtered Reads
                                    │
                    ┌───────────────┴───────────────┐
                    │                               │
              MEGAHIT pre-assembly            SPAdes meta assembly
           (ref for PCR dedup)              (control – no dedup)
                    │                               │
              PCR deduplication                     │
           (samtools markdup)                       │
                    │                               │
              SPAdes meta assembly                  │
              (PCR-assemble)                        │
                    │                               │
                    └───────────┬───────────────────┘
                                ▼
                     contigs ≥ 2 000 bp
                                │
                        Bowtie2 re-mapping
                                │
          ┌──────────────┬──────┴──────┬──────────────────┐
          ▼              ▼             ▼                  ▼
       UniT          COMEBin      MetaDecoder         SemiBin2
          │              │             │                  │
          └──────────────┴──────┬──────┴──────────────────┘
                                ▼
                          RefineM
                      (filter bins)
                                │
                    ┌───────────┴───────────┐
                    ▼                       ▼
              CheckM lineage_wf      GTDB-Tk classify_wf
                    │                       │
                    └───────────┬───────────┘
                                 ▼
                    MIMAG gate (comp ≥ 50, con ≤ 10)
                                 │
                    ┌────────────┴────────────┐
                    ▼                         ▼
         dRep per combination         dRep across combinations
          (within each of the 48)    (all passing MAGs at once)
                    │                         │
                    └────────────┬────────────┘
                                 ▼
                        Dereplicated MAG sets
```

> **Kaiju is disabled.** Profiling on the raw reads is not
> part of this workflow: the Kaiju rule has been **removed from**
> `rules/preprocess.smk` (only an explanatory comment remains), its target is
> absent from `rule all`, and its `kaiju:`/`threads.kaiju` config keys have been
> removed to avoid orphan entries. Rationale: the community-composition results
> for this study come from the genome-resolved MAGs (CheckM + GTDB-Tk), so Kaiju
> was dropped to keep the DAG focused. To reinstate it: (1) restore
> `kaiju: {nodes: "/mnt/storage3/thliao/data/kaiju_db/nodes.dmp", names:
> "/mnt/storage3/thliao/data/kaiju_db/names.dmp", db:
> "/mnt/storage3/thliao/data/kaiju_db/kaiju_db_nr_euk.fmi"}` and
> `threads.kaiju: 64` in `config.yaml`; (2) **rewrite the rule body from scratch
> in `rules/preprocess.smk`** (or recover it from git history) and add its
> outputs to `rule all`; (3) add the three database paths to `check_inputs` so
> they are gated before a long run. If the paper or final report cites Kaiju (or
> `kaiju_summary.xlsx`), it **must** be reinstated.

## Requirements

| Item | Requirement |
|------|-------------|
| OS | Linux (conda environments pinned to `linux-64`) |
| Snakemake | ≥ 7.0 |
| Conda | Miniconda, Anaconda, or Mamba |
| Storage | ≥ 500 GB recommended |

### Conda environments

Dependencies are split into one environment per tool family (the tools require
incompatible Python versions). The rules reference these environments by name;
`run.sh` creates any that are missing from the matching specification under
`envs/`, and reuses the ones that already exist. To create them manually:

```bash
conda env create -f envs/metag.yml
conda env create -f envs/megahit.yml
conda env create -f envs/binning.yml
conda env create -f envs/comebin_env.yml
conda env create -f envs/metadecoder.yml
conda env create -f envs/SemiBin.yml
conda env create -f envs/gtdbtk-2.3.2.yml
conda env create -f envs/drep.yml
```

Because Snakemake reuses these named environments, editing a specification does
not rebuild the environment automatically. Clone before refreshing, so a failed
solve cannot destroy a working environment:

```bash
conda create -n metag_bak --clone metag          # rollback point
conda env update -n metag -f envs/metag.yml      # apply the change
# or: conda env remove -n metag && conda env create -f envs/metag.yml
```

> **The deployed `metag` environment is frozen.** It was built before conda-forge
> published the repodata patches that pin `liblzma`/`xz` together, so it holds a
> combination (`liblzma 5.6.4` + `xz 5.2.6`) that no longer solves from a plain
> `conda env create`. `envs/metag.yml` has had the four conflicting split pins
> (`liblzma`, `liblzma-devel`, `xz-tools`, `xz-gpl-tools`) removed so a *fresh*
> build is solvable on a new machine — which also means the YAML no longer
> describes the deployed environment. Do **not** delete, `conda update`, or
> `conda env update --prune` the live `metag` env: it is the only known-good copy.
> Use this YAML for new machines only, and clone before changing any live env.

| Environment | Steps |
|-------------|-------|
| `metag` | host removal (`bt2`, `coverm`, `fastqs`), fastp, PCR dedup, SPAdes, contig filtering/remapping |
| `megahit` | MEGAHIT pre-assembly |
| `binning` | UniT, RefineM, CheckM, bin preparation |
| `comebin_env` | COMEBin |
| `metadecoder` | MetaDecoder |
| `SemiBin` | SemiBin2 |
| `gtdbtk-2.3.2` | GTDB-Tk |
| `drep` | dRep 3.6.2, dereplication (carries `mash` and `mummer4`/`nucmer`) |

### External tools

FastQ-Screen and SPAdes are invoked from absolute paths configured in
`config.yaml` (`software:`), with the original notebook paths as defaults:

| Tool | Config key |
|------|------------|
| FastQ-Screen 0.14.1 | `software.fastq_screen` |
| SPAdes 3.15.5 | `software.spades` |

> **Fragile external paths.** All three of these live outside lab storage:
> FastQ-Screen under `/home-user/thliao/...`, SPAdes under `/home-user/...`, and
> `checkm_data` **inside the `drep` conda environment** (a `drep` rebuild erases
> it). Do not clean those directories. Preferred fix: copy `checkm_data` to a lab
> path and re-run `checkm data setRoot <new_path>`, relocate FastQ-Screen and
> SPAdes if possible, then update the matching `config.yaml` keys. They are
> already validated by `check_inputs`/`check_qc_inputs`, so a missing path is
> reported before any long job starts. Never "repair" an environment with
> `conda env update --prune` (see the frozen-`metag` note above).

### Reference databases

| Resource | Purpose | Used by | Config key |
|----------|---------|---------|------------|
| Bowtie2 index | Host-read removal | `bt2` | `bowtie2_index` |
| CoverM reference FASTA | Host-read removal | `coverm` | `coverm_ref` |
| FastQ-Screen config | Multi-genome host screening | `fastqs` | `fastq_screen_conf` |
| GTDB-Tk reference data | Taxonomic classification | `gtdbtk` | `gtdbtk_data` |
| CheckM reference data | Completeness / contamination | `checkm` | `checkm_data` |

## Configuration

Edit `config.yaml`:

```yaml
data_dir: "/path/to/raw/data"      # contains *_1.fq.gz / *_2.fq.gz
tmpdir:   "/path/to/tmp"

filter_methods:
  - bt2
  - coverm
  - fastqs

assembly_treats:
  - control_assemble
  - PCR_assemble

binners:
  - unitem
  - comebin
  - metadecoder
  - semibin2_single

groups:
  - control
  - exp

bowtie2_index:
  control: "/path/to/control/control"
  exp:     "/path/to/exp/exp"

coverm_ref:
  control: "/path/to/control.fna"
  exp:     "/path/to/ref_Acropora_allSymb_comb.fa"

fastq_screen_conf:
  control: "/path/to/control.conf"
  exp:     "/path/to/exp.conf"

software:
  fastq_screen: "/path/to/FastQ-Screen-0.14.1/fastq_screen"
  spades:       "/path/to/spades/bin/spades.py"

gtdbtk_data: "/path/to/gtdb/releaseXXX"          # sets GTDBTK_DATA_PATH
checkm_data: "/path/to/checkm/data"              # sets CHECKM_DATA_PATH

threads:
  fastq_screen: 16
  fastp: 16
  bowtie2_map: 32
  # ... see config.yaml for the full list
```

The three component lists are read directly by the `Snakefile`; commenting out
a value prunes that branch of the DAG.

## Running

`run.sh` checks for the eight named environments and creates any that are
missing from `envs/*.yml`, then runs Snakemake with the conda frontend set to
`conda`. Any arguments are passed through to Snakemake. Unless you pass your own,
`run.sh` adds `--rerun-triggers=mtime` and `--keep-going`, so one failed job no
longer cancels the rest of the DAG:

```bash
# Validate the DAG
bash run.sh -n

# Record the planned job count (the dry-run "Job stats" block)
mkdir -p logs
bash run.sh -n --cores 64 2>&1 | tee "logs/run_$(date +%Y%m%d_%H%M).log"

# Single node
bash run.sh --cores 128 --keep-going

# Slurm (Snakemake 7.x)
bash run.sh --jobs 50 --keep-going --restart-times 2 \
    --cluster "sbatch --cpus-per-task={threads} --mem=64G --time=48:00:00" \
    --latency-wait 120

# Resume after a failure
bash run.sh --cores 128 --keep-going --rerun-incomplete

# Rule graph
bash run.sh -n --rulegraph | dot -Tpng > rulegraph.png
```

To call Snakemake directly instead, keep the frontend explicit:
`snakemake --use-conda --conda-frontend conda ...`.

`run.sh` must be run from the repository and needs `conda` on `PATH` (activate
your base or snakemake environment first).

### Threads and memory budget

Several rules declare 48–96 threads (`comebin`/`checkm` 96, `semibin2_single`/
`unitem`/`gtdbtk` 64, `remap_final_contigs` 48). With a smaller `--cores`,
Snakemake scales each rule down and can then run **only one job at a time**, so a
single deep chain monopolises the machine while the other arm of the matrix waits
— the usual cause of "one group finished, the other is empty". Keep the sum of
concurrent rule threads within the machine (raise `--cores`, lower the
`config.yaml` `threads.*`, or submit via the Slurm example above).

Threads and memory can be overridden without editing `config.yaml`:

```bash
bash run.sh --cores 32 --keep-going \
    --resources mem_mb=200000 \
    --set-threads semibin2_single=32 comebin=32 checkm_lineage_wf=32 \
      gtdbtk_classify=32 remap_final_contigs=32
```

Two traps with `--set-threads` (same root cause as `--rerun-triggers`):

- it is also `nargs='+'`, so it **eats positional targets**: put targets
  **before** the option (`bash run.sh output/x.done --set-threads a=8`) or
  separate with `--`, otherwise Snakemake reports `Unparseable value`;
- a **misspelled rule name is silently ignored**, so confirm the change took
  effect by checking that `Rules claiming more threads will be scaled down` no
  longer appears in the run log.

### Rerun and skip behaviour

`run.sh` adds `--rerun-triggers=mtime` unless you pass your own. The value must be
a **single token** (`--rerun-triggers=mtime`): Snakemake declares this option with
`nargs='+'`, so the space-separated form `--rerun-triggers mtime <target>` makes
argparse greedily swallow the following arguments and abort with
`invalid choice`. This disables the
`code`/`params`/`software-env` triggers so that editing rule code or conda
environments does not silently recompute everything; file modification times and
input checksums are still compared as usual. To force a full re-evaluation (for
example after editing a rule), use `--forceall`/`-F`:

```bash
bash run.sh -F --cores 128 --keep-going
```

After applying a change to rule code or to `envs/*.yml`, either run once with
`-F`/`--forceall` or delete the affected `.done` markers; otherwise the mtime-only
default skips the affected jobs and the change silently does not take effect.
(`--rerun-triggers` in Snakemake 7.32.4 only accepts
`mtime,params,input,software-env,code` — there is no `all`.)

Calling `snakemake` directly instead uses the default triggers
(`mtime,params,input,software-env,code`), so a changed rule or a rebuilt conda
environment re-executes the affected jobs. To keep existing results, pass
`--rerun-triggers=mtime`, or drop the provenance of a single output with
`snakemake --cleanup-metadata <output_file>`.

### Resuming after an interruption

An interrupted run (Ctrl-C, node crash, wall-clock kill) leaves partial outputs
marked incomplete, and Snakemake refuses to continue silently. Recover with:

```bash
snakemake --unlock                                   # clear locks left by a hard kill
bash run.sh --rerun-incomplete --cores 64 --keep-going
```

Finished steps are skipped; only the interrupted jobs are recomputed, after
their stale products are removed.

### Regenerating MetaDecoder output after this fix

MetaDecoder releases up to 1.2.1 share their `{kmers,dpgmm}` cache across
combinations, so any MetaDecoder bins produced before this workflow was patched
may be wrong. Remove **all** of them, not just the combinations that failed: the
shared cache could equally have been reused silently and produced a wrong bin
without any error. Then let Snakemake recompute; downstream QC follows
automatically because the input checksums change:

```bash
snakemake --unlock                       # only if a previous run was killed hard
rm -rf output/binning/metadecoder
rm -f contigs_r2000bp.fasta.2000.metadecoder.*   # stray cache at the repo root (outside output/)
bash run.sh --cores 64 --keep-going
```

Selecting a subset is done by editing the lists in `config.yaml`:

```yaml
filter_methods: [bt2, coverm]
assembly_treats: [control_assemble]
binners: [unitem, semibin2_single]
```

## Input data

Place paired-end reads under `data_dir` as `{sample}_1.fq.gz` and
`{sample}_2.fq.gz`. Samples are discovered by globbing `{sample}_1.fq.gz` and
are processed against every `group`. Sample identifiers should not contain
`_1`/`_2` motifs that the pairing convention could misread.

## Output

```
output/
├── fastqs/                 FastQ-Screen + fastp
├── bt2/                    Bowtie2 mapped BAM
├── bt2_unmapped/           Bowtie2 unmapped (-f 4)
├── coverm/                 CoverM mapped BAM
├── coverm_filtered/        CoverM inverse-filtered BAM
├── fq4dep/                 Filtered FASTQ for assembly
├── megahit_pre/            MEGAHIT pre-assembly
├── PCR_free/               Bowtie2 index + BAM vs MEGAHIT contigs
├── PCR_done/               PCR-deduplicated reads
├── assemble/               SPAdes contigs, ≥2000 bp contigs, sorted BAM
                             (SPAdes raw intermediates: .../{sample}/spades_work/)
├── binning/                Binner outputs
├── bins_prepared/          Decompressed / cleaned bins for QC
├── refinem/                RefineM-filtered bins
├── checkm/                 CheckM results
├── gtdb/                   GTDB-Tk results
├── summary/                MIMAG table, per-combination counts, final_mags.tsv
├── drep_per_workflow/      one dRep run per combination, plus the funnel table
└── drep_cross_workflow/    one dRep run over all passing MAGs
```

> **SPAdes intermediates moved.** Each SPAdes job now assembles inside
> `output/assemble/{treat}/{method}/{group}/{sample}/spades_work/` and only
> `contigs.fasta` is moved to its declared path. This stops a SPAdes rerun from
> wiping the downstream `contigs_r2000bp.fasta`, `bt2_index.done` and
> `*_sorted.bam` that live in the same directory. Intermediates left there by a
> pre-patch run (`scaffolds.fasta`, `spades.log`, `K*/`, `misc/`,
> `input_dataset.yaml`) must be deleted **once by hand**, otherwise old and new
> files coexist.

Per-rule logs are written to `logs/` at the repository root, not under `output/`.

### Harvesting results

`scripts/harvest_mags.py` turns the output tree into two flat tables without
running any tool. It is a **separate workflow, never part of `rule all`**:

```bash
snakemake -s rules/report.smk harvest --cores 1
```

- `harvest/combo_status.tsv` — one row per (binner × treat × method × group ×
  sample) combination, with per-stage done/product counts and a `status` column.
  Combinations that were never executed appear as `not_run`, so "half the matrix
  is missing" is visible at a glance. The statuses are defined **exactly** as:

  | status | meaning |
  |--------|---------|
  | `not_run` | no assembly (no `contigs.fasta`, `contigs_r2000bp.fasta` or `*_sorted.bam`) |
  | `partial` | assembly exists but the binner's `.done` marker is missing |
  | `empty` | binner finished but produced 0 bins |
  | `qc_empty` | bins are prepared but QC produced **no rows** (empty `bin_stats_ext.tsv` and empty GTDB summaries) |
  | `ok` | bins prepared **and** CheckM/GTDB-Tk produced numeric rows |

  `qc_empty` does not require a `.done` marker: a finished-but-empty QC run is
  never `ok`, so "QC silently returned nothing" cannot be mistaken for success.
- `harvest/mag_quality.tsv` — one row per MAG, joining CheckM (`completeness`,
  `contamination`, `strain_heterogeneity`, `genome_size`, `gc`, `contigs`,
  `marker_lineage`) with GTDB-Tk (`gtdb_domain`, `gtdb_taxonomy`).

`harvest/` is git-ignored. The tables are only as complete as the run: before any
CheckM/GTDB-Tk job has finished, `mag_quality.tsv` is empty.

Wildcard values:

| Wildcard | Values |
|----------|--------|
| `{group}` | `control`, `exp` |
| `{method}` | `bt2`, `coverm`, `fastqs` |
| `{treat}` | `control_assemble`, `PCR_assemble` |
| `{binner}` | `unitem`, `comebin`, `metadecoder`, `semibin2_single` |

## Workflow details

### Host-read removal

| Method | Tool(s) | Strategy |
|--------|---------|----------|
| `bt2` | Bowtie2 | Map to host reference, keep unmapped (`samtools view -f 4`), convert to FASTQ |
| `coverm` | CoverM | `coverm make` to a directory, then inverse filter (≥75% aligned, ≥95% identity) |
| `fastqs` | FastQ-Screen + fastp | Screen against the multi-genome config (`--nohits`), then fastp QC |

### Assembly

| Treatment | Input | PCR dedup |
|-----------|-------|-----------|
| `control_assemble` | Filtered reads | No |
| `PCR_assemble` | Reads after MEGAHIT pre-assembly + `samtools markdup` | Yes |

MEGAHIT pre-assembly (`--min-contig-len 250`) produces contigs used only as a
reference for PCR duplicate marking; SPAdes (`--meta`,
`-k 21,33,55,77,99,127`) is the final assembler for both treatments. Contigs
shorter than 2000 bp are discarded before binning.

### Binning

| Tool | Environment | Approach |
|------|-------------|----------|
| UniT | `binning` | Ensemble of MetaBAT2 + MaxBin2 (consensus) |
| COMEBin | `comebin_env` | Contrastive multi-view representation learning (device auto-selected via `torch.cuda.is_available()`; `envs/comebin_env.yml` ships CPU PyTorch, so training runs on CPU unless a GPU build is installed) |
| MetaDecoder | `metadecoder` | Two-layer DPGMM + k-mer frequency model |
| SemiBin2 | `SemiBin` | Self-supervised contrastive learning (`single_easy_bin --self-supervised`) with a fixed `--random-seed 1`, so bins are reproducible across runs |

> **MetaDecoder cache warning.** Releases up to 1.2.1 write a
> `<basename(fasta)>.<min_len>.metadecoder.{kmers,dpgmm}` cache in the working
> directory and a `<bam>.index` cache **next to the input BAM**
> (`output/assemble/...`). This workflow runs the `cluster` step in a private
> working directory and removes the BAM index before `coverage`, so caches are
> never reused across combinations. Any MetaDecoder results produced **before**
> this fix must be considered contaminated and regenerated (see below).

### Bin preparation and quality assessment

`prepare_bins` normalises each binner's output into `output/bins_prepared/`:
gzipped bins from UniT (`*.fna.gz`) and SemiBin2 (`*.fa.gz`) are decompressed,
and only bin FASTA files are copied, so `.COVERAGE`/`.SEED`/`*.tsv` files do not
reach downstream tools. Bin extensions: UniT `.fna`, COMEBin `.fa`,
MetaDecoder `.fasta`, SemiBin2 `.fa`.

| Tool | Function |
|------|----------|
| RefineM | Scaffold statistics, outlier detection, bin filtering |
| CheckM | Lineage-specific completeness and contamination (`lineage_wf`) |
| GTDB-Tk | Taxonomic classification (`classify_wf --skip_ani_screen`) |

> **GTDB-Tk memory.** pplacer allocates ~60 GB for internal nodes (observed up to
> 61.4 GB in the logs), so the rule declares `mem_mb=96000` (~1.6× the peak) to
> keep two GTDB-Tk jobs off a node that only fits one. On a single node `mem_mb`
> is only a scheduling declaration and is not enforced; cap the node with
> `--resources mem_mb=<node RAM>` so the scheduler honours it.

### Output hygiene and empty results

Most rules delete their own products before rebuilding (the exception is
`fastq_screen`, which is re-run with `--force`), so a `.done` marker means the
outputs on disk are fresh and complete. A binner that returns no bins is not
an error: it logs a `WARNING`, and RefineM/CheckM/GTDB-Tk are skipped, writing
empty placeholder tables (`scaffold_stats.tsv`, `outliers.tsv`,
`storage/bin_stats_ext.tsv`, `<prefix>.bac120.summary.tsv`). Run with
`--keep-going` so a single empty branch does not stop the rest.

`prepare_bins` still fails hard when the raw binner directory contains files
that match no expected bin extension, so a naming/extension mismatch is never
recorded as "done with zero bins".

## Deployment on a new machine

1. Install Miniconda/Mamba with `conda ≥ 24.7.1` on `PATH`, and Snakemake ≥ 7 (`mamba create -n snakemake -c bioconda -c conda-forge snakemake=7.32.4`).
2. Clone the repository and place the paired-end reads under `data_dir`.
3. Provide the external tools and reference databases, then update every path in `config.yaml` (`data_dir`, `tmpdir`, `bowtie2_index`, `coverm_ref`, `fastq_screen_conf`, `software`, `gtdbtk_data`, `checkm_data`) and adjust `threads`.
4. Create the eight environments (`mamba env create -f envs/<name>.yml`) or let `run.sh` create the missing ones on first use.
5. Validate the DAG (`bash run.sh -n`), then run with `bash run.sh ...`; resume with `--rerun-incomplete`.

The `check_inputs` rule verifies that all configured input paths exist before
any long job starts.

### Troubleshooting

| Symptom | Probable cause | Remedy |
|---|---|---|
| `check_inputs` fails with `no paired-end samples matched` | `data_dir` contains no `{sample}_1.fq.gz`, so `SAMPLES` is empty | correct `data_dir` |
| `FileNotFoundError` for a reference or tool | a path in `config.yaml` was not updated | update `config.yaml` |
| `check_inputs` fails | one of the configured paths does not exist | read the log under `logs/check_inputs.log` |
| CheckM reports the data folder is unset | `CHECKM_DATA_PATH` not configured | set `checkm_data`, or run `checkm data setRoot <dir>` once |
| `The 'mamba' command is not available` | Snakemake's default conda frontend is `mamba` | use `run.sh` (sets `--conda-frontend conda`) or install mamba |
| `Conda must be version 24.7.1 or later` | a Snakemake 8.x install requires a newer conda | pin Snakemake 7.32.4 (see above) or upgrade conda |
| Conda environment creation fails | `defaults` ToS, network, or a missing build | accept the ToS, use a mirror, or relax pins |
| SPAdes exits with a Python error | interpreter mismatch of the external script | invoke via `conda run -n <env>` or install SPAdes into the environment |
| Disk exhaustion | `tmpdir` or `output/` under-provisioned | relocate to larger storage |
| Edited a rule but nothing reran | `--rerun-triggers=mtime` is the default | run with `-F`/`--forceall`, or delete the relevant `.done` |
| `invalid choice: '<target>'` for `--rerun-triggers` | the space-separated form was used; the option is `nargs='+'` and eats the target | write it as a single token: `--rerun-triggers=mtime` |
| `ERROR: fastq_screen conf parsed to zero entries` | the fastq_screen conf is empty or contains only comments | fix the conf; `SKIP_CONF_CHECK=1` does not bypass this, nor a missing `DATABASE` index |
| A binner log shows `WARNING: ... no bins` | that branch recovered no MAGs | expected; downstream QC is skipped with empty placeholders |
| `Directory cannot be locked` | a previous run was killed hard | run `snakemake --unlock`, then rerun with `--rerun-incomplete` |
| Cluster jobs fail on node hiccups | `--restart-times` defaults to 0 | add `--restart-times 2` (see the Slurm example) |

## Tools

| Tool | Version | Environment | Reference |
|------|---------|-------------|-----------|
| MEGAHIT | 1.2.9 | `megahit` | Li et al., 2015 |
| SPAdes | 3.15.5 | external | Bankevich et al., 2012 |
| Bowtie2 | 2.5.4 | `metag` | Langmead & Salzberg, 2012 |
| CoverM | 0.7.0 | `metag` | Aroney et al., 2025 |
| FastQ-Screen | 0.14.1 | external | Wingett & Andrews, 2018 |
| fastp | 1.3.7 | `metag` | Chen et al., 2018 |
| samtools | 1.18 | `metag` | Danecek et al., 2021 |
| UniT | 1.0.2 | `binning` | Parks et al. |
| COMEBin | 1.0.4 | `comebin_env` | Wang et al., 2024 |
| MetaDecoder | 1.2.1 | `metadecoder` | Liu et al., 2022 |
| SemiBin2 | 2.2.0 | `SemiBin` | Pan et al., 2023 |
| RefineM | 0.1.2 | `binning` | Parks et al. |
| CheckM | 1.2.3 | `binning` | Parks et al., 2015 |
| GTDB-Tk | 2.3.2 | `gtdbtk-2.3.2` | Chaumeil et al., 2020 |
| dRep | 3.6.2 | `drep` | Olm et al., 2017 |
| Kaiju | – | disabled | Menzel et al., 2016 |

## dRep layer (MAG dereplication)

After CheckM/GTDB-Tk, a gene-resolution dereplication layer turns each
combination's bins into non-redundant MAG sets. It is a main-DAG extension and
runs entirely in the `drep` environment. Five rules make it up: `mimag_gate`
(`rules/summary.smk`) and `drep_per_workflow`, `drep_cross_workflow`,
`drep_collect`, `drep_taxonomy` (`rules/drep.smk`).

```
output/checkm/**/storage/bin_stats_ext.tsv   (the single quality source)
        │
   MIMAG gate: comp >= 50 and con <= 10      Bowers 2017 Nat Biotechnol, PMID 28787424
        │  -> output/summary/mimag_bins.tsv        (passing MAGs; the only join table)
        │     output/summary/mimag_summary.tsv     (per-combination counts)
        │     output/summary/mimag_excluded.tsv    (rejected MAGs + exclude_reason)
        ▼
   dRep -comp 50 -con 10 -l 0 --genomeInfo  (the same cutoffs as the gate)
        │  -sa 0.95 (~95% ANI; PMID 19855009, 30504855)
        ├─ output/drep_per_workflow/<binner>/<treat>/<method>/<group>/<sample>/
        │     n_rep = files in dereplicated_genomes/  (non-redundant MAG count)
        │     n_merged = n_mimag - n_rep
        │     output/drep_per_workflow/per_workflow.tsv  (funnel: n_total/n_mimag/n_rep/...)
        └─ output/drep_cross_workflow/data_tables/{Bdb,Cdb,Wdb}.csv  (all MIMAG bins, one job)
              Cdb = genome -> secondary_cluster (membership)
              Wdb = genome + cluster + score   (representative / winner)
        ▼
   output/summary/final_mags.tsv  (membership + representative taxonomy propagation)
```

| Item | Value |
|------|-------|
| Quality cutoffs | `comp >= 50` and `con <= 10`, both inclusive; dRep applies the same two values through `-comp`/`-con`, and `rules/drep.smk` asserts that the `mimag` and `drep` config blocks agree |
| Length filter | off (`-l 0`); dRep only applies it when `length > 1` |
| Clustering | `-pa 0.90 -sa 0.95 -nc 0.1`, `--S_algorithm fastANI` (passed explicitly; the 3.6.2 argparse default is `fastANI`) |
| `--genomeInfo` | 4-column CSV `genome,completeness,contamination,strain_heterogeneity`, taken from the MIMAG gate's parse and covering every staged genome (a gap would make dRep run its own CheckM) |
| Winner score | `comp - 5*con + con*(strh/100) + 0.5*log10(N50) + 0*log10(size) + 1*(centrality - 0.95)` |

Per-pool behaviour: 0 MIMAG → `skipped_empty` (no dRep call, empty
`data_tables/*.csv` placeholders); 1 MIMAG → `passthrough_single` (copied into
`dereplicated_genomes/`, no dRep call); ≥2 → `dereplicated`. Each dRep run uses a
private, freshly wiped work directory. Two checks guard the layer: within a pool,
the manifest, the staged files, the `genomeInfo.csv` rows and the `g_source.list`
lines must all have the same size; and after a dereplicated run the log must
report `100.00% of genomes passed checkM filtering`, which is what the shared
cutoffs imply, so a lower number means the gate and dRep disagreed.
`scripts/drep_stage.py` implements staging, verification, summaries, the funnel
table and the L9 join.

> **dRep work directories.**
> (1) dRep's `WorkDirectory.overwrite` is hardcoded to `True`, so dRep itself
> does not protect an existing work directory; the per-pool directory is private
> and wiped before each run, and dRep refuses a `-g` list when `Bdb` already
> exists (`assert wd.hasDb("Bdb") == False`). Do not point `dRep` at a shared or
> already-populated work directory.
> (2) `dRep compare` overwrites `Cdb/Mdb/Ndb/Widb` but not `Bdb/Wdb`, so a
> directory it touches holds tables from two different runs. It must not share
> the cross-workflow output directory; if `compare` is ever needed, drop `-g` so
> it reuses the existing `Bdb`, or use a fresh directory.

> **Why `fastANI`.** The originally published run used dRep's 3.6.2 default,
> which is `fastANI` (the `ANImf = (DEFAULT)` line in `dRep -h` is stale). The
> workflow passes it explicitly so the value cannot drift with the installed
> version. Switch `drep.s_algorithm` (and record it in Methods) if ANImf is
> required; it aligns whole genomes with `nucmer` and is far slower.

> **Regenerating the dRep layer.** Editing any rule/script here does not
> re-trigger jobs by default (`--rerun-triggers=mtime`); run once with `-F`, or
> delete the affected `derep.done` / `output/drep_*` trees. dRep's work directory
> is wiped and rebuilt per pool, so nothing stale survives a rerun.

## License

[GPL v3](LICENSE)
