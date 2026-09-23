# Standardizing a Bioinformatic Analysis Pipeline for Coral Metagenomics

This repository is a Snakemake workflow for genome-resolved metagenomics of
coral holobiont samples. It removes coral host and Symbiodiniaceae reads,
assembles the remaining reads, recovers metagenome-assembled genomes (MAGs)
with four binners plus the metaWRAP ensemble arm, and assesses their quality and
taxonomy. 

Maintainer: `DAI Yifan` (`yidai204@gmail.com`).

## Table of Contents

- [Experimental design](#experimental-design)
- [Pipeline](#pipeline)
- [Requirements](#requirements)
  - [Conda environments](#conda-environments)
  - [External tools](#external-tools)
  - [Reference databases](#reference-databases)
- [Configuration](#configuration)
- [Running](#running)
  - [Threads and memory budget](#threads-and-memory-budget)
  - [Rerun and skip behaviour](#rerun-and-skip-behaviour)
  - [Resuming after an interruption](#resuming-after-an-interruption)
  - [Regenerating MetaDecoder output after this fix](#regenerating-metadecoder-output-after-this-fix)
- [Input data](#input-data)
- [Output](#output)
  - [Harvesting results](#harvesting-results)
- [Workflow details](#workflow-details)
  - [Host-read removal](#host-read-removal)
  - [Assembly](#assembly)
  - [Binning](#binning)
  - [metaWRAP arm](#metawrap-arm)
  - [Bin preparation and quality assessment](#bin-preparation-and-quality-assessment)
  - [Output hygiene and empty results](#output-hygiene-and-empty-results)
- [Reproducing from scratch](#reproducing-from-scratch)
  - [Regression baseline](#regression-baseline)
  - [Troubleshooting](#troubleshooting)
- [Tools](#tools)
- [dRep layer (MAG dereplication)](#drep-layer-mag-dereplication)
- [Citation](#citation)
- [License](#license)

## Experimental design

Each sample corresponds to one wet-lab protocol. The workflow enumerates every
combination of the following factors:

- host-read removal: `bt2`, `coverm`, `fastqs`
- library treatment: `control_assemble` (no PCR deduplication), `PCR_assemble` (PCR deduplication)
- reference database group: `control` (control database), `exp` (custom *Acropora* + Symbiodiniaceae database)
- binner: `unitem`, `comebin`, `metadecoder`, `semibin2_single`, `metawrap`

The first four binners run on every library treatment, giving
`3 × 2 × 2 × 4 = 48` pipeline runs per sample. The fifth entry, `metawrap`, is a
full ensemble sub-pipeline (see [metaWRAP arm](#metawrap-arm)), not a single
binner, and by default also runs on both treatments (its treatment list is
`metawrap_treats` in `config.yaml`, a subset of `assembly_treats`): that adds
`3 × 2 × 2 × 1 = 12` more, i.e. **60 combinations per sample** with the shipped
config (see `BINNER_TREAT_PAIRS`/`COMBOS` in the `Snakefile`). All component
lists are defined in `config.yaml` and read by the `Snakefile`, so any subset can
be selected.

Dereplication (dRep) is part of this workflow, see
[dRep layer](#drep-layer-mag-dereplication) below. Functional annotation
(Prokka/KofamScan) and the statistical comparison of MAG sets are still done
downstream, in a separate notebook.

## Pipeline

```
Raw reads (paired-end: {sample}_1.fq.gz / {sample}_2.fq.gz)

  Host-read removal: three alternative methods (`filter_methods`); each one
  yields its own filtered read set for every downstream combination:
    fastqs   FastQ-Screen → fastp                 → output/fastqs/<group>/<sample>/{sample}_[12]_cleaned.fastq.gz
    bt2      Bowtie2 map → unmapped → FASTQ        → output/fq4dep/bt2/<group>/<sample>_{1,2}.fq.gz
    coverm   CoverM map → inverse filter → FASTQ   → output/fq4dep/coverm/<group>/<sample>_{1,2}.fq.gz
                                  │
                ┌─────────────────┴─────────────────┐
                │                                   │
        control_assemble                      PCR_assemble
        (no deduplication)                    (PCR deduplication)
                │                                   │
                │                         MEGAHIT pre-assembly
                │                         (reference for dedup)
                │                                   │
                │                         Bowtie2 re-map reads to
                │                         MEGAHIT contigs
                │                                   │
                │                         PCR deduplication: samtools
                │                         collate → fixmate → markdup,
                │                         duplicates dropped, BAM → FASTQ
                │                                   │
          SPAdes meta assembly             SPAdes meta assembly
          (--meta --only-assembler)             (PCR-assemble)
                │                                   │
                └─────────────────┬─────────────────┘
                                  ▼
                        contigs > 2 000 bp
                                  │
                  Bowtie2 re-mapping to the final contigs
                  → {sample}_sorted.bam (+ .bai)
                                  │
   ┌──────────────┬───────────────┼────────────────┬───────────────────────┐
   ▼              ▼               ▼                ▼                       ▼
  UniT         COMEBin      MetaDecoder        SemiBin2          metaWRAP arm (5th)
 (+consensus)              (coverage →     (single sample)   bins the same contigs
   │              │         seed → cluster)        │          with its OWN coverage
   │              │               │                │          (bwa + jgi_summarize):
   │              │               │                │          MetaBAT2 + MaxBin2 +
   │              │               │                │          CONCOCT → bin_refinement
   │              │               │                │          → reassemble_bins
   └──────────────┴───────┬───────┴────────────────┘                    │
                          ▼                                             │
                     RefineM                                     no RefineM
                  (filter bins)                    (its bins are per-bin reassemblies
                          │                         whose contigs are absent from
                          │                         contigs_r2000bp.fasta and
                          │                         {sample}_sorted.bam; a passthrough
                          │                         refinem.done is emitted instead)
                          └────────────────────┬────────────────────────┘
                                               ▼
                            ┌──────────────────┴──────────────────┐
                            ▼                                     ▼
                  CheckM lineage_wf (1.2.3)            GTDB-Tk classify_wf
                  + checkm qa → results.tsv
                            │                                     │
                            └──────────────────┬──────────────────┘
                                               ▼
                              MIMAG gate (comp ≥ 50, con ≤ 10)
                                               │
                            ┌──────────────────┴──────────────────┐
                            ▼                                     ▼
               dRep per combination                  dRep across combinations
        (one run per combination: 48 from the              (all passing MAGs
         four binners + 12 from the metaWRAP                    at once)
         arm = 60 per sample)
                            │                                     │
                            └──────────────────┬──────────────────┘
                                               ▼
                                   Dereplicated MAG sets
```

> Kaiju is disabled. Profiling on the raw reads is not
> part of this workflow: the Kaiju rule has been removed from
> `rules/preprocess.smk` (only an explanatory comment remains), its target is
> absent from `rule all`, and its `kaiju:`/`threads.kaiju` config keys have been
> removed to avoid orphan entries. Rationale: the community-composition results
> for this study come from the genome-resolved MAGs (CheckM + GTDB-Tk), so Kaiju
> was dropped to keep the DAG focused. To reinstate it: (1) restore
> `kaiju: {nodes: "/mnt/storage3/thliao/data/kaiju_db/nodes.dmp", names:
> "/mnt/storage3/thliao/data/kaiju_db/names.dmp", db:
> "/mnt/storage3/thliao/data/kaiju_db/kaiju_db_nr_euk.fmi"}` and
> `threads.kaiju: 64` in `config.yaml`; (2) rewrite the rule body from scratch
> in `rules/preprocess.smk` (or recover it from git history) and add its
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
conda env create -f envs/metawrap-env.yml
```

The Snakemake driver itself runs in its own environment (not created by
`run.sh`). `envs/snakemake.yml` is the authoritative freeze of that environment
(`conda env export` from the deployment machine: 191 pins, Snakemake 7.32.4,
pulp 2.7.0, Python 3.11.16); its `prefix:` line was removed so the environment
can be created anywhere:

```bash
mamba env create -f envs/snakemake.yml
```

`run.sh` reads the pinned version out of that file and warns if the active
Snakemake is a different one.

Because Snakemake reuses these named environments, editing a specification does
not rebuild the environment automatically. Clone before refreshing, so a failed
solve cannot destroy a working environment:

```bash
conda create -n metag_bak --clone metag          # rollback point
conda env update -n metag -f envs/metag.yml      # apply the change
# or: conda env remove -n metag && conda env create -f envs/metag.yml
```

> The deployed `metag` environment is frozen. It was built before conda-forge
> published the repodata patches that pin `liblzma`/`xz` together, so it holds a
> combination (`liblzma 5.6.4` + `xz 5.2.6`) that no longer solves from a plain
> `conda env create`. `envs/metag.yml` has had the four conflicting split pins
> (`liblzma`, `liblzma-devel`, `xz-tools`, `xz-gpl-tools`) removed so a *fresh*
> build is solvable on a new machine, which also means the YAML no longer
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
| `metawrap-env` | metaWRAP 1.3.2 ensemble arm (metaBAT2/MaxBin2/CONCOCT, bin_refinement, reassemble_bins) + its bundled CheckM 1.0.12 |

### External tools

FastQ-Screen and SPAdes are invoked from absolute paths configured in
`config.yaml` (`software:`), with the original notebook paths as defaults:

| Tool | Config key |
|------|------------|
| FastQ-Screen 0.14.1 | `software.fastq_screen` |
| SPAdes 3.15.5 | `software.spades` |

> Fragile external paths. All three of these live outside lab storage:
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

Copy the annotated template and fill the 12 machine-specific paths (plus
`metawrap.local_base`):

```bash
cp config.example.yaml config.yaml
```

`config.yaml` ships with the same placeholders as `config.example.yaml`, so
either file can be the starting point. Edit `config.yaml`:

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
  - metawrap          # ensemble arm; requires the two `metawrap:` keys below

# Optional: treatments the metaWRAP arm runs (subset of assembly_treats).
# Remove the key to run it on every treatment.
metawrap_treats:
  - control_assemble
  - PCR_assemble

metawrap:
  env: metawrap-env               # named conda env (pre-created; run.sh reuses it)
  guard: "tools/mw_guard.sh"      # REQUIRED: path relative to the repo root
  guard_version: "1.1"            # documentation; identity is enforced by guard_md5
  guard_md5: "50e6127e06deebd82c751e8ce575e0ab"
  local_base: "/path/to/node-local/scratch"  # REQUIRED: short, node-local scratch dir
  refine_mem_gb: 40
  reassemble_mem_gb: 40
  keep_short_work: 0
  archive_existing: 0

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

`run.sh` checks for the nine named environments and creates any that are
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
your base or snakemake environment first). If conda is not on `PATH`, set
`CONDA_BASE=/path/to/miniconda3`; the script has no hardcoded install path.

### Threads and memory budget

Several rules declare 48–96 threads (`comebin`/`checkm` 96, `semibin2_single`/
`unitem`/`gtdbtk` 64, `remap_final_contigs` 48). With a smaller `--cores`,
Snakemake scales each rule down and can then run **only one job at a time**, so a
single deep chain monopolises the machine while the other arm of the matrix waits,
which is the usual cause of "one group finished, the other is empty". Keep the sum of
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
`mtime,params,input,software-env,code`; there is no `all`.)

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
│                            (metaWRAP arm: metawrap/<...>/04_FINAL_BINS_FOR_GTDB/,
│                             with the guard's 01_/02_/03_ trees alongside)
├── bins_prepared/          Decompressed / cleaned bins for QC
├── refinem/                RefineM-filtered bins (metaWRAP arm: passthrough)
├── checkm/                 CheckM results (checkm.done + qa `results.tsv`)
├── gtdb/                   GTDB-Tk results
├── summary/                MIMAG table, per-combination counts, final_mags.tsv
├── drep_per_workflow/      one dRep run per combination, plus the funnel table
└── drep_cross_workflow/    one dRep run over all passing MAGs
```

> SPAdes intermediates moved. Each SPAdes job now assembles inside
> `output/assemble/{treat}/{method}/{group}/{sample}/spades_work/` and only
> `contigs.fasta` is moved to its declared path. This stops a SPAdes rerun from
> wiping the downstream `contigs_r2000bp.fasta`, `bt2_index.done` and
> `*_sorted.bam` that live in the same directory. Intermediates left there by a
> pre-patch run (`scaffolds.fasta`, `spades.log`, `K*/`, `misc/`,
> `input_dataset.yaml`) must be deleted **once by hand**, otherwise old and new
> files coexist.

Per-rule logs are written to `logs/` at the repository root, not under `output/`.

### Harvesting results

`tools/harvest_mags.py` turns the output tree into two flat tables without
running any tool. It is a **separate workflow, never part of `rule all`**:

```bash
snakemake -s rules/report.smk harvest --cores 1
```

- `harvest/combo_status.tsv`: one row per (binner × treat × method × group ×
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
- `harvest/mag_quality.tsv`: one row per MAG, joining CheckM (`completeness`,
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
| `{binner}` | `unitem`, `comebin`, `metadecoder`, `semibin2_single`, `metawrap` |

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

> MetaDecoder cache warning. Releases up to 1.2.1 write a
> `<basename(fasta)>.<min_len>.metadecoder.{kmers,dpgmm}` cache in the working
> directory and a `<bam>.index` cache **next to the input BAM**
> (`output/assemble/...`). This workflow runs the `cluster` step in a private
> working directory and removes the BAM index before `coverage`, so caches are
> never reused across combinations. Any MetaDecoder results produced **before**
> this fix must be considered contaminated and regenerated (see below).

### metaWRAP arm

The `metawrap` binner is **not a single binner**: it is metaWRAP's ensemble
sub-pipeline: `binning` (metaBAT2 + MaxBin2 + CONCOCT), `bin_refinement`
(cross-set consolidation/de-replication), and `reassemble_bins` (per-bin SPAdes
reassembly), driven by an external guard wrapper configured at
`metawrap.guard` (with a short node-local scratch dir at `metawrap.local_base`).
The guard stages the final bins under
`output/binning/metawrap/<treat>/<method>/<group>/<sample>/04_FINAL_BINS_FOR_GTDB/`
and archives its own CheckM 1.0.12 report under `.../03_BIN_REASSEMBLY/`.

- **Guard identity is pinned.** `config.yaml` ships `guard: "tools/mw_guard.sh"`
  (a path relative to the repository root, resolved against the Snakefile's
  directory) with `guard_version: "1.1"` and
  `guard_md5: "50e6127e06deebd82c751e8ce575e0ab"`. `rules/metawrap.smk` aborts at
  parse time when `metawrap` is in `binners` and the guard is empty, missing,
  unreadable, or its md5 differs from `guard_md5`. Do not edit the guard: its
  md5 is asserted, and any change (including an EOL translation) breaks the run.
  `.gitattributes` marks it `-text` so it is stored byte-for-byte.
- **First use requires a one-time environment patch.** The guard writes into the
  active `metawrap-env` at runtime (not captured by `envs/metawrap-env.yml`):
  `install_rmtree_patch()` appends `${CONDA_PREFIX}/lib/python2.7/site-packages/sitecustomize.py`,
  and `patch-l3` rewrites `.../metawrap-modules/reassemble_bins.sh` leaving a
  `.bak`. Before the first run, apply and verify:
  `bash tools/mw_guard.sh patch-l3 && bash tools/mw_guard.sh verify-l3`.
  (The effect of these patches on *results* has not been tested with a
  positive/negative control; code reading only.)

- **RefineM is skipped for this arm** (`metawrap_ingest` passes the bins through
  to `output/refinem/metawrap/.../`): metaWRAP's reassembled contigs are not in
  `contigs_r2000bp.fasta` / `{sample}_sorted.bam`, which RefineM's
  `scaffold_stats` requires. `metawrap_ingest` also cross-checks that the final
  bin basenames equal the guard's CheckM 1.0.12 key column when that archive is
  present.
- **Reported completeness/contamination come from the shared `checkm_lineage_wf`
  (CheckM 1.2.3)**, so every arm is scored on one ruler and the cross-workflow
  dRep winner is not biased by mixing CheckM versions. metaWRAP's internal
  preselection, however, still uses its bundled CheckM 1.0.12.
- **Three metaWRAP-only filters**: `bin_refinement` silently drops input bins
  outside 50 kb–20 Mb; the `-c/-x` cherry-pick (fed from `mimag`) keeps only
  bins passing the MIMAG cutoffs; `choose_best_bin.py` then picks, per genome,
  `completeness + 5 × (100 − contamination)` (N50 breaks ties).
- **Coverage is computed by metaWRAP itself** (`bwa index`/`bwa mem` +
  `jgi_summarize_bam_contig_depths`), not from the pipeline's bowtie2 BAMs.
- **0 surviving MAGs is a job failure**, not an `empty` status: metaWRAP
  hard-errors when nothing passes.

Set `metawrap_treats` (a subset of `assembly_treats`) to restrict this arm;
remove the key to run it on every treatment. `rules/metawrap.smk` aborts at parse
time when `metawrap.guard` or `metawrap.local_base` is empty, when the guard file
is missing/unreadable, or when its md5 does not match `metawrap.guard_md5`.

### Bin preparation and quality assessment

`prepare_bins` normalises each binner's output into `output/bins_prepared/`:
gzipped bins from UniT (`*.fna.gz`) and SemiBin2 (`*.fa.gz`) are decompressed,
and only bin FASTA files are copied, so `.COVERAGE`/`.SEED`/`*.tsv` files do not
reach downstream tools. Bin extensions: UniT `.fna`, COMEBin `.fa`,
MetaDecoder `.fasta`, SemiBin2 `.fa`, metaWRAP `.fa` (the metaWRAP arm reads the
guard's `04_FINAL_BINS_FOR_GTDB/`, then `metawrap_ingest` passes those bins
through to `output/refinem/metawrap/` without RefineM).

| Tool | Function |
|------|----------|
| RefineM | Scaffold statistics, outlier detection, bin filtering |
| CheckM | Lineage-specific completeness and contamination (`lineage_wf`) |
| GTDB-Tk | Taxonomic classification (`classify_wf --skip_ani_screen`) |

> **GTDB-Tk memory is a per-thread figure.** `classify_wf` passes
> `-j <--pplacer_cpus>`, and when that option is omitted pplacer falls back to
> `--cpus` (`gtdbtk/classify.py`: `self.pplacer_cpus = max(pplacer_cpus if
> pplacer_cpus else cpus, 1)`). pplacer caches the reference likelihoods into a
> `VIRT = RES = <parent footprint>` block and then forks one worker per
> `--pplacer_cpus`; Unix copy-on-write shares those pages, so **the host reports
> `PARENT_MEMORY * (N_CHILDREN + 1)` without using additional physical memory**.
> That *reported* total — not the physical footprint — is what a cgroup or an HPC
> allocation accounts for, which is why the same job can be killed on one node
> and pass on another. Official guidance for hitting this: `--scratch_dir`
> (mmap to disk, trades RAM for speed) and `--pplacer_cpus 1` — see
> <https://ecogenomics.github.io/GTDBTk/faq.html> (*"GTDB-Tk reaches the memory
> limit / pplacer crashes"*).
>
> Three different numbers apply, so quote each with its source. For bacteria with
> the reference data this pipeline pins (`gtdbtk-2.3.2`, r207/r214) GTDB-Tk's own
> threshold is **`PPLACER_MIN_RAM_BAC_SPLIT = 55` GB** (`gtdbtk/config/common.py`;
> ~61 GB observed here), the current published requirement for the bacterial
> reference tree is **~140 GB** (**950 GB** with `--full_tree`,
> <https://ecogenomics.github.io/GTDBTk/installing/index.html#hardware-requirements>),
> and the FAQ's 150 GB is the illustration used in its own example. `--full_tree`
> in 2.3.2 is guarded at `PPLACER_MIN_RAM_BAC_FULL = 320` GB, archaea at
> `PPLACER_MIN_RAM_ARC = 40` GB.
>
> A memory kill is silent by design: the rule fails with a bare
> `PplacerException` while
> `classify/intermediate_results/pplacer/tree_*/pplacer.class_level.<domain>.out`
> stops after `Caching likelihood information ... done.` at `Preparing the edges
> for baseball...` **with no error line** — the process was killed, it did not
> report a failure. Re-run the single job before touching parameters:
>
> ```bash
> conda activate snakemake
> bash run.sh --configfile tools/mini/fixture/config.mini.yaml --cores 16 \
>   output/gtdb/<binner>/<treat>/<method>/<group>/<sample>/gtdb.done
> ```
>
> The rule therefore pins `--pplacer_cpus 1` and `--scratch_dir`, and declares
> `retries: 1` — rule-level retries are Snakemake's mechanism for fallible rules,
> with the `--retries` CLI overriding it globally for unreliable clusters
> (<https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#defining-retries-for-fallible-rules>).
> `resources.gtdbtk_mem_mb` (200000) is the scheduler-facing declaration: it
> covers the reported `2 x 55-61 GB` of the pinned data and still has room for
> `2 x 140 GB` after a reference-data upgrade. On a single node `mem_mb` is only
> a declaration and is not enforced; cap the node with
> `--resources mem_mb=<node RAM>`, and prefer a compute node over a login node.

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

## Reproducing from scratch

1. Install Miniconda/Mamba with `conda ≥ 24.7.1` on `PATH`, then create the
   Snakemake driver environment from the pinned specification and activate it:
   ```bash
   mamba env create -f envs/snakemake.yml
   conda activate snakemake
   ```
   `envs/snakemake.yml` is the `conda env export` of the deployment machine's
   driver environment (Snakemake 7.32.4, pulp 2.7.0, Python 3.11.16); its
   `prefix:` line was stripped so the file works on another machine. The personal
   `ursky` channel was also dropped (verified unused here: no dependency in this
   file is hosted on it; the same check removed it from `envs/drep.yml`, and it is
   required only by `envs/metawrap-env.yml`). The
   minimal equivalent, if you prefer to build it by hand, is
   `mamba create -n snakemake -c bioconda -c conda-forge snakemake=7.32.4 pulp=2.7.0`.
   `run.sh` then checks that the active Snakemake matches this file.
2. Clone the repository and place the paired-end reads under `data_dir` as
   `{sample}_1.fq.gz` / `{sample}_2.fq.gz`.
3. Provide the external assets and reference databases:
   - external tools: FastQ-Screen and SPAdes at absolute paths
     (`software.fastq_screen`, `software.spades`);
   - reference databases: Bowtie2 index (both groups), CoverM reference FASTA,
     fastq_screen configs, GTDB-Tk data (`gtdbtk_data`), CheckM data
     (`checkm_data`, pre-populated with `checkm data setRoot <dir>`).
4. Copy the config template and fill the 12 machine-specific paths plus
   `metawrap.local_base` (must be a **short, node-local** scratch dir):
   ```bash
   cp config.example.yaml config.yaml
   ```
5. Create the nine tool environments (`mamba env create -f envs/<name>.yml`) or
   let `run.sh` create the missing ones on first use. The Snakemake driver env is
   separate and not created by `run.sh`.
6. Apply the metaWRAP guard's one-time environment patch (see
   [metaWRAP arm](#metawrap-arm)):
   ```bash
   bash tools/mw_guard.sh patch-l3 && bash tools/mw_guard.sh verify-l3
   ```
7. Validate the DAG (`bash run.sh -n`), then run with `bash run.sh ...`; resume
   with `--rerun-incomplete`.

The `check_inputs` rule verifies that all configured input paths exist before
any long job starts (the guard is checked earlier, at DAG parse time).

### Regression baseline

`tools/mini/` holds a tiny end-to-end fixture (2 genomes, one sample, one binner arm
plus the metaWRAP arm). `tools/mini/validate_run.sh` clones the repo, generates the
fixture, fingerprints the DAG and asserts the outputs:

```bash
bash tools/mini/validate_run.sh --repo <repo> --dry-only   # [fingerprint] MATCH, total=28
bash tools/mini/validate_run.sh --repo <repo> --cores 16   # 31 PASS / 0 WARN / 0 FAIL + VALIDATE: OK
```

The validated run took **1h39m19s** (28 jobs): `metawrap_bins` 55m48s,
GTDB-Tk ×2 29m11s, CheckM ×2 10m44s, SPAdes 85s, the rest ~2 min. This covers
only one arm/sample; other binners, `PCR_assemble`, the `fastqs` and `coverm`
methods, multiple groups/samples in parallel and real-scale resources are **not**
covered.

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
| `PplacerException` from `gtdbtk_classify`, and `pplacer.class_level.<domain>.out` stops at `Preparing the edges for baseball...` | pplacer was killed for memory: the host accounts `PARENT x (1 + pplacer_cpus)` (~2 x 55-61 GB for the pinned ref data) | re-run the single job first; keep one GTDB-Tk job at a time (`--set-threads gtdbtk_classify=16`), keep the rule's `--pplacer_cpus 1`/`--scratch_dir`, raise `resources.gtdbtk_mem_mb` |
| A job that failed on one node succeeds unchanged on another node | environment-induced kill (memory pressure or a cgroup limit), not a code defect | compare the two nodes' free RAM before touching parameters; run the matrix on a compute node |

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
| CheckM | 1.0.12 | `metawrap-env` | metaWRAP arm's internal preselection only |
| metaWRAP | 1.3.2 | `metawrap-env` | Uritskiy et al., 2018 |
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
output/checkm/**/storage/bin_stats_ext.tsv   (comp/con + genome stats, from CheckM 1.2.3)
output/checkm/**/results.tsv                 (checkm qa -o 2: strain heterogeneity, marker lineage)
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

> `strain_heterogeneity` / `marker_lineage` source. CheckM's
> `storage/bin_stats_ext.tsv` contains neither a strain-heterogeneity key nor a
> title-case marker-lineage key (this holds for CheckM 1.0.12 *and* 1.2.3; the
> qa table's `Marker lineage` header is title-case, the ext dict's key is
> lowercase `marker lineage`). The `checkm_qa` rule therefore re-derives
> `checkm qa -o 2 --tab_table` into `output/checkm/**/results.tsv`, and
> `tools/harvest_mags.py` prefers that table for `strain_heterogeneity` (and
> falls back to it for `marker_lineage`), so the documented winner score's
> `con*(strh/100)` term is actually evaluated. Adding this rule makes one cheap
> `checkm qa` job per combination run on the next invocation (it re-parses
> existing HMMER output; it does not re-run `checkm lineage_wf`).

Per-pool behaviour: 0 MIMAG → `skipped_empty` (no dRep call, empty
`data_tables/*.csv` placeholders); 1 MIMAG → `passthrough_single` (copied into
`dereplicated_genomes/`, no dRep call); ≥2 → `dereplicated`. Each dRep run uses a
private, freshly wiped work directory. Two checks guard the layer: within a pool,
the manifest, the staged files, the `genomeInfo.csv` rows and the `g_source.list`
lines must all have the same size; and after a dereplicated run the log must
report `100.00% of genomes passed checkM filtering`, which is what the shared
cutoffs imply, so a lower number means the gate and dRep disagreed.
`tools/drep_stage.py` implements staging, verification, summaries, the funnel
table and the L9 join.

> dRep work directories.
> (1) dRep's `WorkDirectory.overwrite` is hardcoded to `True`, so dRep itself
> does not protect an existing work directory; the per-pool directory is private
> and wiped before each run, and dRep refuses a `-g` list when `Bdb` already
> exists (`assert wd.hasDb("Bdb") == False`). Do not point `dRep` at a shared or
> already-populated work directory.
> (2) `dRep compare` overwrites `Cdb/Mdb/Ndb/Widb` but not `Bdb/Wdb`, so a
> directory it touches holds tables from two different runs. It must not share
> the cross-workflow output directory; if `compare` is ever needed, drop `-g` so
> it reuses the existing `Bdb`, or use a fresh directory.

> Why `fastANI`. The originally published run used dRep's 3.6.2 default,
> which is `fastANI` (the `ANImf = (DEFAULT)` line in `dRep -h` is stale). The
> workflow passes it explicitly so the value cannot drift with the installed
> version. Switch `drep.s_algorithm` (and record it in Methods) if ANImf is
> required; it aligns whole genomes with `nucmer` and is far slower.

> Regenerating the dRep layer. Editing any rule/script here does not
> re-trigger jobs by default (`--rerun-triggers=mtime`); run once with `-F`, or
> delete the affected `derep.done` / `output/drep_*` trees. dRep's work directory
> is wiped and rebuilt per pool, so nothing stale survives a rerun.

## Citation

If you use this workflow, please cite it. Machine-readable metadata are in
[`CITATION.cff`](CITATION.cff).

> TODO (maintainer): fill the remaining placeholders in `CITATION.cff`
> (`version` is `1.2.2` and `date-released` is set; `repository-code`, `doi`,
> affiliation and ORCID still need to be added), and add the preferred citation
> (paper/DOI) here once one exists.

## License

[GPL v3](LICENSE)
