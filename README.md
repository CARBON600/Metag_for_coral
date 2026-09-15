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

Dereplication (dRep), functional annotation (Prokka/KofamScan) and the
statistical comparison of MAG sets are performed downstream in a separate
notebook and are not part of this workflow.

## Pipeline

```
Raw reads
  → host removal (bt2 | coverm | fastqs)
  → MEGAHIT pre-assembly → samtools markdup → SPAdes   (PCR_assemble)
  → SPAdes                                             (control_assemble)
  → contigs ≥ 2000 bp → Bowtie2 re-mapping
  → binning (UniT | COMEBin | MetaDecoder | SemiBin2)
  → prepare_bins → RefineM → CheckM / GTDB-Tk
```

Kaiju taxonomic profiling is disabled.

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
```

Because Snakemake reuses these named environments, editing a specification does
not rebuild the environment automatically; refresh it explicitly with
`conda env update -n <name> -f envs/<name>.yml` (or remove and recreate it).

| Environment | Steps |
|-------------|-------|
| `metag` | host removal (`bt2`, `coverm`, `fastqs`), fastp, PCR dedup, SPAdes, contig filtering/remapping |
| `megahit` | MEGAHIT pre-assembly |
| `binning` | UniT, RefineM, CheckM, bin preparation |
| `comebin_env` | COMEBin |
| `metadecoder` | MetaDecoder |
| `SemiBin` | SemiBin2 |
| `gtdbtk-2.3.2` | GTDB-Tk |

### External tools

FastQ-Screen and SPAdes are invoked from absolute paths configured in
`config.yaml` (`software:`), with the original notebook paths as defaults:

| Tool | Config key |
|------|------------|
| FastQ-Screen 0.14.1 | `software.fastq_screen` |
| SPAdes 3.15.5 | `software.spades` |

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

`run.sh` checks for the seven named environments and creates any that are
missing from `envs/*.yml`, then runs Snakemake with the conda frontend set to
`conda`. Any arguments are passed through to Snakemake:

```bash
# Validate the DAG
bash run.sh -n

# Single node
bash run.sh --cores 128 --keep-going

# Slurm (Snakemake 7.x)
bash run.sh --jobs 50 --keep-going \
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
├── binning/                Binner outputs
├── bins_prepared/          Decompressed / cleaned bins for QC
├── refinem/                RefineM-filtered bins
├── checkm/                 CheckM results
├── gtdb/                   GTDB-Tk results
└── logs/                   Per-rule logs
```

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
| COMEBin | `comebin_env` | Contrastive multi-view representation learning |
| MetaDecoder | `metadecoder` | Two-layer DPGMM + k-mer frequency model |
| SemiBin2 | `SemiBin` | Self-supervised contrastive learning (`single_easy_bin --self-supervised`) |

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

## Deployment on a new machine

1. Install Miniconda/Mamba and Snakemake ≥ 7 (`mamba create -n snakemake -c bioconda -c conda-forge snakemake=7.32.4`).
2. Clone the repository and place the paired-end reads under `data_dir`.
3. Provide the external tools and reference databases, then update every path in `config.yaml` (`data_dir`, `tmpdir`, `bowtie2_index`, `coverm_ref`, `fastq_screen_conf`, `software`, `gtdbtk_data`, `checkm_data`) and adjust `threads`.
4. Create the seven environments (`mamba env create -f envs/<name>.yml`) or let `run.sh` create the missing ones on first use.
5. Validate the DAG (`bash run.sh -n`), then run with `bash run.sh ...`; resume with `--rerun-incomplete`.

The `check_inputs` rule verifies that all configured input paths exist before
any long job starts.

### Troubleshooting

| Symptom | Probable cause | Remedy |
|---|---|---|
| `Nothing to be done` | `data_dir` contains no `{sample}_1.fq.gz`, so `SAMPLES` is empty | correct `data_dir` |
| `FileNotFoundError` for a reference or tool | a path in `config.yaml` was not updated | update `config.yaml` |
| `check_inputs` fails | one of the configured paths does not exist | read the log under `logs/check_inputs.log` |
| CheckM reports the data folder is unset | `CHECKM_DATA_PATH` not configured | set `checkm_data`, or run `checkm data setRoot <dir>` once |
| `The 'mamba' command is not available` | Snakemake's default conda frontend is `mamba` | use `run.sh` (sets `--conda-frontend conda`) or install mamba |
| `Conda must be version 24.7.1 or later` | a Snakemake 8.x install requires a newer conda | pin Snakemake 7.32.4 (see above) or upgrade conda |
| Conda environment creation fails | `defaults` ToS, network, or a missing build | accept the ToS, use a mirror, or relax pins |
| SPAdes exits with a Python error | interpreter mismatch of the external script | invoke via `conda run -n <env>` or install SPAdes into the environment |
| Disk exhaustion | `tmpdir` or `output/` under-provisioned | relocate to larger storage |

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
| Kaiju | – | disabled | Menzel et al., 2016 |

## License

[GPL v3](LICENSE)
