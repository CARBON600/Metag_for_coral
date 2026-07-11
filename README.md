# Standardizing a Bioinformatic Analysis Pipeline for Coral Metagenomics
> Genome-resolved benchmarking of wet-lab protocols and bioinformatic workflow components to improve MAG recovery from coral holobiont metagenomes.
## Overview
Coral holobiont metagenomes are often dominated by coral host and Symbiodiniaceae DNA, limiting genome-resolved recovery of microbial members.  
This project benchmarks integrated **wet lab + bioinformatic** workflow combinations to optimize **metagenome-assembled genome (MAG)** recovery, quality and downstream functional interpretability.
The workflow is written in [Snakemake](https://snakemake.github.io/) and supports
conda-based dependency management. It is designed for Linux HPC clusters (with or
without Slurm).

## Pipeline Summary

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
                          Kaiju taxonomy
                         (on raw reads)
```

---

## Requirements

| Category | Requirement |
|----------|-------------|
| **OS** | Linux (conda environments pinned to `linux-64`) |
| **Snakemake** | ≥ 7.0 |
| **Conda** | Miniconda or Anaconda |
| **Storage** | ≥ 500 GB recommended |

### Reference Databases

| Resource | Purpose | Used by |
|----------|---------|---------|
| Bowtie2 index | Host-read removal | `bt2` method |
| CoverM reference FASTA | Host-read removal | `coverm` method |
| FastQ-Screen config file | Multi-genome host screening | `fastqs` method |
| Kaiju DB | Taxonomic profiling | `kaiju_raw` (`nodes.dmp`, `names.dmp`, `.fmi`) |

---

## Quick Start

### 1. Clone and install

```bash
git clone <repo-url>
cd Metag_for_coral-main

# Let Snakemake handle environments automatically (recommended), or:
conda env create -f envs/metag.yml
conda env create -f envs/binning.yml
```

### 2. Edit `config.yaml`

Set your paths, reference databases, and choose which tools to run:

```yaml
data_dir: "/path/to/raw/data"      # contains *_1.fq.gz / *_2.fq.gz
tmpdir:   "/path/to/tmp"

# ── Workflow components (comment out to skip) ──
filter_methods:
  - bt2
  - coverm
  - fastqs

assembly_treats:
  - control_assemble       # no PCR dedup
  - PCR_assemble           # with PCR dedup

binners:
  - unitem
  - comebin
  - metadecoder
  - semibin2_single

# ── Reference databases ──
bowtie2_index:
  control: "/path/to/control/control"
  exp:     "/path/to/exp/exp"

coverm_ref:
  control: "/path/to/control.fna"
  exp:     "/path/to/ref_Acropora_allSymb_comb.fa"

fastq_screen_conf:
  control: "/path/to/control.conf"
  exp:     "/path/to/exp.conf"

kaiju:
  nodes: "/path/to/kaiju/nodes.dmp"
  names: "/path/to/kaiju/names.dmp"
  db:    "/path/to/kaiju/kaiju_db_nr_euk.fmi"

# ── Threads (per rule) ──
threads:
  fastq_screen: 16
  fastp: 16
  bowtie2_map: 32
  # ... (see config.yaml for full list)
```

### 3. Dry-run

```bash
snakemake --use-conda -n
```

### 4. Run

```bash
# Single node
snakemake --use-conda --cores 128 --keep-going

# Slurm cluster
snakemake --use-conda --jobs 50 --keep-going \
    --cluster "sbatch --cpus-per-task={threads} --mem=64G --time=48:00:00" \
    --latency-wait 120
```

### 5. Resume after failure

```bash
snakemake --use-conda --cores 128 --keep-going --rerun-incomplete
```

### 6. Visualise DAG

```bash
snakemake --use-conda -n --dag | dot -Tpng > dag.png
```

---

## Selecting Subsets of Tools

Edit the three lists in `config.yaml` to control which components run:

```yaml
# Compare only bt2 vs coverm, two binners, skip PCR dedup
filter_methods:
  - bt2
  - coverm

assembly_treats:
  - control_assemble

binners:
  - unitem
  - semibin2_single
```

Snakemake's demand-driven execution ensures only the required rules run —
irrelevant dependency chains are automatically pruned.

---

## Input Data

Place paired-end FASTQ files under `data_dir`:

```
{SAMPLE}_1.fq.gz
{SAMPLE}_2.fq.gz
```

Samples are discovered by globbing `{sample}_1.fq.gz`. Each sample is processed
through all `groups` using group-specific reference databases.

---

## Output Structure

```
output/
├── fastqs/                 FastQ-Screen + fastp
│   └── {group}/{sample}/
├── bt2/                    Bowtie2 mapped BAM
│   └── {group}/
├── bt2_unmapped/           Bowtie2 unmapped (-f 4)
│   └── {group}/
├── coverm/                 CoverM mapped BAM
│   └── {group}/
├── coverm_filtered/        CoverM inverse-filtered BAM
│   └── {group}/
├── fq4dep/                 Filtered FASTQ for assembly
│   ├── bt2/{group}/{sample}/
│   └── coverm/{group}/{sample}/
├── PCR_free/               MEGAHIT pre-assembly + remapping
│   └── {method}/{group}/{sample}/
├── PCR_done/               PCR-deduplicated reads
│   └── {method}/{group}/{sample}/
├── assemble/               SPAdes contigs + filtered + BAM
│   └── {treat}/{method}/{group}/{sample}/
├── binning/                Binner outputs
│   └── {binner}/{treat}/{method}/{group}/{sample}/
├── refinem/                RefineM-filtered bins
│   └── {binner}/{treat}/{method}/{group}/{sample}/
├── checkm/                 CheckM results
│   └── {binner}/{treat}/{method}/{group}/{sample}/
├── gtdb/                   GTDB-Tk results
│   └── {binner}/{treat}/{method}/{group}/{sample}/
└── kaiju/                  Kaiju taxonomic summaries
    └── raw/{sample}_kaiju_summary.tsv
```

### Wildcard Key

| Wildcard | Values |
|----------|--------|
| `{group}` | `control`, `exp` |
| `{method}` | `bt2`, `coverm`, `fastqs` |
| `{treat}` | `control_assemble`, `PCR_assemble` |
| `{binner}` | `unitem`, `comebin`, `metadecoder`, `semibin2_single` |

---

## Workflow Details

### Stage 1 — Host-Read Removal

| Method | Tool(s) | Strategy |
|--------|---------|----------|
| `bt2` | Bowtie2 | Map to host reference → keep unmapped (`samtools view -f 4`) → FASTQ |
| `coverm` | CoverM | Map to host reference → inverse filter (≥75% aligned, ≥95% identity) → FASTQ |
| `fastqs` | FastQ-Screen + fastp | Screen against multi-genome config (`--nohits --tag`) → fastp QC |

### Stage 2 — Assembly

| Treatment | Input | PCR Dedup |
|-----------|-------|-----------|
| `control_assemble` | Filtered reads | No |
| `PCR_assemble` | Reads after MEGAHIT pre-assembly + `samtools markdup` | Yes |

MEGAHIT pre-assembly (`--min-contig-len 250`) generates contigs used solely as
a reference for PCR duplicate marking. **SPAdes** (`--meta`,
`-k 21,33,55,77,99,127`) is the final assembler for both treatments. Contigs
shorter than 2 000 bp are discarded before binning.

### Stage 3 — Binning

| Tool | Approach |
|------|----------|
| **UniT** | Ensemble of MetaBAT2 + MaxBin2; consensus binning |
| **COMEBin** | Contrastive multi-view representation learning |
| **MetaDecoder** | Two-layer DPGMM + k-mer frequency model (`--disable_gpu`) |
| **SemiBin2** | Self-supervised contrastive learning (`single_easy_bin --self-supervised`) |

### Stage 4 — Quality Assessment

| Tool | Function |
|------|----------|
| **RefineM** | Scaffold statistics → outlier detection → bin filtering |
| **CheckM** | Lineage-specific completeness & contamination (`lineage_wf`) |
| **GTDB-Tk** | Taxonomic classification (`classify_wf --skip_ani_screen`) |

### Stage 5 — Taxonomic Profiling

| Tool | Input | Output |
|------|-------|--------|
| **Kaiju** | Raw reads | Phylum-level summary table (`kaiju2table -r phylum`) |

---

## Design Notes

- The `rule all` target covers the full factorial design:  
  `filter_methods × assembly_treats × groups × binners × samples`
- `--keep-going` ensures one sample's failure does not block others.
- Intermediate large BAM files from PCR deduplication are automatically cleaned.
- Conda environments are split into `metag` (preprocess + assembly) and
  `binning` (binning + QC) to avoid dependency conflicts between tool versions.
- Output integrity is checked with `test -s` after critical steps (assembly,
  remapping, FASTQ conversion).

---

## Tools

| Tool | Version | Reference |
|------|---------|-----------|
| MEGAHIT | 1.2.9 | Li et al., 2015 |
| SPAdes | meta | Bankevich et al., 2012 |
| Bowtie2 | 2.5.4 | Langmead & Salzberg, 2012 |
| CoverM | 0.7.0 | Aroney et al., 2025 |
| FastQ-Screen | – | Wingett & Andrews, 2018 |
| fastp | – | Chen et al., 2018 |
| samtools | 1.18 | Danecek et al., 2021 |
| UniT | 1.0.2 | Parks et al. |
| COMEBin | 1.0.4 | Wang et al., 2024 |
| MetaDecoder | 1.2.1 | Liu et al., 2022 |
| SemiBin2 | 2.2.0 | Pan et al., 2023 |
| RefineM | 0.1.2 | Parks et al. |
| CheckM | 1.2.3 | Parks et al., 2015 |
| GTDB-Tk | – | Chaumeil et al., 2020 |
| Kaiju | – | Menzel et al., 2016 |

---

## License

[GPL v3](LICENSE)
