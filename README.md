# Standardizing a Bioinformatic Analysis Pipeline for Coral Metagenomics
> Genome-resolved benchmarking of wet-lab protocols and bioinformatic workflow components to improve MAG recovery from coral holobiont metagenomes.
## Overview
Coral holobiont metagenomes are often dominated by coral host and Symbiodiniaceae DNA, limiting genome-resolved recovery of microbial members.  
This project benchmarks integrated **wet lab + bioinformatic** workflow combinations to optimize **metagenome-assembled genome (MAG)** recovery, quality and downstream functional interpretability.
The workflow is written in [Snakemake](https://snakemake.github.io/) and supports
conda-based dependency management. It is designed for Linux HPC clusters (with or
without Slurm).

## Pipeline summary

Raw reads (paired-end)
 │
 ├── FastQ-Screen (+ fastp) ────┐
 ├── Bowtie2 → unmapped → FASTQ ─┤
 └── CoverM → inverse filter ────┤
                                  ▼
                        Filtered reads
                                  │
                    ┌─────────────┴──────────────┐
                    │                            │
              MEGAHIT pre-assembly          SPAdes meta assembly
           (reference for PCR dedup)        (control – no dedup)
                    │                            │
              PCR deduplication                  │
              (samtools markdup)                 │
                    │                            │
              SPAdes meta assembly               │
              (PCR-assemble)                     │
                    │                            │
                    └──────────┬─────────────────┘
                               ▼
                    contigs ≥ 2000 bp
                               │
                    Bowtie2 read re-mapping
                               │
               ┌───────────────┼───────────────┐
               ▼               ▼               ▼               ▼
            UniT           COMEbin        MetaDecoder     SemiBin2
               │               │               │               │
               └───────────────┴───────────────┴───────────────┘
                               ▼
                    RefineM (filter bins)
                               │
                    ┌──────────┴──────────┐
                    ▼                     ▼
              CheckM lineage_wf    GTDB-Tk classify_wf
                    │                     │
                    └──────────┬──────────┘
                               ▼
                        Kaiju taxonomy

## Requirements

- **OS:** Linux (the conda environments are pinned to linux-64 packages)
- **Snakemake ≥ 7.0**
- **Conda** (Miniconda / Anaconda)

### Databases

The following reference databases must be prepared before running:

| Resource | Purpose | Notes |
|----------|---------|-------|
| Bowtie2 index | Host-read removal (bt2 method) | E.g. coral host + Symbiodiniaceae genomes |
| CoverM reference FASTA | Host-read removal (coverm method) | Same genomes, concatenated |
| FastQ-Screen config | Host-read screening (fastqs method) | Bowtie2-based multi-genome config file |
| Kaiju DB | Taxonomic profiling | `kaiju_db_nr_euk.fmi` + `nodes.dmp` + `names.dmp` |

## Installation

```bash
# Clone the repository
git clone <repo-url>
cd Metag_for_coral-main

# Create conda environments (or let Snakemake do it automatically)
conda env create -f envs/metag.yml
conda env create -f envs/binning.yml
Configuration
Edit config.yaml to match your environment:
# Paths
data_dir: "/path/to/raw/data"          # contains *_1.fq.gz and *_2.fq.gz
tmpdir:   "/path/to/tmp"               # CoverM temporary directory

# Sample groups (databases)
groups:
  - control
  - exp

# Workflow components to run (comment out any to skip)
filter_methods:
  - bt2
  - coverm
  - fastqs

assembly_treats:
  - control_assemble     # no PCR dedup
  - PCR_assemble         # with PCR dedup

binners:
  - unitem
  - comebin
  - metadecoder
  - semibin2_single

# Reference databases (update to your paths)
fastq_screen_conf:
  control: "/path/to/control.conf"
  exp:     "/path/to/exp.conf"

bowtie2_index:
  control: "/path/to/control/control"
  exp:     "/path/to/exp/exp"

coverm_ref:
  control: "/path/to/control.fna"
  exp:     "/path/to/ref_Acropora_allSymb_comb.fa"

kaiju:
  nodes: "/path/to/kaiju/nodes.dmp"
  names: "/path/to/kaiju/names.dmp"
  db:    "/path/to/kaiju/kaiju_db_nr_euk.fmi"

# Threads (adjust to your hardware)
threads:
  fastq_screen: 16
  fastp: 16
  bowtie2_map: 32
  # ...
Input data naming
Raw paired-end reads must be placed under data_dir with the naming convention:
{SAMPLE}_1.fq.gz
{SAMPLE}_2.fq.gz
Samples are discovered automatically by globbing. Group membership is defined
by groups in the config – each sample is processed through all groups using
group-specific reference databases.
Usage
1. Dry-run (check DAG)
snakemake --use-conda -n
2. Run on a single node
snakemake --use-conda --cores 128 --keep-going
3. Run on a Slurm cluster
snakemake --use-conda --jobs 50 --keep-going \
    --cluster "sbatch --cpus-per-task={threads} --mem=64G --time=48:00:00" \
    --latency-wait 120
4. Run only a subset
Edit config.yaml to disable unwanted components:
# Compare only bt2 vs coverm, only two binners, skip PCR dedup
filter_methods:
  - bt2
  - coverm

assembly_treats:
  - control_assemble

binners:
  - unitem
  - semibin2_single
5. Resume after failure
snakemake --use-conda --cores 128 --keep-going --rerun-incomplete
6. Visualise the DAG
snakemake --use-conda -n --dag | dot -Tpng > dag.png
Output structure
output/
├── fastqs/                    # FastQ-Screen + fastp output
│   └── {group}/{sample}/
├── bt2/                       # Bowtie2 mapped BAM
│   └── {group}/
├── bt2_unmapped/              # Bowtie2 unmapped BAM (-f 4)
│   └── {group}/
├── coverm/                    # CoverM mapped BAM
│   └── {group}/
├── coverm_filtered/           # CoverM inverse-filtered BAM
│   └── {group}/
├── fq4dep/                    # Filtered FASTQ for assembly
│   ├── bt2/{group}/{sample}/
│   └── coverm/{group}/{sample}/
├── PCR_free/                  # MEGAHIT pre-assembly + remapping
│   └── {method}/{group}/{sample}/
├── PCR_done/                  # PCR-deduplicated reads
│   └── {method}/{group}/{sample}/
├── assemble/                  # SPAdes assemblies + filtered contigs + BAM
│   └── {treat}/{method}/{group}/{sample}/
├── binning/                   # Binner outputs
│   └── {binner}/{treat}/{method}/{group}/{sample}/
├── refinem/                   # RefineM-filtered bins
│   └── {binner}/{treat}/{method}/{group}/{sample}/
├── checkm/                    # CheckM lineage_wf results
│   └── {binner}/{treat}/{method}/{group}/{sample}/
├── gtdb/                      # GTDB-Tk classify_wf results
│   └── {binner}/{treat}/{method}/{group}/{sample}/
└── kaiju/                     # Kaiju taxonomic summaries
    └── raw/{sample}_kaiju_summary.tsv
Path key
Variable	Values
{group}	control, exp
{method}	bt2, coverm, fastqs
{treat}	control_assemble, PCR_assemble
{binner}	unitem, comebin, metadecoder, semibin2_single
Workflow details
Host-read removal (3 parallel methods)
Method	Tool	Strategy
bt2	Bowtie2	Map → keep unmapped (samtools view -f 4) → FASTQ
coverm	CoverM	Map → inverse filter (min 75% aligned, 95% identity) → FASTQ
fastqs	FastQ-Screen + fastp	Screen against multi-genome config → fastp QC
Assembly (2 treatments)
Treatment	Input reads	PCR dedup
control_assemble	Filtered reads	No
PCR_assemble	Reads after MEGAHIT pre-assembly + samtools markdup	Yes
MEGAHIT pre-assembly (--min-contig-len 250) produces contigs used only as
a reference for PCR duplicate marking. SPAdes (--meta) is the final assembler
for both treatments.
Binning (4 tools)
Binner	Key features
UniT	Ensemble of MetaBAT2, MaxBin2; consensus binning
COMEBin	Contrastive multi-view representation learning
MetaDecoder	Two-layer DPGMM + k-mer frequency model
SemiBin2	Self-supervised contrastive learning (--self-supervised)
Quality assessment
1. RefineM — scaffold statistics, outlier detection, bin filtering
2. CheckM — lineage-specific completeness & contamination (lineage_wf)
3. GTDB-Tk — taxonomic classification (classify_wf --skip_ani_screen)
Taxonomic profiling
- Kaiju — run on raw reads; produces phylum-level summary tables for all
samples.
Design notes
- The rule all target covers the complete factorial design:
filter_methods × assembly_treats × groups × binners × samples.
- Missing intermediate files cause only the required rules to re-run.
- --keep-going ensures one sample failure does not block others.
- Intermediate large BAM files from PCR deduplication are automatically cleaned.
- Conda environments are split: metag for preprocess/assembly, binning for
binning/QC – avoiding dependency conflicts between newer and older tool
versions.
