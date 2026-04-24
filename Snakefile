import os
from snakemake.io import glob_wildcards

configfile: "config.yaml"

# -------------------------
# Sample discovery
# -------------------------
SAMPLES = sorted(
    glob_wildcards(os.path.join(config["data_dir"], "{sample}_1.fq.gz")).sample
)

GROUPS = config["groups"]                # ["control", "exp"]
FILTER_METHODS = ["bt2", "coverm", "fastqs"]
ASSEMBLY_TREATS = ["control_assemble", "PCR_assemble"]
BINNERS = ["unitem", "comebin", "metadecoder", "semibin2_single"]

def raw_r1(wc):
    return os.path.join(config["data_dir"], f"{wc.sample}_1.fq.gz")
def raw_r2(wc):
    return os.path.join(config["data_dir"], f"{wc.sample}_2.fq.gz")

def filtered_r1(wc):
    if wc.method == "fastqs":
        return f"output/fastqs/{wc.group}/{wc.sample}/{wc.sample}_1_cleaned.fastq.gz"
    elif wc.method == "bt2":
        return f"output/fq4dep/bt2/{wc.group}/{wc.sample}_1.fq.gz"
    elif wc.method == "coverm":
        return f"output/fq4dep/coverm/{wc.group}/{wc.sample}_1.fq.gz"
    raise ValueError(f"Unknown method: {wc.method}")
def filtered_r2(wc):
    if wc.method == "fastqs":
        return f"output/fastqs/{wc.group}/{wc.sample}/{wc.sample}_2_cleaned.fastq.gz"
    elif wc.method == "bt2":
        return f"output/fq4dep/bt2/{wc.group}/{wc.sample}_2.fq.gz"
    elif wc.method == "coverm":
        return f"output/fq4dep/coverm/{wc.group}/{wc.sample}_2.fq.gz"
    raise ValueError(f"Unknown method: {wc.method}")

def assembly_reads_r1(wc):
    if wc.treat == "control_assemble":
        return filtered_r1(wc)
    elif wc.treat == "PCR_assemble":
        return f"output/PCR_done/{wc.method}/{wc.group}/{wc.sample}_1.fq.gz"
    raise ValueError(f"Unknown treat: {wc.treat}")
def assembly_reads_r2(wc):
    if wc.treat == "control_assemble":
        return filtered_r2(wc)
    elif wc.treat == "PCR_assemble":
        return f"output/PCR_done/{wc.method}/{wc.group}/{wc.sample}_2.fq.gz"
    raise ValueError(f"Unknown treat: {wc.treat}")


def megahit_contigs(wc):
    return f"output/megahit_pre/{wc.method}/{wc.group}/{wc.sample}/{wc.sample}.contigs.fa"

def final_contigs(wc):
    return f"output/assemble/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}/contigs.fasta"
def final_contigs_r2000(wc):
    return f"output/assemble/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}/contigs_r2000bp.fasta"
def final_sorted_bam(wc):
    return f"output/assemble/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}/{wc.sample}_sorted.bam"

def binner_done(wc):
    if wc.binner == "unitem":
        return f"output/binning/unitem/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}/consensus.done"
    elif wc.binner == "comebin":
        return f"output/binning/comebin/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}/comebin.done"
    elif wc.binner == "metadecoder":
        return f"output/binning/metadecoder/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}/cluster.done"
    elif wc.binner == "semibin2_single":
        return f"output/binning/semibin2_single/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}/semibin.done"
    raise ValueError(f"Unknown binner: {wc.binner}")

def binner_genome_dir(wc):
    base = f"output/binning/{wc.binner}/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}"
    mapping = {
        "unitem": f"{base}/consensus/bins",
        "comebin": f"{base}/comebin_res/comebin_res_bins",
        "metadecoder": base,
        "semibin2_single": f"{base}/output_bins",
    }
    return mapping[wc.binner]

def binner_extension(wc):
    mapping = {
        "unitem": "fna",
        "comebin": "fa",
        "metadecoder": "fasta",
        "semibin2_single": "fa",
    }
    return mapping[wc.binner]

include: "rules/preprocess.smk"
include: "rules/assembly.smk"
include: "rules/binning.smk"


rule all:
    input:
        # raw read taxonomic profiling
        expand("output/kaiju/raw/{sample}_kaiju_summary.tsv", sample=SAMPLES),

        # final assembly contigs > 2000bp
        expand(
            "output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta",
            treat=ASSEMBLY_TREATS,
            method=FILTER_METHODS,
            group=GROUPS,
            sample=SAMPLES
        ),

        # final sorted BAMs mapped back to contigs
        expand(
            "output/assemble/{treat}/{method}/{group}/{sample}/{sample}_sorted.bam",
            treat=ASSEMBLY_TREATS,
            method=FILTER_METHODS,
            group=GROUPS,
            sample=SAMPLES
        ),

        # RefineM / CheckM / GTDB-Tk completion markers
        expand(
            "output/refinem/{binner}/{treat}/{method}/{group}/{sample}/refinem.done",
            binner=BINNERS,
            treat=ASSEMBLY_TREATS,
            method=FILTER_METHODS,
            group=GROUPS,
            sample=SAMPLES
        ),
        expand(
            "output/checkm/{binner}/{treat}/{method}/{group}/{sample}/checkm.done",
            binner=BINNERS,
            treat=ASSEMBLY_TREATS,
            method=FILTER_METHODS,
            group=GROUPS,
            sample=SAMPLES
        ),
        expand(
            "output/gtdb/{binner}/{treat}/{method}/{group}/{sample}/gtdb.done",
            binner=BINNERS,
            treat=ASSEMBLY_TREATS,
            method=FILTER_METHODS,
            group=GROUPS,
            sample=SAMPLES
        )