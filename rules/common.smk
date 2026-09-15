import os
from snakemake.io import glob_wildcards

SAMPLES = sorted(
    glob_wildcards(os.path.join(config["data_dir"], "{sample}_1.fq.gz")).sample
)

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

def binner_raw_genome_dir(wc):
    base = f"output/binning/{wc.binner}/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}"
    mapping = {
        "unitem": f"{base}/consensus/bins",
        "comebin": f"{base}/comebin_res/comebin_res_bins",
        "metadecoder": base,
        "semibin2_single": f"{base}/output_bins",
    }
    return mapping[wc.binner]

def binner_prepared_dir(wc):
    return f"output/bins_prepared/{wc.binner}/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}"

def binner_prepared_done(wc):
    return f"{binner_prepared_dir(wc)}/prepare.done"

def binner_bins_gzipped(wc):
    return wc.binner in ("unitem", "semibin2_single")

def binner_extension(wc):
    mapping = {
        "unitem": "fna",
        "comebin": "fa",
        "metadecoder": "fasta",
        "semibin2_single": "fa",
    }
    return mapping[wc.binner]

_BT2_SUFFIX = (".1.bt2", ".1.bt2l")

def bowtie2_index_probe(prefix):
    """Return whichever bowtie2 index flavour exists (>4 GB references get .bt2l)."""
    for suffix in _BT2_SUFFIX:
        if os.path.exists(prefix + suffix):
            return prefix + suffix
    return prefix + _BT2_SUFFIX[0]
