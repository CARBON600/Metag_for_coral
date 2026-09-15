import os
from snakemake.utils import min_version

min_version("7.0")

configfile: "config.yaml"

include: "rules/common.smk"

GROUPS = config["groups"]
FILTER_METHODS = config["filter_methods"]
ASSEMBLY_TREATS = config["assembly_treats"]
BINNERS = config["binners"]

include: "rules/preprocess.smk"
include: "rules/assembly.smk"
include: "rules/binning.smk"


rule all:
    input:
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
