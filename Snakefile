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
include: "rules/summary.smk"
include: "rules/drep.smk"


rule all:
    input:
        # fail fast (also when SAMPLES is empty) and gate QC reference data
        "output/check_inputs.done",
        "output/check_qc_inputs.done",

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
        ),

        # MIMAG gate 1 (thin layer over the shared harvest parser)
        "output/summary/mimag_bins.tsv",
        "output/summary/mimag_summary.tsv",
        "output/summary/mimag_excluded.tsv",

        # dRep layer: per-workflow dereplication + cross-workflow + L9
        expand(
            "output/drep_per_workflow/{binner}/{treat}/{method}/{group}/{sample}/derep.done",
            binner=BINNERS,
            treat=ASSEMBLY_TREATS,
            method=FILTER_METHODS,
            group=GROUPS,
            sample=SAMPLES
        ),
        "output/drep_per_workflow/per_workflow.tsv",
        "output/drep_cross_workflow/data_tables/Bdb.csv",
        "output/drep_cross_workflow/data_tables/Cdb.csv",
        "output/drep_cross_workflow/data_tables/Wdb.csv",
        "output/summary/final_mags.tsv"
