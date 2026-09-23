import os
from snakemake.utils import min_version

min_version("7.0")

configfile: "config.yaml"

include: "rules/common.smk"

GROUPS = config["groups"]
FILTER_METHODS = config["filter_methods"]
ASSEMBLY_TREATS = config["assembly_treats"]
BINNERS = config["binners"]

# The metaWRAP arm is a binner whose two assembly treatments are decoupled from
# `assembly_treats`: it can run on a subset (metawrap_treats). Every other
# binner still runs on every assembly treatment.
METAWRAP_BINNER = "metawrap"
MW_BINNERS = [b for b in BINNERS if b == METAWRAP_BINNER]
OTHER_BINNERS = [b for b in BINNERS if b != METAWRAP_BINNER]

# Absent key => the metaWRAP arm runs on every assembly treatment. A
# present-but-empty list is rejected: it would otherwise silently mean "run
# everything", the opposite of what an explicit empty list reads like.
if "metawrap_treats" in config:
    MW_TREATS = list(config["metawrap_treats"])
    if not MW_TREATS:
        raise ValueError(
            "config.yaml: 'metawrap_treats' is empty; remove the key to run the "
            "metaWRAP arm on every assembly treatment.")
else:
    MW_TREATS = list(ASSEMBLY_TREATS)
if not set(MW_TREATS) <= set(ASSEMBLY_TREATS):
    raise ValueError(
        "metawrap_treats {0} must be a subset of assembly_treats {1}".format(
            sorted(MW_TREATS), sorted(ASSEMBLY_TREATS)))

# (binner, treat) pairs. A single expand() cannot express "these binners on all
# treats, that binner on a subset" because Snakemake's zip() requires every
# wildcard list to have the same length. The cartesian product is therefore
# enumerated once here and reused by the rule files via COMBOS / combo_paths().
BINNER_TREAT_PAIRS = (
    [(b, t) for b in OTHER_BINNERS for t in ASSEMBLY_TREATS]
    + [(b, t) for b in MW_BINNERS for t in MW_TREATS]
)
COMBOS = [
    (b, t, m, g, s)
    for (b, t) in BINNER_TREAT_PAIRS
    for m in FILTER_METHODS
    for g in GROUPS
    for s in SAMPLES
]


def combo_paths(prefix, leaf):
    """All paths "output/<prefix>/<binner>/<treat>/<method>/<group>/<sample>/<leaf>"."""
    return [
        "output/{0}/{1}/{2}/{3}/{4}/{5}/{6}".format(prefix, b, t, m, g, s, leaf)
        for (b, t, m, g, s) in COMBOS
    ]


include: "rules/preprocess.smk"
include: "rules/assembly.smk"
include: "rules/binning.smk"
include: "rules/metawrap.smk"
include: "rules/summary.smk"
include: "rules/drep.smk"


# Single flat target list for `rule all`. Built as one expression (no `*`
# unpacking) because Snakemake's rule parser is not a plain Python parser.
ALL_TARGETS = (
    # fail fast (also when SAMPLES is empty) and gate QC reference data
    ["output/check_inputs.done", "output/check_qc_inputs.done"]
    # final assembly contigs > 2000bp
    + expand(
        "output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta",
        treat=ASSEMBLY_TREATS,
        method=FILTER_METHODS,
        group=GROUPS,
        sample=SAMPLES
    )
    # final sorted BAMs mapped back to contigs
    + expand(
        "output/assemble/{treat}/{method}/{group}/{sample}/{sample}_sorted.bam",
        treat=ASSEMBLY_TREATS,
        method=FILTER_METHODS,
        group=GROUPS,
        sample=SAMPLES
    )
    # RefineM / CheckM / GTDB-Tk completion markers
    + combo_paths("refinem", "refinem.done")
    + combo_paths("checkm", "checkm.done")
    + combo_paths("gtdb", "gtdb.done")
    # MIMAG gate 1 (thin layer over the shared harvest parser)
    + [
        "output/summary/mimag_bins.tsv",
        "output/summary/mimag_summary.tsv",
        "output/summary/mimag_excluded.tsv",
    ]
    # dRep layer: per-workflow dereplication + cross-workflow + L9
    + combo_paths("drep_per_workflow", "derep.done")
    + [
        "output/drep_per_workflow/per_workflow.tsv",
        "output/drep_cross_workflow/data_tables/Bdb.csv",
        "output/drep_cross_workflow/data_tables/Cdb.csv",
        "output/drep_cross_workflow/data_tables/Wdb.csv",
        "output/summary/final_mags.tsv",
    ]
)


rule all:
    input:
        ALL_TARGETS
