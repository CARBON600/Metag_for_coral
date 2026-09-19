# dRep layer: dereplication within each workflow and across all of them, plus
# the L9 join. The cutoffs below are the mimag cutoffs again; the two copies are
# asserted equal here so they cannot drift apart.
import os

assert config["drep"]["min_completion"] == config["mimag"]["min_completeness"], \
    "drep.min_completion must equal mimag.min_completeness"
assert config["drep"]["max_contamination"] == config["mimag"]["max_contamination"], \
    "drep.max_contamination must equal mimag.max_contamination"

DREP_ENV = "drep"
DREP_COMP = config["mimag"]["min_completeness"]
DREP_CON = config["mimag"]["max_contamination"]
DREP_LEN = config["drep"]["min_length"]
DREP_PA = config["drep"]["p_ani"]
DREP_SA = config["drep"]["s_ani"]
DREP_NC = config["drep"]["cov_thresh"]
DREP_ALG = config["drep"]["s_algorithm"]


# One dRep run per (binner, treat, method, group, sample) pool.
rule drep_per_workflow:
    input:
        mimag="output/summary/mimag_bins.tsv"
    output:
        done="output/drep_per_workflow/{binner}/{treat}/{method}/{group}/{sample}/derep.done",
        summary="output/drep_per_workflow/{binner}/{treat}/{method}/{group}/{sample}/summary.tsv"
    conda:
        DREP_ENV
    threads: config["threads"]["drep_pw"]
    resources:
        mem_mb=config["resources"]["drep_pw_mem_mb"]
    log:
        "logs/drep_per_workflow/{binner}/{treat}/{method}/{group}/{sample}.log"
    params:
        wd=lambda wc: "output/drep_per_workflow/{0}/{1}/{2}/{3}/{4}".format(
            wc.binner, wc.treat, wc.method, wc.group, wc.sample),
        comp=DREP_COMP,
        con=DREP_CON,
        length=DREP_LEN,
        p_ani=DREP_PA,
        s_ani=DREP_SA,
        cov=DREP_NC,
        alg=DREP_ALG
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        python scripts/drep_stage.py stage --mode per_workflow \
          --mimag-bins {input.mimag} --output-dir output --workdir {params.wd} \
          --binner {wildcards.binner} --treat {wildcards.treat} \
          --method {wildcards.method} --group {wildcards.group} --sample {wildcards.sample}
        if [ "$(cat {params.wd}/status.txt)" = "dereplicated" ]; then
          dRep dereplicate {params.wd} -g {params.wd}/g_source.list \
            -comp {params.comp} -con {params.con} -l {params.length} \
            -pa {params.p_ani} -sa {params.s_ani} -nc {params.cov} \
            --S_algorithm {params.alg} \
            --genomeInfo {params.wd}/genomeInfo.csv -p {threads}
          python scripts/drep_stage.py verify --log {params.wd}/log/logger.log \
            --derep-dir {params.wd}/dereplicated_genomes
        fi
        python scripts/drep_stage.py summary --workdir {params.wd} \
          --mimag-summary output/summary/mimag_summary.tsv \
          --binner {wildcards.binner} --treat {wildcards.treat} \
          --method {wildcards.method} --group {wildcards.group} --sample {wildcards.sample}
        test -s {output.summary}
        touch {output.done}
        """


# One dRep run over every passing MAG from all workflows combined.
rule drep_cross_workflow:
    input:
        mimag="output/summary/mimag_bins.tsv"
    output:
        bdb="output/drep_cross_workflow/data_tables/Bdb.csv",
        cdb="output/drep_cross_workflow/data_tables/Cdb.csv",
        wdb="output/drep_cross_workflow/data_tables/Wdb.csv",
        manifest="output/drep_cross_workflow/manifest.tsv"
    conda:
        DREP_ENV
    threads: config["threads"]["drep_xw"]
    resources:
        mem_mb=config["resources"]["drep_xw_mem_mb"]
    log:
        "logs/drep_cross_workflow.log"
    params:
        # Taken from the declared output so the work directory follows it; a
        # literal would be a hardcoded prefix of the outputs below.
        wd=lambda wc, output: os.path.dirname(os.path.dirname(output[0])),
        comp=DREP_COMP,
        con=DREP_CON,
        length=DREP_LEN,
        p_ani=DREP_PA,
        s_ani=DREP_SA,
        cov=DREP_NC,
        alg=DREP_ALG
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        python scripts/drep_stage.py stage --mode cross_workflow \
          --mimag-bins {input.mimag} --output-dir output --workdir {params.wd}
        if [ "$(cat {params.wd}/status.txt)" = "dereplicated" ]; then
          dRep dereplicate {params.wd} -g {params.wd}/g_source.list \
            -comp {params.comp} -con {params.con} -l {params.length} \
            -pa {params.p_ani} -sa {params.s_ani} -nc {params.cov} \
            --S_algorithm {params.alg} \
            --genomeInfo {params.wd}/genomeInfo.csv -p {threads}
          python scripts/drep_stage.py verify --log {params.wd}/log/logger.log \
            --derep-dir {params.wd}/dereplicated_genomes
          test -s {output.bdb}
          test -s {output.cdb}
          test -s {output.wdb}
        fi
        """


# Funnel table: one row per pool, reconciled against the passing-MAG list.
rule drep_collect:
    input:
        done=expand(
            "output/drep_per_workflow/{binner}/{treat}/{method}/{group}/{sample}/derep.done",
            binner=BINNERS, treat=ASSEMBLY_TREATS,
            method=FILTER_METHODS, group=GROUPS, sample=SAMPLES
        ),
        mimag="output/summary/mimag_bins.tsv"
    output:
        per_workflow="output/drep_per_workflow/per_workflow.tsv"
    conda:
        DREP_ENV
    resources:
        mem_mb=2000
    log:
        "logs/drep_collect.log"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        python scripts/drep_stage.py collect \
          --per-workflow-root output/drep_per_workflow \
          --mimag-bins {input.mimag} --out {output.per_workflow}
        test -s {output.per_workflow}
        """


# L9: cluster membership + representative taxon propagation.
rule drep_taxonomy:
    input:
        mimag="output/summary/mimag_bins.tsv",
        cdb="output/drep_cross_workflow/data_tables/Cdb.csv",
        wdb="output/drep_cross_workflow/data_tables/Wdb.csv",
        manifest="output/drep_cross_workflow/manifest.tsv"
    output:
        final="output/summary/final_mags.tsv"
    conda:
        DREP_ENV
    resources:
        mem_mb=4000
    log:
        "logs/drep_taxonomy.log"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        python scripts/drep_stage.py taxonomy \
          --mimag-bins {input.mimag} --cdb {input.cdb} --wdb {input.wdb} \
          --manifest {input.manifest} --out {output.final}
        test -s {output.final}
        """
