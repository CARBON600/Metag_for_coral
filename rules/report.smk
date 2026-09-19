import os

# Standalone reporting workflow.  NOT part of the main DAG: run it explicitly with
#   snakemake -s rules/report.smk harvest --cores 1
# It only reads output/ and writes harvest/, so it can never perturb results.
_base = workflow.basedir.rstrip("/\\")
if os.path.basename(_base) == "rules":
    _base = os.path.dirname(_base)
_SCRIPT = os.path.join(_base, "scripts", "harvest_mags.py")


rule harvest:
    output:
        combo="harvest/combo_status.tsv",
        quality="harvest/mag_quality.tsv"
    log:
        "logs/harvest.log"
    params:
        script=_SCRIPT,
        output_dir="output",
        harvest_dir="harvest"
    shell:
        r"""
        mkdir -p $(dirname {log}) {params.harvest_dir}
        exec > {log} 2>&1
        python3 {params.script} --output-dir {params.output_dir} --outdir {params.harvest_dir}
        test -s {output.combo}
        """
