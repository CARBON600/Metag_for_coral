# MIMAG quality gate (Bowers 2017 Nat Biotechnol, PMID 28787424).
#
# Runs the CheckM/GTDB parser that tools/harvest_mags.py already contains, so
# there is one parser rather than two; the passing MAGs go to
# output/summary/mimag_bins.tsv, which is what the dRep layer below joins on.
rule mimag_gate:
    input:
        checkm=combo_paths("checkm", "checkm.done"),
        gtdb=combo_paths("gtdb", "gtdb.done"),
        # Force the CheckM QA tables to exist before harvest reads them (they are
        # the only source of `Strain heterogeneity`).
        checkm_qa=combo_paths("checkm", "results.tsv"),
        qc="output/check_qc_inputs.done"
    output:
        bins="output/summary/mimag_bins.tsv",
        summary="output/summary/mimag_summary.tsv",
        excluded="output/summary/mimag_excluded.tsv"
    conda:
        "drep"
    resources:
        mem_mb=2000
    log:
        "logs/mimag_gate.log"
    params:
        min_comp=config["mimag"]["min_completeness"],
        max_con=config["mimag"]["max_contamination"],
        # Rendered as space-separated tokens and passed straight through.
        binners=" ".join(BINNERS),
        treats=" ".join(ASSEMBLY_TREATS),
        methods=" ".join(FILTER_METHODS),
        groups=" ".join(GROUPS),
        samples=" ".join(SAMPLES)
    shell:
        r"""
        mkdir -p $(dirname {log}) output/summary
        exec > {log} 2>&1
        python tools/harvest_mags.py \
          --output-dir output --outdir output/summary --emit-mimag \
          --min-completeness {params.min_comp} --max-contamination {params.max_con} \
          --binners {params.binners} --treats {params.treats} \
          --methods {params.methods} --groups {params.groups} --samples {params.samples}
        # The script always writes headers, so an empty table is still a
        # non-empty file and these checks only catch a failed run.
        test -s {output.bins}
        test -s {output.summary}
        test -s {output.excluded}
        """
