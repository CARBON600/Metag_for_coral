# ---------------------------------------
# UniT: bin
rule unitem_bin:
    input:
        contigs="output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta",
        bam="output/assemble/{treat}/{method}/{group}/{sample}/{sample}_sorted.bam"
    output:
        done="output/binning/unitem/{treat}/{method}/{group}/{sample}/bin.done"
    threads: config["threads"]["unitem"]
    conda:
        "envs/binning.yaml"
    params:
        outdir=lambda wc: f"output/binning/unitem/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}"
    shell:
        r"""
        mkdir -p {params.outdir}
        unitem bin \
          --bam_files {input.bam} \
          --max40 --max107 \
          --mb2 --mb_verysensitive --mb_sensitive --mb_specific --mb_veryspecific --mb_superspecific \
          -c {threads} \
          {params.outdir}/bin
        touch {output.done}
        """

rule unitem_consensus:
    input:
        "output/binning/unitem/{treat}/{method}/{group}/{sample}/bin.done"
    output:
        done="output/binning/unitem/{treat}/{method}/{group}/{sample}/consensus.done"
    threads: config["threads"]["unitem"]
    conda:
        "envs/binning.yaml"
    params:
        outdir=lambda wc: f"output/binning/unitem/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}"
    shell:
        r"""
        unitem profile \
          -f {params.outdir}/bin/bin_dirs.tsv \
          -c {threads} \
          {params.outdir}/profile

        unitem consensus \
          -f {params.outdir}/bin/bin_dirs.tsv \
          {params.outdir}/profile \
          {params.outdir}/consensus

        touch {output.done}
        """

# ---------------------------------------
# COMEbin
rule comebin:
    input:
        contigs="output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta",
        bam="output/assemble/{treat}/{method}/{group}/{sample}/{sample}_sorted.bam"
    output:
        done="output/binning/comebin/{treat}/{method}/{group}/{sample}/comebin.done"
    threads: config["threads"]["comebin"]
    conda:
        "envs/binning.yaml"
    params:
        outdir=lambda wc: f"output/binning/comebin/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}",
        bamdir=lambda wc: f"output/binning/comebin/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}/bam"
    shell:
        r"""
        mkdir -p {params.outdir}
        mkdir -p {params.bamdir}

        ln -sf $(realpath {input.bam}) {params.bamdir}/$(basename {input.bam})
        ln -sf $(realpath {input.bam}.bai) {params.bamdir}/$(basename {input.bam}.bai)

        # NOTE: verify your run_comebin.sh interface on your installation
        run_comebin.sh \
          -a {input.contigs} \
          -o {params.outdir} \
          -p {params.bamdir} \
          -t {threads}

        touch {output.done}
        """

# ---------------------------------------
# MetaDecoder
rule metadecoder_coverage:
    input:
        contigs="output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta",
        bam="output/assemble/{treat}/{method}/{group}/{sample}/{sample}_sorted.bam"
    output:
        cov="output/binning/metadecoder/{treat}/{method}/{group}/{sample}/{sample}.COVERAGE"
    threads: config["threads"]["metadecoder_cov"]
    conda:
        "envs/binning.yaml"
    params:
        outdir=lambda wc: f"output/binning/metadecoder/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}"
    shell:
        r"""
        mkdir -p {params.outdir}
        metadecoder coverage \
          -b {input.bam} \
          -o {output.cov} \
          --threads {threads} \
          --bin_size 500000
        """

rule metadecoder_fix_coverage:
    input:
        "output/binning/metadecoder/{treat}/{method}/{group}/{sample}/{sample}.COVERAGE"
    output:
        "output/binning/metadecoder/{treat}/{method}/{group}/{sample}/{sample}_fixed.COVERAGE"
    conda:
        "envs/binning.yaml"
    shell:
        r"""
        grep -v '^#' {input} | grep -v '^sequence' > {output}
        """

rule metadecoder_seed:
    input:
        contigs="output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta",
        cov="output/binning/metadecoder/{treat}/{method}/{group}/{sample}/{sample}_fixed.COVERAGE"
    output:
        "output/binning/metadecoder/{treat}/{method}/{group}/{sample}/{sample}.SEED"
    threads: config["threads"]["metadecoder_seed"]
    conda:
        "envs/binning.yaml"
    shell:
        r"""
        metadecoder seed \
          --threads {threads} \
          -f {input.contigs} \
          -o {output}
        """

rule metadecoder_cluster:
    input:
        contigs="output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta",
        cov="output/binning/metadecoder/{treat}/{method}/{group}/{sample}/{sample}_fixed.COVERAGE",
        seed="output/binning/metadecoder/{treat}/{method}/{group}/{sample}/{sample}.SEED"
    output:
        done="output/binning/metadecoder/{treat}/{method}/{group}/{sample}/cluster.done"
    threads: config["threads"]["metadecoder_cluster"]
    conda:
        "envs/binning.yaml"
    params:
        outprefix=lambda wc: f"output/binning/metadecoder/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}/{wc.sample}.metadecoder"
    shell:
        r"""
        metadecoder cluster \
          -f {input.contigs} \
          -c {input.cov} \
          -s {input.seed} \
          -o {params.outprefix} \
          --min_sequence_length 2000 \
          --disable_gpu

        touch {output.done}
        """

# ---------------------------------------
# SemiBin2 single-sample mode
rule semibin2_single:
    input:
        contigs="output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta",
        bam="output/assemble/{treat}/{method}/{group}/{sample}/{sample}_sorted.bam"
    output:
        done="output/binning/semibin2_single/{treat}/{method}/{group}/{sample}/semibin.done"
    threads: config["threads"]["semibin"]
    conda:
        "envs/binning.yaml"
    params:
        outdir=lambda wc: f"output/binning/semibin2_single/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}"
    shell:
        r"""
        mkdir -p {params.outdir}
        SemiBin2 single_easy_bin \
          --self-supervised \
          --input-fasta {input.contigs} \
          --input-bam {input.bam} \
          --output {params.outdir} \
          -t {threads}

        touch {output.done}
        """

# ---------------------------------------
# RefineM for all binners
rule refinem_bins:
    input:
        bin_done=binner_done,
        bam="output/assemble/{treat}/{method}/{group}/{sample}/{sample}_sorted.bam"
    output:
        done="output/refinem/{binner}/{treat}/{method}/{group}/{sample}/refinem.done"
    threads: config["threads"]["refinem"]
    conda:
        "envs/binning.yaml"
    params:
        genomes=binner_genome_dir,
        ext=binner_extension,
        outdir=lambda wc: f"output/refinem/{wc.binner}/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}"
    shell:
        r"""
        mkdir -p {params.outdir}

        refinem scaffold_stats \
          -x {params.ext} \
          -c {threads} \
          {params.genomes} \
          {params.outdir} \
          {input.bam}

        refinem outliers \
          {params.outdir}/scaffold_stats.tsv \
          {params.outdir}

        refinem filter_bins \
          -x {params.ext} \
          {params.genomes} \
          {params.outdir}/outliers.tsv \
          {params.outdir}

        touch {output.done}
        """

# ---------------------------------------
# CheckM
rule checkm_lineage_wf:
    input:
        "output/refinem/{binner}/{treat}/{method}/{group}/{sample}/refinem.done"
    output:
        done="output/checkm/{binner}/{treat}/{method}/{group}/{sample}/checkm.done"
    threads: config["threads"]["checkm"]
    conda:
        "envs/binning.yaml"
    params:
        genomes=lambda wc: f"output/refinem/{wc.binner}/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}",
        ext=binner_extension,
        outdir=lambda wc: f"output/checkm/{wc.binner}/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}"
    shell:
        r"""
        mkdir -p {params.outdir}
        checkm lineage_wf \
          -t {threads} \
          -x {params.ext} \
          {params.genomes} \
          {params.outdir}

        touch {output.done}
        """

# ---------------------------------------
# GTDB-Tk
rule gtdbtk_classify:
    input:
        "output/refinem/{binner}/{treat}/{method}/{group}/{sample}/refinem.done"
    output:
        done="output/gtdb/{binner}/{treat}/{method}/{group}/{sample}/gtdb.done"
    threads: config["threads"]["gtdbtk"]
    conda:
        "envs/binning.yaml"
    params:
        genomes=lambda wc: f"output/refinem/{wc.binner}/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}",
        ext=binner_extension,
        outdir=lambda wc: f"output/gtdb/{wc.binner}/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}",
        prefix=lambda wc: f"{wc.sample}_{wc.binner}"
    shell:
        r"""
        mkdir -p {params.outdir}
        gtdbtk classify_wf \
          --genome_dir {params.genomes} \
          --skip_ani_screen \
          --out_dir {params.outdir} \
          --extension {params.ext} \
          --prefix {params.prefix} \
          --cpus {threads}

        touch {output.done}
        """