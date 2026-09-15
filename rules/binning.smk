rule unitem_bin:
    input:
        contigs="output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta",
        bam="output/assemble/{treat}/{method}/{group}/{sample}/{sample}_sorted.bam"
    output:
        done="output/binning/unitem/{treat}/{method}/{group}/{sample}/bin.done"
    threads: config["threads"]["unitem"]
    conda:
        "binning"
    resources:
        mem_mb=32000
    log:
        "logs/unitem_bin/{treat}/{method}/{group}/{sample}.log"
    params:
        outdir=lambda wc: f"output/binning/unitem/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        unitem bin \
          --bam_files {input.bam} \
          --max40 --max107 \
          --mb2 --mb_verysensitive --mb_sensitive --mb_specific --mb_veryspecific --mb_superspecific \
          -c {threads} \
          {input.contigs} \
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
        "binning"
    resources:
        mem_mb=32000
    log:
        "logs/unitem_consensus/{treat}/{method}/{group}/{sample}.log"
    params:
        outdir=lambda wc: f"output/binning/unitem/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
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

rule comebin:
    input:
        contigs="output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta",
        bam="output/assemble/{treat}/{method}/{group}/{sample}/{sample}_sorted.bam"
    output:
        done="output/binning/comebin/{treat}/{method}/{group}/{sample}/comebin.done"
    threads: config["threads"]["comebin"]
    conda:
        "comebin_env"
    resources:
        mem_mb=64000
    log:
        "logs/comebin/{treat}/{method}/{group}/{sample}.log"
    params:
        outdir=lambda wc: f"output/binning/comebin/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}",
        bamdir=lambda wc: f"output/binning/comebin/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}/bam"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        mkdir -p {params.bamdir}

        ln -sf $(realpath {input.bam}) {params.bamdir}/$(basename {input.bam})
        ln -sf $(realpath {input.bam}.bai) {params.bamdir}/$(basename {input.bam}.bai)

        run_comebin.sh \
          -a {input.contigs} \
          -o {params.outdir} \
          -p {params.bamdir} \
          -t {threads}

        touch {output.done}
        """

rule metadecoder_coverage:
    input:
        contigs="output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta",
        bam="output/assemble/{treat}/{method}/{group}/{sample}/{sample}_sorted.bam"
    output:
        cov="output/binning/metadecoder/{treat}/{method}/{group}/{sample}/{sample}.COVERAGE"
    threads: config["threads"]["metadecoder_cov"]
    conda:
        "metadecoder"
    resources:
        mem_mb=16000
    log:
        "logs/metadecoder_coverage/{treat}/{method}/{group}/{sample}.log"
    params:
        outdir=lambda wc: f"output/binning/metadecoder/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        metadecoder coverage \
          -b {input.bam} \
          -o {output.cov} \
          --threads {threads} \
          --bin_size 500000
        """

rule metadecoder_seed:
    input:
        contigs="output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta",
        cov="output/binning/metadecoder/{treat}/{method}/{group}/{sample}/{sample}.COVERAGE"
    output:
        "output/binning/metadecoder/{treat}/{method}/{group}/{sample}/{sample}.SEED"
    threads: config["threads"]["metadecoder_seed"]
    conda:
        "metadecoder"
    resources:
        mem_mb=32000
    log:
        "logs/metadecoder_seed/{treat}/{method}/{group}/{sample}.log"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        metadecoder seed \
          --threads {threads} \
          -f {input.contigs} \
          -o {output}
        """

rule metadecoder_cluster:
    input:
        contigs="output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta",
        cov="output/binning/metadecoder/{treat}/{method}/{group}/{sample}/{sample}.COVERAGE",
        seed="output/binning/metadecoder/{treat}/{method}/{group}/{sample}/{sample}.SEED"
    output:
        done="output/binning/metadecoder/{treat}/{method}/{group}/{sample}/cluster.done"
    threads: config["threads"]["metadecoder_cluster"]
    conda:
        "metadecoder"
    resources:
        mem_mb=64000
    log:
        "logs/metadecoder_cluster/{treat}/{method}/{group}/{sample}.log"
    params:
        outprefix=lambda wc: f"output/binning/metadecoder/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}/{wc.sample}.metadecoder"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        metadecoder cluster \
          -f {input.contigs} \
          -c {input.cov} \
          -s {input.seed} \
          -o {params.outprefix} \
          --min_sequence_length 2000 \
          --disable_gpu

        touch {output.done}
        """

rule semibin2_single:
    input:
        contigs="output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta",
        bam="output/assemble/{treat}/{method}/{group}/{sample}/{sample}_sorted.bam"
    output:
        done="output/binning/semibin2_single/{treat}/{method}/{group}/{sample}/semibin.done"
    threads: config["threads"]["semibin"]
    conda:
        "SemiBin"
    resources:
        mem_mb=32000
    log:
        "logs/semibin2_single/{treat}/{method}/{group}/{sample}.log"
    params:
        outdir=lambda wc: f"output/binning/semibin2_single/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        SemiBin2 single_easy_bin \
          --self-supervised \
          --input-fasta {input.contigs} \
          --input-bam {input.bam} \
          --output {params.outdir} \
          -t {threads}

        touch {output.done}
        """

# Normalise each binner's output (decompress unitem/semibin bins, keep only bin
# FASTA files) so RefineM/CheckM do not see COVERAGE/SEED/txt files.
rule prepare_bins:
    input:
        bin_done=binner_done
    output:
        done="output/bins_prepared/{binner}/{treat}/{method}/{group}/{sample}/prepare.done"
    conda:
        "binning"
    resources:
        mem_mb=2000
    log:
        "logs/prepare_bins/{binner}/{treat}/{method}/{group}/{sample}.log"
    params:
        raw=binner_raw_genome_dir,
        outdir=binner_prepared_dir,
        ext=binner_extension,
        gz=binner_bins_gzipped
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        if [ "{params.gz}" = "True" ]; then
          for f in {params.raw}/*.{params.ext}.gz; do
            [ -e "$f" ] || continue
            gzip -cd "$f" > {params.outdir}/$(basename "$f" .gz)
          done
        else
          for f in {params.raw}/*.{params.ext}; do
            [ -e "$f" ] || continue
            cp "$f" {params.outdir}/
          done
        fi
        test -n "$(ls -A {params.outdir})"
        touch {output.done}
        """

rule refinem_bins:
    input:
        prepared=binner_prepared_done,
        contigs="output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta",
        bam="output/assemble/{treat}/{method}/{group}/{sample}/{sample}_sorted.bam"
    output:
        done="output/refinem/{binner}/{treat}/{method}/{group}/{sample}/refinem.done"
    threads: config["threads"]["refinem"]
    conda:
        "binning"
    resources:
        mem_mb=16000
    log:
        "logs/refinem_bins/{binner}/{treat}/{method}/{group}/{sample}.log"
    params:
        genomes=binner_prepared_dir,
        ext=binner_extension,
        outdir=lambda wc: f"output/refinem/{wc.binner}/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        mkdir -p {params.outdir}

        refinem scaffold_stats \
          -x {params.ext} \
          -c {threads} \
          {input.contigs} \
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

rule checkm_lineage_wf:
    input:
        "output/refinem/{binner}/{treat}/{method}/{group}/{sample}/refinem.done"
    output:
        done="output/checkm/{binner}/{treat}/{method}/{group}/{sample}/checkm.done"
    threads: config["threads"]["checkm"]
    conda:
        "binning"
    params:
        genomes=lambda wc: f"output/refinem/{wc.binner}/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}",
        ext=binner_extension,
        outdir=lambda wc: f"output/checkm/{wc.binner}/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}",
        checkm_data=config["checkm_data"]
    resources:
        mem_mb=32000
    log:
        "logs/checkm_lineage_wf/{binner}/{treat}/{method}/{group}/{sample}.log"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        export CHECKM_DATA_PATH={params.checkm_data}
        mkdir -p {params.outdir}
        checkm lineage_wf \
          -t {threads} \
          -x {params.ext} \
          {params.genomes} \
          {params.outdir}

        touch {output.done}
        """

rule gtdbtk_classify:
    input:
        "output/refinem/{binner}/{treat}/{method}/{group}/{sample}/refinem.done"
    output:
        done="output/gtdb/{binner}/{treat}/{method}/{group}/{sample}/gtdb.done"
    threads: config["threads"]["gtdbtk"]
    conda:
        "gtdbtk-2.3.2"
    params:
        genomes=lambda wc: f"output/refinem/{wc.binner}/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}",
        ext=binner_extension,
        outdir=lambda wc: f"output/gtdb/{wc.binner}/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}",
        prefix=lambda wc: f"{wc.sample}_{wc.binner}",
        gtdbtk_data=config["gtdbtk_data"]
    resources:
        mem_mb=64000
    log:
        "logs/gtdbtk_classify/{binner}/{treat}/{method}/{group}/{sample}.log"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        export GTDBTK_DATA_PATH={params.gtdbtk_data}
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