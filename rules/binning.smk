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
        rm -rf {params.outdir}/bin
        unitem bin \
          --bam_files {input.bam} \
          --max40 --max107 \
          --mb2 --mb_verysensitive --mb_sensitive --mb_specific --mb_veryspecific --mb_superspecific \
          -c {threads} \
          {input.contigs} \
          {params.outdir}/bin
        n=$(ls -A {params.outdir}/bin 2>/dev/null | wc -l || true)
        echo "unitem bin: $n files"
        [ "$n" -eq 0 ] && echo "WARNING: unitem produced no bins (allowed)"
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
        rm -rf {params.outdir}/profile {params.outdir}/consensus
        unitem profile \
          -f {params.outdir}/bin/bin_dirs.tsv \
          -c {threads} \
          {params.outdir}/profile

        unitem consensus \
          -f {params.outdir}/bin/bin_dirs.tsv \
          {params.outdir}/profile \
          {params.outdir}/consensus

        n=$(ls -A {params.outdir}/consensus/bins 2>/dev/null | wc -l || true)
        echo "unitem consensus: $n bins"
        [ "$n" -eq 0 ] && echo "WARNING: 0 consensus bins (allowed)"
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
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        mkdir -p {params.bamdir}

        ln -sf $(realpath {input.bam}) {params.bamdir}/$(basename {input.bam})
        ln -sf $(realpath {input.bam}.bai) {params.bamdir}/$(basename {input.bam}.bai)

        run_comebin.sh \
          -a {input.contigs} \
          -o {params.outdir} \
          -p {params.bamdir} \
          -t {threads}

        n=$(ls -A {params.outdir}/comebin_res/comebin_res_bins 2>/dev/null | wc -l || true)
        echo "comebin: $n bins"
        [ "$n" -eq 0 ] && echo "WARNING: 0 COMEBin bins (allowed)"
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
        # MetaDecoder caches an index next to the BAM (<bam>.index) and reuses it
        # if present; drop it so a regenerated BAM never reuses a stale index.
        rm -f {input.bam}.index
        metadecoder coverage \
          -b {input.bam} \
          -o {output.cov} \
          --threads {threads} \
          --bin_size 500000
        # MetaDecoder writes an index beside the BAM while it runs; remove it so
        # no undeclared file is left under output/assemble/.
        rm -f {input.bam}.index
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
        rm -f {params.outprefix}.*.fasta
        # MetaDecoder 1.2.1 writes <basename(fasta)>.<min_len>.metadecoder.{{kmers,dpgmm}}
        # into the CWD and loads them on a later run if present.  Every combination
        # uses the same input basename (contigs_r2000bp.fasta) and
        # --min_sequence_length 2000, so a shared CWD silently reuses another
        # combination's cache.  Run in a private CWD and address inputs/output by
        # absolute path so the cache is always rebuilt.
        _root=$PWD
        _cwd=$_root/$(dirname {params.outprefix})/_metadecoder_cwd
        rm -rf "$_cwd"
        mkdir -p "$_cwd"
        _out=$_root/{params.outprefix}
        ( cd "$_cwd" && metadecoder cluster \
            -f "$_root/{input.contigs}" \
            -c "$_root/{input.cov}" \
            -s "$_root/{input.seed}" \
            -o "$_out" \
            --min_sequence_length 2000 \
            --disable_gpu )
        rm -rf "$_cwd"

        n=$(ls -1 {params.outprefix}.*.fasta 2>/dev/null | wc -l || true)
        echo "metadecoder cluster: $n bins"
        [ "$n" -eq 0 ] && echo "WARNING: 0 MetaDecoder bins (allowed)"
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
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        # Fixed seed so bins are reproducible across runs (SemiBin2 default is
        # system-seeded, i.e. "results may vary between runs").
        SemiBin2 single_easy_bin \
          --self-supervised \
          --input-fasta {input.contigs} \
          --input-bam {input.bam} \
          --output {params.outdir} \
          --random-seed 1 \
          -t {threads}

        n=$(ls -A {params.outdir}/output_bins 2>/dev/null | wc -l || true)
        echo "semibin2: $n bins"
        [ "$n" -eq 0 ] && echo "WARNING: 0 SemiBin2 bins (allowed)"
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
        n_all=0
        for f in {params.raw}/*; do
          [ -e "$f" ] || continue
          case "$f" in
            *.COVERAGE|*.SEED|*.done) ;;
            *) n_all=$((n_all+1)) ;;
          esac
        done
        n=$(ls -A {params.outdir} 2>/dev/null | wc -l || true)
        echo "prepare_bins: candidates=$n_all prepared=$n"
        if [ "$n_all" -gt 0 ] && [ "$n" -eq 0 ]; then
          echo "ERROR: {params.raw} has $n_all file(s) but none matched extension '{params.ext}'"
          exit 1
        fi
        [ "$n" -eq 0 ] && echo "WARNING: no bins prepared (allowed)"
        touch {output.done}
        """

rule refinem_bins:
    input:
        checked="output/check_qc_inputs.done",
        prepared=binner_prepared_done,
        contigs="output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta",
        bam="output/assemble/{treat}/{method}/{group}/{sample}/{sample}_sorted.bam"
    output:
        done="output/refinem/{binner}/{treat}/{method}/{group}/{sample}/refinem.done"
    # RefineM filters scaffolds within a bin using the original assembly contigs
    # and the remapped BAM. The metaWRAP arm's bins are per-bin SPAdes
    # reassemblies whose contigs are absent from contigs_r2000bp.fasta and
    # {sample}_sorted.bam, so RefineM does not apply; rules/metawrap.smk's
    # metawrap_ingest provides refinem.done for that binner instead.
    wildcard_constraints:
        binner="unitem|comebin|metadecoder|semibin2_single",
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
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        n=$(ls -A {params.genomes} 2>/dev/null | wc -l || true)
        if [ "$n" -eq 0 ]; then
          echo "WARNING: no prepared bins; writing empty RefineM tables"
          printf 'genome\tvalues\n' > {params.outdir}/scaffold_stats.tsv
          printf 'genome\tvalues\n' > {params.outdir}/outliers.tsv
          touch {output.done}
          exit 0
        fi

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
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        n=$(ls {params.genomes}/*.{params.ext} 2>/dev/null | wc -l || true)
        if [ "$n" -eq 0 ]; then
          echo "WARNING: no bins after RefineM; writing empty CheckM tables"
          mkdir -p {params.outdir}/storage
          : > {params.outdir}/storage/bin_stats_ext.tsv
          : > {params.outdir}/results.tsv
          touch {output.done}
          exit 0
        fi
        checkm lineage_wf \
          -t {threads} \
          -x {params.ext} \
          {params.genomes} \
          {params.outdir}

        touch {output.done}
        """


# `checkm lineage_wf` prints its out_format-2 QA table (the only place that
# carries `Strain heterogeneity`, and a title-case `Marker lineage`) to stdout,
# which lands in the rule log but NOT in storage/bin_stats_ext.tsv. This rule
# re-derives that table as a declared output so tools/harvest_mags.py can read
# it. It re-parses existing HMMER output only: no gene calling, no HMMER rerun.
rule checkm_qa:
    input:
        "output/checkm/{binner}/{treat}/{method}/{group}/{sample}/checkm.done"
    output:
        "output/checkm/{binner}/{treat}/{method}/{group}/{sample}/results.tsv"
    threads: config["threads"]["checkm_qa"]
    conda:
        "binning"
    params:
        outdir=lambda wc: f"output/checkm/{wc.binner}/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}",
        checkm_data=config["checkm_data"]
    resources:
        mem_mb=config["resources"]["checkm_qa_mem_mb"]
    log:
        "logs/checkm_qa/{binner}/{treat}/{method}/{group}/{sample}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        export CHECKM_DATA_PATH={params.checkm_data}
        d={params.outdir}
        # lineage_wf writes lineage.ms at the analysis root; accept it under
        # storage/ too in case a CheckM build places it there.
        ms="$(find "$d" -maxdepth 2 -name 'lineage.ms' | head -n 1 || true)"
        if [ -z "$ms" ]; then
          echo "WARNING: no lineage.ms under $d (empty CheckM branch); writing empty results.tsv"
          : > {output}
          exit 0
        fi
        # `checkm qa` calls getBinIdsFromOutDir(), which lists <analyze_dir>/bins;
        # guard against a partially cleaned directory so this fails soft, not hard.
        if [ ! -d "$d/bins" ]; then
          echo "WARNING: no bins/ under $d; writing empty results.tsv"
          : > {output}
          exit 0
        fi
        checkm qa "$ms" "$d" -o 2 --tab_table -f {output} -t {threads}
        test -s {output}
        """


rule gtdbtk_classify:
    input:
        "output/refinem/{binner}/{treat}/{method}/{group}/{sample}/refinem.done"
    output:
        done="output/gtdb/{binner}/{treat}/{method}/{group}/{sample}/gtdb.done"
    threads: config["threads"]["gtdbtk"]
    # pplacer's memory is charged per forked thread: GTDB-Tk passes
    # -j <--pplacer_cpus or --cpus>, and the host reports
    # PARENT_MEMORY * (N_CHILDREN + 1) even though copy-on-write shares the pages.
    # For the pinned ref data (2.3.2, r207/r214) the parent is ~55-61 GB
    # (PPLACER_MIN_RAM_BAC_SPLIT); current ref data is ~140 GB.  A job killed for
    # memory dies with PplacerException while the class-level *.out stops at
    # "Preparing the edges for baseball..." with no error, so it is intermittent
    # and node-dependent; one retry is worth more than a wider declaration.
    retries: 1
    conda:
        "gtdbtk-2.3.2"
    params:
        genomes=lambda wc: f"output/refinem/{wc.binner}/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}",
        ext=binner_extension,
        outdir=lambda wc: f"output/gtdb/{wc.binner}/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}",
        prefix=lambda wc: f"{wc.sample}_{wc.binner}",
        gtdbtk_data=config["gtdbtk_data"]
    resources:
        # pplacer is accounted as PARENT x (1 + pplacer_cpus); the pinned ref
        # data (2.3.2, r207/r214) is ~55-61 GB per fork, current data ~140 GB.
        # Declare ~2x so the scheduler does not co-schedule GTDB-Tk jobs.
        mem_mb=config["resources"]["gtdbtk_mem_mb"]
    log:
        "logs/gtdbtk_classify/{binner}/{treat}/{method}/{group}/{sample}.log"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        export GTDBTK_DATA_PATH={params.gtdbtk_data}
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        n=$(ls {params.genomes}/*.{params.ext} 2>/dev/null | wc -l || true)
        if [ "$n" -eq 0 ]; then
          echo "WARNING: no bins after RefineM; writing empty GTDB-Tk summaries"
          : > {params.outdir}/{params.prefix}.bac120.summary.tsv
          : > {params.outdir}/{params.prefix}.ar122.summary.tsv
          touch {output.done}
          exit 0
        fi
        # --pplacer_cpus 1 pins pplacer to a single fork (the official answer to the
        # per-thread memory blow-up); --scratch_dir trades RAM for disk.
        gtdbtk classify_wf \
          --genome_dir {params.genomes} \
          --skip_ani_screen \
          --out_dir {params.outdir} \
          --extension {params.ext} \
          --prefix {params.prefix} \
          --cpus {threads} \
          --pplacer_cpus 1 \
          --scratch_dir {params.outdir}/scratch

        touch {output.done}
        """