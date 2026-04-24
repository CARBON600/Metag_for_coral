# ---------------------------------------
# MEGAHIT pre-assembly
# method: bt2 / coverm / fastqs
# group : control / exp
rule megahit_pre:
    input:
        r1=filtered_r1,
        r2=filtered_r2
    output:
        directory("output/megahit_pre/{method}/{group}/{sample}")
    threads: config["threads"]["megahit"]
    conda:
        "envs/metag.yaml"
    params:
        prefix=lambda wc: wc.sample
    shell:
        r"""
        mkdir -p {output}
        megahit \
          -1 {input.r1} \
          -2 {input.r2} \
          -o {output} \
          --out-prefix {params.prefix} \
          -t {threads} \
          --min-contig-len 250
        """

rule build_megahit_bt2_index:
    input:
        contigs=megahit_contigs
    output:
        done="output/PCR_free/index/{method}/{group}/{sample}/build.done"
    threads: config["threads"]["bt2_build_contigs"]
    conda:
        "envs/metag.yaml"
    params:
        prefix=lambda wc: f"output/PCR_free/index/{wc.method}/{wc.group}/{wc.sample}/{wc.sample}"
    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        bowtie2-build {input.contigs} {params.prefix}
        touch {output.done}
        """

rule remap_to_megahit_contigs:
    input:
        idx_done="output/PCR_free/index/{method}/{group}/{sample}/build.done",
        r1=filtered_r1,
        r2=filtered_r2
    output:
        bam="output/PCR_free/{method}/{group}/{sample}.bam"
    threads: config["threads"]["remap_contigs"]
    conda:
        "envs/metag.yaml"
    params:
        prefix=lambda wc: f"output/PCR_free/index/{wc.method}/{wc.group}/{wc.sample}/{wc.sample}"
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        bowtie2 -p {threads} -x {params.prefix} -1 {input.r1} -2 {input.r2} \
          | samtools view -@ {threads} -bS - > {output.bam}
        """

# ---------------------------------------
# PCR dedup and export deduplicated reads
# ---------------------------------------
rule pcr_dedup:
    input:
        bam="output/PCR_free/{method}/{group}/{sample}.bam"
    output:
        dedup_bam="output/PCR_done/{method}/{group}/{sample}_dedup.bam",
        r1="output/PCR_done/{method}/{group}/{sample}_1.fq.gz",
        r2="output/PCR_done/{method}/{group}/{sample}_2.fq.gz",
        stats="output/PCR_done/{method}/{group}/{sample}_stats_file.txt"
    threads: config["threads"]["pcr_dedup"]
    conda:
        "envs/metag.yaml"
    params:
        sorted_bam=lambda wc: f"output/PCR_done/{wc.method}/{wc.group}/{wc.sample}_sort.bam",
        markdup_bam=lambda wc: f"output/PCR_done/{wc.method}/{wc.group}/{wc.sample}_markdup.bam",
        name_bam=lambda wc: f"output/PCR_done/{wc.method}/{wc.group}/{wc.sample}_fq.bam"
    shell:
        r"""
        mkdir -p $(dirname {output.dedup_bam})

        samtools sort -@ {threads} -o {params.sorted_bam} {input.bam}
        samtools index {params.sorted_bam}

        samtools view -F 2304 -b {params.sorted_bam} \
          | samtools collate -@ {threads} -O -u - \
          | samtools fixmate -@ {threads} -m -u - - \
          | samtools sort -@ {threads} -u - \
          | samtools markdup -@ {threads} -f {output.stats} - {params.markdup_bam}

        samtools view -@ {threads} -b -F 1024 {params.markdup_bam} > {output.dedup_bam}

        samtools sort -n -@ {threads} -O BAM -o {params.name_bam} {output.dedup_bam}
        samtools fastq -@ {threads} \
            -1 {output.r1} \
            -2 {output.r2} \
            -0 /dev/null -s /dev/null \
            {params.name_bam}
        """

# ---------------------------------------
# Final assembly: control_assemble
# uses filtered reads directly
rule spades_control_assemble:
    input:
        r1=lambda wc: f"output/fastqs/{wc.group}/{wc.sample}/{wc.sample}_1_cleaned.fastq.gz" if wc.method == "fastqs"
                      else f"output/fq4dep/{wc.method}/{wc.group}/{wc.sample}_1.fq.gz",
        r2=lambda wc: f"output/fastqs/{wc.group}/{wc.sample}/{wc.sample}_2_cleaned.fastq.gz" if wc.method == "fastqs"
                      else f"output/fq4dep/{wc.method}/{wc.group}/{wc.sample}_2.fq.gz"
    output:
        directory("output/assemble/control_assemble/{method}/{group}/{sample}")
    threads: config["threads"]["spades_control"]
    conda:
        "envs/metag.yaml"
    shell:
        r"""
        mkdir -p {output}
        spades.py --meta \
          -1 {input.r1} \
          -2 {input.r2} \
          -k 21,33,55,77,99,127 \
          --only-assembler \
          -t {threads} \
          -o {output}
        """

# ---------------------------------------
# Final assembly: PCR_assemble
# uses PCR-deduplicated reads
rule spades_pcr_assemble:
    input:
        r1="output/PCR_done/{method}/{group}/{sample}_1.fq.gz",
        r2="output/PCR_done/{method}/{group}/{sample}_2.fq.gz"
    output:
        directory("output/assemble/PCR_assemble/{method}/{group}/{sample}")
    threads: config["threads"]["spades_pcr"]
    conda:
        "envs/metag.yaml"
    shell:
        r"""
        mkdir -p {output}
        spades.py --meta \
          -1 {input.r1} \
          -2 {input.r2} \
          -k 21,33,55,77,99,127 \
          --only-assembler \
          -t {threads} \
          -o {output}
        """

# ---------------------------------------
# Keep contigs > 2000 bp
# ---------------------------------------
rule filter_contigs_r2000:
    input:
        "output/assemble/{treat}/{method}/{group}/{sample}/contigs.fasta"
    output:
        "output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta"
    conda:
        "envs/metag.yaml"
    run:
        from Bio import SeqIO
        seqs = [rec for rec in SeqIO.parse(input[0], "fasta") if len(rec.seq) > 2000]
        with open(output[0], "w") as fh:
            SeqIO.write(seqs, fh, "fasta-2line")

# ---------------------------------------
# Build Bowtie2 index for final contigs
# ---------------------------------------
rule build_final_contig_index:
    input:
        contigs="output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta"
    output:
        done="output/assemble/{treat}/{method}/{group}/{sample}/bt2_index.done"
    threads: config["threads"]["bt2_build_contigs"]
    conda:
        "envs/metag.yaml"
    params:
        prefix=lambda wc: f"output/assemble/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}/{wc.sample}_bw2"
    shell:
        r"""
        bowtie2-build {input.contigs} {params.prefix}
        touch {output.done}
        """

# ---------------------------------------
# Remap reads back to final contigs for binning
# ---------------------------------------
rule remap_reads_to_final_contigs:
    input:
        idx_done="output/assemble/{treat}/{method}/{group}/{sample}/bt2_index.done",
        r1=assembly_reads_r1,
        r2=assembly_reads_r2
    output:
        bam="output/assemble/{treat}/{method}/{group}/{sample}/{sample}_sorted.bam",
        bai="output/assemble/{treat}/{method}/{group}/{sample}/{sample}_sorted.bam.bai"
    threads: config["threads"]["remap_final_contigs"]
    conda:
        "envs/metag.yaml"
    params:
        prefix=lambda wc: f"output/assemble/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}/{wc.sample}_bw2",
        rawbam=lambda wc: f"output/assemble/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}/{wc.sample}_raw.bam"
    shell:
        r"""
        bowtie2 -p {threads} -x {params.prefix} -1 {input.r1} -2 {input.r2} \
          | samtools view -@ {threads} -bS - > {params.rawbam}

        samtools sort -@ {threads} -o {output.bam} {params.rawbam}
        samtools index {output.bam}
        rm -f {params.rawbam}
        """