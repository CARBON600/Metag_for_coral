# ---------------------------------------
# FASTQ_SCREEN
# ---------------------------------------
rule fastq_screen:
    input:
        r1=raw_r1,
        r2=raw_r2
    output:
        r1="output/fastqs/{group}/{sample}/{sample}_1.tagged_filter.fastq.gz",
        r2="output/fastqs/{group}/{sample}/{sample}_2.tagged_filter.fastq.gz"
    threads: config["threads"]["fastq_screen"]
    conda:
        "envs/metag.yaml"
    params:
        conf=lambda wc: config["fastq_screen_conf"][wc.group],
        outdir=lambda wc: f"output/fastqs/{wc.group}/{wc.sample}"
    shell:
        r"""
        mkdir -p {params.outdir}
        # NOTE:
        fastq_screen \
          --conf {params.conf} \
          --outdir {params.outdir} \
          --nohits \
          --aligner bowtie2 \
          {input.r1} {input.r2}

        test -s {output.r1}
        test -s {output.r2}
        """

# ---------------------------------------
# FASTP clean after fastq_screen
# ---------------------------------------
rule fastp_clean:
    input:
        r1="output/fastqs/{group}/{sample}/{sample}_1.tagged_filter.fastq.gz",
        r2="output/fastqs/{group}/{sample}/{sample}_2.tagged_filter.fastq.gz"
    output:
        r1="output/fastqs/{group}/{sample}/{sample}_1_cleaned.fastq.gz",
        r2="output/fastqs/{group}/{sample}/{sample}_2_cleaned.fastq.gz",
        json="output/fastqs/{group}/{sample}/{sample}_fastp.json",
        html="output/fastqs/{group}/{sample}/{sample}_fastp.html"
    threads: config["threads"]["fastp"]
    conda:
        "envs/metag.yaml"
    shell:
        r"""
        fastp \
          -i {input.r1} -I {input.r2} \
          -o {output.r1} -O {output.r2} \
          -w {threads} \
          -j {output.json} \
          -h {output.html}
        """

# ---------------------------------------
# Bowtie2 mapping against control / exp reference
# ---------------------------------------
rule bowtie2_map:
    input:
        r1=raw_r1,
        r2=raw_r2
    output:
        bam="output/bt2/{group}/{sample}.bam"
    threads: config["threads"]["bowtie2_map"]
    conda:
        "envs/metag.yaml"
    params:
        index=lambda wc: config["bowtie2_index"][wc.group]
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        bowtie2 -p {threads} -x {params.index} -1 {input.r1} -2 {input.r2} \
          | samtools view -@ {threads} -bS - > {output.bam}
        """

rule bowtie2_unmapped:
    input:
        bam="output/bt2/{group}/{sample}.bam"
    output:
        bam="output/bt2_unmapped/{group}/{sample}.bam"
    threads: config["threads"]["bowtie2_filter"]
    conda:
        "envs/metag.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        samtools view -@ {threads} -b -f 4 {input.bam} > {output.bam}
        """

rule bam_to_fastq_bt2:
    input:
        bam="output/bt2_unmapped/{group}/{sample}.bam"
    output:
        r1="output/fq4dep/bt2/{group}/{sample}_1.fq.gz",
        r2="output/fq4dep/bt2/{group}/{sample}_2.fq.gz"
    threads: config["threads"]["bam_to_fastq"]
    conda:
        "envs/metag.yaml"
    params:
        tmpbam=lambda wc: f"output/fq4dep/bt2/{wc.group}/{wc.sample}_name_sorted.bam"
    shell:
        r"""
        mkdir -p $(dirname {output.r1})
        samtools sort -n -@ {threads} -O BAM -o {params.tmpbam} {input.bam}
        samtools fastq -@ {threads} \
            -1 {output.r1} \
            -2 {output.r2} \
            -0 /dev/null -s /dev/null \
            {params.tmpbam}
        rm -f {params.tmpbam}
        """

# ---------------------------------------
# CoverM mapping / filtering
# ---------------------------------------
rule coverm_map:
    input:
        r1=raw_r1,
        r2=raw_r2
    output:
        bam="output/coverm/{group}/{sample}.bam"
    threads: config["threads"]["coverm"]
    conda:
        "envs/metag.yaml"
    params:
        ref=lambda wc: config["coverm_ref"][wc.group],
        tmpdir=config["tmpdir"]
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        export TMPDIR={params.tmpdir}
        coverm make \
            -r {params.ref} \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {output.bam} \
            -t {threads}
        """

rule coverm_filter:
    input:
        bam="output/coverm/{group}/{sample}.bam"
    output:
        bam="output/coverm_filtered/{group}/{sample}_filtered.bam"
    threads: config["threads"]["coverm"]
    conda:
        "envs/metag.yaml"
    params:
        tmpdir=config["tmpdir"]
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        export TMPDIR={params.tmpdir}
        coverm filter \
          -b {input.bam} \
          -o {output.bam} \
          --inverse \
          --min-read-aligned-percent 0.75 \
          --min-read-percent-identity 0.95 \
          --threads {threads}
        """

rule bam_to_fastq_coverm:
    input:
        bam="output/coverm_filtered/{group}/{sample}_filtered.bam"
    output:
        r1="output/fq4dep/coverm/{group}/{sample}_1.fq.gz",
        r2="output/fq4dep/coverm/{group}/{sample}_2.fq.gz"
    threads: config["threads"]["bam_to_fastq"]
    conda:
        "envs/metag.yaml"
    params:
        tmpbam=lambda wc: f"output/fq4dep/coverm/{wc.group}/{wc.sample}_name_sorted.bam"
    shell:
        r"""
        mkdir -p $(dirname {output.r1})
        samtools sort -n -@ {threads} -O BAM -o {params.tmpbam} {input.bam}
        samtools fastq -@ {threads} \
            -1 {output.r1} \
            -2 {output.r2} \
            -0 /dev/null -s /dev/null \
            {params.tmpbam}
        rm -f {params.tmpbam}
        """

# ---------------------------------------
# Raw read taxonomic profiling with Kaiju
# ---------------------------------------
rule kaiju_raw:
    input:
        r1=raw_r1,
        r2=raw_r2
    output:
        out="output/kaiju/raw/{sample}_kaiju.out",
        summary="output/kaiju/raw/{sample}_kaiju_summary.tsv"
    threads: config["threads"]["kaiju"]
    conda:
        "envs/metag.yaml"
    params:
        nodes=config["kaiju"]["nodes"],
        names=config["kaiju"]["names"],
        db=config["kaiju"]["db"]
    shell:
        r"""
        mkdir -p $(dirname {output.out})
        kaiju \
          -t {params.nodes} \
          -f {params.db} \
          -i {input.r1} \
          -j {input.r2} \
          -o {output.out} \
          -z {threads}

        kaiju2table \
          -t {params.nodes} \
          -n {params.names} \
          -o {output.summary} \
          {output.out} \
          -l superkingdom,phylum,class,order,family,genus,species \
          -r phylum
        """
