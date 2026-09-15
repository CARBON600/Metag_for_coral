rule megahit_pre:
    input:
        r1=filtered_r1,
        r2=filtered_r2
    output:
        contigs="output/megahit_pre/{method}/{group}/{sample}/{sample}.contigs.fa"
    threads: config["threads"]["megahit"]
    conda:
        "megahit"
    resources:
        mem_mb=32000
    log:
        "logs/megahit_pre/{method}/{group}/{sample}.log"
    params:
        prefix=lambda wc: wc.sample,
        outdir=lambda wc: f"output/megahit_pre/{wc.method}/{wc.group}/{wc.sample}"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        megahit \
          -1 {input.r1} \
          -2 {input.r2} \
          -o {params.outdir} \
          --out-prefix {params.prefix} \
          -t {threads} \
          --min-contig-len 250
        test -s {output.contigs}
        """

rule build_megahit_bt2_index:
    input:
        contigs=megahit_contigs
    output:
        done="output/PCR_free/index/{method}/{group}/{sample}/build.done"
    threads: config["threads"]["bt2_build_contigs"]
    conda:
        "metag"
    resources:
        mem_mb=8000
    log:
        "logs/build_megahit_bt2_index/{method}/{group}/{sample}.log"
    params:
        prefix=lambda wc: f"output/PCR_free/index/{wc.method}/{wc.group}/{wc.sample}/{wc.sample}"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
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
        "metag"
    resources:
        mem_mb=16000
    log:
        "logs/remap_to_megahit_contigs/{method}/{group}/{sample}.log"
    params:
        prefix=lambda wc: f"output/PCR_free/index/{wc.method}/{wc.group}/{wc.sample}/{wc.sample}"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        mkdir -p $(dirname {output.bam})
        bowtie2 -p {threads} -x {params.prefix} -1 {input.r1} -2 {input.r2} \
          | samtools view -@ {threads} -bS - > {output.bam}
        """

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
        "metag"
    resources:
        mem_mb=16000
    log:
        "logs/pcr_dedup/{method}/{group}/{sample}.log"
    params:
        sorted_bam=lambda wc: f"output/PCR_done/{wc.method}/{wc.group}/{wc.sample}_sort.bam",
        markdup_bam=lambda wc: f"output/PCR_done/{wc.method}/{wc.group}/{wc.sample}_markdup.bam",
        name_bam=lambda wc: f"output/PCR_done/{wc.method}/{wc.group}/{wc.sample}_fq.bam"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
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

        test -s {output.dedup_bam}
        test -s {output.r1}
        test -s {output.r2}

        rm -f {params.sorted_bam} {params.sorted_bam}.bai
        rm -f {params.markdup_bam}
        rm -f {params.name_bam}
        """

rule spades_control_assemble:
    input:
        r1=filtered_r1,
        r2=filtered_r2
    output:
        contigs="output/assemble/control_assemble/{method}/{group}/{sample}/contigs.fasta"
    threads: config["threads"]["spades_control"]
    conda:
        "metag"
    resources:
        mem_mb=64000
    log:
        "logs/spades_control_assemble/{method}/{group}/{sample}.log"
    params:
        spades=config["software"]["spades"],
        outdir=lambda wc: f"output/assemble/control_assemble/{wc.method}/{wc.group}/{wc.sample}"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        {params.spades} --meta \
          -1 {input.r1} \
          -2 {input.r2} \
          -k 21,33,55,77,99,127 \
          --only-assembler \
          -t {threads} \
          -o {params.outdir}
        test -s {output.contigs}
        """

rule spades_pcr_assemble:
    input:
        r1="output/PCR_done/{method}/{group}/{sample}_1.fq.gz",
        r2="output/PCR_done/{method}/{group}/{sample}_2.fq.gz"
    output:
        contigs="output/assemble/PCR_assemble/{method}/{group}/{sample}/contigs.fasta"
    threads: config["threads"]["spades_pcr"]
    conda:
        "metag"
    resources:
        mem_mb=64000
    log:
        "logs/spades_pcr_assemble/{method}/{group}/{sample}.log"
    params:
        spades=config["software"]["spades"],
        outdir=lambda wc: f"output/assemble/PCR_assemble/{wc.method}/{wc.group}/{wc.sample}"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        {params.spades} --meta \
          -1 {input.r1} \
          -2 {input.r2} \
          -k 21,33,55,77,99,127 \
          --only-assembler \
          -t {threads} \
          -o {params.outdir}
        test -s {output.contigs}
        """

rule filter_contigs_r2000:
    input:
        contigs="output/assemble/{treat}/{method}/{group}/{sample}/contigs.fasta"
    output:
        filtered="output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta"
    conda:
        "metag"
    resources:
        mem_mb=4000
    log:
        "logs/filter_contigs_r2000/{treat}/{method}/{group}/{sample}.log"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        python -c "from Bio import SeqIO; seqs = [rec for rec in SeqIO.parse('{input.contigs}', 'fasta') if len(rec.seq) > 2000]; SeqIO.write(seqs, '{output.filtered}', 'fasta-2line')"
        test -s {output.filtered}
        """

rule build_final_contig_index:
    input:
        contigs="output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta"
    output:
        done="output/assemble/{treat}/{method}/{group}/{sample}/bt2_index.done"
    threads: config["threads"]["bt2_build_contigs"]
    conda:
        "metag"
    resources:
        mem_mb=8000
    log:
        "logs/build_final_contig_index/{treat}/{method}/{group}/{sample}.log"
    params:
        prefix=lambda wc: f"output/assemble/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}/{wc.sample}_bw2"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        bowtie2-build {input.contigs} {params.prefix}
        touch {output.done}
        """

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
        "metag"
    resources:
        mem_mb=16000
    log:
        "logs/remap_reads_to_final_contigs/{treat}/{method}/{group}/{sample}.log"
    params:
        prefix=lambda wc: f"output/assemble/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}/{wc.sample}_bw2",
        rawbam=lambda wc: f"output/assemble/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}/{wc.sample}_raw.bam"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        bowtie2 -p {threads} -x {params.prefix} -1 {input.r1} -2 {input.r2} \
          | samtools view -@ {threads} -bS - > {params.rawbam}

        samtools sort -@ {threads} -o {output.bam} {params.rawbam}
        samtools index {output.bam}
        rm -f {params.rawbam}

        test -s {output.bam}
        test -s {output.bai}
        """