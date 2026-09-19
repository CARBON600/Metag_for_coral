# Fail fast, before any long job starts, if a configured input path is absent.
rule check_inputs:
    output:
        marker="output/check_inputs.done"
    resources:
        mem_mb=1000
    log:
        "logs/check_inputs.log"
    params:
        paths=lambda wc: [
            config["data_dir"],
            config["tmpdir"],
            config["software"]["fastq_screen"],
            config["software"]["spades"],
            config["fastq_screen_conf"]["control"],
            config["fastq_screen_conf"]["exp"],
            bowtie2_index_probe(config["bowtie2_index"]["control"]),
            bowtie2_index_probe(config["bowtie2_index"]["exp"]),
            config["coverm_ref"]["control"],
            config["coverm_ref"]["exp"],
        ],
        data_dir=config["data_dir"],
        n_samples=lambda wc: len(SAMPLES),
        tmpdir=config["tmpdir"],
        conf_files=lambda wc: [
            config["fastq_screen_conf"]["control"],
            config["fastq_screen_conf"]["exp"],
        ]
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        missing=""
        for p in {params.paths}; do
            if [ -e "$p" ]; then echo "OK   $p"; else echo "MISS $p"; missing="$missing $p"; fi
        done
        if [ -n "$missing" ]; then
            echo "ERROR: missing required paths:$missing"
            exit 1
        fi
        if [ "{params.n_samples}" -eq 0 ]; then
            echo "ERROR: no paired-end samples matched {params.data_dir}/*_1.fq.gz"
            exit 1
        fi
        if [ ! -d "{params.tmpdir}" ]; then
            echo "ERROR: tmpdir is not a directory (MEGAHIT/SPAdes --tmp-dir need one): {params.tmpdir}"
            exit 1
        fi
        # Content-level check: every aligner and DATABASE index named by a
        # fastq_screen config must actually exist, not just the config file.
        #
        # Parse each line after stripping a trailing '#' comment and take the
        # LAST remaining token as the path.  The official DATABASE form is
        # "DATABASE <name> <index>" (3 fields), but a trailing comment would
        # otherwise occupy field 3 and be mistaken for the index path.  This
        # handles 2-field, 3-field and commented lines without relaxing the
        # "index must exist" requirement.
        #
        # A missing DATABASE index is always fatal: letting it through only
        # defers the failure by hours to the fastq_screen job.  Set
        # SKIP_CONF_CHECK=1 to downgrade a broken third-party ALIGNER line from
        # a hard stop to a warning.
        _skip_conf=""
        [ "${{SKIP_CONF_CHECK:-0}}" = "1" ] && _skip_conf=1
        for conf in {params.conf_files}; do
            if [ ! -e "$conf" ]; then echo "ERROR: fastq_screen conf missing: $conf"; exit 1; fi
            _conf_rows="$(awk '{{ sub(/#.*/, ""); if ($1 == "") next; print $1, $NF }}' "$conf" || true)"
            if [ -z "$_conf_rows" ]; then
                echo "ERROR: fastq_screen conf parsed to zero entries (empty conf or awk failure?): $conf"
                exit 1
            fi
            while read -r _kind _path; do
                case "$_kind" in
                    BOWTIE2|BOWTIE|BWA|BWAMEM|MINIMAP2)
                        if [ ! -x "$_path" ]; then
                            if [ -n "$_skip_conf" ]; then
                                echo "WARN: aligner not executable ($conf): $_path (SKIP_CONF_CHECK=1)"
                            else
                                echo "ERROR: aligner not executable ($conf): $_path"; exit 1
                            fi
                        fi ;;
                    DATABASE)
                        if [ ! -e "$_path" ] && [ ! -e "$_path.1.bt2" ] && [ ! -e "$_path.1.bt2l" ]; then
                            echo "ERROR: DATABASE index not found ($conf): $_path"; exit 1
                        fi ;;
                esac
            done <<< "$_conf_rows"
        done
        mkdir -p $(dirname {output.marker})
        touch {output.marker}
        """

# GTDB-Tk / CheckM reference data is only needed at the very end of the pipeline,
# so it is gated separately and must not block host removal or assembly.
rule check_qc_inputs:
    output:
        marker="output/check_qc_inputs.done"
    resources:
        mem_mb=1000
    log:
        "logs/check_qc_inputs.log"
    params:
        paths=lambda wc: [
            config["gtdbtk_data"],
            config["checkm_data"],
        ]
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        missing=""
        for p in {params.paths}; do
            if [ -e "$p" ]; then echo "OK   $p"; else echo "MISS $p"; missing="$missing $p"; fi
        done
        if [ -n "$missing" ]; then
            echo "ERROR: missing required QC paths:$missing"
            exit 1
        fi
        for p in {params.paths}; do
            if [ ! -d "$p" ] || [ -z "$(ls -A "$p" 2>/dev/null)" ]; then
                echo "ERROR: QC data directory missing or empty: $p"
                exit 1
            fi
        done
        mkdir -p $(dirname {output.marker})
        touch {output.marker}
        """

rule fastq_screen:
    input:
        checked="output/check_inputs.done",
        r1=raw_r1,
        r2=raw_r2
    output:
        r1="output/fastqs/{group}/{sample}/{sample}_1.tagged_filter.fastq.gz",
        r2="output/fastqs/{group}/{sample}/{sample}_2.tagged_filter.fastq.gz"
    threads: config["threads"]["fastq_screen"]
    conda:
        "metag"
    resources:
        mem_mb=8000
    log:
        "logs/fastq_screen/{group}/{sample}.log"
    params:
        conf=lambda wc: config["fastq_screen_conf"][wc.group],
        outdir=lambda wc: f"output/fastqs/{wc.group}/{wc.sample}",
        fastq_screen=config["software"]["fastq_screen"]
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        rm -f {params.outdir}/{wildcards.sample}_1.fq.gz_temp_subset.fastq \
              {params.outdir}/{wildcards.sample}_2.fq.gz_temp_subset.fastq
        {params.fastq_screen} \
          --conf {params.conf} \
          --outdir {params.outdir} \
          --nohits \
          --force \
          --aligner bowtie2 \
          --threads {threads} \
          {input.r1} {input.r2}

        test -s {output.r1}
        test -s {output.r2}
        """

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
        "metag"
    resources:
        mem_mb=4000
    log:
        "logs/fastp_clean/{group}/{sample}.log"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        fastp \
          -i {input.r1} -I {input.r2} \
          -o {output.r1} -O {output.r2} \
          -w {threads} \
          -j {output.json} \
          -h {output.html}
        test -s {output.r1}
        test -s {output.r2}
        """

rule bowtie2_map:
    input:
        checked="output/check_inputs.done",
        r1=raw_r1,
        r2=raw_r2
    output:
        bam="output/bt2/{group}/{sample}.bam"
    threads: config["threads"]["bowtie2_map"]
    conda:
        "metag"
    resources:
        mem_mb=16000
    log:
        "logs/bowtie2_map/{group}/{sample}.log"
    params:
        index=lambda wc: config["bowtie2_index"][wc.group]
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        mkdir -p $(dirname {output.bam})
        bowtie2 -p {threads} -x {params.index} -1 {input.r1} -2 {input.r2} \
          | samtools view -@ {threads} -bS - > {output.bam}
        test -s {output.bam}
        """

rule bowtie2_unmapped:
    input:
        bam="output/bt2/{group}/{sample}.bam"
    output:
        bam="output/bt2_unmapped/{group}/{sample}.bam"
    threads: config["threads"]["bowtie2_filter"]
    conda:
        "metag"
    resources:
        mem_mb=8000
    log:
        "logs/bowtie2_unmapped/{group}/{sample}.log"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        mkdir -p $(dirname {output.bam})
        samtools view -@ {threads} -b -f 4 {input.bam} > {output.bam}
        test -s {output.bam}
        """

rule bam_to_fastq_bt2:
    input:
        bam="output/bt2_unmapped/{group}/{sample}.bam"
    output:
        r1="output/fq4dep/bt2/{group}/{sample}_1.fq.gz",
        r2="output/fq4dep/bt2/{group}/{sample}_2.fq.gz"
    threads: config["threads"]["bam_to_fastq"]
    conda:
        "metag"
    resources:
        mem_mb=8000
    log:
        "logs/bam_to_fastq_bt2/{group}/{sample}.log"
    params:
        tmpbam=lambda wc: f"output/fq4dep/bt2/{wc.group}/{wc.sample}_name_sorted.bam"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        mkdir -p $(dirname {output.r1})
        samtools sort -n -@ {threads} -O BAM -o {params.tmpbam} {input.bam}
        samtools fastq -@ {threads} \
            -1 {output.r1} \
            -2 {output.r2} \
            -0 /dev/null -s /dev/null \
            {params.tmpbam}
        rm -f {params.tmpbam}
        test -s {output.r1}
        test -s {output.r2}
        """

rule coverm_map:
    input:
        checked="output/check_inputs.done",
        r1=raw_r1,
        r2=raw_r2
    output:
        bam="output/coverm/{group}/{sample}.bam"
    threads: config["threads"]["coverm"]
    conda:
        "metag"
    resources:
        mem_mb=16000
    log:
        "logs/coverm_map/{group}/{sample}.log"
    params:
        ref=lambda wc: config["coverm_ref"][wc.group],
        tmpdir=config["tmpdir"],
        outdir=lambda wc: f"output/coverm/{wc.group}/{wc.sample}_raw"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        mkdir -p $(dirname {output.bam})
        rm -rf {params.outdir}
        export TMPDIR={params.tmpdir}
        coverm make \
            -r {params.ref} \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {params.outdir} \
            -t {threads}
        bam=$(find {params.outdir} -maxdepth 1 -name '*.bam' | head -n 1 || true)
        if [ -z "$bam" ]; then
            echo "ERROR: 'coverm make' left no BAM in {params.outdir}; cannot continue (see the coverm output above)."
            exit 1
        fi
        mv "$bam" {output.bam}
        rm -rf {params.outdir}
        test -s {output.bam}
        """

rule coverm_filter:
    input:
        bam="output/coverm/{group}/{sample}.bam"
    output:
        bam="output/coverm_filtered/{group}/{sample}_filtered.bam"
    threads: config["threads"]["coverm"]
    conda:
        "metag"
    resources:
        mem_mb=16000
    log:
        "logs/coverm_filter/{group}/{sample}.log"
    params:
        tmpdir=config["tmpdir"]
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        mkdir -p $(dirname {output.bam})
        export TMPDIR={params.tmpdir}
        coverm filter \
          -b {input.bam} \
          -o {output.bam} \
          --inverse \
          --min-read-aligned-percent 0.75 \
          --min-read-percent-identity 0.95 \
          --threads {threads}
        test -s {output.bam}
        """

rule bam_to_fastq_coverm:
    input:
        bam="output/coverm_filtered/{group}/{sample}_filtered.bam"
    output:
        r1="output/fq4dep/coverm/{group}/{sample}_1.fq.gz",
        r2="output/fq4dep/coverm/{group}/{sample}_2.fq.gz"
    threads: config["threads"]["bam_to_fastq"]
    conda:
        "metag"
    resources:
        mem_mb=8000
    log:
        "logs/bam_to_fastq_coverm/{group}/{sample}.log"
    params:
        tmpbam=lambda wc: f"output/fq4dep/coverm/{wc.group}/{wc.sample}_name_sorted.bam"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        mkdir -p $(dirname {output.r1})
        samtools sort -n -@ {threads} -O BAM -o {params.tmpbam} {input.bam}
        samtools fastq -@ {threads} \
            -1 {output.r1} \
            -2 {output.r2} \
            -0 /dev/null -s /dev/null \
            {params.tmpbam}
        rm -f {params.tmpbam}
        test -s {output.r1}
        test -s {output.r2}
        """

# Kaiju is disabled: no environment is shipped and the target is absent from
# `rule all`. Its former config keys were removed; see the README for the
# rationale and the DB paths needed to reinstate it.
