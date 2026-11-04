###############################################################################
# bwa aln mapping
###############################################################################

rule bwa_aln_paired:
    input:
        ref=ancient(config["reference"]),
        r1=W("fastq/{sample}_R1.paired.fastq.gz"),
        r2=W("fastq/{sample}_R2.paired.fastq.gz")
    output:
        sai1=temp(W("mapping/paired/{sample}_R1.sai")),
        sai2=temp(W("mapping/paired/{sample}_R2.sai"))
    conda:
        "../envs/bwa.yaml"
    threads:
        config["bwa"]["threads"]
    params:
        aln_opts=lambda wildcards: config["bwa"].get("aln_opts", "") \
                if wildcards.sample in ANCIENT_SAMPLES else ""
    wildcard_constraints:
        sample=NORMAL_SAMPLE_REGEX
    shell:
        r"""
        mkdir -p {WORK_DIR}/mapping/paired

        bwa aln -t {threads} {params.aln_opts} {input.ref} {input.r1} > {output.sai1}
        bwa aln -t {threads} {params.aln_opts} {input.ref} {input.r2} > {output.sai2}
        """


rule bwa_sampe_paired:
    input:
        ref=ancient(config["reference"]),
        sai1=W("mapping/paired/{sample}_R1.sai"),
        sai2=W("mapping/paired/{sample}_R2.sai"),
        r1=W("fastq/{sample}_R1.paired.fastq.gz"),
        r2=W("fastq/{sample}_R2.paired.fastq.gz")
    output:
        temp(W("mapping/paired/{sample}.sam"))
    conda:
        "../envs/bwa.yaml"
    wildcard_constraints:
        sample=NORMAL_SAMPLE_REGEX
    shell:
        r"""
        mkdir -p {WORK_DIR}/mapping/paired
        bwa sampe {input.ref} {input.sai1} {input.sai2} {input.r1} {input.r2} > {output}
        """


###############################################################################
# Split large paired-end FASTQ files into two synchronized parts
###############################################################################

rule split_large_paired_fastq:
    input:
        r1=W("fastq/{sample}_R1.paired.fastq.gz"),
        r2=W("fastq/{sample}_R2.paired.fastq.gz")
    output:
        r1_part1=temp(
            W("fastq_split/{sample}/{sample}_R1.part1.fastq.gz")
        ),
        r1_part2=temp(
            W("fastq_split/{sample}/{sample}_R1.part2.fastq.gz")
        ),
        r2_part1=temp(
            W("fastq_split/{sample}/{sample}_R2.part1.fastq.gz")
        ),
        r2_part2=temp(
            W("fastq_split/{sample}/{sample}_R2.part2.fastq.gz")
        )
    wildcard_constraints:
        sample=LARGE_SAMPLE_REGEX
    conda:
        "../envs/python.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/fastq_split/{wildcards.sample}

        python scripts/split_paired_fastq.py \
            --r1 {input.r1} \
            --r2 {input.r2} \
            --r1-part1 {output.r1_part1} \
            --r1-part2 {output.r1_part2} \
            --r2-part1 {output.r2_part1} \
            --r2-part2 {output.r2_part2}
        """

###############################################################################
# bwa aln mapping: split samples
###############################################################################

rule bwa_aln_large_paired:
    input:
        ref=ancient(config["reference"]),
        r1=W(
            "fastq_split/{sample}/{sample}_R1.part{chunk}.fastq.gz"
        ),
        r2=W(
            "fastq_split/{sample}/{sample}_R2.part{chunk}.fastq.gz"
        )
    output:
        sai1=temp(
            W("mapping/paired_split/{sample}/part{chunk}_R1.sai")
        ),
        sai2=temp(
            W("mapping/paired_split/{sample}/part{chunk}_R2.sai")
        )
    wildcard_constraints:
        sample=LARGE_SAMPLE_REGEX,
        chunk="1|2"
    conda:
        "../envs/bwa.yaml"
    threads:
        config["bwa"]["threads"]
    params:
        aln_opts=lambda wildcards: (
            config["bwa"].get("aln_opts", "")
            if wildcards.sample in ANCIENT_SAMPLES
            else ""
        )
    shell:
        r"""
        mkdir -p \
            {WORK_DIR}/mapping/paired_split/{wildcards.sample}

        bwa aln \
            -t {threads} \
            {params.aln_opts} \
            {input.ref} \
            {input.r1} \
            > {output.sai1}

        bwa aln \
            -t {threads} \
            {params.aln_opts} \
            {input.ref} \
            {input.r2} \
            > {output.sai2}
        """

###############################################################################
# bwa sampe mapping: split samples
###############################################################################

rule bwa_sampe_large_paired:
    input:
        ref=ancient(config["reference"]),
        sai1=W(
            "mapping/paired_split/{sample}/part{chunk}_R1.sai"
        ),
        sai2=W(
            "mapping/paired_split/{sample}/part{chunk}_R2.sai"
        ),
        r1=W(
            "fastq_split/{sample}/{sample}_R1.part{chunk}.fastq.gz"
        ),
        r2=W(
            "fastq_split/{sample}/{sample}_R2.part{chunk}.fastq.gz"
        )
    output:
        temp(
            W("mapping/paired_split/{sample}/part{chunk}.sorted.bam")
        )
    wildcard_constraints:
        sample=LARGE_SAMPLE_REGEX,
        chunk="1|2"
    conda:
        "../envs/bwa.yaml"
    threads:
        config["bwa"]["threads"]
    resources:
        mem_mb=config["bwa"].get("split_mem_mb", 16000)
    params:
        sort_mem=config["bwa"].get("sort_mem", "2G")
    shell:
        r"""
        mkdir -p \
            {WORK_DIR}/mapping/paired_split/{wildcards.sample}

        bwa sampe \
            {input.ref} \
            {input.sai1} \
            {input.sai2} \
            {input.r1} \
            {input.r2} \
        | samtools view \
            -@ {threads} \
            -u \
            - \
        | samtools sort \
            -@ {threads} \
            -m {params.sort_mem} \
            -o {output} \
            -
        """

###############################################################################
# Merge split mappings and produce the standard SAM output
###############################################################################

rule merge_large_paired_mapping:
    input:
        part1=W(
            "mapping/paired_split/{sample}/part1.sorted.bam"
        ),
        part2=W(
            "mapping/paired_split/{sample}/part2.sorted.bam"
        )
    output:
        temp(W("mapping/paired/{sample}.sam"))
    wildcard_constraints:
        sample=LARGE_SAMPLE_REGEX
    conda:
        "../envs/samtools.yaml"
    threads:
        config["bwa"]["threads"]
    params:
        merged_bam=lambda wildcards: W(
            f"mapping/paired_split/{wildcards.sample}/merged.sorted.bam"
        )
    shell:
        r"""
        mkdir -p {WORK_DIR}/mapping/paired
        mkdir -p \
            {WORK_DIR}/mapping/paired_split/{wildcards.sample}

        samtools merge \
            -@ {threads} \
            -f \
            {params.merged_bam} \
            {input.part1} \
            {input.part2}

        samtools view \
            -@ {threads} \
            -h \
            -o {output} \
            {params.merged_bam}

        rm -f {params.merged_bam}
        """


rule flagstat_sam_paired:
    input:
        W("mapping/paired/{sample}.sam")
    output:
        W("reports_mapping/paired/{sample}.sam.flagstat.txt")
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/reports_mapping/paired
        samtools flagstat {input} > {output}
        """