###############################################################################
# bwa aln mapping
###############################################################################

rule bwa_aln_paired:
    input:
        ref=config["reference"],
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
    shell:
        r"""
        mkdir -p {WORK_DIR}/mapping/paired

        bwa aln -t {threads} {params.aln_opts} {input.ref} {input.r1} > {output.sai1}
        bwa aln -t {threads} {params.aln_opts} {input.ref} {input.r2} > {output.sai2}
        """


rule bwa_sampe_paired:
    input:
        ref=config["reference"],
        sai1=W("mapping/paired/{sample}_R1.sai"),
        sai2=W("mapping/paired/{sample}_R2.sai"),
        r1=W("fastq/{sample}_R1.paired.fastq.gz"),
        r2=W("fastq/{sample}_R2.paired.fastq.gz")
    output:
        temp(W("mapping/paired/{sample}.sam"))
    conda:
        "../envs/bwa.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/mapping/paired
        bwa sampe {input.ref} {input.sai1} {input.sai2} {input.r1} {input.r2} > {output}
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
