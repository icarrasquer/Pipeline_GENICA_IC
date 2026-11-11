
###############################################################################
# Map collapsed reads (single-end)
###############################################################################

rule bwa_aln_collapsed:
    input:
        ref=lambda wc: REF_DICT[wc.species],
        fq=W("{species}/fastq/{sample}.collapsed.fastq.gz")
    output:
        temp(W("{species}/mapping/collapsed/{sample}.sai"))
    conda:
        "../envs/bwa.yaml"
    threads:
        config["bwa"]["threads"]
    params:
        aln_opts=lambda wc: (
            config["bwa"].get("aln_opts", "")
            if wc.sample in ANCIENT_SAMPLES
            else ""
        )
    shell:
        r"""
        mkdir -p $(dirname {output})
        bwa aln -t {threads} {params.aln_opts} \
            {input.ref} {input.fq} > {output}
        """


rule bwa_samse_collapsed:
    input:
        ref=lambda wc: REF_DICT[wc.species],
        sai=W("{species}/mapping/collapsed/{sample}.sai"),
        fq=W("{species}/fastq/{sample}.collapsed.fastq.gz")
    output:
        temp(W("{species}/mapping/collapsed/{sample}.bam"))
    conda:
        "../envs/allpurpose.yaml"
    threads: 4
    shell:
        r"""
        set -euo pipefail

        bwa samse {input.ref} {input.sai} {input.fq} \
        | samtools view -@ {threads} -b -o {output} -
        """


###############################################################################
# Map uncollapsed R1 and R2 together (paired-end)
###############################################################################

rule bwa_aln_uncollapsed:
    input:
        ref=lambda wc: REF_DICT[wc.species],
        r1=W("{species}/fastq/{sample}.R1_uncoll.fastq.gz"),
        r2=W("{species}/fastq/{sample}.R2_uncoll.fastq.gz")
    output:
        sai1=temp(W("{species}/mapping/paired/{sample}_R1.sai")),
        sai2=temp(W("{species}/mapping/paired/{sample}_R2.sai"))
    conda:
        "../envs/bwa.yaml"
    threads:
        config["bwa"]["threads"]
    params:
        aln_opts=lambda wc: (
            config["bwa"].get("aln_opts", "")
            if wc.sample in ANCIENT_SAMPLES
            else ""
        )
    shell:
        r"""
        mkdir -p $(dirname {output.sai1})

        bwa aln -t {threads} {params.aln_opts} \
            {input.ref} {input.r1} > {output.sai1}

        bwa aln -t {threads} {params.aln_opts} \
            {input.ref} {input.r2} > {output.sai2}
        """


rule bwa_sampe_uncollapsed:
    input:
        ref=lambda wc: REF_DICT[wc.species],
        sai1=W("{species}/mapping/paired/{sample}_R1.sai"),
        sai2=W("{species}/mapping/paired/{sample}_R2.sai"),
        r1=W("{species}/fastq/{sample}.R1_uncoll.fastq.gz"),
        r2=W("{species}/fastq/{sample}.R2_uncoll.fastq.gz")
    output:
        temp(W("{species}/mapping/paired/{sample}.bam"))
    conda:
        "../envs/allpurpose.yaml"
    threads: 4
    shell:
        r"""
        set -euo pipefail

        bwa sampe {input.ref} {input.sai1} {input.sai2} \
            {input.r1} {input.r2} \
        | samtools view -@ {threads} -b -o {output} -
        """


###############################################################################
# Merge collapsed and uncollapsed alignments
###############################################################################

rule merge_mapping:
    input:
        collapsed=W("{species}/mapping/collapsed/{sample}.bam"),
        paired=W("{species}/mapping/paired/{sample}.bam")
    output:
        bam=W("{species}/mapping/collapsed_paired/{sample}.bam"),
        csi=W("{species}/mapping/collapsed_paired/{sample}.bam.csi")
    conda:
        "../envs/allpurpose.yaml"
    threads: 4
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.bam})

        samtools merge -u - \
            {input.collapsed} {input.paired} \
        | samtools sort -@ {threads} \
            -o {output.bam} -

        samtools index -c -@ {threads} {output.bam}
        """

###############################################################################
# Flagstat of merged BAM
###############################################################################

rule flagstat_merged:
    input:
        bam=W("{species}/mapping/collapsed_paired/{sample}.bam")
    output:
        W("{species}/reports_mapping/collapsed_paired/{sample}.flagstat.txt")
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        r"""
        mkdir -p {WORK_DIR}/{wildcards.species}/reports_mapping/collapsed_paired

        samtools flagstat -@ {threads} {input.bam} > {output}
        """