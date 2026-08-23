###############################################################################
# AdapterRemoval collapse
###############################################################################
rule adapterremoval_collapse:
    input:
        r1p=W("fastq/{sample}_R1.paired.fastq.gz"),
        r2p=W("fastq/{sample}_R2.paired.fastq.gz")
    output:
        W("fastq/{sample}.collapsed.fastq.gz")
    conda:
        "../envs/adapterremoval.yaml"
    threads:
        config["adapterremoval"]["threads"]
    params:
        minq=config["adapterremoval"]["minquality"],
        minlen=config["adapterremoval"]["minlength"]
    shell:
        r"""
        AdapterRemoval \
            --file1 {input.r1p} \
            --file2 {input.r2p} \
            --gzip \
            --threads {threads} \
            --trimns \
            --trimqualities \
            --minquality {params.minq} \
            --minlength {params.minlen} \
            --collapse \
            --outputcollapsed {output}
        """

###############################################################################
# post stats per read (again per-sample/read files)
###############################################################################
rule stats_collapsed:
    input:
        W("fastq/{sample}.collapsed.fastq.gz")
    output:
        W("reports_read/{sample}.collapsed.length.tsv")
    conda:
        "../envs/seqkit.yaml"
    shell:
        r"""
        seqkit stats {input} > {output}
        """

###############################################################################
# bwa aln mapping
###############################################################################

rule bwa_aln_collapsed:
    input:
        ref=ancient(config["reference"]),
        fq=W("fastq/{sample}.collapsed.fastq.gz")
    output:
        temp(W("mapping/collapsed/{sample}.sai"))
    conda:
        "../envs/bwa.yaml"
    threads:
        config["bwa"]["threads"]
    params:
        aln_opts=config["bwa"].get("aln_opts", "")
    shell:
        r"""
        mkdir -p {WORK_DIR}/mapping/collapsed
        bwa aln -t {threads} {params.aln_opts} {input.ref} {input.fq} > {output}
        """


rule bwa_samse_collapsed:
    input:
        ref=ancient(config["reference"]),
        sai=W("mapping/collapsed/{sample}.sai"),
        fq=W("fastq/{sample}.collapsed.fastq.gz")
    output:
        temp(W("mapping/collapsed/{sample}.sam"))
    conda:
        "../envs/bwa.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/mapping/collapsed
        bwa samse {input.ref} {input.sai} {input.fq} > {output}
        """


rule flagstat_sam_collapsed:
    input:
        W("mapping/collapsed/{sample}.sam")
    output:
        W("reports_mapping/collapsed/{sample}.sam.flagstat.txt")
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/reports_mapping/collapsed
        samtools flagstat {input} > {output}
        """