import os
from glob import glob
import pandas as pd

#Config file for the run
configfile: "config/config_droc_oldlibprep.yaml"

#Loading samples codes and new name assignation
META = pd.read_csv("config/Droc_batch1_IDs.txt", sep="\t")

OLD_TO_NEW = dict(zip(META.lib_code, META.sample_name))
NEW_TO_OLD = dict(zip(META.sample_name, META.lib_code))

SAMPLES = list(META.sample_name)

#Define R1 and R2
READS = ["R1", "R2"]
wildcard_constraints:
    read="R1|R2"

#Define input and output directories
INPUT_DIR = config["input_dir"]
WORK_DIR = config["work_dir"]
FASTQC_OUT = os.path.join(WORK_DIR, config["fastqc"]["outdir"])

#Find the R1 and R2 to emrg
def raw_fastqs(sample, read):
    old_sample = NEW_TO_OLD[sample]

    if read == "R1":
        pattern = os.path.join(INPUT_DIR, f"*{old_sample}*_R1_*.fastq.gz")
    elif read == "R2":
        pattern = os.path.join(INPUT_DIR, f"*{old_sample}*_R2_*.fastq.gz")
    else:
        raise ValueError(f"Unexpected read: {read}")

    matches = sorted(glob(pattern))

    if not matches:
        raise ValueError(
            f"No raw FASTQs found for sample={sample}, read={read} in {INPUT_DIR}. "
            f"Tried: {pattern}"
        )

    return matches

#Helper for paths in WORK_DIR
def W(*parts):
    return os.path.join(WORK_DIR, *parts)

###############################################################################
# Set necessary ouputs
###############################################################################

rule all:
    input:
        # cleaned gz fastqs
        expand(W("fastq/{sample}_{read}.clean.fastq.gz"), sample=SAMPLES, read=READS),
        # fastqc outputs for cleaned reads
        expand(os.path.join(FASTQC_OUT, "{sample}_{read}.clean_fastqc.html"), sample=SAMPLES, read=READS),
        # collapsed output
        expand(W("fastq/{sample}.collapsed.fastq.gz"), sample=SAMPLES),
        # mapping outputs
        expand(W("mapping/{sample}.sam"), sample=SAMPLES),
        expand(W("mapping/{sample}_F3084.bam"), sample=SAMPLES),
        expand(W("mapping/{sample}_F3084.bam.bai"), sample=SAMPLES),

        # stats
        expand(W("reports_per_read/{sample}_{read}.pre.stats.tsv"), sample=SAMPLES, read=READS),
        expand(W("reports_per_read/{sample}_{read}.pre.length.tsv"), sample=SAMPLES, read=READS),
        expand(W("reports_per_read/{sample}_{read}.post.stats.tsv"), sample=SAMPLES, read=READS),
        expand(W("reports_per_read/{sample}_{read}.post.length.tsv"), sample=SAMPLES, read=READS),
        expand(W("reports_per_read/{sample}_R1.paired.length.tsv"), sample=SAMPLES),
        expand(W("reports_per_read/{sample}.collapsed.length.tsv"), sample=SAMPLES),
        expand(W("reports_mapping/{sample}.sam.flagstat.txt"), sample=SAMPLES)
        

###############################################################################
# Step 1: merge raw reads into a temp fastq (uncompressed)
###############################################################################
rule merge_raw:
    input:
        lambda wc: raw_fastqs(wc.sample, wc.read)
    output:
        temp(W("tmp/{sample}_{read}.merged.fastq"))
    shell:
        r"""
        mkdir -p {WORK_DIR}/tmp
        zcat {input} > {output}
        """

###############################################################################
# Step 2: initial stats (write per-sample/read stats file)
###############################################################################
rule stats_pre:
    input:
        W("tmp/{sample}_{read}.merged.fastq")
    output:
        W("reports_per_read/{sample}_{read}.pre.stats.tsv"),
        W("reports_per_read/{sample}_{read}.pre.length.tsv")
    params:
        a3 = config["adapters"]["a_3prime"],
        g5 = config["adapters"]["g_5prime_loop"]
    conda:
        "envs/seqkit.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/reports_per_read
        nb_reads=$(grep -c "." {input} | awk '{{print $1/4}}')
        nb_Rd1_SP=$(grep -c "{params.g5}" {input} || echo 0)
        nb_Rd2_SP=$(grep -c "{params.a3}" {input} || echo 0)
        echo -e "{wildcards.sample}_{wildcards.read}\tpre_merged\t$nb_reads\t$nb_Rd1_SP\t$nb_Rd2_SP" > {output[0]}
        seqkit stats {input} > {output[1]}
        """

###############################################################################
# Step 3: adapter trimming + looping 5' trimming until motif hits <= threshold
###############################################################################
rule cutadapt_and_loop:
    input:
        W("tmp/{sample}_{read}.merged.fastq")
    output:
        temp(W("tmp/{sample}_{read}.cutadapt1.fastq"))
    conda:
        "envs/cutadapt.yaml"
    threads:
        config["cutadapt"]["threads"]
    params:
        a3 = config["adapters"]["a_3prime"],
        g5 = config["adapters"]["g_5prime_loop"],
        minlen = config["cutadapt"]["minlen"],
        e = config["cutadapt"]["error_rate"],
        maxhits = config["loop_stop"]["max_hits"]
    shell:
        r"""
        mkdir -p {WORK_DIR}/tmp

        cutadapt -j {threads} -e {params.e} -a {params.a3} {input} -o {output} -m {params.minlen}

        num=$(grep -c {params.g5} {output} || echo 0)
        while [ "$num" -gt {params.maxhits} ]; do
            cutadapt -j {threads} -e {params.e} -g {params.g5} {output} -o {WORK_DIR}/tmp/{wildcards.sample}_{wildcards.read}.loop.fastq -m {params.minlen}
            mv {WORK_DIR}/tmp/{wildcards.sample}_{wildcards.read}.loop.fastq {output}
            num=$(grep -c {params.g5} {output} || echo 0)
        done
        """

###############################################################################
# Step 4: polyG/polyC removal via your perl script
###############################################################################
rule drop_poly:
    input:
        W("tmp/{sample}_{read}.cutadapt1.fastq")
    output:
        temp(W("tmp/{sample}_{read}.poly.fastq"))
    shell:
        r"""
        mkdir -p {WORK_DIR}/tmp
        DropBpFastq_polyC.pl {input} {output}
        """

###############################################################################
# Step 5: quality + length filter -> final cleaned fastq.gz
###############################################################################
rule quality_filter:
    input:
        W("tmp/{sample}_{read}.poly.fastq")
    output:
        W("fastq/{sample}_{read}.clean.fastq.gz")
    conda:
        "envs/cutadapt.yaml"
    threads:
        config["cutadapt"]["threads"]
    params:
        q = config["cutadapt"]["qual"],
        minlen = config["cutadapt"]["minlen"]
    shell:
        r"""
        mkdir -p {WORK_DIR}/fastq {WORK_DIR}/tmp
        cutadapt -j {threads} -q {params.q} {input} -o {WORK_DIR}/tmp/{wildcards.sample}_{wildcards.read}.clean.fastq -m {params.minlen}
        gzip -c {WORK_DIR}/tmp/{wildcards.sample}_{wildcards.read}.clean.fastq > {output}
        rm -f {WORK_DIR}/tmp/{wildcards.sample}_{wildcards.read}.clean.fastq
        """

###############################################################################
# Step 6: post stats per read (again per-sample/read files)
###############################################################################
rule stats_post:
    input:
        W("fastq/{sample}_{read}.clean.fastq.gz")
    output:
        W("reports_per_read/{sample}_{read}.post.stats.tsv"),
        W("reports_per_read/{sample}_{read}.post.length.tsv")
    params:
        a3 = config["adapters"]["a_3prime"],
        g5 = config["adapters"]["g_5prime_loop"]
    conda:
        "envs/seqkit.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/reports_per_read
        nb_reads=$(zgrep -c "." {input} | awk '{{print $1/4}}')
        nb_Rd1_SP=$(zgrep -c {params.g5} {input} || echo 0)
        nb_Rd2_SP=$(zgrep -c {params.a3} {input} || echo 0)
        echo -e "{wildcards.sample}_{wildcards.read}\tpost_clean\t$nb_reads\t$nb_Rd1_SP\t$nb_Rd2_SP" > {output[0]}
        seqkit stats {input} > {output[1]}
        """

###############################################################################
# Step 7: FastQC on cleaned reads
###############################################################################
rule fastqc:
    input:
        W("fastq/{sample}_{read}.clean.fastq.gz")
    output:
        html = os.path.join(FASTQC_OUT, "{sample}_{read}.clean_fastqc.html"),
        zip  = os.path.join(FASTQC_OUT, "{sample}_{read}.clean_fastqc.zip")
    conda:
        "envs/fastqc.yaml"
    shell:
        r"""
        mkdir -p {FASTQC_OUT}
        fastqc {input} -o {FASTQC_OUT}
        """

###############################################################################
# Step 8: pair/repair reads (BBMap repair.sh)
###############################################################################
rule repair_reads:
    input:
        r1=W("fastq/{sample}_R1.clean.fastq.gz"),
        r2=W("fastq/{sample}_R2.clean.fastq.gz")
    output:
        r1p=W("fastq/{sample}_R1.paired.fastq.gz"),
        r2p=W("fastq/{sample}_R2.paired.fastq.gz")
    threads: 1
    resources:
        repair_slots=1,
        mem_mb=45000
    params:
        xmx = config["bbmap"]["xmx"],  # e.g. "60g"
    conda:
        "envs/bbmap.yaml"
    shell:
        r"""
        set -euo pipefail
        module load BBMap/38.96-GCC-10.3.0

        repair.sh -Xmx{params.xmx} \
          in1={input.r1} in2={input.r2} \
          out={output.r1p} out2={output.r2p} \
          overwrite=t
        """

###############################################################################
# Step 9: post stats per read (again per-sample/read files)
###############################################################################
rule stats_paired:
    input:
        W("fastq/{sample}_R1.paired.fastq.gz")
    output:
        W("reports_per_read/{sample}_R1.paired.length.tsv")
    conda:
        "envs/seqkit.yaml"
    shell:
        r"""
        seqkit stats {input} > {output}
        """


###############################################################################
# Step 10: AdapterRemoval collapse
###############################################################################
rule adapterremoval_collapse:
    input:
        r1p=W("fastq/{sample}_R1.paired.fastq.gz"),
        r2p=W("fastq/{sample}_R2.paired.fastq.gz")
    output:
        W("fastq/{sample}.collapsed.fastq.gz")
    conda:
        "envs/adapterremoval.yaml"
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
# Step 10: post stats per read (again per-sample/read files)
###############################################################################
rule stats_collapsed:
    input:
        W("fastq/{sample}.collapsed.fastq.gz")
    output:
        W("reports_per_read/{sample}.collapsed.length.tsv")
    conda:
        "envs/seqkit.yaml"
    shell:
        r"""
        seqkit stats {input} > {output}
        """

###############################################################################
# Step 11: bwa aln mapping
###############################################################################

###############################################################################
# Mapping: bwa aln (SE) on collapsed reads -> SAM -> flagstat -> sorted BAM
###############################################################################

rule bwa_aln:
    input:
        ref=config["reference"],
        fq=W("fastq/{sample}.collapsed.fastq.gz")
    output:
        temp(W("mapping/{sample}.sai"))
    conda:
        "envs/bwa.yaml"
    threads:
        config["bwa"]["threads"]
    params:
        aln_opts=config["bwa"].get("aln_opts", "")
    shell:
        r"""
        mkdir -p {WORK_DIR}/mapping
        bwa aln -t {threads} {params.aln_opts} {input.ref} {input.fq} > {output}
        """

rule bwa_samse:
    input:
        ref=config["reference"],
        sai=W("mapping/{sample}.sai"),
        fq=W("fastq/{sample}.collapsed.fastq.gz")
    output:
        temp(W("mapping/{sample}.sam"))
    conda:
        "envs/bwa.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/mapping
        bwa samse {input.ref} {input.sai} {input.fq} > {output}
        """

rule flagstat_sam:
    input:
        W("mapping/{sample}.sam")
    output:
        W("reports_mapping/{sample}.sam.flagstat.txt")
    conda:
        "envs/samtools.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/reports_mapping
        samtools flagstat {input} > {output}
        """

rule sam_to_sorted_bam:
    input:
        W("mapping/{sample}.sam")
    output:
        bam=W("mapping/{sample}_F3084.bam"),
        bai=W("mapping/{sample}_F3084.bam.bai")
    conda:
        "envs/samtools.yaml"
    threads: 4
    shell:
        r"""
        mkdir -p {WORK_DIR}/mapping

        # Convert SAM -> sorted BAM
        samtools view -@ {threads} -bh -F 3084 {input} \
          | samtools sort -@ {threads} -o {output.bam} -

        samtools index -@ {threads} {output.bam}
        """

 

# Next steps -> mapq20, mapdamage
