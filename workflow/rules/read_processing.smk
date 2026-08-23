###############################################################################
# merge raw reads into a temp fastq (uncompressed)
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
# initial stats (write per-sample/read stats file)
###############################################################################
rule stats_pre:
    input:
        W("tmp/{sample}_{read}.merged.fastq")
    output:
        W("reports_read/{sample}_{read}.pre.stats.tsv"),
        W("reports_read/{sample}_{read}.pre.length.tsv")
    params:
        a3 = config["adapters"]["a_3prime"],
        g5 = config["adapters"]["g_5prime_loop"]
    conda:
        "../envs/seqkit.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/reports_read
        nb_reads=$(grep -c "." {input} | awk '{{print $1/4}}')
        nb_Rd1_SP=$(grep -c "{params.g5}" {input} || echo 0)
        nb_Rd2_SP=$(grep -c "{params.a3}" {input} || echo 0)
        echo -e "{wildcards.sample}_{wildcards.read}\tpre_merged\t$nb_reads\t$nb_Rd1_SP\t$nb_Rd2_SP" > {output[0]}
        seqkit stats {input} > {output[1]}
        """

###############################################################################
# adapter trimming + looping 5' trimming until motif hits <= threshold
###############################################################################
rule cutadapt_and_loop:
    input:
        W("tmp/{sample}_{read}.merged.fastq")
    output:
        temp(W("tmp/{sample}_{read}.cutadapt1.fastq"))
    conda:
        "../envs/cutadapt.yaml"
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
# polyG/polyC removal via your perl script
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
# quality + length filter -> final cleaned fastq.gz
###############################################################################
rule quality_filter:
    input:
        W("tmp/{sample}_{read}.poly.fastq")
    output:
        W("fastq/{sample}_{read}.clean.fastq.gz")
    conda:
        "../envs/cutadapt.yaml"
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
# post stats per read (again per-sample/read files)
###############################################################################
rule stats_post:
    input:
        W("fastq/{sample}_{read}.clean.fastq.gz")
    output:
        W("reports_read/{sample}_{read}.post.stats.tsv"),
        W("reports_read/{sample}_{read}.post.length.tsv")
    params:
        a3 = config["adapters"]["a_3prime"],
        g5 = config["adapters"]["g_5prime_loop"]
    conda:
        "../envs/seqkit.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/reports_read
        nb_reads=$(zgrep -c "." {input} | awk '{{print $1/4}}')
        nb_Rd1_SP=$(zgrep -c {params.g5} {input} || echo 0)
        nb_Rd2_SP=$(zgrep -c {params.a3} {input} || echo 0)
        echo -e "{wildcards.sample}_{wildcards.read}\tpost_clean\t$nb_reads\t$nb_Rd1_SP\t$nb_Rd2_SP" > {output[0]}
        seqkit stats {input} > {output[1]}
        """

###############################################################################
# FastQC on cleaned reads
###############################################################################
rule fastqc:
    input:
        W("fastq/{sample}_{read}.clean.fastq.gz")
    output:
        html = os.path.join(FASTQC_OUT, "{sample}_{read}.clean_fastqc.html"),
        zip  = os.path.join(FASTQC_OUT, "{sample}_{read}.clean_fastqc.zip")
    conda:
        "../envs/fastqc.yaml"
    shell:
        r"""
        mkdir -p {FASTQC_OUT}
        fastqc {input} -o {FASTQC_OUT}
        """

###############################################################################
# repair reads
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
        mem_mb=70000
    params:
        xmx = config["bbmap"]["xmx"],  
    conda:
        "../envs/bbmap.yaml"
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
# post stats per read (again per-sample/read files)
###############################################################################
rule stats_paired:
    input:
        W("fastq/{sample}_R1.paired.fastq.gz")
    output:
        W("reports_read/{sample}_R1.paired.length.tsv")
    conda:
        "../envs/seqkit.yaml"
    shell:
        r"""
        seqkit stats {input} > {output}
        """

###############################################################################
# Kraken2 classification
###############################################################################
rule kraken2:
    input:
        r1=W("fastq/{sample}_R1.paired.fastq.gz"),
        r2=W("fastq/{sample}_R2.paired.fastq.gz")
    output:
        report=W("reports_read/kraken/{sample}.report"),
        kraken=W("reports_read/kraken/{sample}.kraken")
    params:
        db=config["kraken"]["DB"]
    conda:
        "../envs/kraken2.yaml"
    threads: 8
    resources:
        mem_mb=450000
    shell:
        r"""
        mkdir -p {WORK_DIR}/kraken

        kraken2 \
            --use-names \
            --threads {threads} \
            --db {params.db} \
            --report {output.report} \
            --paired \
            {input.r1} \
            {input.r2} \
            > {output.kraken}
        """

###############################################################################
# Extract classified read positions
###############################################################################

def existing_raw_bam(wc):
    return W(f"mapping/paired/{wc.sample}_F{FILTER_FLAG}.bam")

rule extract_kraken_classified_positions:
    input:
        kraken=W("reports_read/kraken/{sample}.kraken"),
        bam=existing_raw_bam
    output:
        bed=W("reports_read/kraken/{sample}.classified.positions.bed"),
        tsv=W("reports_read/kraken/{sample}.classified.positions.tsv"),
        stats=W("reports_read/kraken/{sample}.classified.stats.tsv")
    params:
        min_mapq=MAPQ
    conda:
        "../envs/allpurpose.yaml"
    threads: 2
    shell:
        r"""
        /storage/research/ips_pal/GENOMICS/WORK/GENICA/work_Ines/workflow/scripts/script_extract_reads_bed.sh \
            {input.kraken} \
            {input.bam} \
            {output.bed} \
            {output.tsv} \
            {output.stats} \
            {params.min_mapq}
        """

###############################################################################
# Regions covered by Kraken-classified reads in at least N samples
###############################################################################

rule kraken_genomecov:
    input:
        bed=W("reports_read/kraken/{sample}.classified.positions.bed"),
        fai=config["reference"] + ".fai"
    output:
        temp(W("reports_read/kraken/{sample}.classified.genomecov.bed"))
    conda:
        "../envs/allpurpose.yaml"
    shell:
        "bedtools genomecov -i {input.bed} -g {input.fai} -bg > {output}"


rule kraken_filter_coverage:
    input:
        W("reports_read/kraken/{sample}.classified.genomecov.bed")
    output:
        temp(W("reports_read/kraken/{sample}.classified.cov{cov}.bed"))
    conda:
        "../envs/allpurpose.yaml"
    shell:
        r"""awk -v c={wildcards.cov} 'BEGIN{{OFS="\t"}} $4>=c {{print $1,$2,$3}}' {input} | bedtools merge > {output}"""


rule intersect_all_kraken_beds:
    input:
        lambda wc: expand(
            W("reports_read/kraken/{sample}.classified.cov{cov}.bed"),
            sample=SAMPLES,
            cov=wc.cov
        )
    output:
        W("reports_read/kraken/classified.positions.cov{cov}.minsamples{minsamples}.bed")
    conda:
        "../envs/allpurpose.yaml"
    shell:
        r"""bedtools multiinter -i {input} | awk -v n={wildcards.minsamples} 'BEGIN{{OFS="\t"}} $4>=n {{print $1,$2,$3}}' | bedtools merge > {output}"""