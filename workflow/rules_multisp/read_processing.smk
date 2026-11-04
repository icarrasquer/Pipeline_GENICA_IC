###############################################################################
# merge raw reads into a temp fastq (uncompressed)
###############################################################################
rule merge_raw:
    input:
        lambda wc: raw_fastqs(wc.sample, wc.read)
    output:
        temp(W("{species}/tmp/{sample}_{read}.merged.fastq"))
    shell:
        r"""
        mkdir -p {WORK_DIR}/{wildcards.species}/tmp
        zcat {input} > {output}
        """

###############################################################################
# initial stats (write per-sample/read stats file)
###############################################################################
rule stats_pre:
    input:
        W("{species}/tmp/{sample}_{read}.merged.fastq")
    output:
        W("{species}/reports_read/{sample}_{read}.pre.stats.tsv"),
        W("{species}/reports_read/{sample}_{read}.pre.length.tsv")
    params:
        a3 = config["adapters"]["a_3prime"],
        g5 = config["adapters"]["g_5prime_loop"]
    conda:
        "../envs/seqkit.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/{wildcards.species}/reports_read
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
        W("{species}/tmp/{sample}_{read}.merged.fastq")
    output:
        temp(W("{species}/tmp/{sample}_{read}.cutadapt1.fastq"))
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
        mkdir -p {WORK_DIR}/{wildcards.species}/tmp

        cutadapt -j {threads} -e {params.e} -a {params.a3} {input} -o {output} -m {params.minlen}

        num=$(grep -c {params.g5} {output} || echo 0)
        while [ "$num" -gt {params.maxhits} ]; do
            cutadapt -j {threads} -e {params.e} -g {params.g5} {output} -o {WORK_DIR}/{wildcards.species}/tmp/{wildcards.sample}_{wildcards.read}.loop.fastq -m {params.minlen}
            mv {WORK_DIR}/{wildcards.species}/tmp/{wildcards.sample}_{wildcards.read}.loop.fastq {output}
            num=$(grep -c {params.g5} {output} || echo 0)
        done
        """

###############################################################################
# polyG/polyC removal via perl script
###############################################################################
rule drop_poly:
    input:
        W("{species}/tmp/{sample}_{read}.cutadapt1.fastq")
    output:
        temp(W("{species}/tmp/{sample}_{read}.poly.fastq"))
    shell:
        r"""
        mkdir -p {WORK_DIR}/{wildcards.species}/tmp
        DropBpFastq_polyC.pl {input} {output}
        """

###############################################################################
# quality + length filter -> final cleaned fastq.gz
###############################################################################
rule quality_filter:
    input:
        W("{species}/tmp/{sample}_{read}.poly.fastq")
    output:
        W("{species}/fastq/{sample}_{read}.clean.fastq.gz")
    conda:
        "../envs/cutadapt.yaml"
    threads:
        config["cutadapt"]["threads"]
    params:
        q = config["cutadapt"]["qual"],
        minlen = config["cutadapt"]["minlen"]
    shell:
        r"""
        mkdir -p {WORK_DIR}/{wildcards.species}/fastq {WORK_DIR}/{wildcards.species}/tmp

        cutadapt -j {threads} -q {params.q} {input} \
            -o {WORK_DIR}/{wildcards.species}/tmp/{wildcards.sample}_{wildcards.read}.clean.fastq \
            -m {params.minlen}

        gzip -c {WORK_DIR}/{wildcards.species}/tmp/{wildcards.sample}_{wildcards.read}.clean.fastq > {output}

        rm -f {WORK_DIR}/{wildcards.species}/tmp/{wildcards.sample}_{wildcards.read}.clean.fastq
        """

###############################################################################
# post stats per read (again per-sample/read files)
###############################################################################
rule stats_post:
    input:
        W("{species}/fastq/{sample}_{read}.clean.fastq.gz")
    output:
        W("{species}/reports_read/{sample}_{read}.post.stats.tsv"),
        W("{species}/reports_read/{sample}_{read}.post.length.tsv")
    params:
        a3 = config["adapters"]["a_3prime"],
        g5 = config["adapters"]["g_5prime_loop"]
    conda:
        "../envs/seqkit.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/{wildcards.species}/reports_read
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
        W("{species}/fastq/{sample}_{read}.clean.fastq.gz")
    output:
        html=W("{species}/reports_read/fastqc/{sample}_{read}.clean_fastqc.html"),
        zip=W("{species}/reports_read/fastqc/{sample}_{read}.clean_fastqc.zip")
    conda:
        "../envs/fastqc.yaml"
    params:
        outdir=lambda wc: W(f"{wc.species}/reports_read/fastqc")
    shell:
        r"""
        mkdir -p {params.outdir}
        fastqc {input} -o {params.outdir}
        """

###############################################################################
# repair reads
###############################################################################
rule repair_reads:
    input:
        r1=W("{species}/fastq/{sample}_R1.clean.fastq.gz"),
        r2=W("{species}/fastq/{sample}_R2.clean.fastq.gz")
    output:
        r1p=W("{species}/fastq/{sample}_R1.paired.fastq.gz"),
        r2p=W("{species}/fastq/{sample}_R2.paired.fastq.gz")
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
        W("{species}/fastq/{sample}_R1.paired.fastq.gz")
    output:
        W("{species}/reports_read/{sample}_R1.paired.length.tsv")
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
        r1=W("{species}/fastq/{sample}_R1.paired.fastq.gz"),
        r2=W("{species}/fastq/{sample}_R2.paired.fastq.gz")
    output:
        report=W("{species}/reports_read/kraken/{sample}.report"),
        kraken=W("{species}/reports_read/kraken/{sample}.kraken")
    params:
        db=config["kraken"]["DB"]
    conda:
        "../envs/kraken2.yaml"
    threads: 8
    resources:
        mem_mb=450000
    shell:
        r"""
        mkdir -p {WORK_DIR}/{wildcards.species}/reports_read/kraken

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
    return W(f"{wc.species}/mapping/paired/{wc.sample}_F{FILTER_FLAG}.bam")

rule extract_kraken_classified_positions:
    input:
        kraken=W("{species}/reports_read/kraken/{sample}.kraken"),
        bam=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}.bam")
    output:
        bed=W("{species}/reports_read/kraken/{sample}.classified.positions.bed"),
        tsv=W("{species}/reports_read/kraken/{sample}.classified.positions.tsv"),
        stats=W("{species}/reports_read/kraken/{sample}.classified.stats.tsv")
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
        bed=W("{species}/reports_read/kraken/{sample}.classified.positions.bed"),
        fai=lambda wc: REF_DICT[wc.species] + ".fai"
    output:
        temp(W("{species}/reports_read/kraken/{sample}.classified.genomecov.bed"))
    conda:
        "../envs/allpurpose.yaml"
    shell:
        "bedtools genomecov -i {input.bed} -g {input.fai} -bg > {output}"


rule kraken_filter_coverage:
    input:
        W("{species}/reports_read/kraken/{sample}.classified.genomecov.bed")
    output:
        temp(W("{species}/reports_read/kraken/{sample}.classified.cov{cov}.bed"))
    conda:
        "../envs/allpurpose.yaml"
    shell:
        r"""awk -v c={wildcards.cov} 'BEGIN{{OFS="\t"}} $4>=c {{print $1,$2,$3}}' {input} | bedtools merge > {output}"""


rule intersect_all_kraken_beds:
    input:
        lambda wc: expand(
            W("{species}/reports_read/kraken/{sample}.classified.cov{cov}.bed"),
            species=wc.species,
            sample=SAMPLES_BY_SPECIES[wc.species],
            cov=wc.cov
        )
    output:
        W("{species}/reports_read/kraken/classified.positions.cov{cov}.minsamples{minsamples}.bed")
    conda:
        "../envs/allpurpose.yaml"
    shell:
        r"""bedtools multiinter -i {input} | awk -v n={wildcards.minsamples} 'BEGIN{{OFS="\t"}} $4>=n {{print $1,$2,$3}}' | bedtools merge > {output}"""

###############################################################################
# AdapterRemoval merge
###############################################################################
rule adapterremoval_collapse:
    input:
        r1p=W("{species}/fastq/{sample}_R1.paired.fastq.gz"),
        r2p=W("{species}/fastq/{sample}_R2.paired.fastq.gz")
    output:
        merged=W("{species}/fastq/{sample}.collapsed.fastq.gz"),
        r1=W("{species}/fastq/{sample}.R1_uncoll.fastq.gz"),
        r2=W("{species}/fastq/{sample}.R2_uncoll.fastq.gz")
    conda:
        "../envs/adapterremoval3.yaml"
    threads:
        config["adapterremoval"]["threads"]
    params:
        minq=config["adapterremoval"]["minquality"],
        minlen=config["adapterremoval"]["minlength"]
    shell:
        r"""
        adapterremoval \
            --in-file1 {input.r1p} \
            --in-file2 {input.r2p} \
            --out-format fastq.gz \
            --threads {threads} \
            --trimns \
            --trimqualities \
            --minquality {params.minq} \
            --minlength {params.minlen} \
            --merge \
            --out-merged {output.merged} \
            --out-file1 {output.r1} \
            --out-file2 {output.r2} 
        """
###############################################################################
# post stats per read (again per-sample/read files)
###############################################################################
rule stats_collapsed:
    input:
        W("{species}/fastq/{sample}.collapsed.fastq.gz")
    output:
        W("{species}/reports_read/{sample}.collapsed.length.tsv")
    conda:
        "../envs/seqkit.yaml"
    shell:
        r"""
        seqkit stats {input} > {output}
        """
