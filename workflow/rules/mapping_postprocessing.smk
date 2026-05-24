wildcard_constraints:
    branch="collapsed|paired"

###############################################################################
# SAM -> sorted BAM (remove unmapped and reads with several positions)
###############################################################################

rule sam_to_raw_bam:
    input:
        sam=W("mapping/{branch}/{sample}.sam")
    output:
        bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}.bam"),
        bai=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}.bam.bai")
    params:
        filter_flag=FILTER_FLAG,
        tmp=lambda wc: W(f"tmp/{wc.sample}.{wc.branch}.sorttmp")
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/mapping/{wildcards.branch}
        mkdir -p {WORK_DIR}/tmp

        samtools view -@ {threads} -bh -F {params.filter_flag} {input.sam} \
          | samtools sort -@ {threads} -m 1G -T {params.tmp} -o {output.bam} -

        samtools index -@ {threads} {output.bam}
        """

###############################################################################
# divide organelle and nuclear BAMs
###############################################################################

rule split_organelle_nuclear_bam:
    input:
        bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}.bam"),
        bai=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}.bam.bai")
    output:
        organelle_bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_organelle.bam"),
        nuclear_bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear.bam")
    params:
        organelle_contigs=" ".join(config["organelle_contigs"]),
        organelle_pattern="|".join(c.replace(".", r"\.") for c in config["organelle_contigs"])
    conda:
        "../envs/samtools.yaml"
    threads: 2
    shell:
        r"""
        mkdir -p {WORK_DIR}/mapping/{wildcards.branch}

        samtools view -@ {threads} -b {input.bam} {params.organelle_contigs} > {output.organelle_bam}

        samtools idxstats {input.bam} \
        | cut -f1 \
        | grep -v '^\*$' \
        | grep -v -E '^({params.organelle_pattern})$' \
        > {WORK_DIR}/mapping/{wildcards.branch}/{wildcards.sample}.nuclear_contigs.txt

        samtools view -@ {threads} -b {input.bam} $(cat {WORK_DIR}/mapping/{wildcards.branch}/{wildcards.sample}.nuclear_contigs.txt) > {output.nuclear_bam}
        """

###############################################################################
# remove PCR duplicates
###############################################################################

rule markdup:
    input:
        bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear.bam")
    output:
        bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup.bam"),
        stats=W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear.markdup.stats.txt")
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        r"""
        mkdir -p {WORK_DIR}/mapping/{wildcards.branch}
        mkdir -p {WORK_DIR}/reports_mapping/{wildcards.branch}

        set -euo pipefail

        samtools sort -n -@ {threads} -O bam \
            -T {WORK_DIR}/mapping/{wildcards.branch}/{wildcards.sample}.tmp.namesort \
            {input.bam} \
        | samtools fixmate -m -@ {threads} -O bam - - \
        | samtools sort -@ {threads} -O bam \
            -T {WORK_DIR}/mapping/{wildcards.branch}/{wildcards.sample}.tmp.positionsort \
            - \
        | samtools markdup -r -s -@ {threads} - {output.bam} \
            2> {output.stats}
        """


rule flagstat_markdup:
    input:
        W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup.bam")
    output:
        W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup.bam.flagstat.txt")
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/reports_mapping/{wildcards.branch}
        samtools flagstat {input} > {output}
        """

###############################################################################
# MAPQ filter
###############################################################################

rule MAPQ:
    input:
        W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup.bam")
    output:
        bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.bam"),
        bai=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.bam.bai")
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        r"""
        mkdir -p {WORK_DIR}/mapping/{wildcards.branch}

        samtools view -@ {threads} -bh -q {MAPQ} {input} \
          | samtools sort -@ {threads} -o {output.bam} -

        samtools index -@ {threads} {output.bam}
        """


rule flagstat_MAPQ:
    input:
        W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.bam")
    output:
        W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.bam.flagstat.txt")
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/reports_mapping/{wildcards.branch}
        samtools flagstat {input} > {output}
        """

###############################################################################
# Mean coverage
###############################################################################

rule mean_coverage:
    input:
        bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.bam"),
        bai=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.bam.bai")
    output:
        txt=W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.mean_coverage.txt")
    conda:
        "../envs/samtools.yaml"
    threads: 2
    shell:
        r"""
        mkdir -p {WORK_DIR}/reports_mapping/{wildcards.branch}

        samtools depth -aa {input.bam} | \
        awk '{{sum+=$3; n++}} END {{if(n>0) print sum/n; else print 0}}' > {output.txt}
        """


rule mean_read_length:
    input:
        bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.bam")
    output:
        txt=W(f"reports_read/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.bam.mean_read_length.txt")
    conda:
        "../envs/samtools.yaml"
    threads: 2
    shell:
        r"""
        mkdir -p {WORK_DIR}/reports_read/{wildcards.branch}

        samtools view {input.bam} \
        | awk '{{sum+=length($10); n++}} END {{if(n>0) print sum/n; else print 0}}' \
        > {output.txt}
        """

###############################################################################
# Mapdamage
###############################################################################

rule mapdamage:
    input:
        bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.bam"),
        bai=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.bam.bai"),
        ref=config["reference"]
    output:
        result=W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}/results.txt")
    params:
        outdir=W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}"),
        extra=config.get("mapdamage", {}).get("extra", "")
    conda:
        "../envs/mapdamage.yaml"
    threads: 2
    shell:
        r"""
        mkdir -p {params.outdir}

        mapDamage \
            -i {input.bam} \
            -r {input.ref} \
            -d {params.outdir} \
            {params.extra}
        """