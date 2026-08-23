wildcard_constraints:
    branch="collapsed|paired",
    genomic_region="nuclear|organelle"

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
# MAPQ filter
###############################################################################
rule MAPQ:
    input:
        bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}.bam"),
        bai=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}.bam.bai")
    output:
        bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}.bam"),
        bai=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}.bam.bai")
    conda:
        "../envs/allpurpose.yaml"
    threads: 4
    shell:
        r"""
        mkdir -p {WORK_DIR}/mapping/{wildcards.branch}

        samtools view -@ {threads} -bh -q {MAPQ} {input.bam} \
          | samtools sort -@ {threads} -o {output.bam} -

        samtools index -@ {threads} {output.bam}
        """


rule flagstat_MAPQ:
    input:
        bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}.bam")
    output:
        txt=W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}.bam.flagstat.txt")
    conda:
        "../envs/allpurpose.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/reports_mapping/{wildcards.branch}
        samtools flagstat {input} > {output}
        """

###############################################################################
# Removal of Kraken2-classified contaminant reads
###############################################################################

rule remove_kraken_reads:
    input:
        bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}.bam"),
        tsv=W("reports_read/kraken/{sample}.classified.positions.tsv")
    output:
        bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken.bam"),
        bai=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken.bam.bai")
    conda:
        "../envs/allpurpose.yaml"
    threads: 4
    shell:
        r"""
        set -euo pipefail

        mkdir -p {WORK_DIR}/mapping/{wildcards.branch}

        samtools view -@ {threads} -h {input.bam} \
        | awk -v kraken="{input.tsv}" '
            BEGIN {{
                while ((getline line < kraken) > 0) {{
                    split(line, fields, "\t")
                    if (fields[1] != "read_name")
                        bad[fields[1]] = 1
                }}
                close(kraken)
            }}

            /^@/ {{
                print
                next
            }}

            !($1 in bad) {{
                print
            }}
        ' \
        | samtools view -@ {threads} -b -o {output.bam} -

        samtools index -@ {threads} {output.bam}
        """

rule flagstat_kraken:
    input:
        bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken.bam")
    output:
        txt=W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken.bam.flagstat.txt")
    conda:
        "../envs/allpurpose.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/reports_mapping/{wildcards.branch}
        samtools flagstat {input} > {output}
        """

###############################################################################
# remove PCR duplicates
###############################################################################

rule markdup:
    input:
        bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken.bam")
    output:
        bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.bam"),
        bai=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.bam.bai"),
        stats=W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.markdup.stats.txt")
    conda:
        "../envs/allpurpose.yaml"
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
            
        samtools index -@ {threads} {output.bam}
        """


rule flagstat_markdup:
    input:
        W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.bam")
    output:
        W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.flagstat.txt")
    conda:
        "../envs/allpurpose.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/reports_mapping/{wildcards.branch}
        samtools flagstat {input} > {output}
        """

###############################################################################
# Divide organelle and nuclear BAMs
###############################################################################

rule split_organelle_nuclear_bam:
    input:
        bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.bam"),
        bai=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.bam.bai")
    output:
        organelle_bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_organelle.bam"),
        organelle_bai=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_organelle.bam.bai"),
        nuclear_bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_nuclear.bam"),
        nuclear_bai=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_nuclear.bam.bai")
    params:
        organelle_contigs=" ".join(config["organelle_contigs"]),
        organelle_pattern="|".join(c.replace(".", r"\.") for c in config["organelle_contigs"])
    conda:
        "../envs/allpurpose.yaml"
    threads: 2
    shell:
        r"""
        mkdir -p {WORK_DIR}/mapping/{wildcards.branch}

        samtools view -@ {threads} -b {input.bam} {params.organelle_contigs} > {output.organelle_bam}
        samtools index -@ {threads} {output.organelle_bam}

        samtools idxstats {input.bam} \
        | cut -f1 \
        | grep -v '^\*$' \
        | grep -v -E '^({params.organelle_pattern})$' \
        > {WORK_DIR}/mapping/{wildcards.branch}/{wildcards.sample}.nuclear_contigs.txt

        samtools view -@ {threads} -b {input.bam} $(cat {WORK_DIR}/mapping/{wildcards.branch}/{wildcards.sample}.nuclear_contigs.txt) > {output.nuclear_bam}
        samtools index -@ {threads} {output.nuclear_bam}
        """

###############################################################################
# Coverage statistics across unmasked nuclear and organelle regions
###############################################################################

rule coverage_summary:
    input:
        bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_{{genomic_region}}.bam"),
        bai=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_{{genomic_region}}.bam.bai"),
        masked=config["masked_regions_bed"],
        fai=config["reference"] + ".fai"
    output:
        summary=W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_{{genomic_region}}.coverage_summary.tsv"),
        histogram=W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_{{genomic_region}}.coverage_histogram.tsv")
    params:
        prefix=lambda wc: W(f"tmp/{wc.sample}.{wc.branch}.{wc.genomic_region}"),
        organelle_pattern="|".join(c.replace(".", r"\.") for c in config["organelle_contigs"])
    conda:
        "../envs/allpurpose.yaml"
    threads: 2
    shell:
        r"""
        set -euo pipefail
        mkdir -p {WORK_DIR}/reports_mapping/{wildcards.branch} {WORK_DIR}/tmp

        if [ "{wildcards.genomic_region}" = "organelle" ]; then
            cut -f1,2 {input.fai} | grep -E '^({params.organelle_pattern})[[:space:]]' > {params.prefix}.genome
        else
            cut -f1,2 {input.fai} | grep -v -E '^({params.organelle_pattern})[[:space:]]' > {params.prefix}.genome
        fi

        awk 'NR==FNR{{keep[$1]=1; next}} $1 in keep' \
            {params.prefix}.genome {input.masked} \
        | bedtools sort -i - \
        | bedtools complement -i - -g {params.prefix}.genome \
        > {params.prefix}.unmasked.bed

        printf "sample\tbranch\tgenomic_region\tmean_coverage\tunmasked_region_length_bp\tcovered_at_least_1x_bp\tcovered_at_least_5x_bp\tcovered_at_least_1x_percent\tcovered_at_least_5x_percent\n" \
            > {output.summary}

        samtools depth -@ {threads} -aa -s -b {params.prefix}.unmasked.bed {input.bam} \
        | awk -v sample="{wildcards.sample}" \
              -v branch="{wildcards.branch}" \
              -v region="{wildcards.genomic_region}" \
              -v hist="{params.prefix}.hist" '
            {{
                d=$3; sum+=d; n++;
                if(d>=1) c1++;
                if(d>=5) c5++;
                h[d]++;
            }}
            END {{
                mean=(n ? sum/n : 0);
                p1=(n ? 100*c1/n : 0);
                p5=(n ? 100*c5/n : 0);

                printf "%s\t%s\t%s\t%.8f\t%d\t%d\t%d\t%.8f\t%.8f\n",
                    sample,branch,region,mean,n,c1+0,c5+0,p1,p5;

                for(d in h)
                    print d "\t" h[d] > hist;
            }}' >> {output.summary}

        printf "coverage\tunmasked_length_bp\n" > {output.histogram}
        sort -k1,1n {params.prefix}.hist >> {output.histogram}

        rm -f {params.prefix}.genome {params.prefix}.unmasked.bed {params.prefix}.hist
        """

###############################################################################
# Plot coverage distribution
###############################################################################

rule plot_coverage:
    input:
        histogram=W(
            f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_dedup_{{genomic_region}}.coverage_histogram.tsv")
    output:
        pdf=W(
            f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_dedup_{{genomic_region}}.coverage.pdf")
    params:
        title=lambda wc: (
            f"{wc.sample} – {wc.branch} – {wc.genomic_region}"
        ),
        max_depth=lambda wc: (
            config.get("coverage_plot", {})
            .get("max_depth", {})
            .get(wc.genomic_region, 50)
        )
    conda:
        "../envs/python.yaml"
    shell:
        r"""
        set -euo pipefail

        python ../scripts/plot_coverage.py \
            --histogram {input.histogram} \
            --output {output.pdf} \
            --title "{params.title}" \
            --max-depth {params.max_depth}
        """

###############################################################################
# Mapdamage
###############################################################################

rule mapdamage:
    input:
        bam=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.bam"),
        bai=W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.bam.bai"),
        ref=ancient(config["reference"])
    output:
        runtime=W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}/Runtime_log.txt"),
        misincorp=W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}/misincorporation.txt")
    params:
        outdir=W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}"),
        extra=config.get("mapdamage", {}).get("extra", ""),
        envdir="/storage/research/ips_pal/GENOMICS/WORK/GENICA/work_Ines/workflow/.snakemake/conda/7b5576f7f2ceb98ce826747071521825_"
    log:
        W(f"logs/mapdamage/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.log")
    conda:
        "../envs/mapdamage.yaml"
    threads: 2
    shell:
        r"""
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})

        {params.envdir}/bin/python -c "import sys, pysam; print(sys.executable); print(pysam.__file__)" \
            >> {log} 2>&1

        {params.envdir}/bin/python \
            {params.envdir}/bin/mapDamage \
            -i {input.bam} \
            -r {input.ref} \
            -d {params.outdir} \
            {params.extra} \
            >> {log} 2>&1
        """