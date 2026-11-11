wildcard_constraints:
    genomic_region="nuclear|organelle"

###############################################################################
# SAM -> sorted BAM (remove unmapped and reads with several positions)
###############################################################################

rule sam_to_raw_bam:
    input:
        bam=W("{species}/mapping/collapsed_paired/{sample}.bam")
    output:
        bam=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}.bam"),
        csi=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}.bam.csi")
    params:
        filter_flag=FILTER_FLAG,
        tmp=lambda wc: W(f"{wc.species}/tmp/{wc.sample}.collapsed_paired.sorttmp")
    threads: 4
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/{wildcards.species}/mapping/collapsed_paired
        mkdir -p {WORK_DIR}/{wildcards.species}/tmp

        samtools view -@ {threads} -bh -F {params.filter_flag} {input.bam} \
          -e '!(flag & 1) || (flag & 2)' | samtools sort -@ {threads} -m 1G -T {params.tmp} -o {output.bam} -

        samtools index -c -@ {threads} {output.bam}
        """

###############################################################################
# MAPQ filter
###############################################################################
rule MAPQ:
    input:
        bam=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}.bam"),
        csi=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}.bam.csi")
    output:
        bam=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}.bam"),
        csi=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}.bam.csi")
    conda:
        "../envs/allpurpose.yaml"
    threads: 4
    shell:
        r"""
        mkdir -p {WORK_DIR}/{wildcards.species}/mapping/collapsed_paired

        samtools view -@ {threads} -bh -q {MAPQ} {input.bam} \
          | samtools sort -@ {threads} -o {output.bam} -

        samtools index -c -@ {threads} {output.bam}
        """


rule flagstat_MAPQ:
    input:
        bam=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}.bam")
    output:
        txt=W(f"{{species}}/reports_mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}.flagstat.txt")
    conda:
        "../envs/allpurpose.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/{wildcards.species}/reports_mapping/collapsed_paired
        samtools flagstat {input} > {output}
        """

###############################################################################
# Removal of Kraken2-classified contaminant reads
###############################################################################

rule remove_kraken_reads:
    input:
        bam=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}.bam"),
        tsv=W("{species}/reports_read/kraken/{sample}.classified.positions.tsv")
    output:
        bam=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken.bam"),
        csi=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken.bam.csi")
    conda:
        "../envs/allpurpose.yaml"
    threads: 4
    shell:
        r"""
        set -euo pipefail

        mkdir -p {WORK_DIR}/{wildcards.species}/mapping/collapsed_paired

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

        samtools index -c -@ {threads} {output.bam}
        """

rule flagstat_kraken:
    input:
        bam=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken.bam")
    output:
        txt=W(f"{{species}}/reports_mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken.flagstat.txt")
    conda:
        "../envs/allpurpose.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/{wildcards.species}/reports_mapping/collapsed_paired
        samtools flagstat {input} > {output}
        """

###############################################################################
# Remove PCR duplicates
###############################################################################

rule markdup:
    input:
        bam=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken.bam")
    output:
        bam=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.bam"),
        csi=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.bam.csi"),
        stats=W(f"{{species}}/reports_mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.markdup.stats.txt")
    conda:
        "../envs/allpurpose.yaml"
    threads: 4
    shell:
        r"""
        mkdir -p {WORK_DIR}/{wildcards.species}/mapping/collapsed_paired
        mkdir -p {WORK_DIR}/{wildcards.species}/reports_mapping/collapsed_paired

        set -euo pipefail

        samtools sort -n -@ {threads} -O bam \
            -T {WORK_DIR}/{wildcards.species}/mapping/collapsed_paired/{wildcards.sample}.tmp.namesort \
            {input.bam} \
        | samtools fixmate -m -@ {threads} -O bam - - \
        | samtools sort -@ {threads} -O bam \
            -T {WORK_DIR}/{wildcards.species}/mapping/collapsed_paired/{wildcards.sample}.tmp.positionsort \
            - \
        | samtools markdup -r -s -@ {threads} - {output.bam} \
            2> {output.stats}
            
        samtools index -c -@ {threads} {output.bam}
        """


rule flagstat_markdup:
    input:
        W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.bam")
    output:
        W(f"{{species}}/reports_mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.flagstat.txt")
    conda:
        "../envs/allpurpose.yaml"
    shell:
        r"""
        mkdir -p {WORK_DIR}/{wildcards.species}/reports_mapping/collapsed_paired
        samtools flagstat {input} > {output}
        """

###############################################################################
# Divide organelle and nuclear BAMs
###############################################################################

rule split_organelle_nuclear_bam:
    input:
        bam=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.bam"),
        csi=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.bam.csi")
    output:
        organelle_bam=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_organelle.bam"),
        organelle_csi=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_organelle.bam.csi"),
        nuclear_bam=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_nuclear.bam"),
        nuclear_csi=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_nuclear.bam.csi")
    params:
        organelle_contigs=lambda wc: " ".join(
            get_organelle_contigs(ref_species(wc.species))
        ),
        organelle_pattern=lambda wc: "|".join(
            re.escape(c)
            for c in get_organelle_contigs(ref_species(wc.species))
        )
    conda:
        "../envs/allpurpose.yaml"
    threads: 2
    shell:
        r"""
        mkdir -p {WORK_DIR}/{wildcards.species}/mapping/collapsed_paired

        samtools view -@ {threads} -b {input.bam} {params.organelle_contigs} > {output.organelle_bam}
        samtools index -c -@ {threads} {output.organelle_bam}

        samtools idxstats {input.bam} \
        | cut -f1 \
        | grep -v '^\*$' \
        | grep -v -E '^({params.organelle_pattern})$' \
        > {WORK_DIR}/{wildcards.species}/mapping/collapsed_paired/{wildcards.sample}.nuclear_contigs.txt

        samtools view -@ {threads} -b {input.bam} $(cat {WORK_DIR}/{wildcards.species}/mapping/collapsed_paired/{wildcards.sample}.nuclear_contigs.txt) > {output.nuclear_bam}
        samtools index -c -@ {threads} {output.nuclear_bam}
        """

###############################################################################
# Coverage statistics across unmasked nuclear and organelle regions
###############################################################################

rule coverage_summary:
    input:
        bam=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_{{genomic_region}}.bam"),
        csi=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_{{genomic_region}}.bam.csi"),
        regions=lambda wc: UNMASKED_BED_DICT[ref_species(wc.species)]
    output:
        summary=W(f"{{species}}/reports_mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_{{genomic_region}}.coverage_summary.tsv"),
        histogram=W(f"{{species}}/reports_mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_{{genomic_region}}.coverage_histogram.tsv")
    params:
        prefix=lambda wc: W(f"{wc.species}/tmp/{wc.sample}.collapsed_paired.{wc.genomic_region}"),
        organelle_pattern=lambda wc: "|".join(
            re.escape(c)
            for c in get_organelle_contigs(ref_species(wc.species))
        )
    conda:
        "../envs/allpurpose.yaml"
    threads: 2
    shell:
        r"""
        set -euo pipefail
        mkdir -p {WORK_DIR}/{wildcards.species}/reports_mapping/collapsed_paired {WORK_DIR}/{wildcards.species}/tmp

        if [ "{wildcards.genomic_region}" = "organelle" ]; then
            grep -E '^({params.organelle_pattern})[[:space:]]' \
                {input.regions} > {params.prefix}.regions.bed
        else
            grep -v -E '^({params.organelle_pattern})[[:space:]]' \
                {input.regions} > {params.prefix}.regions.bed
        fi

        printf "sample\tbranch\tgenomic_region\tmean_coverage\tunmasked_region_length_bp\tcovered_at_least_1x_bp\tcovered_at_least_5x_bp\tcovered_at_least_1x_percent\tcovered_at_least_5x_percent\n" \
            > {output.summary}

        samtools depth -@ {threads} -a -s \
            -b {params.prefix}.regions.bed \
            {input.bam} \
        | awk \
            -v sample="{wildcards.sample}" \
            -v branch="collapsed_paired" \
            -v region="{wildcards.genomic_region}" \
            -v hist="{params.prefix}.hist" '
            {{
                d=$3
                sum+=d
                n++
                if(d>=1) c1++
                if(d>=5) c5++
                h[d]++
            }}
            END {{
                mean=(n ? sum/n : 0)
                p1=(n ? 100*c1/n : 0)
                p5=(n ? 100*c5/n : 0)

                printf "%s\t%s\t%s\t%.8f\t%d\t%d\t%d\t%.8f\t%.8f\n",
                    sample,branch,region,mean,n,c1+0,c5+0,p1,p5

                for(d in h)
                    print d "\t" h[d] > hist
            }}' >> {output.summary}

        printf "coverage\tunmasked_length_bp\n" > {output.histogram}
        sort -k1,1n {params.prefix}.hist >> {output.histogram}

        rm -f {params.prefix}.regions.bed {params.prefix}.hist
        """

###############################################################################
# Plot coverage distribution
###############################################################################

rule plot_coverage:
    input:
        histogram=W(
            f"{{species}}/reports_mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_{{genomic_region}}.coverage_histogram.tsv"
        )
    output:
        pdf=W(
            f"{{species}}/reports_mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_{{genomic_region}}.coverage.pdf"
        )
    params:
        title=lambda wc: f"{wc.sample} – collapsed_paired – {wc.genomic_region}",
        max_depth=lambda wc: (
            config.get("coverage_plot", {})
            .get("max_depth", {})
            .get(wc.genomic_region, 50)
        )
    conda:
        "../envs/allpurpose.yaml"
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
# Mean mapped fragment length and fragment length distribution
###############################################################################

rule mapped_fragment_length:
    input:
        bam=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.bam"),
        csi=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.bam.csi")
    output:
        summary=W(f"{{species}}/reports_mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.fragment_length_summary.tsv"),
        distribution=W(f"{{species}}/reports_mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.fragment_length_distribution.tsv")
    conda:
        "../envs/allpurpose.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {WORK_DIR}/{wildcards.species}/reports_mapping/collapsed_paired

        samtools view -F 4 {input.bam} \
        | awk '
            and($2, 2304) == 0 && and($2, 1) != 0 && and($2, 8) == 0 && $9 > 0 {{print $9}}
            and($2, 2304) == 0 && and($2, 1) == 0 {{print length($10)}}
        ' > {output.distribution}.tmp


        printf "sample\tbranch\tn_fragments\tmean_fragment_length\n" > {output.summary}
        awk -v s="{wildcards.sample}" -v b="collapsed_paired" \
            '{{sum+=$1; n++}} END{{printf "%s\t%s\t%d\t%.8f\n",s,b,n,(n?sum/n:0)}}' \
            {output.distribution}.tmp >> {output.summary}

        printf "fragment_length_bp\tcount\n" > {output.distribution}
        sort -n {output.distribution}.tmp | uniq -c \
        | awk '{{print $2 "\t" $1}}' >> {output.distribution}

        rm -f {output.distribution}.tmp
        """
        
###############################################################################
# Mapdamage
###############################################################################

rule mapdamage:
    input:
        bam=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_{{genomic_region}}.bam"),
        csi=W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_{{genomic_region}}.bam.csi"),
        ref=lambda wc: REF_DICT[ref_species(wc.species)]
    output:
        runtime=W(f"{{species}}/reports_mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_{{genomic_region}}/Runtime_log.txt"),
        misincorp=W(f"{{species}}/reports_mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_{{genomic_region}}/misincorporation.txt")
    params:
        outdir=W(f"{{species}}/reports_mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_{{genomic_region}}"),
        extra=config.get("mapdamage", {}).get("extra", ""),
        envdir="/storage/research/ips_pal/GENOMICS/WORK/GENICA/work_Ines/workflow/.snakemake/conda/7b5576f7f2ceb98ce826747071521825_"
    log:
        W(f"{{species}}/logs/mapdamage/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_{{genomic_region}}.log")
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