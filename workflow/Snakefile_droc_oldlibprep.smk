import os
from glob import glob
import pandas as pd

configfile: "config/config_droc_oldlibprep.yaml"

META = pd.read_csv(config["metadata_to_run"], sep="\t")
#Make sure the metadata to run the analysis contains
# id_code seq_id cal_age folder

OLD_TO_NEW = dict(zip(META["seq_id"], META["id_code"]))
NEW_TO_OLD = dict(zip(META["id_code"], META["seq_id"]))
SAMPLE_TO_FOLDER = dict(zip(META["id_code"], META["folder"]))

SAMPLES = list(META["id_code"])
META["cal_age"] = pd.to_numeric(META["cal_age"], errors="coerce")

MODERN_SAMPLES = list(META.loc[META["cal_age"] == 0, "id_code"])
ANCIENT_SAMPLES = list(META.loc[META["cal_age"] > 0, "id_code"])

SAMPLES = MODERN_SAMPLES + ANCIENT_SAMPLES

READS = ["R1", "R2"]
BRANCHES = ["collapsed", "paired"]
ANCIENT_BRANCHES = ["collapsed", "paired"]
MODERN_BRANCHES = ["paired"]

wildcard_constraints:
    read="R1|R2",
    branch="collapsed|paired"

INPUT_DIR = config["input_dir"]
WORK_DIR = config["work_dir"]
FASTQC_OUT = os.path.join(WORK_DIR, "reports_read/fastqc")

FILTER_FLAG = config["samtools"]["filter_flag"]
MAPQ = config["samtools"]["mapq"]


def W(*parts):
    return os.path.join(WORK_DIR, *parts)


def raw_fastqs(sample, read):
    old_sample = NEW_TO_OLD[sample]
    folder = SAMPLE_TO_FOLDER[sample]
    sample_dir = os.path.join(INPUT_DIR, folder)

    if read == "R1":
        pattern = os.path.join(sample_dir, f"*{old_sample}*_R1_*.fastq.gz")
    elif read == "R2":
        pattern = os.path.join(sample_dir, f"*{old_sample}*_R2_*.fastq.gz")
    else:
        raise ValueError(f"Unexpected read: {read}")

    matches = sorted(glob(pattern))

    if not matches:
        raise ValueError(
            f"No raw FASTQs found for sample={sample}, read={read} in {sample_dir}. "
            f"Tried: {pattern}"
        )

    return matches


###############################################################################
# Final outputs
###############################################################################

ALL_TARGETS = []

# read stats
ALL_TARGETS += expand(W("reports_read/{sample}_{read}.pre.length.tsv"), sample=SAMPLES, read=READS)
ALL_TARGETS += expand(W("reports_read/{sample}_{read}.post.length.tsv"), sample=SAMPLES, read=READS)
ALL_TARGETS += expand(W("reports_read/{sample}_R1.paired.length.tsv"), sample=SAMPLES)
ALL_TARGETS += expand(W("reports_read/{sample}.collapsed.length.tsv"), sample=ANCIENT_SAMPLES)

#######
#Ancient DNA
#######
# SAM flagstats
# ALL_TARGETS += expand(W("reports_mapping/{branch}/{sample}.sam.flagstat.txt"),
#                     branch=ANCIENT_BRANCHES, sample=ANCIENT_SAMPLES)
# final BAMs
ALL_TARGETS += expand(W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear.bam"), 
                    branch=ANCIENT_BRANCHES, sample=ANCIENT_SAMPLES)

ALL_TARGETS += expand(W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.bam"),
                    branch=ANCIENT_BRANCHES, sample=ANCIENT_SAMPLES)

# dedup flagstats
ALL_TARGETS += expand(W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup.bam.flagstat.txt"),
                    branch=ANCIENT_BRANCHES, sample=ANCIENT_SAMPLES)

# MAPQ flagstats
ALL_TARGETS += expand(W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.bam.flagstat.txt"),
                    branch=ANCIENT_BRANCHES, sample=ANCIENT_SAMPLES)

# coverage
ALL_TARGETS += expand(W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.mean_coverage.txt"),
                    branch=ANCIENT_BRANCHES, sample=ANCIENT_SAMPLES)

# mean read length
ALL_TARGETS += expand(W(f"reports_read/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.bam.mean_read_length.txt"),
                    branch=ANCIENT_BRANCHES, sample=ANCIENT_SAMPLES)

# # mapDamage
# ALL_TARGETS += expand(W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}/results.txt"),
#                     branch=ANCIENT_BRANCHES, sample=ANCIENT_SAMPLES)

#######
#Modern DNA
#######
# SAM flagstats
# ALL_TARGETS += expand(W("reports_mapping/{branch}/{sample}.sam.flagstat.txt"),
#                     branch=MODERN_BRANCHES, sample=MODERN_SAMPLES)
# final BAMs
ALL_TARGETS += expand(W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear.bam"), 
                    branch=MODERN_BRANCHES, sample=MODERN_SAMPLES)

ALL_TARGETS += expand(W(f"mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.bam"),
                    branch=MODERN_BRANCHES, sample=MODERN_SAMPLES)

# dedup flagstats
ALL_TARGETS += expand(W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup.bam.flagstat.txt"),
                      branch=MODERN_BRANCHES, sample=MODERN_SAMPLES)

# MAPQ flagstats
ALL_TARGETS += expand(W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.bam.flagstat.txt"),
                      branch=MODERN_BRANCHES, sample=MODERN_SAMPLES)

# coverage
ALL_TARGETS += expand(W(f"reports_mapping/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.mean_coverage.txt"),
                      branch=MODERN_BRANCHES, sample=MODERN_SAMPLES)

# mean read length
ALL_TARGETS += expand(W(f"reports_read/{{branch}}/{{sample}}_F{FILTER_FLAG}_nuclear_dedup_q{MAPQ}.bam.mean_read_length.txt"),
                      branch=MODERN_BRANCHES, sample=MODERN_SAMPLES)

rule all:
    input:
        ALL_TARGETS


###############################################################################
# Rule files
###############################################################################

include: "rules/read_processing.smk"
include: "rules/mapping_collapsed.smk"
include: "rules/mapping_paired.smk"
include: "rules/mapping_postprocessing.smk"