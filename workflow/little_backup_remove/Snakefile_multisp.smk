import os
import re
from glob import glob
import pandas as pd

configfile: "config/config_lowcov_xgen.yaml"

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


SAMPLE_TO_SPECIES = dict(
    zip(META["id_code"], META["sp"])
)

SPECIES = META["sp"].unique().tolist()
SAMPLE_SPECIES = META[["id_code", "sp"]].drop_duplicates()

SAMPLES_BY_SPECIES = {
    sp: group["id_code"].tolist()
    for sp, group in META.groupby("sp")
}

REF_META = pd.read_csv(config["references"], sep=r"\s+")

REF_DICT = dict(zip(REF_META["sp"], REF_META["reference"]))

ORGANELLE_META = pd.read_csv(
    config["organelle_contigs_file"],
    sep="\t"
)

def get_organelle_contigs(species):
    return ORGANELLE_META.loc[
        ORGANELLE_META["species"] == species,
        "contig"
    ].dropna().tolist()

UNMASKED_META = pd.read_csv(
    config["unmasked_regions_file"],
    sep="\t"
)

UNMASKED_BED_DICT = dict(
    zip(UNMASKED_META["species"], UNMASKED_META["bed"])
)

READS = ["R1", "R2"]

wildcard_constraints:
    read="R1|R2",
    branch="collapsed|paired",
    sample="|".join(re.escape(s) for s in SAMPLES)

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


print("SAMPLE TEST:", SAMPLES[:5])
print("SAMPLE REGEX:", "|".join(re.escape(s) for s in SAMPLES[:5]))
###############################################################################
# Final outputs
###############################################################################

ALL_TARGETS = []

for sp, samples in SAMPLES_BY_SPECIES.items():

    # Read processing
    ALL_TARGETS += expand(
        W("{species}/reports_read/{sample}_{read}.pre.length.tsv"),
        species=sp, sample=samples, read=READS
    )

    ALL_TARGETS += expand(
        W("{species}/reports_read/{sample}_{read}.post.length.tsv"),
        species=sp, sample=samples, read=READS
    )

    ALL_TARGETS += expand(
        W("{species}/fastq/{sample}_{read}.clean.fastq.gz"),
        species=sp, sample=samples, read=READS
    )

    ALL_TARGETS += expand(
        W("{species}/reports_read/{sample}_R1.paired.length.tsv"),
        species=sp, sample=samples
    )

    ALL_TARGETS += expand(
        W("{species}/reports_read/fastqc/{sample}_{read}.clean_fastqc.html"),
        species=sp, sample=samples, read=READS
    )

    ALL_TARGETS += expand(
        W("{species}/reports_read/{sample}.collapsed.length.tsv"),
        species=sp, sample=samples
    )

    # Merged BAM
    ALL_TARGETS += expand(
        W("{species}/mapping/collapsed_paired/{sample}.bam"),
        species=sp, sample=samples
    )

    # Original mapping flagstat
    ALL_TARGETS += expand(
        W("{species}/reports_mapping/collapsed_paired/{sample}.flagstat.txt"),
        species=sp, sample=samples
    )

    # Kraken classification and positions
    ALL_TARGETS += expand(
        W("{species}/reports_read/kraken/{sample}.classified.positions.{ext}"),
        species=sp, sample=samples, ext=["bed", "tsv"]
    )

    ALL_TARGETS += expand(
        W("{species}/reports_read/kraken/{sample}.classified.stats.tsv"),
        species=sp, sample=samples
    )

    # Mapping statistics after MAPQ filtering
    ALL_TARGETS += expand(
        W(f"{{species}}/reports_mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}.flagstat.txt"),
        species=sp, sample=samples
    )

    # Mapping statistics after Kraken removal
    ALL_TARGETS += expand(
        W(f"{{species}}/reports_mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken.flagstat.txt"),
        species=sp, sample=samples
    )

    # Deduplicated BAM and statistics
    ALL_TARGETS += expand(
        W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.bam"),
        species=sp, sample=samples
    )

    ALL_TARGETS += expand(
        W(f"{{species}}/reports_mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.flagstat.txt"),
        species=sp, sample=samples
    )

    # # Nuclear and organelle BAMs
    # ALL_TARGETS += expand(
    #     W(f"{{species}}/mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_{{genomic_region}}.bam"),
    #     species=sp, sample=samples,
    #     genomic_region=["nuclear", "organelle"]
    # )

    # # Coverage statistics
    # ALL_TARGETS += expand(
    #     W(f"{{species}}/reports_mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_{{genomic_region}}.coverage_summary.tsv"),
    #     species=sp, sample=samples,
    #     genomic_region=["nuclear", "organelle"]
    # )

    # # Coverage distributions
    # ALL_TARGETS += expand(
    #     W(f"{{species}}/reports_mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup_{{genomic_region}}.coverage_histogram.tsv"),
    #     species=sp, sample=samples,
    #     genomic_region=["nuclear", "organelle"]
    # )

    # Mean mapped fragment length and distribution
    ALL_TARGETS += expand(
        W(f"{{species}}/reports_mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.fragment_length_summary.tsv"),
        species=sp, sample=samples
    )

    ALL_TARGETS += expand(
        W(f"{{species}}/reports_mapping/collapsed_paired/{{sample}}_F{FILTER_FLAG}_q{MAPQ}_rmkraken_dedup.fragment_length_distribution.tsv"),
        species=sp, sample=samples
    )


rule all:
    input:
        ALL_TARGETS


###############################################################################
# Rule files
###############################################################################

include: "rules_multisp/read_processing.smk"
include: "rules_multisp/mapping.smk"
include: "rules_multisp/mapping_postprocessing.smk"