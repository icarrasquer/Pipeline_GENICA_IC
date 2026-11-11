import os
import re
from glob import glob
import pandas as pd

configfile: "config/config_lowcov_xgen.yaml"

META = pd.read_csv(config["metadata_to_run"], sep="\t")
for col in ["id_code", "seq_id", "sp", "folder"]: META[col] = META[col].astype(str).str.strip()
META["read_sample"] = META["id_code"].map(lambda s: re.sub(r"_[^_]+_wgs$", "", s) if s.startswith("EBC_") else s)
META["read_group"] = META.apply(lambda r: "EBC" if r.id_code.startswith("EBC_") else r.sp, axis=1)
META["mapping_group"] = META.apply(lambda r: f"EBC/{r.sp}" if r.id_code.startswith("EBC_") else r.sp, axis=1)
READ_META = META.drop_duplicates("read_sample")
NEW_TO_OLD = dict(zip(READ_META["read_sample"], READ_META["seq_id"]))
OLD_TO_NEW = dict(zip(READ_META["seq_id"], READ_META["read_sample"]))
SAMPLE_TO_FOLDER = dict(zip(READ_META["read_sample"], READ_META["folder"]))
READ_SAMPLE = dict(zip(META["id_code"], META["read_sample"]))
READ_SPECIES = dict(zip(META["id_code"], META["read_group"]))
def read_path(sample, filename): return W(READ_SPECIES[sample], filename.format(sample=READ_SAMPLE[sample]))
SAMPLES = META["id_code"].tolist()
META["cal_age"] = pd.to_numeric(META["cal_age"], errors="coerce")
MODERN_SAMPLES = META.loc[META["cal_age"] == 0, "id_code"].tolist()
ANCIENT_SAMPLES = META.loc[META["cal_age"] > 0, "id_code"].tolist()
SAMPLE_TO_SPECIES = dict(zip(META["id_code"], META["sp"]))
SPECIES = META["sp"].unique().tolist()
SAMPLE_SPECIES = META[["id_code", "sp"]].drop_duplicates()
SAMPLES_BY_SPECIES = {sp: g["id_code"].tolist() for sp, g in META.groupby("mapping_group")}
READ_SAMPLES_BY_SPECIES = {sp: g["read_sample"].tolist() for sp, g in READ_META.groupby("read_group")}
def ref_species(species): return species.split("/")[-1]

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
    sample="|".join(re.escape(s) for s in SAMPLES + READ_META["read_sample"].tolist()),
    species="|".join(re.escape(s) for s in list(SAMPLES_BY_SPECIES) + list(READ_SAMPLES_BY_SPECIES))

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

for sp, read_samples in READ_SAMPLES_BY_SPECIES.items():
    # Read processing
    ALL_TARGETS += expand(
        W("{species}/reports_read/{sample}_{read}.pre.length.tsv"),
        species=sp, sample=read_samples, read=READS
    )

    ALL_TARGETS += expand(
        W("{species}/reports_read/{sample}_{read}.post.length.tsv"),
        species=sp, sample=read_samples, read=READS
    )

    ALL_TARGETS += expand(
        W("{species}/fastq/{sample}_{read}.clean.fastq.gz"),
        species=sp, sample=read_samples, read=READS
    )

    ALL_TARGETS += expand(
        W("{species}/reports_read/{sample}_R1.paired.length.tsv"),
        species=sp, sample=read_samples
    )

    ALL_TARGETS += expand(
        W("{species}/reports_read/fastqc/{sample}_{read}.clean_fastqc.html"),
        species=sp, sample=read_samples, read=READS
    )

    ALL_TARGETS += expand(
        W("{species}/reports_read/{sample}.collapsed.length.tsv"),
        species=sp, sample=read_samples
    )


for sp, samples in SAMPLES_BY_SPECIES.items():

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