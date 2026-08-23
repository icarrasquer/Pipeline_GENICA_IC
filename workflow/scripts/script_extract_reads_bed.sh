#!/bin/bash
set -euo pipefail

KRAKEN="$1"
BAM="$2"
OUTPUT_BED="$3"
OUTPUT_TSV="$4"
OUTPUT_STATS="$5"
MIN_MAPQ="$6"

TMP_PREFIX="${OUTPUT_BED%.bed}"
READ_NAMES="${TMP_PREFIX}.classified.names.txt"
TMP_BAM="${TMP_PREFIX}.classified.tmp.bam"

mkdir -p "$(dirname "$OUTPUT_BED")"
mkdir -p "$(dirname "$OUTPUT_TSV")"
mkdir -p "$(dirname "$OUTPUT_STATS")"

###############################################################################
# Extract names of read pairs classified by Kraken
###############################################################################

awk '$1 == "C" {print $2}' "$KRAKEN" \
    | sort -u \
    > "$READ_NAMES"

###############################################################################
# Extract Kraken-classified alignments with MAPQ >= threshold
###############################################################################

samtools view \
    -b \
    -q "$MIN_MAPQ" \
    -N "$READ_NAMES" \
    "$BAM" \
    > "$TMP_BAM"

###############################################################################
# Create BED file with the mapped positions
###############################################################################

bedtools bamtobed -i "$TMP_BAM" \
    | sort -k1,1 -k2,2n \
    > "$OUTPUT_BED"

###############################################################################
# Create read-level alignment table
###############################################################################

samtools view "$TMP_BAM" \
    | awk '
        BEGIN {
            OFS = "\t"
            print "read_name", "chromosome", "position", "mapq", "cigar"
        }
        {
            print $1, $3, $4, $5, $6
        }
    ' \
    > "$OUTPUT_TSV"

###############################################################################
# Count read pairs considered by Kraken
#
# Each Kraken output line represents one paired-end fragment.
###############################################################################

TOTAL_KRAKEN_PAIRS=$(
    awk 'END {print NR + 0}' "$KRAKEN"
)

TOTAL_KRAKEN_READS=$((TOTAL_KRAKEN_PAIRS * 2))

###############################################################################
# Count Kraken-classified pairs
###############################################################################

KRAKEN_CLASSIFIED_PAIRS=$(
    awk '$1 == "C" {count++} END {print count + 0}' "$KRAKEN"
)

###############################################################################
# Count individual R1/R2 alignments retained at the selected MAPQ
###############################################################################

CLASSIFIED_MAPPED_READS=$(
    samtools view -c "$TMP_BAM"
)

###############################################################################
# Write statistics
###############################################################################

awk \
    -v total_pairs="$TOTAL_KRAKEN_PAIRS" \
    -v total_reads="$TOTAL_KRAKEN_READS" \
    -v classified_pairs="$KRAKEN_CLASSIFIED_PAIRS" \
    -v mapped_reads="$CLASSIFIED_MAPPED_READS" \
    -v min_mapq="$MIN_MAPQ" \
    '
    BEGIN {
        OFS = "\t"

        if (total_pairs > 0) {
            classified_pairs_percent = 100 * classified_pairs / total_pairs
        } else {
            classified_pairs_percent = 0
        }

        if (total_reads > 0) {
            mapped_reads_percent = 100 * mapped_reads / total_reads
        } else {
            mapped_reads_percent = 0
        }

        print \
            "total_kraken_pairs", \
            "total_kraken_reads", \
            "kraken_classified_pairs", \
            "kraken_classified_pairs_percent", \
            "classified_mapped_reads_min_mapq", \
            "classified_mapped_reads_min_mapq_percent", \
            "min_mapq"

        printf \
            "%d\t%d\t%d\t%.8f\t%d\t%.8f\t%d\n", \
            total_pairs, \
            total_reads, \
            classified_pairs, \
            classified_pairs_percent, \
            mapped_reads, \
            mapped_reads_percent, \
            min_mapq
    }
    ' \
    > "$OUTPUT_STATS"

###############################################################################
# Remove temporary files
###############################################################################

rm -f "$READ_NAMES" "$TMP_BAM"