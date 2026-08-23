#!/usr/bin/env python3

import argparse
import gzip
import sys
from pathlib import Path
from typing import TextIO


def open_fastq(path: str, mode: str) -> TextIO:
    """Open a FASTQ file, compressed or uncompressed."""
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode, encoding="utf-8")


def read_record(handle: TextIO):
    """Read one four-line FASTQ record."""
    record = [handle.readline() for _ in range(4)]

    if record[0] == "":
        return None

    if any(line == "" for line in record):
        raise ValueError("Incomplete FASTQ record detected.")

    return record


def normalise_read_name(header: str) -> str:
    """
    Return the read identifier without /1, /2 or whitespace-delimited metadata.
    """
    name = header.strip().split()[0]

    if name.endswith("/1") or name.endswith("/2"):
        name = name[:-2]

    return name


def write_record(handle: TextIO, record):
    handle.writelines(record)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Split paired FASTQ files into two synchronized parts. "
            "Read pairs are distributed alternately between part 1 and part 2."
        )
    )

    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--r1-part1", required=True)
    parser.add_argument("--r1-part2", required=True)
    parser.add_argument("--r2-part1", required=True)
    parser.add_argument("--r2-part2", required=True)

    args = parser.parse_args()

    for output_path in (
        args.r1_part1,
        args.r1_part2,
        args.r2_part1,
        args.r2_part2,
    ):
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    pair_count = 0

    with (
        open_fastq(args.r1, "r") as r1_handle,
        open_fastq(args.r2, "r") as r2_handle,
        open_fastq(args.r1_part1, "w") as r1_part1,
        open_fastq(args.r1_part2, "w") as r1_part2,
        open_fastq(args.r2_part1, "w") as r2_part1,
        open_fastq(args.r2_part2, "w") as r2_part2,
    ):
        while True:
            r1_record = read_record(r1_handle)
            r2_record = read_record(r2_handle)

            if r1_record is None and r2_record is None:
                break

            if r1_record is None or r2_record is None:
                raise ValueError(
                    "R1 and R2 contain different numbers of FASTQ records."
                )

            r1_name = normalise_read_name(r1_record[0])
            r2_name = normalise_read_name(r2_record[0])

            if r1_name != r2_name:
                raise ValueError(
                    f"Paired-read names do not match:\n"
                    f"R1: {r1_record[0].strip()}\n"
                    f"R2: {r2_record[0].strip()}"
                )

            if pair_count % 2 == 0:
                write_record(r1_part1, r1_record)
                write_record(r2_part1, r2_record)
            else:
                write_record(r1_part2, r1_record)
                write_record(r2_part2, r2_record)

            pair_count += 1

    if pair_count == 0:
        raise ValueError("The input FASTQ files contain no read pairs.")

    print(
        f"Successfully split {pair_count:,} read pairs:\n"
        f"  part 1: {(pair_count + 1) // 2:,} pairs\n"
        f"  part 2: {pair_count // 2:,} pairs",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()