#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot the percentage of the unmasked genomic region "
            "covered at or above each sequencing depth."
        )
    )

    parser.add_argument(
        "--histogram",
        required=True,
        help="Coverage histogram produced by coverage_summary.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output PDF path.",
    )

    parser.add_argument(
        "--title",
        required=True,
        help="Plot title.",
    )

    parser.add_argument(
        "--max-depth",
        type=int,
        default=50,
        help="Maximum depth displayed on the x-axis.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    histogram_path = Path(args.histogram)
    output_path = Path(args.output)

    if not histogram_path.exists():
        raise FileNotFoundError(
            f"Histogram file does not exist: {histogram_path}"
        )

    df = pd.read_csv(histogram_path, sep="\t")

    required_columns = {
        "coverage",
        "unmasked_length_bp",
    }

    if not required_columns.issubset(df.columns):
        raise ValueError(
            "Histogram must contain coverage and "
            "unmasked_length_bp columns."
        )

    df["coverage"] = pd.to_numeric(
        df["coverage"],
        errors="raise",
    )

    df["unmasked_length_bp"] = pd.to_numeric(
        df["unmasked_length_bp"],
        errors="raise",
    )

    df = (
        df.groupby("coverage", as_index=False)[
            "unmasked_length_bp"
        ]
        .sum()
        .sort_values("coverage")
    )

    total_unmasked_length = df["unmasked_length_bp"].sum()

    if total_unmasked_length <= 0:
        raise ValueError(
            "The total unmasked region length is zero."
        )

    df["covered_at_least_depth_bp"] = (
        df["unmasked_length_bp"]
        .iloc[::-1]
        .cumsum()
        .iloc[::-1]
    )

    df["covered_at_least_depth_percent"] = (
        100
        * df["covered_at_least_depth_bp"]
        / total_unmasked_length
    )

    max_depth = max(1, args.max_depth)

    plot_df = df[df["coverage"] <= max_depth].copy()

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.step(
        plot_df["coverage"],
        plot_df["covered_at_least_depth_percent"],
        where="post",
    )

    ax.axvline(
        1,
        linestyle="--",
        linewidth=1,
    )

    ax.axvline(
        5,
        linestyle="--",
        linewidth=1,
    )

    ax.set_title(args.title)
    ax.set_xlabel("Minimum coverage depth")
    ax.set_ylabel(
        "Unmasked region covered at or above depth (%)"
    )

    ax.set_xlim(0, max_depth)
    ax.set_ylim(0, 100)

    ax.grid(
        axis="y",
        alpha=0.3,
    )

    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


if __name__ == "__main__":
    main()