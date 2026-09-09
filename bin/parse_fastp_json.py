#!/usr/bin/env python3
"""
Parse fastp JSON reports into a single tab-separated summary table.

fastp writes a rich JSON report per sample, but JSON is awkward to consume
downstream. This script flattens the metrics that matter for metagenomics QC
into one row per sample:

  raw / filtered read counts and bases, retention rate, Q20/Q30 rates,
  GC content, duplication rate, adapter trimming, read length.

Usage:
    parse_fastp_json.py S01.fastp.json S02.fastp.json > fastp_summary.tsv
    parse_fastp_json.py *.fastp.json --output fastp_summary.tsv
"""

import sys
import csv
import json
import argparse
from pathlib import Path

# Column order of the output table. Kept explicit rather than derived from the
# data so the schema stays stable even if fastp adds fields between versions.
COLUMNS = [
    "sample",
    "raw_reads",
    "raw_bases",
    "filtered_reads",
    "filtered_bases",
    "reads_retained_pct",
    "bases_retained_pct",
    "raw_q20_rate",
    "raw_q30_rate",
    "filtered_q20_rate",
    "filtered_q30_rate",
    "raw_gc_content",
    "filtered_gc_content",
    "duplication_rate",
    "adapter_trimmed_reads",
    "adapter_trimmed_bases",
    "insert_size_peak",
    "reads_failed_quality",
    "reads_failed_length",
    "reads_failed_n_bases",
    "reads_too_short",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarise fastp JSON reports into one TSV table"
    )
    parser.add_argument(
        "json_files",
        type=Path,
        nargs="+",
        help="One or more fastp JSON reports",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output TSV path (default: stdout)",
    )
    return parser.parse_args()


def sample_name_from_path(path):
    """Derive the sample ID from a '<sample>.fastp.json' filename."""
    name = path.name
    for suffix in (".fastp.json", ".json"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def pct(numerator, denominator):
    """Percentage, guarding against division by zero on empty inputs."""
    if not denominator:
        return "NA"
    return f"{100.0 * numerator / denominator:.2f}"


def rate(value):
    """fastp reports rates as 0-1 floats; render as a fixed-precision string."""
    if value is None:
        return "NA"
    return f"{float(value):.4f}"


def extract_metrics(report, sample):
    """Flatten one parsed fastp JSON document into a summary row."""
    before = report.get("summary", {}).get("before_filtering", {})
    after = report.get("summary", {}).get("after_filtering", {})
    filtering = report.get("filtering_result", {})
    adapter = report.get("adapter_cutting", {})
    duplication = report.get("duplication", {})

    raw_reads = before.get("total_reads", 0)
    raw_bases = before.get("total_bases", 0)
    filtered_reads = after.get("total_reads", 0)
    filtered_bases = after.get("total_bases", 0)

    return {
        "sample": sample,
        "raw_reads": raw_reads,
        "raw_bases": raw_bases,
        "filtered_reads": filtered_reads,
        "filtered_bases": filtered_bases,
        "reads_retained_pct": pct(filtered_reads, raw_reads),
        "bases_retained_pct": pct(filtered_bases, raw_bases),
        "raw_q20_rate": rate(before.get("q20_rate")),
        "raw_q30_rate": rate(before.get("q30_rate")),
        "filtered_q20_rate": rate(after.get("q20_rate")),
        "filtered_q30_rate": rate(after.get("q30_rate")),
        "raw_gc_content": rate(before.get("gc_content")),
        "filtered_gc_content": rate(after.get("gc_content")),
        # duplication is absent when fastp runs without duplication analysis
        "duplication_rate": rate(duplication.get("rate")) if duplication else "NA",
        "adapter_trimmed_reads": adapter.get("adapter_trimmed_reads", 0),
        "adapter_trimmed_bases": adapter.get("adapter_trimmed_bases", 0),
        "insert_size_peak": report.get("insert_size", {}).get("peak", "NA"),
        "reads_failed_quality": filtering.get("low_quality_reads", 0),
        "reads_failed_length": filtering.get("too_long_reads", 0),
        "reads_failed_n_bases": filtering.get("too_many_N_reads", 0),
        "reads_too_short": filtering.get("too_short_reads", 0),
    }


def main():
    args = parse_args()

    rows = []
    for json_path in args.json_files:
        if not json_path.exists():
            sys.stderr.write(f"ERROR: fastp JSON not found: {json_path}\n")
            sys.exit(1)

        try:
            with open(json_path) as handle:
                report = json.load(handle)
        except json.JSONDecodeError as exc:
            sys.stderr.write(f"ERROR: {json_path} is not valid JSON: {exc}\n")
            sys.exit(1)

        rows.append(extract_metrics(report, sample_name_from_path(json_path)))

    rows.sort(key=lambda row: row["sample"])

    handle = open(args.output, "w", newline="") if args.output else sys.stdout
    try:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    finally:
        if args.output:
            handle.close()


if __name__ == "__main__":
    main()
