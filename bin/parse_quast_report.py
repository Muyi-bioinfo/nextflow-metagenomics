#!/usr/bin/env python3
"""
Parse QUAST report.tsv files into a single tab-separated assembly summary.

QUAST writes one report per assembly, which is fine for reading but useless for
comparing assemblies. This script flattens them into one row per assembly:

  contig counts, total length, largest contig, N50/N75/N90, L50/L75/L90, auN,
  GC content, contig counts and bases above each length threshold, N density.

Two things make this more than a transpose:

  * QUAST metric names are prose ("# contigs (>= 1000 bp)") and vary with the
    options used. They are mapped onto a fixed column schema so the output stays
    stable across QUAST versions, with "NA" where a metric was not reported.
    N90/L90/auN in particular only appear with --report-all-metrics (QUAST 5.2+).

  * The assembly label carries both identities the pipeline cares about —
    "<sample>.<assembler>" — so sample and assembler are split back out into
    their own columns. With --assembler both, a sample yields two rows that can
    be compared directly.

Usage:
    parse_quast_report.py S01.megahit.quast.report.tsv > assembly_summary.tsv
    parse_quast_report.py *.quast.report.tsv --output assembly_summary.tsv
"""

import sys
import csv
import argparse
from pathlib import Path

# Assemblers the pipeline can run. Used to split "<sample>.<assembler>" labels;
# sample IDs may themselves contain dots, so the suffix is matched explicitly
# rather than splitting on the last separator.
KNOWN_ASSEMBLERS = ("megahit", "metaspades")

# Column order of the output table. Kept explicit rather than derived from the
# reports so the schema stays stable even when QUAST's metric set changes.
COLUMNS = [
    "sample",
    "assembler",
    "n_contigs",
    "total_length",
    "largest_contig",
    "gc_percent",
    "n50",
    "n75",
    "n90",
    "l50",
    "l75",
    "l90",
    "aun",
    "n_contigs_ge_1kb",
    "n_contigs_ge_5kb",
    "n_contigs_ge_10kb",
    "n_contigs_ge_25kb",
    "n_contigs_ge_50kb",
    "total_length_ge_1kb",
    "total_length_ge_10kb",
    "total_length_ge_50kb",
    "ns_per_100kbp",
]

# QUAST metric name -> output column.
#
# Note the difference between the two families of metrics:
#   "# contigs"            counts contigs at or above --min-contig (default 500)
#   "# contigs (>= 0 bp)"  counts every contig the assembler emitted
# The unqualified names are the ones that describe the assembly as it will be
# used downstream, so those are what land in n_contigs / total_length.
METRIC_MAP = {
    "# contigs": "n_contigs",
    "Total length": "total_length",
    "Largest contig": "largest_contig",
    "GC (%)": "gc_percent",
    "N50": "n50",
    "N75": "n75",
    "N90": "n90",
    "L50": "l50",
    "L75": "l75",
    "L90": "l90",
    "auN": "aun",
    "# contigs (>= 1000 bp)": "n_contigs_ge_1kb",
    "# contigs (>= 5000 bp)": "n_contigs_ge_5kb",
    "# contigs (>= 10000 bp)": "n_contigs_ge_10kb",
    "# contigs (>= 25000 bp)": "n_contigs_ge_25kb",
    "# contigs (>= 50000 bp)": "n_contigs_ge_50kb",
    "Total length (>= 1000 bp)": "total_length_ge_1kb",
    "Total length (>= 10000 bp)": "total_length_ge_10kb",
    "Total length (>= 50000 bp)": "total_length_ge_50kb",
    "# N's per 100 kbp": "ns_per_100kbp",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarise QUAST report.tsv files into one TSV table"
    )
    parser.add_argument(
        "report_files",
        type=Path,
        nargs="+",
        help="One or more QUAST report.tsv files",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output TSV path (default: stdout)",
    )
    return parser.parse_args()


def split_label(label):
    """
    Split a '<sample>.<assembler>' assembly label into its two identities.

    Falls back to (label, "NA") for labels this pipeline did not generate, so a
    hand-run QUAST report still produces a usable row instead of a crash.
    """
    for assembler in KNOWN_ASSEMBLERS:
        suffix = f".{assembler}"
        if label.endswith(suffix):
            return label[: -len(suffix)], assembler
    return label, "NA"


def read_report(path):
    """
    Read one QUAST report.tsv into [(assembly_label, {metric: value}), ...].

    The file is metric-per-row, assembly-per-column. One column is the normal
    case here (one QUAST run per assembly), but several columns are handled too
    so a manually combined report still parses.
    """
    with open(path, newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))

    if not rows:
        sys.stderr.write(f"ERROR: QUAST report is empty: {path}\n")
        sys.exit(1)

    header = rows[0]
    if not header or header[0] != "Assembly":
        sys.stderr.write(
            f"ERROR: {path} does not look like a QUAST report.tsv "
            f"(first cell is {header[0]!r}, expected 'Assembly')\n"
        )
        sys.exit(1)

    labels = header[1:]
    metrics = [{} for _ in labels]

    for row in rows[1:]:
        if not row:
            continue
        name = row[0]
        for index in range(len(labels)):
            # Guard against ragged rows: QUAST leaves a cell empty when a metric
            # does not apply to one of several assemblies in a combined report.
            value = row[index + 1] if index + 1 < len(row) else ""
            metrics[index][name] = value

    return list(zip(labels, metrics))


def build_row(label, metrics):
    """Map one assembly's QUAST metrics onto the fixed output schema."""
    sample, assembler = split_label(label)

    row = {column: "NA" for column in COLUMNS}
    row["sample"] = sample
    row["assembler"] = assembler

    for quast_name, column in METRIC_MAP.items():
        value = metrics.get(quast_name, "")
        # QUAST writes "-" for metrics it computed but cannot express
        if value not in ("", "-"):
            row[column] = value

    return row


def main():
    args = parse_args()

    rows = []
    for report_path in args.report_files:
        if not report_path.exists():
            sys.stderr.write(f"ERROR: QUAST report not found: {report_path}\n")
            sys.exit(1)

        for label, metrics in read_report(report_path):
            rows.append(build_row(label, metrics))

    if not rows:
        sys.stderr.write("ERROR: no assemblies found in the given QUAST reports\n")
        sys.exit(1)

    missing_n90 = [r["sample"] for r in rows if r["n90"] == "NA"]
    if missing_n90:
        sys.stderr.write(
            "WARNING: N90 absent for: "
            + ", ".join(missing_n90)
            + " — QUAST reports it only with --report-all-metrics (QUAST >= 5.2)\n"
        )

    rows.sort(key=lambda row: (row["sample"], row["assembler"]))

    handle = open(args.output, "w", newline="") if args.output else sys.stdout
    try:
        writer = csv.DictWriter(
            handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    finally:
        if args.output:
            handle.close()


if __name__ == "__main__":
    main()
