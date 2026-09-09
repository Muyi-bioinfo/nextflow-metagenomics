#!/usr/bin/env python3
"""
Validate and parse samplesheet for nextflow-metagenomics pipeline.

Expected columns: sample,fastq_1,fastq_2,group,batch,host

fastq_2 may be left empty for single-end data; single_end=True is written
to the validated output in that case.

Checks:
- File existence (fastq_1 always; fastq_2 only when non-empty)
- Duplicate sample IDs
- Missing required metadata
- File extensions (.fastq.gz, .fq.gz)

Outputs validated samplesheet to stdout (adds single_end column).
"""

import sys
import csv
import argparse
from pathlib import Path
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate samplesheet for metagenomics pipeline"
    )
    parser.add_argument(
        "samplesheet",
        type=Path,
        help="Input samplesheet CSV file"
    )
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=None,
        help="Base directory for resolving relative paths (default: current directory)"
    )
    return parser.parse_args()


def check_file_exists(filepath, sample_id, read_type, base_dir=None):
    """Check if FASTQ file exists."""
    if not filepath:
        return f"Sample '{sample_id}': {read_type} file path is empty"

    path = Path(filepath)

    # If relative path and base_dir provided, resolve relative to base_dir
    if base_dir and not path.is_absolute():
        path = base_dir / path

    if not path.exists():
        return f"Sample '{sample_id}': {read_type} file not found: {filepath}"

    if not path.is_file():
        return f"Sample '{sample_id}': {read_type} is not a file: {filepath}"

    return None


def check_file_extension(filepath, sample_id, read_type):
    """Check if file has valid FASTQ extension."""
    valid_extensions = ('.fastq.gz', '.fq.gz', '.fastq', '.fq')
    path = Path(filepath)

    if not any(str(path).endswith(ext) for ext in valid_extensions):
        return (f"Sample '{sample_id}': {read_type} has invalid extension. "
                f"Expected: {', '.join(valid_extensions)}")

    return None


def validate_samplesheet(samplesheet_path, base_dir=None):
    """
    Validate samplesheet format and content.

    Returns:
        tuple: (is_valid, errors, rows)
    """
    errors = []
    rows = []
    sample_ids = set()

    # Check if samplesheet exists
    if not samplesheet_path.exists():
        return False, [f"Samplesheet not found: {samplesheet_path}"], []

    # Required columns (fastq_2 is present as a column but may be empty for SE)
    required_columns = {'sample', 'fastq_1', 'group', 'batch', 'host'}

    try:
        with open(samplesheet_path, 'r') as f:
            reader = csv.DictReader(f)

            # Check header
            if not reader.fieldnames:
                return False, ["Samplesheet is empty"], []

            missing_cols = required_columns - set(reader.fieldnames)
            if missing_cols:
                return False, [f"Missing required columns: {', '.join(missing_cols)}"], []

            # Validate each row
            for line_num, row in enumerate(reader, start=2):  # start=2 因为第1行是header
                sample_id = row.get('sample', '').strip()
                fastq_1 = row.get('fastq_1', '').strip()
                fastq_2 = row.get('fastq_2', '').strip()
                group = row.get('group', '').strip()
                batch = row.get('batch', '').strip()
                host = row.get('host', '').strip()

                # Check for empty sample ID
                if not sample_id:
                    errors.append(f"Line {line_num}: Sample ID is empty")
                    continue

                # Check for duplicate sample IDs
                if sample_id in sample_ids:
                    errors.append(f"Line {line_num}: Duplicate sample ID: {sample_id}")
                    continue

                sample_ids.add(sample_id)

                # Check required metadata
                if not group:
                    errors.append(f"Sample '{sample_id}': Missing 'group' metadata")
                if not batch:
                    errors.append(f"Sample '{sample_id}': Missing 'batch' metadata")
                if not host:
                    errors.append(f"Sample '{sample_id}': Missing 'host' metadata")

                # Check FASTQ files
                if not fastq_1:
                    errors.append(f"Sample '{sample_id}': Missing fastq_1 path")
                else:
                    # Check file existence
                    err = check_file_exists(fastq_1, sample_id, "fastq_1", base_dir)
                    if err:
                        errors.append(err)
                    else:
                        # Check extension
                        err = check_file_extension(fastq_1, sample_id, "fastq_1")
                        if err:
                            errors.append(err)

                # fastq_2 is optional: empty means single-end
                single_end = not bool(fastq_2)
                if fastq_2:
                    err = check_file_exists(fastq_2, sample_id, "fastq_2", base_dir)
                    if err:
                        errors.append(err)
                    else:
                        err = check_file_extension(fastq_2, sample_id, "fastq_2")
                        if err:
                            errors.append(err)

                # Store validated row
                rows.append({
                    'sample':      sample_id,
                    'fastq_1':     fastq_1,
                    'fastq_2':     fastq_2,
                    'single_end':  str(single_end),
                    'group':       group,
                    'batch':       batch,
                    'host':        host
                })

    except Exception as e:
        return False, [f"Error reading samplesheet: {str(e)}"], []

    # Final check
    if not rows:
        errors.append("No valid samples found in samplesheet")

    is_valid = len(errors) == 0
    return is_valid, errors, rows


def main():
    args = parse_args()

    # Validate samplesheet
    is_valid, errors, rows = validate_samplesheet(args.samplesheet, args.base_dir)

    if not is_valid:
        sys.stderr.write("ERROR: Samplesheet validation failed\n\n")
        for error in errors:
            sys.stderr.write(f"  - {error}\n")
        sys.stderr.write("\n")
        sys.exit(1)

    # Output validated samplesheet to stdout
    if rows:
        writer = csv.DictWriter(
            sys.stdout,
            fieldnames=['sample', 'fastq_1', 'fastq_2', 'single_end', 'group', 'batch', 'host']
        )
        writer.writeheader()
        writer.writerows(rows)

    sys.exit(0)


if __name__ == "__main__":
    main()
