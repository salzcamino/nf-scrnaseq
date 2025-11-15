#!/usr/bin/env python3

"""
Validate and convert samplesheet CSV to the format expected by the pipeline.

Expected columns:
- sample: Sample ID
- fastq_1: Path to read 1 FASTQ file
- fastq_2: Path to read 2 FASTQ file (optional for single-end)
- protocol: scRNA-seq protocol (optional, defaults to pipeline parameter)
"""

import os
import sys
import csv
import argparse
from pathlib import Path


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and reformat samplesheet for nf-scrnaseq pipeline"
    )
    parser.add_argument(
        "input",
        metavar="INPUT",
        help="Input samplesheet CSV file"
    )
    parser.add_argument(
        "output",
        metavar="OUTPUT",
        help="Output validated samplesheet CSV file"
    )
    return parser.parse_args()


def check_samplesheet(file_in, file_out):
    """
    Validate samplesheet and write out validated version.

    Args:
        file_in: Input samplesheet path
        file_out: Output validated samplesheet path
    """
    REQUIRED_COLUMNS = ['sample', 'fastq_1']
    OPTIONAL_COLUMNS = ['fastq_2', 'protocol']

    # Check file exists
    if not Path(file_in).is_file():
        print(f"ERROR: Samplesheet file does not exist: {file_in}", file=sys.stderr)
        sys.exit(1)

    # Read and validate samplesheet
    samples = {}
    with open(file_in, 'r') as fin:
        reader = csv.DictReader(fin)

        # Check required columns
        if not reader.fieldnames:
            print("ERROR: Samplesheet is empty", file=sys.stderr)
            sys.exit(1)

        for col in REQUIRED_COLUMNS:
            if col not in reader.fieldnames:
                print(f"ERROR: Missing required column: {col}", file=sys.stderr)
                sys.exit(1)

        # Process rows
        for line_num, row in enumerate(reader, start=2):
            sample = row['sample'].strip()
            fastq_1 = row['fastq_1'].strip()
            fastq_2 = row.get('fastq_2', '').strip()
            protocol = row.get('protocol', '').strip()

            # Validate sample ID
            if not sample:
                print(f"ERROR: Line {line_num}: Sample ID cannot be empty", file=sys.stderr)
                sys.exit(1)

            # Check for duplicate samples
            if sample in samples:
                print(f"ERROR: Line {line_num}: Duplicate sample ID: {sample}", file=sys.stderr)
                sys.exit(1)

            # Validate FASTQ files
            if not fastq_1:
                print(f"ERROR: Line {line_num}: fastq_1 cannot be empty for sample {sample}", file=sys.stderr)
                sys.exit(1)

            # Check file extensions
            valid_extensions = ['.fastq.gz', '.fq.gz', '.fastq', '.fq']
            if not any(fastq_1.endswith(ext) for ext in valid_extensions):
                print(f"WARNING: Line {line_num}: fastq_1 does not have a standard FASTQ extension: {fastq_1}", file=sys.stderr)

            if fastq_2 and not any(fastq_2.endswith(ext) for ext in valid_extensions):
                print(f"WARNING: Line {line_num}: fastq_2 does not have a standard FASTQ extension: {fastq_2}", file=sys.stderr)

            # Validate protocol if provided
            valid_protocols = ['10x_3prime', '10x_5prime', 'dropseq', 'smartseq2', 'smartseq3', 'indrop']
            if protocol and protocol not in valid_protocols:
                print(f"WARNING: Line {line_num}: Unknown protocol '{protocol}'. Valid options: {', '.join(valid_protocols)}", file=sys.stderr)

            # Store validated sample
            samples[sample] = {
                'sample': sample,
                'fastq_1': fastq_1,
                'fastq_2': fastq_2,
                'protocol': protocol
            }

    # Check we have at least one sample
    if not samples:
        print("ERROR: No valid samples found in samplesheet", file=sys.stderr)
        sys.exit(1)

    # Write validated samplesheet
    with open(file_out, 'w', newline='') as fout:
        fieldnames = ['sample', 'fastq_1', 'fastq_2', 'protocol']
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()
        for sample_data in samples.values():
            writer.writerow(sample_data)

    print(f"Successfully validated {len(samples)} samples", file=sys.stderr)


def main():
    """Main function."""
    args = parse_args()
    check_samplesheet(args.input, args.output)


if __name__ == '__main__':
    main()
