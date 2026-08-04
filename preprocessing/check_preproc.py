#!/usr/bin/env python3

"""
Check that preprocessing output files exist (and are non-empty) for every
sample subfolder in a run folder.

Expected output files per sample <SAMPLE>/ subfolder:
    <SAMPLE>_R1.fastq.bz2
    <SAMPLE>_R2.fastq.bz2
    <SAMPLE>_UN.fastq.bz2
    <SAMPLE>_summary.stats

A sample is marked as FAILED if any of these files is missing, or if a
.bz2 file has size <= 14 bytes (the size of an empty bz2 archive), or if
the .stats file is 0 bytes.

Usage:
    python check_preproc_outputs.py -d /path/to/run_folder
"""

import argparse
from pathlib import Path

EMPTY_BZ2_SIZE = 14  # bytes: size of a bz2 archive with no content


def required_files(sample_dir: Path):
    sample = sample_dir.name
    return {
        f"{sample}_R1.fastq.bz2": EMPTY_BZ2_SIZE,
        f"{sample}_R2.fastq.bz2": EMPTY_BZ2_SIZE,
        f"{sample}_UN.fastq.bz2": EMPTY_BZ2_SIZE,
        f"{sample}_summary.stats": 0,
    }


def check_sample(sample_dir: Path):
    """Return True if sample passed, False otherwise."""
    for fname, min_size in required_files(sample_dir).items():
        fpath = sample_dir / fname
        if not fpath.is_file() or fpath.stat().st_size <= min_size:
            return False
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Check preprocessing outputs for each sample in a run folder."
    )
    parser.add_argument("-d", "--folder", required=True, help="Run folder")
    args = parser.parse_args()

    run_folder = Path(args.folder).resolve()
    sample_dirs = sorted(p for p in run_folder.iterdir() if p.is_dir())

    failed = [d.name for d in sample_dirs if not check_sample(d)]

    print(f"Checked {len(sample_dirs)} samples in {run_folder}")
    print(f"Passed: {len(sample_dirs) - len(failed)}  Failed: {len(failed)}")

    if failed:
        out_path = run_folder.parent / f"{run_folder.name}_samples_failed_preproc.txt"
        out_path.write_text("\n".join(failed) + "\n")
        print(f"Failed samples written to: {out_path}")
    else:
        print("All samples passed.")


if __name__ == "__main__":
    main()
