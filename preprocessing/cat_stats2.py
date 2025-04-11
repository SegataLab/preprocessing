#!/usr/bin/env python

__author__ = ('Francesco Asnicar (f.asnicar@unitn.it), '
              'Paolo Manghi (paolo.manghi@unitn.it)')
__version__ = '0.4'
__date__ = '28 January 2025'

import sys
import os
import argparse
import numpy as np
import logging

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

def read_params():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-i', '--input_files', required=True, nargs='+', type=str,
                   help='List of input files to concatenate')
    p.add_argument('-o', '--output_file', required=True, default=None, type=str,
                   help='The output file to write')
    p.add_argument('--inp_sep', required=False, default="\t", type=str,
                   help='Define the field-separator character for the input files')
    p.add_argument('--out_sep', required=False, default="\t", type=str,
                   help='Define the field-separator character for the output file')
    return p.parse_args()

def check_params(args):
    missing_files = [f for f in args.input_files if not os.path.isfile(f)]
    if missing_files:
        logger.error(f"The following input files are missing: {', '.join(missing_files)}")
        sys.exit(1)

    if len(args.inp_sep) > 1:
        logger.error(f"Input field-separator must be a single character: '{args.inp_sep}'")
        sys.exit(1)

    if len(args.out_sep) > 1:
        logger.error(f"Output field-separator must be a single character: '{args.out_sep}'")
        sys.exit(1)

def cat(input_files, inp_sep, output_file):
    header = None
    data = []

    for stat_file in input_files:
        if os.path.basename(stat_file) == os.path.basename(output_file):
            logger.warning(f"Skipping output file '{output_file}' from input list.")
            continue

        logger.info(f"Processing file: {stat_file}")
        with open(stat_file) as f:
            for row in f:
                if not header and row.startswith('#'):
                    header = row.strip().split(inp_sep)

                if not row.startswith('#'):
                    data.append([os.path.basename(row.strip().split(inp_sep)[0])] + row.strip().split(inp_sep)[1:])

    if header:
        data = [header] + data

    if not data:
        logger.error("No valid data found in input files.")
        sys.exit(1)

    return data

def line_of_interest(sample_id):
    return (sample_id.endswith('_R1.fastq.bz2') or sample_id.endswith('_R2.fastq.bz2') or
            sample_id.endswith('_UN.fastq.bz2') or sample_id.endswith('_UP.fastq.bz2'))

def line_to_str_or_float(line):
    return [float(f) if f.replace('.', '', 1).isdigit() else str(f) for f in line]

def summarize(data):
    funcs = {'n_of_bases': sum,
             'n_of_reads': sum,
             'min_read_len': min,
             'median_read_len': np.median,
             'mean_read_len': np.mean,
             'max_read_len': max}
    convs = {'n_of_bases': int,
             'n_of_reads': int,
             'min_read_len': int,
             'median_read_len': float,
             'mean_read_len': float,
             'max_read_len': int}

    logger.info("Summarizing data.")
    interesting_lines = [line_to_str_or_float(line) for line in data if line_of_interest(line[0])]

    if not interesting_lines:
        logger.error("No lines of interest found in the input data.")
        sys.exit(1)

    summary = []
    for pos, fld in enumerate(data[0]):
        if pos == 0:
            summary.append('_'.join(interesting_lines[0][0].split('_')[:-1]))
        else:
            try:
                summary.append(convs[fld](funcs[fld]([field[pos] for field in interesting_lines])))
            except KeyError:
                logger.warning(f"Field '{fld}' not recognized for summarization.")
                summary.append('NA')

    return data + [summary]

def write_output_file(output_file, out_sep, data):
    if os.path.isfile(output_file):
        logger.error(f"Output file '{output_file}' already exists.")
        sys.exit(1)

    logger.info(f"Writing output to file: {output_file}")
    with open(output_file, 'w') as f:
        f.write('\n'.join([out_sep.join(map(str, a)) for a in data]) + '\n')

if __name__ == "__main__":
    args = read_params()
    check_params(args)

    logger.info("Starting file concatenation and summarization.")
    data = summarize(cat(args.input_files, args.inp_sep, args.output_file))
    write_output_file(args.output_file, args.out_sep, data)

    logger.info("Process completed successfully.")
    sys.exit(0)
