#!/usr/bin/env python


__author__ = ('Francesco Asnicar (f.asnicar@unitn.it), '
              'Paolo Manghi (paolo.manghi@unitn.it)')
__version__ = '0.3'
__date__ = '2 May 2018'


import sys
import os
import glob
import argparse
import numpy as np


def info(s, init_new_line=False, exit=False, exit_value=0):
    if init_new_line:
        sys.stdout.write('\n')

    sys.stdout.write('{}'.format(s))
    sys.stdout.flush()

    if exit:
        sys.exit(exit_value)


def error(s, init_new_line=False, exit=False, exit_value=1):
    if init_new_line:
        sys.stderr.write('\n')

    sys.stderr.write('[e] {}\n'.format(s))
    sys.stderr.flush()

    if exit:
        sys.exit(exit_value)


def read_params():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-i', '--input_dir', required=True, default=None, type=str,
                   help='Input folder containing the files to cat')
    p.add_argument('-e', '--extension', required=False, default=".stats", type=str,
                   help='The extension of the files to cat')
    p.add_argument('-o', '--output_file', required=True, default=None, type=str,
                   help='The output file to write')
    p.add_argument('--inp_sep', required=False, default="\t", type=str,
                   help='Define the field-separator character for the input files')
    p.add_argument('--out_sep', required=False, default="\t", type=str,
                   help='Define the field-separator character for the output file')

    return p.parse_args()


def check_params(args):
    if not args.input_dir.endswith('/'):
        args.input_dir += '/'

    if not args.extension.startswith('.'):
        args.extension = '.' + args.extension

    if len(args.inp_sep) > 1:
        error('input field-separator is more than one characther, "{}"'.format(args.inp_sep), exit=True)

    if len(args.out_sep) > 1:
        error('output field-separator is more than one characther, "{}"'.format(args.out_sep), exit=True)


def cat(input_dir, ext, inp_sep, output_file):
    header = None
    data = []

    for stat_file in sorted(glob.glob(input_dir + '*' + ext)):
        if os.path.basename(stat_file) == os.path.basename(output_file):  # should check not to include the output file
            continue

        with open(stat_file) as f:
            for row in f:
                if not header and row.startswith('#'):
                    header = row.strip().split(inp_sep)

                if not row.startswith('#'):
                    data.append(row.strip().split(inp_sep))

    if header:
        data = [header] + data

    return data


def line_of_interest(sample_id):
    return (sample_id.endswith('_R1.fastq') or sample_id.endswith('_R2.fastq') or
            sample_id.endswith('_UN.fastq') or sample_id.endswith('_UP.fastq'))


def line_to_str_or_float(line):
    return [float(f) if str(f)[-1].isdigit() else str(f) for f in line]


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

    interesting_lines = [line_to_str_or_float(line) for line in data if line_of_interest(line[0])]

    for pos, fld in enumerate(data[0]):
        if pos == 0:
            summary = ['_'.join(interesting_lines[0][0].split('_')[:-1])]
        else:
            summary += [convs[fld](funcs[fld]([field[pos] for field in interesting_lines]))]

    return data + [summary]


def write_output_file(output_file, out_sep, data):
    if os.path.isfile(output_file):
        error('output file "{}" already exists'.format(output_file), exit=True)

    with open(output_file, 'w') as f:
        f.write('\n'.join([out_sep.join(map(str, a)) for a in data]) + '\n')


if __name__ == "__main__":
    args = read_params()
    check_params(args)
    data = summarize(cat(args.input_dir, args.extension, args.inp_sep, args.output_file))
    write_output_file(args.output_file, args.out_sep, data)
    sys.exit(0)
