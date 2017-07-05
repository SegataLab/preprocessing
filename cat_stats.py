#!/usr/bin/env python


__author__ = 'Francesco Asnicar (f.asnicar@unitn.it)'
__version__ = '0.1'
__date__ = '5 July 2017'


import sys
import os
import glob
import argparse


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
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input_dir', required=True, default=None, type=str, help='Input folder containing the files to cat')
    p.add_argument('-e', '--extension', required=False, default=".stats", type=str, help='The extension of the files to cat, default ".stats"')
    p.add_argument('-o', '--output_file', required=True, default=None, type=str, help='The output file to write')
    p.add_argument('--inp_sep', required=False, default="\t", type=str, help='Define the field-separator character for the input files, default "\\t"')
    p.add_argument('--out_sep', required=False, default="\t", type=str, help='Define the field-separator character for the output file, default "\\t"')

    return p.parse_args()


def check_params(args):
    if not args.input_dir.endswith('/'):
        args.input_dir += '/'

    if not args.extension.startswith('.'):
        args.extension = '.'+args.extension

    if len(args.inp_sep) > 1:
        error('input field-separator is more than one characther, "{}"'.format(args.inp_sep), exit=True)

    if len(args.out_sep) > 1:
        error('output field-separator is more than one characther, "{}"'.format(args.out_sep), exit=True)


def cat(input_dir, ext, inp_sep):
    header = None
    data = []

    for stat_file in sorted(glob.glob(input_dir+'*'+ext)):
        with open(stat_file) as f:
            for row in f:
                if not header and row.startswith('#'):
                    header = row.strip().split(inp_sep)

                if not row.startswith('#'):
                    data.append(row.strip().split(inp_sep))

    if header:
        data = [header]+data

    return data


def write_output_file(output_file, out_sep, data):
    if os.path.isfile(output_file):
        error('output file "{}" already exists'.format(output_file), exit=True)

    with open(output_file, 'w') as f:
        f.write('\n'.join([out_sep.join(a) for a in data]) + '\n')


if __name__ == "__main__":
    args = read_params()
    check_params(args)
    data = cat(args.input_dir, args.extension, args.inp_sep)
    write_output_file(args.output_file, args.out_sep, data)
    sys.exit(0)
