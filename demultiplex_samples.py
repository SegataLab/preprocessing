#!/usr/bin/env python


__author__ = 'Francesco Asnicar (f.asnicar@unitn.it)'
__version__ = '0.3'
__date__ = '4 April 2018'


import sys
import os
import glob
import argparse
import shutil


if sys.version_info[0] < 3:
    raise Exception("Not running Python3")


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
    p.add_argument('-i', '--input_folder', required=True, type=str,
                   help="Input folder containing the samples' folders to rename")
    p.add_argument('-o', '--output_folder', required=True, type=str,
                   help="Output folder where move the renamed folders/files")
    p.add_argument('-m', '--mapping_file', required=True, type=str,
                   help='The mapping file containing the sample name and the forward and '
                        'reverse indexes associated')
    p.add_argument('-s', '--separator', default='\t', type=str,
                   help='The column separator in the mapping file')
    p.add_argument('-1', '--index_primer_1', default=2, type=int,
                   help='The column of the Index primer 1')
    p.add_argument('-2', '--index_primer_2', default=1, type=int,
                   help='The column of the Index primer 2')
    p.add_argument('-n', '--sample_name', default=0, type=int,
                   help='The column of the Sample IDs')

    return p.parse_args()


def check_params(args):
    if not os.path.isdir(args.input_folder):
        error('', exit=True)

    if not os.path.isfile(args.mapping_file):
        error('', exit=True)

    if len(args.separator) > 1:
        error('', exit=True)


def read_mapping_file(map_file, sep, sn, i1, i2):
    mapp = {}
    sample_names = []

    with open(map_file) as f:
        for row in f:
            if row.startswith('#'):
                continue

            row_clean = row.strip().split(sep)

            if (row_clean[i1], row_clean[i2]) not in mapp:
                mapp[(row_clean[i1], row_clean[i2])] = row_clean[sn].replace(' ', '_')
                sample_names.append(row_clean[sn].replace(' ', '_'))
            else:
                error('indexes "({},{})" already present, duplicated?'
                      .format(row_clean[i1], row_clean[i2]), exit=True)

    if len(sample_names) != len(set(sample_names)):
        error('There are duplicated sample names, check your mapping file', exit=True)

    return mapp


def demultiplex(input_folder, output_folder, mapping):
    for fld in glob.iglob(os.path.join(input_folder, '*')):
        new_name = None

        for k, v in mapping.items():
            if ('_'.join(k) in fld) or ('_'.join(k[::-1]) in fld):
                new_name = v
                break

        if not new_name:
            error('skipping "{}", no indexes found in mapping'.format(fld))
            continue

        output = os.path.join(output_folder, new_name)

        if os.path.isdir(output):
            error('skipping "{}", folder already exists'.format(output))
            continue
        else:
            info('Creating folder "{}"\n'.format(output))
            os.mkdir(output, mode=0o775)

        for src in (glob.glob(os.path.join(fld, '*_R1.fastq.bz2')) +
                    glob.glob(os.path.join(fld, '*_R2.fastq.bz2')) +
                    glob.glob(os.path.join(fld, '*_UN.fastq.bz2')) +
                    glob.glob(os.path.join(fld, '*_summary.stats'))):
            suffix = src.split('_')[-1]
            dest = os.path.join(output, new_name + '_' + suffix)

            info('Copying "{}" into "{}"\n'.format(src, dest))
            shutil.copy2(src, dest)


if __name__ == "__main__":
    args = read_params()
    check_params(args)
    mapp = read_mapping_file(args.mapping_file, args.separator, args.sample_name,
                             args.index_primer_1, args.index_primer_2)
    demultiplex(args.input_folder, args.output_folder, mapp)
    sys.exit(0)
