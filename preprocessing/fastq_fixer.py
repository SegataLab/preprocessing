#!/usr/bin/env python3


__author__ = 'Francesco Asnicar (f.asnicar@unitn.it)'
__version__ = '0.01'
__date__ = '7 June 2019'


import sys
import bz2
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import argparse as ap
import time


if sys.version_info[0] < 3:
    raise Exception("fastq_fixer requires Python 3, your current Python version is {}.{}.{}"
                    .format(sys.version_info[0], sys.version_info[1], sys.version_info[2]))


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
    p = ap.ArgumentParser(formatter_class=ap.ArgumentDefaultsHelpFormatter)

    p.add_argument('-i', '--input', type=str, required=True, help="")
    p.add_argument('-o', '--output', type=str, required=False, default=None, help="")
    p.add_argument('--verbose', action='store_true', default=False, help="Prints more stuff")
    p.add_argument('--overwrite', action='store_true', default=False, help="Overwrite output file if exists")
    p.add_argument('-v', '--version', action='version',
                   version='fastq_fixer.py version {} ({})'.format(__version__, __date__),
                   help="Prints the current fastq_fixer.py version and exit")
    return p.parse_args()


def check_params(args, verbose=False):
    if not os.path.isfile(args.input):
        error('input file "{}" does not exists'.format(args.input), exit=True)

    if os.path.isfile(args.output) and (not args.overwrite):
        output_prefix, ext = os.path.splitext(args.output)
        timestamp = str(datetime.datetime.today().strftime('%Y%m%d%H%M%S'))
        output = output_prefix + "_" + timestamp + ext

        if verbose:
            info('Output file "{}" exists, writing output to "{}"'.format(args.output, output), exit=True)

        args.output = output

    if verbose:
        info('Arguments: {}\n'.format(vars(args)))


def fastq_fixer():
    args = read_params()

    if args.verbose:
        info('fastq_fixer.py version {} ({})\n'.format(__version__, __date__))
        info('Command line: {}\n\n'.format(' '.join(sys.argv)), init_new_line=True)

    check_params(args, verbose=args.verbose)

    for r in SeqIO.parse(args.input):
        pass

    SeqIO.write(seqr, args.output, 'fastq')


if __name__ == '__main__':
    t0 = time.time()
    fastq_fixer()
    t1 = time.time()
    info('Total elapsed time {}s\n'.format(int(t1 - t0)), init_new_line=True)
    sys.exit(0)
