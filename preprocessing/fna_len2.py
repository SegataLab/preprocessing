#!/usr/bin/env python

import sys
import argparse
from Bio import SeqIO
import numpy as np
import gzip
import bz2
import io

def openr(fn, mode="rt"):  # Force text mode
    if fn is None:
        return sys.stdin
    if fn.endswith(".bz2"):
        return bz2.open(fn, mode)
    elif fn.endswith(".gz"):
        return gzip.open(fn, mode)
    else:
        return open(fn, mode)

def openw(fn, mode="wt"): # Force text mode
    if fn is None:
        return sys.stdout
    if fn.endswith(".bz2"):
        return bz2.open(fn, mode)
    elif fn.endswith(".gz"):
        return gzip.open(fn, mode)
    else:
        return open(fn, mode)

def read_params():
    p = argparse.ArgumentParser(description='Display the length of each fasta entry')
    p.add_argument('inp_f', metavar='INPUT_FILE', nargs='?', default=None, type=str,
        help="the input fna file [stdin if not present]")
    p.add_argument('out_f', metavar='OUTPUT_FILE', nargs='?', default=None, type=str,
        help="the output txt file [stdout if not present]")
    p.add_argument('-q', action='store_true', help="set this for fastq")
    p.add_argument('-t', '--total', action='store_true', default=False,
        help="Print only the sum of the length of all sequences\n")
    p.add_argument('-s', '--stat', action='store_true', default=False,
        help="Print only the statistics about the length of sequences\n")
    p.add_argument('-g', '--gaps', action='store_true', default=False,
        help="Print in addition to the length of the entries also the number of gaps\n")
    return vars(p.parse_args())

if __name__ == '__main__':
    par = read_params()
    lvec = []

    with openr(par['inp_f']) as inf, openw(par['out_f']) as outf:
      if par['stat'] or par['total']:
        for record in SeqIO.parse(inf, "fastq" if par['q'] else "fasta"):
            lvec.append(len(record.seq))

        if not lvec:
            sys.stderr.write(f'[e] "{par["inp_f"]}" no values found!\n')
            sys.exit(0)

        if par['total']:
            outf.write(str(sum(lvec)) + "\n")

        if par['stat']:
            lvec = np.array(lvec)
            outf.write("#samplename\tn_of_bases\tn_of_reads\tmin_read_len\tmedian_read_len\tmean_read_len\tmax_read_len\n")
            outf.write(f"{par['inp_f']}\t{np.sum(lvec)}\t{lvec.size}\t{np.min(lvec)}\t{np.median(lvec)}\t{np.mean(lvec)}\t{np.max(lvec)}\n")
      else:
        if par['gaps']:
          for record in SeqIO.parse(inf, "fastq" if par['q'] else "fasta"):
            outf.write(f"{record.id}\t{len(record.seq)}\t{str(record.seq).count('-')}\n")
        else:
          for record in SeqIO.parse(inf, "fastq" if par['q'] else "fasta"):
            outf.write(f"{record.id}\t{len(record.seq)}\n")
