#!/usr/bin/env python


from __future__ import with_statement
import sys
import argparse
from Bio import SeqIO
import numpy as np
import bz2
import gzip


def openr(fn, mode="r"):
    if fn is None:
        return sys.stdin

    if fn.endswith(".bz2"):
        if sys.version_info[0] < 3:
            return bz2.BZ2File(fn)
        else:
            return bz2.open(fn, 'rt')
    elif fn.endswith(".gz"):
        if sys.version_info[0] < 3:
            return None  # need to check if gzip is different in Python2
        else:
            return gzip.open(fn, 'rt')
    else:
        return open(fn, mode)


def openw(fn):
    if fn is None:
        return sys.stdout

    if fn.endswith(".bz2"):
        if sys.version_info[0] < 3:
            return bz2.BZ2File(fn, "w")
        else:
            return bz2.open(fn, 'wt')
    elif fn.endswith(".gz"):
        if sys.version_info[0] < 3:
            return None  # need to check if gzip is different in Python2
        else:
            return gzip.open(fn, 'wt')
    else:
        return open(fn, "w")


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
    samplename = None

    if par['inp_f']:
        samplename = par['inp_f']

        if '/' in samplename:
            samplename = samplename[samplename.rfind('/') + 1:]

        if '.' in samplename:
            samplename = samplename[:samplename.rfind('.')]
    else:
        samplename = 'stdin'
        samplename += '_fastq' if par['q'] else '_fasta'

    with openw(par['out_f']) as outf:
        for r in SeqIO.parse(openr(par['inp_f']), "fastq" if par['q'] else "fasta"):
            lenn = len(r.seq)

            if par['stat'] or par['total']:
                lvec.append(lenn)
            else:
                if par['gaps']:
                    outf.write("\t".join([r.id, str(lenn), str(str(r.seq).count('-'))]) + "\n")
                else:
                    outf.write("\t".join([r.id, str(lenn)]) + "\n")

        if par['stat'] or par['total']:
            if not lvec:
                sys.stderr.write('[e] "{}" no values found!\n'.format(samplename))
                sys.exit(0)

        if par['total']:
            outf.write(str(np.sum(lvec)) + "\n")

        if par['stat']:
            outf.write('\t'.join(["#samplename", "n_of_bases", "n_of_reads", "min_read_len", "median_read_len", "mean_read_len", "max_read_len"]) + '\n')
            outf.write('\t'.join([str(a) for a in [samplename, np.sum(lvec), len(lvec), np.min(lvec), np.median(lvec), np.mean(lvec), np.max(lvec)]]) + '\n')
