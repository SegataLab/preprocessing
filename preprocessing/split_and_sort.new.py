#!/usr/bin/env python


from Bio import SeqIO
import argparse
import bz2


parser = argparse.ArgumentParser()
parser.add_argument("--R1", help="R1 in", required=True)
parser.add_argument("--R2", help="R2 in", required=True)
parser.add_argument("-p", '--prefix', help="Destination prefix", default='out')
parser.add_argument("-u", '--unpaired', help="File with unpaired read IDs")
args = parser.parse_args()

filR1 = bz2.open(args.prefix + '_R1.fastq.bz2', 'wt')
filR2 = bz2.open(args.prefix + '_R2.fastq.bz2', 'wt')
filUN = bz2.open(args.prefix + '_UN.fastq.bz2', 'wt')

R1_sortedIndex = SeqIO.index(args.R1, "fastq")
R2_sortedIndex = SeqIO.index(args.R2, "fastq")
unpaired_reads = set([i.strip() for i in bz2.open(args.unpaired, 'rt')])

for r1 in iter(sorted(R1_sortedIndex)):
    if r1 in R2_sortedIndex:
        SeqIO.write([R1_sortedIndex[r1]], filR1, 'fastq')
    elif r1 in unpaired_reads:
        SeqIO.write([R1_sortedIndex[r1]], filUN, 'fastq')

for r2 in sorted(R2_sortedIndex):
    if r2 in R1_sortedIndex:
        SeqIO.write([R2_sortedIndex[r2]], filR2, 'fastq')
    elif r2 in unpaired_reads:
        SeqIO.write([R2_sortedIndex[r2]], filUN, 'fastq')

filR1.close()
filR2.close()
filUN.close()
