#!/usr/bin/env python


from Bio import SeqIO
import argparse
import bz2


parser = argparse.ArgumentParser()

parser.add_argument("--R1", help="R1 in", required=True)
parser.add_argument("--R2", help="R2 in", required=True)
parser.add_argument("-p",'--prefix', help="Destination prefix", default='out')

args = parser.parse_args()

filR1 = bz2.BZ2File(args.prefix+'.R1.fastq.bz2', 'w')
filR2 = bz2.BZ2File(args.prefix+'.R2.fastq.bz2', 'w')
filUP = bz2.BZ2File(args.prefix+'.UP.fastq.bz2', 'w')

R1_sortedIndex = SeqIO.index(args.R1, "fastq")
R2_sortedIndex = SeqIO.index(args.R2, "fastq")

for r1 in iter(sorted(R1_sortedIndex)):
	if r1 not in R2_sortedIndex:
		SeqIO.write([R1_sortedIndex[r1]], filUP, 'fastq')		
	else: 
		SeqIO.write([R1_sortedIndex[r1]], filR1, 'fastq')

for r2 in sorted(R2_sortedIndex):
	if r2 not in R1_sortedIndex:
		SeqIO.write([R2_sortedIndex[r2]], filUP, 'fastq')
	else: 
		SeqIO.write([R2_sortedIndex[r2]], filR2, 'fastq')

filR1.close()
filR2.close()
filUP.close()
