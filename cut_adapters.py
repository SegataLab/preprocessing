#!/usr/bin/env python


"""
https://www.ncbi.nlm.nih.gov/tools/vecscreen/about/#Parameters
Vector contamination usually occurs at the beginning or end of a sequence; therefore,
different criteria are applied for terminal and internal matches
"""

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

inp1 = sys.argv[1]  # input folder containing NCBI errors
inp2 = sys.argv[2]  # input folder containing the assemblies
outf = sys.argv[3]  # output folder where saving the cleaned assemblies

for f in os.listdir(inp1):
    print(f)

    exclude = []
    trim = {}

    with open(os.path.join(inp1, f)) as g:
        for r in g:
            if "PROKKA_contig" in r:
                if len(r.strip().split('\t')) == 3:
                    exclude.append(r.strip().split('\t')[0])
                elif len(r.strip().split('\t')) == 4:
                    trim[r.strip().split('\t')[0]] = [(int(ii.split('..')[0]), int(ii.split('..')[-1])) for ii in r.strip().split('\t')[-2].split(',')]

    # print('exclude')
    # print('\n'.join(exclude))
    # print('\ntrim')
    # print('\n'.join([str(ll) + ': ' + str(trim[ll]) for ll in trim]))
    # print()

    out_clean = []

    for r in SeqIO.parse(os.path.join(inp2, f.split('.')[0] + '.fasta'), "fasta"):
        if [r.id for cex in exclude if cex in r.id]:
            # print('{}\texclude'.format(r.id))
            continue

        if [r.id for ctr in trim if ctr in r.id]:
            # print('{}\ttrim'.format(r.id))
            multiple_trim = {'before': [], 'after': []}
            seq = r.seq
            offset = 0

            for s, e in trim[r.id.split('|')[-1]]:
                before = seq[:s - offset]
                after = seq[e - offset + 1:]

                if multiple_trim['after']:  # if there are, remove the last one added
                    multiple_trim['after'] = multiple_trim['after'][:-1]

                if len(before) >= 1000:
                    multiple_trim['before'].append(SeqRecord(before, id="{}___before[:{}]".format(r.id, s), description=""))

                if len(after) >= 1000:
                    multiple_trim['after'].append(SeqRecord(after, id="{}___after[{}+1:]".format(r.id, e), description=""))

                offset += e
                seq = seq[e - offset + 1:]

            out_clean += multiple_trim['before']
            out_clean += multiple_trim['after']
            continue

        # print(r.id)
        out_clean.append(r)

    # print()
    # print('\n'.join([r.id for r in out_clean]))
    # print()

    with open(os.path.join(outf, f.split('.')[0] + '.fasta'), 'w') as f:
        SeqIO.write(out_clean, f, "fasta")
