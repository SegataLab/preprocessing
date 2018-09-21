#!/usr/bin/env python3


from Bio import SeqIO
from Bio import pairwise2


couples = [('N702', 'S506'),
           ('N702', 'S508'),
           ('N703', 'S508'),
           ('N705', 'S502'),
           ('N706', 'S502'),
           ('N707', 'S502'),
           ('N708', 'S502'),
           ('N706', 'S508'),
           ('N707', 'S508'),
           ('N708', 'S508'),
           ('N709', 'S508'),
           ('N710', 'S508'),
           ('N711', 'S508'),
           ('N705', 'S508'),
           ('N712', 'S508')]
seq2idx = {'CTAGTACG': 'N702', 'CGTACTAG': 'N702',
           'TTCTGCCT': 'N703', 'AGGCAGAA': 'N703',
           'AGGAGTCC': 'N705', 'GGACTCCT': 'N705',
           'CATGCCTA': 'N706', 'TAGGCATG': 'N706',
           'GTAGAGAG': 'N707', 'CTCTCTAC': 'N707',
           'CCTCTCTG': 'N708', 'CAGAGAGG': 'N708',
           'AGCGTAGC': 'N709', 'GCTACGCT': 'N709',
           'CAGCCTCG': 'N710', 'CGAGGCTG': 'N710',
           'TGCCTCTT': 'N711', 'AAGAGGCA': 'N711',
           'TCCTCTAC': 'N712', 'GTAGAGGA': 'N712',
           'CTCTCTAT': 'S502',
           'ACTGCATA': 'S506',
           'CTAAGCCT': 'S508'}
output = {}

for r in SeqIO.parse('Undetermined_S0_L007_R1_001.fastq', 'fastq'):
# for r in SeqIO.parse('Undetermined_S0_L007_R2_001.fastq', 'fastq'):
    n7xx = None
    s5xx = None
    s1, s2 = r.description.split(' ')[-1].split(':')[-1].split('+')
    # print(r.description, r.description.split(' ')[-1].split(':')[-1].split('+'))

    for seq, idx in seq2idx.items():
        m1 = [a[0] for a in pairwise2.align.globalms(seq, s1, 9, -7, -9, -9) if (a[-3] >= 40) and (a[-1] == 8)]
        m2 = [a[0] for a in pairwise2.align.globalms(seq, s2, 9, -7, -9, -9) if (a[-3] >= 40) and (a[-1] == 8)]

        if m1:
            if len(m1) > 1:
                print('"{}" more than one match for n7xx'.format(r.id))
                break

            if not n7xx:
                n7xx = m1[0]
            else:
                print('"{}" new match found for n7xx'.format(r.id))
                break

        if m2:
            if len(m2) > 1:
                print('"{}" more than one match for s5xx'.format(r.id))
                break

            if not s5xx:
                s5xx = m2[0]
            else:
                print('"{}" new match found for s5xx'.format(r.id))
                break

        if n7xx and s5xx:
            if (seq2idx[n7xx], seq2idx[s5xx]) in couples:
                idd = '{}_{}_R1'.format(seq2idx[n7xx], seq2idx[s5xx])

                if idd in output:
                    output[idd].append(r)
                else:
                    output[idd] = [r]

                break

for fn, seqr in output.items():
    print('writing {}.fastq'.format(fn))
    SeqIO.write(seqr, '{}.fastq'.format(fn), 'fastq')
