#!/usr/bin/env python


"""
https://www.ncbi.nlm.nih.gov/tools/vecscreen/about/#Parameters
Vector contamination usually occurs at the beginning or end of a sequence; therefore,
different criteria are applied for terminal and internal matches
"""

import sys
from Bio import SeqIO


"""
$ blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true \
         -evalue 700 -searchsp 1750000000000 -outfmt 6 -word_size 11                \
         -db /shares/CIBIO-Storage/CM/mir/databases/univec/UniVec_Core              \
         -query input_assembly
"""


def already_matched(contig_match, hits_found, contig_id):
    cm_ss, cm_ee = contig_match

    if contig_id in hits_found:
        hf_ss, hf_ee = hits_found[contig_id]

        if (cm_ss >= hf_ss) and (cm_ee <= hf_ee):
            return True

    return False


fwrv = {'gnl|uv|NGB00727.1:1-51': 51,  # forward
        'gnl|uv|NGB00728.1:1-51': 51,  # forward
        'gnl|uv|NGB00729.1:1-51': 51,  # forward
        'gnl|uv|NGB00730.1:1-51': 51,  # forward
        'gnl|uv|NGB00734.1:1-51': 51,  # forward
        'gnl|uv|NGB00735.1:1-47': 47,  # reverse
        'gnl|uv|NGB00736.1:1-47': 47,  # reverse
        'gnl|uv|NGB00738.1:1-47': 47,  # reverse
        'gnl|uv|NGB00739.1:1-47': 47,  # reverse
        'gnl|uv|NGB00740.1:1-47': 47,  # reverse
        'gnl|uv|NGB00741.1:1-47': 47,  # reverse
        'gnl|uv|NGB00742.1:1-47': 47,  # reverse
        'gnl|uv|NGB00743.1:1-47': 47,  # reverse
        'gnl|uv|NGB00744.1:1-47': 47,  # reverse
        'gnl|uv|NGB00745.1:1-47': 47,  # reverse
        'gnl|uv|NGB00746.1:1-47': 47}  # reverse
univec_db = '/shares/CIBIO-Storage/CM/mir/databases/univec/UniVec_Core'
input_assembly = sys.argv[1]
offset = 25
contig_lens = dict([(r.id, len(r)) for r in SeqIO.parse(input_assembly, "fasta")])
adapter_lens = dict([(r.id, len(r)) for r in SeqIO.parse(univec_db, "fasta")])
output_b6o = sys.argv[2]
hits = [(l.strip().split('\t')[0],  # contig_id
         l.strip().split('\t')[1],  # adapter_id
         float(l.strip().split('\t')[2]),  # percentage_identity
         int(l.strip().split('\t')[3]),  # len_match
         int(l.strip().split('\t')[5]),  # gaps
         float(l.strip().split('\t')[10]),  # e_value
         (int(l.strip().split('\t')[6]), int(l.strip().split('\t')[7])),  # contig_start, contig_end
         (int(l.strip().split('\t')[8]), int(l.strip().split('\t')[9])))  # adapter_start, adapter_end
        for l in open(output_b6o) if l.strip().split('\t')[1] in fwrv]
hits_found = {}

for entry in hits:
    contig_id, adapter_id, pid, len_match, gaps, evalue, contig_match, adapter_match = entry

    if (pid >= 100) and ((len_match / adapter_lens[adapter_id]) >= 35):
        if not already_matched(contig_match, hits_found, contig_id):
            hits_found[contig_id] = contig_match
            print('M1', entry)
        elif len_match >= adapter_lens[adapter_id]:
            print('M1', entry, 'multiple')

    if (pid >= 95) and (len_match >= adapter_lens[adapter_id]) and (gaps <= 0):
        if not already_matched(contig_match, hits_found, contig_id):
            hits_found[contig_id] = contig_match
            print('M3', entry)
        elif len_match >= adapter_lens[adapter_id]:
            print('M3', entry, 'multiple')

    if (pid >= 90) and (evalue <= 0.01):
        if not already_matched(contig_match, hits_found, contig_id):
            hits_found[contig_id] = contig_match
            print('M4', entry)
        elif len_match >= adapter_lens[adapter_id]:
            print('M4', entry, 'multiple')

    if (pid >= 100) and (fwrv[adapter_id] in adapter_match) and (evalue < 10):
        if not already_matched(contig_match, hits_found, contig_id):
            hits_found[contig_id] = contig_match
            print('M2', entry)
        elif len_match >= adapter_lens[adapter_id]:
            print('M2', entry, 'multiple')
