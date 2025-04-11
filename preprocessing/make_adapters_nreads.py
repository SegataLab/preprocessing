#!/usr/bin/env python


import os
import pandas as pd
import numpy as np
import argparse
import sys


parser = argparse.ArgumentParser()
parser.add_argument("-d", '--folder', required=True, help="Run folder from which load the FastQC results")
args = parser.parse_args()

folder = args.folder

if not os.path.isdir(folder):
    print(f"ERROR: Run folder \"{folder}\" does not exists!")
    sys.exit(1)

samples_to_skip = set(['.ipynb_checkpoints', '.empty', 'Stats', 'Reports'])

## ADAPTERS
adapters = {}

for sample in os.listdir(folder):
    folder_sample = os.path.join(folder, sample)

    if (sample in samples_to_skip) or (not os.path.isdir(folder_sample)):
        continue

    adapters[sample] = 0.0
    n = 0

    for lane in os.listdir(os.path.join(folder_sample, 'fastqc_raw')):
        folder_sample_lane = os.path.join(folder_sample, 'fastqc_raw', lane)

        if not os.path.isdir(folder_sample_lane):
            continue

        save = False
        tmp_list = []

        with open(os.path.join(folder_sample_lane, 'fastqc_data.txt')) as f:
            for r in f:
                if r.startswith('>>Adapter Content'):
                    save = True

                if save and r.startswith('>>END_MODULE'):
                    save = False
                    break

                if save:
                    tmp_list.append(r.strip().split('\t'))

        adapters[sample] += sum(map(float, tmp_list[-1][1:]))
        n += 1

    adapters[sample] /= n

## NREADS
nreads = {}

for sample in os.listdir(folder):
    folder_sample = os.path.join(folder, sample)

    if (sample in samples_to_skip) or (not os.path.isdir(folder_sample)):
        continue

    nreads[sample] = 0

    for lane in os.listdir(os.path.join(folder_sample, 'fastqc_raw')):
        folder_sample_lane = os.path.join(folder_sample, 'fastqc_raw', lane)
        
        if not os.path.isdir(folder_sample_lane):
            continue
        
        with open(os.path.join(folder_sample_lane, 'fastqc_data.txt')) as f:
            for r in f:
                if r.startswith('Total Sequences'):
                    nreads[sample] += int(r.strip().split('\t')[1])
                    break

## WRITING OUTPUT
if not len(set(adapters.keys()) - set(nreads.keys())):
    print('writing {}_adapters_nreads.tsv'.format(os.path.dirname(folder)))
    
    with open('{}_adapters_nreads.tsv'.format(os.path.dirname(folder)), 'w') as f:
        f.write('#sample\tadapters%\tnreads\n')
        f.write('\n'.join(['{}\t{}\t{}'.format(sample, adapters[sample], nreads[sample]) for sample in nreads.keys()]) + '\n')
else:
    print('ERROR: adapters and nreads keys not identical!')
    sys.exit(1)

sys.exit()
