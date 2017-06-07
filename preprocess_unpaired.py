#!/usr/bin/env python
#Author: Duy Tin Truong (duytin.truong@unitn.it)
#        at CIBIO, University of Trento, Italy


__author__ = 'Duy Tin Truong (duytin.truong@unitn.it), Francesco Asnicar (f.asnicar@unitn.it)'
__version__ = '0.4'
__date__ = '7 June 2017'


import sys
import os
import glob
import argparse
from doit_loader import DoitLoader


def read_params():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input_dir', required=True, default=None, type=str,
        help='The input directory, path must be absolute!')
    p.add_argument('--extension', required=False, default=".fastq.gz", type=str,
        help='The extension of the raw input files. Default ".fastq.gz"')
    p.add_argument('--use_threads', required=False, action='store_true',
        help='Use threads')
    p.add_argument('--nprocs_main', required=False, default=1, type=int,
        help='Number of tasks to run in parallel')
    p.add_argument('--nprocs_bowtie2', required=False, default=1, type=int,
        help='Number of bowtie2 processors')
    p.add_argument('--remove_ribosomes', required=False, default=False, action='store_true',
        help='Remove ribosomes for mRNA datasets')
    p.add_argument('--clean', required=False, default=False, action='store_true',
        help='clean')
    p.add_argument('--keep_intermediate', required=False, default=False, action='store_true',
        help="If specified the script won't remove intermediate files")
    p.add_argument('--bowtie2_indexes', required=False, default='/CM/databases/bowtie2_indexes', type=str,
        help='Folder containing the bowtie2 indexes of the genomes to be removed from the samples. Default "/CM/databases/bowtie2_indexes"')

    args = p.parse_args()

    if not args.input_dir.startswith('/'):
        print 'ERROR, input directory is not an absolute path: "', args.input_dir, '"'
        sys.exit(1)

    if not args.input_dir.endswith('/'):
        args.input_dir += '/'

    return args


def get_inputs(input_dir, ext):
    inputs = dict()
    R1 = sorted(glob.glob('{}*R1*{}'.format(input_dir, ext)))
    R2 = sorted(glob.glob('{}*R2*{}'.format(input_dir, ext)))

    if not (len(R1) + len(R2)):
        for folder in os.listdir(input_dir):
            if not folder.endswith('/'):
                folder += '/'

            R1 = sorted(glob.glob('{}{}*R1*{}'.format(input_dir, folder, ext)))
            R2 = sorted(glob.glob('{}{}*R2*{}'.format(input_dir, folder, ext)))

            if R1 and R2:
                inputs[folder] = ([a[len(input_dir)+len(folder):] for a in R1], [a[len(input_dir)+len(folder):] for a in R2])
    else:
        folder = input_dir.split('/')[-2]
        input_dir = '/'.join(input_dir.split('/')[:-2])

        if not folder.endswith('/'):
            folder += '/'

        if not input_dir.endswith('/'):
            input_dir += '/'

        inputs[folder] = ([a[len(input_dir)+len(folder):] for a in R1], [a[len(input_dir)+len(folder):] for a in R2])

    if input_dir and inputs:
        return input_dir, inputs
    else:
        return None, None


def concatenate_reads(input_dir, inputs):
    merged = dict()

    for folder, (R1s, R2s) in inputs.iteritems():
        out_prefix = set(['_'.join(a.split('_')[:3]) for a in R1s + R2s])
        merged[folder] = list()

        if len(out_prefix) == 1:
            out_prefix = list(out_prefix)[0]
        else:
            print "\nERROR, concatenate_reads() cannot finds common filename!"
            print 'out_prefix', out_prefix
            print 'folder', folder
            print 'R1s', R1s
            print 'R2s', R2s
            print "\n"
            sys.exit(1)

        for Rs, oR in zip([R1s, R2s], [out_prefix+'.R1', out_prefix+'.R2']):
            cmd = 'zcat {}'.format(' '.join([input_dir+folder+a for a in Rs]))
            DoitLoader.add_task([input_dir+folder+oR+'.fastq'], [input_dir+folder+a for a in Rs], [cmd, {'stdout': input_dir+folder+oR+'.fastq'}])

            cmd = 'fna_len.py {} {} -q --stat'.format(input_dir+folder+oR+'.fastq', input_dir+folder+oR+'.stats')
            DoitLoader.add_task([input_dir+folder+oR+'.stats'], [input_dir+folder+oR+'.fastq'], [cmd])

            merged[folder].append(oR+'.fastq')

    return merged


def quality_control(input_dir, merged, keep_intermediate):
    qc = dict()
    files_to_remove = []

    for folder, Rs in merged.iteritems():
        qc[folder] = list()

        for R in Rs:
            oR = R[:R.rfind('.')]

            cmd = 'trim_galore --nextera --stringency 5 --length 75 --quality 20 --max_n 2 --trim-n --dont_gzip --no_report_file --suppress_warn --output_dir {} {}'.format(input_dir+folder, input_dir+folder+R)
            DoitLoader.add_task([input_dir+folder+oR+'_trimmed.fq'], [input_dir+folder+R], [cmd])

            if not keep_intermediate:
                files_to_remove.append(((input_dir+folder+oR+'_trimmed.fq', ), input_dir+folder+R))

            cmd = 'fna_len.py {} {} -q --stat'.format(input_dir+folder+oR+'_trimmed.fq', input_dir+folder+oR+'_trimmed.stats')
            DoitLoader.add_task([input_dir+folder+oR+'_trimmed.stats'], [input_dir+folder+oR+'_trimmed.fq'], [cmd])

            qc[folder].append(oR+'_trimmed.fq')

    if files_to_remove:
        for dep, to_rm in files_to_remove:
            DoitLoader.add_task([], list(dep), ['rm {}'.format(to_rm)], uptodate=[False])

    return qc


def screen_contaminating_dnas(input_dir, qc, bowtie2_indexes, keep_intermediate, remove_ribosomes=False, nprocs_bowtie2=1):
    files_to_remove = []
    screened = dict()
    cont_dnas = ['hg19', 'phiX174']

    if remove_ribosomes:
        cont_dnas += ['SILVA_119.1_SSURef_Nr99_tax_silva', 'SILVA_119_LSURef_tax_silva']

    for folder, Rs in qc.iteritems():
        screened[folder] = tuple()

        for R in Rs:
            outf = R[:R.rfind('.')]
            Rext = R[R.rfind('.'):]
            final = None

            for cont_dna in cont_dnas:
                iR = outf
                suffix = '_{}'.format(cont_dna.replace('_', '-').replace('.', '-'))
                outf += suffix

                cmd = 'bowtie2 -x {} -U {} -S {} -p {} --sensitive-local --un {}'.format(cont_dna, input_dir+folder+iR+Rext, input_dir+folder+outf+'.sam', nprocs_bowtie2, input_dir+folder+outf+'.fastq')
                DoitLoader.add_task([input_dir+folder+outf+'.fastq', input_dir+folder+outf+'.sam'], [input_dir+folder+iR+Rext], [cmd, {'env': {'BOWTIE2_INDEXES': bowtie2_indexes}}])

                if not keep_intermediate:
                    files_to_remove.append(((input_dir+folder+outf+'.fastq', ), ' '.join([input_dir+folder+iR+Rext, input_dir+folder+outf+'.sam'])))

                Rext = '.fastq'
                final = outf+'.fastq'

                cmd = 'fna_len.py {} {} -q --stat'.format(input_dir+folder+outf+'.fastq', input_dir+folder+outf+'.stats')
                DoitLoader.add_task([input_dir+folder+outf+'.stats'], [input_dir+folder+outf+'.fastq'], [cmd])

            screened[folder] += tuple([final])

    if files_to_remove:
        for dep, to_rm in files_to_remove:
            DoitLoader.add_task([], list(dep), ['rm {}'.format(to_rm)], uptodate=[False])

    return screened


def split_and_sort(input_dir, screened, keep_intermediate):
    files_to_remove = []

    for folder, (R1, R2) in screened.iteritems():
        out = R1[:R1.find('.')]
        put = R1[R1.rfind('R1'):R1.rfind('.')].replace('R1', '')

        if (out != R2[:R2.find('.')]) or (put != R2[R2.rfind('R2'):R2.rfind('.')].replace('R2', '')):
            print "\nERROR, split_and_sort() cannot finds common filename!"
            print 'R1', out, put
            print 'R2', R2[:R2.find('.')], R2[R2.rfind('R2'):R2.rfind('.')].replace('R2', '')
            print "\n"
            sys.exit(1)

        cmd = 'split_and_sort.py --R1 {} --R2 {} --prefix {}'.format(input_dir+folder+R1, input_dir+folder+R2, input_dir+folder+out+put)
        DoitLoader.add_task([input_dir+folder+out+put+'.R1.fastq.bz2', input_dir+folder+out+put+'.R2.fastq.bz2', input_dir+folder+out+put+'.UP.fastq.bz2'], [input_dir+folder+R1, input_dir+folder+R2], [cmd])

        if not keep_intermediate:
            files_to_remove.append(((input_dir+folder+out+put+'.R1.fastq.bz2', input_dir+folder+out+put+'.R2.fastq.bz2', input_dir+folder+out+put+'.UP.fastq.bz2'), ' '.join([input_dir+folder+R1, input_dir+folder+R2])))

    if files_to_remove:
        for dep, to_rm in files_to_remove:
            DoitLoader.add_task([], list(dep), ['rm {}'.format(to_rm)], uptodate=[False])


if __name__ == "__main__":
    args = read_params()
    answer = None
    code = 1

    input_dir, inputs = get_inputs(args.input_dir, args.extension) # get input files
    merged = concatenate_reads(input_dir, inputs) # concatenate reads
    qc = quality_control(input_dir, merged, args.keep_intermediate) # trim-galore quality control
    screened = screen_contaminating_dnas(input_dir, qc, args.bowtie2_indexes, args.keep_intermediate, args.remove_ribosomes, args.nprocs_bowtie2) # bowtie2 remove contaminating DNAs
    split_and_sort(input_dir, screened, args.keep_intermediate)

    if args.clean:
        doit_args = ['clean']
        answer = raw_input('Are you sure that you want to clean all intermediate files (Y, N)? ')
    else:
        doit_args = ['-n', str(args.nprocs_main), '--db-file', os.path.join(args.input_dir, '.doit.db')]

        if args.use_threads:
            doit_args += ['-P', 'thread']

    if (answer is None) or (answer.upper()[0] == 'Y'):
        code = DoitLoader.run(doit_args)

    sys.exit(code)
