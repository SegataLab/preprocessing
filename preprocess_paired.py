#!/usr/bin/env python
#Author: Duy Tin Truong (duytin.truong@unitn.it)
#        at CIBIO, University of Trento, Italy


__author__ = 'Duy Tin Truong (duytin.truong@unitn.it)'
__version__ = '0.2'
__date__ = '28 February 2017'


import sys
import os
import glob
import argparse
from doit_loader import DoitLoader


def read_params():
    p = argparse.ArgumentParser()
    p.add_argument('--nprocs_main', required=False, default=1, type=int,
        help='Number of tasks to run in parallel')
    p.add_argument('--nprocs_bowtie2', required=False, default=1, type=int,
        help='Number of bowtie2 processors')
    p.add_argument('-i', '--input_dir', required=True, default=None, type=str,
        help='The input directory, path must be absolute!')
    p.add_argument('--remove_ribosomes', required=False, default=False, action='store_true',
        help='Remove ribosomes for mRNA datasets.')
    p.add_argument('--clean', required=False, default=False, action='store_true',
        help='clean')
    p.add_argument('--keep_intermediate', required=False, default=False, action='store_true',
        help="If specified the script won't remove intermediate files")
    p.add_argument('--use_threads', required=False, action='store_true',
        help='Use threads.')
    p.add_argument('--bowtie2_indexes', required=False, default='/CM/databases/bowtie2_indexes', type=str,
        help='Folder containing the bowtie2 indexes of the genomes to be removed from the samples.')

    args = p.parse_args()

    if not args.input_dir.startswith('/'):
        print 'ERROR, input directory is not an absolute path: "', args.input_dir, '"'
        sys.exit(1)

    if not args.input_dir.endswith('/'):
        args.input_dir += '/'

    return args


def get_inputs(input_dir, ext):
    inputs = dict()
    R1s = sorted(glob.glob('{}*R1*{}'.format(input_dir, ext)))
    R2s = sorted(glob.glob('{}*R2*{}'.format(input_dir, ext)))

    if not (len(R1s) and len(R2s)):
        for folder in os.listdir(input_dir):
            if not folder.endswith('/'):
                folder += '/'

            R1s = sorted(glob.glob('{}{}*R1*{}'.format(input_dir, folder, ext)))
            R2s = sorted(glob.glob('{}{}*R2*{}'.format(input_dir, folder, ext)))
            inputs[folder] = ([a[len(input_dir)+len(folder):] for a in R1s], [a[len(input_dir)+len(folder):] for a in R2s])
    else:
        folder = input_dir.split('/')[-2]

        if not folder.endswith('/'):
                folder += '/'

        inputs[folder] = ([a[len(input_dir):] for a in R1s], [a[len(input_dir):] for a in R2s])
        input_dir = '/'.join(input_dir.split('/')[:-2])

    if not input_dir.endswith('/'):
        input_dir += '/'

    return input_dir, inputs


def concatenate_reads(input_dir, inputs):
    merged = dict()

    for folder, (R1s, R2s) in inputs.iteritems():
        if not folder.endswith('/'):
            folder += '/'

        oR1 = list(set(['_'.join(a.split('R1')[0].split('_')[:-2]) for a in R1s]))[0]+'_R1'
        oR2 = list(set(['_'.join(a.split('R2')[0].split('_')[:-2]) for a in R2s]))[0]+'_R2'

        cmd = 'zcat {}'.format(' '.join([input_dir+folder+a for a in R1s]))
        DoitLoader.add_task([input_dir+folder+oR1+'.fastq'], [input_dir+folder+a for a in R1s], [cmd, {'stdout': input_dir+folder+oR1+'.fastq'}])

        cmd = 'fna_len.py {} {} -q --stat'.format(input_dir+folder+oR1+'.fastq', input_dir+folder+oR1+'.stats')
        DoitLoader.add_task([input_dir+folder+oR1+'.stats'], [input_dir+folder+oR1+'.fastq'], [cmd])

        cmd = 'zcat {}'.format(' '.join([input_dir+folder+a for a in R2s]))
        DoitLoader.add_task([input_dir+folder+oR2+'.fastq'], [input_dir+folder+a for a in R2s], [cmd, {'stdout': input_dir+folder+oR2+'.fastq'}])

        cmd = 'fna_len.py {} {} -q --stat'.format(input_dir+folder+oR2+'.fastq', input_dir+folder+oR2+'.stats')
        DoitLoader.add_task([input_dir+folder+oR2+'.stats'], [input_dir+folder+oR2+'.fastq'], [cmd])

        merged[folder] = (oR1+'.fastq', oR2+'.fastq')

    return merged


def quality_control(input_dir, merged, keep_intermediate):
    qc = dict()

    for folder, (R1, R2) in merged.iteritems():
        if not folder.endswith('/'):
            folder += '/'

        outR1 = '.'.join(R1.split('.')[:-1])
        outR2 = '.'.join(R2.split('.')[:-1])

        cmd = 'trim_galore --nextera --stringency 5 --paired --length 75 --quality 20 --max_n 2 --trim-n --dont_gzip --no_report_file --suppress_warn --output_dir {} {} {}'.format(input_dir+folder, input_dir+folder+R1, input_dir+folder+R2)
        DoitLoader.add_task([outR1+'_val_1.fq', outR2+'_val_2.fq'], [input_dir+folder+R1, input_dir+folder+R2], [cmd])

        if not keep_intermediate:
            cmd = 'rm -f {} {}'.format(input_dir+folder+R1, input_dir+folder+R2)
            DoitLoader.add_task([], [], [cmd])

        cmd = 'fna_len.py {} {} -q --stat'.format(input_dir+folder+outR1+'_val_1.fq', input_dir+folder+outR1+'_val_1.stats')
        DoitLoader.add_task([input_dir+folder+outR1+'_val_1.stats'], [input_dir+folder+outR1+'_val_1.fq'], [cmd])

        cmd = 'fna_len.py {} {} -q --stat'.format(input_dir+folder+outR2+'_val_2.fq', input_dir+folder+outR2+'_val_2.stats')
        DoitLoader.add_task([input_dir+folder+outR2+'_val_2.stats'], [input_dir+folder+outR2+'_val_2.fq'], [cmd])

        qc[folder] = (outR1+'_val_1.fq', outR2+'_val_2.fq')

    return qc


def screen_contaminating_dnas(input_dir, qc, bowtie2_indexes, keep_intermediate, remove_ribosomes=False, nprocs_bowtie2=1):
    files_to_remove = []
    files_to_compress = []
    cont_dnas = ['hg19', 'phiX174']

    if remove_ribosomes:
        cont_dnas += ['SILVA_119.1_SSURef_Nr99_tax_silva', 'SILVA_119_LSURef_tax_silva']

    for folder, (R1, R2) in qc.iteritems():
        if not folder.endswith('/'):
            folder += '/'

        outR1 = R1[:R1.rfind('.')]
        R1ext = R1[R1.rfind('.'):]
        outR2 = R2[:R2.rfind('.')]
        R2ext = R2[R2.rfind('.'):]

        outf = '_'.join(outR1.split('_')[:-3])

        if outf != '_'.join(outR2.split('_')[:-3]):
            print "\nERROR, wrong filename"
            print outR1
            print outR2
            print outf
            print "\n"
            sys.exit(1)

        for cont_dna in cont_dnas:
            inR1, inR2 = outR1, outR2
            suffix = '_{}'.format(cont_dna.replace('_', '-').replace('.', '-'))
            outf += suffix
            outR1 = outf+'.1'
            outR2 = outf+'.2'

            cmd = 'bowtie2 -x {} -1 {} -2 {} -S {} --sensitive-local --un-conc {}'.format(cont_dna, input_dir+folder+inR1+R1ext, input_dir+folder+inR2+R2ext, input_dir+folder+outf+'.sam', input_dir+folder+outf+'.fastq')

            if nprocs_bowtie2 > 1:
                cmd += ' -p {}'.format(nprocs_bowtie2)

            DoitLoader.add_task([input_dir+folder+outR1+R1ext, input_dir+folder+outR2+R2ext,input_dir+folder+outf+'.sam'], [input_dir+folder+inR1+R1ext, input_dir+folder+inR2+R2ext], [cmd, {'env': {'BOWTIE2_INDEXES': bowtie2_indexes}}])

            if not keep_intermediate:
                files_to_remove.append(input_dir+folder+inR1+R1ext)
                files_to_remove.append(input_dir+folder+inR2+R2ext)
                files_to_remove.append(input_dir+folder+outf+'.sam')

            files_to_compress = [input_dir+folder+outR1+R1ext, input_dir+folder+outR2+R2ext]

            R1ext = '.fastq'
            R2ext = '.fastq'

            cmd = 'fna_len.py {} {} -q --stat'.format(input_dir+folder+outR1+R1ext, input_dir+folder+outR1+'.stats')
            DoitLoader.add_task([input_dir+folder+outR1+'.stats'], [input_dir+folder+outR1+R1ext], [cmd])

            cmd = 'fna_len.py {} {} -q --stat'.format(input_dir+folder+outR2+R2ext, input_dir+folder+outR2+'.stats')
            DoitLoader.add_task([input_dir+folder+outR2+'.stats'], [input_dir+folder+outR2+R2ext], [cmd])

    if keep_intermediate:
        cmd = 'rm -f {}'.format(' '.join(files_to_remove))
        DoitLoader.add_task([], [], [cmd])

    cmd = 'bzip2 -z {}'.format(files_to_compress[0])
    DoitLoader.add_task([files_to_compress[0][:files_to_compress[0].rfind('.')]+'.bz2'], [], [cmd])

    cmd = 'bzip2 -z {}'.format(files_to_compress[1])
    DoitLoader.add_task([files_to_compress[1][:files_to_compress[1].rfind('.')]+'.bz2'], [], [cmd])


def main(args):
    ext = '.fastq.gz'
    answer = None
    code = 1

    input_dir, inputs = get_inputs(args.input_dir, ext) # get input files
    merged = concatenate_reads(input_dir, inputs) # concatenate reads
    qc = quality_control(input_dir, merged, args.keep_intermediate) # trim-galore quality control
    screen_contaminating_dnas(input_dir, qc, args.bowtie2_indexes, args.keep_intermediate, args.remove_ribosomes, args.nprocs_bowtie2) # bowtie2 remove contaminating DNAs

    if args.clean:
        doit_args = ['clean']
        answer = raw_input('Are you sure that you want to clean all intermediate files (Y, N)? ')
    else:
        doit_args = ['-n', str(args.nprocs_main), '--db-file', os.path.join(args.input_dir, '.doit.db')]

        if args.use_threads:
            doit_args += ['-P', 'thread']

    if (answer is None) or (answer.upper()[0] == 'Y'):
        code = DoitLoader.run(doit_args)

    return code


if __name__ == "__main__":
    args = read_params()
    sys.exit(main(args))
