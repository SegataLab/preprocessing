#!/usr/bin/env python


__author__ = 'Duy Tin Truong (duytin.truong@unitn.it), Francesco Asnicar (f.asnicar@unitn.it)'
__version__ = '0.6'
__date__ = '4 July 2017'


import sys
import os
import glob
import argparse
from doit_loader import DoitLoader


def info(s, init_new_line=False, exit=False, exit_value=0):
    if init_new_line:
        sys.stdout.write('\n')

    sys.stdout.write('{}'.format(s))
    sys.stdout.flush()

    if exit:
        sys.exit(exit_value)


def error(s, init_new_line=False, exit=False, exit_value=1):
    if init_new_line:
        sys.stderr.write('\n')

    sys.stderr.write('[e] {}\n'.format(s))
    sys.stderr.flush()

    if exit:
        sys.exit(exit_value)


def read_params():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input_dir', required=True, default=None, type=str, help='The input directory, path must be absolute!')
    p.add_argument('-e', '--extension', required=False, default=".fastq.gz", type=str, help='The extension of the raw input files. Default ".fastq.gz"')
    p.add_argument('-t', '--use_threads', required=False, action='store_true', help='Use threads')
    p.add_argument('-m', '--nprocs_main', required=False, default=1, type=int, help='Number of tasks to run in parallel')
    p.add_argument('-b', '--nprocs_bowtie2', required=False, default=1, type=int, help='Number of bowtie2 processors')
    p.add_argument('-r', '--remove_ribosomes', required=False, default=False, action='store_true', help='Remove ribosomes for mRNA datasets')
    p.add_argument('-c', '--clean', required=False, default=False, action='store_true', help='clean')
    p.add_argument('-k', '--keep_intermediate', required=False, default=False, action='store_true', help="If specified the script won't remove intermediate files")
    p.add_argument('-x', '--bowtie2_indexes', required=False, default='/CM/databases/bowtie2_indexes', type=str, help='Folder containing the bowtie2 indexes of the genomes to be removed from the samples. Default "/CM/databases/bowtie2_indexes"')

    return p.parse_args()


def check_params(args):
    if not args.input_dir.startswith('/'):
        error('input directory is not an absolute path: "{}"'.format(args.input_dir), exit=True)

    if not args.input_dir.endswith('/'):
        args.input_dir += '/'


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
        merged[folder] = list()
        out_prefix = set(['_'.join(a.split('_')[:3]) for a in R1s + R2s])

        if len(out_prefix) == 1:
            out_prefix = list(out_prefix)[0]
        else:
            out_prefix = set(['_'.join(a.split('_')[:2]) for a in R1s + R2s])

            if len(out_prefix) == 1:
                out_prefix = list(out_prefix)[0]
            else:
                error('concatenate_reads() cannot finds common filename!\n    {}'.format('\n    '.join('{} "{}"'.format(a, b) for a, b in zip(['out_prefix', 'folder', 'R1s', 'R2s'], [out_prefix, folder, R1s, R2s]))), exit=True)

        # info('    folder: {}\n'.format(folder))
        # info('out_prefix: {}\n\n'.format(out_prefix))

        for Rs, oR in zip([R1s, R2s], [out_prefix+'.R1', out_prefix+'.R2']):
            cmd = 'zcat {}'.format(' '.join([input_dir+folder+a for a in Rs]))
            DoitLoader.add_task([input_dir+folder+oR+'.fastq'], [input_dir+folder+a for a in Rs], [cmd, {'stdout': input_dir+folder+oR+'.fastq'}])

            cmd = 'fna_len.py {} {} -q --stat'.format(input_dir+folder+oR+'.fastq', input_dir+folder+oR+'.stats')
            DoitLoader.add_task([input_dir+folder+oR+'.stats'], [input_dir+folder+oR+'.fastq'], [cmd])

            merged[folder].append(oR+'.fastq')

    return merged


def quality_control(input_dir, merged, keep_intermediate):
    qc = dict()

    for folder, Rs in merged.iteritems():
        qc[folder] = list()

        for R in Rs:
            oR = R[:R.rfind('.')]

            cmd = 'trim_galore --nextera --stringency 5 --length 75 --quality 20 --max_n 2 --trim-n --dont_gzip --no_report_file --suppress_warn --output_dir {} {}'.format(input_dir+folder, input_dir+folder+R)
            DoitLoader.add_task([input_dir+folder+oR+'_trimmed.fq'], [input_dir+folder+R], [cmd, {'stdout': '/dev/null'}])

            cmd = 'fna_len.py {} {} -q --stat'.format(input_dir+folder+oR+'_trimmed.fq', input_dir+folder+oR+'_trimmed.stats')
            DoitLoader.add_task([input_dir+folder+oR+'_trimmed.stats'], [input_dir+folder+oR+'_trimmed.fq'], [cmd])

            if not keep_intermediate:
                DoitLoader.add_task([], [input_dir+folder+oR+'_trimmed.fq', input_dir+folder+oR+'_trimmed.stats'], ['rm {}'.format(input_dir+folder+R)])

            qc[folder].append(oR+'_trimmed.fq')

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
            to_remove_dep = []
            to_remove = []
            outf = R[:R.rfind('.')]
            Rext = R[R.rfind('.'):]
            final = None

            for cont_dna in cont_dnas:
                iR = outf
                suffix = '_{}'.format(cont_dna.replace('_', '-').replace('.', '-'))
                outf += suffix

                cmd = 'bowtie2 -x {} -U {} -S {} -p {} --sensitive-local --un {}'.format(cont_dna, input_dir+folder+iR+Rext, input_dir+folder+outf+'.sam', nprocs_bowtie2, input_dir+folder+outf+'.fastq')
                DoitLoader.add_task([input_dir+folder+outf+'.fastq', input_dir+folder+outf+'.sam'], [input_dir+folder+iR+Rext], [cmd, {'stdout': '/dev/null', 'env': {'BOWTIE2_INDEXES': bowtie2_indexes}}])

                # info('bowtie2\n')
                # info('     input: {}\n'.format(input_dir+folder+iR+Rext))
                # info('    output: {}\n'.format(input_dir+folder+outf+'.fastq'))
                # info('            {}\n\n'.format(input_dir+folder+outf+'.sam'))

                cmd = 'fna_len.py {} {} -q --stat'.format(input_dir+folder+outf+'.fastq', input_dir+folder+outf+'.stats')
                DoitLoader.add_task([input_dir+folder+outf+'.stats'], [input_dir+folder+outf+'.fastq'], [cmd])

                # info('fna_len.py\n')
                # info('     input: {}\n'.format(input_dir+folder+outf+'.fastq'))
                # info('    output: {}\n\n'.format(input_dir+folder+outf+'.stats'))

                if not keep_intermediate:
                    # DoitLoader.add_task([], [input_dir+folder+outf+'.fastq', input_dir+folder+outf+'.stats'], ['rm {} {}'.format(input_dir+folder+iR+Rext, input_dir+folder+outf+'.sam')])

                    to_remove_dep.append(input_dir+folder+outf+'.stats')
                    to_remove.append(input_dir+folder+iR+Rext)
                    to_remove.append(input_dir+folder+outf+'.sam')

                    # info('rm\n')
                    # info('    remove: {}\n'.format(input_dir+folder+iR+Rext))
                    # info('            {}\n'.format(input_dir+folder+outf+'.sam'))
                    # info('       dep: {}\n'.format(input_dir+folder+outf+'.fastq'))
                    # info('            {}\n\n'.format(input_dir+folder+outf+'.stats'))

                Rext = '.fastq'
                final = outf+'.fastq'

            if not keep_intermediate:
                to_remove_dep.append(input_dir+folder+outf+'.fastq')

            if (not keep_intermediate) and to_remove and to_remove_dep:
                DoitLoader.add_task([], to_remove_dep, ['rm {}'.format(' '.join(to_remove))])

                # info('rm\n')
                # info('        to_remove: {}\n'.format('\n                   '.join(to_remove)))
                # info('    to_remove_dep: {}\n\n'.format('\n                   '.join(to_remove_dep)))

            screened[folder] += tuple([final])

    return screened


def split_and_sort(input_dir, screened, keep_intermediate):
    files_to_remove = []

    for folder, (R1, R2) in screened.iteritems():
        out = R1[:R1.find('.')]
        put = R1[R1.rfind('R1'):R1.rfind('.')].replace('R1', '')

        if (out != R2[:R2.find('.')]) or (put != R2[R2.rfind('R2'):R2.rfind('.')].replace('R2', '')):
            error('split_and_sort() cannot finds common filename!\n    {}'.format('\n    '.join(['{} "{}"'.format(a, b) for a, b in zip(['R1', 'R2'], [out+put, R2[:R2.find('.')]+R2[R2.rfind('R2'):R2.rfind('.')].replace('R2', '')])])), exit=True)

        cmd = 'split_and_sort.py --R1 {} --R2 {} --prefix {}'.format(input_dir+folder+R1, input_dir+folder+R2, input_dir+folder+out+put)
        DoitLoader.add_task([input_dir+folder+out+put+'.R1.fastq.bz2', input_dir+folder+out+put+'.R2.fastq.bz2', input_dir+folder+out+put+'.UP.fastq.bz2'], [input_dir+folder+R1, input_dir+folder+R2], [cmd])

        cmd = 'fna_len.py {} {} -q --stat'.format(input_dir+folder+out+put+'.R1.fastq.bz2', input_dir+folder+out+put+'.R1.stats')
        DoitLoader.add_task([input_dir+folder+out+put+'.R1.stats'], [input_dir+folder+out+put+'.R1.fastq.bz2'], [cmd])

        cmd = 'fna_len.py {} {} -q --stat'.format(input_dir+folder+out+put+'.R2.fastq.bz2', input_dir+folder+out+put+'.R2.stats')
        DoitLoader.add_task([input_dir+folder+out+put+'.R2.stats'], [input_dir+folder+out+put+'.R2.fastq.bz2'], [cmd])

        cmd = 'fna_len.py {} {} -q --stat'.format(input_dir+folder+out+put+'.UP.fastq.bz2', input_dir+folder+out+put+'.UP.stats')
        DoitLoader.add_task([input_dir+folder+out+put+'.UP.stats'], [input_dir+folder+out+put+'.UP.fastq.bz2'], [cmd])

        if not keep_intermediate:
            DoitLoader.add_task([], [input_dir+folder+out+put+'.R1.fastq.bz2', input_dir+folder+out+put+'.R2.fastq.bz2', input_dir+folder+out+put+'.UP.fastq.bz2', input_dir+folder+out+put+'.R1.stats', input_dir+folder+out+put+'.R2.stats', input_dir+folder+out+put+'.UP.stats'], ['rm {} {}'.format(input_dir+folder+R1, input_dir+folder+R2)])

        cmd = 'cat_stats.py -i {} -o {}'.format(input_dir+folder, input_dir+folder+out+put+'_summary.stats')
        DoitLoader.add_task([input_dir+folder+out+put+'.R1.stats', input_dir+folder+out+put+'.R2.stats', input_dir+folder+out+put+'.UP.stats'], [input_dir+folder+out+put+'_summary.stats'], [cmd])


if __name__ == "__main__":
    args = read_params()
    check_params(args)
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
