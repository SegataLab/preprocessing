#!/usr/bin/env python


__author__ = 'Francesco Asnicar (f.asnicar@unitn.it)'
__version__ = '0.1'
__date__ = '21 Sep 2018'


import os
import sys
import bz2
import glob
import gzip
import argparse
import subprocess as sb


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
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument('-i', '--input_dir', required=True, type=str, help="Path to input directory")
    p.add_argument('-e', '--extension', required=False, default=".fastq.gz", type=str, help="The extension of the raw input files")

    procs = p.add_argument_group("Params for the number of processors to use")
    procs.add_argument('-b', '--nprocs_bowtie2', required=False, default=1, type=int, help="Number of bowtie2 processors")

    rm = p.add_argument_group('Params for what contaminats should be removed')
    rm.add_argument('--rm_hsap', required=False, default=False, action='store_true', help="Remove H. sapiens genome")
    rm.add_argument('--rm_mmus', required=False, default=False, action='store_true', help="Remove M. musculus genome")
    rm.add_argument('--rm_rrna', required=False, default=False, action='store_true', help="Remove rRNA (for mRNA datasets)")

    p.add_argument('-k', '--keep_intermediate', required=False, default=False, action='store_true',
                   help="If specified the script won't remove intermediate files")
    p.add_argument('-x', '--bowtie2_indexes', required=False, default='/shares/CIBIO-Storage/CM/mir/databases/bowtie2_indexes',
                   type=str, help="Folder containing the bowtie2 indexes of the genomes to be removed from the samples")
    p.add_argument('--dry_run', required=False, default=False, action='store_true', help="Print commands do not execute them")
    return p.parse_args()


def check_params(args):
    if not os.path.isdir(args.input_dir):
        error('input folder "{}" does not exists'.format(args.input_dir), exit=True)

    if args.input_dir.endswith('/'):
        args.input_dir = args.input_dir[:-1]


def preflight_check(dry_run=False):
    if dry_run:
        info('preflight_check()\n', init_new_line=True)

    cmds = ['zcat -h', 'fna_len.py -h', 'trim_galore -h', 'bowtie2 -h', 'split_and_sort.py -h', 'cat_stats.py -h']

    for cmd in cmds:
        if dry_run:
            info('{}\n'.format(cmd))
            continue

        try:
            with open(os.devnull, 'w') as devnull:
                sb.check_call(cmd.split(' '), stdout=devnull, stderr=devnull)
        except Exception as e:
            error('preflight_check()\n{}\n{}'.format(cmd, e), exit=True)


def get_inputs(input_dir, ext):
    R1 = sorted(glob.glob(os.path.join(input_dir, '*R1*{}'.format(ext))))
    R2 = sorted(glob.glob(os.path.join(input_dir, '*R2*{}'.format(ext))))

    return (R1, R2)


def concatenate_reads(input_dir, inputs_r1s_r2s, dry_run=False):
    if dry_run:
        info('concatenate_reads()\n', init_new_line=True)

    out_prefix = os.path.basename(input_dir)
    R1s, R2s = inputs_r1s_r2s
    outR1fastq = '{}.R1.fastq'.format(os.path.join(input_dir, out_prefix))
    outR2fastq = '{}.R2.fastq'.format(os.path.join(input_dir, out_prefix))
    outR1stats = '{}.R1.stats'.format(os.path.join(input_dir, out_prefix))
    outR2stats = '{}.R2.stats'.format(os.path.join(input_dir, out_prefix))

    # R1 decompress
    if not os.path.isfile(outR1fastq):
        if dry_run:
            info('{} > {}\n'.format(' '.join(R1s), outR1fastq))
        else:
            g = open(outR1fastq, 'w')

            for inpR in R1s:
                # decompress input file
                if inpR.endswith('.bz2'):
                    with bz2.open(inpR, 'rt') as f:
                        g.write(f.read())
                elif inpR.endswith('.gz'):
                    with gzip.open(inpR, 'rt') as f:
                        g.write(f.read())

            g.close()

    # R1 stats
    if not os.path.isfile(outR1stats):
        cmd = 'fna_len.py -q --stat {} {}'.format(outR1fastq, outR1stats)

        if dry_run:
            info('{}\n'.format(cmd))
        else:
            try:
                sb.check_call(cmd.split(' '))
            except Exception as e:
                if os.path.exists(outR1stats):
                    os.remove(outR1stats)

                error('concatenate_reads()\n{}\n{}'.format(cmd, e), exit=True)

    # R2 decompress
    if not os.path.isfile(outR2fastq):
        if dry_run:
            info('{} > {}\n'.format(' '.join(R2s), outR2fastq))
        else:
            g = open(outR2fastq, 'w')

            for inpR in R2s:
                # decompress input file
                if inpR.endswith('.bz2'):
                    with bz2.open(inpR, 'rt') as f:
                        g.write(f.read())
                elif inpR.endswith('.gz'):
                    with gzip.open(inpR, 'rt') as f:
                        g.write(f.read())

            g.close()

    # R2 stats
    if not os.path.isfile(outR2stats):
        cmd = 'fna_len.py -q --stat {} {}'.format(outR2fastq, outR2stats)

        if dry_run:
            info('{}\n'.format(cmd))
        else:
            try:
                sb.check_call(cmd.split(' '))
            except Exception as e:
                if os.path.exists(outR2stats):
                    os.remove(outR2stats)

                error('concatenate_reads()\n{}\n{}'.format(cmd, e), exit=True)

    return (os.path.basename(outR1fastq), os.path.basename(outR2fastq))


def quality_control(input_dir, merged_r1_r2, keep_intermediate, dry_run=False):
    if dry_run:
        info('quality_control()\n', init_new_line=True)

    qc = []

    for R in merged_r1_r2:
        oR = R[:R.rfind('.')]

        if not os.path.isfile('{}_trimmed.fq'.format(oR)):
            cmd = ('trim_galore --nextera --stringency 5 --length 75 --quality 20 --max_n 2 --trim-n --dont_gzip '
                   '--no_report_file --suppress_warn --output_dir {} {}').format(input_dir, os.path.join(input_dir, R))

            if dry_run:
                info('{}\n'.format(cmd))
            else:
                try:
                    with open(os.devnull, 'w') as devnull:
                        sb.check_call(cmd.split(' '), stdout=devnull)
                except Exception as e:
                    if os.path.exists('{}_trimmed.fq'.format(oR)):
                        os.remove('{}_trimmed.fq'.format(oR))

                    error('quality_control()\n{}\n{}'.format(cmd, e), exit=True)

        if not os.path.isfile('{}_trimmed.stats'.format(oR)):
            cmd = 'fna_len.py -q --stat {0}_trimmed.fq {0}_trimmed.stats'.format(os.path.join(input_dir, oR))

            if dry_run:
                info('{}\n'.format(cmd))
            else:
                try:
                    sb.check_call(cmd.split(' '))
                except Exception as e:
                    if os.path.exists('{}_trimmed.stats'.format(oR)):
                        os.remove('{}_trimmed.stats'.format(oR))

                    error('concatenate_reads()\n{}\n{}'.format(cmd, e), exit=True)

        if not keep_intermediate:
            if os.path.isfile(R):
                if dry_run:
                    info('rm {}\n'.format(R))
                else:
                    os.remove(R)

        qc.append('{}_trimmed.fq'.format(oR))

    return tuple(qc)


def screen_contaminating_dnas(input_dir, qced_r1_r2, bowtie2_indexes, keep_intermediate, rm_hsap, rm_rrna, rm_mmus,
                              nprocs_bowtie2=1, dry_run=False):
    if dry_run:
        info('screen_contaminating_dnas()\n', init_new_line=True)

    screened = []
    to_removes = []
    cont_dnas = ['phiX174']

    if rm_hsap:
        cont_dnas += ['hg19']

    if rm_mmus:
        cont_dnas += ['mmusculus_black6_GCA_000001635_8']

    if rm_rrna:
        cont_dnas += ['SILVA_132_SSURef_Nr99_tax_silva', 'SILVA_132_LSURef_tax_silva']

    for R in qced_r1_r2:
        outf = R[:R.rfind('.')]
        Rext = R[R.rfind('.'):]
        final = None

        for cont_dna in cont_dnas:
            iR = outf
            suffix = '_{}'.format(cont_dna.replace('_', '-').replace('.', '-'))
            outf += suffix

            if not os.path.isfile('{}.fastq'.format(outf)):
                cmd = 'bowtie2 -x {} -U {} -S {}.sam -p {} --sensitive-local --un {}.fastq'.format(cont_dna,
                                                                                                   os.path.join(input_dir, iR + Rext),
                                                                                                   os.path.join(input_dir, outf),
                                                                                                   nprocs_bowtie2,
                                                                                                   os.path.join(input_dir, outf))

                if dry_run:
                    info('{}\n'.format(cmd))
                else:
                    try:
                        with open(os.devnull, 'w') as devnull:
                            sb.check_call(cmd.split(' '), stdout=devnull, env={'BOWTIE2_INDEXES': bowtie2_indexes})
                    except Exception as e:
                        for i in '{0}.sam {0}.fastq'.format(os.path.join(input_dir, outf)).split(' '):
                            if os.path.exists(i):
                                os.reomve(i)

                        error('quality_control()\n{}\n{}'.format(cmd, e), exit=True)

            if not os.path.isfile('{}.stats'.format(outf)):
                cmd = 'fna_len.py -q --stat {0}.fastq {0}.stats'.format(os.path.join(input_dir, outf))

                if dry_run:
                    info('{}\n'.format(cmd))
                else:
                    try:
                        sb.check_call(cmd.split(' '))
                    except Exception as e:
                        if os.path.exists('{}.stats'.format(outf)):
                            os.remove('{}.stats'.format(outf))

                        error('concatenate_reads()\n{}\n{}'.format(cmd, e), exit=True)

            if not keep_intermediate:
                to_removes.append(iR + Rext)
                to_removes.append('{}.sam'.format(outf))

            Rext = '.fastq'
            final = outf + '.fastq'

        screened.append(final)

    if to_removes:
        for to_remove in to_removes:
            if os.path.isfile(to_remove):
                if dry_run:
                    info('rm {}'.format(to_remove))
                else:
                    os.remove(to_remove)

    return tuple(screened)


def split_and_sort(input_dir, screened_r1_r2, keep_intermediate, dry_run=False):
    if dry_run:
        info('split_and_sort()\n', init_new_line=True)

    R1, R2 = screened_r1_r2
    out = R1[:R1.find('.')]
    put = R1[R1.rfind('R1'):R1.rfind('.')].replace('R1', '')

    if (out != R2[:R2.find('.')]) or (put != R2[R2.rfind('R2'):R2.rfind('.')].replace('R2', '')):
        error('split_and_sort() ::: cannot finds common filename!\n    R1: {}\n    R2: {}\n   out: {}\n   put: {}'
              .format(R1, R2, out, put), exit=True)

    if not (os.path.isfile(out + put + '_R1.fastq.bz2') and
            os.path.isfile(out + put + '_R2.fastq.bz2') and
            os.path.isfile(out + put + '_UN.fastq.bz2')):
        cmd = 'split_and_sort.py --R1 {} --R2 {} --prefix {}'.format(os.path.join(input_dir, R1),
                                                                     os.path.join(input_dir, R2),
                                                                     os.path.join(input_dir, out + put))

        if dry_run:
            info('{}\n'.format(cmd))
        else:
            try:
                sb.check_call(cmd.split(' '))
            except Exception as e:
                for i in [os.path.isfile(out + put + '_R1.fastq.bz2'),
                          os.path.isfile(out + put + '_R2.fastq.bz2'),
                          os.path.isfile(out + put + '_UN.fastq.bz2')]:
                    if os.path.exists(i):
                        os.remove(i)

                error('split_and_sort()\n{}\n{}'.format(cmd, e), exit=True)

    if not os.path.isfile(out + put + '_R1.stats'):
        cmd = 'fna_len.py -q --stat {0}_R1.fastq.bz2 {0}_R1.stats'.format(os.path.join(input_dir, out + put))

        if dry_run:
            info('{}\n'.format(cmd))
        else:
            try:
                sb.check_call(cmd.split(' '))
            except Exception as e:
                if os.path.exists(out + put + '_R1.stats'):
                    os.remove(out + put + '_R1.stats')

                error('split_and_sort()\n{}\n{}'.format(cmd, e), exit=True)

    if not os.path.isfile(out + put + '_R2.stats'):
        cmd = 'fna_len.py -q --stat {0}_R2.fastq.bz2 {0}_R2.stats'.format(os.path.join(input_dir, out + put))

        if dry_run:
            info('{}\n'.format(cmd))
        else:
            try:
                sb.check_call(cmd.split(' '))
            except Exception as e:
                if os.path.exists(out + put + '_R2.stats'):
                    os.remove(out + put + '_R2.stats')

                error('split_and_sort()\n{}\n{}'.format(cmd, e), exit=True)

    if not os.path.isfile(out + put + '_UN.stats'):
        cmd = 'fna_len.py -q --stat {0}_UN.fastq.bz2 {0}_UN.stats'.format(os.path.join(input_dir, out + put))

        if dry_run:
            info('{}\n'.format(cmd))
        else:
            try:
                sb.check_call(cmd.split(' '))
            except Exception as e:
                if os.path.exists(out + put + '_UN.stats'):
                    os.remove(out + put + '_UN.stats')

                error('split_and_sort()\n{}\n{}'.format(cmd, e), exit=True)

    if not keep_intermediate:
        for R in screened_r1_r2:
            if os.path.isfile(R):
                if dry_run:
                    info('rm {}\n'.format(R))
                else:
                    os.remove(R)

    if not os.path.isfile(out + put + '_summary.stats'):
        cmd = 'cat_stats.py -i {} -o {}'.format(input_dir, os.path.join(input_dir, out + put + '_summary.stats'))

        if dry_run:
            info('{}\n'.format(cmd))
        else:
            try:
                sb.check_call(cmd.split(' '))
            except Exception as e:
                if os.path.exists(out + put + '_summary.stats'):
                    os.remove(out + put + '_summary.stats')

                error('split_and_sort()\n{}\n{}'.format(cmd, e), exit=True)


if __name__ == "__main__":
    args = read_params()
    check_params(args)

    preflight_check(dry_run=args.dry_run)  # check that all the needed software are available

    inputs_r1s_r2s = get_inputs(args.input_dir, args.extension)  # get input files

    if args.dry_run:
        info('inputs_r1s_r2s: {}\n'.format(inputs_r1s_r2s), init_new_line=True)

    merged_r1_r2 = concatenate_reads(args.input_dir, inputs_r1s_r2s, dry_run=args.dry_run)  # concatenate reads

    if args.dry_run:
        info('merged_r1_r2: {}\n'.format(merged_r1_r2))

    qced_r1_r2 = quality_control(args.input_dir, merged_r1_r2, args.keep_intermediate, dry_run=args.dry_run)  # quality control

    if args.dry_run:
        info('qced_r1_r2: {}\n'.format(qced_r1_r2))

    # bowtie2 remove contaminating DNAs
    screened_r1_r2 = screen_contaminating_dnas(args.input_dir, qced_r1_r2, args.bowtie2_indexes, args.keep_intermediate,
                                               args.rm_hsap, args.rm_rrna, args.rm_mmus,
                                               nprocs_bowtie2=args.nprocs_bowtie2, dry_run=args.dry_run)

    if args.dry_run:
        info('screened_r1_r2: {}\n'.format(screened_r1_r2))

    split_and_sort(args.input_dir, screened_r1_r2, args.keep_intermediate, dry_run=args.dry_run)

    sys.exit(0)
