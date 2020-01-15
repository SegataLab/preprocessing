#!/usr/bin/env python


__author__ = 'Francesco Asnicar (f.asnicar@unitn.it)'
__version__ = '0.2.4'
__date__ = '30 December 2019'


import os
import sys
import bz2
import glob
import gzip
import argparse
import subprocess as sb
import multiprocessing as mp
import time
import shutil


if sys.version_info[0] < 3:
    raise Exception("Preprocessing requires Python 3, your current Python version is {}.{}.{}"
                    .format(sys.version_info[0], sys.version_info[1], sys.version_info[2]))


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
    p.add_argument('-e', '--extension', required=False, default=".fastq.gz",
                   choices=[".fastq.gz", ".fq.gz", ".fastq.bz2", ".fq.bz2"], help="The extension of the raw input files")
    p.add_argument('-s', '--samplename', required=False, default="", help="Specify the sample name")

    p.add_argument('-f', '--forward', required=False, default="R1",
                   help="Identifier to distinguish forward reads in the input folder")
    p.add_argument('-r', '--reverse', required=False, default="R2",
                   help="Identifier to distinguish reverse reads in the input folder")

    procs = p.add_argument_group("Params for the number of processors to use")
    procs.add_argument('-n', '--nproc', required=False, default=2, type=int, help="Number of threads to use")
    procs.add_argument('-b', '--nproc_bowtie2', required=False, default=2, type=int, help="Number of bowtie2 processors")

    rm = p.add_argument_group('Params for what contaminats should be removed')
    rm.add_argument('--rm_hsap', required=False, default=False, action='store_true', help="Remove H. sapiens genome")
    rm.add_argument('--rm_mmus', required=False, default=False, action='store_true', help="Remove M. musculus genome")
    rm.add_argument('--rm_rrna', required=False, default=False, action='store_true', help="Remove rRNA (for mRNA datasets)")

    p.add_argument('-k', '--keep_intermediate', required=False, default=False, action='store_true',
                   help="If specified the script won't remove intermediate files")
    p.add_argument('-x', '--bowtie2_indexes', required=False, default='/shares/CIBIO-Storage/CM/mir/databases/bowtie2_indexes',
                   type=str, help="Folder containing the bowtie2 indexes of the genomes to be removed from the samples")
    p.add_argument('--dry_run', required=False, default=False, action='store_true', help="Print commands do not execute them")
    p.add_argument('--verbose', required=False, default=False, action='store_true', help="Makes preprocessing verbose")
    p.add_argument('-v', '--version', action='version', version='Preprocessing version {} ({})'.format(__version__, __date__),
                   help="Prints the current Preprocessing version and exit")
    return p.parse_args()


def check_params(args):
    if not os.path.isdir(args.input_dir):
        error('input folder "{}" does not exists'.format(args.input_dir), exit=True)

    if args.input_dir.endswith('/'):
        args.input_dir = args.input_dir[:-1]


def initt(terminating_):
    # This places terminating in the global namespace of the worker subprocesses.
    # This allows the worker function to access `terminating` even though it is
    # not passed as an argument to the function.
    global terminating
    terminating = terminating_


def preflight_check(dry_run=False, verbose=False):
    if dry_run or verbose:
        info('preflight_check()\n', init_new_line=True)

    cmds = ['fna_len.py -h', 'trim_galore --help', 'bowtie2 -h',
            'split_and_sort.new.py -h',  # 'split_and_sort.py -h',
            'cat_stats.py -h']

    for cmd in cmds:
        if dry_run or verbose:
            info('{}\n'.format(cmd))

        if dry_run:
            continue

        try:
            with open(os.devnull, 'w') as devnull:
                sb.check_call(cmd.split(' '), stdout=devnull, stderr=devnull)
        except Exception as e:
            error('preflight_check()\n{}\n{}'.format(cmd, e), exit=True)


def get_inputs(input_dir, fwd, rev, sn, ext, verbose=False):
    if verbose:
        info('get_inputs()\n', init_new_line=True)

    R1 = sorted([os.path.join(input_dir, i) for i in os.listdir(input_dir) if (fwd in i.replace(sn, '')) and i.endswith(ext)])
    R2 = sorted([os.path.join(input_dir, i) for i in os.listdir(input_dir) if (rev in i.replace(sn, '')) and i.endswith(ext)])

    return (R1, R2)


def concatenate_reads(input_dir, inputs_r1s_r2s, nproc=1, dry_run=False, verbose=False):
    if dry_run or verbose:
        info('concatenate_reads()\n', init_new_line=True)

    out_prefix = os.path.basename(input_dir)
    R1s, R2s = inputs_r1s_r2s
    outR1fastq = '{}.R1.fastq'.format(os.path.join(input_dir, out_prefix))
    outR2fastq = '{}.R2.fastq'.format(os.path.join(input_dir, out_prefix))
    outR1stats = '{}.R1.stats'.format(os.path.join(input_dir, out_prefix))
    outR2stats = '{}.R2.stats'.format(os.path.join(input_dir, out_prefix))
    tasks = [(R1s, outR1fastq, outR1stats, dry_run, verbose),
             (R2s, outR2fastq, outR2stats, dry_run, verbose)]
    terminating = mp.Event()

    with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
        try:
            fastqs = [a for a in pool.imap_unordered(concatenate_reads_mp, tasks, chunksize=1)]
        except Exception as e:
            error('concatenate_reads()\n    tasks: {}\n    e: {}'.format(tasks, str(e)), init_new_line=True, exit=True)

    return tuple(fastqs)


def concatenate_reads_mp(x):
    if not terminating.is_set():
        try:
            inps, out_fastq, out_stats, dry_run, verbose = x

            # decompress
            if not os.path.isfile(out_fastq):
                if dry_run or verbose:
                    info('cat {} > {}\n'.format(' '.join(inps), out_fastq))

                if not dry_run:
                    g = open(out_fastq, 'wb')

                    for inpR in inps:
                        # decompress input file
                        if inpR.endswith('.bz2'):
                            with bz2.open(inpR, 'rb') as f:
                                # g.write(f.read())
                                shutil.copyfileobj(f, g)
                        elif inpR.endswith('.gz'):
                            with gzip.open(inpR, 'rb') as f:
                                # g.write(f.read())
                                shutil.copyfileobj(f, g)

                    g.close()

            # stats
            if not os.path.isfile(out_stats):
                cmd = 'fna_len.py -q --stat {} {}'.format(out_fastq, out_stats)

                if dry_run or verbose:
                    info('{}\n'.format(cmd))

                if not dry_run:
                    sb.check_call(cmd.split(' '))

            return os.path.basename(out_fastq)
        except Exception as e:
            terminating.set()

            for i in [out_fastq, out_stats]:
                if os.path.exists(i):
                    os.remove(i)

            error('concatenate_reads_mp()\n    x: {}\n    e: {}'.format(x, str(e)), init_new_line=True)
            raise
    else:
        terminating.set()


def quality_control(input_dir, merged_r1_r2, keep_intermediate, sn, nproc=1, dry_run=False, verbose=False):
    if dry_run or verbose:
        info('quality_control()\n', init_new_line=True)

    tasks = zip(merged_r1_r2, [input_dir] * len(merged_r1_r2),
                [keep_intermediate] * len(merged_r1_r2),
                [dry_run] * len(merged_r1_r2),
                [verbose] * len(merged_r1_r2))
    terminating = mp.Event()

    with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
        try:
            qc = [a for a in pool.imap_unordered(quality_control_mp, tasks, chunksize=1)]
        except Exception as e:
            error('quality_control()\ntasks: {}\n    e: {}'.format(tasks, e), init_new_line=True, exit=True)

    r1 = [i for i in qc if "R1" in i.replace(sn, '')]
    r2 = [i for i in qc if "R2" in i.replace(sn, '')]

    if len(r1) > 1:
        error('quality_control(): more than one R1 detected: [{}]'.format(', '.join(r1)), exit=True)

    if len(r2) > 1:
        error('quality_control(): more than one R2 detected: [{}]'.format(', '.join(r2)), exit=True)

    return tuple([r1[0], r2[0]])


def quality_control_mp(x):
    if not terminating.is_set():
        try:
            R, input_dir, keep_intermediate, dry_run, verbose = x

            oR = R[:R.rfind('.')]

            if not os.path.isfile('{}_trimmed.fq'.format(os.path.join(input_dir, oR))):
                cmd = ('trim_galore --nextera --stringency 5 --length 75 --quality 20 --max_n 2 --trim-n --dont_gzip '
                       '--no_report_file --suppress_warn --output_dir {} {}').format(input_dir, os.path.join(input_dir, R))

                # command for Moreno, no --nextera
                #cmd = ('trim_galore --stringency 5 --length 75 --quality 20 --max_n 2 --trim-n --dont_gzip '
                #       '--no_report_file --suppress_warn --output_dir {} {}').format(input_dir, os.path.join(input_dir, R))

                if dry_run or verbose:
                    info('{}\n'.format(cmd))

                if not dry_run:
                    with open(os.devnull, 'w') as devnull:
                        sb.check_call(cmd.split(' '), stdout=devnull)

            if not os.path.isfile('{}_trimmed.stats'.format(os.path.join(input_dir, oR))):
                cmd = 'fna_len.py -q --stat {0}_trimmed.fq {0}_trimmed.stats'.format(os.path.join(input_dir, oR))

                if dry_run or verbose:
                    info('{}\n'.format(cmd))

                if not dry_run:
                    sb.check_call(cmd.split(' '))

            # remove([R], keep_intermediate, folder=input_dir, dry_run=dry_run, verbose=verbose)
            return '{}_trimmed.fq'.format(oR)
        except Exception as e:
            terminating.set()

            for i in ['{}_trimmed.fq'.format(os.path.join(input_dir, oR)), '{}_trimmed.stats'.format(os.path.join(input_dir, oR))]:
                if os.path.exists(i):
                    os.remove(i)

            error('quality_control_mp()\n    x: {}\n    e: {}'.format(x, e), init_new_line=True)
            raise
    else:
        terminating.set()


def screen_contaminating_dnas(input_dir, qced_r1_r2, bowtie2_indexes, keep_intermediate, rm_hsap, rm_rrna, rm_mmus,
                              nprocs_bowtie2=1, dry_run=False, verbose=False):
    if dry_run or verbose:
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

            if not os.path.isfile('{}.fastq'.format(os.path.join(input_dir, outf))):
                cmd = ('bowtie2 -x {} -U {} -p {} --sensitive-local --un {}.fastq'
                       .format(os.path.join(bowtie2_indexes, cont_dna),
                               os.path.join(input_dir, iR + Rext),
                               nprocs_bowtie2,
                               os.path.join(input_dir, outf)))

                if dry_run or verbose:
                    info('{}\n'.format(cmd))

                if not dry_run:
                    try:
                        with open(os.devnull, 'w') as devnull:
                            sb.check_call(cmd.split(' '), stdout=devnull, stderr=devnull)
                    except Exception as e:
                        if os.path.exists('{0}.fastq'.format(os.path.join(input_dir, outf))):
                            os.remove('{0}.fastq'.format(os.path.join(input_dir, outf)))

                        error('screen_contaminating_dnas()\n{}\n{}'.format(cmd, e), exit=True)

            if not os.path.isfile('{}.stats'.format(os.path.join(input_dir, outf))):
                cmd = 'fna_len.py -q --stat {0}.fastq {0}.stats'.format(os.path.join(input_dir, outf))

                if dry_run or verbose:
                    info('{}\n'.format(cmd))

                if not dry_run:
                    try:
                        sb.check_call(cmd.split(' '))
                    except Exception as e:
                        if os.path.exists('{}.stats'.format(outf)):
                            os.remove('{}.stats'.format(outf))

                        error('concatenate_reads()\n{}\n{}'.format(cmd, e), exit=True)

            if not keep_intermediate:
                to_removes.append(iR + Rext)

            Rext = '.fastq'
            final = outf + '.fastq'

        screened.append(final)

    remove(to_removes, keep_intermediate, folder=input_dir, dry_run=dry_run, verbose=verbose)
    return tuple(screened)


def split_and_sort(input_dir, screened_r1_r2, keep_intermediate, nproc=1, dry_run=False, verbose=False):
    if dry_run or verbose:
        info('split_and_sort()\n', init_new_line=True)

    R1, R2 = screened_r1_r2
    out = R1[:R1.find('.')]
    put = '_'.join([a for a in
                    R1[R1.rfind('R1'):R1.rfind('.')].replace('R1', '').replace('trimmed', '').replace('phiX174', '').replace('hg19', '').split('_')
                    if a])
    outR2 = R2[:R2.find('.')]
    putR2 = '_'.join([a for a in
                      R2[R2.rfind('R2'):R2.rfind('.')].replace('R2', '').replace('trimmed', '').replace('phiX174', '').replace('hg19', '').split('_')
                      if a])

    if (out != outR2) or (put != putR2):
        error('split_and_sort() ::: cannot finds common filename!\n    R1: {}\n    R2: {}\n   out: {}\n   put: {}'
              .format(R1, R2, out, put), exit=True)

    if not (os.path.isfile(out + put + '_R1.fastq.bz2') and
            os.path.isfile(out + put + '_R2.fastq.bz2') and
            os.path.isfile(out + put + '_UN.fastq.bz2')):
        # cmd = 'split_and_sort.py --R1 {} --R2 {} --prefix {}'.format(os.path.join(input_dir, R1),
        #                                                              os.path.join(input_dir, R2),
        #                                                              os.path.join(input_dir, out + put))
        cmd = 'split_and_sort.new.py --R1 {} --R2 {} --prefix {}'.format(os.path.join(input_dir, R1),
                                                                         os.path.join(input_dir, R2),
                                                                         os.path.join(input_dir, out + put))

        if dry_run or verbose:
            info('{}\n'.format(cmd))

        if not dry_run:
            try:
                sb.check_call(cmd.split(' '))
            except Exception as e:
                for i in [os.path.isfile(out + put + '_R1.fastq.bz2'),
                          os.path.isfile(out + put + '_R2.fastq.bz2'),
                          os.path.isfile(out + put + '_UN.fastq.bz2')]:
                    if os.path.exists(i):
                        os.remove(i)

                error('split_and_sort()\n{}\n{}'.format(cmd, e), exit=True)

    tasks = []

    if not os.path.isfile(out + put + '_R1.stats'):
        tasks.append(('fna_len.py -q --stat {0}_R1.fastq.bz2 {0}_R1.stats'.format(os.path.join(input_dir, out + put)),
                      '{0}_R1.stats'.format(os.path.join(input_dir, out + put)), dry_run, verbose))

    if not os.path.isfile(out + put + '_R2.stats'):
        tasks.append(('fna_len.py -q --stat {0}_R2.fastq.bz2 {0}_R2.stats'.format(os.path.join(input_dir, out + put)),
                      '{0}_R2.stats'.format(os.path.join(input_dir, out + put)), dry_run, verbose))

    if not os.path.isfile(out + put + '_UN.stats'):
        tasks.append(('fna_len.py -q --stat {0}_UN.fastq.bz2 {0}_UN.stats'.format(os.path.join(input_dir, out + put)),
                      '{0}_UN.stats'.format(os.path.join(input_dir, out + put)), dry_run, verbose))

    terminating = mp.Event()

    with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
        try:
            [_ for _ in pool.imap_unordered(split_and_sort_mp, tasks, chunksize=1)]
        except Exception as e:
            error('split_and_sort()\ntasks: {}\n    e: {}'.format(tasks, e), init_new_line=True, exit=True)

    if not os.path.isfile(out + put + '_summary.stats'):
        cmd = 'cat_stats.py -i {} -o {}'.format(input_dir, os.path.join(input_dir, out + put + '_summary.stats'))

        if dry_run or verbose:
            info('{}\n'.format(cmd))

        if not dry_run:
            try:
                sb.check_call(cmd.split(' '))
            except Exception as e:
                if os.path.exists(out + put + '_summary.stats'):
                    os.remove(out + put + '_summary.stats')

                error('split_and_sort()\n{}\n{}'.format(cmd, e), exit=True)

    remove(screened_r1_r2, keep_intermediate, folder=input_dir, dry_run=dry_run, verbose=verbose)
    return (out + put + '_R1.fastq.bz2', out + put + '_R2.fastq.bz2', out + put + '_UN.fastq.bz2')


def split_and_sort_mp(x):
    if not terminating.is_set():
        try:
            cmd, output, dry_run, verbose = x

            if dry_run or verbose:
                info('{}\n'.format(cmd))

            if not dry_run:
                sb.check_call(cmd.split(' '))
        except Exception as e:
            terminating.set()

            if os.path.exists(output):
                os.remove(output)

            error('split_and_sort_mp()\n    x: {}\n    e: {}'.format(x, e), init_new_line=True)
            raise
    else:
        terminating.set()


def remove(to_remove, keep_intermediate, folder=None, dry_run=False, verbose=False):
    if verbose:
        info('remove()\n', init_new_line=True)

    if not args.keep_intermediate:
        for r in to_remove:
            rf = os.path.join(folder, r) if folder else r

            if os.path.isfile(rf):
                if args.verbose:
                    info('rm {}\n'.format(rf))

                if not args.dry_run:
                    os.remove(rf)


if __name__ == "__main__":
    t0 = time.time()
    args = read_params()

    if args.verbose:
        info('Preprocessing version {} ({})\n'.format(__version__, __date__))
        info('Command line: {}\n'.format(' '.join(sys.argv)), init_new_line=True)

    check_params(args)
    preflight_check(dry_run=args.dry_run, verbose=args.verbose)
    inputs_r1s_r2s = get_inputs(args.input_dir, args.forward, args.reverse, args.samplename, args.extension, verbose=args.verbose)

    if (len(inputs_r1s_r2s[0]) == 0) or (len(inputs_r1s_r2s[1]) == 0):
        error('No input files detected!\nR1s: {}\nR2s: {}'.format(inputs_r1s_r2s[0], inputs_r1s_r2s[1]), exit=True)

    if args.dry_run or args.verbose:
        info('inputs_r1s: {}\n'.format('\n            '.join(inputs_r1s_r2s[0])), init_new_line=True)
        info('inputs_r2s: {}\n'.format('\n            '.join(inputs_r1s_r2s[1])))

    merged_r1_r2 = concatenate_reads(args.input_dir, inputs_r1s_r2s, nproc=args.nproc, dry_run=args.dry_run, verbose=args.verbose)

    if args.dry_run or args.verbose:
        info('merged_r1: {}\n'.format(merged_r1_r2[0]), init_new_line=True)
        info('merged_r2: {}\n'.format(merged_r1_r2[1]))

    qced_r1_r2 = quality_control(args.input_dir, merged_r1_r2, args.keep_intermediate, args.samplename,
                                 nproc=args.nproc, dry_run=args.dry_run, verbose=args.verbose)
    remove(merged_r1_r2, args.keep_intermediate, folder=args.input_dir, dry_run=args.dry_run, verbose=args.verbose)

    if args.dry_run or args.verbose:
        info('qced_r1: {}\n'.format(qced_r1_r2[0]), init_new_line=True)
        info('qced_r2: {}\n'.format(qced_r1_r2[1]))

    screened_r1_r2 = screen_contaminating_dnas(args.input_dir, qced_r1_r2, args.bowtie2_indexes, args.keep_intermediate,
                                               args.rm_hsap, args.rm_rrna, args.rm_mmus,
                                               nprocs_bowtie2=args.nproc_bowtie2 if args.nproc_bowtie2 > args.nproc else args.nproc,
                                               dry_run=args.dry_run, verbose=args.verbose)
    remove(qced_r1_r2, args.keep_intermediate, folder=args.input_dir, dry_run=args.dry_run, verbose=args.verbose)

    if args.dry_run or args.verbose:
        info('screened_r1: {}\n'.format(screened_r1_r2[0]), init_new_line=True)
        info('screened_r2: {}\n'.format(screened_r1_r2[1]))

    splitted_and_sorted = split_and_sort(args.input_dir, screened_r1_r2, args.keep_intermediate,
                                         nproc=args.nproc, dry_run=args.dry_run, verbose=args.verbose)
    remove(screened_r1_r2, args.keep_intermediate, folder=args.input_dir, dry_run=args.dry_run, verbose=args.verbose)

    if args.dry_run or args.verbose:
        info('splitted_and_sorted: {}\n'.format(splitted_and_sorted[0]), init_new_line=True)
        info('                     {}\n'.format(splitted_and_sorted[1]))
        info('                     {}\n'.format(splitted_and_sorted[2]))

    if args.verbose:
        info('time elapsed: {} s\n'.format(int(time.time() - t0)), init_new_line=True)

    sys.exit(0)
