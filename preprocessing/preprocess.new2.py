#!/usr/bin/env python

__author__ = 'Francesco Asnicar (f.asnicar@unitn.it)'
__version__ = '0.3.8'
__date__ = '18 April 2025'

import os
import sys
import bz2
import gzip
import glob
import argparse
import subprocess as sb
import multiprocessing as mp
import logging
import shutil
from Bio import SeqIO

if sys.version_info < (3, 6):
    raise Exception("Preprocessing requires Python 3.6 or higher. Your current Python version is {}.{}.{}".format(
        sys.version_info[0], sys.version_info[1], sys.version_info[2]))

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

def log_and_exit(message, code=1):
    logger.error(message)
    sys.exit(code)

def read_params():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input_dir', required=True, help="Path to input directory")
    parser.add_argument('-e', '--extension', default=".fastq.gz", choices=[".fastq.gz", ".fq.gz", ".fastq.bz2", ".fq.bz2"],
        help="Extension of the raw input files")
    parser.add_argument('-s', '--samplename', default="", help="Sample name")
    parser.add_argument('-f', '--forward', default="_R1_", help="Identifier for forward reads")
    parser.add_argument('-r', '--reverse', default="_R2_", help="Identifier for reverse reads")
    parser.add_argument('-n', '--nproc', type=int, default=2,
        help="Number of threads, used in some parts (i.e., decompressing and concatenating and Bowtie2)")
    parser.add_argument('-p', '--paired_end', action='store_true', help="Indicates paired-end sequencing")
    parser.add_argument('-k', '--keep_intermediate', action='store_true', help="Keep intermediate files")
    parser.add_argument('--dry_run', action='store_true', help="Print commands without executing them")
    parser.add_argument('--verbose', action='store_true', help="Enable verbose logging")
    parser.add_argument('-x', '--bowtie2_indexes', default='/shares/CIBIO-Storage/CM/scratch/databases/bowtie2_indexes',
        help="Path to Bowtie2 indexes")
    parser.add_argument('--rm_hsap', action='store_true', help="Remove Homo sapiens genome (hg19)")
    parser.add_argument('--rm_GRCh38', action='store_true', help="Remove Homo sapiens genome (GRCh38.p14 from GCF_000001405.40)")
    parser.add_argument('--rm_mmus', action='store_true', help="Remove Mus musculus genome")
    parser.add_argument('--rm_rrna', action='store_true', help="Remove rRNA sequences")
    parser.add_argument('--rm_pcin', action='store_true', help="Remove Phascolarctos cinereus genome")
    parser.add_argument('--rm_pcoq', action='store_true', help="Remove Propithecus coquereli genome")
    parser.add_argument('--rm_mmur', action='store_true', help="Remove Microcebus murinus genome")
    parser.add_argument('--rm_mmul', action='store_true', help="Remove Macaca mulatta genome")
    parser.add_argument('--rm_ptro', action='store_true', help="Remove Pan troglodytes genome")
    parser.add_argument('--rm_sbol', action='store_true', help="Remove Saimiri boliviensis genome")
    parser.add_argument('--rm_vvar', action='store_true', help="Remove Varecia variegata genome")
    parser.add_argument('--rm_clup', action='store_true', help="Remove Canis lupus familiaris genome")
    parser.add_argument('--rm_cjac', action='store_true', help="Remove Callithrix jacchus genome")
    parser.add_argument('--rm_agig', action='store_true', help="Remove Aldabrachelys gigantea genome")
    parser.add_argument('--rm_alho', action='store_true', help="Remove Allochrocebus lhoesti genome")
    parser.add_argument('--rm_soed', action='store_true', help="Remove Saguinus oedipus genome")
    parser.add_argument('--rm_dmel', action='store_true', help="Remove Drosophila melanogaster genome")
    return parser.parse_args()

def validate_input_dir(input_dir):
    logger.debug(f"Validating input directory: {input_dir}")
    if not os.path.isdir(input_dir):
        log_and_exit(f"Input folder '{input_dir}' does not exist.")
    return os.path.abspath(input_dir)

def execute_command(cmd, dry_run=False):
    logger.info(f"Executing: {cmd}")
    logger.debug(f"Dry run mode: {dry_run}")
    if dry_run:
        logger.info("Dry run enabled. Command not executed.")
        return
    try:
        sb.check_call(cmd, shell=True, stderr=sys.stderr)
    except sb.CalledProcessError as e:
        log_and_exit(f"Command failed: {cmd}\nError: {str(e)}")

def remove(to_remove, keep_intermediate, folder=None, dry_run=False):
    logger.info("Executing remove function.")
    if keep_intermediate:
        logger.info("Skipping removal of intermediate files (keep_intermediate=True).")
        return

    for file in to_remove:
        file_path = os.path.join(folder, file) if folder else file

        if os.path.isfile(file_path):
            logger.info(f"Removing file: {file_path}")
            if not dry_run:
                try:
                    os.remove(file_path)
                except Exception as e:
                    logger.error(f"Failed to remove {file_path}: {str(e)}")
        else:
            logger.warning(f"File not found: {file_path}")

def run_fna_len(input_files, output_files, dry_run=False, nproc=1):
    if len(input_files) != len(output_files):
        log_and_exit("Error: input_files and output_files lists must have the same length.")

    logger.info(f"Running fna_len.py on {','.join(input_files)}")
    tasks = list(zip(input_files, output_files, [dry_run] * len(input_files)))

    try:
        with mp.Pool(processes=nproc) as pool:
            pool.map(process_fna_len, tasks)
    except OSError as e:
        log_and_exit(f"fna_len failed: {str(e)}")

    logger.debug(f"Completed fna_len.py processing.")

def process_fna_len(task):
    input_file, output_file, dry_run = task
    logger.info(f"Processing {input_file} -> {output_file}")
    cmd = f"fna_len2.py -q --stat {input_file} > {output_file}"
    execute_command(cmd, dry_run)

def run_cat_stats(output_dir, stat_files, final_output_file, dry_run=False):
    logger.info(f"Running cat_stats.py to consolidate statistics.")
    input_files = " ".join(stat_files)
    cmd = f"cat_stats2.py -i {input_files} -o {final_output_file}"
    execute_command(cmd, dry_run)

def get_files(input_dir, extension, forward, reverse):
    logger.info("Retrieving input files.")
    logger.debug(f"Input directory: {input_dir}, Extension: {extension}, Forward: {forward}, Reverse: {reverse}")
    forward_files = sorted(glob.glob(os.path.join(input_dir, f"*{forward}*{extension}")))
    reverse_files = sorted(glob.glob(os.path.join(input_dir, f"*{reverse}*{extension}")))
    if not forward_files or not reverse_files:
        log_and_exit("No forward or reverse read files found in the input directory.")
    logger.debug(f"Forward files: {forward_files}, Reverse files: {reverse_files}")
    return forward_files, reverse_files

def concatenate_reads(forward_file_list, reverse_file_list, input_dir, dry_run=False, nproc=1):
    logger.info("Concatenating reads in parallel.")

    forward_output_file = os.path.join(input_dir, "concatenated_R1.fastq")
    reverse_output_file = os.path.join(input_dir, "concatenated_R2.fastq")

    logger.debug(f"Forward file list: {forward_file_list}, Output file: {forward_output_file}")
    logger.debug(f"Reverse file list: {reverse_file_list}, Output file: {reverse_output_file}")
    logger.debug(f"Input directory: {input_dir}, dry run: {dry_run}, nproc: {nproc}")

    outputs = []
    if os.path.isfile(forward_output_file) and os.path.isfile(reverse_output_file):
        logger.info("Both output files already exist. Skipping concatenation.")
        return [forward_output_file, reverse_output_file]

    if dry_run:
        logger.info("Dry run enabled. Simulated concatenation.")
        return [forward_output_file, reverse_output_file]

    tasks = [(forward_file_list, forward_output_file), (reverse_file_list, reverse_output_file)]

    try:
        with mp.Pool(processes=nproc) as pool:
            outputs = pool.map(decompress_and_concatenate, tasks)

        logger.debug(f"Concatenated output files: {outputs}")
    except OSError as e:
        log_and_exit(f"Failed to concatenate reads: {str(e)}")

    return outputs

def decompress_and_concatenate(task):
    file_list, output_file = task
    logger.info(f"Concatenating {','.join(file_list)} into {output_file}")

    try:
        with open(output_file, 'wb') as g:
            for file in file_list:
                logger.debug(f"Processing: {file}")
                with (gzip.open if file.endswith('.gz') else bz2.open if file.endswith('.bz2') else open)(file, 'rb') as f:
                    shutil.copyfileobj(f, g)

        return output_file
    except Exception as e:
        log_and_exit(f"Failed to concatenate {output_file}: {str(e)}")

# DO WE WANT TO DO PAIRED-END TRIMMING?
def quality_control(input_files, output_dir, dry_run=False, nproc=1):
    logger.info(f"Performing quality control on {','.join(input_files)}")
    tasks = [(input_file, output_dir, dry_run) for input_file in input_files]

    try:
        with mp.Pool(processes=nproc) as pool:
            trimmed_files = pool.map(process_quality_control, tasks)
    except OSError as e:
        log_and_exit(f"quality_control failed: {str(e)}")

    logger.debug(f"Completed quality control.")
    return trimmed_files

def process_quality_control(task):
    input_file, output_dir, dry_run = task
    logger.info(f"Processing {input_file}")
    cmd = (
        "trim_galore --nextera --length 75 --2colour 20 --max_n 2 --trim-n -j 1 "
        f"--no_report_file --suppress_warn --output_dir {output_dir} {input_file}"
    )
    execute_command(cmd, dry_run)
    trimmed_file = os.path.join(output_dir, os.path.basename(input_file).replace('.fastq', '_trimmed.fq'))
    if not os.path.isfile(trimmed_file):
        log_and_exit(f"Quality control failed for {input_file}.")
    return trimmed_file

def get_unpaired(input_dir, r1, r2, samplename):
    logger.info("Identifying unpaired reads.")
    logger.debug(f"Input R1: {r1}, Input R2: {r2}")
    
    r1_path = os.path.join(input_dir, r1)
    r2_path = os.path.join(input_dir, r2)
    unpaired_file = os.path.join(input_dir, f"{samplename}_unpaired.txt.bz2")

    r1_index = SeqIO.index(r1_path, "fastq")
    r2_index = SeqIO.index(r2_path, "fastq")

    unpaired_r1 = (i for i in r1_index if i not in r2_index)
    unpaired_r2 = (i for i in r2_index if i not in r1_index)

    with bz2.open(unpaired_file, 'wt') as f:
        f.write('\n'.join(unpaired_r1) + '\n' + '\n'.join(unpaired_r2) + '\n')

    logger.info(f"Unpaired reads written to: {unpaired_file}")
    return unpaired_file

def remove_contaminants(input_file, bowtie2_index_dir, output_dir, contaminant_list, nproc=1, dry_run=False):
    logger.info("Removing contaminants.")
    current_input = input_file
    to_remove = []
    fna_len_lists = {'inp': [], 'out': []}
    rx = 'R1' if '_R1' in input_file else 'R2' if '_R2' in input_file else None  # MAYBE WE SHOULD USE THE FWD/REV PARAMS?

    if not rx:
        log_and_exit(f"Failed to identify forward (R1) or reverse (R2) from input file {input_file}.")

    for contaminant in contaminant_list:
        logger.info(f"Removing contaminant: {contaminant}")
        output_file = os.path.join(output_dir, os.path.basename(current_input).replace('.fq', f'_{contaminant}.fq'))

        cmd = (
            f"bowtie2 -x {os.path.join(bowtie2_index_dir, contaminant)} "
            f"-U {current_input} -p {nproc} --sensitive-local "
            f"--un {output_file} > /dev/null 2>&1"
        )

        execute_command(cmd, dry_run)

        if not os.path.isfile(output_file):
            log_and_exit(f"Contaminant removal failed for {current_input} with contaminant {contaminant}.")

        fna_len_lists['inp'].append(output_file)
        fna_len_lists['out'].append(os.path.join(input_dir, f'stats_screened_{rx}_{contaminant}.txt'))
        to_remove.append(current_input)
        current_input = output_file

    return current_input, fna_len_lists, to_remove

def split_and_sort(input_dir, r1, r2, unpaired_file, samplename, dry_run=False):
    logger.info("Splitting and sorting reads.")
    output_r1 = os.path.join(input_dir, f"{samplename}_R1.fastq.bz2")
    output_r2 = os.path.join(input_dir, f"{samplename}_R2.fastq.bz2")
    output_un = os.path.join(input_dir, f"{samplename}_UN.fastq.bz2")
    output_unpaired = os.path.join(input_dir, unpaired_file)

    if all(os.path.isfile(f) for f in [output_r1, output_r2, output_un]):
        logger.info("Split and sort output files already exist. Skipping step.")
        return output_r1, output_r2, output_un

    cmd = f"split_and_sort.new2.py --R1 {r1} --R2 {r2} --prefix {os.path.join(input_dir, samplename)}"
    cmd += f" --unpaired {output_unpaired}" if unpaired_file else ""
    execute_command(cmd, dry_run)

    return output_r1, output_r2, output_un

if __name__ == "__main__":
    args = read_params()

    if args.verbose:
        logger.setLevel(logging.DEBUG)
        logger.info(f'Preprocessing version {__version__} ({__date__})')
        logger.info(f"Command line: {' '.join(sys.argv)}")

    logger.info("Validating input directory.")
    input_dir = validate_input_dir(args.input_dir)

    logger.info("Starting preprocessing pipeline.")

    forward_files, reverse_files = get_files(input_dir, args.extension, args.forward, args.reverse)
    logger.info(f"Found {len(forward_files)} forward files and {len(reverse_files)} reverse files.")

    concatenated_forward, concatenated_reverse = concatenate_reads(forward_files, reverse_files, input_dir, args.dry_run, args.nproc)
    stats_concatenated_R1 = os.path.join(input_dir, "stats_concatenated_R1.txt")
    stats_concatenated_R2 = os.path.join(input_dir, "stats_concatenated_R2.txt")
    stats_forward = [stats_concatenated_R1]
    stats_reverse = [stats_concatenated_R2]

    trimmed_forward, trimmed_reverse = quality_control([concatenated_forward, concatenated_reverse], input_dir, args.dry_run, args.nproc)
    stats_trimmed_R1 = os.path.join(input_dir, "stats_trimmed_R1.txt")
    stats_trimmed_R2 = os.path.join(input_dir, "stats_trimmed_R2.txt")
    stats_forward.append(stats_trimmed_R1)
    stats_reverse.append(stats_trimmed_R2)

    run_fna_len([concatenated_forward, concatenated_reverse, trimmed_forward, trimmed_reverse],
                [stats_concatenated_R1, stats_concatenated_R2, stats_trimmed_R1, stats_trimmed_R2],
                args.dry_run, args.nproc)
    remove([concatenated_forward, concatenated_reverse], args.keep_intermediate, input_dir, args.dry_run)

    unpaired_file = get_unpaired(input_dir, trimmed_forward, trimmed_reverse, args.samplename) if args.paired_end else None

    contaminants = ['phiX174']
    contaminants.append("hg19") if args.rm_hsap else None
    contaminants.append("GRCh38p14") if args.rm_GRCh38 else None
    contaminants.append("mmusculus_black6_GCA_000001635_8") if args.rm_mmus else None
    contaminants.extend(["SILVA_132_SSURef_Nr99_tax_silva", "SILVA_132_LSURef_tax_silva"]) if args.rm_rrna else None
    contaminants.append("Phascolarctos_cinereus__GCA_900166895.1__tgac_v2.0") if args.rm_pcin else None
    contaminants.append("Propithecus_coquereli_GCF_000956105.1") if args.rm_pcoq else None
    contaminants.append("Microcebus_murinus_GCF_000165445.2") if args.rm_mmur else None
    contaminants.append("mmulatta_GCF_003339765.1") if args.rm_mmul else None
    contaminants.append("ptroglodytes_GCF_028858775.1") if args.rm_ptro else None
    contaminants.append("sboliviensis_GCF_016699345.2") if args.rm_sbol else None
    contaminants.append("vvariegata_GCA_028533085.1") if args.rm_vvar else None
    contaminants.append("clfamiliaris_GCF_000002285.5") if args.rm_clup else None
    contaminants.append("cjacchus_GCF_011100555.1") if args.rm_cjac else None
    contaminants.append("agigantea_GCA_026122505.1") if args.rm_agig else None
    contaminants.append("GCA_963574325") if args.rm_alho else None
    contaminants.append("GCA_031835075") if args.rm_soed else None
    contaminants.append("Drosophila_melanogaster_GCF_000001215") if args.rm_dmel else None

    logger.info("Screening for contaminants.")
    screened_forward, fna_len_dict_forward, to_remove_forward = remove_contaminants(trimmed_forward, args.bowtie2_indexes,
                                                                                    input_dir, contaminants, args.nproc, args.dry_run)
    stats_forward.extend(fna_len_dict_forward['out'])

    screened_reverse, fna_len_dict_reverse, to_remove_reverse = remove_contaminants(trimmed_reverse, args.bowtie2_indexes,
                                                                                    input_dir, contaminants, args.nproc, args.dry_run)
    stats_reverse.extend(fna_len_dict_reverse['out'])

    run_fna_len(fna_len_dict_forward['inp'] + fna_len_dict_reverse['inp'],
                fna_len_dict_forward['out'] + fna_len_dict_reverse['out'],
                args.dry_run, args.nproc)
    remove(to_remove_reverse + to_remove_forward, args.keep_intermediate, input_dir, args.dry_run)

    split_outputs = split_and_sort(input_dir, screened_forward, screened_reverse, unpaired_file, args.samplename, args.dry_run)
    stats_split_R1 = os.path.join(input_dir, "stats_split_R1.txt")
    stats_split_R2 = os.path.join(input_dir, "stats_split_R2.txt")
    stats_split_UN = os.path.join(input_dir, "stats_split_UN.txt")
    run_fna_len(split_outputs, [stats_split_R1, stats_split_R2, stats_split_UN], args.dry_run, args.nproc)
    remove([screened_forward, screened_reverse, unpaired_file], args.keep_intermediate, input_dir, args.dry_run)

    logger.info("Summarize statistics using cat_stats.py.")
    stat_files = stats_forward + stats_reverse + [stats_split_R1, stats_split_R2, stats_split_UN]
    summary_stats_file = os.path.join(input_dir, f"{args.samplename}_summary.stats")
    run_cat_stats(input_dir, stat_files, summary_stats_file, args.dry_run)

    remove(stat_files, args.keep_intermediate, input_dir, args.dry_run)

    logger.info("Preprocessing pipeline completed successfully.")
