#!/usr/bin/env python

from Bio import SeqIO
import argparse
import bz2
import itertools
import multiprocessing as mp
import logging

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

def index_fastq(file_path):
    logger.info(f"Indexing FASTQ file: {file_path}")
    return SeqIO.index(file_path, "fastq")

def write_fastq_bz2(task):
    read_iterator, output_file = task
    logger.info(f"Writing {output_file}")
    with bz2.open(output_file, 'wt') as f:
        SeqIO.write(read_iterator, f, "fastq")
    logger.info(f"Finished writing")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--R1", help="R1 input file", required=True)
    parser.add_argument("--R2", help="R2 input file", required=True)
    parser.add_argument("-p", '--prefix', help="Output prefix", default='out')
    parser.add_argument("-u", '--unpaired', help="File with unpaired read IDs")
    args = parser.parse_args()

    logger.info("Starting split_and_sort")
    logger.info(f"Arguments: {args}")

    R1_index = index_fastq(args.R1)
    R2_index = index_fastq(args.R2)

    logger.info("Creating read iterators")
    paired_r1_id_sorted = sorted((r1_id for r1_id in R1_index if r1_id in R2_index))
    paired_r2_id_sorted = sorted((r2_id for r2_id in R2_index if r2_id in R1_index))
    paired_r1_iterator = (R1_index[r1_id] for r1_id in paired_r1_id_sorted)
    paired_r2_iterator = (R2_index[r2_id] for r2_id in paired_r2_id_sorted)

    tasks = [(paired_r1_iterator, args.prefix + '_R1.fastq.bz2'),
             (paired_r2_iterator, args.prefix + '_R2.fastq.bz2')]

    if args.unpaired:
        logger.info(f"Reading unpaired reads")
        unpaired_reads = set(line.strip() for line in bz2.open(args.unpaired, 'rt'))
        logger.info(f"Creating unpaired iterator")
        unpaired_reads_iter = itertools.chain(
            (R1_index[r1_id] for r1_id in R1_index if r1_id in unpaired_reads),
            (R2_index[r2_id] for r2_id in R2_index if r2_id in unpaired_reads)
        )
        tasks.append((unpaired_reads_iter, args.prefix + '_UN.fastq.bz2'))

    logger.info("Starting writing output files")
    for task in tasks:
        write_fastq_bz2(task)

    logger.info("split_and_sort completed successfully")

if __name__ == "__main__":
    main()
