#!/bin/bash

# Function to calculate statistics
calculate_stats() {
    local FASTQ_FILE=$1
    local OUTPUT_FILE=$2

    # Initialize variables
    total_bases=0
    total_reads=0
    min_length=-1
    max_length=0
    lengths_file=$(mktemp)

    # Read the FASTQ file
    gawk 'NR%4==2 {
        read_length=length($0);
        total_bases+=read_length;
        total_reads++;
        if (min_length == -1 || read_length < min_length) min_length=read_length;
        if (read_length > max_length) max_length=read_length;
        print read_length >> "'$lengths_file'"
    }
    END {
        if (total_reads > 0) {
            mean_length=total_bases/total_reads;
        } else {
            min_length=0;
            max_length=0;
            mean_length=0;
        }
        print total_bases, total_reads, min_length, mean_length, max_length;
    }' min_length=-1 $FASTQ_FILE > ${FASTQ_FILE}_stats.txt

    # Extract values from temp file
    read total_bases total_reads min_length mean_length max_length < ${FASTQ_FILE}_stats.txt

    # Calculate median read length
    median_length=$(sort -n $lengths_file | awk '{a[i++]=$1} END {if (i%2==0) print (a[i/2-1]+a[i/2])/2; else print a[i/2]}')

    # Clean up temp files
    rm ${FASTQ_FILE}_stats.txt
    rm $lengths_file

    # Print results to the output file
    echo -e "#samplename\t n_of_bases\t n_of_reads\t min_read_len\t median_read_len\t mean_read_len\t max_read_len" >> "$OUTPUT_FILE"
    echo -e "$(basename $FASTQ_FILE)\t $total_bases\t $total_reads\t $min_length\t $median_length\t $mean_length\t $max_length" >> "$OUTPUT_FILE"
}

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ] && [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_fastq1> [<input_fastq2>] <output_file1> [<output_file2>]"
    exit 1
fi

# Single-end or paired-end
if [ "$#" -eq 2 ]; then
    INPUT_FASTQ1=$1
    OUTPUT_FILE1=$2

    # Calculate stats for single-end FASTQ
    calculate_stats "$INPUT_FASTQ1" "$OUTPUT_FILE1"

elif [ "$#" -eq 4 ]; then
    INPUT_FASTQ1=$1
    INPUT_FASTQ2=$2
    OUTPUT_FILE1=$3
    OUTPUT_FILE2=$4

    # Calculate stats for paired-end FASTQ
    calculate_stats "$INPUT_FASTQ1" "$OUTPUT_FILE1"
    calculate_stats "$INPUT_FASTQ2" "$OUTPUT_FILE2"
fi
