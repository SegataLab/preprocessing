#!/bin/bash


if [ $# -ne 2 ]; then
    echo "No arguments";
    echo "Usage:";
    echo "./make_sample_list.sh <output_file> <input_folder>";
    exit 1;
fi;

if [[ -z $1 ]]; then
    echo "Output file parameter is empty";
    exit 1;
fi;

if [[ -z $2 ]]; then
    echo "Input folder parameter is empty";
    exit 1;
fi;

if [ -f $1 ]; then
    echo "Output file already exists";
    exit 1;
fi;

if [ ! -d $2 ]; then
    echo "Input folder does not exists";
    exit 1;
fi;

echo -n "Writing output file: $1 ... ";

echo -e "#sample\tn_cleaned_reads" > $1;
tail -n1 $2/*/*_summary.stats | grep -v "^$" | grep -v "^==>" | cut -f1,3 >> $1;

echo "Done!";
