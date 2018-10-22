#!/bin/bash


# binaries="preprocess.py split_and_sort.py cat_stats.py fna_len.py"
binaries="preprocess.new.py split_and_sort.new.py cat_stats.py fna_len.py"

mkdir -p $PREFIX/bin

for i in $binaries; do
    cp preprocessing/$i $PREFIX/bin;
done
