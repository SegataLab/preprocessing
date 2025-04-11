#!/bin/bash


# conda deactivate && conda activate fastqc
# parallel -j10 "../fastqc.sh {}" ::: `ls`

echo $1;

mkdir -p $1/fastqc_raw;

for i in $(ls $1/*.fastq.gz); do
    echo "    $i";
    gzip -dk $i;
    fastqc $i -o $1/fastqc_raw/ -t 8 -f fastq --extract -q;
    rm ${i%.gz};
done;

# conda deactivate && conda activate
