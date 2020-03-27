#!/bin/bash


# conda deactivate && conda activate fastqc

echo $1;

for i in $(ls $1/*.fastq.gz); do
    echo "    $i";
    gzip -dk $i;
    fastqc $i -o $1/fastqc_raw/ -t 8 -f fastq -q;
    rm ${i%.gz};
done;

# conda deactivate && conda activate
