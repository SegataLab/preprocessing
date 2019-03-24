#!/bin/bash


# safely copying metagenomes to be preprocessed

w='/shares/CIBIO-Storage/CM/tmp/fasnicar_preproc'  # working directory
d='run24_pool2'  # name of folder in working directory
p='_pool2_'  # pool to be considered from the seq_internal folder

# assumes to be inside the run folder in "/shares/CIBIO-Storage/CM/mir/seq_internal"
# assumes that the run folder contains files and not folders

for i in $(ls -d *$p*.fastq.gz | cut -f-4 -d'_' | sort | uniq); do
    if [ ! -d $w/$d/$i ]; then
        echo "mkdir -p $w/$d/$i"
        mkdir -p $w/$d/$i;
    fi;
    for j in $(ls $i*.fastq.gz); do
        if [ ! -f $w/$d/$i/$j ]; then
            echo "cp -an $j $w/$d/$i/"
            cp -an $j $w/$d/$i/;
        else
            f=`echo $j | rev | cut -f2- -d'_' | rev`;
            n=`echo $j | rev | cut -f1 -d'_' | rev | cut -f1 -d'.'`;
            e=`echo $j | rev | cut -f1 -d'_' | rev | cut -f2- -d'.'`;
            printf -v m "%03d" $((n+1));
            echo "cp -an $j $w/$d/$i/${f}_${m}.${e}"
            cp -an $j $w/$d/$i/${f}_${m}.${e};
        fi;
    done;
done;
