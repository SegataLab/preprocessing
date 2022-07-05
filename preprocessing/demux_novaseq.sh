# conda deactivate && conda activate fasnicar_bcl2fastq; 

ulimit -n 999999; 
bcl2fastq --barcode-mismatches 1 -r 16 -w 16 -p 32 -R RUN_FOLDER --sample-sheet SAMPLE_SHEET.CSV -o OUTPUT_FOLDER

#rm -r Reports/ Stats/ Undetermined_S0_L00*.fastq.gz
#for i in $(ls *.fastq.gz | rev | cut -f5- -d'_' | rev | sort -r | uniq); do mkdir ${i}; echo $i; mv ${i}_*.fastq.gz ${i}/; done
