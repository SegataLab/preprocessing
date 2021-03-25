# conda deactivate && conda activate fasnicar_bcl2fastq; 

ulimit -n 999999; 
bcl2fastq --barcode-mismatches 1 -r 16 -w 16 -p 32 -R RUN_FOLDER --sample-sheet SAMPLE_SHEET.CSV -o OUTPUT_FOLDER

#sudo cp -a 201103_A01083_0039_AHYNF2DSXX/ ../seq_internal/NovaSeq_fastqs/
#rm -r Reports/ Stats/ Undetermined_S0_L00*.fastq.gz
#for i in $(ls | rev | cut -f5- -d'_' | rev | sort | uniq); do mkdir ${i}; echo $i; mv ${i}_*.fastq.gz ${i}/; done
