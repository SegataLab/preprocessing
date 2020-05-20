# conda deactivate && conda activate fasnicar_bcl2fastq; 

ulimit -n 999999; 
bcl2fastq --barcode-mismatches 1 -r 16 -w 16 -p 32 -R RUN_FOLDER --sample-sheet SAMPLE_SHEET.CSV -o OUTPUT_FOLDER
