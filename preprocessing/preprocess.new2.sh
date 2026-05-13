export PATH=/shares/CIBIO-Storage/CM/scratch/users/f.asnicar/git/preprocessing/preprocessing:$PATH  # cat_stats.py & fna_len.py
export PATH=/shares/CIBIO-Storage/CM/scratch/users/f.asnicar/00_scripts/bbmap:$PATH  # repair.sh

python /shares/CIBIO-Storage/CM/scratch/users/f.asnicar/git/preprocessing/preprocessing/preprocess.new2.py $*
