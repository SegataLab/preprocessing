#export PATH=/shares/CIBIO-Storage/CM/mir/tools/bin:$PATH
export PATH=/shares/CIBIO-Storage/CM/scratch/users/f.asnicar/git/preprocessing/preprocessing:$PATH  # split_and_sort.py & cat_stats.py & fna_len.py
#export PATH=/shares/CIBIO-Storage/CM/scratch/users/f.asnicar/git/pyphlan:$PATH  # fna_len.py

python /shares/CIBIO-Storage/CM/scratch/users/f.asnicar/git/preprocessing/preprocessing/preprocess.new2.py $*
