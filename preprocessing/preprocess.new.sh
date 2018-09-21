export PATH=/shares/CIBIO-Storage/CM/mir/tools/bin:$PATH
export PATH=/shares/CIBIO-Storage/CM/scratch/users/f.asnicar/hg/preprocessing:$PATH  # split_and_sort.py & cat_stats.py
export PATH=/shares/CIBIO-Storage/CM/scratch/users/f.asnicar/hg/pyphlan:$PATH  # fna_len.py

python preprocess.new.py $*
