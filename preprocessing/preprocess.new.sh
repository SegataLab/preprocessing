export PATH=/shares/CIBIO-Storage/CM/mir/tools/bin:$PATH
export PATH=/shares/CIBIO-Storage/CM/news/users/f.asnicar/git/preprocessing/preprocessing:$PATH  # split_and_sort.py & cat_stats.py
export PATH=/shares/CIBIO-Storage/CM/news/users/f.asnicar/git/pyphlan:$PATH  # fna_len.py

python /shares/CIBIO-Storage/CM/news/users/f.asnicar/git/preprocessing/preprocessing/preprocess.new.py $*
