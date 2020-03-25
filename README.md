## How to install it

```
conda install preprocessing -c fasnicar
```


## How to use it

Example command line:

```
parallel -j NCPU 'preprocess.sh -i {} [other params]' ::: `ls input_folder`
```

where:

- `preprocess.sh` takes one parameter which is the input folder containing the raw reads
- input folder should contains the raw reads


## Building conda package

```
conda build recipe --python 3 -c bioconda
```
