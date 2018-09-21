## How to use it


Example command line:


```
$ parallel -j NCPU 'preprocess.sh -i {} [other params]' ::: `ls input_folder`
```


Where:
 * `preprocess.sh` takes one parameter which is the input folder containing the raw reads
 * input folder content should look like:

    > ...R1_L001.fastq.gz
      ...R2_L001.fastq.gz
      ...
      ...R1_L008.fastq.gz
      ...R2_L008.fastq.gz
