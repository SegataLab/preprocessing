#!/bin/bash

# Function to display help message
show_help() {
    echo "Usage: $0 -f <folder_path> -m <mapper_type> -r <reference_genome> -p <threads> -o <output_dir> [-h]"
    echo
    echo "Required arguments:"
    echo "  -f  Sample Folder path containing the input fastq files the sample."
    echo "  -m  Mapper type ('bowtie' or 'kraken'). Default is 'bowtie'."
    echo "  -r  Reference genome (e.g., 'hg19')."
    echo "  -p  Number of threads to use. Default is 2."
    echo
    echo "Optional arguments:"
    echo "  -o  Output directory where results will be saved. Defaults to the input folder path."
    echo "  -h  Show this help message and exit."
}

# Activate the Conda environment
. /shares/CIBIO-Storage/CM/scratch/prebiomics/tools/anaconda3/etc/profile.d/conda.sh
conda activate preprocessing

# Default values
mapper='bowtie'
threads=2
output_dir=""

# Parse command-line arguments
while getopts ":f:m:r:p:o:h" opt; do
    case ${opt} in
        f )
            FOLDER_PATH=$OPTARG
            ;;
        m )
            mapper=$OPTARG
            ;;
        r )
            ref=$OPTARG
            ;;
        p )
            threads=$OPTARG
            ;;
        o )
            output_dir=$OPTARG
            ;;
        h )
            show_help
            exit 0
            ;;
        \? )
            echo "Invalid option: -$OPTARG" >&2
            show_help
            exit 1
            ;;
        : )
            echo "Option -$OPTARG requires an argument." >&2
            show_help
            exit 1
            ;;
    esac
done

# Check for required arguments
if [ -z "$FOLDER_PATH" ] || [ -z "$ref" ]; then
    echo "Error: Folder path and reference genome are required."
    show_help
    exit 1
fi

# If output directory is not specified, use the input folder path
if [ -z "$output_dir" ]; then
    output_dir="$FOLDER_PATH"
fi

# Validate mapper type
if [[ "$mapper" != "bowtie" && "$mapper" != "kraken" ]]; then
    echo "Error: <mapper_type> must be either 'bowtie' or 'kraken'"
    exit 1
fi

# Check for the presence of .bz2 files in the input folder
if find "$FOLDER_PATH" -type f -name '*.bz2' | grep -q '.'; then
    echo "Error: It seems the preprocessing has already been done for this sample!"
    exit 1
fi

# Check for the presence of fastq files in the input folder
if ! find "$FOLDER_PATH" -type f -name '*fastq*' | grep -q '.'; then
    echo "Error: There is no fastq file inside the sample folder provided!"
    exit 1
fi

# Extract the sample name from the folder name
SAMPLE_NAME=$(basename "$FOLDER_PATH")

# Find all forward (R1) and reverse (R2) files
R1_FILES=($(find "$FOLDER_PATH" -type f -name '*R1*fastq.gz'))
R2_FILES=($(find "$FOLDER_PATH" -type f -name '*R2*fastq.gz'))

####################################################################### Concatenating

echo "Concatenating started!"
start_time=$(date +%s)
# Function to concatenate files
concatenate_files() {
    local FILE_ARRAY=("${!1}")
    local SAMPLE_NAME=$2
    local FILE_SUFFIX=$3
    local OUTPUT_DIR=$4

    if [ "${#FILE_ARRAY[@]}" -eq 1 ]; then
        # If there is only one file, just decompress it
        pigz -d -p 8 -c "${FILE_ARRAY[0]}" > "${OUTPUT_DIR}/${SAMPLE_NAME}_${FILE_SUFFIX}.fq"
    else
        # If there are multiple files, decompress and concatenate them
        for FILE in "${FILE_ARRAY[@]}"; do
            pigz -d -p 8 -c "$FILE" >> "${OUTPUT_DIR}/${SAMPLE_NAME}_${FILE_SUFFIX}.fq"
        done
    fi
}

export -f concatenate_files

# Run the concatenation in parallel using background jobs
concatenate_files R1_FILES[@] "$SAMPLE_NAME" "R1" "$output_dir" &
pid1=$!
concatenate_files R2_FILES[@] "$SAMPLE_NAME" "R2" "$output_dir" &
pid2=$!

# Wait for both background jobs to complete
wait $pid1
wait $pid2

end_time=$(date +%s)
execution_time=$((end_time - start_time))
echo "Concatenating ended! (time of execution: ${execution_time}s)"

echo

######################################################################### Quality control

echo "Quality control started!"
start_time=$(date +%s)

values=("${output_dir}/${SAMPLE_NAME}_R1.fq" "${output_dir}/${SAMPLE_NAME}_R2.fq")
quality_control() {
    local input=$1
    trim_galore --nextera --stringency 5 --length 75 --2colour 20 --max_n 2 --trim-n -j 1 --dont_gzip \
    --no_report_file --suppress_warn \
    --output_dir $2 ${input} > /dev/null 2>&1
}

export -f quality_control

# Use GNU Parallel to run the process_value function in parallel
parallel quality_control ::: ${values[@]} ::: ${output_dir}

end_time=$(date +%s)
execution_time=$((end_time - start_time))
echo "Quality control ended! (time of execution: ${execution_time}s)"

echo

######################################################################### Screen_contaminating_dnas

echo "Screen contaminating dnas started!"
start_time=$(date +%s)

ref_bowtie_address='/shares/CIBIO-Storage/CM/scratch/databases/bowtie2_indexes'
ref_kraken_address='/shares/CIBIO-Storage/CM/scratch/databases/kraken2_databases'
first_bowtie_genome=${ref_bowtie_address}/phiX174
second_bowtie_genome=${ref_bowtie_address}/${ref}
first_kraken_genome=${ref_kraken_address}/phiX174
second_kraken_genome=${ref_kraken_address}/${ref}

values=("${output_dir}/${SAMPLE_NAME}_R1_trimmed.fq" "${output_dir}/${SAMPLE_NAME}_R2_trimmed.fq")

mapping() {
    local input=$1
    local first_bowtie_genome=$2
    local second_bowtie_genome=$3
    local first_kraken_genome=$4
    local second_kraken_genome=$5
    local mapper=$6
    local ref=$7

    output="${input%.fq}_phiX174.fq"

    if [ "$mapper" == "bowtie" ]; then
        bowtie2 -x "$first_bowtie_genome"  -U "$input"  -p "$threads"  --sensitive-local  --un "${output}" > /dev/null 2>&1
    elif [ "$mapper" == "kraken" ]; then
        kraken2 --db "$first_kraken_genome" --confidence 0.05 --minimum-hit-groups 3 --threads "$threads" --unclassified-out "${output}" "$input" > /dev/null 2>&1
    fi

    input="${output}"
    output="${input%.fq}_$ref.fq"

    if [ "$mapper" == "bowtie" ]; then
        bowtie2 -x "$second_bowtie_genome" -U "$input" -p "$threads" --sensitive-local --un "${output}" > /dev/null 2>&1
    elif [ "$mapper" == "kraken" ]; then
        kraken2 --db "$second_kraken_genome" --confidence 0.05 --minimum-hit-groups 3 --threads "$threads" --unclassified-out "${output}" "$input" > /dev/null 2>&1
    fi
}

export -f mapping

# Use GNU Parallel to run the process_value function in parallel
parallel mapping ::: "${values[@]}" ::: "$first_bowtie_genome" ::: "$second_bowtie_genome" ::: "$first_kraken_genome" ::: "$second_kraken_genome" ::: "$mapper" ::: "$ref"

end_time=$(date +%s)
execution_time=$((end_time - start_time))
echo "Screen contaminating dnas ended! (time of execution: ${execution_time}s)"

echo

############################################################################################ Getting unpaired indices

echo "Getting unpaired indices analysis started!"
start_time=$(date +%s)

values=("${output_dir}/${SAMPLE_NAME}_R1_trimmed.fq" "${output_dir}/${SAMPLE_NAME}_R2_trimmed.fq")
indexing() {
    local input=$1
    filename="${input##*/}"
    filename_without_extension="${filename%.fq}"
    name="index_${filename_without_extension}"

    seqkit seq -n ${input} | cut -d' ' -f1 | sort > $2/${name}
}

export -f indexing

# Use GNU Parallel to run the process_value function in parallel
parallel indexing ::: ${values[@]} ::: ${output_dir}

# Find unique identifiers (only in one of the files)
comm -3 ${output_dir}/index_${SAMPLE_NAME}_R1_trimmed  ${output_dir}/index_${SAMPLE_NAME}_R2_trimmed | gawk '{print $1}' > ${output_dir}/${SAMPLE_NAME}_unpaired.txt

# Delete the intermediate files
rm ${output_dir}/index_${SAMPLE_NAME}_R1_trimmed ${output_dir}/index_${SAMPLE_NAME}_R2_trimmed

end_time=$(date +%s)
execution_time=$((end_time - start_time))
echo "Getting unpaired indices analysis ended! (time of execution: ${execution_time}s)"

echo

############################################################################################# Split and Sorting

echo "Split and sorting started!"
start_time=$(date +%s)

values=("${FOLDER_PATH}/${SAMPLE_NAME}_R1_trimmed_phiX174_${ref}.fq" "${FOLDER_PATH}/${SAMPLE_NAME}_R2_trimmed_phiX174_${ref}.fq")
indexing_step2() {
    local input=$1
    filename="${input##*/}"
    filename_without_extension="${filename%.fq}"
    name="index_${filename_without_extension}"
    seqkit seq -n ${input} | cut -d' ' -f1 | sort > $2/${name}
}

export -f indexing_step2

# Use GNU Parallel to run the process_value function in parallel
parallel indexing_step2 ::: ${values[@]} ::: ${FOLDER_PATH}

# Indices in R1 or R2
comm -23 ${FOLDER_PATH}/index_${SAMPLE_NAME}_R1_trimmed_phiX174_${ref}  ${FOLDER_PATH}/index_${SAMPLE_NAME}_R2_trimmed_phiX174_${ref}  > ${FOLDER_PATH}/only_R1
comm -23 ${FOLDER_PATH}/index_${SAMPLE_NAME}_R2_trimmed_phiX174_${ref}  ${FOLDER_PATH}/index_${SAMPLE_NAME}_R1_trimmed_phiX174_${ref}  > ${FOLDER_PATH}/only_R2

# Indices in R1 and R2
comm -12 ${FOLDER_PATH}/index_${SAMPLE_NAME}_R1_trimmed_phiX174_${ref}  ${FOLDER_PATH}/index_${SAMPLE_NAME}_R2_trimmed_phiX174_${ref}  > ${FOLDER_PATH}/R1_and_R2

# Find the common indices of unpaired file and only_R1 or only_R2 files
comm -12 ${FOLDER_PATH}/${SAMPLE_NAME}_unpaired.txt ${FOLDER_PATH}/only_R1 > ${FOLDER_PATH}/unpaired_index_1
comm -12 ${FOLDER_PATH}/${SAMPLE_NAME}_unpaired.txt ${FOLDER_PATH}/only_R2 > ${FOLDER_PATH}/unpaired_index_2

# Building fastq files
seqtk subseq ${FOLDER_PATH}/${SAMPLE_NAME}_R1_trimmed_phiX174_${ref}.fq  ${FOLDER_PATH}/R1_and_R2  >  ${FOLDER_PATH}/${SAMPLE_NAME}_R1.fastq
seqtk subseq ${FOLDER_PATH}/${SAMPLE_NAME}_R2_trimmed_phiX174_${ref}.fq  ${FOLDER_PATH}/R1_and_R2  >  ${FOLDER_PATH}/${SAMPLE_NAME}_R2.fastq


# Building fastq files
seqtk subseq ${FOLDER_PATH}/${SAMPLE_NAME}_R1_trimmed_phiX174_${ref}.fq  ${FOLDER_PATH}/R1_and_R2  >  ${FOLDER_PATH}/${SAMPLE_NAME}_R1.fastq
seqtk subseq ${FOLDER_PATH}/${SAMPLE_NAME}_R2_trimmed_phiX174_${ref}.fq  ${FOLDER_PATH}/R1_and_R2  >  ${FOLDER_PATH}/${SAMPLE_NAME}_R2.fastq

seqtk subseq ${FOLDER_PATH}/${SAMPLE_NAME}_R1_trimmed.fq ${FOLDER_PATH}/unpaired_index_1 > ${FOLDER_PATH}/UN_R1.fastq
seqtk subseq ${FOLDER_PATH}/${SAMPLE_NAME}_R2_trimmed.fq ${FOLDER_PATH}/unpaired_index_2 > ${FOLDER_PATH}/UN_R2.fastq

# Merge two reverse and forward unpaired files together
paste -d '\n' --serial ${FOLDER_PATH}/UN_R1.fastq ${FOLDER_PATH}/UN_R2.fastq  > ${FOLDER_PATH}/${SAMPLE_NAME}_UN.fastq

# Doing the statistics
bash statistics.sh  ${FOLDER_PATH}/${SAMPLE_NAME}_UN.fastq  ${FOLDER_PATH}/${SAMPLE_NAME}_UN.fastq.stats
#./statistics.sh ${FOLDER_PATH}/${SAMPLE_NAME}_R1.fastq  ${FOLDER_PATH}/${SAMPLE_NAME}_R2.fastq  ${FOLDER_PATH}/${SAMPLE_NAME}_R1.fastq.stats  ${FOLDER_PATH}/${SAMPLE_NAME}_R2.fastq.stats

end_time=$(date +%s)
execution_time=$((end_time - start_time))
echo "Split and sorting ended! (time of execution: ${execution_time}s)"

echo

############################################################################################# Statistics

echo "Computing statistics started!"
start_time=$(date +%s)

list1=("${FOLDER_PATH}/${SAMPLE_NAME}_R1.fq" "${FOLDER_PATH}/${SAMPLE_NAME}_R1_trimmed.fq")
list2=("${FOLDER_PATH}/${SAMPLE_NAME}_R1_trimmed_phiX174.fq" "${FOLDER_PATH}/${SAMPLE_NAME}_R1_trimmed_phiX174_${ref}.fq" "${FOLDER_PATH}/${SAMPLE_NAME}_R1.fastq")
statistics() {
    local input1=$1
    input2="${input1/_R1/_R2}"
    output1="${input1%.fq}.stats"
    output2="${input2%.fq}.stats"
    bash statistics.sh ${input1} ${input2} ${output1} ${output2}
}

export -f statistics

# Use GNU Parallel to run the process_value function in parallel
parallel statistics ::: ${list1[@]}
parallel statistics ::: ${list2[@]}

end_time=$(date +%s)
execution_time=$((end_time - start_time))
echo "Computing statistics ended! (time of execution: ${execution_time}s)"

echo

############################################################################################## Merging statistics

echo "Merging statistics started!"
start_time=$(date +%s)

# Stats files to be merged
files=(${FOLDER_PATH}/${SAMPLE_NAME}_R1.stats ${FOLDER_PATH}/${SAMPLE_NAME}_R1_trimmed.stats ${FOLDER_PATH}/${SAMPLE_NAME}_R1_trimmed_phiX174.stats ${FOLDER_PATH}/${SAMPLE_NAME}_R1_trimmed_phiX174_${ref}.stats ${FOLDER_PATH}/${SAMPLE_NAME}_R2.stats ${FOLDER_PATH}/${SAMPLE_NAME}_R2_trimmed.stats ${FOLDER_PATH}/${SAMPLE_NAME}_R2_trimmed_phiX174.stats ${FOLDER_PATH}/${SAMPLE_NAME}_R2_trimmed_phiX174_${ref}.stats ${FOLDER_PATH}/${SAMPLE_NAME}_R1.fastq.stats ${FOLDER_PATH}/${SAMPLE_NAME}_R2.fastq.stats ${FOLDER_PATH}/${SAMPLE_NAME}_UN.fastq.stats)

# Output file
output_file="${FOLDER_PATH}/${SAMPLE_NAME}_summary.stats"

# Check if there are any files provided
if [ ${#files[@]} -eq 0 ]; then
    echo "No files provided."
    exit 1
fi

# Extract header from the first file and write it to the output file
head -n 1 "${files[0]}" > "$output_file"

# Loop through all files and append their content (excluding the header) to the output file
for file in "${files[@]}"; do
    tail -n +2 "$file" >> "$output_file"
done

# Extract the last three lines
last_three_lines=$(tail -n 3 "$output_file")

# Extract each line separately
line1=$(echo "$last_three_lines" | sed -n '1p')
line2=$(echo "$last_three_lines" | sed -n '2p')
line3=$(echo "$last_three_lines" | sed -n '3p')

# Extract individual columns from each line
file1=$(echo "$line1" | gawk '{print $1}')
file2=$(echo "$line2" | gawk '{print $1}')
file3=$(echo "$line3" | gawk '{print $1}')

col2_1=$(echo "$line1" | gawk '{print $2}')
col2_2=$(echo "$line2" | gawk '{print $2}')
col2_3=$(echo "$line3" | gawk '{print $2}')

col3_1=$(echo "$line1" | gawk '{print $3}')
col3_2=$(echo "$line2" | gawk '{print $3}')
col3_3=$(echo "$line3" | gawk '{print $3}')

col4_1=$(echo "$line1" | gawk '{print $4}')
col4_2=$(echo "$line2" | gawk '{print $4}')
col4_3=$(echo "$line3" | gawk '{print $4}')

col5_1=$(echo "$line1" | gawk '{print $5}')
col5_2=$(echo "$line2" | gawk '{print $5}')
col5_3=$(echo "$line3" | gawk '{print $5}')

col7_1=$(echo "$line1" | gawk '{print $7}')
col7_2=$(echo "$line2" | gawk '{print $7}')
col7_3=$(echo "$line3" | gawk '{print $7}')

# Calculate the new values
name=$(echo "$file3" | sed 's/_UN\.fastq$//')

col2_sum=$((col2_1 + col2_2 + col2_3))
col3_sum=$((col3_1 + col3_2 + col3_3))

col4_min=$(echo "$col4_1" "$col4_2" "$col4_3" | gawk '{min = $1; for (i=1; i<=NF; i++) if ($i < min) min = $i; print min}')
col5=$(echo "$col5_1")

# Use `bc` to maintain 10 decimal places for the average
col6=$(echo "scale=10; $col2_sum / $col3_sum" | bc)

col7_max=$(echo "$col7_1" "$col7_2" "$col7_3" | gawk '{max = $1; for (i=1; i<=NF; i++) if ($i > max) max = $i; print max}')

# Convert only col6 to 10 decimal places
col6=$(printf "%.10f" "$col6")

# Create the new line
new_line="$name $col2_sum $col3_sum $col4_min $col5 $col6 $col7_max"

# Append the new line to the file
echo "$new_line" >> "$output_file"

# Remove ".fq" from the file
sed -i 's/.fq//g' "$output_file"

end_time=$(date +%s)
execution_time=$((end_time - start_time))
echo "Merging statistics ended! (time of execution: ${execution_time}s)"

echo

############################################################################################### Compressing the outputs

echo "Compressing the outputs started!"
start_time=$(date +%s)

# Compressing the files
files=("${FOLDER_PATH}/${SAMPLE_NAME}_R1.fastq" "${FOLDER_PATH}/${SAMPLE_NAME}_R2.fastq" "${FOLDER_PATH}/${SAMPLE_NAME}_UN.fastq" "${FOLDER_PATH}/${SAMPLE_NAME}_unpaired.txt")
zipping() {
    local input=$1
    pbzip2 -f -p2 ${input}
}

export -f zipping

# Use GNU Parallel to run the zipping function in parallel
parallel zipping ::: ${files[@]}

end_time=$(date +%s)
execution_time=$((end_time - start_time))
echo "Compressing the outputs ended! (time of execution: ${execution_time}s)"

echo

############################################################################################# Removing unneccessary files

echo "Removing intermediate files!"

rm ${FOLDER_PATH}/index_${SAMPLE_NAME}_R1_trimmed_phiX174_${ref} ${FOLDER_PATH}/index_${SAMPLE_NAME}_R2_trimmed_phiX174_${ref}  ${FOLDER_PATH}/only_R1 ${FOLDER_PATH}/only_R2 ${FOLDER_PATH}/R1_and_R2 ${FOLDER_PATH}/unpaired_index_1 ${FOLDER_PATH}/unpaired_index_2 ${FOLDER_PATH}/UN_R1.fastq ${FOLDER_PATH}/UN_R2.fastq
rm ${FOLDER_PATH}/${SAMPLE_NAME}*fq ${FOLDER_PATH}/${SAMPLE_NAME}_R*.stats ${FOLDER_PATH}/${SAMPLE_NAME}_UN*.stats
