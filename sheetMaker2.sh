#!/bin/bash

# Open output directory.
cd $3/sheetMaker

# Create csv file and add header.
touch "$2".csv
echo "sample,fastq_1,fastq_2,strandedness" > "$2".csv

# Return to working directory.
cd $1

# Dictionaries for checking if both R1 and R2 filepaths are found before echoing to the csv.
declare -A fastq_1_files
declare -A fastq_2_files


# Rename all files to end with fq.fz
find "$1" -type f -regex ".*\.unclassified\.gz" | sort | while read -r file; do
    # Rename the file.
    newname=$(echo $file | sed 's/\(.*\)\.fq\.unclassified\.gz/\1.unclassified.fq.gz/')
    mv $file $newname
done

# Find all files in the given filepath that are of filetype .unclassified.fq.gz and loop through each in a sorted order.
find "$1" -type f -regex ".*\.unclassified\.fq\.gz" | sort | while read -r file; do

    # Extract the patient sample identifier and whether reverse or not (5' to 3' or 3' to 5') by cutting.
    sample=$(basename $file | cut -d'_' -f1-4)
    direction=$(basename $file | cut -d'_' -f5)

    # Add the filepath to the appropriate dictionary of filepaths to temporarily store it until I encounter the next filepath in the pair.
    if [[ $direction == "1" ]]; then
        fastq_1_files["$sample"]="$file"
    else
        fastq_2_files["$sample"]="$file"
    fi

    # Wait until there exists an entry for both filepaths in the dictionaries before echoing to the csv. Then clear the dictionaries.
    if [ -n "${fastq_1_files[$sample]}" ] && [ -n "${fastq_2_files[$sample]}" ]; then
        echo "${sample},${fastq_1_files[$sample]},${fastq_2_files[$sample]},unstranded" >> "$3"/sheetMaker/"$2".csv
        fastq_1_files["$sample"]=""
        fastq_2_files["$sample"]=""
    fi

done

echo "$3"/sheetMaker/"$2".csv

