#!/usr/bin/env bash

rnaseqOutputDir=$1
nextflowCommand1=$2
nextflowCommand2=$3

# Create work directory for temp files.
cd $rnaseqOutputDir
mkdir work

# Run nextflow (only FastQC + Trimming + FastQC) with the csv filepath from sheetMaker.sh and outdir argument.
$nextflowCommand1

# Create directory for trimmed files to run ribodetector on
cd trimgalore
mkdir fastq_files_for_ribodetector

echo "<------------------------------------------------------------------------------------------------ Step1"

# Find all fq.gz files made by trimgalore and sort them away from the txt files in a seperate folder.
find . -type f -regex ".*\.fq\.gz" -exec mv {} fastq_files_for_ribodetector/ \;

# Enter the folder with only the fq.gz files.
cd fastq_files_for_ribodetector
echo "<------------------------------------------------------------------------------------------------ Step2"

declare -A seenSamples

find . -type f -regex ".*\.fq\.gz" | sort | while read -r file; do
	#echo $file
	sample=$(basename $file | cut -d'_' -f1-4)
	#echo $sample
	if [ -z ${seenSamples[$sample]} ]; then
		file1=$rnaseqOutputDir/trimgalore/fastq_files_for_ribodetector/$sample"_1_val_1.fq.gz"
		file2=$rnaseqOutputDir/trimgalore/fastq_files_for_ribodetector/$sample"_2_val_2.fq.gz" 
		ribodetector -t 45 \
		 -l 151 \
		  -i $file1 $file2 \
		  -e both \
		  -m 22 \
		  --chunk_size 1028 \
		  -o ${file1}".nonrrna.1.fq" ${file2}".nonrrna.2.fq"
	fi
	seenSamples[$sample]=1
done

echo "<------------------------------------------------------------------------------------------------ Step3"

# Do MultiQC on output from ribodetector.
multiqc "$rnaseqOutputDir/trimgalore/fastq_files_for_ribodetector --outdir $rnaseqOutputDir/multiqc"

# Run second sheetMaker OBS! will break if the app doesnt exist in this specific directory, hardcoded for now!
cd "/media/baseripper/KINGSTON/KuzeyK/DGE_Analysis_App_ACT/DGE_Analysis_App"
./sheetMaker2.sh "$rnaseqOutputDir/trimgalore/fastq_files_for_ribodetector DGEA_projSheet2 $rnaseqOutputDir"
cd $rnaseqOutputDir

# Run the rest of the pipeline.
$nextflowCommand2
