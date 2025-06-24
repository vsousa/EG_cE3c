#!/bin/bash

# S. L. Mendes 06/05/2021
# Script to run fastp on all fq.gz files from the squalius low coverage project 

#S. L. Mendes edited on 08.07.2022 and again on 08.05.2023

#This is to be run in FASTPFOLDER (see bellow)

# SETTINGS (change here the file names and folder names)
#folder where we want to save the new clean fq.gz files and where we also have the individuals/fasta files
FASTPFOLDER=path/to/folder/where/you/want/to/save/output/files;
# path to the folder where we stored the fastq files 
FQ_FILES=path/to/folder/where/you/have/your/fq.gz/files;
# path to the list of file names (each name corresponds to two files - _1 and _2, each containing one of the paired ends, see script below)
#This list is a simple txt file where each line contains 1 file name (e.g. SASx1_EKDN220048020-1A_HKN5GDSX5_L2) minus the _1 or _2
#Note the script below - there is no need to write a file in the list twice, as the script is writen to get both file_name_1 and file_name_2 for each file_name 
FILELIST=/path/to/filelist/file_names_125inds.txt;


# read the POPMAP file and for each individual ID, 
# search the corresponding FASTQ files (_1 and _2)
# and run fastp
while read -r indID; 
do
	# run fastp to get rid of the PolyG tails and keep only reads with 150bp
	fastp -i ${FQ_FILES}/${indID}_1.fq.gz -I ${FQ_FILES}/${indID}_2.fq.gz --dont_overwrite -o ${FASTPFOLDER}/${indID}_1.fq.gz -O ${FASTPFOLDER}/${indID}_2.fq.gz -g --thread 6 -q 20 -j ${FASTPFOLDER}/${indID}.json -h ${FASTPFOLDER}/${indID}.html
done < ${FILELIST}

