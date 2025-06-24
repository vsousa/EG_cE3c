#!/bin/bash

#Sofia L. Mendes 

#edited on 23/03/2023 to merge the bam files of the individuals that were sequenced twice by Novogene and thus have, after mapping
#against the S. cephalus reference genome, two separate BAM files. These will be merged before marking duplicates with Picard because there may be 
#duplicates between the two files that we would otherwise miss. 

#edited on 25.05.2023 for the 125 individuals of the lcWGS dataset

# Script to merge BAM files

#This is to be run in BAMFOLDER (see bellow)

# SETTINGS (change here the file names and folder names)
#folder where we want to save the new bam files and where we also have the list of files to merge
BAMFOLDER=/path/to/folder/where/you/want/to/save/your/merged/bam/files;
# folder where we have the unmerged bam files after mapping with BWA to S. cephalus reference genome
UNMERGED=/path/to/folder/where/you/have/your/bam/files;
#The two folders above can be the same folder, depends on how one is organizing the pipeline in different folders
# path to the list of files to merge and names of the output merged file
#this should be a text file with three tab delimited columns
#columns 1 and 2 are the names of the two bam files to be merged
#column 3 is the output bam file name
BAMLIST=/path/to/file/with/list/of/bam/to/merge/bam_to_merge_2files_lcWGS_125inds_scephalus.txt;

# read the BAMLIST file and for each individual, 
# merge the files on column 1 (fileA) and column 2 (fileB) using samtools
# then sort the output with samtools
# and save it with the name fileOUT
while read -r fileA fileB fileOUT; 
do
	# merge and save the output to a temporary file
	samtools merge -@ 5 ${BAMFOLDER}/${fileOUT}.bam ${UNMERGED}/sorted_${fileA}.bam ${UNMERGED}/sorted_${fileB}.bam
	# sort the temporary file and save the output in a final bam file
	samtools sort -@ 5 -o ${BAMFOLDER}/sorted_${fileOUT}.bam ${BAMFOLDER}/${fileOUT}.bam # sort bam file, save in this (-o) file
	# remove the intermediate temporary bam file
	rm ${BAMFOLDER}/${fileOUT}.bam
	# remove the two initial bam files
	#rm ${UNMERGED}/sorted_${fileA}.bam
	#rm ${UNMERGED}/sorted_${fileB}.bam
done < ${BAMLIST}




