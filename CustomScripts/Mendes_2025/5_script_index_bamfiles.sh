#!/bin/bash

#Sofia L. Mendes edited on 31/05/2022 and again on 24.03.2023

# Script to index all bam files for ANGSD

#This is to be run in BAMFOLDER (see bellow)

# SETTINGS (change here the file names and folder names)
#folder where we have the new bam files that need indexing
BAMFOLDER=/path/to/folder/with/bam/files/to/index;
#popmap with the individual ID and the population the individual belongs to.
#The file should have two columns, one with the ID and another with the population
#this is the same popmap from the Picard mark duplicates step
POPMAP=/path/to/popmap/popmap_125inds_lcWGS_scephalus_merged.txt;

# read the POPMAP file and for each individual (ind), 
# index the bam file with samtools
while read -r ind pop; 
do
	# index the bam file with samtools
	samtools index -@ 4 ${BAMFOLDER}/sorted_markeddup_${ind}.bam
done < ${POPMAP}

