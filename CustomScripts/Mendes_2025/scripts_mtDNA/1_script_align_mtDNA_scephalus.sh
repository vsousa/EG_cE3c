#!/bin/bash

# Author: Vitor C. Sousa
# Date: 11.09.2018

# Author: Vitor C. Sousa
# Date: 11.09.2018

#Sofia L. Mendes edited on 27/06/2023 to map the mtDNA of the 125 individuals against the squalius cephalus mitocondrial genome
#unlike in the nuclear genome, unmapped reads will be removed 

# Script to run BWA and Picard to create BAM files with only sorted mapped reads after removal of unmapped reads and duplicates

#This is to be run in BWAFOLDER (see bellow)

# SETTINGS (change here the file names and folder names)
#folder where we want to save the bam files and where we also have the popmap
BWAFOLDER=path/to/folder;
# path to the reference genome
REFERENCE=/path/to/reference/genome/file/NC_031540.1.fasta;
# path to the fastq files
FQ_FILES=/pat/to/fastq/files/after/fastp;
# path to the popmap file with the individual ID and populations
POPMAP=/path/to/popmap.txt;


# read the POPMAP file and for each individual ID, 
# search the corresponding FASTQ file, uncompress it,
# and align it using BWA and create the BAM files
while read -r indID popID rgID; 
do
	# use bwa to align and samtools to output the bamfile
	# the RG is used to add sample information to each bam file
	# this sample information is required for freeBayes to distinguish 
	# bam files from different individuals (i.e. different samples ID and SM)
	# align and create bam file
	GROUPINFO='@RG\tID:'${rgID}'\tSM:'${rgID}'\tLB:lib1';
	bwa mem -t 8 -R ${GROUPINFO} ${REFERENCE} ${FQ_FILES}/${indID}_1.fq.gz ${FQ_FILES}/${indID}_2.fq.gz | \
	samtools view -bS -@ 8 | samtools sort -@ 8 -o ${BWAFOLDER}/sorted_${indID}.bam # sort bam file, save in this (-o) file
	# remove unmapped reads and sort again because the vast majority of reads do not map to the Squalius cephalus mitochondrial genome
	# Note that for the nuclear data we do not do this, that is, we do not remove the unmapped reads
	# The -F option means "Do not output alignments with any bits set in INT present in the FLAG field
	# and the 4 FLAG means "unmapped". This means we will not keep unmapped reads.
	samtools view -b -F 4 -@ 8 ${BWAFOLDER}/sorted_${indID}.bam | samtools sort -@ 8 -o ${BWAFOLDER}/mapped_sorted_${indID}.bam
	# remove the bam file that still had all the unmapped reads
	rm ${BWAFOLDER}/sorted_${indID}.bam
	# Verify mate-pair information between mates and fix if needed with Picard 
	picard FixMateInformation \
		I=${BWAFOLDER}/mapped_sorted_${indID}.bam \
		O=${BWAFOLDER}/n_sorted_fixedmate_${indID}.bam \
		TMP_DIR=/media/shared/smendes/11_MTDNA_tree/1_bam_files_mtcephalus_all_inds/temp_dir \
		SORT_ORDER=queryname >> ${BWAFOLDER}/n_fixedmate_${indID}.log 2>&1
	#remove the mappped_sorted bam file (prior to FixMateInformation)
	rm ${BWAFOLDER}/mapped_sorted_${indID}.bam
	# Mark duplicates
	picard MarkDuplicates \
		I=${BWAFOLDER}/n_sorted_fixedmate_${indID}.bam \
		O=${BWAFOLDER}/n_tempmarkeddup_${indID}.bam \
		M=${BWAFOLDER}/marked_dup_metrics_${indID}.txt \
		REMOVE_DUPLICATES=true \
		TMP_DIR=/media/shared/smendes/11_MTDNA_tree/1_bam_files_mtcephalus_all_inds/temp_dir \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 >> ${BWAFOLDER}/n_tempmarkeddup_${indID}.log 2>&1
	#remove the name sorted output of FixMateInformation
	rm ${BWAFOLDER}/n_sorted_fixedmate_${indID}.bam
	#sort by coordinates the temporary marked duplicates bam file with samtools v1.19
	samtools sort -@ 6 ${BWAFOLDER}/n_tempmarkeddup_${indID}.bam -o ${BWAFOLDER}/mapped_sorted_markeddup_${indID}.bam
	#remove the temporary marked duplicates unsorted file
	rm ${BWAFOLDER}/n_tempmarkeddup_${indID}.bam
	# check file
	# samtools view -H ${BWAFOLDER}/mapped_sorted_markeddup_${indID}.bam  | grep '@RG'
done < ${POPMAP}

