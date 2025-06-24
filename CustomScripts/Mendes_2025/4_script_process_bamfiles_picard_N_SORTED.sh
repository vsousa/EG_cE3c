#!/bin/bash

#Sofia L. Mendes

#Process the lcWGS alignments against S. cephalus

#This is to be run in PICARDFOLDER (see bellow)

# SETTINGS (change here the file names and folder names)
#BWAFOLDER is where the bam files of the alignments against the s. cephalus reference genome are saved
BWAFOLDER=/path/to/folder/containing/your/sorted/bam/files;
# path to the reference genome
REFERENCE=/media/shared/smendes/SCEPHALUS_ref_genome/ncbi-genomes-2022-05-31/GCA_022829025.1_ASM2282902v1_genomic.fna;
#PICARDFOLDER is the folder where we want to save the new bam files
PICARDFOLDER=/path/to/output/folder;
# path to the popmap file with the individual ID and populations
POPMAP=/path/to/popmap/file/popmap_125inds_lcWGS_scephalus_merged.txt;
#popmap should be a two column tab separated file where in the first column you have the file names (indID) and in the second column the sampling location they belong to
#note some bam files names have changes compared to the popmap used for mapping with bwa,
#as bam files from individuals sequenced in multiple runs have been merged
#see example below (from the popmap actually used)
#AT18549_EKDN230012706-1A_HV33YDSX5_L4	oitaven
#AT18551_merged	oitaven
#AT17974_EDSW200006862-1a_HCLVMDSXY_L2	minho_arnoia
#AT17975_merged	minho_arnoia
#SCPeg24_FDSW210081762-1r_HVV2WDSXY_L3	pego
#SCPeg29_FDSW210081763-1r_HVV2WDSXY_L3	pego
#AT18766_EDSW200006864-1a_HCLVMDSXY_L2	douro_tera
#AT18768_EDSW200006863-1a_HCLVMDSXY_L2	douro_tera
#AT16757_EKDN230012707-1A_HV33YDSX5_L2	douro_curueno
#AT16758_merged	douro_curueno


while read -r indID popID; 
do
	# ID of library in $indID
	echo "**********************\nProcessing sample sorted_${indID}\n***************************\n"
		
	# Verify mate-pair information between mates and fix if needed with Picard 
	# Also sort the output by queryname so that Markduplicates can mark and remove all duplicates
	picard FixMateInformation \
		I=${BWAFOLDER}/sorted_${indID}.bam \
		O=${PICARDFOLDER}/n_sorted_fixedmate_${indID}.bam \
		TMP_DIR=/media/shared/smendes/6_process_bam_lcWGS_125inds_scephalus/temp_dir \
		SORT_ORDER=queryname >> ${PICARDFOLDER}/n_fixedmate_${indID}.log 2>&1
	# Mark duplicates
	picard MarkDuplicates \
		I=${PICARDFOLDER}/n_sorted_fixedmate_${indID}.bam \
		O=${PICARDFOLDER}/n_tempmarkeddup_${indID}.bam \
		M=${PICARDFOLDER}/marked_dup_metrics_${indID}.txt \
		REMOVE_DUPLICATES=true \
		TMP_DIR=/media/shared/smendes/6_process_bam_lcWGS_125inds_scephalus/temp_dir \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 >> ${PICARDFOLDER}/n_tempmarkeddup_${indID}.log 2>&1
	#remove the name sorted output of FixMateInformation
	rm ${PICARDFOLDER}/n_sorted_fixedmate_${indID}.bam
	#sort by coordinates the temporary marked duplicates bam file with samtools v1.19
	samtools sort -@ 6 ${PICARDFOLDER}/n_tempmarkeddup_${indID}.bam -o ${PICARDFOLDER}/sorted_markeddup_${indID}.bam
	#remove the temporary marked duplicates unsorted file
	rm ${PICARDFOLDER}/n_tempmarkeddup_${indID}.bam
done < ${POPMAP}

