#!/bin/bash

#Sofia L. Mendes
#31/05/2022
#map the squalius lcWGS against the squalius cephalus reference genome
#edited on 16/03/2023 to map new low coverage individuals against the squalius cephalus reference genome


# Script to run BWA and samtools to map against the genome and create BAM files

#This is to be run in BWAFOLDER (see bellow)

# SETTINGS (change here the file names and folder names)
#folder where we want to save the bam files and where we also have the popmap
BWAFOLDER=/path/to/where/you/want/to/save/bamfiles;
# path to the reference genome
REFERENCE=/media/shared/smendes/SCEPHALUS_ref_genome/ncbi-genomes-2022-05-31/GCA_022829025.1_ASM2282902v1_genomic.fna;
# path to the fastq files after fastp
FQ_FILES=/path/to/folder/where/you/have/your/fastp/outputs;
# path to the popmap file with the individual ID and populations
POPMAP=/path/to/popmap/file/popmap_125inds_all_files.txt;
#popmap should be a two column tab separated file where in the first column you have the file names (indID) and in the second column the sampling location they belong to
#see example below (a portion of the actual popmap used)
#AT18549_EKDN230012706-1A_HV33YDSX5_L4	oitaven
#AT18551_EKDN220047963-1A_HLWCNDSX5_L2	oitaven
#AT18551_EKDN220047963-1A_HTWYYDSX5_L1	oitaven
#AT18551_EKDN220047963-1A_HYCHVDSX5_L1	oitaven
#AT17974_EDSW200006862-1a_HCLVMDSXY_L2	minho_arnoia
#AT17975_FDSW210081761-1r_HVV2WDSXY_L3	minho_arnoia
#AT17975_FDSW210081761-1r_HWLYCDSXY_L3	minho_arnoia
#SCPeg24_FDSW210081762-1r_HVV2WDSXY_L3	pego
#SCPeg29_FDSW210081763-1r_HVV2WDSXY_L3	pego
#AT18766_EDSW200006864-1a_HCLVMDSXY_L2	douro_tera
#AT18768_EDSW200006863-1a_HCLVMDSXY_L2	douro_tera
#AT16757_EKDN230012707-1A_HV33YDSX5_L2	douro_curueno
#AT16758_EKDN220047965-1A_HLWCNDSX5_L2	douro_curueno
#AT16758_EKDN220047965-1A_HTWYYDSX5_L1	douro_curueno
#AT16758_EKDN220047965-1A_HYCHVDSX5_L1	douro_curueno


# read the POPMAP file and for each individual ID, 
# search the corresponding FASTQ file, uncompress it,
# and align it using BWA and create the BAM files
while read -r indID popID; 
do
	# use bwa to align and samtools to output sorted bamfiles
	bwa mem -t 6 ${REFERENCE} ${FQ_FILES}/${indID}_1.fq.gz ${FQ_FILES}/${indID}_2.fq.gz | \
	samtools view -b -@ 6 | samtools sort -@ 6 -o ${BWAFOLDER}/sorted_${indID}.bam # sort bam file by coordinate and save in this (-o) file
done < ${POPMAP}

