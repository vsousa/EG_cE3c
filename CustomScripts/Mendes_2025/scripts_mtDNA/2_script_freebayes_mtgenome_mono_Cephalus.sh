#!/bin/bash


# Sofia L. Mendes
# 14.06.2021

#Edited by Catarina Bernardo -18/08/23

#Script to perform genotype calling on the bam files of the alignments against the Squalius cephalus mitochondrial genome

# SETTINGS (change here the file names and folder names)
FREEBAYESFOLDER=/folder/where/you/want/to/save/your/vcf;
# path to the Squalius cephalus mitochondrial genome
REFERENCE=/path/to/reference/genome/NC_031540.1.fasta;
# path to aligned bam files
BAMFOLDER=/path/to/folder/where/you/have/your/bam/files;
# path to the popmap file with the individual ID and populations
POPMAP=/path/to/popmap_freebayes.txt;


# read the POPMAP file and for each individual ID, 
# create a file with the bam file names indID.bam
while read -r indID popID; 
do
	# >> means that we keep appending this to the file 
	echo ${BAMFOLDER}/${indID}.bam >> ${FREEBAYESFOLDER}/listBamFiles; 
done < ${POPMAP}

#FreeBayes crashed due to limited memory. Use command line below to correct that.
# ulimit -S -s 131072
# call freebayes
freebayes -f ${REFERENCE} -L ${FREEBAYESFOLDER}/listBamFiles -p 1 --min-mapping-quality 30 --min-base-quality 20 --report-monomorphic --vcf ${FREEBAYESFOLDER}/mtdna_data_withmono_freebayes.vcf
#bgzip ${FREEBAYESFOLDER}/mtdna_data_withmono_freebayes.vcf
