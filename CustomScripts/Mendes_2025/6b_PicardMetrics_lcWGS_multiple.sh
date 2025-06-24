#!/bin/bash

#Based on #get_multistats.sh" by Vitor C. Sousa, adapted for the Squalius GBS data by Sofia L. Mendes on 22.09.2020

#Edited for lcWGS alignments against S. cephalus after processing the bam files with Picard on 07.06.2022
#Edited again on 25.03.2023 for the new 2023 batch of samples
#Edited on 30.05.2023 to run on the final 125 individuals of the low coverage dataset

# reference, that is, the SQUALIUS CEPHALUS reference genome
REFERENCE="/media/shared/smendes/SCEPHALUS_ref_genome/ncbi-genomes-2022-05-31/GCA_022829025.1_ASM2282902v1_genomic.fna";
echo "Squalius cephalus reference genome is file=${REFERENCE}";

#folder where the script should be run and where we want to keep the output files
PICARDFOLDER="/path/to/folder/where/you/want/to/save/output";
echo ${PICARDFOLDER};

#folder where the bam files are
BAMFOLDER="/path/to/final/processed/bam/files";
echo ${BAMFOLDER};

#popmap with the sample/individual id (sampleID) and the population the individual belongs to
#The file should have two columns, one with the ID and another with the population
#this is the same popmap from the Picard mark duplicates step
POPMAP=/path/to/popmap/popmap_125inds_lcWGS_scephalus_merged.txt;


# read the file with the sample id and perform the analyses for each sample
while read -r sampleID popID; do		
	# ID of library in $sampleID
	echo "**********************\nProcessing sample sorted_markeddup_${sampleID}\n***************************\n"
		
	# get statistics about mapping quality 
	picard CollectMultipleMetrics \
		I=${BAMFOLDER}/sorted_markeddup_${sampleID}.bam \
		O=${PICARDFOLDER}/sorted_markeddup_${sampleID}.metrics \
		R=${REFERENCE} \
		PROGRAM=CollectSequencingArtifactMetrics \
		PROGRAM=CollectAlignmentSummaryMetrics \
		PROGRAM=CollectInsertSizeMetrics \
		PROGRAM=CollectGcBiasMetrics \
		PROGRAM=CollectQualityYieldMetrics >> ${PICARDFOLDER}/mapstatsorted_markeddup_${sampleID}.log 2>&1
		
done < ${POPMAP}
