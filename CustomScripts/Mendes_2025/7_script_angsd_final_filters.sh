#!/bin/bash

#Sofia L. Mendes
#15.06.2021

# Script to obtain genotype likelihoods from the bam files with ANGSD


#This is to be run in ANGSDFOLDER (see bellow)

# SETTINGS (change here the file names and folder names)
#folder where we want to save the outputs
ANGSDFOLDER=/path/to/angsd/output/folder;
# folder where we have the final indexed bam files
BAMFOLDER=/path/to/folder/with/indexed/bam/files;
# path to the reference genome
REFERENCE=/media/shared/smendes/SCEPHALUS_ref_genome/ncbi-genomes-2022-05-31/GCA_022829025.1_ASM2282902v1_genomic.fna;
# path to the popmap file with the individual ID and populations
#this is the same popmap from the Picard mark duplicates step
POPMAP=/path/to/popamp/popmap_125inds_lcWGS_scephalus_merged.txt;
# file with the names of the chromossomes/ scaffolds of the reference genome we want to analyse
# For S. cephalus, we will look at sites on the 25 long scaffolds that correspond to the 25 chromossomes
CHROM=/path/to/text/file/with/chromossome/names/chromossome_names_scephalus.txt;

# read the POPMAP file and create a file with the list of all the bam files
# for the individuals in out POPMAP
while read -r indID popID; 
do
	# >> means that we keep appending this to the file 
	echo ${BAMFOLDER}/sorted_markeddup_${indID}.bam >> ${ANGSDFOLDER}/listBamFiles.txt; 
done < ${POPMAP}


# obtain genotype likelihoods with ANGSD
# -b/-bam is the list of bam files
# -ref indicates the reference genome
# -out is the location (${ANGSDFOLDER}) and prefix (lcwgs_squalius_74inds) of all output files
# -nThreads is the number of threads
# -doVcf - if set to one, output in vcf format as well 
# -uniqueOnly - if set to 1 it means "retain only uniquely mapping reads"
# -remove_bads - 
# -only_proper_pairs - If set to 1, consider only proper pairs. If set to 0, use reads even if they are not in proper pairs
# -trim - option for trimming at the tips of reads. If set to 0, no trimming.
# -minMapQ - minimum mapping quality
# -minQ - minimum base quality
# -doCounts - 
# -minInd - Only keep sites with at least minIndDepth (default is 1) from at least [int] individuals
# -setMinDepthInd - Define minIndDepth for the previous option, i.e., discard individual if sequencing depth for an individual is below [int]. This filter is only applied to analysis which are based on counts of alleles i.e. analysis that uses -doCounts
# -setMinDepth - Discard site if total sequencing depth (all individuals added together) is below [int].
# -setMaxDepthInd - Discard individual if sequencing depth for an individual is above [int], i.e., discard individual if sequencing depth for an individual is above [int]. This filter is only applied to analysis which are based on counts of alleles i.e. analysis that uses -doCounts
# -GL - genotype likelihood model (1 means Samtools, 2 means GATK, check the manual for more options http://www.popgen.dk/angsd/index.php/Genotype_Likelihoods)
# -doGlf - 
# -doMajorMinor (if 1, Infer major and minor from GL (both alleles are inferred from genotype likelihoods) ; If 2, Infer major and minor from allele counts; for more options see manual http://www.popgen.dk/angsd/index.php/Major_Minor)
# -skipTriallelic
# -dosnpstat - if set to 1 it will calculate statistics for the SNPs
# -doMaf - 
# -doPost - 
# -SNP_pval - 
# -doHWE - if set to 1, estimate the divination from HWE for each site 


#Run ANGSD per chromossome
while read -r chromID; 
do
angsd -b ${ANGSDFOLDER}/listBamFiles.txt -ref ${REFERENCE} -r ${chromID} -out ${ANGSDFOLDER}/lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_${chromID}_sc -nThreads 4 \
   -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -minMapQ 30 -minQ 20 -doCounts 1 -minInd 124 -setMinDepthInd 2 -setMaxDepthInd 14 -setMaxDepth 1082 \
   -GL 1 -doGlf 2 -doMajorMinor 1 -skipTriallelic 1 -dosnpstat 1 -doMaf 1 -doPost 1 -SNP_pval 1e-6 -doHWE 1 -sb_pval 0.05 -qscore_pval 0.05 -edge_pval 0.05 -mapq_pval 0.05 
done < ${CHROM}

