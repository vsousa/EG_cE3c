#!/bin/bash
# Vitor Sousa, last updated 13/01/2020 
# Please e-mail me if you find any problems
# vmsousa[at]fc.ul.pt

# REQUIRED PROGRAMS
# - htslib, samtools, vcftools, bcftools, R, tabix

# REQUIRED FILES
# - vcffile
# - indpopinfo text file with two columns with individual ID and corresponding POP
#        NOTE: indpopinfo file defines the order of individuals and populations.
#              Thus, this file can have just a sub-set of individuals from the VCF,
#              for which you want to obtain the SFS.

################################################
# 1 - SETTINGS 
################################################

# As input you should provide a filtered VCF file
# Filtered based on DP (try to have DP>10)
# Filter based on missing data per site, per individual.
# NOTE: Do not use a VCF with a MAF filter applied, 
#       as the SFS is affected by the MAF filters!

# tag for VCF file (vcf file with format "vcffile".vcf)
vcffile="final_NO_maf_merged_soaradeMS_DPfiltered_inds_less_50missing";
# tag for the resulting files
vcftag=final_${vcffile};
# tag for file with ind pop info (i.e. information about the individuals and corresponding pop)
# this can have just a sub-set of individuals and populations of the VCF.
indidpop="popmap_no_out_soaradeMS_50missing_carol_ocrezacanha_pyrsul.txt";
# tag for folder with output data
outSFSfolder="block_SFS";
# block size in bp
block_length=200;
# minimum sample size per pop (must be given as a string separated by comma)
# order of pops must be the same as in indidpop file!
ind_threshold="2,3,2";
# minimum median distance between consecutive SNPs in a good block
# blocks with a median distance equal or lower than this are discarded.
dist_threshold=2;
# path to folder with Rscripts
scriptsRfolder="./Scripts_VCFtoSFS";
# random sampling of individuals? FALSE for deterministic sampling (always selecting the top individuals ranked according to content of data)
randomInd=F;
# seed for random sampling
seed=6126151;

############################################################################
# 2. Get the genotype matrix and other relevant information from the VCF
############################################################################

# Check that your vcf is indexed with tabix and compressed with bgzip
#bgzip -c ${vcffile}.vcf > ${vcffile}.vcf.gz;
tabix -p vcf ${vcffile}.vcf.gz;

# 2.1. Get the genotypes for each individual at each site (GT field)
bcftools="/home/ALUNOSFC/fc43159/programas/bcftools/bcftools";
${bcftools} query -f '[%GT\t]\n' ${vcffile}.vcf.gz > ${vcftag}.GT;
# Use sed to re-code genotypes:
# Replace 0/0 or 0|0 by 0 (homozygote reference)
# Replace 0/1 or 0|1 by 1 (heterozygote)
# Replace 1/0 or 1|0 by 1 (heterozygote)
# Replace 1/1 or 1|1 by 2 (homozygote alternative)
# Replace ./. or .|. by -1 (missing data)
# Replace . by -1 (missing data)
sed -i "s/0\/0/0/g;s/0\/1/1/g;s/1\/0/1/g;s/1\/1/2/g;s/0|0/0/g;s/0|1/1/g;s/1|0/1/g;s/1|1/2/g;s/\.\/\./-1/g;s/\.|\./-1/g;s/\./-1/g;" ${vcftag}.GT;

# 2.2. Get depth of coverage per site per individual (DP field)
${bcftools} query -f '[%DP\t]\n' ${vcffile}.vcf.gz > ${vcftag}.DP;
sed -i "s/\.\t/0\t/g" ${vcftag}.DP;

# 2.3. Get the chromosome (or scaffold) code and position (CHROM POS fields)
# Get the list of positions in the VCF file
${bcftools} query -f '%CHROM\t%POS\n' ${vcffile}.vcf.gz > ${vcftag}.chrpos;

# 2.4. Use vcftools to get several statitiscs
# --site-mean-depth: mean depth per site
# --het: heterozygosity per individual 
# --missing-indv: missing data per individual
# --missing-site: missing data per site
# Get the mean and variance depth of coverage per site
vcftools --gzvcf ${vcffile}.vcf.gz --site-mean-depth --out ${vcftag};
vcftools --gzvcf ${vcffile}.vcf.gz --het --out ${vcftag};
vcftools --gzvcf ${vcffile}.vcf.gz --missing-indv --out ${vcftag};
vcftools --gzvcf ${vcffile}.vcf.gz --missing-site --out ${vcftag};

# 2.5. Get the list of individuals in the VCF file
# create a file with the individuals in the VCF file
zgrep -n -m 1 ^#CHROM ${vcffile}.vcf.gz > indsIn${vcftag};
# remove the first 9 elements of the resulting file, since those are not individual IDs.
cut -f10- indsIn${vcftag} > tmp_indsIn${vcftag};
mv tmp_indsIn${vcftag} indsIn${vcftag};


############################################################################
# 3. Get SFS by resampling individuals with blocks without missing data
############################################################################

# The command line has the following arguments:
# - filename_vcf: filename tag for the file with genotypes (filename_vcf.GT), chromosome position info (filename_vcf.CHRMPOS), etc.
# - indpopinfofilename: file with the ID of individuals and corresponding population
#                        we want to look at. This file must have two columns, 1st indID, 2nd Pop ID.
#                        NOTE: This defines the order of individuals and populations!
#                        This file can have just a sub-set of individuals from the VCF.
# - foldertag: tag for the folder with output SFS
# - block_length: define the length (distance among positions for defining a block)
# - ind_threshold: string with number of diploid individuals per pop (P1, P2, P3), e.g. "5,2,3"
# - dist_threshold: minimum mean distance between consecutive SNPs in a good block

Rscript ${scriptsRfolder}/BlockSFS_jan2020.r ${vcftag} ${indidpop} ${outSFSfolder} ${block_length} ${ind_threshold} ${dist_threshold} ${randomInd} ${seed};


##################################################################
# 4. Generate block-bootstrap replicates
##################################################################

# The command line has the following arguments:
# -foldertag: tag for the folder with SFS
# -block_length: length of blocks used
# -ind_threshold: string with number of diploid individuals per pop (P1, P2, P3), e.g. "5,2,3"
# -dist_threshold: minimum mean distance between consecutive SNPs in a good block
# -numboot number: of bootstrap replicates
# -mafOrDaf string: with "m" for MAF (minor allele frequency) and "d" for DAF (derived allele frequency)
# -randomSNP: T for sampling 1 random SNP per block. F to sample all linked SNPs per block.
# -seed seed of random number generator
#Example of use:
# Rscript blockBootstrap_SFS.r block_SFS 10000 "10,10,10" 2 10 m F 626325
nboot=1000;
mafdaf="m";
randomSNP="F";
Rscript ${scriptsRfolder}/blockBootstrap_SFS.r block_SFS ${block_length} ${ind_threshold} ${dist_threshold} ${nboot} ${mafdaf} ${randomSNP} ${seed};


##################################################################
# 5. Generate block-bootstrap replicates with 1SNP per block
##################################################################

randomSNP="T";
Rscript ${scriptsRfolder}/blockBootstrap_SFS.r block_SFS ${block_length} ${ind_threshold} ${dist_threshold} ${nboot} ${mafdaf} ${randomSNP} ${seed};

