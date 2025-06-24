#!/bin/bash

#Sofia L. Mendes
#02.11.2023

# Script to obtain NJ tree using PCAngsd


#Folder where we want to run PCAngsd to obtain Neightbour-joining trees
#Script should be run in this folder
NJFOLDER=/path/to/folder
#Name of the beagle.gz input file and the output files
BFILE="merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc"
#MAF filter to apply
MAF=0.01
#MAF term to add to output file name
MAFTERM=MAF001
#File with sample names to add on the tree (same order as in the beagle file, of course)
TREESAMPLES=/media/shared/smendes/14_NJtree_PCAngsd/treesamples_125inds_lcWGS_scephalus_merged.txt


#Path to PCAngsd
PCANGSD=/home/smendes/programs/pcangsd/pcangsd.py


#Perform a PCA with PCAngsd without specifying the number of eigenvalues
python $PCANGSD -beagle ${NJFOLDER}/${BFILE}.beagle.gz -o ${NJFOLDER}/${BFILE}_${MAFTERM} -minMaf ${MAF} -sites_save -tree -tree_samples ${TREESAMPLES} -threads 8 > NJtree_PCAngsd_${BFILE}.txt
