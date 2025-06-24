#!/bin/bash

#Sofia L. Mendes
#16.06.2021

# Script to obtain perform PCA using PCAngsd and run NGSadmix

#edit for 125 inds on 20.06.2023

#Folder where we want to perform the PCA and where we have the input files
#Script should be run in this folder
PCAFOLDER=/path/to/folder/where/we/will/save/output
#Name of the beagle.gz input file and the output files
BFILE="merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc"
#MAF filter to apply
MAF=0.01
#MAF term to add to output file name
MAFTERM=MAF001


#Path to PCAngsd
PCANGSD=/home/smendes/programs/pcangsd/pcangsd.py


#Perform a PCA with PCAngsd without specifying the number of eigenvalues
#python $PCANGSD -beagle ${PCAFOLDER}/${BFILE}.beagle.gz -o ${PCAFOLDER}/${BFILE}_${MAFTERM} -minMaf ${MAF} -threads 8 -admix > PCAlog_${BFILE}.txt

#Perform a PCA with PCAngsd based on various numbers of eigenvalues
for E in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16;
do 
python $PCANGSD -beagle ${PCAFOLDER}/${BFILE}.beagle.gz -o ${PCAFOLDER}/${BFILE}_${MAFTERM} -minMaf ${MAF} -threads 5 -admix -e ${E} > PCAlog_${BFILE}_e${E}.txt
done


#Perform a NGSadmix for various K values
for K in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17;
do 
NGSadmix -likes ${PCAFOLDER}/${BFILE}.beagle.gz -K ${K} -o ${PCAFOLDER}/NGSadmix_${BFILE}_${MAFTERM}_${K} -P 5 -minMaf ${MAF} > NGSAdmixlog_${BFILE}_K${K}.txt
done
