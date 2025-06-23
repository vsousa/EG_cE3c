#!/bin/bash
#SBATCH --job-name=saf
#SBATCH --time=04:00:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p fct
#SBATCH -q cpca097952021
#SBATCH --array=0-59 # This will compute for two populations


#Sofia L. Mendes
#10.10.2023
#Vitor Sousa 13.10.2023 - adjusted to add SLURM commands
# Script to run ANGSD to all the files ending with ".filelist"
# The array length should be the 
# Obtains the SAF file for each population
# Example of usage:
#       sbatch saf_test_slurm.sh

#Folder where we want to save the outputs of the FST 
SAFFOLDER=./SAF

#Folder where the bam files are
BAMFOLDER=./7_indexed_bam_lcWGS_125inds_scephalus;

#Ancestral FASTA file (Squalius cephalus)
ANC=./SCEPHALUS_ref_genome/GCA_022829025.1_ASM2282902v1_genomic.fna

#Chromossomes to use from the reference genome (Squalius cephalus)
#This will allow us to restrict the analysis to the 25 well assembled chromossomes from the Squalius cephalus genome
CHRM=./SCEPHALUS_ref_genome/chromossome_names_scephalus.txt

#Filters we want to apply
#Mapping quality, base quality, uniquely mapping reads etc (does not vary between populations)
FILTERS="-minMapQ 30 -minQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0"

# get the list of files with list of bamfiles for each population
# i.e. files that end with .filelist
FILES=(*.filelist)
   
# vector with pop names, need to loop through file names
for index in ${!FILES[@]}
do
    POPNAME[$index]=${FILES[$index]%.*} # remove the extension of each file
    echo ${POPNAME[$index]}
done

# Get the number of individuals per pop (equal to the number of lines in each file of FILES)
#Number of individuals required to have data
for index in ${!FILES[@]}
do
    # get the number of lines, removing the first line that has the pop name
    let nlines=$(wc -l < ${FILES[$index]})-1;
    # save the number of lines into vector
    minIND[$index]="-minInd ${nlines}"
    echo ${minIND[$index]}
done

# Print the list of filters used
echo "Running ANGGSD with the following filters: \n ${FILTERS}"

# Running 
echo "Runing the pop ${POPNAME[$SLURM_ARRAY_TASK_ID]} of task id $SLURM_ARRAY_TASK_ID"

#Step 1 - OBTAIN SAF FILE FOR EACH POP
# Loop through all the populations using the ARRAY of SLURM - variable SLURM_ARRAY_TASK_ID
angsd -b ${FILES[$SLURM_ARRAY_TASK_ID]} -anc $ANC -rf $CHRM -out ${SAFFOLDER}/${POPNAME[$SLURM_ARRAY_TASK_ID]} $FILTERS ${minIND[$SLURM_ARRAY_TASK_ID]} -dosaf 1 -gl 1

