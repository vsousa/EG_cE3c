#!/bin/bash
#SBATCH --job-name=SFS2D
#SBATCH --mem=5G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -p fct
#SBATCH -q cpca158582022
#SBATCH --array=0-1 # Specify here the number of pairwise combinations


#Sofia L. Mendes
#10.10.2023
#Vitor Sousa 18.10.2023 - adjusted to add SLURM commands
# Script to run ANGSD to get the pairwise 2D-SFS for all pairwise combinations of saf files in a folder
# The array length should be the number of pairwise comparisons
# Requires the SAF file for each populations
# Example of usage:
#       sbatch sfs2d_test_slurm.sh

# Get the 2D-SFS for all pairs of populations

#Folder where we want to save the outputs of the calculation
SFSFOLDER=./SFS

#Folder where the saf files are located
SAFFOLDER=./SAF;

#Folder where we want to save the outputs of the FST 
FSTFOLDER=./FST

# Check that folder SAFFOLDER exist, and create it if it does not
if [ -d "${SAFFOLDER}" ] 
then
    echo "Computing 2D-SFS. Looking for *.saf files in folder: ${SAFFOLDER}"
else
    echo "ERROR: Folder ${SAFFOLDER} with *.saf files not found."
fi

# Check that the SFSFOLDER exists
if [ -d "${SFSFOLDER}" ] 
then
    echo "Saving SFS files in folder: ${SFSFOLDER}"
else
    mkdir ${SFSFOLDER}
    echo "Saving SFS files in folder: ${SFSFOLDER}"
fi

# Check that folder FSTFOLDER exists
if [ -d "${FSTFOLDER}" ] 
then
    echo "FST results will be saved in folder: ${FSTFOLDER}"
else
    echo "ERROR: FST folder ${FSTFOLDER} not found."
fi

# get the list of files with list of bamfiles for each population
# i.e. files that end with .filelist
SAFFILES=(./${SAFFOLDER}/*.saf.gz)

# get vector with pop names, need to loop through file names
for index in ${!SAFFILES[@]}
do
  POPNAME[$index]=$(basename ${SAFFILES[$index]} .saf.gz) # use basename to remove the .saf.gz file extension
  echo ${POPNAME[$index]}
done

# Create variables that contain the combination of populations
COUNT=0
let ENDi=${#POPNAME[@]}-2
let ENDj=${#POPNAME[@]}-1
for i in $(seq 0 ${ENDi})
do
  let iplus=${i}+1
  for j in $(seq ${iplus} ${ENDj})
  do
      echo ${i} ${j} # combination of index of populations
      POPi[$COUNT]=${POPNAME[${i}]} # save the population 1 in POPi
      POPj[$COUNT]=${POPNAME[${j}]} # save the population 2 in POPj
      ((COUNT++)) # increase the counter
  done
done

# STEP 2 - Get the 2D-SFS for each combination of populations, running each combination in parallel
realSFS ./${SAFFOLDER}/${POPi[$SLURM_ARRAY_TASK_ID]}.saf.idx ./${SAFFOLDER}/${POPj[$SLURM_ARRAY_TASK_ID]}.saf.idx -P 1 > ./${SFSFOLDER}/${POPi[$SLURM_ARRAY_TASK_ID]}_${POPj[$SLURM_ARRAY_TASK_ID]}.ml

# STEP 3 - prepare the fst for easy window analysis etc
realSFS fst index ${SAFFOLDER}/${POPi[$SLURM_ARRAY_TASK_ID]}.saf.idx ${SAFFOLDER}/${POPj[$SLURM_ARRAY_TASK_ID]}.saf.idx -sfs ${SFSFOLDER}/${POPi[$SLURM_ARRAY_TASK_ID]}_${POPj[$SLURM_ARRAY_TASK_ID]}.ml -P 1 -fstout ${FSTFOLDER}/${POPi[$SLURM_ARRAY_TASK_ID]}_${POPj[$SLURM_ARRAY_TASK_ID]}

# STEP 4 - get mean FST for each pair of populations
realSFS fst stats ${FSTFOLDER}/${POPi[$SLURM_ARRAY_TASK_ID]}_${POPj[$SLURM_ARRAY_TASK_ID]}.fst.idx -P 1 > ${FSTFOLDER}/${POPi[$SLURM_ARRAY_TASK_ID]}_${POPj[$SLURM_ARRAY_TASK_ID]}.meanfst


