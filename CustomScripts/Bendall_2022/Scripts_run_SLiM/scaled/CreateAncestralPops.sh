#!/bin/bash

# SLURM THINGS HERE

# SETTINGS
file="sim_createAncestral_template.s" # name of file with template of model
folder="Ancestral_500Kb" # name of the folder where all results are saved
tag="anc" # tag added to beggining of each file
nruns=1000;  # Number of runs - change it later

# current folder
currentfolder="$PWD";

# check that folder exists and create folder in current directory
if [ -d ${folder} ]; then
	echo "Folder ${folder} already exists.";
else 
	mkdir ${folder};
fi

# go to correct folder and copy relevant files
cd ${folder}
cp ${currentfolder}/slim  .
cp ${currentfolder}/$file .

# Loop through all the parameter combinations
for recrate in 2.5e-8 2.5e-7 # recombination rate
do
	for chrm in X A # chromosome (A-autosome, X-x chromosome)
	do
		for sexratio in 0.5 # sex-ratio
		do
			for timeend in 10000 # time of end of simulations
				
				# file name tag
				# create a file with this name replacing the parameters of the template by the current combination of parameters
				outfile="${tag}_${chrm}_sr${sexratio}_r${recrate}_${timeend}"
				echo "Running ${outfile}"							
				if [ -e $file ]; then
				# check if this is the first run or other files are present
					if [ `ls -1 ${outfile}* 2>/dev/null | wc -l ` -eq 0 ]; then
					    echo "First run. Not deleting pre-existing files!"
					else
					    echo "Deleting existing files..."
					    rm $outfile*;
					    rm ms_$outfile*;
					fi
	
					# fixed params
					popsize=1000; # population size in number of individuals
					              # note that for X 1000 individuals correspond to 1500 gene copies
					tend=${timeend}; # time of split
					mutrate="2.5e-7"; # mutation rate

					# replace parameter values in the template file
					sed "s/TAGFILE/\"$outfile\"/g" $file > $outfile.s;
					sed -i "s/CHR/\"${chrm}\"/g" $outfile.s;
					sed -i "s/SEXRATIO/${sexratio}/g" $outfile.s;													
					sed -i "s/MUTRATE/${mutrate}/g" $outfile.s;
					sed -i "s/RECRATE/${recrate}/g" $outfile.s;
					sed -i "s/POPSIZEANC/${popsize}/g" $outfile.s;
					sed -i "s/ENDGEN/${tend}/g" $outfile.s;
	
					# Run SLiM to perform simulations in parallel
					seq 1 $nruns | parallel -j 6 ./slim -s 9{}51{}2 -d "rep={}" ${outfile}.s > ${outfile}.log
	
				fi
			done
		done
	done
done
# go back to original folder
#cd ${currentfolder}/${folder}
# move all the files to current folder
#cp -r ${path_torun}${folder}/* .
rm slim
#rm -r ${path_torun}${folder}




