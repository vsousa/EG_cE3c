#!/bin/bash

# SLURM THINGS HERE

# SETTINGS
file="sim_divsel3_scaled_template.s" # name of file with template of model
folder="divsel3_500Kb_scaled_extra" # name of the folder where all results are saved
tag="divsel" # tag added to beggining of each file
nruns=1000;  # Number of runs - change it later

# current folder
currentfolder="$PWD";

# check that folder exists and create folder in current directory
if [ -d ${folder} ]; then
	echo "Folder ${folder} already exists.";
else 
	mkdir ${folder};
fi

# go to folder and copy relevant files
cd ${folder}
cp ${currentfolder}/slim  .
cp ${currentfolder}/$file .

# Loop through all the parameter combinations
for recrate in 2.5e-7 2.5e-8 # recombination rate
do
	for chrm in A X # chromosome (A-autosome, X-x chromosome)
	do
		for freqmut in 0.1 0.5 # initial frequency in number of chromosomes
		do 
			for migrate in 0.000 0.00034 0.0017 0.0034 # migration rate
			do 
				for selmutben in 0.000 0.0067 0.01333 0.02667 0.05333 0.0667 0.1334 # selective coefficient
				do 
					for selmutdel in 0.000 # selective of deleterious muts
					do
						for sexratio in 0.5 # sex-ratio
						do
							for dominanceben in 0.01 0.5 # dominance of beneficial mutations
							do 						
								for dominancedel in 0.5 # dominance of deleterious mutations
								do 
									for timeendrec in 500 # time split recent t1
									do						
										for timeendold in 2000 # time split t2
										do						
											# file name tag				
											# create a file with this name replacing the parameters of the template by the current combination of parameters				
											outfile="${tag}_${chrm}_h${dominanceben}_${dominancedel}_m${migrate}_sr${sexratio}_s${selmutben}_${selmutdel}_r${recrate}_f${freqmut}_${timeendold}"
											out_recent="${tag}_${chrm}_h${dominanceben}_${dominancedel}_m${migrate}_sr${sexratio}_s${selmutben}_${selmutdel}_r${recrate}_f${freqmut}_${timeendrec}"
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
												tendr=${timeendrec}; # time of split time 1
												tendo=${timeendold}; # time of split time 2
												mutrate="2.5e-7"; # mutation rate

												# tag of file with ancestral data
												ancestral="anc_${chrm}_sr${sexratio}_r${recrate}_10000";

												# replace parameter values in the template file
												sed "s/TAGFILEEND/\"$outfile\"/g" $file > $outfile.s;
												sed -i "s/TAGFILESAVE/\"${out_recent}\"/g" $outfile.s;
												sed -i "s/ANCFILE/\"${ancestral}\"/g" $outfile.s;
												sed -i "s/CHR/\"${chrm}\"/g" $outfile.s;
												sed -i "s/SEXRATIO/${sexratio}/g" $outfile.s;
												sed -i "s/MIGRATE/${migrate}/g" $outfile.s;
												sed -i "s/SELMUTDEL/${selmutdel}/g" $outfile.s;
												sed -i "s/SELMUTBEN/${selmutben}/g" $outfile.s;
												sed -i "s/DOMINANCEDEL/${dominancedel}/g" $outfile.s;
												sed -i "s/DOMINANCEBEN/${dominanceben}/g" $outfile.s;
												sed -i "s/FREQ/${freqmut}/g" $outfile.s;									
												sed -i "s/MUTRATE/${mutrate}/g" $outfile.s;
												sed -i "s/RECRATE/${recrate}/g" $outfile.s;
												sed -i "s/POPSIZEANC/${popsize}/g" $outfile.s;
												sed -i "s/POPSIZEBOT/${popsize}/g" $outfile.s;
												sed -i "s/POPSIZEDESC/${popsize}/g" $outfile.s;
												sed -i "s/SAVEGEN/${tendr}/g" $outfile.s;
												sed -i "s/ENDGEN/${tendo}/g" $outfile.s;
	
												# Run SLiM to perform simulations in parallel
												seq 1 $nruns | parallel -j 3 ./slim -s 9{}51{}2 -d "rep={}" ${outfile}.s > ${outfile}.log
	
												# Read ms files and put them together into a single one
												filems_old="ms_${outfile}_all.txt";
												filems_rec="ms_${out_recent}_all.txt";
												# go through all the runs
												for((rep=1; rep<$nruns+1; rep++))
												do
													# run slim in a for loop
													#./slim -s 9${rep}51${rep}2 -d "rep=${rep}" ${outfile}.s > ${outfile}.log
												
													# put all ms files together for old divergence time
													if [ -e "ms_${outfile}_${rep}.txt" ]; then
											    			cat "ms_${outfile}_${rep}.txt" >> $filems_old;
														rm "ms_${outfile}_${rep}.txt";
													else 
														echo "File ms_${outfile}_${rep}.txt does not exist"
													fi 
													# put all ms files together for recent divergence time
													if [ -e "ms_${out_recent}_${rep}.txt" ]; then
											    			cat "ms_${out_recent}_${rep}.txt" >> $filems_rec;
														rm "ms_${out_recent}_${rep}.txt";
													else 
														echo "File ms_${out_recent}_${rep}.txt does not exist"
													fi 
												done
												# Read trajectory files and put them together into a single one
												file_traj="traj_${outfile}.benmut"
												for((rep=1; rep<$nruns+1; rep++))
												do
													# put all files with male data together
													if [ -e "${outfile}_${rep}.benmuts" ]; then
											    			cat "${outfile}_${rep}.benmuts" >> $file_traj;
														rm "${outfile}_${rep}.benmuts";
													else 
														echo "File ${outfile}_${rep}.benmuts not found!"
													fi 
												done
												# Read mean fitness trajectory files and put them together into a single one
												file_traj="traj_${outfile}.fitness"
												for((rep=1; rep<$nruns+1; rep++))
												do
													# put all files with male data together
													if [ -e "${outfile}_${rep}.meanfit" ]; then
											    			cat "${outfile}_${rep}.meanfit" >> $file_traj;
														rm "${outfile}_${rep}.meanfit";
													else 
														echo "File ${outfile}_${rep}.meanfit not found!"
													fi 
												done
												# Read fitness landscape files and put them together into a single one
												file_traj="traj_${outfile}.fitlandscape"
												for((rep=1; rep<$nruns+1; rep++))
												do
													# put all files with male data together
													if [ -e "${outfile}_${rep}.fitlandscape" ]; then
											    			cat "${outfile}_${rep}.fitlandscape" >> $file_traj;
														rm "${outfile}_${rep}.fitlandscape";
													else 
														echo "File ${outfile}_${rep}.fitlandscape not found!"
													fi 
												done
											fi
										done
									done
								done
							done 
						done
					done
				done
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




