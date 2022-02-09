#!/bin/bash

# SETTINGS
file="sim_divsel3_scaled_template.s" # name of file with template of model
folder="sawfly_divsel3_500Kb_scaled" # name of the folder where all results are saved
tag="sawfly_divsel3" # tag added to beggining of each file
nruns=1000;  # Number of runs - change it later

# info about ancestral pop
tag_anc="anc_saw"; # tag added to beggining of each file
folderanc=Ancestral_500Kb; # folder of ancestral pop

# info about genome
seqsize=499999; # Number of sites in sequence
posmut=249999; # position of beneficial mutation 

# info about demography
pop1size=226; # size of population 1
pop2size=752; # size of population 2
mutrate="3.5e-7"; # mutation rate
ratiom01m10=0.0468; # ratio of migration m01 to m10

# current folder
currentfolder="$PWD";

# check that folder exists and create folder in current directory
if [ -d ${folder} ]; then
	echo "Folder ${folder} already exists.";
else 
	mkdir ${folder};
fi

cd ${folder}
cp ${currentfolder}/slim  .
cp ${currentfolder}/$file .

# All the parameter values are given for the X chromosome (haplodiploids)
for recrate in 3.5e-7 1.05e-6 # recombination rate
do
	for chrm in A X # chromosome (A-autosome, X-x chromosome)
	do
		for freqmut in 0.01 0.1 # initial frequency in number of chromosomes
		do 
			for migrate in 3.65e-4 # migration rate
			do 
				for selmutben in 0.000 0.0002 0.0004 0.0009 0.0030 0.0042 0.0058 0.0081 0.0113 0.0156 0.0217 0.0302 0.0420 0.0583 0.0810 0.113 0.156 0.217 0.302 # selective coefficient (as a function of the Ne of N. pinetum)				
				do 
					for selmutdel in 0.000 # selective of deleterious muts
					do
						for sexratio in 0.3 # sex-ratio
						do
							for dominanceben in 0.01 0.50 # dominance of beneficial mutations
							do 						
								for dominancedel in 0.5 # dominance of deleterious mutations
								do 
									for timeendrec in 500 # time of end of simulations (split time=tend-8000)
									do						
										for timeendold in 1549 # time of end of simulations (split time=tend-8000)
										do						
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
												tendr=${timeendrec};
												tendo=${timeendold};
												
												# tag of file with ancestral data
												ancestral="${tag_anc}_${chrm}_sr${sexratio}_r${recrate}_10000";

												# change settings here
												sed "s/TAGFILEEND/\"$outfile\"/g" $file > $outfile.s;
												sed -i "s/TAGFILESAVE/\"${out_recent}\"/g" $outfile.s;
												sed -i "s/ANCFILE/\"${ancestral}\"/g" $outfile.s;
												sed -i "s/FOLDERANC/${folderanc}/g" $outfile.s;
												sed -i "s/CHR/\"${chrm}\"/g" $outfile.s;
												sed -i "s/SEXRATIO/${sexratio}/g" $outfile.s;
												sed -i "s/MIGRATE/${migrate}/g" $outfile.s;
												sed -i "s/SELMUTDEL/${selmutdel}/g" $outfile.s;
												sed -i "s/SELMUTBEN/${selmutben}/g" $outfile.s;
												sed -i "s/DOMINANCEDEL/${dominancedel}/g" $outfile.s;
												sed -i "s/DOMINANCEBEN/${dominanceben}/g" $outfile.s;
												sed -i "s/FREQ/${freqmut}/g" $outfile.s;
	

												# change these less often...									
												sed -i "s/MUTRATE/${mutrate}/g" $outfile.s;
												sed -i "s/RECRATE/${recrate}/g" $outfile.s;
												sed -i "s/POP1SIZE/${pop1size}/g" $outfile.s;
												sed -i "s/POP2SIZE/${pop2size}/g" $outfile.s;
												sed -i "s/RATIOM01M10/${ratiom01m10}/g" $outfile.s;
												sed -i "s/SAVEGEN/${tendr}/g" $outfile.s;
												sed -i "s/ENDGEN/${tendo}/g" $outfile.s;
												sed -i "s/POSMUT/${posmut}/g" $outfile.s;
												sed -i "s/SEQSIZE/${seqsize}/g" $outfile.s;
												
	
												# Run SLIM2 to perform simulations in parallel
												#seq 10 | parallel -n0 -q ls ${outfile}.s
												seq 1 $nruns | parallel -j 6 ./slim -s 9{}51{}2 -d "rep={}" ${outfile}.s > ${outfile}.log
	
												# Read ms files and put them together into a single one
												filems_old="ms_${outfile}_all.txt";
												filems_rec="ms_${out_recent}_all.txt";
												#filems_female="ms_${outfile}_neutral.txt"
												for((rep=1; rep<$nruns+1; rep++))
												do
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
												# Read fitness files and put them together into a single one
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
												# Read fitness files and put them together into a single one
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




