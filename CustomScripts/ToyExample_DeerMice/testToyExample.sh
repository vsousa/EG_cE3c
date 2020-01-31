#!/bin/bash
# Vitor Sousa, last updated 24/01/2020 

echo "IMPORTANT: Run line by line! Do not run as ./testToyExample.sh";

##############################################################
##############################################################
# 1. Go from VCF to SFS
##############################################################
##############################################################

# Copy the files from the repository
cp -r ../Fastsimcoal_VCFtoSFS/* .

##############################################################
# 1.1 Modify the ProcessVCF.sh file with the following settings
##############################################################

# ***************************************************************** 
# EDIT MANUALLY ProcessVCF.sh 
# Copy the following
# ***************************************************************** 
# tag for VCF file (vcf file with format "vcffile".vcf)
vcffile="filtered_3scaf";
# tag for the resulting files
vcftag=final_${vcffile};
# tag for file with ind pop info (i.e. information about the individuals and corresponding pop)
# this can have just a sub-set of individuals and populations of the VCF.
indidpop="IndPopInfo_OffOnSH.txt";
# tag for folder with output data
outSFSfolder="block_SFS";
# block size in bp
block_length=1000;
# minimum sample size per pop (must be given as a string separated by comma)
# order of pops must be the same as in indidpop file!
ind_threshold="1,2,3";
# minimum median distance between consecutive SNPs in a good block
# blocks with a median distance equal or lower than this are discarded.
dist_threshold=2;
# path to folder with Rscripts
scriptsRfolder="./Scripts_VCFtoSFS";
# random sampling of individuals? FALSE for deterministic sampling (always selecting the top individuals ranked according to content of data)
randomInd=F;
# seed for random sampling
seed=6126151;
# *****************************************************************

# Make the bash script executable
chmod +x ProcessVCF.sh

##################################################
# 1.2. Run the script to get the SFS based on VCF
##################################################
./ProcessVCF.sh

# tag for files and folders
tagfilefolder=${block_length}bp_dist${dist_threshold}_ind${ind_threshold//\,/\_};

# This will create the folder with the observed SFS
ls block_SFS_${tagfilefolder}

##############################################################
##############################################################
# 2. Run fastsimcoal 2
##############################################################
##############################################################

# Copy folder to run fastsimcoal2
cp -r ../Fastsimcoal_Runs/runFsc.sh .
chmod +x *.sh

# Copy the observed SFS with all linked sites 
# This is used to estimate parameters and find the expected SFS for each model
cp ./block_SFS_${tagfilefolder}/SFS_${tagfilefolder}/${tagfilefolder}_MAF_nomon_DSFS.obs .

# Make fastsimcoal2 and scripts (ending in *.sh) executable
chmod +x fsc26

#################################################################################
# 2.1. Modify the file runFsc.sh to have the general settings defined for any model.
#################################################################################

# ***************************************************************** 
# EDIT MANUALLY runFsc.sh 
# Copy the following
# ***************************************************************** 
# ##--- Name of the fastsimcoal2 executable
fsc=fsc26 
echo "Fastsimcoal version ${fsc}"

##--- Name of folder where the error and output log files will be saved
msgs=conOutput
echo "Error and log files saved in folder ${msgs}"

##--- Number of runs (this can be use to keep adding runs for each model)
# NOTE: You should have 100 runs!! Here 3 just to test.
runBase=1   # initial run number
numRuns=3 # final run number

##--- Fastsimcoal related parameters
# -n option
# NOTE: You should use 100000 to 1000000. Here 10000 just to be quick to test.
maxNumSims=10000	
# -L option
# NOTE: You should have at least 40 cycles. Here 10 just to be quick to test.
maxNumLoopsInBrentOptimization=10
# -C option
minValidSFSEntry=1
# -c option (Number of cores)
numCores=1  

#--Derived allele frequency sprectrum or Minor allele frequency spectrum?
# For minor allele frequency spectrum use "-m", 
# For derived allele frequency spectrum use "-d"
SFStype="-m"
#SFStype="-d"

#-- Monomorphic sites?
#useMonoSites=""    #Uncomment this line to use monomorphic sites
useMonoSites="-0" #Uncomment this line NOT to use monomorphic sites
quiet="-q" #-q option
	
#--multiSFS?
#multiSFS="" # Uncomment this line and comment next if you do not use the --multiSFS option
multiSFS="--multiSFS" #--multiSFS
# ***************************************************************** 



##############################################################
# 2.2. runFsc.sh
##############################################################
# It requires the following arguments
# - poptag: tag of the populations analysed
# - tplEstTag: tag of the model, i.e. tag of the EST and TPL file 
#   NOTE: the resulting folders and file will be named tplEstTag
# - obsSFSfileTag: tag of the name of the file with the observed SFS 
#                  This can be anything, it does not need to be in the name format required by fastsimcoal2)
#                  NOTE: if you have multiple pairwise 2D SFS files, they are required to have 1_0, 2_0, 2_1, etc. indicating the pairwise comparison
# - obsFileEnding: tag for the ending of the SFS file according to fastsimcoal2 requirements:
#                  "DSFS.obs" for multiSFS derived allele
#                  "MSFS.obs" for multiSFS MAF
#                  "jointDAFpop1_0.obs" for 2D derived allele
#                  "jointMAFpop1_0.obs" for 2D MAF
#                  "DAFpop0.obs" for 1D derived allele
#                  "MAFpop0.obs" for 1D MAF
./runFsc.sh NCS nomig_S2B ${tagfilefolder}_MAF_nomon_DSFS.obs MSFS.obs

# This will create the following folders
ls NCS-nomig_S2B
ls conOutput_NCS-nomig_S2B

# Check if the error files are zero-byte file, i.e. that they are empty
ls -alth conOutput_NCS-nomig_S2B/*.err

##############################################################
##############################################################
# 3. Process the output of fastsimcoal2
##############################################################
##############################################################

# copy the scripts to analyse output of fastsimcoal2
mkdir Scripts_AnalyseFsc
cp -r ../Fastsimcoal_ProcessOutput/Scripts_AnalyseFsc/* ./Scripts_AnalyseFsc
ls Scripts_AnalyseFsc
cd Scripts_AnalyseFsc

####################################################
# 3.1 Modify settings of file: AnalyseFscResults.r
####################################################

# ***************************************************************** 
# EDIT MANUALLY AnalyseFscResults.r
# Copy the following
# ***************************************************************** 
settings <- list()
# population tag 
settings$poptag <- "NCS"
# model tag
settings$modeltag <- "nomig_S2B" 
# population names according to order in Obs SFS
settings$pop.names <- c("SouthOff", "OnSH", "NorthOff")
# path to folder with results
# here if your results are in a server you could use
# username@server.address:/folderInServer
settings$pathtofolder <- paste("../",settings$poptag,"-",settings$modeltag,sep="")
# observed SFS file name and path used to get maximum likelihood estimates
# this can contain linked SNPs
settings$obsfilename <- paste("../block_SFS_1000bp_dist2_ind1_2_3/SFS_1000bp_dist2_ind1_2_3/1000bp_dist2_ind1_2_3_MAF_nomon_DSFS.obs",sep="")
# observed SFS with only independent SNPs
settings$obsfilename_unlinkedSNPs <- paste("../block_SFS_1000bp_dist2_ind1_2_3/SFS_1000bp_dist2_ind1_2_3/1000bp_dist2_ind1_2_3_1SNP_MSFS.obs",sep="") 
# need an option for multi-SFS
settings$multiSFS <- TRUE
# -C option with minimum SFS counts. All entries with less than -C are pooled together
settings$minentry <- 1
# ***************************************************************** 

# 3.2. Run the script AnalyseFscResults.r or run the Script line by line in RStudio
Rscript AnalyseFscResults.r

######################################################################################
######################################################################################
# 4. Check the output and files and delete everything to keep just the required files
######################################################################################
######################################################################################

#This folder should contain the following files:
#- filtered_3scaf.vcf: VCF file with deer mice data from Pfeifer et al (2018) MBE.
#- fsc26: fastsimcoal2 executable (Linux)
#- IndPopInfo_OffOnSH.txt: file with the list of individuals in the order we want to get in the SFS
#- nomig_S2B.est: example EST file for fastsimcoal2 run
#- nomig_S2B.tpl: example TPL file for fastsimcoal2 run
#- testToyExample.sh: bash script to run
#- README_ToyExample.txt: this file
