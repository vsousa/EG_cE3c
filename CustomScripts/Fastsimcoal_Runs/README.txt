# Vitor Sousa
# 20.01.2020 -----------------

# This script will use the script runFsc.sh to call several runs of fastsimcoal.
# REQUIREMENTS:
# - fsc26 : fastsimcoal2 executable
# - SFS.obs : observed SFS
# - model*.tpl : TPL file for a model (or many models)
# - model*.est : EST file for a model (or many models)

# You need to modify runFsc.sh to adjust the fastsimcoal2 settings.
# You will need to specify the following fastsimcoal2 options:
# -n option, maxNumSims=1000	
# -L option, maxNumLoopsInBrentOptimization=10
# -C option, minValidSFSEntry=1
# -c option (Number of cores), numCores=1  
#--Derived allele frequency sprectrum or Minor allele frequency spectrum?
# For minor allele frequency spectrum use "-m", 
# For derived allele frequency spectrum use "-d"
#-- Monomorphic sites?
#useMonoSites=""    #Uncomment this line to use monomorphic sites
#useMonoSites="-0" #Uncomment this line NOT to use monomorphic sites
#--multiSFS?
#multiSFS="" # Uncomment this line and comment next if you do not use the --multiSFS option
# multiSFS="--multiSFS" #--multiSFS
# Once you change these settings, you can analyse your data with
# different models, all using the same settings.

# Make fastsimcoal2 and scripts (ending in *.sh) executable
chmod +x fsc26
chmod +x *.sh

# runFsc.sh requires the following arguments:
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
./runFsc.sh NCS mig_S2B 1000bp_dist2_ind3_4_5_MAF_nomon_DSFS.obs MSFS.obs