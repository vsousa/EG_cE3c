Vitor Sousa
20.01.2020 -----------------

# Make fastsimcoal2 and scripts (ending in *.sh) executable
chmod +x fsc26
chmod +x *.sh

# Modify the file runFsc.sh to have the general settings defined for any model.
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

./runFsc.sh NCS nomig_S2B 1000bp_dist2_ind3_4_5_MAF_nomon_DSFS.obs MSFS.obs