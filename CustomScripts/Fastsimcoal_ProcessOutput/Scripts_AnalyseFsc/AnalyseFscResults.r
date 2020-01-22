# Vitor Sousa
# 21/01/2020
# This script will automatically process the output 
# of several runs of fastsimcoal2 for a given model.
# You need to specify the settings, where you specify the paths to the results.
# This will work even if your files are in a remote server.

# Load functions
source("utilFscOutput.r")
source("ParFileInterpreter_VS.r")

####################
## Settings       ##
####################

settings <- list()
# population tag 
settings$poptag <- "NCS"
# model tag
settings$modeltag <- "nomig_S2B" 
# population names according to order in Obs SFS
settings$pop.names <- c("NorthOff", "OnSH", "SouthOff")
# path to folder with results
# here if your results are in a server you could use
# username@server.address:/folderInServer
settings$pathtofolder <- paste("../RunFastsimcoal/",settings$poptag,"-",settings$modeltag,sep="")
# observed SFS file name and path used to get maximum likelihood estimates
# this can contain linked SNPs
settings$obsfilename <- paste("../RunFastsimcoal/1000bp_dist2_ind3_4_5_MAF_nomon_DSFS.obs",sep="")
# observed SFS with only independent SNPs
settings$obsfilename_unlinkedSNPs <- paste("../RunFastsimcoal/1000bp_dist2_ind3_4_5_1SNP_MSFS.obs",sep="") 
# need an option for multi-SFS
settings$multiSFS <- TRUE
# -C option with minimum SFS counts. All entries with less than -C are pooled together
settings$minentry <- 1


# MAIN STARTS HERE - NO NEED TO CHANGE SETTINGS AFTER THIS

# 1. Create folder to save the runs
folderName <- paste(settings$poptag,settings$modeltag, sep="-")
dir.create(folderName)

# 2. Get the header and parameters with the bestlikelihood
getALLFilesServer(settings)

# 3. Get the bestlikelihoods into correct folder
setwd(paste("./", folderName ,sep=""))
results <- getbestparam(settings)
# look at the best parameter estimates
results$bestlik

# 4. Get the best run files from server into correct folder
setwd("../")
getRUNFilesServer(settings, results)

# 5. Get fit of obs to expected SFS
pdf(file=paste(folderName,"/FitSFS.pdf",sep=""), width=8, height=5)
getfitobsexp(settings, results)
dev.off()

# 6. Call script to make the plot of the model and re-scale parameters
pdf(file=paste(folderName,"/ModelParam.pdf",sep=""), width=6, height=6)
# pdf(file=paste(settings$poptag,"-",settings$modeltag,"_model.pdf", sep=""), width=8, height=10)
maxparfile <- paste("./", folderName ,"/run", results$run2, "/", folderName, "/", folderName, "_maxL", sep="")
parFileInterpreter(args=maxparfile, pop.names=settings$pop.names, gentime=1, printPDF=FALSE)
dev.off()

# 7. Get the likelihood and AIC for dataset with unlinked sites
if(settings$multiSFS) {
  aic_lhood <- computelhoodAIC(settings, results)
  aic_lhood
  
} else {
  stop("Re-computing the likelihood for unlinked sites currently only workds for multidimensional SFS.")
  # TO DO!!!
}



