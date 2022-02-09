# Vitor Sousa
# updated 01/02/2021 
# removed comments and cleaned.
# 20/07/2020
# This script will compute the likelihood and AIC
# of several models, based on the multiSFS of independent SNPs.
# Note that if you are working of pairs of 2D-SFS, 
# even if using just independent SNPs,
# the likelihood will still be a composite likelihood 
# because the pairs of 2D SFS are not independent.
# NOTE 1: Thus, this only works for joint multiSFS.
# NOTE 2: This only works for the observed SFS without monomorphic sites.
# You need to specify the settings, where you specify:
# - the paths to the maxPar files of the different models.
# - the paths to the bootstrap observed multiSFS replicates.
# The script will get the maxPar files for different models.
# For each model, it will obtain the expected SFS.
# Then, it computes the likelihood for the different bootstrap SFS
# based on the expected SFS.
# This will work even if your files are in a remote server.
# Note that this only works

# Load functions
source("./Scripts_VCFtoSFS/utilFscOutput.r")
source("./Scripts_VCFtoSFS/ParFileInterpreter_VS.r")

####################
## Settings       ##
####################

settings <- list()
# model tags
settings$modeltag <- c("D_Hybrid","E_SC_TG","F_SC_CT")
# model population tag 
settings$poptag <- c("COCRCANPs","COCRCANPs","COCRCANPs")
# model code - the name of each model is given by poptag-code 
settings$code <- c("a_ncm_hs","a_sc_ct","a_sc_gt")
# number of parameters of each model
settings$nparam <- c(9, 9, 9)
# individual sample size for each pop (in number of individuals)
settings$ss <- c(1,2,1)
# path to maxPar files
# here if your results are in a server you could use
# username@server.address:/folderInServer
settings$pathmaxpar <- c("../../FSC_results_manuscript_08_07_2020/best_runs_02_06_2020/D_HybridSpeciation/run38",
                         "../../FSC_results_manuscript_08_07_2020/best_runs_02_06_2020/E_CT_SecContact_TG/run39",
                         "../../FSC_results_manuscript_08_07_2020/best_runs_02_06_2020/F_TG_SecContact_CT/run15")
# path to observed SFS bootstrap replicates 
settings$bootpath <- paste("./block_SFS_200bp_dist2_ind1_2_1/Bootstrap_1SNP_200bp_dist2_ind1_2_1/jsfs",sep="")
# tag for observed SFS file 
settings$bootfile <- paste("sfs_1snp_200bp_dist2_ind1_2_1_boot_nm_MSFS",sep="")
# number of bootstrap replicates
settings$nboot <- 1000
# path to observed SFS multisfs
settings$obsfile <- paste("./block_SFS_200bp_dist2_ind1_2_1/SFS_200bp_dist2_ind1_2_1/200bp_dist2_ind1_2_1_MAF_nomon_DSFS.obs")
# -C option with minimum SFS counts. All entries with less than -C are pooled together
settings$minentry <- 1
# path to faststimcoal2 executable
settings$fsc <- paste("./fsc26")
# number of coalescent simulations to get the expected SFS
settings$ncoalsim <- 1e5
# number of replicates to get the variation in the approximation of expected SFS
settings$nreplicate <- 10
# maf or daf tag "d" for daf, "m" for maf
settings$mdaf <- "m"

##############################################################
# MAIN
##############################################################
# you do not need to change anything else after this line.

# create folder to run the simulations
folderSaveSFS <- "ExpectedMultiSFS"
dir.create(folderSaveSFS)
setwd(paste("./",folderSaveSFS,sep=""))

# Create terminal to run fastsimcoal2 from RStudio session (this works for windows)
termID <- rstudioapi::terminalCreate("fsc", shellType = "win-wsl-bash")

############################################################################
# 1. Simulate the expected SFS for each model according to maxL.par file
############################################################################

for(i in 1:length(settings$pathmaxpar)) {
  print(paste("parameters",i))
  # if the terminal is not running, start the new iteration of the for loop
  folder_expsfs <- paste(settings$poptag[i],"-", settings$code[i], "_maxL",sep="")
  # copy observed file name to folder where simulations are done
  # changing the file name
  if( settings$mdaf == "m") {
    cmd_copyobs <- paste("cd ", folderSaveSFS, ";\n cp .", settings$obsfile, " ", settings$poptag[i],"-", settings$code[i], "_maxL_MSFS.obs;\n", sep="")
    file_expsfs <- paste(settings$poptag[i],"-", settings$code[i],"_maxL_MSFS.txt", sep="")  
    save_expfsf <- paste(settings$poptag[i],"-", settings$code[i], "_MSFS.txt", sep="")
  } else {
    cmd_copyobs <- paste("cd ", folderSaveSFS, ";\n cp .", settings$obsfile, " ", settings$poptag[i],"-", settings$code[i], "_maxL_DSFS.obs;\n", sep="")
    file_expsfs <- paste(settings$poptag[i],"-", settings$code[i],"_maxL_DSFS.txt", sep="")  
    save_expfsf <- paste(settings$poptag[i],"-", settings$code[i], "_DSFS.txt", sep="")
  }
  # run the command lines
  rstudioapi::terminalSend(termID, cmd_copyobs)
  # fastsimcoal2 options to run
  cmd_runfsc <- paste("for j in $(seq 1 ", settings$nreplicate,"); \n
                       do \n
                          echo Simulating expected SFS $j;\n
                          ../fsc26 -i ../", settings$pathmaxpar[i],"/", settings$poptag[i],"-", settings$code[i], "/", settings$poptag[i],"-", settings$code[i],"_maxL.par -q --multiSFS -c2 -B2 -n ", format(settings$ncoalsim, scientific=FALSE)," -", settings$mdaf, " ;\n
                       if [ $j -eq 1 ];\n
                          then\n
                               cat ./", folder_expsfs, "/", file_expsfs, " > ", save_expfsf, ";\n
                          else \n
                               tail -n 1 ./", folder_expsfs, "/", file_expsfs, " >> ", save_expfsf, ";\n
                          fi  \n
                       done \n
                      rm -r ", folder_expsfs, "; \n
                      cd .. \n", sep="")  
  #print(cmd_runfsc)
  # run the command lines
  rstudioapi::terminalSend(termID, cmd_runfsc)  
}

# kill the terminal
rstudioapi::terminalKill(termID)



############################################################################
# 2. Read the exppected SFS for each model 
############################################################################

# read expSFS
expSFS <- lapply(c(1:length(settings$pathmaxpar)), function(i) {
  # read the expected SFS of the first model
  #folder_expsfs <- paste(settings$poptag[i],"-", settings$code[i], "_maxL",sep="")
  if(settings$mdaf=="m") {
    file_expsfs <- paste(settings$poptag[i],"-", settings$code[i],"_MSFS.txt", sep="")  
  } else {
    file_expsfs <- paste(settings$poptag[i],"-", settings$code[i],"_DSFS.txt", sep="")
  }
  # read the expected SFS
  matrix(scan(paste(file_expsfs, sep="/"), skip=2, nlines=settings$nreplicate), ncol = settings$nreplicate)
})
setwd("../")

# check that all expSFS have the same size and that it is the same as prod of sample size
sfssize <- unlist(lapply(expSFS, function(x) nrow(x)))
if(length(unique(sfssize))!=1 ) {
  stop("Not all expected SFS have the same size! Check the maxL.para file for each model.")
}

# check that all expSFS are equal to the 
if(sum(sfssize!=prod((settings$ss*2)+1))>0) {
  stop("Expected SFS do not have the expected size! Check the sample size (settings$ss) and number of lineages in maxL.par file.")
}

################################################################################################
# 3. Read the 1SNP sfs bootstrap replicates and compute corresponding likelihood for each model
################################################################################################

# read the observed SFS
obsboot <- sapply(1:settings$nboot, function(i) {
  #print(paste("reading boot ", i))
  # read the bootstrap obs SFS
  scan(paste(settings$bootpath,"/",settings$bootfile, i, ".obs",sep=""), skip=2)
})

summary(colSums(obsboot))

# distribution of the number of SNPs per entry
hist(obsboot[obsboot>0], xlim=c(0,1000), breaks=1000)

# 40 seems to recap the results with all the SNPs
settings$minentry <- 1
hist(apply(obsboot, 2, function(x) sum(x>settings$minentry)))

# compute the likelihood
lhood <- sapply(1:settings$nboot, function(i) {
  
  # # compute the likelihood according to the models
  # model_lhood_tmp <- sapply(expSFS, function(x) {
  #   apply(x,2,function(col) computelhood(obsboot[,i], col, settings$minentry))
  # })
  # 
  # # compute the mean likelihood for each model (among the replicates performed)
  # model_lhood <- apply(model_lhood_tmp, 2, mean)
  # 
  # option 2
  # compute the mean SFS and get the likelihood based on mean SFS
  mean_expSFS <- sapply(expSFS, function(x) {
    rowMeans(x)
  }) 
  model_lhood <- apply(mean_expSFS, 2, function(col) {
    computelhood(obsboot[,i], col, settings$minentry)
  })
  
  # compute the delta likelihood
  maxobs_lhood <- computelhood(obsboot[,i], obsboot[,i]/sum(obsboot[,i]), settings$minentry)
  delta_lhood <- model_lhood - maxobs_lhood
  
  # output a vector with the likelihood and delta_likelihood
  c(model_lhood,delta_lhood)
})

boxplot(t(lhood[9:16,]))

# compute the AIC for each model
aic <- apply(lhood, 2, function(col) {
  AIC(col[1:length(settings$nparam)],settings$nparam)
})

# compute the relative likelihood
# accoding to formula of Excoffier et al 2013
rellhood <- apply(aic, 2, function(col) {
  delta_aic <- col-min(col)
  lhood <- exp(-0.5*delta_aic)
  denominator <- sum(lhood)
  lhood/denominator
})

boxplot(t(rellhood))

table(unlist(apply(rellhood, 2, function(col) {which(col==max(col))})))

tmp <- unlist(apply(rellhood, 2, function(col) {which(col>0.95)}))
str(tmp)
barplot(table(tmp))

tmp <- data.frame(model=rep(settings$modeltag, each=settings$nboot), 
                  rellhood=as.numeric(t(rellhood)),
                  sim=rep(1:settings$nboot, times=length(settings$modeltag)))
library(ggplot2)

# pdf(file=paste("relLhood_1snp_",paste(settings$ss,collapse=""),"_nboot",settings$nboot,"_minentry",settings$minentry,".pdf",sep=""), width=10, height=6)
pdf(file=paste("relLhood_3modelsA-C_1snp_",paste(settings$ss,collapse=""),"_nboot",settings$nboot,"_minentry",settings$minentry,".pdf",sep=""), width=6, height=4)
ggplot(data=tmp, aes(x = model, y = rellhood)) +  
  geom_boxplot(aes(fill = model), alpha = 0.5, outlier.alpha = 0.05, outlier.colour = NULL, outlier.fill = NULL) 
dev.off()  

# TO DO:
# ggplot with distribution of delta likelihood
# repeat just for the 3 more complex models. 
# it is clear that all models lead to overlaping lhoods, except for model 3.


