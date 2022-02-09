# SCRIPT TO PROCESS SUMMARY FILES OF MS OUTPUT
# creates summary files with average and quantile of summary statistics 
# (e.g. FST, PI, etc.) along the chromosome on a genome scan at windows of 20Kb analysis.
# These summary files can then be used for functions to plot results

rm(list=ls()) # clean memory
# load required packages
library(RColorBrewer)
library(doParallel)
library(ggplot2)
library(MBA)
library(lubridate)
library(reshape2)
library(colorRamps)
library(scales)
library(grid)
library(gridExtra)
library(gridGraphics)
# register the number of cores for parallel
registerDoParallel(cores=5)

# SETTINGS
# number of simulations
nsims <- 1000
# type of mutations sampled
muttype <- "all"
# length of simulated chromosomes
seqlength <- 500000
# sample size from each population (in number of genome copies, e.g. 20 indicates 10 diploid inds)
ss <- 20
# number of populations
npop <- 2
# sliding window analysis
window.size=20000
slide.size=20000
# label for the runs
extrachrm <- rep(c("A","HD"), times=1)
# TRUE if the name of files indicate the initial frequency
freqTag <- TRUE 


# Get the file tag and create the legend text for each scenario
freqTag <- TRUE
tag <- list()
legtext <- list()
scenario <- 1


# COMBINATION OF PARAMETERS
# migration rate
migrate_v <-c("0.000", "0.00034", "0.0017", "0.0034")
# selection of mutation under divergent selection
selmutben_v <- c("0.000", "0.0067", "0.01333" ,"0.02667", "0.05333", "0.0667", "0.1334")
# recombination rate
recrate_v <- c("2.5e-7")
# initial frequency of allele a at locus under divergent selection
freqmut_v <- c("0.001","0.01","0.1","0.5") # done for a given frequency at a time
# dominance of allele a at locus under divergent selection
dominance_v <- c("0.01","0.5")
# time of split of the two populations
time_v <- c("2000")
# folder can be different for different sims. Use this vector for those cases
folder_v <- c("divsel3_500Kb_scaled_merged",
              "divsel3_500Kb_scaled",
              "divsel3_500Kb_scaled_extra_neutral")

# Get the number of parameter combinations
combnum <- length(migrate_v)*length(selmutben_v)*length(recrate_v)*length(dominance_v)*length(freqmut_v)*length(time_v)*2 

# tag to name of output files
#savetag <- paste("FitM2_h_r",recrate_v[1], "_f",freqmut_v[1],sep="")
savetag <- paste("FitM2_h_r",recrate_v[1],sep="")

# for debug
# uncomment and do not run the for loops
# chrm <- "A"
# migrate <- migrate_v[3]
# selmutben <- selmutben_v[6]
# selmutdel <- c("0.000")
# recrate <- recrate_v[1]
# sexratio <- c("0.5")
# dominanceben <- c("0.01")
# dominancedel <- c("0.5")
# freqmut <- freqmut_v[1]
# timeend <- 2000

# Analyse the summary files for:
# - "sumstat" all sims, including those where the allele a was lost
# - "sumstat_cb" only sims conditional on keeping the beneficial allele a in pop1
# initialize a list that will save the results
meanquantilestat <- list()
count <- 1
for(cond_folder in c("sumstat","sumstat_cb")) {
  # Auxiliary files to save the file tag and legend text for each scenario
  tag <- list()
  legtext <- list()
  
  # Auxiliary files to save the file tag and legend text for each scenario
  tag <- list()
  legtext <- list()
  
  # Go through all combinations of parameters
  scenario <- 1
  for(migrate in migrate_v) {
    for(selmutben in selmutben_v) {
      for(selmutdel in c("0.000")) {
        for(recrate in recrate_v) {
          for(sexratio in c("0.5")) {
            for(dominanceben in dominance_v) {
              for(dominancedel in c("0.5")) {
                for(freqmut in freqmut_v) {
                  for(timeend in time_v) {
                    for(chrm in c("A","X")) {
                      
                      # Get the file tag for the input files
                      if(freqTag) {
                        tag[[scenario]]=paste("divsel_", chrm, "_h", dominanceben,"_", dominancedel,"_m", migrate ,"_sr", sexratio, "_s", selmutben, "_", selmutdel ,"_r", recrate, "_f", freqmut, "_", timeend, sep="")    
                      } else {
                        tag[[scenario]]=paste("divsel_", chrm, "_h", dominanceben,"_", dominancedel,"_m", migrate ,"_sr", sexratio, "_s", selmutben, "_", selmutdel ,"_r", recrate, sep="")  
                      }
                      print(tag[[scenario]])
                      
                      # Create the legend text for each
                      if(chrm == "X") {
                        extrachrm <- "HD"
                      } else if(chrm == "A") {
                        extrachrm <- "D"
                      }
                      
                      # save legend text for each scenario
                      legtext[[scenario]] <- paste(extrachrm, "_m", migrate, "_r=", recrate, "_s=", selmutben, "_h=", dominanceben, "_f=", freqmut ,sep="")
                      print(legtext[[scenario]])
                      
                      scenario <- scenario + 1 
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  # Get the parameters for each scenario
  migrate <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("m", a)],"m"))[2])}))
  recrate <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("^r", a)],"r"))[2])}))
  dominance <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("^h", a)],"h"))[2])}))
  chrm <- unlist(lapply(legtext, function(x) {a <- unlist(strsplit(x, "_"))[1]}))
  selection <- unlist(lapply(legtext, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("^s", a)],"s="))[2])}))
  freq <- unlist(lapply(legtext, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("f", a)],"f="))[2])}))
  
  # get the index of combinations for Autosome (D - Diploid)
  aindex <- which(chrm=="D")
  
  # get the index of combinations for X-chromosome (HD - Hemizygous)
  xindex <- which(chrm=="HD")
  
  # for the neutral case the dominance does not matter and with freq 0.001 and 0.01
  # the neutral sims were only done for for h=0.01 and not h=0.5
  # get the correct file name for the neutral cases with initial freq less than 0.1
  eval_neut <-(selection==0.000 & freq==0.001) | ((selection==0.000 & freq==0.01))
  tag[eval_neut] <- sub(pattern="h0.5", replace="h0.01", x=tag[eval_neut])

  # Define the folder for a given combination of parameters
  folder <- rep(folder_v[2], times=combnum)
  folder[freq==0.1 | freq==0.5] <- folder_v[1]
  folder[eval_neut] <- folder_v[3]
  
  # check that all the files exist
  eval <- numeric(length(tag))
  {
    for(scenario in 1:length(tag)) {
      eval[scenario] <- file.exists(paste("./",folder[scenario], "/", cond_folder,"/window_sumstat_", window.size,"_", tag[scenario],sep=""))
    }
  }
  
  # 1. Read the data
  # if files do not exist show their path and file name
  if(sum(eval==0)>0) {
    for(i in which(eval==0)) {
      print(paste("File does not exist: ./",folder[i], "/", cond_folder,"/window_sumstat_", window.size,"_", tag[i],sep=""))
    }
    stop("Some files do not exist!")
  } 
  
  # summary statistics recorded for each window
  list_sumstat <- c("pi1","pi2","dxy","fst")
  
  # read the data
  {
    # read chromosome wide results
    sumstat_window <- foreach(scenario=1:length(tag)) %dopar% {
      read.table(file=paste("./",folder[scenario], "/", cond_folder,"/window_sumstat_", window.size,"_", tag[scenario],sep=""), stringsAsFactors = FALSE, header=T)
    }
  }
  
  # check that all the datasets have the correct number of simulations
  if(cond_folder=="sumstat") {
    if(sum((unlist(lapply(sumstat_window, ncol))/length(list_sumstat))!=nsims)>0) {
      stop("Not all scenarios have the correct number of sims!")
    }  
  }
  
  # Plot the statistics as a function of the migration rate
  migrate <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("m", a)],"m"))[2])}))
  recrate <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("^r", a)],"r"))[2])}))
  dominance <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("^h", a)],"h"))[2])}))
  chrm <- unlist(lapply(legtext, function(x) {a <- unlist(strsplit(x, "_"))[2]}))
  selection <- unlist(lapply(legtext, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("^s", a)],"s="))[2])}))
  freq <- unlist(lapply(legtext, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("f", a)],"f="))[3])}))
  
  
  # Get a summary of the mean and quantiles for each summary statistic 
  meanquantilestat[[count]] <- foreach(i=1:length(list_sumstat)) %dopar% {
    meanstat <- matrix(unlist(lapply(sumstat_window, function(x) rowMeans(x[,seq(i,ncol(x),by=4),drop=FALSE]))), nrow=nrow(sumstat_window[[1]]))
    quantile05 <- matrix(unlist(lapply(sumstat_window, function(x) apply(x[,seq(i,ncol(x),by=4),drop=FALSE], 1, function(z) {quantile(z,0.05)}))), nrow=nrow(sumstat_window[[1]]))
    quantile25 <- matrix(unlist(lapply(sumstat_window, function(x) apply(x[,seq(i,ncol(x),by=4),drop=FALSE], 1, function(z) {quantile(z,0.25)}))), nrow=nrow(sumstat_window[[1]]))
    quantile75 <- matrix(unlist(lapply(sumstat_window, function(x) apply(x[,seq(i,ncol(x),by=4),drop=FALSE], 1, function(z) {quantile(z,0.75)}))), nrow=nrow(sumstat_window[[1]]))
    quantile95 <- matrix(unlist(lapply(sumstat_window, function(x) apply(x[,seq(i,ncol(x),by=4),drop=FALSE], 1, function(z) {quantile(z,0.95)}))), nrow=nrow(sumstat_window[[1]]))
    list(mean=meanstat, q05=quantile05, q25=quantile25, q75=quantile75, q95=quantile95)
  }

  count <- count + 1
}

#Save a summary file with these statistics
saveRDS(list(stats=meanquantilestat, tag=unlist(tag)), file = paste(savetag, "_summary_window.rds",sep=""))
# tmp <- readRDS(file = paste(savetag, "_summary_window.rds",sep=""))
