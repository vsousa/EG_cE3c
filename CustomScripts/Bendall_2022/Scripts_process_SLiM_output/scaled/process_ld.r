# SCRIPT TO READ MS FILES TO COMPUTE LD
# PLOT THE LD PATTERNS
# LD is measured with statistic r^2

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

# COMBINATIONS OF PARAMETERS TESTED
# migration rate
migrate_v <-c("0.000", "0.00034", "0.0017", "0.0034")
# selection of mutation under divergent selection
#selmutben_v <- c("0.000", "0.0067", "0.01333" ,"0.02667", "0.05333", "0.0667", "0.1334")
selmutben_v <- c("0.000")
# recombination rate
recrate_v <- c("2.5e-7")
# initial frequency of allele a at locus under divergent selection
#freqmut_v <- c("0.001", "0.01", "0.1", "0.5")
freqmut_v <- c("0.001")
# dominance of allele a at locus under divergent selection
#dominance_v <- c("0.01","0.5")
dominance_v <- c("0.01")
# time of split of the two populations
time_v <- c("2000")
# TRUE if the name of files indicate the initial frequency
freqTag <- TRUE

# Get the number of parameter combinations
combnum <- length(migrate_v)*length(selmutben_v)*length(recrate_v)*length(dominance_v)*length(freqmut_v)*length(time_v)*2 

# tag to name of output files
savetag <- paste("neutral_ld_f", freqmut_v[1], sep="")

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
cond_folder <- "sumstat"
# Auxiliary files to save the file tag and legend text for each scenario
tag <- list()
legtext <- list()

# define folder with the results of the simulations
folder_v <- c("divsel3_500Kb_scaled_extra_neutral")

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
                    
                    # get the tag name of the output files
                    if(freqTag) {
                      tag[[scenario]]=paste("divsel_", chrm, "_h", dominanceben,"_", dominancedel,"_m", migrate ,"_sr", sexratio, "_s", selmutben, "_", selmutdel ,"_r", recrate, "_f", freqmut, "_", timeend ,sep="")  
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

# Define the folder for a given combination of parameters
folder <- rep(folder_v[1], times=combnum)

# check that all the files exist
eval <- numeric(length(tag))
{
  for(scenario in 1:length(tag)) {
    eval[scenario] <- file.exists(paste("./",folder[scenario], "/", cond_folder,"/r2_mean_sd", tag[scenario],sep=""))
  }
}

# 1. Read the data
# if files do not exist show their path and file name
if(sum(eval==0)>0) {
  for(i in which(eval==0)) {
    print(paste("File does not exist: ./",folder[i], "/", cond_folder,"/r2_mean_sd", tag[i],sep=""))
  }
  stop("Some files do not exist!")
} 

# read the data
{
  # read chromosome wide results
  meanr2 <- foreach(scenario=1:length(tag)) %dopar% {
    readRDS(file=paste("./",folder[scenario], "/", cond_folder,"/r2_mean_sd", tag[scenario],sep=""))
  }
}

# check that all the datasets have the correct number of simulations
if(cond_folder=="sumstat") {
  if(sum(unlist(lapply(meanr2, ncol))!=nsims)>0) {
    stop("Not all scenarios have the correct number of sims!")
  }  
}


# for each combination of parameters save the summary of the r2 across the 25 windows of 20Kb
summaryr2 <- lapply(meanr2, function(x) {
  # merge results for pop1 (x[[1]]) and pop2 (x[[2]])
  # compute mean and quantile for each
  c(mean(c(x[[1]],x[[2]]), na.rm=TRUE),
    quantile(c(x[[1]],x[[2]]), c(0.05,0.25,0.5,0.75,0.95), na.rm=TRUE))
})
# merge all results into matrix
summaryr2 <- do.call(rbind, summaryr2)

# save data.frame with summary of results
df <- data.frame(chrm=chrm, m=migrate, s=selection, f=freq, r=recrate,
           meanr2=summaryr2[,1],
           q05=summaryr2[,2],
           q25=summaryr2[,3],
           q50=summaryr2[,4],
           q75=summaryr2[,5],
           q95=summaryr2[,6])
write.table(df, file="./Linked_Sellocus/neutral_meanr2.txt", row.names=FALSE, quote=FALSE)