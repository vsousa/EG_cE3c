# SCRIPT TO PROCESS SUMMARY FILES OF MS OUTPUT
# creates summary files with average and quantile of summary statistics 
# (e.g. FST, PI, SFS) that can then be used for functions to plot results
# from SLiM simulations performed according to the Neodiprion demograhic history IM model

# clean memory
rm(list=ls())

# load required packages
library(RColorBrewer)
library(doParallel)
library(ggplot2)
registerDoParallel(cores=6)

# SETTINGS
# number of simulations
nsims <- 1000
# type of mutations sampled
muttype <- "all"
# length of simulated chromosomes
seqlength <- 500000
# position of locus with mutation under divergent selection
positionben <- 250000
# sample size from each population (in number of genome copies, e.g. 20 indicates 10 diploid inds)
ss <- c(8,12) # sample size for pop1 and pop2
# number of populations
npop <- 2
# sliding window analysis
window.size <- 20000 # window size
slide.size  <- 20000 # slide size

# COMBINATIONS OF PARAMETERS TESTED
# migration rate
migrate_v <-c("3.65e-4")
# selective coefficient (as a function of the Ne of N. pinetum))
selmutben_v <- c("0.000","0.0002","0.0004","0.0009","0.0030","0.0042","0.0058","0.0081","0.0113","0.0217","0.0302","0.0420","0.0583","0.0810","0.113","0.156","0.217","0.302","0.0305", "0.0609", "0.1218", "0.2437", "0.3046", "0.6092") # selective coefficient (as a function of the Ne of N. pinetum))
#selmutben_v <- c("0.000","0.0305", "0.0609", "0.1218", "0.2437", "0.3046", "0.6092","0.0002","0.0004","0.0009","0.0030","0.0042","0.0058","0.0081","0.0113","0.0217","0.0302","0.0420","0.0583","0.0810","0.113","0.156","0.217","0.302") # selective coefficient (as a function of the Ne of N. pinetum))
# recombination rate
recrate_v <- c("1.05e-6")
# initial frequency of allele a at locus under divergent selection
freqmut_v <- c("0.1")
# dominance of allele a at locus under divergent selection
dominance_v <- c("0.01","0.5")
# time of split of the two populations
time_v <- c(1549)
# TRUE if the name of files indicate the initial frequency
freqTag <- TRUE

# define folder with the results of the simulations
folder_v <- c("sawfly_divsel3_500Kb_scaled_rec2_3","sawfly_divsel3_500Kb_scaled_rec2_3")

# Get the number of parameter combinations
combnum <- length(migrate_v)*length(selmutben_v)*length(recrate_v)*length(dominance_v)*length(freqmut_v)*length(time_v)*2 

# tag to name of output files
savetag <- paste("sawfly_divsel3_t1549_scaled_r2_3_", recrate_v[1],"_w",window.size,"_f",freqmut_v[1],sep="") # hrec stands for high recombination  
if(recrate_v[1]=="3.5e-7") {
  savetag <- paste("sawfly_divsel3_t1549_scaled_r2_3_w",window.size,sep="") 
} 


# for debug
# chrm <- "A"
# migrate <- migrate_v[1]
# selmutben <- selmutben_v[1]
# selmutdel <- c("0.000")
# recrate <- recrate_v[1]
# sexratio <- c("0.5")
# dominanceben <- c("0.01")
# dominancedel <- c("0.5")
# freqmut <- freqmut_v[1]
# cond_folder <- "sumstat"

# Analyse the summary files for:
# - "sumstat" all sims, including those where the allele a was lost
# - "sumstat_cb" only sims conditional on keeping the beneficial allele a in pop1
# initialize a list that will save the results of window based analysis
meanquantilestat <- list()
count <- 1
for(cond_folder in c("sumstat","sumstat_cb")) {
  # Auxiliary files to save the file tag and legend text for each scenario
  tag <- list()
  legtext <- list()
  
  # Go through all combinations of parameters
  scenario <- 1
  for(migrate in migrate_v) {
    for(selmutben in selmutben_v) {
      for(selmutdel in c("0.000")) {
        for(recrate in recrate_v) {
          for(sexratio in c("0.3")) {
            for(dominanceben in dominance_v) {
              for(dominancedel in c("0.5")) {
                for(freqmut in freqmut_v) {
                  for(timeend in time_v) {
                    for(chrm in c("A","X")) {
                      # Get the file tag for the input files
                      if(freqTag) {
                        tag[[scenario]]=paste("sawfly_divsel3_", chrm, "_h", dominanceben,"_", dominancedel,"_m", migrate ,"_sr", sexratio, "_s", selmutben, "_", selmutdel ,"_r", recrate, "_f", freqmut, "_", timeend, sep="")  
                      } else {
                        tag[[scenario]]=paste("sawfly_divsel3_", chrm, "_h", dominanceben,"_", dominancedel,"_m", migrate ,"_sr", sexratio, "_s", selmutben, "_", selmutdel ,"_r", recrate, sep="")  
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
  folder[aindex] <- folder_v[2]
  
  # check that all the files exist
  eval <- numeric(length(tag))
  for(scenario in 1:length(tag)) {
      eval[scenario] <- file.exists(paste("./",folder[scenario], "/", cond_folder,"/genome_sumstat_", tag[scenario],sep=""))
  }
  
  # 1. Read the data
  # if files do not exist show their path and file name
  if(sum(eval==0)>0) {
    for(i in which(eval==0)) {
      print(paste("File does not exist: ./",folder[i], "/", cond_folder,"/genome_sumstat_", tag[i],sep=""))
    }
    stop("Some files do not exist!")
  } 
  
  # Read the data
  {
    # read genome wide results
    sumstat_genome <- foreach(scenario=1:length(tag)) %dopar% {
      read.table(file=paste("./",folder[scenario], "/", cond_folder,"/genome_sumstat_", tag[scenario],sep=""), stringsAsFactors = FALSE, header=T)  
    }
    
    sumstat_window <- foreach(scenario=1:length(tag)) %dopar% {
      read.table(file=paste("./",folder[scenario], "/", cond_folder,"/window_sumstat_", tag[scenario],sep=""), stringsAsFactors = FALSE, header=T) 
    }
    
    if(cond_folder=="sumstat") {
      sfs2d <- foreach(scenario=1:length(tag)) %dopar% {
        read.table(file=paste("./",folder[scenario],"/", cond_folder,"/sfs2d_", tag[scenario],sep=""), stringsAsFactors = FALSE, header=F)
      }    
    }
  }

  # check that all the datasets have the correct number of simulations
  if(cond_folder=="sumstat") {
    if(sum(unlist(lapply(sumstat_genome, nrow))!=nsims)>0) {
      stop("Not all scenarios have the correct number of sims!")
    } 
    # NEUTRAL SFS for Diploid and Haplodiploid
    # Save the mean SFS for the neutral case, summing over all the 1000 simulations
    sfs2dsum <- lapply(sfs2d[which(selection==0)], function(x) {matrix(colSums(x), nrow=ss+1)}) 
    # get the sum of the SFS for the diploid case
    sfs2dsum_d <- Reduce(f="+", sfs2dsum[chrm[which(selection==0)]=="D"])
    # get the sum of the SFS for the haplodiploid case
    sfs2dsum_hd <- Reduce(f="+", sfs2dsum[chrm[which(selection==0)]=="HD"])
    # Save a file with the list of SFS
    saveRDS(list(neutral_sfs_d=sfs2dsum_d, neutral_sfs_hd=sfs2dsum_hd), file = paste("./SummaryFiles/",savetag, "_neutral_sfs.rds",sep=""))
  }  

  # COMPUTE THE MEAN AND QUANTILES OF SEVERAL STATISTICS
  mean_stats <- foreach(i=1:ncol(sumstat_genome[[1]])) %dopar% {
    as.numeric(unlist(lapply(sumstat_genome, function(x) mean(x[,i]))))
  }
  # mean_stats <- matrix(unlist(mean_stats), ncol=4)
  quantile_stats <- foreach(i=1:ncol(sumstat_genome[[1]])) %dopar% {
    matrix(as.numeric(unlist(lapply(sumstat_genome, function(x) quantile(x[,i], c(0.05, 0.25, 0.5, 0.75, 0.95))))), ncol=5, byrow=T)
  }
  # save a summary of the results into a table
  savetext <- c("pi1","pi2","dxy","fst")
  savetextaux <- c("mean","quantile0.05","quantile0.25","quantile0.50","quantile0.75","quantile0.95")
  tmp <- data.frame(m = migrate[aindex], 
                    r = recrate[aindex], 
                    d = dominance[aindex],
                    s = selection[aindex],
                    f = freq[aindex])
  # save the results into one file for each stat
  for(i in 1:length(savetext)) {
    aux <- cbind(tmp, mean_stats[[i]][aindex])
    aux <- cbind(aux, quantile_stats[[i]][aindex,])
    aux <- cbind(aux, mean_stats[[i]][xindex])
    aux <- cbind(aux, quantile_stats[[i]][xindex,])
    colnames(aux) <- c(names(tmp),paste("D", savetext[i], savetextaux, sep="_"),paste("HD", savetext[i], savetextaux, sep="_"))
    write.table(aux, file=paste("./SummaryFiles/", savetext[i],"_", savetag, cond_folder,".txt",sep=""), quote=FALSE, row.names = FALSE)
  }
  
  # SAVE SUMMARIES ALONG THE CHROMOSOME
  # summary statistics recorded for each window
  list_sumstat <- c("pi1","pi2","dxy","fst")
  # check that all the datasets have the correct number of simulations
  if(cond_folder=="sumstat") {
    if(sum((unlist(lapply(sumstat_window, ncol))/length(list_sumstat))!=nsims)>0) {
      stop("Not all scenarios have the correct number of sims!")
    }  
  }
  
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
  
# Save a summary file with these statistics
saveRDS(list(stats=meanquantilestat, tag=unlist(tag)), file = paste("./SummaryFiles/", savetag, "_summary_window.rds",sep=""))
# tmp <- readRDS(file = paste(savetag, "_summary_window.rds",sep=""))
