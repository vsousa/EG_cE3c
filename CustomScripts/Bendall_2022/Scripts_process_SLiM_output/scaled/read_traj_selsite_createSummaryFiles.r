# SCRIPT PROCESS ALLELE FREQUENCY TRAJECTORIES
# creates summary files with average and quantile of allele frequency trajectories
# and files with summary statistics (e.g. FST, probability of keeping allele a, etc.)

# clean memory - remove all variables
rm(list=ls())
# load required packages
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
registerDoParallel(cores=2)

# load files with functions to compute summary statistics
source("./Rscripts/computeStats.r")

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
# label for the runs
extrachrm <- rep(c("A","HD"), times=1)
# minimum number of runs required to compute some statistics
minruns <- 10 
# TRUE if the name of files indicate the initial frequency
freqTag <- TRUE 

# COMBINATIONS OF PARAMETERS
# migration rate
migrate_v <-c("0.000", "0.00034", "0.0017", "0.0034")
# selection of mutation under divergent selection
selmutben_v <- c("0.000", "0.0067", "0.01333" ,"0.02667", "0.05333", "0.0667", "0.1334")
# recombination rate
recrate_v <- c("2.5e-7")
# initial frequency of allele a at locus under divergent selection
freqmut_v <- c("0.001", "0.01", "0.1", "0.5")
# dominance of allele a at locus under divergent selection
dominance_v <- c("0.01","0.5")
# time of split of the two populations
time_v <- c("2000")
# folder can be different for different sims. Use this vector for those cases
folder_v <- c("divsel3_500Kb_scaled",
              "divsel3_500Kb_scaled_merged",
              "divsel3_500Kb_scaled_extra_neutral")

# get the number of combinations of parameters
combnum <- length(migrate_v)*length(selmutben_v)*length(recrate_v)*length(dominance_v)*length(freqmut_v)*length(time_v)*2 

# tag to name of output files
savetag <- paste("Traj_t",time_v[1],"_r",recrate_v[1], sep="")

# Auxiliary files to save the file tag and legend text for each scenario
tag <- list()
legtext <- list()

# DEBUG
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
# type <- 1

# Loop through all combinations of parameters
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
                      tag[[scenario]] <- paste("traj_divsel_", chrm, "_h", dominanceben,"_", dominancedel,"_m", migrate ,"_sr", sexratio, "_s", selmutben, "_", selmutdel ,"_r", recrate, "_f", freqmut, "_", timeend, sep="")    
                    } else {
                      tag[[scenario]] <- paste("traj_divsel_", chrm, "_h", dominanceben,"_", dominancedel,"_m", migrate ,"_sr", sexratio, "_s", selmutben, "_", selmutdel ,"_r", recrate, "_", timeend, sep="")    
                    }
                    
                    # Create the legend text for each
                    if(chrm == "X") {
                      extrachrm <- "H"
                    } else if(chrm == "A") {
                      extrachrm <- "D"
                    }
                    
                    # save legend text for each scenario
                    legtext[[scenario]] <- paste(extrachrm, "_m", migrate, "_r=", recrate, "_s=", selmutben, "_h=", dominanceben, "_f=", freqmut ,sep="")
                    #print(legtext[[scenario]])
                    
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

# get the index of combinations for X-chromosome (H - Hemizygous)
xindex <- which(chrm=="H")

# for the neutral case the dominance does not matter and with freq 0.001 and 0.01 
# the neutral sims were only done for for h=0.01 and not h=0.5
# get the correct file name for the neutral cases with initial freq less than 0.1
eval_neut <-(selection==0.000 & freq==0.001) | ((selection==0.000 & freq==0.01))
tag[eval_neut] <- sub(pattern="h0.5", replace="h0.01", x=tag[eval_neut])    

# Define the folder for a given combination of parameters
folder <- rep(folder_v[1], times=combnum)
folder[freq==0.1 | freq==0.5] <- folder_v[2]
folder[eval_neut] <- folder_v[3]

# check that all the files exist
eval <- numeric(length(tag))
{
  for(scenario in 1:length(tag)) {
    tmpfile <- paste("./",folder[scenario], "/", tag[scenario],".benmut",sep="")
    eval[scenario] <- file.exists(tmpfile)
  }
}

# 1. Read the data
# if files do not exist show their path and file name
if(sum(eval==0)>0) {
  for(i in which(eval==0)) {
    print(paste("File does not exist: ./",folder[i], "/", tag[i],".benmut",sep=""))
  }
  stop("Some files do not exist!")
} 
# read trajectory for all simulations of a given parameter combination
traj <- foreach(scenario=1:length(tag)) %dopar% {
  tmpfile <- paste("./",folder[scenario], "/", tag[scenario],".benmut",sep="")
  if(file.exists(tmpfile)) {
    res <- matrix(scan(file=tmpfile), ncol=6, byrow = T)  
  } else {
    stop(paste("File:", tmpfile, " not found", sep=""))
  }
  res
}

# Define last generation - this can be adjusted if we want to look at other times
ngenfinal <- time_v[1]

# 2. Remove duplicated trajectories from traj 
# (i.e. cases where the trajectory was exactly the same)
# due to appending of sims issues with the merged cases
traj_original <- traj
# traj <- lapply(traj, function(x) {
#   # to find the duplicated, we need to find cases where the run and gen are repeated
#   dup_tf <- duplicated(apply(x, 1, function(row) {paste(row[1:2],collapse = "_")}))
#   x[!dup_tf,]
# })
# same code using parallel package
traj <- foreach(scenario=1:length(tag)) %dopar% {
  # variable x refers to the trajectory of each parameter combination
  x <-  traj[[scenario]]
  # to find the duplicated, we need to find cases where the run and gen are repeated
  dup_tf <- duplicated(apply(x, 1, function(row) {paste(row[1:2],collapse = "_")}))
  # output the trajectories without the duplicated cases
  x[!dup_tf,]
}
rm(traj_original) 

# 3. keep only the last generation
# lastgen <- lapply(traj, function(x) {
#   # get the index corresponding to the final generation for each simulation
#   lastgen_i <- which(x[,2]==ngenfinal)
#   # output matrix just with the frequency at each generation for each simulation
#   x[lastgen_i,]
# })
# same code using parallel package
lastgen <- foreach(scenario=1:length(tag)) %dopar% {
  # variable x refers to the trajectory of each parameter combination
  x <-  traj[[scenario]]
  # get the index corresponding to the final generation for each simulation
  lastgen_i <- which(x[,2]==ngenfinal)
  # output matrix just with the frequency at each generation for each simulation
  x[lastgen_i,]
}
#str(lastgen) 

# check that all cases have the correct number of sims
if(sum(sapply(lastgen, nrow)!=nsims)>0) {
  stop("Number of values is different from number of sims!")
}
# get a matrix of booleans with TRUE/FALSE to evaluate the sims 
# where the allele beneficial in p1 was lost
lost_sim <- sapply(lastgen, function(x) {x[1:nsims,3]==0}, simplify=TRUE)
# str(lost_sim)

# check that all the other params are the same for A and X
if(sum(migrate[aindex]!=migrate[xindex])>0 |
   sum(recrate[aindex]!=recrate[xindex])>0 |
   sum(dominance[aindex]!=dominance[xindex])>0 |
   sum(selection[aindex]!=selection[xindex])>0 |
   sum(freq[aindex]!=freq[xindex])>0) {
    stop("Error: different order of params for A and X!")
}

# 4. Save mean trajectories and quantiles after removing duplicates and corresponding tag into a file
# loop through each scenario to compute the mean and quantiles of allele frequencies across all runs
summary_traj <- foreach(scenario=1:length(tag)) %dopar% {
  # temporary variable that saves the trajectory for a given scenario
  d <- traj[[scenario]]
  # mean trajectory and quantiles of allele a in p1
  # get the mean at each generation at pop1
  tmp_mean_p1 <- as.numeric(tapply(d[,3], INDEX=d[,2], FUN=mean, simplify = TRUE))
  # get the quantile at each generation at pop1
  tmp_quantile_p1 <- matrix(as.numeric(unlist(tapply(d[,3], INDEX=d[,2], FUN=function(x) {quantile(x, c(0.05, 0.25,0.5,0.75,0.95))}, simplify = TRUE))),ncol=5, byrow=T)
  # mean trajectory and quantiles of allele a in p2
  # get the mean at each generation at pop2
  tmp_mean_p2 <- as.numeric(tapply(d[,4], INDEX=d[,2], FUN=mean, simplify = TRUE))
  # get the quantile at each generation at pop2
  tmp_quantile_p2 <- matrix(as.numeric(unlist(tapply(d[,4], INDEX=d[,2], FUN=function(x) {quantile(x, c(0.05, 0.25,0.5,0.75,0.95))}, simplify = TRUE))),ncol=5, byrow=T)
  # merge into a summary matrix
  cbind(d[1:length(tmp_mean_p1),2], 
        tmp_mean_p1,tmp_quantile_p1,
        tmp_mean_p2,tmp_quantile_p2)
}

# CONDITIONAL ON KEEPING ALLELE
# loop through each scenario to compute the mean and quantiles of allele frequencies across runs that maintained the allele a in pop1
# define the time points at which we have the frequency data
gens <- sort(unique(traj[[1]][,2]))
summary_traj_cb <- foreach(scenario=1:length(tag)) %dopar% {
  # condition of maintaining the allele a
  cond <- rep(!lost_sim[,scenario], each=length(gens))
  if(sum(cond)>0) {
    # temporary variable that saves the trajectory for a given scenario
    # conditional on maintaining the allele a
    d <- traj[[scenario]][cond,]
    # mean trajectory and quantiles of allele a in p1
    # get the mean at each generation at pop1
    tmp_mean_p1 <- as.numeric(tapply(d[,3], INDEX=d[,2], FUN=mean, simplify = TRUE))
    # get the quantile at each generation at pop1
    tmp_quantile_p1 <- matrix(as.numeric(unlist(tapply(d[,3], INDEX=d[,2], FUN=function(x) {quantile(x, c(0.05, 0.25,0.5,0.75,0.95))}, simplify = TRUE))),ncol=5, byrow=T)
    # mean trajectory and quantiles of allele a in p2
    # get the mean at each generation at pop2
    tmp_mean_p2 <- as.numeric(tapply(d[,4], INDEX=d[,2], FUN=mean, simplify = TRUE))
    # get the quantile at each generation at pop2
    tmp_quantile_p2 <- matrix(as.numeric(unlist(tapply(d[,4], INDEX=d[,2], FUN=function(x) {quantile(x, c(0.05, 0.25,0.5,0.75,0.95))}, simplify = TRUE))),ncol=5, byrow=T)
  } else {
    tmp_mean_p1 <- NA
    tmp_quantile_p1 <- rep(NA, times=5)
    tmp_mean_p2 <- NA
    tmp_quantile_p2 <- rep(NA, times=5)
  }
  # merge into a summary matrix
  cbind(gens, 
        tmp_mean_p1,tmp_quantile_p1,
        tmp_mean_p2,tmp_quantile_p2)
}

saveRDS(list(summary_traj=summary_traj, summary_traj_cb=summary_traj_cb), file = paste(savetag, "_summary_traj.rds",sep=""))
# tmp <- readRDS(file = paste(savetag, "_summary.rds",sep=""))

#########################################################
# Compute and save statistics 
#########################################################

# 1. MEAN ALLELE FREQUENCY
# get the mean allele frequency for pop1 and pop2 at last generation
af <- sapply(lastgen, function(x) {colMeans(x[,3:4])})
# str(af)

# 2. PROBABILITY OF "FIXATION"- REACHING DETERMINISTIC EQUILIBRIUM FREQS IN BOTH POPS
#    Define a threshold for each migration rate
#    threshold based on mean allele freq from deterministic sims
#    for selection between 0.005<s<0.1
#    and migration rates 0.0000001000 0.0003162278 0.0010000000 0.0031622777
migdetfixP1 <- c(1.0000000, 0.9706445, 0.9198146, 0.8240068)
# migration threshold for the pop2 frequency
migdetfixP2 <- c(1.539097e-06, 8.593618e-02, 1.698009e-01, 3.019988e-01)
# migration threshold allele freq difference between pop1 and pop2
migdetfixdiff <- c(0.9999983, 0.8847083, 0.7500137, 0.5220079)
# round the thresholds based on deterministic sims
migdetfixP1 <- round(migdetfixP1, digits=3)
migdetfixP2 <- round(migdetfixP2, digits=3)
migdetfixdiff <- round(migdetfixdiff, digits=3)
# find the threshold for each scenario
pthreshold <- numeric(combnum)
for(i in 1:length(migrate_v)) {
  pthreshold[migrate==as.numeric(migrate_v[i])] <- migdetfixdiff[i]
}

# pfix = probability that frequency in last generation is larger than pthreshold
pfix <- unlist(sapply(1:combnum, function(i) {
                        x <- lastgen[[i]];
                        # get number of runs higher than threshold and divide by those
                        # where mutation was not lost (i.e. with a freq > 0)
                        sum((x[,3]-x[,4]) >= pthreshold[i])}))
# number of runs where the allele was maintained in pop1
nkeep <- unlist(sapply(1:combnum, function(i) {
  x <- lastgen[[i]];
  # get number of runs higher than threshold and divide by those
  # where mutation was not lost (i.e. with a freq > 0)
  sum(x[,3] > 0)}))

# DEAL WITH POTENTIAL ERRORS DUE TO SMALL NUMBER OF RUNS
# cases in which the proportion of runs here freq in last generation is less than 1/nsims is considered error
error_pfix <- 1/nsims
pfix[pfix<=error_pfix] <- 0

# 3. AVERAGE TIME TO REACH "FIXATION"
#   Find the average time for reavhing a frequency larger than the deterministic migration-selection equilibrium
mean_tfix <- foreach(i=1:combnum) %dopar% {
  # get the index of lines where it reached "fixation"
  # defining fixation as the difference between the allele freqs at P1 and P2
  fixed_i <- which((traj[[i]][,3]- traj[[i]][,4]) >= pthreshold[i]);
  if(length(fixed_i) > 0) {
    # get just the run and generation info for those that fixed
    fixed <- traj[[i]][fixed_i,1:2, drop=FALSE]
    # get index of runs that reached fixation
    fruns_i <- unique(fixed[,1])
    if(length(fruns_i)>minruns) {
      # go through each run and get the minimum generation for each run
      mtime_runs <- sapply(fruns_i, function(j) {min(fixed[fixed[,1]==j,2])})
    } else {
      mtime_runs <- NA;
    }   
  } else {
    mtime_runs <- NA;
  }
  mean(mtime_runs)
}
mean_tfix <- unlist(mean_tfix)


# 4. PROBABILITY OF KEEPING THE ALLELE
#    Probability of keeping beneficial allele in p1
ploss <- unlist(lapply(lastgen, function(x) sum(x[,3]==0)))/nsims
pkeep <- 1-ploss


# 5. FST BASED ON POPULATION ALLELE FREQUENCIES 
#    FST computed based on Lasne et al. 2017 equation based on population frequencies (rather than sample frequencies)
# DEAL WITH POTENTIAL ERRORS DUE TO SMALL VALUES
# minimum FST defined as 0.0001
minfst <- 1e-4
# compute fst for each sim and then make the mean of FST across sims
fst_sim <- sapply(lastgen, function(x) {
    vapply(1:nsims, function(j) {
      getFst_freq_lasne(x[j,3], x[j,4])}, FUN.VALUE = numeric(1))
    }, simplify = TRUE)
#str(fst_sim)
# get the mean fst per simulation
fst_eq_sim <- colMeans(fst_sim)
fst_eq_sim[fst_eq_sim<minfst] <- minfst
# get the mean fst conditional on keeping the a allele
fst_sim[lost_sim] <- NA
fst_eq_sim_cb <- colMeans(fst_sim, na.rm=TRUE)
# replace cases with very low FST by minfst
fst_eq_sim_cb[fst_eq_sim_cb<minfst] <- minfst
# replace the cases with NA by minfst
fst_eq_sim_cb[is.na(fst_eq_sim_cb)] <- minfst
# replace cases with less then this number of runs by missing data
notlost <- colSums(!lost_sim)
fst_eq_sim_cb[notlost <= minruns] <- NA


# 6. MEAN ALLELE FREQUENCIES CONDITIONAL ON KEEPING ALLELE
#    This is equivalent to the area under the curve.
#    the higher the area under the curve the quicker it reaches fixation.
#    average allele freq across all generations for each sims
p1a <- sapply(traj, function(x) {
  tapply(x[,3], INDEX=x[,1], FUN=mean, simplify = TRUE)[1:nsims]
}, simplify=TRUE)
#str(p1a)
# allele freq of allele a in pop2 for all sims
p2a <- sapply(traj, function(x) {
  tapply(x[,4], INDEX=x[,1], FUN=mean, simplify = TRUE)[1:nsims]
}, simplify=TRUE)
#str(p2a)

# replace lost sims by missing data
p1a_cb <- p1a
p1a_cb[lost_sim] <- NA
p2a_cb <- p2a
p2a_cb[lost_sim] <- NA

# compute average sum of allele frequencies across sims
p1a_mean <- colMeans(p1a, na.rm=TRUE)
p2a_mean <- colMeans(p2a, na.rm=TRUE)

# compute the average allele frequencies only for not lost cases
p1a_cb_mean <- colMeans(p1a_cb, na.rm=TRUE)
p2a_cb_mean <- colMeans(p2a_cb, na.rm=TRUE)

# SAVE THE RESULTS INTO TABLES
# Save a summary of the key summary results
tmp <- data.frame(m=migrate,rec=recrate, h=dominance,chrm=chrm,s=selection,f=freq,
                  pfix=pfix,mtfix=mean_tfix,pkeep=pkeep,fstall=fst_eq_sim,fstcb=fst_eq_sim_cb,
                  meangen_p1all=p1a_mean,meangen_p2all=p2a_mean,meangen_p1cb=p1a_cb_mean,meangen_2cb=p2a_cb_mean)
write.table(tmp, file=paste("./Linked_Sellocus/Summary_", savetag, ".txt",sep=""), quote=FALSE, row.names = FALSE)
