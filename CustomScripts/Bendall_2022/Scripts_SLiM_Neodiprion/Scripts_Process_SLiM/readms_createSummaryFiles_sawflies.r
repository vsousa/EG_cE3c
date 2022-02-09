# PROCESS ms files from SLiM simulations
# according to demographic model inferred for sawflied
# Neodiprion lecontei and N. lecontei.
# read ms files and output summary files 

# clean memory: remove all variables from memory
rm(list=ls())
# load used packages
library(stringr)
library(RColorBrewer)
library(doParallel)
registerDoParallel(cores=5)

# load files with functions to compute summary statistics
source("computeStats_SFS.r")
# source("ldfunctions.r")

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

# define folder with the results of the simulations
# For diploid (A) use folder sawfly_divsel3_500Kb_scaled_rec2_3 and chrm_v <- c("A")
# For haplodiploid (X) use folder sawfly_divsel3_500Kb_scaled and chrm_v <- c("X")
folder <- rep("sawfly_divsel3_500Kb_scaled_rec2_3",times=1000)
chrm_v <- c("A","X")
# folder <- rep("sawfly_divsel3_500Kb_scaled",times=1000)
# chrm_v <- c("X")

scenario <- 1

# COMBINATIONS OF PARAMETERS TESTED
# migration rate
migrate_v <-c("3.65e-4")
# selective coefficient (as a function of the Ne of N. pinetum))
selmutben_v <- c("0.000","0.0305", "0.0609", "0.1218", "0.2437", "0.3046", "0.6092","0.0002","0.0004","0.0009","0.0030","0.0042","0.0058","0.0081","0.0113","0.0217","0.0302","0.0420","0.0583","0.0810","0.113","0.156","0.217","0.302") 
# recombination rate
recrate_v <- c("1.05e-6")
# initial frequency of allele a at locus under divergent selection
freqmut_v <- c("0.1")
# dominance of allele a at locus under divergent selection
dominanceben_v <- c("0.01","0.5")
# time of split of the two populations
timeend_v <- c(1549)
# TRUE if the name of files indicate the initial frequency
freqTag <- TRUE
# sex-ratio defined as number of males/(number males + number of females)
sexratio_v <- 0.3

# total sample size
tss <- sum(ss)

# determine the window size around the selected site we want to look at
# assuming that we will want to have the same number of sites as in the observed data
# adjust theta to have the same expected number of SNPs as in observed data
# mean number of SNPs depending on window size
# 100Kb - 19.02 
# 50kb  -  9.96
# 25kb  -  5.93
# 10kb  -  3.92
# vector with the number of snps and corresponding window size
obs_window_snp <- c(19.02, 9.96, 5.93, 3.92) 
obs_window_size <- c(100000,50000,25000,10000)     

# the scaling was done such that 500Kb windows in SLIM would have the same number of SNPs as 
# 50000 callable sites, concatenating all the RAD loci together.
# The total number of callable sites in the observed RAD dataset was 252277 sites
# and there were 9994 SNPs out of the 252277 callable sites in the observed dataset
# read the observed SFS with the total number of SNPs and callable sites
# read observed sfs for intergenic regions
obs <- read.table("../SFS/excluded_1kb_Neo_L_P_2DSFS.joint2DMSFS", skip=1, header=TRUE, row.names = 1)
# remove the monomorphic sites, 
# setting the entry with fixed ancestral in both pops to zero, i.e. entry (0,0) 
# that in matrices of R corresponding to entry [1,1]
obs_nomon <- obs
obs_nomon[1,1] <- 0
# callable sites
ncallable <- sum(obs)
nsnps <- sum(obs_nomon)
# read the results from neutral sims to get the number of SNPs out of 500Kb in sims according to Neodiprion sawfly sims.
tag <- paste("sawfly_divsel3_X_h0.01_0.5_m", migrate_v[1] ,"_sr", sexratio_v[1], "_s0.000_0.000_r", recrate_v[1], "_f", freqmut_v[1], "_", timeend_v[1] ,sep="")  
if(file.exists(paste("./sawfly_divsel3_500Kb_scaled_rec2_3/ms_", tag,"_", muttype,".txt",sep=""))) {
  # read the ms files to get the average number of SNPs
  ms_filename <- paste("./sawfly_divsel3_500Kb_scaled_rec2_3/ms_", tag,"_", muttype,".txt",sep="")
  ms <- readLines(ms_filename)
  positions_file <- lapply(lapply(lapply(str_split(grep('positions:', ms, value=T), pattern = "positions:"), function(x) x[2]), function(x) as.numeric(str_split(x, " ", simplify=T))[-1]), function(x) round(x*seqlength))
  
  # Get the ms as a matrix and positions as a vector, discarding duplicated positions
  l <- foreach(sim=1:nsims) %dopar% {
    # read only the sequence data
    start <- (3*sim)+((sim-1)*tss)+1
    end <- start+tss-1
    ms_original <- ms[start:end]
    # 2nd line, positions
    positions_tmp <- positions_file[[sim]]
    
    # Delete monomorphic sites in the resampled ms
    ms_matrix_tmp <- strsplit(ms_original, split="")
    ms_matrix_tmp <- t(matrix(as.numeric(unlist(ms_matrix_tmp)), ncol=tss))
    # Get the duplicated entries
    # dups <- which(duplicated(positions_tmp))
    # dups_sites <- as.numeric(matrix(c(dups-1,dups), nrow=2, byrow = TRUE))
    
    # remove the sites that are monomorphic
    freqSites_pop <- colSums(ms_matrix_tmp)
    mono <- which(freqSites_pop==0 | freqSites_pop== tss) 
    
    # sites to remove
    #rm_sites <- union(mono, dups)
    rm_sites <- mono
    
    # remove the sites with two mutations, keeping only one of them
    if(length(rm_sites)>0) {
      return(list(p=positions_tmp[-rm_sites], m=ms_matrix_tmp[,-rm_sites]))  
    } else { 
      return(list(p=positions_tmp, m=ms_matrix_tmp))  
    }
    #list(p=positions_tmp, m=ms_matrix_tmp)
  }
  
  # get the positions and ms_matrix as a list
  positions_original <- lapply(l, function(x) x$p)
} else {
  stop("could not open an example file to determine the average number of SNPs under neutrality.")
}
# get the number of SNPs per sim
nsnps_neutral <- mean(unlist(lapply(positions_original, length)))

# get the window size we will investigate, with the same number of SNPs as in real data using a simple proportion
# if 500,000bp (seqlength) in SLIM sims correspond to nsnsp_neutral
# then x bp in SLIM correspond to obs_window_snp
window.size <- round((seqlength*obs_window_snp)/nsnps_neutral)

# we will need to have windows where the selected site are at the middle
# adding the above window.add to below and upper the selected site
window.add <- round(window.size/2)


# # for debug
# chrm <- "A"
# migrate <- migrate_v[1]
# selmutben <- selmutben_v[19]
# selmutdel <- c("0.000")
# recrate <- recrate_v[1]
# sexratio <- sexratio_v[1]
# dominanceben <- dominanceben_v[1]
# dominancedel <- c("0.5")
# timeend <- timeend_v[1]
# freqTag <- TRUE
# freqmut <- freqmut_v[1]

# loop through all combination of parameters
sumstat <- list()
sumstat_cb <- list()
scenario <- 1
for(migrate in migrate_v) {
  for(selmutben in selmutben_v) {
    for(selmutdel in c("0.000")) {
      for(recrate in recrate_v) {
        for(sexratio in sexratio_v[1]) {
          for(dominanceben in dominanceben_v) {
            for(dominancedel in c("0.5")) {
              for(freqmut in freqmut_v) {
                for(timeend in timeend_v) {
                  for(chrm in chrm_v) {
                    if(freqTag) {
                      tag=paste("sawfly_divsel3_", chrm, "_h", dominanceben,"_", dominancedel,"_m", migrate ,"_sr", sexratio, "_s", selmutben, "_", selmutdel ,"_r", recrate, "_f", freqmut, "_", timeend ,sep="")  
                    } else {
                      tag=paste("sawfly_divsel3_", chrm, "_h", dominanceben,"_", dominancedel,"_m", migrate ,"_sr", sexratio, "_s", selmutben, "_", selmutdel ,"_r", recrate, sep="")  
                    }
                    
                    print(tag)
                    # Create folder to save the sumstat if it is the first or if folder does not exist
                    if(scenario==1 | !dir.exists(paste("./",folder[scenario], "/sumstat",sep=""))) {
                      dir.create(paste("./",folder[scenario], "/sumstat",sep=""))
                      dir.create(paste("./",folder[scenario], "/sumstat_cb",sep=""))
                    }
                    
                    # check if file exist
                    if( file.exists(paste("./", folder[scenario],"/ms_", tag,"_", muttype,".txt",sep="")) ) {
                      # read the ms files and perform overall analysis
                      # testms <- readMS(paste("./", folder[scenario],"/ms_", tag,"_", muttype,".txt",sep=""), big.data = T)
                      ms_filename <- paste("./", folder[scenario],"/ms_", tag,"_", muttype,".txt",sep="")
                      ms <- readLines(ms_filename)
                      positions_file <- lapply(lapply(lapply(str_split(grep('positions:', ms, value=T), pattern = "positions:"), function(x) x[2]), function(x) as.numeric(str_split(x, " ", simplify=T))[-1]), function(x) round(x*seqlength))
                      
                      # Get the ms as a matrix and positions as a vector, discarding duplicated positions
                      l <- foreach(sim=1:nsims) %dopar% {
                        # read only the sequence data
                        start <- (3*sim)+((sim-1)*tss)+1
                        end <- start+tss-1
                        ms_original <- ms[start:end]
                        # 2nd line, positions
                        positions_tmp <- positions_file[[sim]]
                        
                        # Delete monomorphic sites in the resampled ms
                        ms_matrix_tmp <- strsplit(ms_original, split="")
                        ms_matrix_tmp <- t(matrix(as.numeric(unlist(ms_matrix_tmp)), ncol=tss))
                        # Get the duplicated entries
                        # dups <- which(duplicated(positions_tmp))
                        # dups_sites <- as.numeric(matrix(c(dups-1,dups), nrow=2, byrow = TRUE))
                        
                        # remove the sites that are monomorphic
                        freqSites_pop <- colSums(ms_matrix_tmp)
                        mono <- which(freqSites_pop==0 | freqSites_pop== tss) 
                        
                        # sites to remove
                        #rm_sites <- union(mono, dups)
                        rm_sites <- mono
                        
                        # remove the sites with two mutations, keeping only one of them
                        if(length(rm_sites)>0) {
                          return(list(p=positions_tmp[-rm_sites], m=ms_matrix_tmp[,-rm_sites]))  
                        } else {
                          return(list(p=positions_tmp, m=ms_matrix_tmp))  
                        }
                        #list(p=positions_tmp, m=ms_matrix_tmp)
                      }
                      
                      # get the positions and ms_matrix as a list
                      positions_original_all <- lapply(l, function(x) x$p)
                      ms_matrix_all <- lapply(l, function(x) x$m)
                      rm(l)
                      # aux for index of pop
                      aux_pop <- list(c(1:ss[1]),c((ss[1]+1):tss)) 
                      
                      # positions of the list where the mutation was not lost
                      benmut <- which(unlist(lapply(positions_original_all, function(x) sum(x==positionben)))==2)
                      
                      # get the sites just at the selected window size
                      # initialize variables to save the mean_sumstat (all sims), mean_sumstat_cb (conditional on keeping beneficial)
                      # get the sites just at the selected window size
                      # initialize variables to save the mean_sumstat (all sims), mean_sumstat_cb (conditional on keeping beneficial), 
                      mean_sumstat <- matrix(NA, ncol=4, nrow=length(obs_window_size)) # each column corresponds to pi1, pi2, dxy, fst
                      colnames(mean_sumstat) <- c("pi_pop1","pi_pop2","dxy","fst")
                      mean_sumstat_cb <- matrix(NA, ncol=4, nrow=length(obs_window_size)) # each column corresponds to pi1, pi2, dxy, fst
                      colnames(mean_sumstat_cb) <- c("pi_pop1","pi_pop2","dxy","fst")
                      # quantiles fst
                      quantiles_fst <- matrix(NA, ncol=9, nrow=length(obs_window_size)) # each column corresponds to pi1, pi2, dxy, fst
                      colnames(quantiles_fst) <- c("fst0.025","fst0.05","fst0.10","fst0.25","fst0.5","fst0.75","fst0.90","fst0.95","fst0.975")
                      quantiles_fst_cb <- matrix(NA, ncol=9, nrow=length(obs_window_size)) # each column corresponds to pi1, pi2, dxy, fst
                      colnames(quantiles_fst_cb) <- c("fst0.025","fst0.05","fst0.10","fst0.25","fst0.5","fst0.75","fst0.90","fst0.95","fst0.975")
                       
                      
                      # variable saving the combination of parameters
                      paramcomb <- matrix(c(rep(c(chrm, dominanceben, migrate ,sexratio, selmutben, recrate, freqmut, timeend), each=length(obs_window_size)),obs_window_size),nrow=length(obs_window_size))
                      paramcomb <- data.frame(chr=paramcomb[,1],h=paramcomb[,2],m=paramcomb[,3],sr=paramcomb[,4],s=paramcomb[,5],r=paramcomb[,6],
                                              f=paramcomb[,7],t=paramcomb[,8], window=paramcomb[,9])  
                      
                      # FOR LOOP through window sizes here
                      for(w in 1:length(obs_window_size)) {
                        # auxiliary variable defining the windows that depend on the window.add
                        seq_windows <- c(positionben-window.add[w],positionben+window.add[w])
                        
                        # get the index of SNPs to keep
                        snps_to_keep_i <- lapply(positions_original_all, function(x) which(x >= seq_windows[1] & x <= seq_windows[2]))
                        # get a subset of ms_matrix_all and positions_original_all just with sites in the target window
                        positions_original <- foreach(sim=1:nsims) %dopar% {
                          positions_original_all[[sim]][snps_to_keep_i[[sim]]]
                        }
                        ms_matrix <- foreach(sim=1:nsims) %dopar% {
                          ms_matrix_all[[sim]][,snps_to_keep_i[[sim]],drop=FALSE]
                        }
                        # check the the mean number of snps in positions is the same as in ms_matrix
                        if( mean(unlist(lapply(positions_original, length))) != mean(unlist(lapply(ms_matrix, ncol))) ) {
                          stop("Error in the selection of positions and ms_matrix!")
                        }
                        
                        # check that the mean number of SNPs is within 0.9 to 1.1 of the observed in real data under neutrality
                        if(selmutben=="0.000") {
                          mean_snps_window_sim <- mean(unlist(lapply(snps_to_keep_i, length)))
                          if(mean_snps_window_sim<0.75*obs_window_snp[w] | mean_snps_window_sim>1.25*obs_window_snp[w]) {
                            stop("The number of SNPs in the window is not the same as in sims")
                          }
                        }
                        
                        ############################################################
                        # Compute average pi, dxy and pairwise FST for each simulation
                        mean_pi_dxy_fst <- foreach(sim=1:nsims) %dopar% {
                          getMeanFst_Pi_Dxy(t(ms_matrix[[sim]]), ss)
                        }
                        mean_sumstat_sims <- matrix(unlist(mean_pi_dxy_fst), ncol=4, byrow = TRUE)
                        
                        # get the mean across sims
                        mean_sumstat[w,] <- colMeans(mean_sumstat_sims)
                        mean_sumstat_cb[w,] <- colMeans(mean_sumstat_sims[benmut,,drop=FALSE])
                        quantiles_fst[w,] <- quantile(mean_sumstat_sims[,4], c(0.025,0.05,0.10,0.25,0.5,0.75,0.9,0.95,0.975))
                        quantiles_fst_cb[w,] <- quantile(mean_sumstat_sims[benmut,4], c(0.025,0.05,0.10,0.25,0.5,0.75,0.9,0.95,0.975))
                      }
                      sumstat[[scenario]] <- cbind(paramcomb,mean_sumstat,quantiles_fst)
                      sumstat_cb[[scenario]] <- cbind(paramcomb,mean_sumstat_cb,quantiles_fst_cb)
                    } else {
                      stop("File not found!")
                    } 
                    
                    # go to next simulation!
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

length(sumstat)

# convert into a matrix
sumstat_matrix <- do.call(rbind, sumstat)
sumstat_cb_matrix <- do.call(rbind, sumstat_cb)

# merge the stats and r2 and save into file
write.table(sumstat_matrix, file=paste("./",folder[1], "/sumstat/sumstat_X_A_fst_targetwindow_f0.1.txt",sep=""), 
            row.names = F, col.names = TRUE, quote = FALSE)
# conditional on keeping benmut
write.table(sumstat_cb_matrix, file=paste("./",folder[1], "/sumstat_cb/sumstat_cb_X_A_fst_targetwindow_f0.1.txt",sep=""), 
            row.names = F, col.names = TRUE, quote = FALSE)

# remove variables
# rm(mean_sumstat,mean_pi_dxy_fst,mean_r2_pop,mean_sumstat_cb,mean_r2_pop_cb)
