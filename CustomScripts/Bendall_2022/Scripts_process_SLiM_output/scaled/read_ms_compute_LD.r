# PROCESS ms files from SLiM simulations
# read ms files and output summary files 
# computing the mean r2 (LD measure) 
# between all pairs of SNPs within a given window

# clean memory: remove all variables from memory
rm(list=ls())
# load used packages
library(stringr)
library(RColorBrewer)
library(doParallel)
registerDoParallel(cores=5)

# load files with functions to compute summary statistics
source("./Rscripts/computeStats.r")
source("./Rscripts/ldfunctions.r")

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
ss <- 20
# number of populations
npop <- 2
# sliding window analysis
window.size <- 20000 # window size
slide.size  <- 20000 # slide size

# define folder with the results of the simulations
folder <- rep("divsel3_500Kb_scaled_extra_neutral", times=10000)
scenario <- 1

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
#dominanceben_v <- c("0.01","0.5")
dominanceben_v <- c("0.01")
# time of split of the two populations
time_v <- c("2000")
# TRUE if the name of files indicate the initial frequency
freqTag <- TRUE

# loop through all combination of parameters
for(migrate in migrate_v) {
  for(selmutben in selmutben_v) {
    for(selmutdel in c("0.000")) {
      for(recrate in recrate_v) {
        for(sexratio in c("0.5")) {
          for(dominanceben in dominanceben_v) {
            for(dominancedel in c("0.5")) {
              for(freqmut in freqmut_v) {
                for(timeend in time_v) {
                  for(chrm in c("A","X")) {
                    # get the tag name of the output files
                    if(freqTag) {
                      tag=paste("divsel_", chrm, "_h", dominanceben,"_", dominancedel,"_m", migrate ,"_sr", sexratio, "_s", selmutben, "_", selmutdel ,"_r", recrate, "_f", freqmut, "_", timeend ,sep="")  
                    } else {
                      tag=paste("divsel_", chrm, "_h", dominanceben,"_", dominancedel,"_m", migrate ,"_sr", sexratio, "_s", selmutben, "_", selmutdel ,"_r", recrate, sep="")  
                    }
                    print(tag)
                    
                    # Create folder to save the summary statistics
                    # folder sumstat_cb is for results conditional on retaining allele a (cb stands for conditional beneficial)
                    if(scenario==1 | !dir.exists(paste("./",folder[scenario], "/sumstat_cb",sep=""))) {
                      dir.create(paste("./",folder[scenario], "/sumstat_cb",sep=""))
                    }
                    # folder sumstat is for results of all simulations
                    if(scenario==1 | !dir.exists(paste("./",folder[scenario], "/sumstat",sep=""))) {
                      dir.create(paste("./",folder[scenario], "/sumstat",sep=""))
                    }    
                    
                    # check if file exist
                    if( file.exists(paste("./", folder[scenario],"/ms_", tag,"_", muttype,".txt",sep="")) ) {
                      # read the ms files and perform overall analysis
                      ms_filename <- paste("./", folder[scenario],"/ms_", tag,"_", muttype,".txt",sep="")
                      ms <- readLines(ms_filename)
                      positions_file <- lapply(lapply(lapply(str_split(grep('positions:', ms, value=T), pattern = "positions:"), function(x) x[2]), function(x) as.numeric(str_split(x, " ", simplify=T))[-1]), function(x) round(x*seqlength))
                      
                      # auxiliary variable defining the windows
                      seq_windows <- seq(0,seqlength,by=window.size)
                      # Get the ms as a matrix and positions as a vector, discarding duplicated positions
                      # note: duplicated positions can occur due to neutral mutations in the same position
                      l <- foreach(sim=1:nsims) %dopar% {
                        # read only the sequence data for each sim
                        # for each sim we can get the line in the ms file corresponding to the start and end of the simulated haplotypes
                        start <- (3*sim)+((sim-1)*(ss*npop))+1
                        end <- start+(ss*npop)-1
                        ms_original <- ms[start:end]
                        # 2nd line, positions
                        positions_tmp <- positions_file[[sim]]
                        
                        # Delete monomorphic sites in the resampled ms
                        ms_matrix_tmp <- strsplit(ms_original, split="")
                        ms_matrix_tmp <- t(matrix(as.numeric(unlist(ms_matrix_tmp)), ncol=ss*npop))
                        # Get the duplicated entries
                        dups <- which(duplicated(positions_tmp))
                        
                        # remove the sites that are monomorphic
                        freqSites_pop <- colSums(ms_matrix_tmp)
                        mono <- which(freqSites_pop==0 | freqSites_pop== (ss*npop)) 
                        
                        # sites to remove
                        rm_sites <- union(mono, dups)
                        
                        # remove the sites with two mutations, keeping only one of them
                        list(p=positions_tmp[-rm_sites], m=ms_matrix_tmp[,-rm_sites])
                      }
                      
                      # get the positions and ms_matrix as a list
                      positions_original <- lapply(l, function(x) x$p)
                      ms_matrix <- lapply(l, function(x) x$m)
                      # aux for index of pop
                      aux_pop <- list(c(1:ss),c((ss+1):(2*ss))) 

                      ###########################################################
                      # Estimate LD for each population, window and simulation
                      mean_r2_pop <- list()
                      # go for each population
                      for(pop in 1:npop) {
                        # tic()
                        mean_r2 <- foreach(sim=1:nsims) %dopar% {
                          # Get the frequency
                          freqSites_pop <- colSums(ms_matrix[[sim]][aux_pop[[pop]],, drop=FALSE])
                          # sum(is.na(freqSites))
                          
                          # get the sites to remove (those that have freq=0 or freq=1)
                          rm_sites <- which(freqSites_pop==0 | freqSites_pop==ss)
                          # get the ms without these positions
                          ms_matrix_pop <- ms_matrix[[sim]][aux_pop[[pop]],-rm_sites, drop=FALSE]
                          # str(ms_matrix_pop)
                          # discard the monomorphic positions and sites with more than two mutations
                          positions_pop <- positions_original[[sim]][-rm_sites]
                          # str(positions_pop)
                          
                          # check that the length of positions and that the number of columns of ms_matrix is the same
                          if( (ncol(ms_matrix_pop) != length(positions_pop)) ) {
                            stop("Error: positions has a different size than ms matrix.")
                          }
                          
                          # check that there are no duplicates left
                          dups <- duplicated(positions_pop)
                          if( sum(dups) > 0) {
                            positions_pop <- positions_pop[!dups]
                            if(sum(duplicated(positions_pop))>0) {
                              stop("Error: there are still duplicated positions.")
                            }
                          }
                          
                          # Get the index of the window of each site
                          index_window <- findInterval(positions_pop, seq_windows)
                          # Vector with the nuumber of windows for each case
                          window_id <- 1:(length(seq_windows)-1)
                          
                          # initialize the variable to save the mean_r2
                          vapply(window_id,FUN=function(x) { get_meanr2(ms_matrix_pop,index_window,x) }, FUN.VALUE = numeric(1))
                        }
                        # toc()
                        # save the results into a list
                        mean_r2_pop[[pop]] <- matrix(unlist(mean_r2), ncol=nsims)
                      }
                      
                      # Save RDS file with mean_r2 for each pop
                      saveRDS(mean_r2_pop, file=paste("./",folder[scenario], "/sumstat/r2_mean_sd", tag,sep=""))
                      
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
