# PROCESS ms files from SLiM simulations
# read ms files and output summary files 

# clean memory: remove all variables from memory
rm(list=ls())
# load used packages
library(stringr)
library(RColorBrewer)
library(doParallel)
registerDoParallel(cores=5)

# load files with functions to compute summary statistics
source("./Rscripts/computeStats.r")

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
folder <- rep("divsel3_500Kb_scaled", times=10000)
scenario <- 1

# COMBINATIONS OF PARAMETERS TESTED
# migration rate
migrate_v <-c("0.000", "0.00034", "0.0017", "0.0034")
# selection of mutation under divergent selection
selmutben_v <- c("0.000","0.0067", "0.01333" ,"0.02667", "0.05333", "0.0667", "0.1334")
# recombination rate
recrate_v <- c("2.5e-7")
# initial frequency of allele a at locus under divergent selection
freqmut_v <- c("0.001", "0.01", "0.1", "0.5")
# dominance of allele a at locus under divergent selection
dominanceben_v <- c("0.01","0.5")
# time of split of the two populations
time_v <- c("2000")
# TRUE if the name of files indicate the initial frequency
freqTag <- TRUE

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
                    if(scenario==1 & !dir.exists(paste("./",folder[scenario], "/sumstat_cb",sep=""))) {
                      dir.create(paste("./",folder[scenario], "/sumstat_cb",sep=""))
                    }
                    # folder sumstat is for results of all simulations
                    if(scenario==1 & !dir.exists(paste("./",folder[scenario], "/sumstat",sep=""))) {
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
                      aux_pop <- matrix(c(1:ss,(ss+1):(2*ss)), ncol=2) 
                      
                      # positions of the list where the mutation was not lost
                      benmut <- which(unlist(lapply(positions_original, function(x) sum(x==positionben)))==1)
                      nbenmuts <- length(benmut)
                      
                      # Save the position of the sims where the beneficial mutation was not lost
                      # write.table(c(length(benmut)/nsims, benmut), 
                      #             file=paste("./",folder[scenario], "/sumstat_cb/benmut", tag,sep=""), col.names = FALSE, row.names = FALSE)
                      
                      
                      # ############################################################
                      # Compute genome-wide average pi, dxy and pairwise FST for each simulation
                      # conditional on retaining selected allele a (i.e. only computing based on benmut)
                      mean_pi_dxy_fst <- foreach(sim=benmut) %dopar% {
                        getMeanFst_Pi_Dxy(t(ms_matrix[[sim]]), ss)
                      }
                      mean_sumstat <- matrix(unlist(mean_pi_dxy_fst), ncol=4, byrow = TRUE)
                      write.table(mean_sumstat, file=paste("./",folder[scenario], "/sumstat_cb/genome_sumstat_", tag,sep=""),
                                  row.names = F, col.names = c("pi_pop1","pi_pop2","dxy","fst"))
                      # remove variables
                      rm(mean_sumstat,mean_pi_dxy_fst)

                      # Compute genome-wide average pi, dxy and pairwise FST for each simulation
                      # all sims, irrespective of retaining beneficial allele
                      mean_pi_dxy_fst <- foreach(sim=1:nsims) %dopar% {
                        getMeanFst_Pi_Dxy(t(ms_matrix[[sim]]), ss)
                      }
                      mean_sumstat <- matrix(unlist(mean_pi_dxy_fst), ncol=4, byrow = TRUE)
                      write.table(mean_sumstat, file=paste("./",folder[scenario], "/sumstat/genome_sumstat_", tag,sep=""),
                                  row.names = F, col.names = c("pi_pop1","pi_pop2","dxy","fst"))
                      # remove variables
                      rm(mean_sumstat,mean_pi_dxy_fst)

                      
                      ############################################################
                      # Window estimates of pi, dxy and pairwise FST for each window and simulation
                      # conditional on retaining selected allele a (i.e. only computing based on benmut)
                      window_pi_dxy_fst <- foreach(sim=benmut) %dopar% {
                        # Get the index of the window of each site
                        index_window <- findInterval(positions_original[[sim]], seq_windows)
                        # Vector with index of windows
                        window_id <- 1:(length(seq_windows)-1)
                        
                        # initialize the variable to save the mean_r2
                        vapply(window_id,FUN=function(x) {getMeanFst_Pi_Dxy(t(ms_matrix[[sim]][,index_window==x]), ss) }, FUN.VALUE = numeric(4))
                      }
                      # create position for each window
                      aux <- seq_windows[1:(length(seq_windows)-1)]
                      xaxis <- aux[1:(length(aux))]+(aux[2]/2)
                      window_sumstat <- do.call("rbind",window_pi_dxy_fst)
                      write.table(t(window_sumstat), file=paste("./",folder[scenario],"/sumstat_cb/window_sumstat_", window.size, "_", tag,sep=""), 
                                  col.names = rep(c("pi_pop1","pi_pop2","dxy","fst"), times=nbenmuts), row.names = xaxis, quote = FALSE)
                      # remove variables
                      rm(window_sumstat, window_pi_dxy_fst)
                      
                      # Window estimates of pi, dxy and pairwise FST for each window and simulation
                      # all sims, irrespective of retaining beneficial allele
                      window_pi_dxy_fst <- foreach(sim=1:nsims) %dopar% {
                        # Get the index of the window of each site
                        index_window <- findInterval(positions_original[[sim]], seq_windows)
                        # Vector with index of windows
                        window_id <- 1:(length(seq_windows)-1)
                        
                        # initialize the variable to save the mean_r2
                        vapply(window_id,FUN=function(x) {getMeanFst_Pi_Dxy(t(ms_matrix[[sim]][,index_window==x]), ss) }, FUN.VALUE = numeric(4))
                      }
                      # create position for each window
                      aux <- seq_windows[1:(length(seq_windows)-1)]
                      xaxis <- aux[1:(length(aux))]+(aux[2]/2)
                      window_sumstat <- do.call("rbind",window_pi_dxy_fst)
                      write.table(t(window_sumstat), file=paste("./",folder[scenario],"/sumstat/window_sumstat_", window.size, "_", tag,sep=""), 
                                  col.names = rep(c("pi_pop1","pi_pop2","dxy","fst"), times=nsims), row.names = xaxis, quote = FALSE)
                      # remove variables
                      rm(window_sumstat, window_pi_dxy_fst)
                      
                    }  else {
                      stop(paste("File:", "./", folder[scenario],"/ms_", tag,"_", muttype,".txt not found!",sep=""))
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
