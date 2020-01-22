# Functions to random sample individuals from blocks
# Vitor Sousa 19/01/2020

#' GETSORTEDINDS_RANDOM
#' sample individuals at random from a block of SNPs without missing data for a given pop
#' INPUT
#'    @param eval_inds vector with index of individuals belonging to a given pop (NOTE: assuming that the index of individual is the same across all blocks)
#'    @param missingdata_ind vector with the proportion of missing data for each individual
#'    @param ind_threshold number of individuals to sample
#'    @return vector with index of individuals sampled
getsortedinds_random <- function(eval_inds, missingdata_ind, ind_threshold) {
  # small value to avoid having individuals with zero probability of being sampled
  nonzero <- 0.00001
  # transform missing data, such that individuals are sampled with prob 1-missingdata
  transform_md <- missingdata_ind
  #transform_md[missingdata_ind>0] <- 0.001 
  # check there are enough individuals
  if(length(eval_inds)<ind_threshold) {
    stop("Not enough individuals to sample!")
  } else {
    # sample individuals proportionally to its missing data
    sampled_inds <- sample(eval_inds, size=ind_threshold, replace = F, prob=(1-(transform_md[eval_inds]))+nonzero)  
  }
  sampled_inds
}


#' GETSORTEDINDS_DETERMINISTIC
#' check if all populations pass the threshold of individuals with data
#' INPUT
#'    @param eval_inds vector with index of individuals belonging to a given pop (NOTE: assuming that the index of individual is the same across all blocks)
#'    @param missingdata_ind vector with the proportion of missing data for each individual
#'    @param ind_threshold number of individuals to sample
#'    @return vector with index of individuals sampled
getsortedinds_deterministic <- function(eval_inds, missingdata_ind, ind_threshold) {
  # get index of individuals that belong to target pop
  eval_inds[order(missingdata_ind[eval_inds])[1:ind_threshold]]
}

#' RESAMPLE_BLOCK_DETERMINISTIC
#' resample deterministically a given number of individuals without missing data from each block.
#' INPUT
#'    @param geno_block nsites x nind genotype matrix with genotypes coded as 0, 1, 2, NA.
#'    @param ind_threshold vector of size npop with the number of individuals to sample from each pop without missing data
#'    @param ind_inpop_i vector of size nind with the pop index of each individual, e.g. indpop_i=c(1,1,1,2,2) mean that inds 1-3 belong to pop1, inds 4-5 belong to pop2
#'    @param pop_ind_i  list of size npop, each with a vector with index of individuals belonging to a given pop (NOTE: assuming that the index of individual is the same across all blocks)
#'    @param npop number of pops
#' OUTPUT
#'    @return resampled_geno_pop list of size npop with matrix of genotypes with resampled_genotypes
#'            npol_snps: number of snps
resample_block_deterministic <- function(geno_block, ind_threshold, ind_inpop_i, pop_ind_i, npop) {
  # initilize output
  result <- list(geno=NA, snp=NA)
  
  # get the missing data per individual at these sites
  missingdata_ind <- colSums(is.na(geno_block))/nrow(geno_block)
  # for each pop, get the index of ind_threshold individuals sorted according to missing data
  sortedinds <-sapply(1:npop, function(i) {
      getsortedinds_deterministic(pop_ind_i[[i]], missingdata_ind, ind_threshold[i])})
  # check that there is no missing data in the index of sorted individuals
  if(sum(is.na(unlist(sortedinds)))>0) {
    stop("Error sampling individuals! sortedinds with missing data!")
  }
  
  # get the missind data per site
  missingdata_site <- rowSums(is.na(geno_block[,unlist(sortedinds),drop=F]))/ncol(geno_block)
  
  # get the sites with no missing data
  sites_no_missingdata <- which(missingdata_site==0)
  
  # check that there is at least one site without missing data for the sampled individuals
  if(length(sites_no_missingdata)>0) {
    # check that there is no missing data in the selected individuals
    if(sum(is.na(geno_block[sites_no_missingdata,unlist(sortedinds),drop=F]))>0) {
      stop("Error! missing data still found!")
    } else {
      # get the genotypes for each pop for sites without missing data
      resampled_geno_pop <- lapply(sortedinds, function(x) geno_block[sites_no_missingdata, x, drop=FALSE])
      # number of polymorphic sites
      npol_snps <- length(sites_no_missingdata)
      # check that the number of sites is correct
      if(nrow(resampled_geno_pop[[1]])!=npol_snps) {stop("Wrong number of SNPs after resampling!")}
    }
    result <- list(geno=resampled_geno_pop, snp=npol_snps)
  } 
  result
}

#' RESAMPLE_BLOCK_RANDOM
#' resample randomly a given number of individuals without missing data from a block
#' INPUT
#'    @param geno_block nsites x nind genotype matrix with genotypes coded as 0, 1, 2, NA.
#'    @param ind_threshold vector of size npop with the number of individuals to sample from each pop without missing data
#'    @param ind_inpop_i vector of size nind with the pop index of each individual, e.g. indpop_i=c(1,1,1,2,2) mean that inds 1-3 belong to pop1, inds 4-5 belong to pop2
#'    @param pop_ind_i  list of size npop, each with a vector with index of individuals belonging to a given pop (NOTE: assuming that the index of individual is the same across all blocks)
#'    @param npop number of pops
#' OUTPUT
#'    @return resampled_geno_pop list of size npop with matrix of genotypes with resampled_genotypes
#'            npol_snps: number of snps
resample_block_random <- function(geno_block, ind_threshold, ind_inpop_i, pop_ind_i, npop) {
  # initilize output
  result <- list(geno=NA, snp=NA)
  
  # get the missing data per individual at these sites
  missingdata_ind <- colSums(is.na(geno_block))/nrow(geno_block)
  # for each pop, get the index of ind_threshold individuals sorted according to missing data
  sortedinds <- sapply(1:npop, function(i) {
    getsortedinds_random(eval_inds=pop_ind_i[[i]], missingdata_ind, ind_threshold[i])})
  # check that there is no missing data in the index of sorted individuals
  if(sum(is.na(unlist(sortedinds)))>0) {
    stop("Error sampling individuals! sortedinds with missing data!")
  }
  
  # get the missind data per site
  missingdata_site <- rowSums(is.na(geno_block[,unlist(sortedinds),drop=F]))/ncol(geno_block)
  
  # get the sites with no missing data
  sites_no_missingdata <- which(missingdata_site==0)
  
  # check that there is at least one site without missing data for the sampled individuals
  if(length(sites_no_missingdata)>0) {
    # check that there is no missing data in the selected individuals
    if(sum(is.na(geno_block[sites_no_missingdata,unlist(sortedinds),drop=F]))>0) {
      stop("Error! missing data still found!")
    } else {
      # get the genotypes for each pop for sites without missing data
      resampled_geno_pop <- lapply(sortedinds, function(x) geno_block[sites_no_missingdata, x, drop=FALSE])
      # number of polymorphic sites
      npol_snps <- length(sites_no_missingdata)
      # check that the number of sites is correct
      if(nrow(resampled_geno_pop[[1]])!=npol_snps) {stop("Wrong number of SNPs after resampling!")}
    }
    result <- list(geno=resampled_geno_pop, snp=npol_snps)
  } 
  result
}


#' RESAMPLE_SCAFFOLD
#' function that given data from a scaffold, will split it into blocks, and for each block resample genotypes
#' INPUT:
#'    @param geno_sc : nsites in scaffold x nind genotype matrix with genotypes coded as 0, 1, 2, NA (NOTE: individuals must be sorted according to populations)
#'    @param snp_position : vector of size nsites with the position of sites within the scaffold
#'    @param randomInd TRUE if sampling of individuals is random. FALSE if sampling is deterministic, sampling always the individuals with less missing data from each block.  
#'    @param block_length : number with the length of blocks. The scaffold is divided into blocks of block length.
#'    @param dist_treshold : number with threshold distance. Only blocks with median distance above this between consecutive SNPs are kept.
#'    @param ind_threshold : vector of size npop with the number of individuals to sample from each pop without missing data
#'    @param ind_inpop_i : vector of size nind with the pop index of each individual, e.g. indpop_i=c(1,1,1,2,2) mean that inds 1-3 belong to pop1, inds 4-5 belong to pop2
#'    @param pop_ind_i   : list of size npop, each with a vector with index of individuals belonging to a given pop (NOTE: assuming that the index of individual is the same across all blocks)
#'    @param npop : number of pops
#' OUTPUT
#'    @return resampled: list of size nblocks with SNPs with. Each element of the list has:
#'    @return resampled_geno_pop: list of size npop with matrix of genotypes with resampled_genotypes
#'    @return npol_snps: vector with number of snps for each block
#'    @return distsnp: vector with number of snps for each block
#'    NOTE: The resampled list contains blocks that do not pass the dist_threshold. Those have 0 SNPs.
resample_scaffold <- function(geno_sc, snp_position, randomInd, block_length, dist_threshold, ind_threshold, ind_inpop_i, pop_ind_i, npop) {
  
  # Need be sure that positions are sorted ascendingly
  # order the geno_sc by increasing snp_position
  # it is possible that users give an unsorted matrix
  sorted_i <- order(snp_position)
  snp_position <- snp_position[sorted_i]
  geno_sc <- geno_sc[sorted_i,]
  
  # Divide the positions into blocks of block_length
  # aux contains the sequence of the initial position of each block
  aux <- floor(seq(min(snp_position),max(snp_position)+block_length,by=block_length))
  num_blocks_tmp <- length(aux)
  
  # check that the aux variable has a difference equal to block_length between all consecutive elements
  if(sum(aux[2:num_blocks_tmp]-aux[1:(num_blocks_tmp-1)]!=block_length)>0) {
    stop("Error! The length of the block is not correctly defined in aux!")
  }
  
  # Find blocks with SNPs
  block_bins_snps <- findInterval(snp_position, aux, rightmost.closed=T)
  snps_block <- table(block_bins_snps)
  id_blocks_snps <- as.numeric(names(snps_block))
  # Get a list of index of sites from each block
  selected_pos_sc <- sapply(1:length(id_blocks_snps), function(i) {which(block_bins_snps==id_blocks_snps[i])})
  # Get the number of snps in each block
  num_snps_block_sc <- unlist(lapply(selected_pos_sc, length))
  
  # Go through each block with SNPs
  # for each block, resample the individuals deterministically
  # The output of getdist_resample_block for each block is
  # - list with resampled genotypes res$geno
  # - list with number of snps
  # - list with dist_consecutiveSNPs
  resampled <- sapply(1:length(snps_block), function(i) {
                      getdist_resample_block(geno_block=geno_sc[selected_pos_sc[[i]],,drop=F], 
                                             snp_position_block=snp_position[selected_pos_sc[[i]]], 
                                             num_snps_block=num_snps_block_sc[i], 
                                             randomInd=randomInd,
                                             dist_threshold, ind_threshold, ind_inpop_i, pop_ind_i, npop)
                      
  }, simplify = FALSE)
  
  resampled
}


#' GETDIST_RESAMPLE_BLOCK
#' function that given data from a a blocks will get the pairwise distance between consecutive SNPs
#' If a block is considered good, it resamples genotypes from the individuals without missing data
#' maximizing the number of snps with data across all those individuals.
#' INPUT:
#' @param geno_block nsites in block x nind genotype matrix with genotypes in scaffold coded as 0, 1, 2, NA. 
#' @param snp_position_block vector of size nsites in block with the position of sites within the scaffold
#' @param num_snps_block number of snps in the block
#' @param randomInd TRUE if sampling of individuals is random. FALSE if sampling is deterministic, sampling always the individuals with less missing data from each block.  
#' @param dist_threshold number with threshold distance. Only blocks with median distance above this between consecutive SNPs are kept.
#' @param ind_threshold vector of size npop with the number of individuals to sample from each pop without missing data
#' @param ind_inpop_i vector of size nind with the pop index of each individual, e.g. indpop_i=c(1,1,1,2,2) mean that inds 1-3 belong to pop1, inds 4-5 belong to pop2
#' @param pop_ind_i list of size npop, each with a vector with index of individuals belonging to a given pop (NOTE: assuming that the index of individual is the same across all blocks)
#' @param npop number of pops
# OUTPUT
#' @return  resampled_geno_pop: list of size npop with matrix of genotypes with resampled_genotypes
#'          npol_snps: vector with number of snps for each block
getdist_resample_block <- function(geno_block, snp_position_block, num_snps_block, randomInd, dist_threshold, ind_threshold, ind_inpop_i, pop_ind_i, npop) {
  
  # initialize variables to NA. This is the output if block has no sites resampled
  dist_consecutiveSNPs <- NA
  res <- list(geno=NA, snp=NA) 
  
  # if there are polymorphic sites in this block
  if(num_snps_block>0) {
    # if there is at least one pairs of SNPS compute the distance among them
    if(num_snps_block>1) {
      # get the distance among consecutive SNPs
      dist_consecutiveSNPs <- snp_position_block[2:num_snps_block]-snp_position_block[1:(num_snps_block-1)]
    } else if (num_snps_block==1){ # set the distance into a very large value
      dist_consecutiveSNPs <- 99999
    }
    
    # Select only blocks where the median distance among SNPs is larger than dist_theshold
    if(sum(is.na(dist_consecutiveSNPs))==0) {
      if(median(dist_consecutiveSNPs)>dist_threshold) {
        # resample the data from the block
        if(!randomInd) {
          res <- resample_block_deterministic(geno_block, ind_threshold, ind_inpop_i, pop_ind_i, npop)  
        } else {
          res <- resample_block_random(geno_block, ind_threshold, ind_inpop_i, pop_ind_i, npop)  
        }
      }
    } 
  }
  list(geno=res$geno, snp=res$snp, distsnp=dist_consecutiveSNPs)
}

#' GETSFS_SINGLESNPBLOCK
#' obtain genotypes after sampling 1 random SNP per block, 
#' @param resampled_geno_block list of size nblocks. Each entry of the list has a list with three entries.
#'                 Each entre has a matrix with nsites x nind_i, 
#'                 where nind_i is the number of individuals in pop i. 
#'                 Note that nsites is equal in all matrices of the list, but nind_i can vary across populations.
#'                 NOTE: NO MISSING DATA ALLOWED!
#' @param nindpops : vector of size npop with the number of individuals per pop
#' @return multidimensional array with the derived joint SFS
getsfs_singlesnpblock <- function(resampled_geno_block, nindpops) {
  # number of pops
  npop <- length(nindpops)
  # initialize the joint SFS as a multi-dimensional array with zeros
  jointsfs <- array(0, dim=(nindpops*2)+1)
  # number of blocks
  nblock <- length(resampled_geno_block)
    
  # check that there is no missing data
  if(sum(is.na(resampled_geno_block))>0) {
    stop("There is missing data in the genotype matrix used to compute SFS!")
  } else {
    
  # get a matrix with ncol=npop and nrow=nblocks  
  allfreq_single <- sapply(resampled_geno_block, function(block) {
    
      # get the allfreq across pops
      allfreq <- sapply(block, function(x) rowSums(x[,,drop=FALSE]))
      # check that allfreq is a matrix and not a vector
      if(is.null(dim(allfreq))) {
        allfreq <- matrix(allfreq, ncol=npop)
      }
      # check that allfreq is a matrix with npop columns
      if(ncol(allfreq)!=npop) {
         stop("Error in computing allele freqs while sampling 1 SNP per block.")
      }

      # discard monomorphic sites - get the freq across pops
      freqPops <- rowSums(allfreq)
      fixed_derived <- sum((nindpops*2))
      snps_sample <- which(freqPops > 0 & freqPops < fixed_derived) 
      nsnps <- length(snps_sample)
      if(nsnps>1) {
        snp_i <- sample(snps_sample, size=1)  
      } else if(nsnps==1) {
        snp_i <- snps_sample
      } else {
        snp_i <- NA
      }
      if(is.na(snp_i)) {
        return(rep(NA, times=npop))
      } else {
        return(allfreq[snp_i,,drop=FALSE])   
      }
    })
    
    # check that the number of columns is the same as number of blocks
    if(!(ncol(allfreq_single)==nblock)) {
      stop("Error while resampling 1 SNP per block. Number of sites is different from number of blocks.")
    } 
  
    # remove the missing data
    rm_i <- which(colSums(is.na(allfreq_single))==npop)
    if(length(rm_i)>0) {
      allfreq_single <- allfreq_single[,-rm_i]  
    }
    print(paste("Number of blocks with 1 SNP after removing missing data is ", ncol(allfreq_single), sep=""))
    
    # compute the joint SFS
    # add 1 to allfreq, since allfreq=0 corresponds to the entry 1 of a muldimentsional array
    # transpose because each entry is a row
    allfreq_i <- t(allfreq_single+1)
  
    # Using table to get the number of snps with a given frequency
    # transform this into a string with the coordinates, separated by comma
    allfreq_s <- apply(allfreq_i, 1, function(x) paste(x, collapse = ","))
    # count the frequency of each string with the coordinates
    freq_spectrum <- table(allfreq_s)
    # get the index as strings
    index_s <- names(freq_spectrum) 
    # transform strings into matrix of indices
    index_i <- t(sapply(strsplit(index_s,split=","), as.numeric))
    
    # allfreq contains the index of the entries in the jointsfs array we need to update
    # NOTE: multidimensional arrays can be access with matrices with coordinates
    #       each row corresponds to an entry, each column to the entry in 1st, 2nd, 3rd... dimension
    jointsfs[index_i] <- freq_spectrum
    
    # NOTE: no need to add the monomorphic sites, because by definition only SNPs are sampled.
  }
  jointsfs
}




#' GETJOINTSFS
#' obtain the joint SFS given a list of genotype matrices, 
#' each entry in the list corresponding to the genotype data from a given population.
#' @param genopops list of size npop. Each entry of the list has a matrix with nsites x nind_i, 
#'                 where nind_i is the number of individuals in pop i. 
#'                 Note that nsites is equal in all matrices of the list, but nind_i can vary across populations.
#'                 NOTE: NO MISSING DATA ALLOWED!
#' @param nindpops : vector of size npop with the number of individuals per pop
getjointsfs <- function(genopops, nindpops) {
  # number of pops
  npop <- length(genopops)
  # initialize the joint SFS as a multi-dimensional array with zeros
  jointsfs <- array(0, dim=(nindpops*2)+1)
  # check number of populations
  if(npop!=length(nindpops)) {
    stop("Number of populations is incorrect in genotype matrix or in nindpops!")
  }
  # check that there is no missing data
  if(sum(is.na(genopops))>0) {
    stop("There is missing data in the genotype matrix used to compute SFS!")
  } else {
    # get the allele frequency as the sum across individuals at each site per pop
    allfreq <- sapply(genopops, function(x) rowSums(x[,,drop=FALSE]))  
    # add 1 to allfreq, since allfreq=0 corresponds to the entry 1 of a muldimentsional array
    allfreq <- allfreq+1
    
    # Using table to get the number of snps with a given frequency
    # transform this into a string with the coordinates, separated by comma
    allfreq_s <- apply(allfreq, 1, function(x) paste(x, collapse = ","))
    # count the frequency of each string with the coordinates
    freq_spectrum <- table(allfreq_s)
    # get the index as strings
    index_s <- names(freq_spectrum) 
    # transform strings into matrix of indices
    index_i <- t(sapply(strsplit(index_s,split=","), as.numeric))
    
    # allfreq contains the index of the entries in the jointsfs array we need to update
    # NOTE: multidimensional arrays can be access with matrices with coordinates
    #       each row corresponds to an entry, each column to the entry in 1st, 2nd, 3rd... dimension
    jointsfs[index_i] <- freq_spectrum
    
    # add the last entry (corresponding to fixed derived - monomorphic)
    # to the entry 0,0,0,...0 that corresponds to fixed ancestral - monomorphic
    entries_mono <- matrix(c(rep(1,times=npop),(nindpops*2)+1), nrow=2, byrow=TRUE)
    jointsfs[entries_mono[1,,drop=FALSE]] <- sum(jointsfs[entries_mono])
    jointsfs[entries_mono[2,,drop=FALSE]] <- 0
    
  }
  jointsfs
}

#' GETMARGINAL_PAIRWISE_2DSFS
#' obtain the marginal pairwise 2D-SFS given a joint multidimensional SFS.
#' @param jointsfs multidimensional array with the joint SFS
getmarginal_pairwise_2dsfs <- function(jointsfs) {
  # sample size per pop
  samplesize <- dim(jointsfs)
  # number of pops
  npops <- length(samplesize)
  
  # get index of pairwsise combinations
  pairwise_i <- combn(npops,2)
  # get the 2D sfs for each pairwise_i combination
  # using apply to a multidimensional array, summing the marginal of 2 dimensions (given by each column in pairwise_i)
  # where the rows refer to pop1 and cols to pop2
  # hence need to transpose the resulting matrix
  marg2dsfs <- apply(pairwise_i, 2, function(col) {t(apply(jointsfs,col,sum))})
  # output a list with the 2D-SFS and corresponding pairwise comparisons
  list(sfs2d=marg2dsfs, pairwise=pairwise_i)
}

#' GETMARGINAL_1DSFS
#' obtain the marginal 1D-SFS given a joint multidimensional SFS.
#' @param jointsfs multidimensional array with the joint SFS
getmarginal_1dsfs <- function(jointsfs) {
  # sample size per pop
  samplesize <- dim(jointsfs)
  # number of pops
  npops <- length(samplesize)
  # get the 1D sfs for each population
  marg1dsfs <- sapply(1:npops, function(i) {apply(jointsfs,i,sum)})
  # output a list with the 2D-SFS and corresponding pairwise comparisons
  marg1dsfs
}


# Function to add matrices using the appropriate sum operator
addMatrix <- function(x) Reduce("+", x);

#' WRITE_MULTISFS
#' save the multidimensional SFS in fastsimcoal2 input file obs SFS format
#' works for DAF and MAF
#' @param jointsfs multidimensional array with the joint SFS
#' @param outfilename string with the output file name
write_multisfs <- function(jointsfs, outfilename) {
  # sample size per pop
  samplesize <- dim(jointsfs)
  # number of pops
  npops <- length(samplesize)
  # convert multidimensional array to a single vector
  # 1st. need to permute the order of dimensions, as R starts writing to file looping from the first
  #      and for fastsimcoal2 input we need "0,0,0" "0,0,1" ... "0,0,n3" "0,1,0" "0,1,1" ... "0,1,n3" ... "n1,n2,n3" 
  multisfs <- aperm(jointsfs, c(npops:1))
  multisfs <- as.vector(multisfs)
  # write the input file for fastsimcoal
  write(paste("1 observation. No. of demes and sample sizes are on next line
              ", npops, " ", paste(samplesize-1, collapse=" "), sep=" "), file=outfilename)
  write(multisfs, file=outfilename, ncolumns = length(multisfs), append = T)
}

#' WRITE_2DSFS
#' save the 2D-SFS in fastsimcoal2 input file obs SFS format
#' works for DAF and MAF
#' @param sfs2d matrix with pairwise 2D-SFS, with rows corresponding to pop2 and columns to pop1
#' @param outfilename string with the output file name
write_2dsfs <- function(sfs2d, outfilename) {
  # sample size per pop (pop2, pop1)
  samplesize <- dim(sfs2d)
  # check that it is a matrix
  if(!(length(samplesize)==2)) {
    stop("function write_2dsfs only works for 2D matrices as input.")
  }
  # get the col and row names
  # note that the samplesize of pop1 is the ncol
  # note that the samplesize of pop2 is the nrow
  code_pop1 <- paste("d0_", 0:(samplesize[2]-1), sep="")
  code_pop2 <- paste("d1_", 0:(samplesize[1]-1), sep="")
    
  # write the input file for fastsimcoal
  write(paste("1 observation."), file=outfilename)
  write.table(sfs2d, file=outfilename, col.names = code_pop1, row.names = code_pop2, append = T, quote=FALSE)
}

#' WRITE_1DSFS
#' save the 1D-SFS in fastsimcoal2 input file obs SFS format
#' works for DAF and MAF
#' @param sfs1d vector with 1D-SFS
#' @param outfilename string with the output file name
write_1dsfs <- function(sfs1d, outfilename) {
  # sample size per pop
  samplesize <- dim(t(sfs1d)) # need to use t() because dim of vector is NULL
  # check that it is a vector
  if(!(sum(samplesize)==(length(sfs1d)+1))) {
    stop("function write_1dsfs only works for 1D vectors as input.")
  }
  # get the col names
  # note that the samplesize of pop1 is the ncol
  code_pop1 <- paste("d0_", 0:(samplesize[2]-1), sep="")
  
  # write the input file for fastsimcoal
  write(paste("1 observation."), file=outfilename)
  write.table(t(sfs1d), file=outfilename, col.names = code_pop1, row.names = F, append = T, quote=FALSE)
}


#' WRITE.SFS
#' save the joint, marginal 2D and marginal 1D SFS
#' @param jointsfs multidimensional array with the joint SFS
#' @param outfiletag_sfs string with the output file tag name
write.sfs <- function(jointsfs, outfiletag_sfs) {
  # save multi dimensional SFS
  outfilename <- paste(outfiletag_sfs, "_DSFS.obs",sep="")
  write_multisfs(jointsfs, outfilename)
  # get marginal 2D SFS
  marg_2dsfs <- getmarginal_pairwise_2dsfs(jointsfs)
  for(i in 1:length(marg_2dsfs$sfs2d)) {
    outfilename <- paste(outfiletag_sfs, "_jointDAFpop", marg_2dsfs$pairwise[2,i]-1,"_", marg_2dsfs$pairwise[1,i]-1,".obs",sep="")
    write_2dsfs(marg_2dsfs$sfs2d[[i]], outfilename)
  }
  marg_1dsfs <- getmarginal_1dsfs(jointsfs)
  for(i in 1:length(marg_1dsfs)) {
    outfilename <- paste(outfiletag_sfs, "_DAFpop", i-1,".obs",sep="")
    write_1dsfs(marg_1dsfs[[i]], outfilename)
  }
}


#' DERIVEDTOMAF
#' transforms derived joint multidimensional SFS into MAF SFS.
#' NOTE: minor allele frequency is computed across populations. 
#' The MAF is not defined within each pop nor pairwise comparison.
#' INPUT:
#'   @param jointsfs multidimensional array with the joint SFS
#' RETURN
#'   @return jointmafsfs multidimensional array with the joint MAF SFS
derivedtomaf <- function(jointsfs) {
  # get sample sizes
  n <- dim(jointsfs)
  
  # initilize with zeros a multidimensional array to save the folded SFS
  mafdsfs <- array(0,dim=n)
  
  # compute threshold frequency, corresponding to 50% of frequency across all samples
  # sum(n) gives the total sample size (-length(n) is used because we discard the zero entries, 1 for each pop)
  threshold_freq <- 0.5*(sum(n)-length(n))
  
  # get all possible number of combinations of frequencies across pops
  # this is done by geting a sequence from 1,...,n[1]; 1,...,n[2], etc.
  # then we apply expand.grid to that list
  comb <- as.matrix(expand.grid(sapply(n, function(x) 1:x)))
  # check
  if(!(length(jointsfs)==nrow(comb))) {
    stop("jointsfs has a different size from possible combinations")
  }
  
  # for each combination apply the function to get the MAF
  # note that we treat the combination as the vector with the index of the entries in the multidimensional array
  # NOTE: we need to use drop=FALSE because to get the index we need to use a matrix
  for(i in 1:nrow(comb)) {
    # get the frequency, which is obtained by substracting 1 because the first index is 1 and not zero
    freq <- sum(comb[i,]-1)
    if(freq < threshold_freq) {
      minor <- comb[i,,drop=FALSE]
      mafdsfs[minor] <- mafdsfs[minor]+jointsfs[minor]  
    } else if(freq==threshold_freq) {
      # if freq is 50% then we consider 50% minor and 50% major allele
      minor <- comb[i,,drop=FALSE]
      major <- n - comb[i,,drop=FALSE] + 1
      mafdsfs[minor] <- mafdsfs[minor]+(0.5*jointsfs[minor]+0.5*jointsfs[major]) 
    } else {
      major <- comb[i,,drop=FALSE]
      minor <- n - comb[i,,drop=FALSE] + 1
      mafdsfs[minor] <- mafdsfs[minor]+jointsfs[major]  
    }
  }
  mafdsfs
}

# Genotype to sequence
gentoseq <- function(x) {
  if(is.na(x)) {
    seq <- c(9,9)    
  } else {
    if(x==0) {
      seq <- c(0,0)
    } else if(x==1) {
      seq <- c(0,1)
    } else if(x==2) {
      seq <- c(1,1)
    }
  }  
  seq
}


#############################################################################
# Set of functions used to read and process VCF to SFS
# Vitor Sousa 17/01/2020
#############################################################################


####################################################
# READ FILES 
####################################################

#' READ_IND_POP_INFO
#' reads a file with the IDs, sampling location and distance for each individual
#' INPUT:
#'   @param filename : string with the name of the file
#' OUTPUT:
#'   @return data frame with nind x 4 columns with  indID, popID, OnOff, distance
read_ind_pop_info <- function(filename, phenotype) {
  indpopinfo <- read.table(paste(filename,sep=""), colClasses=c("character","character","numeric","numeric"), header=T) 
  indpopinfo
}


#' READ_IND_ID_VCF
#' read the individual ID labels in the VCF file
#' INPUT:
#'   @param filename : string with the name of the file
#' OUTPUT:
#'   @return data frame with nind x 3 columns with  indID, N_SITES, MEAN_DEPTH
read_ind_id_vcf <- function(filename) {
  indidvcf <- read.table(paste("./data/",filename,".idepth",sep=""), colClasses=c("character", "numeric", "numeric"), header=T)  
  indidvcf
}

#' GET_INDPOPINFO_INDVCF
#' get index of individuals for which we have pop info
#' this is for cases where the VCF file contains a different number of individuals tha the pop info
#' we only keep the individuals that have the same ID in the ind pop info and in the VCF file.
#' INPUT:
#'   @param indpopinfo : data frame with ind and pop info information - only requirement first column must have individual ID
#'   @param indidvcf   : vector with ind ID labels in VCF file - must have individual ID
#' OUTPUT:
#'   @return list with the $indexinpopinfo and $indexretainvcf
#'   with the index of individuals in both files for the indpopinfo and VCF file respectively.
#'   It prints the number of individuals in common.
get_indpopinfo_indvcf <- function(indpopinfo, indidvcf) {
  count <- 0
  indexinpopinfo <- numeric() # index of the indid individuals in indpopinfo
  indexretainvcf <- numeric() # index of individuals to retain from indid 
  for(i in 1:length(indidvcf)) {
    if (sum(indpopinfo[,1]==indidvcf[i]) > 0) {
      indexinpopinfo <- c(indexinpopinfo, which(indpopinfo[,1]==indidvcf[i]))
      indexretainvcf <- c(indexretainvcf, i)
      count <- count + 1 
    }
  }
  print(count)
  
  # Check that the order of individuals is the same
  if(sum(indpopinfo[indexinpopinfo,1]!=indidvcf[indexretainvcf])>0) {
    stop("Error!! In function get_indpopinfo_indvcf: indpopinfo ind id is not the same as indidvcf in the vcf file")    
  }
  
  list(indexinpopinfo=indexinpopinfo,indexretainvcf=indexretainvcf)  
}


#' Function to plot color bar
#' Copied from http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
#' example of call:
#' color.bar(colorRampPalette(c("light green", "yellow", "orange", "red"))(100), -1)
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}



#' CHECK_NA_NAN_INF
#' prints the number of NA, NAN and Inf values
#' INPUT
#'   @param x : matrix or numeric array to check
#' OUTPUT
#'  @return prints the number of NA, NAN and Inf values
check_na_nan_inf <- function(x) {
  print(paste("#NA=",sum(is.na(x))))
  print(paste("#NAN=",sum(is.nan(x))))
  print(paste("#Inf=",sum(is.infinite(x))))
}




#' SAMPLE_IND
#' sample individuals
#' INPUT:
#'   @param indgenotypes : vector of size nind with the genotypes of each individual for a given site (coded as 0,1,2,NA)
#'   @param numselect : integer with number of individuals to sample from
#' RETURN
#'   @return index of selected individuals (sampled from among non missing data)
sample_ind <- function(indgenotypes, numselect) {
  sample(which(!is.na(indgenotypes)), size=numselect, replace=FALSE)
}


#' GET_FIRSTINDPOPS
#' returns vector with index of individuals when pop changes, e.g. if pop of individuals is c(1,1,1,2,2,3), this function returns c(1,4,6)
#' INPUT
#'	@param popidind : numeric vector of size nind with the corresponding pop for each individual (sorted such that all inds from the same pop are consecutive)
#'	@param ninds : number of individuals
#' RETURN
#'	@return numeric vector of size npop with the index of the first individual from each population
get_firstindpops <- function(popidind, ninds) {	
  resindex <- 1 # the individual of first pop is by definition at index 1
  for(i in 2:ninds) {
    if(popidind[i]!=popidind[i-1]) {
      resindex <- c(resindex, i)
    }		
  }
  resindex
}



#' GETREPEATEDSITESINDEX
#' Function to removed repeated SNPs that might be in the VCF file. Outputs index of sites that are not repeated.
#' INPUT
#'   @param siteinfo : matrix with nsnps x 2 (or more cols), where the first column is the contig ID and the second the position
#' RETURNS
#'   @return numeric vector of size nsnpsToKeep with the index of the sites to keep
getrepeatedsitesindex <- function(siteinfo) {
  contigs <- unique(siteinfo[,1])  
  indexres <- c()  
  #   for(i in 1:length(contigs)) {
  #     # For each contig get the number of sites
  #     evaluate <- which(sites[,1]==contigs[i]) 
  #     evaluate <- evaluate[duplicated(sites[evaluate,2])]    
  #   }  
  indexres <- lapply(contigs, FUN=function(x){evaluate <- which(siteinfo[,1]==x); evaluate[duplicated(siteinfo[evaluate,2])]})
  unlist(indexres)
}


#This function creates a color scale for use with the image()
#function. Input parameters should be consistent with those
#used in the corresponding image plot. The "horiz" argument
#defines whether the scale is horizonal(=TRUE) or vertical(=FALSE).
image.scale <- function(z, zlim, col = rainbow(12), breaks, horiz=TRUE, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){ylim<-c(0,1); xlim<-range(breaks)}
  if(!horiz){ylim<-range(breaks); xlim<-c(0,1)}
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}


