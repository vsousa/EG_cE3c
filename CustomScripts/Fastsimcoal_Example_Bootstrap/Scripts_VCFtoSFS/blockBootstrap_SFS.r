# BOOTSTRAP_RESAMPLING 
#' Script that reads the resampled_geno_block with genotype info
#' and builds resampled bootstrap replicates, by resampling blocks with replacement

#' The command line has the following arguments:
#' @param foldertag tag for the folder with SFS
#' @param block_length length of blocks used
#' @param ind_threshold string with number of diploid individuals per pop (P1, P2, P3), e.g. "5,2,3"
#' @param dist_threshold minimum mean distance between consecutive SNPs in a good block
#' @param numboot number of bootstrap replicates
#' @param mafOrDaf string with "m" for MAF (minor allele frequency) and "d" for DAF (derived allele frequency)
#' @param randomSNP T for sampling 1 random SNP per block. F to sample all linked SNPs per block.
#' @param seed seed of random number generator
#' Example of use:
#' Rscript blockBootstrap_SFS.r block_SFS 10000 "10,10,10" 2 10 m F 626325
# REQUIRES THE FOLLOWING FILES
# - functions_resample.r with definition of functions
# read the command line arguments

args <- commandArgs(trailingOnly = TRUE)

foldertag <- as.character(args[1]) # tag for the folder with the resulting resampled block SFS
block_length <- as.numeric(args[2]) # define the length (distance among positions for defining a block)
ind_threshold_s <- as.character(args[3]) # string with number of diploid individuals per pop "(P1, P2, P3)"
dist_threshold <- as.numeric(args[4]) # minimum mean distance between consecutive SNPs in a good block
numboot <- as.numeric(args[5]) # number of bootstrap replicates
mafOrDaf <- as.character(args[6]) # string with "m" for MAF (minor allele frequency) and "d" for DAF (derived allele frequency)
randomSNP <- as.logical(args[7]) # T for sampling 1 random SNP per block. F to sample all linked SNPs per block.
seed <- as.numeric(args[8]) # seed of random number generator

# check if all arguments are passed to the function
if(length(args)!=8) {
  stop("Not all necessary arguments given as input. The command line needs to have the following arguments in the following order:
       foldertag: tag for the folder with SFS
       block_length: length of blocks used
       ind_threshold: string with number of diploid individuals per pop (P1, P2, P3), e.g. \"5,2,3\"
       dist_threshold: minimum mean distance between consecutive SNPs in a good block
       numboot: number of bootstrap replicates
       mafOrDaf: string with m for MAF (minor allele frequency) and \"d\" for DAF (derived allele frequency)
       randomSNP: T for sampling 1 random SNP per block. F to sample all linked SNPs per block.
       seed: seed of random number generator
       Example of use:
       Rscript blockBootstrap_SFS.r block_SFS 10000 \"3,1,7\" 2 10 m F 626325
       ")
}


# DEBUG
# tag of vcf file name
#foldertag <- "block_SFS" # tag for the folder with the resulting resampled block SFS
#block_length <- 200 # define the length (distance among positions for defining a block)
#ind_threshold_s <- "3,4,3" # number of diploid individuals per pop (P1, P2, P3)
#dist_threshold <- 2 # minimum mean distance between consecutive SNPs in a good block
#numboot <- 10 # number of bootstrap replicates
#mafOrDaf <- "d" # string with "m" for MAF (minor allele frequency) and "d" for DAF (derived allele frequency)
#randomSNP <- T # T for sampling 1 random SNP per block. F to sample all linked SNPs per block.
#seed <- 3483476

# Print setting
print("Performing Bootstrap resampling with the following arguments")
print(paste("Folder with SFS resampled files:", foldertag))
print(paste("Block length:", block_length))
print(paste("Minimum number of individuals per pop:", ind_threshold_s))
print(paste("Minimal distance between consecutive SNPs:", dist_threshold))
print(paste("Number of bootstrap replicates:", numboot))
print(paste("MAF (m) or DAF (d)? ", mafOrDaf))
print(paste("Sample 1 random SNP per block? ", randomSNP))
print(paste("Seed for random number generator: ", seed))

# set the seed of random number generator
set.seed(seed)

# get the individual threshold with sample sizes per pop as a vector
ind_threshold <- as.numeric(unlist(strsplit(ind_threshold_s, ",")))
npop <- length(ind_threshold)

# load the functions used in this script
source("./Scripts_VCFtoSFS/functions_resample.r")
options(warn=-1)

# FILES AND FOLDERS -------------------------------
# filename output file tag
filename <- paste(block_length, "bp_dist", dist_threshold,"_ind",paste(ind_threshold,collapse="_"),sep="")
folderblock_name <- paste(foldertag,"_",filename,sep="")
setwd(folderblock_name)
# Define the bootstrap folder
if(randomSNP==TRUE) {
  bootstrapfolder <- paste("Bootstrap_1SNP_",filename,sep="") 
} else {
  bootstrapfolder <- paste("Bootstrap_",filename,sep="") 
}
# create folder to save the bootstrap SFS
dir.create(paste("./", bootstrapfolder, sep=""))

# LOAD THE RESAMPLED_GENO_BLOCK
bootstrap_filename <- paste(filename,"_originaldata", sep="")
# load variable resampled_geno_block
load(file=bootstrap_filename)
print(paste("Loaded file with resampled data from file ", folderblock_name, "/",bootstrap_filename, 
            "with data from ", length(resampled_geno_block), " blocks.",sep=""))
# check
if(!exists("resampled_geno_block")) {
  stop("Could not load file with resampled genotype data used to build original SFS.")
}
 
# number of blocks
numblocks <- length(resampled_geno_block)
totalsites <- numeric(numboot) 

# For each bootstrap replicate
setwd(bootstrapfolder)
for(boot in 1:numboot) {
  
  if(randomSNP==TRUE) {
    # filename of boot tag
    filename_boot <- paste("sfs_1snp_",filename, "_boot", sep="")
    # sample one SNP per bootstrap replicate
    # and get the jointsfs
    jointsfs <- getsfs_singlesnpblock(resampled_geno_block, ind_threshold)
  } else {
    # sample the scaffold index with replacement
    boot_sc <- sample(c(1:numblocks), size=numblocks, replace=T)
    # find the blocks such that the number of polymorphic sites is between 0.99 and 1.01 of the total length
    #eval1 <- cumsum(selected_length_sc[boot_sc])<(sum(selected_length_sc)*1.01)
    #eval2 <- !cumsum(selected_length_sc[boot_sc])>(sum(selected_length_sc)*0.99)
    #index_last_sc <- mean(c(sum(eval1),sum(eval2)))
    #boot_sc <- boot_sc[1:index_last_sc]
    boot_sc <- sort(boot_sc)
    # get the bootstrap resampled geno
    boot_geno <- resampled_geno_block[boot_sc]
    # filename of boot tag
    filename_boot <- paste("sfs_",filename, "_boot", sep="")
    # if keeping all sites from each block
    # for each list, get the genotypes for pop i and add rbind then into single matrix
    boot_geno_merged <- sapply(c(1:npop), function(i) {do.call(rbind,lapply(boot_geno, function(x) x[[i]]))})
    # get the joint sfs as a multidimensional array
    jointsfs <- getjointsfs(genopops=boot_geno_merged, nindpops=ind_threshold)
  }
  
  # get total number of sites (TO CHECK: this could include monomorphic sites!!)
  totalsites[boot] <- sum(jointsfs)
  #str(boot_geno_merged)
  print(paste("Boot ", boot,": number of blocks ", numblocks, ", Number of sites (can still include monomorphic sites) is ", totalsites[boot], sep=""))
  
  # check the dimensions of the joint SFS
  if(!(sum(dim(jointsfs)==((ind_threshold*2)+1))==npop)) {
    stop("Incorrect dimensions of joint SFS after resampling.")
  }
  
  # Discard monomorphic sites
  jointsfs_nomon <- jointsfs
  entries_mono <- matrix(c(rep(1,times=npop),(ind_threshold*2)+1), nrow=2, byrow=TRUE)
  jointsfs_nomon[entries_mono] <- 0
  
  # Get the MAF (folded SFS)
  mafsfs <- derivedtomaf(jointsfs)
  mafsfs_nomon <- derivedtomaf(jointsfs_nomon)
  
  # save multi dimensional SFS
  dir.create("jsfs")
  if(mafOrDaf=="d") {
    # with monomorphic sites
    outfilename <- paste("./jsfs/",filename_boot, "_DSFS", boot,".obs",sep="")
    write_multisfs(jointsfs, outfilename)
    # without monomorphic sites
    outfilename <- paste("./jsfs/",filename_boot, "_nm_DSFS", boot,".obs",sep="")
    write_multisfs(jointsfs_nomon, outfilename)
  } else if(mafOrDaf=="m") {
    # with monomorphic sites
    outfilename <- paste("./jsfs/",filename_boot, "_MSFS", boot,".obs",sep="")
    write_multisfs(mafsfs, outfilename)
    # without monomorphic sites
    outfilename <- paste("./jsfs/",filename_boot, "_nm_MSFS", boot,".obs",sep="")
    write_multisfs(mafsfs_nomon, outfilename)
  }
  
  # save the 2D pairwise-SFS
  dir.create("pair2D")
  if(mafOrDaf=="d") {
    # with monomorphic sites
    fileend <- "_jointDAFpop"
    marg_2dsfs <- getmarginal_pairwise_2dsfs(jointsfs)
    for(i in 1:length(marg_2dsfs$sfs2d)) {
      outfilename <- paste("./pair2D/",filename_boot, boot, fileend, marg_2dsfs$pairwise[[i]][2]-1,"_", marg_2dsfs$pairwise[[i]][1]-1,".obs",sep="")
      write_2dsfs(marg_2dsfs$sfs2d[[i]], outfilename)
    }
    # without monomorphic sites
    fileend <- "_nm_jointDAFpop"
    marg_2dsfs_nomon <- getmarginal_pairwise_2dsfs(jointsfs_nomon)
    for(i in 1:length(marg_2dsfs_nomon$sfs2d)) {
      outfilename <- paste("./pair2D/",filename_boot, boot, fileend, marg_2dsfs_nomon$pairwise[[i]][2]-1,"_", marg_2dsfs_nomon$pairwise[[i]][1]-1,".obs",sep="")
      write_2dsfs(marg_2dsfs_nomon$sfs2d[[i]], outfilename)
    }
  } else {
    # with monomorphic sites
    fileend <- "_jointMAFpop"
    marg_2dsfs <- getmarginal_pairwise_2dsfs(mafsfs)
    for(i in 1:length(marg_2dsfs$sfs2d)) {
      outfilename <- paste("./pair2D/",filename_boot, boot, fileend, marg_2dsfs$pairwise[[i]][2]-1,"_", marg_2dsfs$pairwise[[i]][1]-1,".obs",sep="")
      write_2dsfs(marg_2dsfs$sfs2d[[i]], outfilename)
    }
    # without monomorphic sites
    fileend <- "_nm_jointMAFpop"
    marg_2dsfs_nomon <- getmarginal_pairwise_2dsfs(mafsfs_nomon)
    for(i in 1:length(marg_2dsfs_nomon$sfs2d)) {
      outfilename <- paste("./pair2D/",filename_boot, boot, fileend, marg_2dsfs_nomon$pairwise[[i]][2]-1,"_", marg_2dsfs_nomon$pairwise[[i]][1]-1,".obs",sep="")
      write_2dsfs(marg_2dsfs_nomon$sfs2d[[i]], outfilename)
    }
  }
    
  # save marginal 1D SFS 
  dir.create("sfs1D")
  if(mafOrDaf=="d") {
    # with monomorphic sites
    fileend <- "_DAFpop"
    marg_1dsfs <- getmarginal_1dsfs(jointsfs)
    for(i in 1:length(marg_1dsfs)) {
      outfilename <- paste("./sfs1D/",filename_boot, boot, fileend, i-1,".obs",sep="")
      write_1dsfs(marg_1dsfs[[i]], outfilename)
    }
    # without monomorphic sites
    fileend <- "_nm_DAFpop"
    marg_1dsfs_nomon <- getmarginal_1dsfs(jointsfs_nomon)
    for(i in 1:length(marg_1dsfs)) {
      outfilename <- paste("./sfs1D/",filename_boot, boot, fileend, i-1,".obs",sep="")
      write_1dsfs(marg_1dsfs[[i]], outfilename)
    }
  } else if(mafOrDaf=="m") {
    # with monomorphic sites
    fileend <- "_MAFpop"
    marg_1dsfs <- getmarginal_1dsfs(mafsfs)
    for(i in 1:length(marg_1dsfs)) {
      outfilename <- paste("./sfs1D/",filename_boot, boot, fileend, i-1,".obs",sep="")
      write_1dsfs(marg_1dsfs[[i]], outfilename)
    }
    # without monomorphic sites
    fileend <- "_nm_MAFpop"
    marg_1dsfs_nomon <- getmarginal_1dsfs(mafsfs_nomon)
    for(i in 1:length(marg_1dsfs)) {
      outfilename <- paste("./sfs1D/",filename_boot, boot, fileend, i-1,".obs",sep="")
      write_1dsfs(marg_1dsfs[[i]], outfilename)
    }
  }
  
}
setwd("../../")
