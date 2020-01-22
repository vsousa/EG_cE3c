#' Script that reads the genotype (GT) and (CHRMPOS) fields from a VCF
#' and builds the SFS, based on a block resmplaing approach

 #' Example of use:
#' Rscript BlockSFS_jan2020.r myvcf_file indpopinfo.txt block_SFS 10000 "10,10,10" 2 T 162512;
# REQUIRES THE FOLLOWING FILES
# - functions_resample.r with definition of functions
# read the command line arguments
args <- commandArgs(trailingOnly = TRUE)

filename_vcf <- as.character(args[1]) # tag for the file with genotypes (filename_vcf.GT), chromosome position info (filename_vcf.CHRMPOS), etc.
indpopinfofilename <- as.character(args[2]) # filename of the file with IndPopInfo (1st column indID, 2nd column popID)
foldertag <- as.character(args[3]) # tag for the folder with the resulting resampled block SFS
block_length <- as.numeric(args[4]) # define the length (distance among positions for defining a block)
ind_threshold_s <- as.character(args[5]) # string with number of diploid individuals per pop "(P1, P2, P3)"
dist_threshold <- as.numeric(args[6]) # minimum mean distance between consecutive SNPs in a good block
randomS <- as.logical(args[7]) # string with TRUE for randomly sampling individuals. FALSE for sampling at each block always the individuals with less missing data.
seed <- as.numeric(args[8]) # number with seed of random number generator

# check if all arguments are passed to the function
if(length(args)!=8) {
  stop("Not all necessary arguments given as input. The command line needs to have the following arguments in the following order:
       filename_vcf: filename tag for the file with genotypes (filename_vcf.GT), chromosome position info (filename_vcf.CHRMPOS), etc.
       indpopinfofilename: file with the ID of individuals and corresponding population
       foldertag: tag for the folder with output SFS
       block_length: define the length (distance among positions for defining a block)
       ind_threshold: string with number of diploid individuals per pop (P1, P2, P3), e.g. \"5,2,3\"
       dist_threshold: minimum mean distance between consecutive SNPs in a good block
       randomS: boolean. TRUE for randomly sampling individuals. FALSE for sampling at each block always the individuals with less missing data.
       seed: number with seed of random number generator.
       Example of use:
       Rscript BlockSFS_jan2020.r myvcf_file indpopinfo.txt block_SFS 10000 \"4,2,3\" 2 T 162512
       ")
}

# load the functions used in this script
source("./Scripts_VCFtoSFS/functions_resample.r")
options(warn=-1)

# DEBUG
# tag of vcf file name
# filename_vcf <- "final_filtered_3scaf" # tag for the file with genotypes (filename_vcf.GT), chromosome position info (filename_vcf.CHRMPOS), etc.
# indpopinfofilename <- "IndPopInfo_OffOnSH.txt" # filename of the file with IndPopInfo (1st column indID, 2nd column popID)
# foldertag <- "block_SFS" # tag for the folder with the resulting resampled block SFS
# block_length <- 1000 # define the length (distance among positions for defining a block)
# ind_threshold_s <- "1,2,3" # number of diploid individuals per pop (P1, P2, P3)
# dist_threshold <- 2 # minimum mean distance between consecutive SNPs in a good block
# randomS <- TRUE # TRUE for randomly sampling individuals. FALSE for sampling at each block always the individuals with less missing data.
# seed <- 236235  # seed for random number generator

# Print setting
print("Performing resampling with the following arguments")
print(paste("Tag for VCF file:",filename_vcf))
print(paste("File with indID and pop info:", indpopinfofilename))
print(paste("Folder where resulting files are saved:", foldertag))
print(paste("Block length:", block_length))
print(paste("Minimum number of individuals per pop:", ind_threshold_s))
print(paste("Minimal distance between consecutive SNPs:", dist_threshold))
print(paste("Sampling individuals randomly? (if FALSE individuals with less missing data are always chosen)", randomS))
print(paste("Seed for random number generator: ", seed))


# get the individual threshold with sample sizes per pop as a vector
ind_threshold <- as.numeric(unlist(strsplit(ind_threshold_s, ",")))

# set the seed of random number generator
set.seed(seed)

################################################################
# Read info about individual ID, Sampling location, phenotypes
################################################################
{# Read the pop map, file where first column is the ind ID and 2nd the POP
indpopinfo <- read.table(indpopinfofilename, stringsAsFactors = FALSE, header=TRUE)  
# Read individual IDs as in VCF file
inds_file <- paste("indsIn",filename_vcf, sep="");
indidvcf <- scan(inds_file, what="character");
# number of individuals in vcf file
nindvcf  <- length(indidvcf) 

# Get index of individuals common in vcf and for which we have pop info (i.e. individuals found in both files)
# commonind has the index
# $indexinpopinfo : index in indpopinfo of the individuals found in both the VCF and indpopinfo matrix
# $indexretainvcf : index in indidvcf of the individuals found in both the VCF and indpopinfo matrix
commonind <- get_indpopinfo_indvcf(indpopinfo, indidvcf)

# vcf.indpopinfo keeps the information about the individuals in the vcf file
# that are sorted in the same order as they are in the file indpopinfo
vcf.indpopinfo <- indpopinfo[commonind$indexinpopinfo, ]

# check that the individuals in vcf.indpopinfo have the same order as in indpopinfo
if(sum(vcf.indpopinfo[,1]!=indidvcf[commonind$indexretainvcf])>0) {
  stop("Error! The order of individuals is not the same in indpopinfo and indidvcf!")
}

# number of individuals
nind <- length(vcf.indpopinfo[,1])

# Get number of populations and pop names
npop <- length(unique(vcf.indpopinfo[,2]))
pop.names <- unique(vcf.indpopinfo[,2])

# Print info for user
print(paste("Read ", indpopinfofilename, ", with ", nrow(indpopinfo), " individuals, and ", npop, " populations.",sep=""))
print(paste("Read ", inds_file, ", with ", nindvcf, " individuals.",sep=""))
print(paste("Found ", nind, " individuals in common between ", indpopinfofilename, " and ", inds_file, sep=""))
print(paste("NOTE: individuals in VCF will be sorted according to order in ", indpopinfofilename, sep=""))
}

##################################################
# Read data files
##################################################
{# Read the genotype data
gt_file <-  paste(filename_vcf ,".GT", sep="");
nind_tmp <- length(scan(gt_file, nlines=1));
tgendata_original <- matrix(scan(gt_file, na.strings="-1"), ncol=nind_tmp, byrow=T);
# check
if(nind != ncol(tgendata_original)) {
  stop("Error! Wrong number of individuals in GT file!")
};
nsites <- nrow(tgendata_original)
print(paste("Read genotypes from ", gt_file, ", with ", nsites, " sites, and ", nind_tmp, " individuals.",sep=""))
rm(nind_tmp)

# Get a matrix with the missing data (i.e.genotypes that did not pass filters)
# this is important because VCF tools only updates the genotypes and not the associated DP, AD and PL
# missingdata <- is.na(tgendata_original)  
# ## Read the depth
# dp_file <- paste(filename_vcf ,".DP", sep="");
# nind_tmp <- length(scan(dp_file, nlines=1));
# dp <- matrix(scan(dp_file), ncol=nind_tmp, byrow=T);
# # check
# if(nind != ncol(dp)) {
#   stop("Error! Wrong number of individuals in DP file!")
# };
# print(paste("Read depth of coverage ", dp_file, ", with ", nrow(dp), "sites, and ", nind_tmp, " individuals.",sep=""))
# 
# # Replace the DP outside filters by missing data
# dp[missingdata] <- NA
# # check the missing data is the same
# if(sum(is.na(dp)) != sum(is.na(tgendata_original))) {
#   stop("Error! Missing data entries differ between GT and DP fields!")
# }

# sort the tgendata keeping only the selected individuals, sorted according to order in indpopinfo
tgendata <- tgendata_original[,commonind$indexretainvcf]
rm(tgendata_original)
#dp <- dp[,commonind$indexretainvcf]

# Check that we only have entries 0,1,2 and NA
checkentries <- sum(tgendata==0,na.rm=T)+sum(tgendata==1,na.rm=T)+sum(tgendata==2,na.rm=T)+sum(is.na(tgendata))
if(checkentries!= ncol(tgendata)*nrow(tgendata)) {
  stop("Error in genotype GT file! There are entries different from 0,1,2,NA!")
}

# read the chromosome and position
chrmposfile <- paste(filename_vcf,".chrpos",sep="")
chrmpos <- matrix(scan(chrmposfile, what="character"), ncol=2, byrow=T)
# check that the number of sites is the same as in GT file
if(nrow(chrmpos)!=nsites) {
  stop("Error in chrmpos file! There are different number of sites in CHRM_POP and genotype GT file!")
}


# get the distance among SNPs
scaffolds <- unique(chrmpos[,1])
#str(scaffolds)
# print info
print(paste("Read CHROM POS file ", chrmposfile, ", with ", nrow(chrmpos), " sites. Found ", length(scaffolds), " chromosomes/scaffolds/contigs.",sep=""))

}


# ***************************************************************#
# Resample without missing data from target populations ##########
# ***************************************************************#

# filename output file tag
filename <- paste(block_length, "bp_dist", dist_threshold,"_ind",paste(ind_threshold,collapse="_"),sep="")

# create folder to save the results
folderblock_name <- paste(foldertag,"_",filename,sep="")
dir.create(folderblock_name)
setwd(folderblock_name)

# Assuming that all blocks have the same order of individuals, which must be the case
npop <- length(pop.names)
# Get a list of index of individuals from each pop
pop_ind_i <- sapply(1:npop, function(i) {which(indpopinfo[,2]==pop.names[i])})
# Get a vector with the pop of each individual
ind_inpop_i <- rep(1:npop, times=unlist(lapply(pop_ind_i, length)))
print(paste("Resampling without missing data: blocks_length=", block_length, ", distance_threshold=", dist_threshold, ", individuals per pop (", paste(ind_threshold, collapse = ","), ")",sep=""))

# Go through each scaffold to define blocks and select SNP with individuals without missing data at those blocks
# NOTE: Assuming that blocks are independent, different individuals are sampled in different blocks.
res_scaf <- sapply(1:length(scaffolds), function(i) {
  # define number of blocks
  select_sc <- chrmpos[,1]==scaffolds[i]
  # partition the scaffold into contiguous blocks with block_length based on polymorphic sites
  snp_position <- as.numeric(chrmpos[select_sc,2]) 
  # if there are more than one callable site per scaffold perform analysis
  if(length(snp_position)>0) {
    # select the sites from a given scaffold
    geno_sc <- tgendata[select_sc,,drop=F]
    # resample the genotypes per scaffold
    res <- resample_scaffold(geno_sc=geno_sc, 
                             snp_position=snp_position, 
                             randomInd=randomS, 
                             block_length=block_length, 
                             dist_threshold, ind_threshold, ind_inpop_i, pop_ind_i, npop)
  } else {
    res <- list(geno=NA,snp=0,distsnp=NA)
  } 
  res
})

# res_scaf is a list of lists()
if(length(res_scaf) != length(scaffolds)) {
  stop("Error! After resampling, the number of scaffolds is incorrect!")
}
print(paste("Number of scaffolds after resampling is", length(res_scaf)))

# merge the scaffolds into a single list
res_scaf <- unlist(res_scaf, recursive = FALSE)
print(paste("Number of blocks (across all scaffolds) after resampling is", length(res_scaf)))

# remove blocks that did not pass distance threshold between consecutive SNPs, i.e. those with distsnp=NA
# or that still contained missing data
bad_blocks <- sapply(res_scaf, function(x) {is.na(x$snp)}, simplify=TRUE)
if(sum(bad_blocks)>0) {
  # discard blocks without enough SNPs
  res_scaf <- res_scaf[-which(bad_blocks)]
}
print(paste("Number of BAD blocks removed:", sum(bad_blocks)))

# Get the multidimenstional SFS across all blocks
# 1st. Merge the geno list from each block
resampled_geno_block <- sapply(res_scaf, function(x) {x$geno}, simplify=FALSE)
# check if there is missing data
if(sum(is.na(resampled_geno_block))>0) {
  stop("Missing data found in genotypes after resampling!")
}
# for each list, get the genotypes for pop i and add rbind then into single matrix
resampled_geno <- sapply(c(1:npop), function(i) {do.call(rbind,lapply(resampled_geno_block, function(x) x[[i]]))})
print(paste("After resampling: Number of blocks is ", length(resampled_geno_block), ", Number of sites (can still include monomorphic sites) is ", nrow(resampled_geno[[1]]), sep=""))
# (TO CHECK: this could include monomorphic sites!! Discard monomorphic sites here!)

# 2.1 Get the joint sfs as a multidimensional array
jointsfs <- getjointsfs(genopops=resampled_geno, nindpops=ind_threshold)
# check the dimensions of the joint SFS
if(!(sum(dim(jointsfs)==((ind_threshold*2)+1))==npop)) {
  stop("Incorrect dimensions of joint SFS after resampling.")
}

# 2.2. Get the joint SFS sampling 1 SNP per block
jointsfs_1snp <- getsfs_singlesnpblock(resampled_geno_block, ind_threshold)
# check the dimensions of the joint SFS
if(!(sum(dim(jointsfs_1snp)==((ind_threshold*2)+1))==npop)) {
  stop("Incorrect dimensions of joint SFS after resampling 1 SNP per block.")
}

# 3rd. Discard monomorphic sites
jointsfs_nomon <- jointsfs
entries_mono <- matrix(c(rep(1,times=npop),(ind_threshold*2)+1), nrow=2, byrow=TRUE)
jointsfs_nomon[entries_mono] <- 0

# 4th. Get the MAF (folded SFS)
mafsfs <- derivedtomaf(jointsfs)
mafsfs_nomon <- derivedtomaf(jointsfs_nomon)
# 4.2. Get the MAF for the sfs with 1 random SNP per block
mafsfs_1snp <- derivedtomaf(jointsfs_1snp)

# SAVE FILES
foldersfs_name <- paste("SFS",filename,sep="_")
dir.create(foldersfs_name)

# a) Save the resampled geno block. This is used to generate the non-parametric bootstrap files
bootstrap_filename <- paste(filename,"_originaldata", sep="")
save(resampled_geno_block, file=bootstrap_filename)
print(paste("Saved file with resampled data (NEEDED FOR NON-PARAMETRIC BOOTSTRAP) in file ", folderblock_name, "/",bootstrap_filename, sep=""))

# b) Save the joint DSFS, and marginal 2D and 1D SFS with and without monomorphic sites
tagsfs <- paste(foldersfs_name,"/",filename,sep="")
write.sfs(jointsfs,tagsfs)
outfilename_1snp <- paste(tagsfs, "_1SNP",sep="")
write_multisfs(jointsfs_1snp, paste(outfilename_1snp, "_DSFS.obs", sep=""))

tagsfs_nomon <- paste(foldersfs_name,"/",filename,"_nomon",sep="")
write.sfs(jointsfs_nomon, tagsfs_nomon)

# c) Save the joint MAF, MAF marginal 1D SFS and MAF pairwise 2D SFS, with and without monomorphic sites
tagsfs <- paste(foldersfs_name,"/",filename,"_MAF",sep="")
write.sfs(mafsfs,tagsfs)
write_multisfs(mafsfs_1snp, paste(outfilename_1snp, "_MSFS.obs", sep=""))
tagsfs_nomon <- paste(foldersfs_name,"/",filename,"_MAF_nomon",sep="")
write.sfs(mafsfs_nomon, tagsfs_nomon)
print(paste("Saved SFS files in folder ", foldersfs_name, sep=""))

setwd("../")

