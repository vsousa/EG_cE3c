#Vitor Sousa - 11/10/2023

#Noticed that the indels are not aligned. Easier to use a package

# install.packages('vcfR')
#install.packages("insect")
#install.packages("ape")
library(vcfR)
library(ape)
library(seqinr)
# library(insect)

# read VCF file
vcf <- read.vcfR("mtdna_data_withmono_freebayes.vcf", verbose = FALSE )
# print info about the vcf
vcf

# Show the meta information
head(vcf@meta)

# Fix region (info for each site and sample)
head(getFIX(vcf))
# What info do we have for each site?

# Get DP field
dp <- extract.gt(vcf, element = "DP", as.numeric = TRUE)

#catarina bernardo - 21/12/2023
# calculating mean dp per individual
individual_codes <- read.table("C:/Users/cayna/OneDrive - Universidade de Lisboa/CatarinaBernardo/vcf_vstest/individual_codes.txt", quote="\"", comment.char="")

ind_names_all <- colnames(dp1)
inds_to_keep <- individual_codes$V1
dp1 <- as.matrix(dp)
dp1
evalinds <- ind_names_all %in% inds_to_keep 
dpind<- colMeans(dp, na.rm = T)

dpind_to_keep <- dpind[evalinds]
mean(dpind_to_keep)
# distribution of DP across sites for individual 2
hist(dp[,2], xlab="DP depth of coverage", main=colnames(dp)[2])
# plot boxplot distribution of DP for each individuals
boxplot(dp, outline=FALSE)

# MERGE THE REF AND ALT ALLELES
# for each position, get the reference and alternative alleles
ref <- getREF(vcf)
alt <- getALT(vcf)

# for the alt, split the alleles according to the comma
alt_all <- strsplit(alt, split=",")
str(alt_all)

# CREATE A LIST COMBINING THE REF AND ALT ALLELES
alleles_nogaps <- list()
for(i in 1:length(ref)) {
  alleles_nogaps[[i]] <- c(ref[[i]], alt_all[[i]])  
}

# get the gt command
gt <- extract.gt(vcf, element="GT", as.numeric = FALSE)
str(gt)

# check that the number of positions in alleles_nogaps is the same as gt
if(nrow(gt)!=length(alleles_nogaps)) {
  stop("Error in length of gt or alleles_nogaps!")
}

# Go through each individual and get the corresponding fasta file 
seq_ind <- matrix(NA,ncol=ncol(gt),nrow=nrow(gt))
for(i in 1:ncol(gt)) { # loop for each individual
  for(j in 1:nrow(gt)) { # loop for each site
    # the genotypes are coded as 0, 1, 2, ... - add 1 such that ref is 1, alt is 2, and so on 
    all_index <- as.numeric(gt[j,i])+1
    seq_ind[j,i] <- alleles_nogaps[[j]][all_index]
  }
}  

# filter each individual, replacing sites with less than dp by "n"
str(seq_ind)
str(dp)

min_dpfilter <- 30
seq_ind[dp<min_dpfilter] <- "N"
seq_ind[is.na(seq_ind)] <- "N"

# the sequence of each individuals will be the concatenation of each column
dna <- apply(seq_ind, MARGIN = 2, function(col) paste(col, collapse=""))
str(dna)
# open files
output <- file("test_vsscript.fa", "w")
# name of individuals
indlabel <- colnames(gt)
# add lines for each individual
for(i in 1:ncol(gt)) {
  # write individual ID
  writeLines(paste(">", indlabel[i], sep=""), con=output)
  # write the DNA sequence
  writeLines(dna[i], con=output)
}
# close the file
close(output)


seqinr::write.fasta(dna, file.out="test_vsscript_string.fa", names=colnames(gt))

