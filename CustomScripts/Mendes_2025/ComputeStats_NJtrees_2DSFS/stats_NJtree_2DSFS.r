# Script to:
# 1. compute FST, dxy and pi from 2D-SFS obtained with ANGSD, discarding sites with frequency of 0.5 in both populations
# 2. estimate a NJ tree from pairwise dxy distance matrix
# 3. estimate bootstrap support for each node of the NJ tree based on dxy, generating resampled 2D-SFS assuming all sites are independent
# The statistics are computed by masking (removing) sites of the 2D-SFS
# corresponding to the entry where all individuals in both populations are heterozygotes (i.e., entry with half the sample size in each pop)
# as this might represent either missing data or mapping errors.
# 
# Author: Vitor Sousa
# Last update: 23rd May 2025
# requires the following script:
# - computeStats_SFS.r
# requires the following packages:
# - xlsx
# - stringdist
# - ape
# - TreeTools
# - phangorn
# - ggplot2

# Requites the following files:
# - ordem_pops_coordenadas.xlsx excel file with the population order sorted according to latitude

# load the required scripts
source("computeStats_SFS.r")

# load required packages
library("xlsx")
library("stringdist")
library("ape")
library("TreeTools")
library("phangorn")
library("ggplot2")
# read the order of populations according to latitude, and info about the number of individuals per pop
popsample <- read.xlsx("ordem_pops_coordenadas.xlsx",sheetIndex = "ordem_pops")
str(popsample)

# DEBUG
# popsample <- read.table("popsample_debug.txt", header=TRUE)
# str(popsample)

# number of populations
npop <- nrow(popsample)

# number of bootstrap replicates
nboot <- 1000

# go through each pairwise combination and compute the statistics based on the 2D-SFS
# index of population comparisons
index <- combn(1:nrow(popsample),2)
npaircomparisons <- ncol(index) # number of pairwise comparisons

# initialize the variables to save the pi, dxy and FST
# statistics based on all sites
pi1 <- matrix(NA, nrow=npop, ncol=npop, dimnames = list(popsample[,1],popsample[,1])) # pi - average number of pairwise differences Pop1
pi2 <- matrix(NA, nrow=npop, ncol=npop, dimnames = list(popsample[,1],popsample[,1])) # pi - average number of pairwise differences Pop2
dxy <- matrix(NA, nrow=npop, ncol=npop, dimnames = list(popsample[,1],popsample[,1])) # dxy - average number of pairwise differences Pop1 and Pop2
fst <- matrix(NA, nrow=npop, ncol=npop, dimnames = list(popsample[,1],popsample[,1])) # fst between Pop1 and Pop2

# masked statistics - discarding the entry corresponding to sites where all individuals are heterozygotes - could be because of missing data or wrong mapping
pi1_masked <- matrix(NA, nrow=npop, ncol=npop, dimnames = list(popsample[,1],popsample[,1])) # pi - average number of pairwise differences Pop1
pi2_masked <- matrix(NA, nrow=npop, ncol=npop, dimnames = list(popsample[,1],popsample[,1])) # pi - average number of pairwise differences Pop2
dxy_masked <- matrix(NA, nrow=npop, ncol=npop, dimnames = list(popsample[,1],popsample[,1])) # dxy - average number of pairwise differences Pop1 and Pop2
fst_masked <- matrix(NA, nrow=npop, ncol=npop, dimnames = list(popsample[,1],popsample[,1])) # fst between Pop1 and Pop2
fst_angsd <- matrix(NA, nrow=npop, ncol=npop, dimnames = list(popsample[,1],popsample[,1])) # fst between Pop1 and Pop2, computed by ANGSD

# pi masked 95% CI
pi1_lo_masked <- matrix(NA, nrow=npop, ncol=npop, dimnames = list(popsample[,1],popsample[,1])) # 2.5% quantile of pi1 distribution
pi1_up_masked <- matrix(NA, nrow=npop, ncol=npop, dimnames = list(popsample[,1],popsample[,1])) # 97.5% quantile of pi1 distribution
pi2_lo_masked <- matrix(NA, nrow=npop, ncol=npop, dimnames = list(popsample[,1],popsample[,1])) # 2.5% quantile of pi2 distribution
pi2_up_masked <- matrix(NA, nrow=npop, ncol=npop, dimnames = list(popsample[,1],popsample[,1])) # 97.5% quantile of pi2 distribution

# dxy masked 95% CI
dxy_lo_masked <- matrix(NA, nrow=npop, ncol=npop, dimnames = list(popsample[,1],popsample[,1])) # 2.5% quantile of dxy distribution
dxy_up_masked <- matrix(NA, nrow=npop, ncol=npop, dimnames = list(popsample[,1],popsample[,1])) # 97.5% quantile of dxy distribution
# fst masked 95% CI
fst_lo_masked <- matrix(NA, nrow=npop, ncol=npop, dimnames = list(popsample[,1],popsample[,1])) # 2.5% quantile of fst distribution
fst_up_masked <- matrix(NA, nrow=npop, ncol=npop, dimnames = list(popsample[,1],popsample[,1])) # 97.5% quantile of fst distribution

# array to save the dxy of resampled SFS
dxy_resampled <- array(NA, dim=c(npop,npop,nboot))


# loop through each pairwise comparison
for(i in 957:ncol(index)) {
# for(i in 1:10) {
  pop1 <- popsample[index[1,i],1]
  pop2 <- popsample[index[2,i],1]

  # get the statistics from the 2D-SFS
  # get the file names
  filenametag1 <- paste(pop1, "_", pop2, sep="")
  filenametag2 <- paste(pop2, "_", pop1, sep="")
  # check which file exist as the name can be Pop1_Pop2.ml or Pop2_Pop1.ml
  if(file.exists(paste("./MLfiles/",filenametag1, ".ml",sep="")) & file.size(paste("./MLfiles/",filenametag1, ".ml",sep=""))>0) {
    filenametag <- filenametag1
    # indices of each pop - if the file name is Pop1_Pop2.ml, then pop1 is the first and pop2 is the second
    pop1_i <- index[1,i]
    pop2_i <- index[2,i]
  } else if(file.exists(paste("./MLfiles/",filenametag2, ".ml",sep="")) & file.size(paste("./MLfiles/",filenametag2, ".ml",sep=""))>0) {
    filenametag <- filenametag2
    # indices of each pop - if file name is Pop2_Pop1.ml, then pop1 is second and pop2 is the first
    pop1_i <- index[2,i]
    pop2_i <- index[1,i]
  }
  # print(filenametag)
  
  # Read the 2D SFS estimated with ANGSD as a vector
  sfs <- scan(paste("./MLfiles/",filenametag, ".ml",sep=""))
  # sfs
  
  # Read the FST computed by ANGSD
  fst_mean_angsd <- scan(paste("./FST/",filenametag, ".meanfst",sep=""))
  # fst_mean_angsd
  
  # get the sample size of each pop
  nind1 <- popsample[pop1_i,3]
  nind2 <- popsample[pop2_i,3]
  
  # convert to 2D-SFS matrix - rows correspond to pop1, but this depends on the filename
  # if(filenametag==filenametag1) {
  #   sfs2d <- matrix(sfs, nrow=(nind1*2)+1, byrow=T)
  # } else if(filenametag==filenametag2) {
  #   sfs2d <- matrix(sfs, nrow=(nind1*2)+1, byrow=T)
  # }
  # #sfs2d <- matrix(sfs, ncol=(nind1*2)+1, byrow=F)  
  sfs2d <- matrix(sfs, nrow=(nind1*2)+1, byrow=T)
  # print(sfs2d)
  
  # get the stats
  res <- get_pidxyfst_2dsfs_rel(sfs2d)
  # print(res)
  # mask the sites corresponding to all individuals treated as heterozygotes in both populations
  # as this might represent either missing data or mapping errors 
  sfs2d[nind1+1,nind2+1] <- 0
  res_masked <- get_pidxyfst_2dsfs_rel(sfs2d)
  
  # save the results of non-masked statistics
  pi1[pop1_i,pop2_i] <- res$pi[1]
  pi1[pop2_i,pop1_i] <- res$pi[2]
  
  pi2[pop2_i,pop1_i] <- res$pi[2]
  pi2[pop1_i,pop2_i] <- res$pi[1]
  
  dxy[pop1_i,pop2_i] <- res$dxy
  dxy[pop2_i,pop1_i] <- res$dxy
  
  fst[pop1_i,pop2_i] <- res$fst
  fst[pop2_i,pop1_i] <- res$fst
  
  # save the results of masked statistics
  pi1_masked[pop1_i,pop2_i] <- res_masked$pi[1]
  pi1_masked[pop2_i,pop1_i] <- res_masked$pi[2]
  
  pi2_masked[pop2_i,pop1_i] <- res_masked$pi[2]
  pi2_masked[pop1_i,pop2_i] <- res_masked$pi[1]
  
  dxy_masked[pop1_i,pop2_i] <- res_masked$dxy
  dxy_masked[pop2_i,pop1_i] <- res_masked$dxy
  
  fst_masked[pop1_i,pop2_i] <- res_masked$fst
  fst_masked[pop2_i,pop1_i] <- res_masked$fst
  
  # save the results of FST computed with ANGSD
  fst_angsd[pop1_i,pop2_i] <- fst_mean_angsd[2]
  fst_angsd[pop2_i,pop1_i] <- fst_mean_angsd[2]
  
  # Perform a bootstrap of sites from the 2D-SFS 
  # resampling sites to obtain CI for the different statistics
  # assuming that all sites are independent
  nsites <- sum(sfs)
  resampled_sfs <- rmultinom(n=nboot, size=nsites, prob=sfs)
  
  # compute the stats for each resampled 2D-SFS
  resampled_stats <- apply(resampled_sfs, 2, function(col) {
    # convert to 2D-SFS matrix
    # if(filenametag==filenametag1) {
    #   tmp2d <- matrix(col, nrow=(nind1*2)+1, byrow=T)
    # } else if(filenametag==filenametag2) {
    #   tmp2d <- matrix(col, nrow=(nind1*2)+1, byrow=T)
    # }
    tmp2d <- matrix(col, nrow=(nind1*2)+1, byrow=T)
    # mask entries 
    tmp2d[nind1+1,nind2+1] <- 0
    # get the stats
    res_tmp <- get_pidxyfst_2dsfs_rel(tmp2d)
    # output pi1, pi2, dxy, fst
    c(res_tmp$pi, res_tmp$dxy, res_tmp$fst)
  })
  
  # get the 0.025 and 0.975 quantiles for each stat, representing the 95% CI
  stat_quantiles <- apply(resampled_stats, 1, function(row) quantile(row,c(0.025,0.975)))
  
  # save the quantiles for the 95% CI
  pi1_lo_masked[pop1_i,pop2_i] <- stat_quantiles[1,1]
  pi1_lo_masked[pop2_i,pop1_i] <- stat_quantiles[1,2]
  
  pi1_up_masked[pop1_i,pop2_i] <- stat_quantiles[2,1]
  pi1_up_masked[pop2_i,pop1_i] <- stat_quantiles[2,2]
  
  pi2_lo_masked[pop2_i,pop1_i] <- stat_quantiles[1,2]
  pi2_lo_masked[pop1_i,pop2_i] <- stat_quantiles[1,1]
  
  pi2_up_masked[pop2_i,pop1_i] <- stat_quantiles[2,2]
  pi2_up_masked[pop1_i,pop2_i] <- stat_quantiles[2,1]
  
  dxy_lo_masked[pop1_i,pop2_i] <- stat_quantiles[1,3]
  dxy_up_masked[pop2_i,pop1_i] <- stat_quantiles[2,3]
  
  fst_lo_masked[pop1_i,pop2_i] <- stat_quantiles[1,4]
  fst_up_masked[pop2_i,pop1_i] <- stat_quantiles[2,4]
  
  # save the resampled dxy values into an array
  dxy_resampled[pop1_i,pop2_i,] <- resampled_stats[3,]
  dxy_resampled[pop2_i,pop1_i,] <- resampled_stats[3,]
  
}

# save the files
write.table(pi1, file="./PI1_allcomparisons.txt", quote=FALSE)
write.table(pi1_masked, file="./PI1_masked_allcomparisons.txt", quote=FALSE)

write.table(pi2, file="./PI2_allcomparisons.txt", quote=FALSE)
write.table(pi2_masked, file="./PI2_masked_allcomparisons.txt", quote=FALSE)

write.table(dxy, file="./DXY_allcomparisons.txt", quote=FALSE)
write.table(dxy_masked, file="./DXY_masked_allcomparisons.txt", quote=FALSE)

write.table(fst, file="./FST_allcomparisons.txt", quote=FALSE)
write.table(fst_masked, file="./FST_masked_allcomparisons.txt", quote=FALSE)

write.table(fst_angsd, file="./FST_ANGSD_allcomparisons.txt", quote=FALSE)

# save the 95% CI files
tmp_tosave <- list(pi1_lo_masked,pi1_up_masked,
                   pi2_lo_masked,pi2_up_masked,
                   dxy_lo_masked,dxy_up_masked,
                   fst_lo_masked,fst_up_masked)
save(tmp_tosave, file="ConfidenceIntervals.data")

# get the average PI per pop
comp1 <- rowMeans(pi1_masked, na.rm=TRUE)
comp2 <- rowMeans(pi2_masked, na.rm=TRUE)
pimean <- rowMeans(cbind(comp1,comp2), na.rm=TRUE)
write.table(pimean, file="./PI_masked_perpop_mean_allcomparisons.txt", quote=FALSE)

# write dxy into file in phylip format
writeDist(dxy_masked, file = "phylypdxy.phylip", format = "phylip")

# saved the dxy_resampled array with resampled matrices of dxy
save(dxy_resampled, file="dxy_resampled.data")

###################################
# NJ tree based on dxy ############
###################################

# read the dxy_masked and dxy_resampled
dxy_masked <- as.matrix(read.table("DXY_masked_allcomparisons.txt"))
str(dxy_masked)

# load the dxy_resampled array
load("dxy_resampled.data")
str(dxy_resampled)

# obtain NJ tree using fastME OLS
tree <- fastme.ols(dxy_masked)
plot(tree, type = "unrooted", show.tip = TRUE, show.node.label=FALSE)
# obtain the root using midpoint rooting
tree_rooted <- midpoint(tree)
plot(tree_rooted)

# write the tree in newick format
write.tree(tree_rooted, file="treefastMEOLS.newick")

# number of bootstrap replicates is the same as the dimension of array
nboot <- dim(dxy_resampled)[3]

# get a tree for each bootstrap replicate
# 1. get matrix of dxy based on the resampled dxy
# 2. get a NJ tree using fastme.ols method for each case
resampled_trees <- apply(dxy_resampled, 3, function(tmp) {
  # rename columns and rows of matrix, the same as in dxy_masked
  colnames(tmp) <- colnames(dxy_masked)
  rownames(tmp) <- rownames(dxy_masked)
  # get tree for each bootstrap matrix of dxy
  test <- fastme.ols(tmp)
  # root the tree using midpoint rooting
  test_rooted <- midpoint(test)})

# compute the proportion of times a given node appears in the bootstraped trees
clad <- prop.clades(tree_rooted, resampled_trees, rooted = TRUE)
clad
# compute the propotion of times a given edge appears in the bootstraped trees
boot <- prop.clades(tree_rooted, resampled_trees)
boot

# save file with tree drawing
pdf(file="NJ_tree.pdf", width=20, height=12)
plot(tree_rooted)
nodelabels(clad/nboot*100, cex=0.5)
# drawSupportOnEdges(boot/nboot*100)
dev.off()

# save tree in nexus format with bootstrap values
tree_rooted$node.label <- clad/nboot*100
write.tree(tree_rooted, file = "treeDxy_fastMEOLS.nwk")


###################################################################################################
# Look at the Pi per population with confidence intervals based on quantiles of resampling 2D-SFS #
###################################################################################################

load("ConfidenceIntervals.data")
str(tmp_tosave)

pi1_lo_masked <- tmp_tosave[[1]]
str(pi1_lo_masked)

pi1_up_masked <- tmp_tosave[[2]]
str(pi1_up_masked)

pi2_lo_masked <- tmp_tosave[[3]]
str(pi2_lo_masked)

pi2_up_masked <- tmp_tosave[[4]]
str(pi2_up_masked)

# read the mean pi per pop
pimean <- read.table("PI_masked_perpop_mean_allcomparisons.txt", col.names="pimean")

# get the minimum PI per pop for lower 95% CI 
lo1 <- apply(pi1_lo_masked, 1, function(row) min(row, na.rm=TRUE))
lo2 <- apply(pi2_lo_masked, 1, function(row) min(row, na.rm=TRUE))
pi_lo <- apply(cbind(lo1, lo2), 1, function(row) min(row, na.rm=TRUE))
# get the maximum PI per pop for upper 95% CI 
up1 <- apply(pi1_up_masked, 1, function(row) max(row, na.rm=TRUE))
up2 <- apply(pi2_up_masked, 1, function(row) max(row, na.rm=TRUE))
pi_up <- apply(cbind(up1, up2), 1, function(row) max(row, na.rm=TRUE))

# plot using ggplot
pipops <- data.frame(pop=1:length(pi_lo),pimean=pimean,pi_lo,pi_up)
ggplot(data=pipops, aes(x=pop, y=pimean)) +
  geom_point(aes(x=pop, y=pimean)) +
  geom_pointrange(aes(ymin=pi_lo, ymax=pi_up), color="#07bbc1")

# plot using plot functions
plot(pipops$pimean, ylim=c(0,0.0055))
segments(x0=1:length(pipops$pimean),y0=pipops$pi_lo, y1=pipops$pi_up)

# save file
write.table(tmp, file="PI_mean_95CI.txt", quote=F)


####################################################################
# GET ESTIMATE OF TMRCA FOR COMPARISONS WITHIN AND BETWEEN SPECIES #
####################################################################

# index of columns within species comparisons
caro <- c(1:15,17) # carolitertii
pyre <- c(16,18:31) # pyrenaicus
caet <- c(32:34) # caetobrigus
tart <- c(35:45) # tartessicus
torg <- c(46) # torgalensis
arad <- c(47:48) # aradensis
vale <- c(49) # valentinus
mala <- c(50:52) # malacitanus

# intra-specific comparisons
# remove the species with just one pop
index_pops_intra <- list(caro = c(1:15,17), # carolitertii
                   pyre = c(16,18:31), # pyrenaicus
                   caet = c(32:34), # caetobrigus
                   tart = c(35:45), # tartessicus
                   #torg = c(46), # torgalensis
                   arad = c(47:48), # aradensis
                   # vale = c(49), # valentinus
                   mala = c(50,52)) # malacitanus) # excluding JucarBV
                   # mala = c(50:52)) # malacitanus)

dxymat <- dxy_masked
# get the minimum- warning because of species just with 1 pop
min_intra <- sapply(index_pops_intra, function(x) min(dxymat[x,x], na.rm=TRUE))
# get the maximim - warning because of species just with 1 pop
max_intra <- sapply(index_pops_intra, function(x) max(dxymat[x,x], na.rm=TRUE))
# get the maximim - warning because of species just with 1 pop
mean_intra <- sapply(index_pops_intra, function(x) mean(dxymat[x,x], na.rm=TRUE))

# mean and range of dxy values of pairwise comparisons within species
dxyintra <- data.frame(pop=names(mean_intra), intra_mean_dxy=mean_intra, intra_min_dxy=min_intra, intra_max_dxy=max_intra)

# inter specific comparisons
index_pops <- list(caro = c(1:15,17), # carolitertii
                   pyre = c(16,18:31), # pyrenaicus
                   caet = c(32:34), # caetobrigus
                   tart = c(35:45), # tartessicus
                   torg = c(46), # torgalensis
                   arad = c(47:48), # aradensis
                   vale = c(49), # valentinus
                   mala = c(50,52)) # malacitanus) # excluding JucarBV

sp_pairs <- t(combn(1:length(index_pops),2)) # get the pairs of indices of pops to compare
min_inter <- apply(sp_pairs, 1, function(row) {min(dxymat[index_pops[[row[2]]],index_pops[[row[1]]]], na.rm=TRUE)})
max_inter <- apply(sp_pairs, 1, function(row) {max(dxymat[index_pops[[row[2]]],index_pops[[row[1]]]], na.rm=TRUE)})
mean_inter <- apply(sp_pairs, 1, function(row) {mean(dxymat[index_pops[[row[2]]],index_pops[[row[1]]]], na.rm=TRUE)})

# data frame with comparisons between pops of different species
pop.names <- c("caro","pyre","caet","tart","torg","arad","vale","mala")
dxyinter <- data.frame(pop1=pop.names[sp_pairs[,1]], pop2=pop.names[sp_pairs[,2]], mean_dxy=mean_inter, min_dxy=min_inter, max_dxy=max_inter)


# convert into estimates of TMRCA in years
gentime <- 3 # generation time in years
mutrate <- 5.97e-9 # mutation rate in mutations per site per generation
tmrca <- (dxyinter[,3:5]/(2*mutrate))
tmrca/1e6
colnames(tmrca) <- c("mean tmrca","min tmrca","max tmrca")

# tmrca between samples of two populations from different species in Million of generations ago
tmrca_inter_Mgen <- cbind(dxyinter[,1:2], round(tmrca/1e6, digits=2))

# tmrca between samples of two populations within the same species in Million of generations ago
tmrca_intra_Mgen <- cbind(dxyintra[,1], round((dxyintra[,2:4]/(2*mutrate))/1e6, digits=2))

# tmrca between two samples from the same population within a species in Millions of generations ago
tmrca_withinpop_Mgen <- cbind(pipops[,1], round((pipops[,2:4]/(2*mutrate))/1e6, digits=2))

# save plot with the distribution of pi, dxy and tmrca
pdf(file="Distribution_dxy_tmrca.pdf", width=8, height=6)
# create a data.frame to represent the data in ggplot
tmp <- data.frame(class=rep(c("withinPop","intraSpecies","interSpecies"), times=c(nrow(pipops),nrow(dxyintra),nrow(dxyinter))),
                  divergence=c(pipops$pimean,dxyintra$intra_mean_dxy,dxyinter$mean_dxy))
ggplot(data=tmp, aes(x=divergence, fill=class, y = after_stat(density))) +
  geom_histogram(alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#b2df8a","#a6cee3", "#1f78b4")) +
  labs(fill="")

# create a data.frame to represent TMRCA dates in ggplot
tmp <- data.frame(class=rep(c("withinPop","intraSpecies","interSpecies"), times=c(nrow(pipops),nrow(dxyintra),nrow(dxyinter))),
                  tmrca=c(tmrca_withinpop_Mgen$pimean,tmrca_intra_Mgen$intra_mean_dxy,tmrca_inter_Mgen$`mean tmrca`))
ggplot(data=tmp, aes(x=tmrca, fill=class, y = after_stat(density))) +
  geom_histogram(alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#b2df8a","#a6cee3", "#1f78b4")) +
  labs(fill="")
dev.of




# Compare the dxy with fst
# read the dxy_masked and dxy_resampled
fst_masked <- as.matrix(read.table("FST_masked_allcomparisons.txt"))
str(fst_masked)

# get the 95% CI limits for dxy anf fst
dxy_lo_masked <- tmp_tosave[[5]]
dxy_up_masked <- tmp_tosave[[6]]
fst_lo_masked <- tmp_tosave[[7]]
fst_up_masked <- tmp_tosave[[8]]

# fill the pairwise matrices of lo and up 
# a <- matrix(c(1,2,3,4,2,2,5,7,3,5,3,9,4,7,9,4), ncol=4, nrow=4, byrow=T)
# a[1,1] <- NA
# a[2,2] <- NA
# a[3,3] <- NA
# a[4,4] <- NA
# a[2,1] <- NA
# a[1,3] <- NA
# a[1,4] <- NA
# a[3,2] <- NA
# a[4,2] <- NA
# a[3,4] <- NA
# a
# 
# uppera <- a[upper.tri(a)]
# lowa <- t(a)[upper.tri(a)]
# uppera[is.na(uppera)] <- lowa[is.na(uppera)]

# fill the NA with the data from each lower or upper triangular matrix
dxy_lo <-  dxy_lo_masked[upper.tri(dxy_lo_masked)]
dxy_lo_lower <- t(dxy_lo_masked)[upper.tri(dxy_lo_masked)]
dxy_lo[is.na(dxy_lo)] <- dxy_lo_lower[is.na(dxy_lo)]
if(sum(is.na(dxy_lo))!=0) {
  stop("error: still missing data!")
}

dxy_up <-  dxy_up_masked[upper.tri(dxy_up_masked)]
dxy_up_lower <- t(dxy_up_masked)[upper.tri(dxy_up_masked)]
dxy_up[is.na(dxy_up)] <- dxy_up_lower[is.na(dxy_up)]
if(sum(is.na(dxy_up))!=0) {
  stop("error: still missing data!")
}

fst_lo <-  fst_lo_masked[upper.tri(fst_lo_masked)]
fst_lo_lower <- t(fst_lo_masked)[upper.tri(fst_lo_masked)]
fst_lo[is.na(fst_lo)] <- fst_lo_lower[is.na(fst_lo)]
if(sum(is.na(fst_lo))!=0) {
  stop("error: still missing data!")
}

fst_up <-  fst_up_masked[upper.tri(fst_up_masked)]
fst_up_lower <- t(fst_up_masked)[upper.tri(fst_up_masked)]
fst_up[is.na(fst_up)] <- fst_up_lower[is.na(fst_up)]
if(sum(is.na(fst_up))!=0) {
  stop("error: still missing data!")
}

# get a vector where each element corresponds to the population
pop_index <- numeric(ncol(dxy_masked))
for(i in 1:length(index_pops)) {
  pop_index[ index_pops[[i]] ] <- i
}

# create a matrix with the pop_index
index_mat_row <- matrix(pop_index, ncol=ncol(dxy_masked), nrow=ncol(dxy_masked), byrow=T)
index_mat_col <- matrix(pop_index, ncol=ncol(dxy_masked), nrow=ncol(dxy_masked), byrow=F)

# get the upper triangular matrix for each index of row and column
index_mat_row_up <- index_mat_row[upper.tri(index_mat_row)]
index_mat_col_up <- index_mat_col[upper.tri(index_mat_col)]

# get true or false for comparisons of the same species
intra_specific <- character() # initialize as a vector of strings
eval <- index_mat_row_up == index_mat_col_up
intra_specific[eval] <- "intra"
intra_specific[!eval] <- "inter"

# create data.frame to plot using ggplot
tmp <- data.frame(dxy=dxy_masked[upper.tri(dxy_masked)],
                  dxy_lo, dxy_up, 
                  fst=fst_masked[upper.tri(fst_masked)],
                  fst_lo, fst_up,
                  comparison=as.factor(intra_specific))

ggplot(data=tmp, aes(x=dxy, y=fst)) +
  geom_pointrange(aes(ymin=fst_lo, ymax=fst_up, colour=comparison), shape=21, fatten=.005, size=1) +
  geom_pointrange(aes(xmin=dxy_lo, xmax=dxy_up, colour=comparison), shape=21, fatten=.005, size=1) +
  scale_x_log10()
  
condition_i <- which((tmp$dxy<0.005 & tmp$fst<0.50) & tmp$comparison=="inter")

# Compare the da with fst
aux <- matrix(NA, ncol=ncol(dxy_masked), nrow=ncol(dxy_masked))
aux[upper.tri(aux)] <- 1:nrow(tmp)

row_col_i <- sapply(condition_i, function(x) which(aux==x, arr.ind=TRUE))

pops_condition <- apply(row_col_i, 2, function(col) colnames(dxy_masked)[col])

# Check da vs fst
pi1_masked <- as.matrix(read.table("PI1_masked_allcomparisons.txt"))
pi2_masked <- as.matrix(read.table("PI2_masked_allcomparisons.txt"))

da <- dxy_masked-((pi1_masked+pi2_masked)/2)
da1 <- dxy_masked-pi1_masked
da2 <- dxy_masked-pi1_masked
# min value
# replace negative values by a small da value
da[da<0] <- min(da[da>0], na.rm=T)/2

tmp <- data.frame(dxy=dxy_masked[upper.tri(dxy_masked)],
                  da=da[upper.tri(da)], 
                  fst=fst_masked[upper.tri(fst_masked)],
                  comparison=as.factor(intra_specific))

# proxy for migration rate up to a constant (mutation rate)
meanpi_tmp <- ((pi1_masked+pi2_masked)/2)
mutrate <- 5.97e-9
tmp$fst[tmp$fst<0] <- min(tmp$fst[tmp$fst>0])/10
m <- (mutrate*(1-tmp$fst))/(meanpi_tmp[upper.tri(meanpi_tmp)]*tmp$fst)
m1 <- (mutrate*(1-tmp$fst))/(pi1_masked[upper.tri(pi1_masked)]*tmp$fst)
m2 <- (mutrate*(1-tmp$fst))/(pi2_masked[upper.tri(pi2_masked)]*tmp$fst)

tmp <- data.frame(dxy=dxy_masked[upper.tri(dxy_masked)],
                  da=da[upper.tri(da)], 
                  da1=da1[upper.tri(da1)],
                  da2=da2[upper.tri(da2)],
                  fst=fst_masked[upper.tri(fst_masked)],
                  comparison=as.factor(intra_specific), 
                  m=m, 
                  m1=m1,
                  m2=m2)

pdf(file="dxy_da_fst.pdf", width=10, height=10)
ggplot(data=tmp, aes(x=dxy, y=fst)) +
  geom_point(aes(x=dxy, y=fst, colour = comparison)) +
  scale_x_log10()

ggplot(data=tmp, aes(x=dxy, y=1-fst)) +
  geom_point(aes(x=dxy, y=1-fst, colour = comparison)) +
  scale_x_log10()


# ggplot(data=tmp, aes(x=da, y=1-fst)) +
  # geom_point(aes(x=da, y=1-fst, colour = comparison)) +
  # scale_x_log10() + geom_vline(xintercept = c(0.005, 0.02))

ggplot(data=tmp, aes(x=da, y=(1-fst)/fst)) +
  geom_point(aes(x=da, y=(1-fst)/fst, colour = comparison)) +
  scale_x_log10() + scale_y_log10() #+geom_vline(xintercept = c(0.005, 0.02))

#+

ggplot(data=tmp, aes(x=da, y=m)) +
  geom_point(aes(x=da, y=m, colour = comparison)) +
  scale_x_log10() + scale_y_log10() #+ geom_vline(xintercept = c(0.005, 0.02))

ggplot(data=tmp, aes(x=da1, y=m1)) +
  geom_point(aes(x=da1, y=m1, colour = comparison)) +
  scale_x_log10() + scale_y_log10() #+ geom_vline(xintercept = c(0.005, 0.02))

ggplot(data=tmp, aes(x=da2, y=m2)) +
  geom_point(aes(x=da2, y=m2, colour = comparison)) +
  scale_x_log10() + scale_y_log10() #+ geom_vline(xintercept = c(0.005, 0.02))

dev.off()


# Get the TMRCA
mutrate <- 5.97e-9 # mutation rate in mutations per site per generation
tmp$tmrca <- tmp$dxy/(2*mutrate)

ggplot(data=tmp, aes(x=tmrca)) +
  geom_histogram(aes(x=tmrca, bg = comparison)) #+
  # scale_x_log10() + scale_y_log10() #+ geom_vline(xintercept = c(0.005, 0.02))
