# PLOTS OF THE MANUSCRIPT FOR THE RESULTS OF SELECTED SITE AND LINKED SELECTION

# clean memory and load required packages
rm(list=ls())

##########################################################################
# Supplementary Table S5
# summary mean ratio linked sites, scaled
##########################################################################

# read the data
sumstat <- read.table("SummaryFiles/fst_FitM2_h_r2.5e-7sumstat.txt", header=TRUE, na.strings ="NA")
sumstat_cb <-read.table("SummaryFiles/fst_FitM2_h_r2.5e-7sumstat_cb.txt", header=TRUE, na.strings ="NA")

# check that both files have the same order of combination of parameters
# check that the combination of parameters is the same for H and D
if(sum(sumstat[,c(1:5)]!=sumstat_cb[,c(1:5)])>0) {
  stop("Wrong order of combination of parameters!")
}

# combination of parameters
freq_v <- unique(sumstat$f)
dom_v <- unique(sumstat$d)

# Create a matrix with 6 columns and 4 rows for each case
fstsummary <- matrix(NA,ncol=6,nrow=length(freq_v)*3*2)
# fill the fstsummary for all sims
count <- 1
{
  # s>0, m>0
  eval <- sumstat$m>0 & sumstat$s>0
  for(i in 1:length(freq_v)) {
    
    # recessive
    # get index of correct combination of cases
    eval_i <- which(eval & sumstat$f==freq_v[i] & sumstat$d==0.01)
    # compute the ratio
    rec_ratio <- sumstat$HD_fst_mean[eval_i]/sumstat$D_fst_mean[eval_i]
    
    # codominant
    eval_i <- which(eval & sumstat$f==freq_v[i] & sumstat$d==0.5)
    # compute the ratio
    cod_ratio <- sumstat$HD_fst_mean[eval_i]/sumstat$D_fst_mean[eval_i]
    
    # get the mean, minimum and maximum
    fstsummary[count,] <- c(mean(rec_ratio), min(rec_ratio), max(rec_ratio),
                            mean(cod_ratio), min(cod_ratio), max(cod_ratio))
    
    count <- count + 1
  }
  # s>0, m=0
  eval <- sumstat$m==0 & sumstat$s>0
  for(i in 1:length(freq_v)) {
    
    # recessive
    # get index of correct combination of cases
    eval_i <- which(eval & sumstat$f==freq_v[i] & sumstat$d==0.01)
    # compute the ratio
    rec_ratio <- sumstat$HD_fst_mean[eval_i]/sumstat$D_fst_mean[eval_i]
    
    # codominant
    eval_i <- which(eval & sumstat$f==freq_v[i] & sumstat$d==0.5)
    # compute the ratio
    cod_ratio <- sumstat$HD_fst_mean[eval_i]/sumstat$D_fst_mean[eval_i]
    
    # get the mean, minimum and maximum
    fstsummary[count,] <- c(mean(rec_ratio), min(rec_ratio), max(rec_ratio),
                            mean(cod_ratio), min(cod_ratio), max(cod_ratio))
    
    count <- count + 1
  }
  # s=0
  eval <- sumstat$s==0
  for(i in 1:length(freq_v)) {
    
    # recessive
    # get index of correct combination of cases
    eval_i <- which(eval & sumstat$f==freq_v[i] & sumstat$d==0.01)
    # compute the ratio
    rec_ratio <- sumstat$HD_fst_mean[eval_i]/sumstat$D_fst_mean[eval_i]
    
    # codominant
    eval_i <- which(eval & sumstat$f==freq_v[i] & sumstat$d==0.5)
    # compute the ratio
    cod_ratio <- sumstat$HD_fst_mean[eval_i]/sumstat$D_fst_mean[eval_i]
    
    # get the mean, minimum and maximum
    fstsummary[count,] <- c(mean(rec_ratio), min(rec_ratio), max(rec_ratio),
                            mean(cod_ratio), min(cod_ratio), max(cod_ratio))
    
    count <- count + 1
  }
}
# fill for sims conditional on beneficial
{
  # s>0, m>0
  eval <- sumstat_cb$m>0 & sumstat_cb$s>0
  for(i in 1:length(freq_v)) {
    
    # recessive
    # get index of correct combination of cases
    eval_i <- which(eval & sumstat_cb$f==freq_v[i] & sumstat_cb$d==0.01)
    # compute the ratio
    rec_ratio <- sumstat_cb$HD_fst_mean[eval_i]/sumstat_cb$D_fst_mean[eval_i]
    
    # codominant
    eval_i <- which(eval & sumstat_cb$f==freq_v[i] & sumstat_cb$d==0.5)
    # compute the ratio
    cod_ratio <- sumstat_cb$HD_fst_mean[eval_i]/sumstat_cb$D_fst_mean[eval_i]
    
    # get the mean, minimum and maximum
    fstsummary[count,] <- c(mean(rec_ratio), min(rec_ratio), max(rec_ratio),
                            mean(cod_ratio), min(cod_ratio), max(cod_ratio))
    
    count <- count + 1
  }
  # s>0, m=0
  eval <- sumstat_cb$m==0 & sumstat_cb$s>0
  for(i in 1:length(freq_v)) {
    
    # recessive
    # get index of correct combination of cases
    eval_i <- which(eval & sumstat_cb$f==freq_v[i] & sumstat_cb$d==0.01)
    # compute the ratio
    rec_ratio <- sumstat_cb$HD_fst_mean[eval_i]/sumstat_cb$D_fst_mean[eval_i]
    
    # codominant
    eval_i <- which(eval & sumstat_cb$f==freq_v[i] & sumstat_cb$d==0.5)
    # compute the ratio
    cod_ratio <- sumstat_cb$HD_fst_mean[eval_i]/sumstat_cb$D_fst_mean[eval_i]
    
    # get the mean, minimum and maximum
    fstsummary[count,] <- c(mean(rec_ratio), min(rec_ratio), max(rec_ratio),
                            mean(cod_ratio), min(cod_ratio), max(cod_ratio))
    
    count <- count + 1
  }
  # s=0
  eval <- sumstat_cb$s==0
  for(i in 1:length(freq_v)) {
    
    # recessive
    # get index of correct combination of cases
    eval_i <- which(eval & sumstat_cb$f==freq_v[i] & sumstat_cb$d==0.01)
    # compute the ratio
    rec_ratio <- sumstat_cb$HD_fst_mean[eval_i]/sumstat_cb$D_fst_mean[eval_i]
    
    # codominant
    eval_i <- which(eval & sumstat_cb$f==freq_v[i] & sumstat_cb$d==0.5)
    # compute the ratio
    cod_ratio <- sumstat_cb$HD_fst_mean[eval_i]/sumstat_cb$D_fst_mean[eval_i]
    
    # get the mean, minimum and maximum
    fstsummary[count,] <- c(mean(rec_ratio), min(rec_ratio), max(rec_ratio),
                            mean(cod_ratio), min(cod_ratio), max(cod_ratio))
    
    count <- count + 1
  }
}
write.table(fstsummary, file="./Figures/SupTable3_linked_scaled.txt", quote=FALSE,
            col.names = c("mean_h0.01","min_h0.01","max_h0.01",
                          "mean_h0.50","min_h0.50","max_h0.50"), row.names = rep(freq_v, times=6))
