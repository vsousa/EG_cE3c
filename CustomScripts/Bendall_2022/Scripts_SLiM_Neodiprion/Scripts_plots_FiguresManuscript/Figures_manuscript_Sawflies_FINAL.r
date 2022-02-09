# Script to read the observed and simulated SFS
# according to the Neodiprion pinetum and N. lecontei demographic history
# Plot figures of manuscript.

# load required packages
library(svglite)
library(ggsci)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(gridExtra)
library(gridGraphics)
library(ggpubr)
library(ggtext)

# load functions
source("computeStats_SFS.r") # functions to compute stats and SFS
source("utilFscOutput.r") # functions to process fastsimcoal2 results

# create folder to save the figures
dir.create("Figures_manuscript_sawflies")

###################################################################
# Supplementary Figure S6 
# Compare the expected 2D-SFS obtained for fastsimcoal2 
# with the observed 2D-SFS of Neodipion sawflies
# and with the neutral SFS done with the SLiM simulations for
# diploid and haplodiploid case. If the scaling of parameters
# is correct, then the neutral SFS of SLiM should 
# be identical to the fastsimcoal2 expected SFS.
###################################################################

# read the observed SFS
# read observed sfs for intergenic regions
obs <- read.table("../SFS/excluded_1kb_Neo_L_P_2DSFS.joint2DMSFS", skip=1, header=TRUE, row.names = 1)
# remove the monomorphic sites, 
# setting the entry with fixed ancestral in both pops to zero, i.e. entry (0,0) 
# that in matrices of R corresponding to entry [1,1]
obs_nomon <- obs
obs_nomon[1,1] <- 0

# Obtain the summary statistics (pi, dxy, FST) 
# for the observed SFS for intergenic regions
obs_exc1kb <- get_pidxyfst_2dsfs(obs_nomon)

# Save the FST values of the observed sfs 
fstmeanobs <- obs_exc1kb$fst

# read expected 2D-SFS based on fastsimcoal2 IM model
expsfs <- read.table("../SFS/expectedSFS_mig_jointDAFpop1_0.txt", header=TRUE, row.names = 1)
# compute the statistics of the expected SFS of fastsimcoal2
expsfs_stats <- get_pidxyfst_2dsfs(expsfs)

# read the expected SFS for neutral case obtained with SLiM 
#slim_sfs <- readRDS(file="../SummaryFiles/sawfly_divsel3_t1549_scaled_r2_3_neutral_sfs.rds")
slim_sfs <- readRDS(file="../SummaryFiles/sawfly_divsel3_t1549_scaled_r2_3_1.05e-6_w20000_f0.1_neutral_sfs.rds")
# diploid neutral SFS (sum of SNPs across runs)
sfs2dsum_d <- slim_sfs$neutral_sfs_d
# haplodiploid neutral SFS (sum of SNPs across runs)
sfs2dsum_hd <- slim_sfs$neutral_sfs_hd

# get the stats for the SLiM neutral SFS for diploid and haplodiploid
stats_sfs_hd <- get_pidxyfst_2dsfs(sfs2dsum_hd)
stats_sfs_d <- get_pidxyfst_2dsfs(sfs2dsum_d)

# compare sumstats
# Compare the pairwise differences and dxy 
# since the SLiM neutral SFS is a mixture of either 2000 or 4000 sims (1000 for each init freq and dominance)
# and from recessive and co-dominant case
# the theta must be divided by 2000 (or 4000), and then multiplied by the inverse of the ratio
# ratio of expected neutral sfs to observed sfs in numbers of sites
nsites_exp <- 50000    # assumed number of sites when performing the scaling of params into SLiM 
                       # (check file DemographicHistory_ParamEstimates_SLIM_rescaling.xlsx, tab: ScaledByThetaNeodiprion)
nsites_obs <- sum(obs) # sum of sites in the observed 2D sfs of intergenic regions (including monomorphic sites)
r_eo <- nsites_exp/nsites_obs # ratio of expected/observed number of sites

# PI - pairwise differences
divideby <- 2000 # 4000 for cases where we analysed 2 initial freq and 2 dominance (1000 sims for each - see above)
pi_obs <- obs_exc1kb$pi[2:1]*r_eo # observed PI for intergenic regions
pi_fsc26 <- expsfs_stats$pi[2:1]*sum(obs_nomon)*r_eo # expected SFS PI based on fastsimcoal2
pi_hd <- (stats_sfs_hd$pi / divideby) # expected neutral PI for haplodiploid based on SLiM sims (need to change the order of pops)
pi_d <- (stats_sfs_d$pi / divideby) # expected neutral PI for diploid based on SLiM sims (need to change the order of pops)

# DXY - number of differences between pop1 and pop2
dxy_obs <- obs_exc1kb$dxy*r_eo # observed Dxyfor intergenic regions
dxy_fsc26 <- expsfs_stats$dxy*sum(obs_nomon)*r_eo # expected SFS Dxy based on fastsimcoal2
dxy_hd <- (stats_sfs_hd$dxy / divideby) # expected neutral Dxy for haplodiploid based on SLiM sims (need to change the order of pops)
dxy_d <- (stats_sfs_d$dxy / divideby) # expected neutral Dxy for diploid based on SLiM sims (need to change the order of pops)

# FST
# matrix with the FST weighted across the entire SFS
# and per site for the 4 conditions
# observed, fastsimcoal2, slim haplodiplod, slim diploid
fstmat <- t(matrix(c(obs_exc1kb$fst_site,obs_exc1kb$fst,
                     expsfs_stats$fst_site,expsfs_stats$fst,
                     stats_sfs_hd$fst_site,stats_sfs_hd$fst, 
                     stats_sfs_d$fst_site, stats_sfs_d$fst), ncol=4))

# MARGINAL 1D-SFS
# marginal 1D N. lecontei
marg_lec <- matrix(c(log10(rowSums(obs_nomon)), 
                     log10(rowSums(expsfs*sum(obs_nomon))),
                     log10(rowSums(t(sfs2dsum_hd/sum(sfs2dsum_hd))*sum(obs_nomon))),
                     log10(rowSums(t(sfs2dsum_d/sum(sfs2dsum_d))*sum(obs_nomon)))), ncol=4)
# marginal 1D N. pinetum
marg_pin <- matrix(c(log10(colSums(obs_nomon)), 
                     log10(colSums(expsfs*sum(obs_nomon))),
                     log10(colSums(t(sfs2dsum_hd/sum(sfs2dsum_hd))*sum(obs_nomon))),
                     log10(colSums(t(sfs2dsum_d/sum(sfs2dsum_d))*sum(obs_nomon)))), ncol=4)


# Set the figure layout with 8 panels corresponding to
# a) Fit of PI and DXY
# b) Fit of FST
# c) Fit of 1D marginal SFS of N. pinetum
# d) Fit of 1D marginal SFS of N. lecontei
# e) Observed 2D-SFS
# f) Expected 2D-SFS obtained with fastsimcoal2 under IM model
# g) Neutral 2D-SFS obtained with SLiM for haplodiploid 
# h) Neutral 2D-SFS obtained with SLiM for diploid
# NOTE: each 2D SFS plot actually contains the plot and the scale.
# need to account for that in the layout.
# pdf("./Figures_manuscript_sawflies/SupFig_S6.pdf", width = 18, height=7.2)
svglite("./Figures_manuscript_sawflies/SupFig_S6.svg", width = 18, height=7.2)
{
  layout(matrix(c(1,3,5:8,2,4,9:12), nrow=2, byrow=T), 
         widths = c(0.2, 0.25, 0.2, 0.05, 0.2, 0.05))
  par(mar=c(6,6,2,2))
  # a) barplot comparing PI and DXY
  # labels of barplot
  bar.names <-c(expression(paste(pi, " N. lec.")), expression(paste(pi, " N. pin.")),"Dxy")
  barplot(t(matrix(c(pi_obs,dxy_obs,pi_fsc26,dxy_fsc26,pi_hd,dxy_hd,pi_d,dxy_d), ncol=4)), beside=TRUE,
          names.arg = bar.names, 
          legend.text =c("obs.","exp. fsc","SLiM HD", "SLiM D"), 
          args.legend = list(x = "topleft",cex=1.5, bty="n"), xlab = "", 
          ylab="#pairwise diff. (50Kb)", main="",
          cex.axis = 1.5, cex.names = 2, cex.lab=1.5)  
  
  
  # b) barplot comparing FST
  bar.names <-c(expression(paste('F'['ST'], " weighted")), expression(paste('F'['ST'], " per site")))
  barplot(fstmat, beside=TRUE,
          names.arg = bar.names, 
          legend.text =c("obs.","exp. fsc","slim HD", "slim D"), 
          xlab = "", ylab=expression('F'['ST']), ylim=c(0,1), main="",
          args.legend = list(x = "topleft",cex=1.5, bty="n"),
          cex.axis = 1.5, cex.names = 2, cex.lab=1.5)  
  
  
  # c) marginal 1D SFS
  barplot(t(marg_lec),
          beside = TRUE, names.arg = 0:(nrow(obs_nomon)-1), 
          legend.text =c("obs","exp fsc","slim HD", "slim D"), 
          xlab = "N. lecontei", ylab="log10(#SNPs)", 
          main="", ylim=c(2.2,3.8), xpd = FALSE, args.legend = list(x = "top",cex=1.5, bty="n"),
          cex.axis = 1.5, cex.names = 2, cex.lab=1.5)
  barplot(t(marg_pin),
          beside = TRUE, names.arg = 0:(ncol(obs_nomon)-1), 
          legend.text =c("obs","exp fsc","slim HD", "slim D"), 
          xlab = "N. pinetum", ylab="log10(#SNPs)", main="", 
          ylim=c(2.2,3.8), xpd = FALSE, args.legend = list(x = "top",cex=1.5, bty="n"),
          cex.axis = 1.5, cex.names = 2, cex.lab=1.5)
  
  # Plot the 2D SFS and compare the SFS observed, fastsimcoal2, Slim HD and Slim D
  par(mar=c(5,5,3,1))
  plot2dSFS_single(obs_nomon, maintext = "Observed 2D-SFS", xtag = "N. lecontei", ytag = "N. pinetum", minentry = 3)
  plot2dSFS_single(expsfs*sum(obs_nomon), maintext = "Expected SFS fsc", xtag = "N. lecontei", ytag = "N. pinetum", minentry = 3)
  plot2dSFS_single(t(sfs2dsum_hd/sum(sfs2dsum_hd))*sum(obs_nomon), maintext = "SLiM SFS HD", xtag = "N. lecontei", ytag = "N. pinetum", minentry = 3)
  plot2dSFS_single(t(sfs2dsum_d/sum(sfs2dsum_d))*sum(obs_nomon), maintext = "SLiM SFS D", xtag = "N. lecontei", ytag = "N. pinetum", minentry = 3)
  dev.off()
}



#######################################################
# Figure 8 MAIN
# Compare the mean FST for simulations of haplodiploid
# and diploid according to sawfly demographic history
# with different selective coefficients, 
# initial frequencies and dominance
########################################################

# vector with the number of snps and corresponding window size
obs_window_snp <- c(19.02, 9.96, 5.93, 3.92) 
window_size <- c(100000,50000,25000,10000) 


# read the file with the sumstat results, including pi and fst
target_freq <- 0.1
sumstat <- read.table(paste("../SummaryFiles/sumstat_X_A_fst_targetwindow_f",target_freq,".txt",sep=""), header=TRUE, na.strings = "NA")
str(sumstat)

fst_genome <- read.table(paste("../SummaryFiles/fst_sawfly_divsel3_t1549_scaled_r2_3_1.05e-6_w20000_f",target_freq,"sumstat.txt",sep=""), header=TRUE, na.strings = "NA")
str(fst_genome)

# add columns of fst_genome to the sumstat data.frame
tmp_all <- data.frame(chr=c(rep("A", nrow(fst_genome)),rep("x", nrow(fst_genome))),
                 h=c(fst_genome$d,fst_genome$d),
                 m=c(fst_genome$m,fst_genome$m),
                 sr=rep(0.3,times=nrow(fst_genome)*2),
                 s=c(fst_genome$s,fst_genome$s),
                 r=c(fst_genome$r,fst_genome$r),
                 f=c(fst_genome$f,fst_genome$f),
                 t=rep(1549,times=nrow(fst_genome)*2),
                 window=rep(1e7, times=nrow(fst_genome)*2),
                 fst=c(fst_genome$D_fst_mean, fst_genome$HD_fst_mean),
                 fst0.05=c(fst_genome$D_fst_quantile0.05, fst_genome$HD_fst_quantile0.05),
                 fst0.25=c(fst_genome$D_fst_quantile0.25, fst_genome$HD_fst_quantile0.25),
                 fst0.5=c(fst_genome$D_fst_quantile0.50, fst_genome$HD_fst_quantile0.50),
                 fst0.75=c(fst_genome$D_fst_quantile0.75, fst_genome$HD_fst_quantile0.75),
                 fst0.95=c(fst_genome$D_fst_quantile0.95, fst_genome$HD_fst_quantile0.95))

# color scheme for H and D
mycols <- rev(brewer.pal(n=10,"RdBu"))[c(2,9)]
colvs <- mycols

# merge the data.frames just with fst values
sumstat_tmp <- sumstat[,c(1:9,13,15,17:19,21)]
str(sumstat_tmp)

# check that the names are the same
if(sum(colnames(sumstat_tmp) != colnames(tmp_all))>0) {
  stop("wrong column names in the two data.frames")
}

# merge the two data.frames just with info about fst
sumstat_fst <- rbind(tmp_all, sumstat_tmp)

# select only selection coefficients less than 0.3
row_i <- which(sumstat_fst$s <= 0.31)
sumstat_fst_row <- sumstat_fst[row_i,]

# list to save figures for different window sizes
figwindow <- list()
count <- 1

# vector with the number of snps and corresponding window size
# obs_window_snp <- c(19.02, 9.96, 5.93, 3.92) 
window_size <- c(50000,100000, 1e7) 
maxfst <- c(0.9,0.9,0.775)       # for 50kb, 100kb, 10Mb
minfst <- c(0.375,0.375,0.575)   # for 50kb, 100kb, 10Mb
pointsize <- 0.04
for(target_h in c(0.01,0.5)) {
  for(target_f in c(target_freq)) {
    # define the target window size
    for(i in 1:length(window_size)) {
      target.window <- window_size[i]
      # get the index corresponding to the target window size
      window_i <- which(target.window==sumstat_fst_row$window)
      
      # define a subset of results just for the selected sumstat 
      sumstat_sel <- sumstat_fst_row[window_i,]
      
      eval <- sumstat_sel$h==target_h & sumstat_sel$f==target_f
      sum(eval)
      sumstat_sel_hf <- sumstat_sel[eval,]
      
      # make the plot
      figwindow[[count]] <- ggplot(data=sumstat_sel_hf, aes(y=fst, x=s, fill=chr, group=chr, colour=chr)) +
        geom_line(size=1.2) + ylim(minfst[i],maxfst[i]) +
        geom_ribbon(aes(ymin=fst0.25, ymax=fst0.75, x=s, fill=chr, group=chr, colour=chr), alpha=0.1, size=0) +
        ggtitle(paste("h",target_h,"_f",target_f ,"_r",sumstat$r[1],"_w", target.window,sep="")) +
        theme_classic() +
        labs(y="Differentiation (F<sub>ST</sub>)", x = "Selection coefficient (s)") +
        theme(plot.title = element_text(size=8), 
              axis.text=element_text(size=10), 
              axis.title.x=element_markdown(size=12), 
              axis.title.y=element_markdown(size=12),
              legend.text = element_text(size=12))   +
        # geom_hline(yintercept=fstmeanobs[1], linetype = 'dashed') +
        geom_vline(xintercept = 0.156, linetype = 'dashed') +
        scale_color_manual(aesthetics = c("colour", "fill"), values=colvs) 
      count <- count + 1
    }
  }
}

# FIGURE OF GENOME SCAN
# positionben <- 250000/1000
# window.size <- 20000/1000
# slide.size <- 20000/1000
# seqlength <- 500000/1000
# seq_windows <- c(rev(seq(positionben-(window.size/2),0,by=-slide.size)),
#                  seq(positionben+(window.size/2),seqlength,by=slide.size))
# midpoint <- seq_windows[2:length(seq_windows)]-(slide.size/2)
positionben <- 5e6/1e6
window.size <- ((20000*1e7)/5e5)/1e6
slide.size <- ((20000*1e7)/5e5)/1e6
seqlength <- 1e7/1e6
seq_windows <- c(rev(seq(positionben-(window.size/2),0,by=-slide.size)),
                 seq(positionben+(window.size/2),seqlength,by=slide.size))
midpoint <- seq_windows[2:length(seq_windows)]-(slide.size/2)

# color scheme for H and D
mycols <- rev(brewer.pal(n=10,"RdBu"))[c(2,9)]
colvs <- mycols

# read the files with summary of window analyses
window.size <- 20
genscan <- readRDS(paste("../SummaryFiles/sawfly_divsel3_t1549_scaled_r2_3_1.05e-6_w",window.size*1000 ,"_f",target_freq,"_summary_window.rds",sep=""))
# get the combination of parameters for each case in window summary
tag <- genscan$tag
freq_scan <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("^f0.", a)],"f"))[2])}))
# get the results for FST
# the index 2 corresponds to results of sumstat_cb (1 for sumstat, 2 for sumstat_cb)
# the index 4 corresponds to fst (1 for pi pop1, 2 for pi pop2, 3 for dxy, 4 for fst)
fst_scan <- genscan$stats[[1]][[4]]
# fst_scan is a list of 5 matrices corresponding to the mean and quantiles 0.05, 0.25, 0.75 and 0.95
# Each matrix contains:
# - 25 rows corresponding to the 25 windows 20Kb 
# - columns corresponding to the parameter combinations in tag                                

# get the combination of parameters for the genome scan results
migrate_scan <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("m", a)],"m"))[2])}))
recrate_scan <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("^r", a)],"r"))[2])}))
dominance_scan <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("^h", a)],"h"))[2])}))
chrm_scan <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_"))[3]}))
selection_scan <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("^s0.", a)],"s"))[2])}))
freq_scan <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("^f0.", a)],"f"))[2])}))

# plot the genome scan along the genome corresponding to cases 1 and 2
case_sel <- unique(selection_scan)
#case_sel <- unique(selection_scan)[c(length(case_sel)-2,length(case_sel))]
case_sel <- unique(selection_scan)[c(length(case_sel)-2)]
# maxfst <- rep(0.775, times=length(case_sel)*2)
# minfst <- rep(0.65, times=length(case_sel)*2)
maxfst <- rep(0.815, times=length(case_sel)*2)
minfst <- rep(0.575, times=length(case_sel)*2)
# target_freq <- 0.1
plotscan <- list()
count <- 1
# loop through target cases
for(i in 1:length(case_sel)) {
  # loop through the dominance
  for(h in c(0.01, 0.5)) {
    # select the trajectory for D and H that correspond to the selected case
    eval_d <- which(freq_scan==target_freq & selection_scan==case_sel[i] & dominance_scan==h & chrm_scan=="A")
    eval_h <- which(freq_scan==target_freq & selection_scan==case_sel[i] & dominance_scan==h & chrm_scan=="X")
    
    # create a data.frame with the required info
    # diploid
    dfd <- data.frame(chrm=rep("D", nrow(fst_scan$mean)), 
                      pos= midpoint, 
                      fst_mean=fst_scan$mean[,eval_d], 
                      fst_25=fst_scan$q25[,eval_d],
                      fst_75=fst_scan$q75[,eval_d])
    # hemizygous 
    dfh <- data.frame(chrm=rep("H", nrow(fst_scan$mean)), 
                      pos=midpoint, 
                      fst_mean=fst_scan$mean[,eval_h], 
                      fst_25=fst_scan$q25[,eval_h],
                      fst_75=fst_scan$q75[,eval_h])
    # merge the two data.frames
    df <- rbind(dfd, dfh)

    # plot the FST along the chromosome
    plotscan[[count]] <- ggplot(data=df, aes(x=pos, y=fst_mean, group=chrm, colour=chrm)) + 
      geom_line(size=1.2)  + 
      # geom_ribbon(aes(ymin=fst_25, ymax=fst_75, group=chrm, fill=chrm), alpha=0.2) +
      geom_ribbon(aes(ymin=fst_25, ymax=fst_75, fill=chrm, group=chrm, colour=chrm), alpha=0.1, size=0) +
      theme_classic() +
      scale_color_manual(aesthetics = c("colour", "fill"), values=colvs) + 
      scale_y_continuous(limits = c(minfst[i], maxfst[i])) + 
      scale_x_continuous(limits = c(0, 10)) + 
      labs(y="Differentiation (F<sub>ST</sub>)", x = "Position (Mb)") +
      labs(colour = "") + 
      theme(axis.text=element_text(size=10), 
            axis.title.x=element_markdown(size=12), 
            axis.title.y=element_markdown(size=12),
            legend.text = element_text(size=12))   +
      # geom_hline(yintercept=fstmeanobs[1], linetype = 'dashed') +
      ggtitle(paste("s",case_sel[i],"_h",h,"_f",target_freq ,"_r",recrate_scan[1],"_w", slide.size,sep="")) +
      theme(plot.title = element_text(size=8))
    
    # add the expected value under neutrality
    # neutral <- getexpectedneutral_stats(N=1000, sr=0.5, U=2.5e-7, L=5e5, window.size=2e4, tsplit=2000, m=case_mig[i])$fst
    
    # plotscan[[count]] <- plotscan[[count]] + 
    #   geom_hline(yintercept=neutral, linetype="dotted", color = "black")
    # 
    # increase the count
    count <- count + 1
  } 
}

# Figure with only 100Kb and 10Mb
fig_v1 <- ggarrange(
  plotlist=c(figwindow,plotscan)[c(2:3,7,5:6,8)],
  nrow=2, ncol=3, common.legend = TRUE, legend="right", widths = c(0.3,0.3,0.4),
  labels = "auto")
fig_v1

# save into files
pdf(paste("./Figures_manuscript_sawflies/Fig8.pdf",sep=""), width=10, height=6)
fig_v1
dev.off()

svglite(paste("./Figures_manuscript_sawflies/Fig8.svg",sep=""), width=10, height=6)
fig_v1
dev.off()



