# PLOTS OF THE MANUSCRIPT FOR THE RESULTS OF SELECTED SITE AND LINKED SELECTION

# clean memory and load required packages
rm(list=ls())
library(doParallel)
library(ggplot2)
library(MBA)
library(lubridate)
library(reshape2)
library(colorRamps)
library(scales)
library(grid)
library(gridExtra)
library(gridGraphics)
library(ggpubr)
library(ggtext)
library(svglite)
library(RColorBrewer)
library(patchwork)
library(DescTools)
registerDoParallel(cores=2)

# load the functions to plot the heatmaps
source("./Rscripts/functions_plotHeatMap.r")
# load functions to compute expected value of stats under neutrality
source("../Scripts_process_SLiM_output/scaled/Rscripts/computeStats.r")

# COMBINATION OF PARAMETERS
# migration rate
migrate_v <-c("0.000", "0.00034", "0.0017", "0.0034")
# selection of mutation under divergent selection
selmutben_v <- c("0.000", "0.0067", "0.01333" ,"0.02667", "0.05333", "0.0667", "0.1334")
# recombination rate
recrate_v <- c("2.5e-7")
# initial frequency of allele a at locus under divergent selection
freqmut_v <- c("0.001", "0.01", "0.1", "0.5")
# dominance of allele a at locus under divergent selection
dominance_v <- c("0.01","0.5")
# time of split of the two populations
time_v <- c("2000")

# Get the number of parameter combinations
combnum <- length(migrate_v)*length(selmutben_v)*length(recrate_v)*length(dominance_v)*length(freqmut_v)*length(time_v)*2 


################################################
# SETTINGS AND READ DATA FOR FIGS 2-5
################################################

# SETTINGS
savetag <- "same_m"
# initial allele frequency 
target_freq <- 0.1
# cases combination of parameters for selection
case_sel <- as.numeric(selmutben_v[c(4,7)])
case_mig <- rep(0.0034, times=length(case_sel)) # for same_m
# number of generations to plot in trajectories for the cases
ngen <- c(2000, 300)
# for same_m
maxfst <- c(0.3,0.5) # scale used when only means ploted
# define the values for scaled by Ne  selection and migration in log10 scale
# effective size (Ne) that corresponds to a 2Ne=1500
# of a population of 1000 individuals at a hemizygous locus 
# (500 females with two copies and 500 males with one copy)
Ne <- 750

# all values larger than maxval replaced by maxval
maxval <- 1.5
# minimum and maximum values in plot
ylim_min <- 0
ylim_max <- maxval
# to avoid divisions by zero, replace 0 by small value
zeroval <- 1e-3
# title for plot
ytext <- "FST"
tagstat <- paste("H/D", ytext)
# define a relative error, 
# such that all values between 1*(1-relerror) and 1*relerror are set to 1
relerror <- 0.05

# filenames with results from simulations
filename_site <- "./SummaryFiles/Summary_Traj_t2000_r2.5e-7.txt"
filename_traj <- "./SummaryFiles/Traj_t2000_r2.5e-7_summary_traj.rds"
filename_linked <- "./SummaryFiles/fst_FitM2_h_r2.5e-7sumstat.txt"
filename_linked_cb <- "./SummaryFiles/fst_FitM2_h_r2.5e-7sumstat_cb.txt"
filename_scan <- "./SummaryFiles/FitM2_h_r2.5e-7_f0.1_summary_window.rds"

# index of window with selected site
window_selsite <- 13

# color scheme for H and D
mycols <- rev(brewer.pal(n=10,"RdBu"))[c(2,9)]

# color for trajectory for the D and H
colvs <- mycols
# END SETTINGS

######################################
# READ AND PROCESS RESULTS
######################################

  # selected site
  # read the summary table for the results at the selected site
  site_sumstat <- read.table(filename_site, header=TRUE, na.strings="NA")
  str(site_sumstat)
  # read the allele trajectory files for selected site
  traj <- readRDS(filename_traj)
  # keep only the trajectories for the conditional on keeping allele a
  traj <- traj$summary_traj_cb
  # assume that the order of combination of parameters is the same as in site_sumstat
  # check that at least the number of elements is the same
  if(length(traj)!=nrow(site_sumstat)) {
    stop("The site_sumstat and traj have different number of parameter combinations.")
  }

  #################
  # linked sites
  sumstat <- read.table(file=paste(filename_linked,sep=""), header=TRUE, na.strings = "NA")
  str(sumstat)
  summary(sumstat)
  sumstat_cb <- read.table(file=paste(filename_linked_cb,sep=""), header=TRUE, na.strings = "NA")
  str(sumstat_cb)
  summary(sumstat_cb)
  # check that sumstat and sumstat_cb have the same combination of parameters
  if(sum(sumstat[,1:5]!=sumstat_cb[,1:5])>0) {
    stop("The sumstat and sumstat_cb have different parameter combinations.")
  }
  
  # read the files with summary of window analyses
  genscan <- readRDS(filename_scan)
  # get the combination of parameters for each case in window summary
  tag <- genscan$tag
  freq_scan <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("^f0.", a)],"f"))[2])}))
  # get the results for FST
  # the index 2 corresponds to results of sumstat_cb (1 for sumstat, 2 for sumstat_cb)
  # the index 4 corresponds to fst (1 for pi pop1, 2 for pi pop2, 3 for dxy, 4 for fst)
  fst_scan <- genscan$stats[[2]][[4]]
  # fst_scan is a list of 5 matrices corresponding to the mean and quantiles 0.05, 0.25, 0.75 and 0.95
  # Each matrix contains:
  # - 25 rows corresponding to the 25 windows 20Kb 
  # - columns corresponding to the parameter combinations in tag      
  
  # FST at window of selected site across all sims (without condition on runs that kept beneficial allele)
  selectedsite_fst_scan <- genscan$stats[[1]][[4]] # all sims
  selectedsite_fst_scan_cb <- genscan$stats[[2]][[4]] # conditional (_cb)
  
  # keep only the sims with freq=target_freq
  sumstat <- sumstat[sumstat$f==target_freq,]
  sumstat_cb <- sumstat_cb[sumstat_cb$f==target_freq,]
  traj <- traj[which(site_sumstat$f==target_freq)]
  site_sumstat <- site_sumstat[site_sumstat$f==target_freq,]
  for(i in 1:length(fst_scan)) {
    fst_scan[[i]] <- fst_scan[[i]][,freq_scan==target_freq]  
    selectedsite_fst_scan[[i]] <- selectedsite_fst_scan[[i]][,freq_scan==target_freq]
    selectedsite_fst_scan_cb[[i]] <- selectedsite_fst_scan_cb[[i]][,freq_scan==target_freq]
  }
  
  # check that the combination of parameters is the same for selected site and linked sites
  deval <- which(site_sumstat$chrm=="D")
  heval <- which(site_sumstat$chrm=="H")
  if(sum(site_sumstat[deval,c(1:3,5:6)]!=sumstat[,1:5])>0) {
    stop("The selsite and linked sites files do not have the same combination of parameters.")
  }
  
  # get the combination of parameters for each case
  migrate <- sumstat$m
  recrate <- sumstat$r
  dominance <- sumstat$d
  selection <- sumstat$s
  freq <- sumstat$f
  # get the combination of parameters for the genome scan results
  migrate_scan <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("m", a)],"m"))[2])}))
  recrate_scan <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("^r", a)],"r"))[2])}))
  dominance_scan <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("^h", a)],"h"))[2])}))
  chrm_scan <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_"))[2]}))
  selection_scan <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("^s0.", a)],"s"))[2])}))
  
  # select the trajectory for D and H that correspond to the selected case
  eval_d_scan <- which(chrm_scan=="A")
  eval_h_scan <- which(chrm_scan=="X")
  
  # check that the combination of parameters are the same
  tmp <- data.frame(m=migrate_scan[eval_d_scan], r=recrate_scan[eval_d_scan],d=dominance_scan[eval_d_scan],s=selection_scan[eval_d_scan])
  if(sum(tmp != sumstat[,1:4])>0) {
    stop("The genome scan and linked sites files do not have the same combination of parameters.")
  }
  
  # Get only results for the window with selected site, 
  # to avoid division of small values average across all windows for neutral case when m=0.0034
  for(i in 1:length(selectedsite_fst_scan)) {
    # for the neutral case with high migration
    selneutral <- selection_scan==0 & migrate_scan==0.0034
    # average across all windows, which can be done as all windows are neutral
    # this reduces the stochastic variation
    neutral_highmig <- colMeans(selectedsite_fst_scan[[i]][,selneutral])
    neutral_highmig_cb <- colMeans(selectedsite_fst_scan_cb[[i]][,selneutral])
    # get only the values at the window with selected site    
    selectedsite_fst_scan[[i]] <- selectedsite_fst_scan[[i]][window_selsite,]
    selectedsite_fst_scan_cb[[i]] <- selectedsite_fst_scan_cb[[i]][window_selsite,]
    # replace the neutral cases
    selectedsite_fst_scan[[i]][selneutral] <- neutral_highmig
    selectedsite_fst_scan_cb[[i]][selneutral] <- neutral_highmig_cb
  }
  
  # define the selection and migration rate
  # transformation to selection values
  selplot <- round(2*Ne*selection)
  # transformation to migration values
  migplot <- round(2*Ne*migrate, digits=1)
  
  # indicate if ratio of stats is ploted in logscale
  logscale <- F
  # indicate if ratio of stats have reversed colours
  revcol <- F
  
  #############################################################
  #############################################################
  # FIGURE 2
  #############################################################
  #############################################################
  
  # define the breaks of the scale
  breaks <- c(0.5,0.75,0.95,1.05,1.25,1.50)
  
  # ALL SIMS - NON CONDITIONAL
  # Plot the heatmaps of H/D statistic
  # at window with selected site
  selsite_window <- plotheat_rel_nointerpol(dstat=selectedsite_fst_scan$mean[eval_d_scan],
                                    hdstat=selectedsite_fst_scan$mean[eval_h_scan],
                                    tagstat=tagstat,
                                    zeroval=zeroval,
                                    breaks=breaks,
                                    dominance_v=as.numeric(dominance_v),
                                    freqmut_v=target_freq,
                                    dom=dominance,
                                    freq=freq,
                                    sel=selplot,
                                    mig=migplot,
                                    revcol = revcol,
                                    log10scale=logscale,
                                    relerror=relerror,
                                    xlabel="Selection (2*Ns*)", ylabel="Migration (2*Nm*)")
  # linked sites
  linked <- plotheat_rel_nointerpol(dstat=sumstat$D_fst_mean,
                                   hdstat=sumstat$HD_fst_mean,
                                   tagstat=tagstat,
                                   zeroval=zeroval,
                                   breaks=breaks,
                                   dominance_v=as.numeric(dominance_v),
                                   freqmut_v=target_freq,
                                   dom=dominance,
                                   freq=freq,
                                   sel=selplot,
                                   mig=migplot,
                                   revcol = revcol,
                                   log10scale=logscale,
                                   relerror=relerror,
                                   xlabel="Selection (2*Ns*)", ylabel="Migration (2*Nm*)")
  
  # ALL SIMS
  # window with selected site
  fstmax_plot <- c(0.6,0.6)
  pointsize <- 0.4
  fst_selwindow_points <- plot_fst_points_v1(mstat=selectedsite_fst_scan$mean, 
                                             q25stat=selectedsite_fst_scan$q25, 
                                             q75stat=selectedsite_fst_scan$q75,
                                             tagstat="", maxval=fstmax_plot[1], minval=0, 
                                             dominance_v=as.numeric(dominance_v), 
                                             dom=dominance_scan, mig=migrate_scan, sel=selection_scan, 
                                             chrm=chrm_scan, case_mig=case_mig[1], colsvs=colsvs, pointsize = pointsize)
  
  # Linked sites 
  fst_linked_points <- plot_fst_points_v2(mstatA=sumstat$D_fst_mean, 
                                          q25statA=sumstat$D_fst_quantile0.25, 
                                          q75statA=sumstat$D_fst_quantile0.75, 
                                          mstatX=sumstat$HD_fst_mean, 
                                          q25statX=sumstat$HD_fst_quantile0.25, 
                                          q75statX=sumstat$HD_fst_quantile0.75, 
                                          tagstat="", maxval=fstmax_plot[2], minval=0, 
                                          dominance_v=as.numeric(dominance_v), 
                                          dom=dominance, mig=migrate, sel=selection, 
                                          case_mig=case_mig[1], colsvs=colsvs, pointsize = pointsize)
  # Figure 2 with points 
  fig2_points <- ggarrange(
    plotlist=c(fst_selwindow_points, selsite_window$figs, fst_linked_points, linked$figs)[c(1,3,5,7,2,4,6,8)],
    nrow=2, ncol=4, common.legend = TRUE, legend="right", labels = "auto", vjust=0.5) +
    theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm")) 
    
  fig2_points
  
  
  # save figures in PDF and SVG formats
  # PDF
  pdf(file=paste("./Figures/Fig2_main.pdf",sep=""), width=10, height=4.8)
  fig2_points
  dev.off()
  
  # SVG 
  svglite(file=paste("./Figures/Fig2_main.svg",sep=""), width=10, height=4.8)
  fig2_points
  dev.off()
  
  
  #############################################################
  #############################################################
  # FIG. 3 Probability of keeping beneficial allele
  # PROBABILITY OF KEEPING THE ALLELE
  #############################################################
  #############################################################
  
  # plot the probability of keeping the alleles for the two extreme migration rates
  # Probability of keep
  nsims <- 1000
  pkeep <- site_sumstat$pkeep # probability of retaining derived allele
  nkeep <- pkeep*nsims # number of sims that kept the derived allele
  
  # get the confidence interval assuming binomial CI
  tmp <- BinomCI(x=nkeep, n=nsims, conf.level = 0.95, method = c("clopper-pearson"), rand = 123)
  low_ci <- tmp[,2]
  upper_ci <- tmp[,3]
  rm(tmp)
  
  # auxiliary variable to define the line type
  aux <- c(1,1,1,2)
  pointsize <- 0.2
    
  # create data.frame to use in ggplot
  # using only the no mig and mig
  selplot_site <- round(2*Ne*site_sumstat$s)
  migplot_site <- round(2*Ne*site_sumstat$m, digits=1) 
  
  pkeep_fig <- list()
  inset_fig <- list()
  
  # Loop through dominance coefficients
  for(h in 1:length(dominance_v)) {
    # get the runs that have the selected freq
    eval_dom_mig <- site_sumstat$h==as.numeric(dominance_v[h]) & (site_sumstat$m==as.numeric(migrate_v[1]) | site_sumstat$m==as.numeric(migrate_v[4]))
    df_keep <- data.frame(chr=site_sumstat$chrm[eval_dom_mig], s=selplot_site[eval_dom_mig], m=as.factor(migplot_site[eval_dom_mig]), 
                          pkeep=pkeep[eval_dom_mig], lci=low_ci[eval_dom_mig], uci=upper_ci[eval_dom_mig])
    
    pkeep_fig[[h]] <- ggplot(data=df_keep, aes(x=s, y=pkeep, group=interaction(m, chr), color=chr)) + 
      geom_line(aes(linetype=m),size=1.1, show.legend = FALSE) + 
      # geom_point(aes(shape=chr), size=pointsize, show.legend = FALSE) +
      geom_pointrange(aes(ymin=lci, ymax=uci, group=interaction(m, chr), fill=chr), position=position_dodge(0.1), show.legend = FALSE, size=pointsize)  +
      theme_classic() +
      xlab("Selection (2*Ns*)") + ylab("Probability of retention") +
      # scale_x_log10() +
      scale_y_continuous(limits = c(0, 1)) + 
      scale_linetype_manual(values=c("solid", "dashed"))+
      scale_color_manual(values=mycols) + 
      theme(axis.title.x =element_markdown(size=16),
            axis.text.x = element_text(size = 14),
            axis.title.y = element_markdown(size=16),
            axis.text.y = element_text(size = 14)) 
    
    inset_fig[[h]] <- pkeep_fig[[h]]  + 
      scale_x_continuous(limits = c(0, 41)) + xlab("") + ylab("") + 
      theme_minimal() +
      theme(axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10)) 
  }
  
  # without inset
  fig3_pkeep <- ggarrange(
    plotlist=pkeep_fig, 
    nrow=1, ncol=2, common.legend = TRUE, legend="right")
  # with inset
  pkeep_fig[[2]] <- pkeep_fig[[2]] + inset_element(inset_fig[[2]], left = 0.25, bottom = -0.05, right = 1.05, top = 0.8)
  fig3_pkeep_inset <- ggarrange(
    plotlist=pkeep_fig, 
    nrow=1, ncol=2, common.legend = TRUE, legend="right")
  
  # Save figure 3 files
  pdf(file=paste("./Figures/Fig3_inset.pdf",sep=""), width=9, height=4)
  fig3_pkeep_inset
  dev.off()
  
  svglite(file=paste("./Figures/Fig3_inset.svg",sep=""), width=9, height=4)
  fig3_pkeep_inset
  dev.off()
  
  
  #############################################################
  #############################################################
  # FIGURE 4 and FIGURE 5
  # CONDITIONAL on keeping beneficial allele
  # used for Fig. 4 and Fig. 5
  #############################################################
  #############################################################
  
  ##############
  # HEATMAPS
  ##############

  # Plot the heatmaps of H/D statistic
  # at window with selected site
  selsite_window <- plotheat_rel_nointerpol(dstat=selectedsite_fst_scan_cb$mean[eval_d_scan],
                                            hdstat=selectedsite_fst_scan_cb$mean[eval_h_scan],
                                            tagstat=tagstat,
                                            zeroval=zeroval,
                                            breaks=breaks,
                                            dominance_v=as.numeric(dominance_v),
                                            freqmut_v=target_freq,
                                            dom=dominance,
                                            freq=freq,
                                            sel=selplot,
                                            mig=migplot,
                                            revcol = revcol,
                                            log10scale=logscale,
                                            relerror=relerror,
                                            xlabel="Selection (2*Ns*)", ylabel="Migration (2*Nm*)")
  
  # linked sites conditional on keeping allele
  linked <- plotheat_rel_nointerpol(dstat=sumstat_cb$D_fst_mean,
                                    hdstat=sumstat_cb$HD_fst_mean,
                                    tagstat=tagstat,
                                    zeroval=zeroval,
                                    breaks=breaks,
                                    dominance_v=as.numeric(dominance_v),
                                    freqmut_v=target_freq,
                                    dom=dominance,
                                    freq=freq,
                                    sel=selplot,
                                    mig=migplot,
                                    revcol = revcol,
                                    log10scale=logscale,
                                    relerror=relerror,
                                    xlabel="Selection (2*Ns*)", ylabel="Migration (2*Nm*)")
  
  #################################################
  # ALLELE TRAJECTORY AND FST SCAN
  #################################################
  
  # plot the allele trajectories for the cases selected
  plottraj <- list()
  count <- 1
  # loop through target cases
  for(i in 1:length(case_sel)) {
    # loop through the dominance
    for(h in c(0.01, 0.5)) {
      # select the trajectory for D and H that correspond to the selected case
      eval_i <- which(migrate==case_mig[i] & selection==case_sel[i] & dominance==h)
      # create a data.frame with the required info
      # diploid
      tmpd <- traj[deval[eval_i]][[1]]
      dfd <- data.frame(chrm=rep("D", nrow(tmpd)+1), 
                        gen=c(0,tmpd[,colnames(tmpd)=="gens"]), 
                        p1=c(target_freq, tmpd[,colnames(tmpd)=="tmp_mean_p1"]),
                        p1_25=c(target_freq, tmpd[,which(colnames(tmpd)=="tmp_mean_p1")+2]),
                        p1_75=c(target_freq, tmpd[,which(colnames(tmpd)=="tmp_mean_p1")+4]),
                        p2=c(target_freq, tmpd[,colnames(tmpd)=="tmp_mean_p2"]),
                        p1_25=c(target_freq, tmpd[,which(colnames(tmpd)=="tmp_mean_p2")+2]),
                        p1_75=c(target_freq, tmpd[,which(colnames(tmpd)=="tmp_mean_p2")+4]))
      # hemizygous 
      tmph <- traj[heval[eval_i]][[1]]
      dfh <- data.frame(chrm=rep("H", nrow(tmph)+1), 
                        gen=c(0,tmph[,colnames(tmph)=="gens"]), 
                        p1=c(target_freq, tmph[,colnames(tmph)=="tmp_mean_p1"]),
                        p1_25=c(target_freq, tmph[,which(colnames(tmph)=="tmp_mean_p1")+2]),
                        p1_75=c(target_freq, tmph[,which(colnames(tmph)=="tmp_mean_p1")+4]),
                        p2=c(target_freq, tmph[,colnames(tmph)=="tmp_mean_p2"]),
                        p1_25=c(target_freq, tmph[,which(colnames(tmph)=="tmp_mean_p2")+2]),
                        p1_75=c(target_freq, tmph[,which(colnames(tmph)=="tmp_mean_p2")+4]))
      # merge the two data.frames
      df <- rbind(dfd, dfh)
      
      # plot the trajectories
      plottraj[[count]] <- ggplot(data=df, aes(x=gen, y=p1, group=chrm, colour=chrm)) + 
        geom_line(size=1)  + 
        # geom_ribbon(aes(ymin=p1_25, ymax=p1_75, group=chrm, fill=chrm), alpha=0.2) +
        theme_classic() +
        scale_color_manual(values=colvs) + 
        scale_y_continuous(limits = c(0, 1)) + 
        scale_x_continuous(limits = c(0, ngen[i])) + 
        # geom_line(data=df, aes(x=gen, y=p2, group=chrm, colour=chrm), size=2, linetype="dashed") +
        geom_line(data=df, aes(x=gen, y=p2, group=chrm, colour=chrm), size=1) +
        labs(y="Frequency of allele *a* (*q*)", x = "Time (generations)") +
        labs(colour = "") + 
        theme(axis.text=element_text(size=10), axis.title=element_text(size=12), 
              axis.title.y=element_markdown(size=12),
              legend.text = element_text(size=12))   
      
      # increase the count
      count <- count + 1
    }
  }
  
  # plot the genome scan along the genome corresponding to cases 1 and 2
  plotscan <- list()
  count <- 1
  # loop through target cases
  for(i in 1:length(case_sel)) {
    # loop through the dominance
    for(h in c(0.01, 0.5)) {
      # select the trajectory for D and H that correspond to the selected case
      eval_d <- which(migrate_scan==case_mig[i] & selection_scan==case_sel[i] & dominance_scan==h & chrm_scan=="A")
      eval_h <- which(migrate_scan==case_mig[i] & selection_scan==case_sel[i] & dominance_scan==h & chrm_scan=="X")
      
      # create a data.frame with the required info
      # diploid
      dfd <- data.frame(chrm=rep("D", nrow(fst_scan$mean)), 
                        pos=seq(10,500,by=20), 
                        fst_mean=fst_scan$mean[,eval_d], 
                        fst_25=fst_scan$q25[,eval_d],
                        fst_75=fst_scan$q75[,eval_d])
      # hemizygous 
      dfh <- data.frame(chrm=rep("H", nrow(fst_scan$mean)), 
                        pos=seq(10,500,by=20), 
                        fst_mean=fst_scan$mean[,eval_h], 
                        fst_25=fst_scan$q25[,eval_h],
                        fst_75=fst_scan$q75[,eval_h])
      # merge the two data.frames
      df <- rbind(dfd, dfh)
      
      # plot the FST along the chromosome
      plotscan[[count]] <- ggplot(data=df, aes(x=pos, y=fst_mean, group=chrm, colour=chrm)) + 
        geom_line(size=1)  + 
        # geom_ribbon(aes(ymin=fst_25, ymax=fst_75, group=chrm, fill=chrm), alpha=0.2) +
        theme_classic() +
        scale_color_manual(values=colvs) + 
        scale_y_continuous(limits = c(0, maxfst[i])) + 
        scale_x_continuous(limits = c(0, 500)) + 
        labs(y="Differentiation (F<sub>ST</sub>)", x = "Position (Kb)") +
        labs(colour = "") + 
        theme(axis.text=element_text(size=10), axis.title.y=element_markdown(size=12), 
              axis.title.x=element_markdown(size=12), legend.text = element_text(size=12))   
      
      # add the expected value under neutrality
      neutral <- getexpectedneutral_stats(N=1000, sr=0.5, U=2.5e-7, L=5e5, window.size=2e4, tsplit=2000, m=case_mig[i])$fst
      
      plotscan[[count]] <- plotscan[[count]] + 
        geom_hline(yintercept=neutral, linetype="dotted", color = "black")
      
      # increase the count
      count <- count + 1
    } 
  }
 
 
  ###############################################################
  # FST 
  ###############################################################
 
  fstmax_plot <- c(0.6,0.3)
  pointsize <- 0.4
  fst_selwindow_points <- plot_fst_points_v1(mstat=selectedsite_fst_scan_cb$mean, 
                                             q25stat=selectedsite_fst_scan_cb$q25, 
                                             q75stat=selectedsite_fst_scan_cb$q75,
                                             tagstat="", maxval=fstmax_plot[1], minval=0, 
                                             dominance_v=as.numeric(dominance_v), 
                                             dom=dominance_scan, mig=migrate_scan, sel=selection_scan, 
                                             chrm=chrm_scan, case_mig=case_mig[1], colsvs=colsvs, pointsize)

  # Linked sites 
  fst_linked_points <- plot_fst_points_v2(mstatA=sumstat_cb$D_fst_mean, 
                                          q25statA=sumstat_cb$D_fst_quantile0.25, 
                                          q75statA=sumstat_cb$D_fst_quantile0.75, 
                                          mstatX=sumstat_cb$HD_fst_mean, 
                                          q25statX=sumstat_cb$HD_fst_quantile0.25, 
                                          q75statX=sumstat_cb$HD_fst_quantile0.75, 
                                          tagstat="", maxval=fstmax_plot[2], minval=0, 
                                          dominance_v=as.numeric(dominance_v), 
                                          dom=dominance, mig=migrate, sel=selection, 
                                          case_mig=case_mig[1], colsvs=colsvs, pointsize)
  
  ###################################
  # Fig 4. RECESSIVE
  ###################################
  
  # change the outer margins of figures to increase white space between them
  for(i in 1:length(selsite_window$figs)) {
    selsite_window$figs[[i]] <-  selsite_window$figs[[i]] + theme(plot.margin=unit(c(1,2,4,1),"lines")) 
    linked$figs[[i]] <- linked$figs[[i]] + theme(plot.margin=unit(c(1,2,4,1),"lines")) 
    fst_selwindow_points[[i]] <- fst_selwindow_points[[i]] + theme(plot.margin=unit(c(1,2,2,1),"lines")) 
    fst_linked_points[[i]] <- fst_linked_points[[i]] + theme(plot.margin=unit(c(1,2,2,1),"lines")) 
  }
  
  for(i in 1:length(plottraj)) {
    plottraj[[i]] <-  plottraj[[i]] + theme(plot.margin=unit(c(1,4,2,1),"lines")) 
    plotscan[[i]] <-  plotscan[[i]] + theme(plot.margin=unit(c(1,2,2,1),"lines")) 
  }
  
  
  # with heatmap, scan and trajectory
  index <- c(1,3)
  # with heatmap, fst points, scan and trajectory
  figlist <- c(selsite_window$figs[1], fst_selwindow_points[1], plottraj[index],
               linked$figs[1],fst_linked_points[1], plotscan[index])
  fig4_v3 <- ggarrange(
    plotlist=figlist,
    nrow=2, ncol=4, common.legend = FALSE, legend="none", labels = "auto")
  
  # save fig4 as PDF and as SVG
  pdf(file=paste("./Figures/Fig4_main_m", case_mig[1],".pdf",sep=""), width=12.6, height=5.4)
  fig4_v3
  dev.off()
  
  svglite(file=paste("./Figures/Fig4_main_m", case_mig[1],".svg",sep=""), width=12.6, height=5.4)
  fig4_v3
  dev.off()
  
  
  ###################################
  # Fig 5. CODOMINANT
  ###################################
  
  # with heatmap, scan and trajectory
  index <- c(2,4)
  figlist <- c(selsite_window$figs[2], plottraj[index],
               linked$figs[2], plotscan[index])

  # with heatmap, fst points, scan and trajectory
  figlist <- c(selsite_window$figs[2], fst_selwindow_points[2], plottraj[index],
               linked$figs[2],fst_linked_points[2], plotscan[index])
  fig5_v3 <- ggarrange(
    plotlist=figlist,
    nrow=2, ncol=4, common.legend = FALSE, legend="none", labels = "auto")

  # save fig5 as PDF and as svg
  pdf(file=paste("./Figures/Fig5_main_m", case_mig[1],".pdf",sep=""), width=12.6, height=5.4)
  fig5_v3
  dev.off()
  
  svglite(file=paste("./Figures/Fig5_main_m", case_mig[1],".svg",sep=""), width=12.6, height=5.4)
  fig5_v3
  dev.off()
  
##################################################################
##################################################################
# SUPPLEMENTARY FIGURES 
##################################################################
##################################################################

  ##########################################################################
  # Sup Fig. 1 
  # check that the neutral statistics for the neutral case are correct
  ##########################################################################
  
  # read the statistics and the quantiles for pi, dxy, fst
  # FST - all simulations - not conditioning on retaining the allele a
  scaled_fst <- read.table(file=paste("./SummaryFiles/fst_FitM2_h_r2.5e-7sumstat.txt",sep=""), header=TRUE, na.strings = "NA")
  str(scaled_fst)
  summary(scaled_fst)
  
  # PI - all simulations - not conditioning on retaining the allele a
  scaled_pi <- read.table(file=paste("./SummaryFiles/pi1_FitM2_h_r2.5e-7sumstat.txt",sep=""), header=TRUE, na.strings = "NA")
  str(scaled_pi)
  summary(scaled_pi)
  
  # DXY - all simulations - not conditioning on retaining the allele a
  scaled_dxy <- read.table(file=paste("./SummaryFiles/dxy_FitM2_h_r2.5e-7sumstat.txt",sep=""), header=TRUE, na.strings = "NA")
  str(scaled_dxy)
  summary(scaled_dxy)
  
  # LD patterns - r^2
  summarymeanr2 <- read.table(file="./SummaryFiles/neutral_meanr2.txt", header=TRUE)
  str(summarymeanr2)
  summary(summarymeanr2)
  
  # get the expected values for the statistics
  expected <- list()
  for(i in 1:length(migrate_v)) {
    expected[[i]] <- unlist(getexpectedneutral_stats(N=1000, sr=0.5, U=2.5e-7, L=5e5, window.size=2e4, tsplit=2e3, m=as.numeric(migrate_v[i]))  )
  }
  expected <- do.call(rbind, expected)
  
  # FST
  {
    # for through migration rates
    df <- list()
    for(i in 1:length(migrate_v)) {
      target_m <- as.numeric(migrate_v[i])
      eval <- scaled_fst$m==target_m & scaled_fst$s==0
      tmp <- scaled_fst[eval,]
      
      # create a data.frame with mean, quantile0.25, quantile0.75
      df[[i]] <- data.frame(medianstat=c(mean(tmp$D_fst_quantile0.50),mean(tmp$HD_fst_quantile0.50)),
                            q005=c(mean(tmp$D_fst_quantile0.05),mean(tmp$HD_fst_quantile0.05)),
                            q025=c(mean(tmp$D_fst_quantile0.25),mean(tmp$HD_fst_quantile0.25)),
                            q075=c(mean(tmp$D_fst_quantile0.75),mean(tmp$HD_fst_quantile0.75)),
                            q095=c(mean(tmp$D_fst_quantile0.95),mean(tmp$HD_fst_quantile0.95)),
                            group=rep(c("D","H"), each=1), mig=rep(target_m, times=2))
    }
    # merge list into a single data frame
    df <- do.call(rbind, df)
    # define colors
    colvs <- c(gray(0.5),1)
    # make boxplot
    fst <- ggplot(df, aes(x=factor(mig), colour=group)) +
      geom_boxplot(aes(x=factor(mig), colour=group, ymin=q005, lower=q025, middle=medianstat, upper=q075, ymax=q095), stat="identity") + 
      xlab("Migration rate") +
      ylab(expression(F[ST])) +
      scale_color_manual(values=colvs) +
      # scale_color_grey() +
      theme_classic() +
      geom_hline(yintercept = expected[,5], linetype="dashed", color = gray(0.6))
  }
  
  # PI
  {
    # for through migration rates
    df <- list()
    for(i in 1:length(migrate_v)) {
      target_m <- as.numeric(migrate_v[i])
      eval <- scaled_pi$m==target_m & scaled_pi$s==0
      tmp <- scaled_pi[eval,]
      
      # create a data.frame with mean, quantile0.25, quantile0.75
      df[[i]] <- data.frame(medianstat=c(mean(tmp$D_pi1_quantile0.50),mean(tmp$HD_pi1_quantile0.50)),
                            q005=c(mean(tmp$D_pi1_quantile0.05),mean(tmp$HD_pi1_quantile0.05)),
                            q025=c(mean(tmp$D_pi1_quantile0.25),mean(tmp$HD_pi1_quantile0.25)),
                            q075=c(mean(tmp$D_pi1_quantile0.75),mean(tmp$HD_pi1_quantile0.75)),
                            q095=c(mean(tmp$D_pi1_quantile0.95),mean(tmp$HD_pi1_quantile0.95)),
                            group=rep(c("D","H"), each=1), mig=rep(target_m, times=2))
    }
    # merge list into a single data frame
    df <- do.call(rbind, df)
    # define colors
    colvs <- c(gray(0.5),1)
    # make boxplot
    pi <- ggplot(df, aes(x=factor(mig), colour=group)) +
      geom_boxplot(aes(x=factor(mig), colour=group, ymin=q005, lower=q025, middle=medianstat, upper=q075, ymax=q095), stat="identity") + 
      xlab("Migration rate") +
      ylab("Pairwise diff. P1") +
      scale_color_manual(values=colvs) +
      # scale_color_grey() +
      theme_classic() +
      geom_hline(yintercept = expected[,1], linetype="dashed", color = gray(0.6))
  }
  
  # DXY
  {
    # for through migration rates
    df <- list()
    for(i in 1:length(migrate_v)) {
      target_m <- as.numeric(migrate_v[i])
      eval <- scaled_dxy$m==target_m & scaled_dxy$s==0
      tmp <- scaled_dxy[eval,]
      
      # create a data.frame with mean, quantile0.25, quantile0.75
      df[[i]] <- data.frame(medianstat=c(mean(tmp$D_dxy_quantile0.50),mean(tmp$HD_dxy_quantile0.50)),
                            q005=c(mean(tmp$D_dxy_quantile0.05),mean(tmp$HD_dxy_quantile0.05)),
                            q025=c(mean(tmp$D_dxy_quantile0.25),mean(tmp$HD_dxy_quantile0.25)),
                            q075=c(mean(tmp$D_dxy_quantile0.75),mean(tmp$HD_dxy_quantile0.75)),
                            q095=c(mean(tmp$D_dxy_quantile0.95),mean(tmp$HD_dxy_quantile0.95)),
                            group=rep(c("D","H"), each=1), mig=rep(target_m, times=2))
    }
    # merge list into a single data frame
    df <- do.call(rbind, df)
    # define colors
    colvs <- c(gray(0.5),1)
    # make boxplot
    dxy <- ggplot(df, aes(x=factor(mig), colour=group)) +
      geom_boxplot(aes(x=factor(mig), colour=group, ymin=q005, lower=q025, middle=medianstat, upper=q075, ymax=q095), stat="identity") + 
      xlab("Migration rate") +
      ylab("Dxy") +
      scale_color_manual(values=colvs) +
      # scale_color_grey() +
      theme_classic() +
      geom_hline(yintercept = expected[,3], linetype="dashed", color = gray(0.6))
  }
  
  # R2
  {
    
    # merge list into a single data frame
    df <- summarymeanr2
    # define colors
    colvs <- c(gray(0.5),1)
    # make boxplot
    r2 <- ggplot(df, aes(x=factor(m), colour=chrm)) +
      geom_boxplot(aes(x=factor(m), colour=chrm, ymin=q05, lower=q25, middle=q50, upper=q75, ymax=q95), stat="identity") + 
      xlab("Migration rate") +
      ylab(expression(LD(r^2))) +
      scale_color_manual(values=colvs) +
      # scale_color_grey() +
      theme_classic() 
  }
  
  # Arrange the plots and plot them to file
  supfig1 <- ggarrange(plotlist=list(pi, dxy, fst, r2), 
                       nrow=2, ncol=2,  labels = c("a", "b", "c", "d"),
                       common.legend = TRUE, legend="right")
  # save SVG file that can be opened with Inkscape
  svglite(file=paste("./Figures/SupFigS1.svg",sep=""), width=7, height=7)
  supfig1
  dev.off()
  # save pdf
  pdf(file=paste("./Figures/SupFigS1.pdf",sep=""), width=7, height=7)
  supfig1
  dev.off()
  
  
##########################################################################
# Sup. Fig. S2  and
# Sup. Fig. S4
##########################################################################

# SETTINGS 
# filename with results for 500-kb
filename_linked <- "./SummaryFiles/fst_FitM2_h_r2.5e-7"
# define the values for scaled by Ne  selection and migration in log10 scale
# effective size (Ne) that corresponds to a 2Ne=1500
# of a population of 1000 individuals at a hemizygous locus 
# (500 females with two copies and 500 males with one copy)
Ne <- 750
nsims <- 1000 # number of simulations
threshold_nsims <- 10 # minimum number of sims kept to be considered
# all values larger than maxval replaced by maxval
maxval <- 1.5
# minimum and maximum values in plot
ylim_min <- 0
ylim_max <- maxval
# to avoid divisions by zero, replace 0 by small value
zeroval <- 1e-3
# define a relative error, 
# such that all values between 1*(1-relerror) and 1*relerror are set to 1
relerror <- 0.05
# to avoid numeric issues due to ratios of small values
# replace values smaller than a given threshold by the threshold value
min_fst_threshold <- 0.02
# indicate if ratio of stats is plotted in logscale
logscale <- F
# indicate if ratio of stats have reversed colors
revcol <- F
# define the breaks of the scale
breaks <- c(0.5,0.75,0.95,1.05,1.25,1.50)
filename_scan <- "./SummaryFiles/FitM2_h_r2.5e-7_summary_window.rds"
filename_site <- "./SummaryFiles/Summary_Traj_t2000_r2.5e-7.txt"
window_selsite <- 13 # index of window with selected site
freqmut_v <- c(0.001, 0.01, 0.1, 0.5) # target frequencies
# END OF SETTINGS

# selected site
# read the summary table for the results at the selected site
site_sumstat <- read.table(filename_site, header=TRUE, na.strings="NA")
str(site_sumstat)

# all simulations - not conditioning on retaining the allele a
sumstat <- read.table(file=paste(filename_linked, "sumstat.txt",sep=""), header=TRUE, na.strings = "NA")
str(sumstat)
summary(sumstat)

# conditional on retaining the allele a
sumstat_cb <- read.table(file=paste(filename_linked,"sumstat_cb.txt",sep=""), header=TRUE, na.strings = "NA")
str(sumstat_cb)
summary(sumstat_cb)

# check that the combination of parameters is the same for selected site and linked sites
deval <- which(site_sumstat$chrm=="D")
heval <- which(site_sumstat$chrm=="H")
if(sum(site_sumstat[deval,c(1:3,5:6)]!=sumstat[,1:5])>0) {
  stop("The selsite and linked sites files do not have the same combination of parameters.")
}


# RATIO of FST at window with selectes site
  # read the files with summary of window analyses
  genscan <- readRDS(filename_scan)
  # get the combination of parameters for each case in window summary
  tag <- genscan$tag
  # freq_scan <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("^f0.", a)],"f"))[2])}))
  # get the results for FST
  # the index 2 corresponds to results of sumstat_cb (1 for sumstat, 2 for sumstat_cb)
  # the index 4 corresponds to fst (1 for pi pop1, 2 for pi pop2, 3 for dxy, 4 for fst)
  # fst_scan <- genscan$stats[[1]][[4]]
  # fst_scan is a list of 5 matrices corresponding to the mean and quantiles 0.05, 0.25, 0.75 and 0.95
  # Each matrix contains:
  # - 25 rows corresponding to the 25 windows 20Kb 
  # - columns corresponding to the parameter combinations in tag      
  
  # FST at window of selected site across all sims (without condition on runs that kept beneficial allele)
  selectedsite_fst_scan <- genscan$stats[[1]][[4]] # all sims
  selectedsite_fst_scan_cb <- genscan$stats[[2]][[4]] # conditional (_cb)
 
  # get the combination of parameters for each case
  migrate <- sumstat$m
  recrate <- sumstat$r
  dominance <- sumstat$d
  selection <- sumstat$s
  freq <- sumstat$f
  # get the combination of parameters for the genome scan results
  migrate_scan <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("m", a)],"m"))[2])}))
  recrate_scan <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("^r", a)],"r"))[2])}))
  dominance_scan <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("^h", a)],"h"))[2])}))
  chrm_scan <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_"))[2]}))
  selection_scan <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("^s0.", a)],"s"))[2])}))
  freq_scan <- unlist(lapply(tag, function(x) {a <- unlist(strsplit(x, "_")); as.numeric(unlist(strsplit(a[grep("^f0.", a)],"f"))[2])}))
  # select the trajectory for D and H that correspond to the selected case
  eval_d_scan <- which(chrm_scan=="A")
  eval_h_scan <- which(chrm_scan=="X")
  
  # check that the combination of parameters are the same
  tmp <- data.frame(m=migrate_scan[eval_d_scan], r=recrate_scan[eval_d_scan],d=dominance_scan[eval_d_scan],s=selection_scan[eval_d_scan],f=freq_scan[eval_d_scan])
  if(sum(tmp[tmp$s>0,] != sumstat[sumstat$s>0,1:5])>0) {
    stop("The genome scan and linked sites files do not have the same combination of parameters.")
  }
  
  # Get only results for the window with selected site, 
  # to avoid division of small values average across all windows for neutral case when m=0.0034
  for(i in 1:length(selectedsite_fst_scan)) {
    # for the neutral case with high migration
    selneutral <- selection_scan==0 & migrate_scan==0.0034
    # average across all windows, which can be done as all windows are neutral
    # this reduces the stochastic variation
    neutral_highmig <- colMeans(selectedsite_fst_scan[[i]][,selneutral])
    neutral_highmig_cb <- colMeans(selectedsite_fst_scan_cb[[i]][,selneutral])
    # get only the values at the window with selected site    
    selectedsite_fst_scan[[i]] <- selectedsite_fst_scan[[i]][window_selsite,]
    selectedsite_fst_scan_cb[[i]] <- selectedsite_fst_scan_cb[[i]][window_selsite,]
    # replace the neutral cases
    selectedsite_fst_scan[[i]][selneutral] <- neutral_highmig
    selectedsite_fst_scan_cb[[i]][selneutral] <- neutral_highmig_cb
  }
  
  # get the number of runs per case
  nkeep <- site_sumstat$pkeep*nsims
  str(nkeep) # size of nkeep is twice the number of combinations
  # replace the cases with less than threshold_nsims by NA
  # get the nkeep for D and for H
  d_nkeep <- nkeep[site_sumstat$chrm=="D"]
  h_nkeep <- nkeep[site_sumstat$chrm=="H"]
  selectedsite_fst_scan_cb$mean[nkeep <= threshold_nsims] <- NA
  sumstat_cb$D_fst_mean[which(d_nkeep <= threshold_nsims)] <- NA
  sumstat_cb$HD_fst_mean[which(h_nkeep <= threshold_nsims)] <- NA
  
  # transformation to selection values
  selplot <- as.numeric(format(2*Ne*selection, digits = 0))
  # transformation to migration values
  migplot <- as.numeric(format(2*Ne*migrate, digits = 1))
  
  # title for plot
  ytext <- "FST"
  tagstat <- paste("H/D", ytext)
  
  # all simulations
  # at window with selected site
  selsite_window <- plotheat_rel_nointerpol(dstat=selectedsite_fst_scan$mean[eval_d_scan],
                                            hdstat=selectedsite_fst_scan$mean[eval_h_scan],
                                            tagstat=tagstat,
                                            zeroval=zeroval,
                                            breaks=breaks,
                                            dominance_v=as.numeric(dominance_v),
                                            freqmut_v=freqmut_v,
                                            dom=dominance,
                                            freq=freq,
                                            sel=selplot,
                                            mig=migplot,
                                            revcol = revcol,
                                            log10scale=logscale,
                                            relerror=relerror,
                                            xlabel="2*Ns*", ylabel="2*Nm*")
  
  # conditional on allele retention
  selsite_window_cb <- plotheat_rel_nointerpol(dstat=selectedsite_fst_scan_cb$mean[eval_d_scan],
                                            hdstat=selectedsite_fst_scan_cb$mean[eval_h_scan],
                                            tagstat=tagstat,
                                            zeroval=zeroval,
                                            breaks=breaks,
                                            dominance_v=as.numeric(dominance_v),
                                            freqmut_v=freqmut_v,
                                            dom=dominance,
                                            freq=freq,
                                            sel=selplot,
                                            mig=migplot,
                                            revcol = revcol,
                                            log10scale=logscale,
                                            relerror=relerror,
                                            xlabel="2*Ns*", ylabel="2*Nm*")
  
  # Write the Supplementary Table with the Summary of FST values for windowd with selected site
  tmp <- data.frame(h=dominance_scan[eval_d_scan],
                    f=freq_scan[eval_d_scan],
                    m=migrate_scan[eval_d_scan], 
                    s=selection_scan[eval_d_scan],
                    fstall_d=selectedsite_fst_scan$mean[eval_d_scan],
                    fstall_h=selectedsite_fst_scan$mean[eval_h_scan],
                    fstcb_d=selectedsite_fst_scan_cb$mean[eval_d_scan],
                    fstcb_h=selectedsite_fst_scan_cb$mean[eval_h_scan])
  # write.table(tmp, file="./Figures/Summary_selsite_window_meanFST.txt", row.names = F, quote = FALSE)


# RATIO of FST at linked sites in 500-kb window

    # check that the combination of parameters is the same in both
    # sumstat and sumstat_cb, i.e that the first 5 columns are the same
    if(sum(sumstat[,1:5]!=sumstat_cb[,1:5])>0) {
      stop("The two files do not have the same combination of parameters.")
    }
    
    # Plot the heatmaps of H/D statistic
    # all simulations
    relstat_sumstat_linked <- plotheat_rel_nointerpol(dstat=sumstat$D_fst_mean,
                                              hdstat=sumstat$HD_fst_mean,
                                              tagstat=tagstat,
                                              zeroval=zeroval,
                                              breaks=breaks,
                                              dominance_v=as.numeric(dominance_v),
                                              freqmut_v=as.numeric(freqmut_v),
                                              dom=dominance,
                                              freq=freq,
                                              sel=selplot,
                                              mig=migplot,
                                              revcol = revcol,
                                              log10scale=logscale,
                                              relerror=relerror,
                                              xlabel="2*Ns*", ylabel="2*Nm*")
    # conditional on allele retention
    relstat_sumstat_linked_cb <- plotheat_rel_nointerpol(dstat=sumstat_cb$D_fst_mean,
                                                 hdstat=sumstat_cb$HD_fst_mean,
                                                 tagstat=tagstat,
                                                 zeroval=zeroval,
                                                 breaks=breaks,
                                                 dominance_v=as.numeric(dominance_v),
                                                 freqmut_v=as.numeric(freqmut_v),
                                                 dom=dominance,
                                                 freq=freq,
                                                 sel=selplot,
                                                 mig=migplot,
                                                 revcol = revcol,
                                                 log10scale=logscale,
                                                 relerror=relerror,
                                                 xlabel="2*Ns*", ylabel="2*Nm*")
    
# Arrange the plots and plot them to file
# ratio at selected site, all sims and conditional
# for different initial frequencies

# figure with 4 columns corresponding to q=1/2N, 0.01, 0.1, 0.5
# and 4 rows corresponding to 20kb-rec, 500kb-rec, 20kb-cod, 500kb-cod
supfig2 <- ggarrange(plotlist=c(selsite_window$figs[c(1:4)], 
                                   relstat_sumstat_linked$figs[c(1:4)],
                                   selsite_window$figs[c(5:8)], 
                                   relstat_sumstat_linked$figs[c(5:8)]), 
                        nrow=4, ncol=length(freqmut_v), 
                        common.legend = TRUE, legend="right")
supfig2



# Arrange the plots and plot them to file
# ratio at selected site, all sims and conditional
# for different initial frequencies
# figure with 4 columns corresponding to q=1/2N, 0.01, 0.1, 0.5
# and 4 rows corresponding to 20kb-rec, 500kb-rec, 20kb-cod, 500kb-cod
supfig4 <- ggarrange(plotlist=c(selsite_window_cb$figs[c(1:4)], 
                                   relstat_sumstat_linked_cb$figs[c(1:4)],
                                   selsite_window_cb$figs[c(5:8)], 
                                   relstat_sumstat_linked_cb$figs[c(5:8)]), 
                        nrow=4, ncol=length(freqmut_v), 
                        common.legend = TRUE, legend="right")
supfig4


# save SVG file that can be opened with Inkscape
svglite(file=paste("./Figures/SupFigS2.svg",sep=""), width=12, height=8)
supfig2
dev.off()

# save PDF
pdf(file=paste("./Figures/SupFigS2.pdf",sep=""), width=12, height=8)
supfig2
dev.off()

# save SVG file that can be opened with Inkscape
svglite(file=paste("./Figures/SupFigS4.svg",sep=""), width=12, height=8)
supfig4
dev.off()

# save PDF
pdf(file=paste("./Figures/SupFigS4.pdf",sep=""), width=12, height=8)
supfig4
dev.off()



##########################################################################################
# Supplementary Fig. S3 
# probability of keeping the allele 
# at two extreme migration rates
##########################################################################################

# SETTINGS
# define the colors for D and H
colvs <- c(gray(0.7), 1)
filename <- "./SummaryFiles/Summary_Traj_t2000_r2.5e-7.txt"
# define the values for scaled by Ne  selection and migration in log10 scale
# effective size (Ne) that corresponds to a 2Ne=1500
# of a population of 1000 individuals at a hemizygous locus 
# (500 females with two copies and 500 males with one copy)
Ne <- 750
# To avoid the log10(0) issue to represent neutral case 
# replace the zero migration and zero selection arbitrary values
nonzero_sel <- 10^0.5
# frequency values in figure
freqmut_v <- freqmut_v <- c("0.001","0.01","0.1", "0.5")

# END OF SETTINGS

pdfFile=FALSE # TRUE for PDF, FALSE for SVG
{
  
  
  # read the summary table for the results at the selected site
  site_sumstat <- read.table(filename, header=TRUE, na.strings="NA")
  str(site_sumstat)
  
  # Get the combination of parameters for each case
  migrate <- site_sumstat$m
  chrm <- site_sumstat$chrm
  dominance <- site_sumstat$h
  selection <- site_sumstat$s
  freq <- site_sumstat$f
  
  # transformation to selection values
  selplot <- (2*Ne*selection)
  # transformation to migration values
  selplot[selplot==0] <- nonzero_sel
  # transform using log10 scale
  selplot <- log10(selplot)
  unique(selplot)
  
  # Probability of keep
  pkeep <- site_sumstat$pkeep
  
  # labels of chromosomes
  dip <- unique(chrm)
  
  # auxiliary variable to define the line type
  aux <- c(1,1,1,2)
  
  # save SVG file that can be opened with Inkscape
  if(pdfFile) {
    pdf(file=paste("./Figures/SupFigS3.pdf",sep=""), width=10, height=5)
  } else {
    svglite(file=paste("./Figures/SupFigS3.svg",sep=""), width=10, height=5)  
  }
  
  # figure with 2 rows and 4 columns
  par(mfrow=c(length(dominance_v),length(freqmut_v)), mar=c(5,5,1,1))
  
  # Loop through dominance coefficients
  for(h in 1:length(dominance_v)) {
    # loop through initial frequency considered (1 plot for each initial frequency)
    for(f in 1:length(freqmut_v)) {
      # freq value to plot
      freqplotval <- as.numeric(freqmut_v[f]) 
      # get the runs that have the selected freq
      freqcond <- freq==freqplotval 
      
      # start empty plot
      plot(-10, -10, type="b", pch=16, col=colvs[1],
           xlab="log(2Ns)", ylab="Prob(keep a)", ylim=c(0,1), xlim=c(min(selplot)*0.95,max(selplot)*1.05), xaxt='n', cex.lab=1.6, cex.axis=1.2)
      
      # loop through the lowest and highest migration rate
      for(m in c(1,4)) {
        # loop through the Diploid and Hemizygous chromosomes
        for(i in 1:length(dip)) {
          # condition
          eval <- chrm==dip[i] & 
            dominance==as.numeric(dominance_v[h]) & 
            migrate==as.numeric(migrate_v[m]) & 
            freqcond 
          # lines for cases with 2Ns>0
          lines(selplot[eval & selection > 0], 
                pkeep[eval & selection > 0], 
                type="b", pch=14+i, col=colvs[i], lwd=2, lty=aux[m])
          # points for case of 2Ns=0
          points(selplot[eval & selection == 0], 
                 pkeep[eval & selection == 0], 
                 type="p", pch=14+i, col=colvs[i], lwd=2, lty=aux[m])
        }
      }
      axis(side=1, at=c(0.5,1,1.5,2,2.5), labels = c("s=0",1.0,1.5,2.0,2.5), cex.axis=1.2)
      # legend("topleft", c("H","D",paste("f=",freqmut_v, sep="")), pch=c(-1,-1,14+c(1:length(freqmut_v))), lty=c(1,1,rep(0, times=length(freqmut_v))), col=c(colvs[2],colvs[1],rep(1,times=length(freqmut_v))))
    }
  }
  dev.off()
}


##########################################################################
# Sup. Fig. S5 
# comparison of H/D FST for linked sites 
# for different recombination rates
##########################################################################

# COMBINATION OF PARAMETERS
# migration rate
migrate_v <-c("0.000", "0.00034", "0.0017", "0.0034")
# selection of mutation under divergent selection
selmutben_v <- c("0.000", "0.0067", "0.01333" ,"0.02667", "0.05333", "0.0667", "0.1334")
# recombination rate
recrate_v <- c("2.5e-7","2.5e-8")
# initial frequency of allele a at locus under divergent selection
freqmut_v <- c("0.1", "0.5")
# dominance of allele a at locus under divergent selection
dominance_v <- c("0.01","0.5")
# time of split of the two populations
time_v <- c("2000")
# effective size (Ne) that corresponds to a 2Ne=1500
# of a population of 1000 individuals at a hemizygous locus 
# (500 females with two copies and 500 males with one copy)
Ne <- 750
# To avoid the log10(0) issue to represent neutral case and no-migration
# replace the zero migration and zero selection arbitrary values
nonzero_mig <- 10^(-0.7)
nonzero_sel <- 10^0.5
# all values larger than maxval replaced by maxval
maxval <- 1.5
# minimum and maximum values in plot
ylim_min <- 0
ylim_max <- maxval
# to avoid divisions by zero, replace 0 by small value
zeroval <- 1e-3
# title for plot
ytext <- "FST"
tagstat <- paste("H/D", ytext)
# define a relative error, 
# such that all values between 1*(1-relerror) and 1*relerror are set to 1
relerror <- 0.05
# END SETTINGS

  # Lower recombination
  # all simulations - not conditioning on retaining the allele a
  sumstat_lowr <- read.table(file=paste("./SummaryFiles/fst_FitM2_h_r2.5e-8sumstat.txt",sep=""), header=TRUE, na.strings = "NA")
  str(sumstat_lowr)
  summary(sumstat_lowr)
  
  # conditional on retaining the allele a
  sumstat_lowr_cb <- read.table(file=paste("./SummaryFiles/fst_FitM2_h_r2.5e-8sumstat_cb.txt",sep=""), header=TRUE, na.strings = "NA")
  str(sumstat_lowr_cb)
  summary(sumstat_lowr_cb)
  
  # check that the combination of parameters is the same in both
  # sumstat and sumstat_cb, i.e that the first 5 columns are the same
  if(sum(sumstat_lowr[,1:5]!=sumstat_lowr_cb[,1:5])) {
    stop("The two files do not have the same combination of parameters.")
  }
  
  
  # Higher recombination
  # all simulations - not conditioning on retaining the allele a
  sumstat <- read.table(file=paste("./SummaryFiles/fst_FitM2_h_r2.5e-7sumstat.txt",sep=""), header=TRUE, na.strings = "NA")
  str(sumstat)
  summary(sumstat)
  
  # conditional on retaining the allele a
  sumstat_cb <- read.table(file=paste("./SummaryFiles/fst_FitM2_h_r2.5e-7sumstat_cb.txt",sep=""), header=TRUE, na.strings = "NA")
  str(sumstat_cb)
  summary(sumstat_cb)
  
  # check that the combination of parameters is the same in both
  # sumstat and sumstat_cb, i.e that the first 5 columns are the same
  if(sum(sumstat[,1:5]!=sumstat_cb[,1:5])>0) {
    stop("The two files do not have the same combination of parameters.")
  }
  
  # Get only the cases with initial allele freq of 0.1 and 0.5
  eval <- sumstat$f==0.1 | sumstat$f==0.5
  sumstat <- sumstat[eval,]
  sumstat_cb <- sumstat_cb[eval,]
  
  # check that the combination of parameters is the same for low and high rec
  if(sum(sumstat_lowr[,c(1,3:5)]!=sumstat[,c(1,3:5)])>0) {
    stop("The two files do not have the same combination of parameters.")
  }
  
  
  # Get the combination of parameters for each case
  migrate <- sumstat$m
  dominance <- sumstat$d
  selection <- sumstat$s
  freq <- sumstat$f
  
  # define the values for scaled by Ne  selection and migration in log10 scale
  # transformation to selection values
  selplot <- as.numeric(format(2*Ne*selection, digits = 0))
  # transformation to migration values
  migplot <- as.numeric(format(2*Ne*migrate, digits=1))
  
  unique(selplot)
  unique(migplot)
  
  # indicate if ratio of stats is ploted in logscale
  logscale <- F
  # indicate if ratio of stats have reversed colours
  revcol <- F
  
  # define the breaks of the scale
  breaks <- c(0.5,0.75,0.95,1.05,1.25,1.50)
  
  # Plot the heatmaps of H/D statistic
  # higher recombination rate
  # all simulations
  relstat_sumstat <- plotheat_rel_nointerpol(dstat=sumstat$D_fst_mean,
                                                    hdstat=sumstat$HD_fst_mean,
                                                    tagstat=tagstat,
                                                    zeroval=zeroval,
                                                    breaks=breaks,
                                                    dominance_v=as.numeric(dominance_v),
                                                    freqmut_v=as.numeric(freqmut_v),
                                                    dom=dominance,
                                                    freq=freq,
                                                    sel=selplot,
                                                    mig=migplot,
                                                    revcol = revcol,
                                                    log10scale=logscale,
                                                    relerror=relerror,
                                                    xlabel="2*Ns*", ylabel="2*Nm*")
  
  
  # conditional on allele retention
  relstat_sumstat_cb <- plotheat_rel_nointerpol(dstat=sumstat_cb$D_fst_mean,
                                               hdstat=sumstat_cb$HD_fst_mean,
                                               tagstat=tagstat,
                                               zeroval=zeroval,
                                               breaks=breaks,
                                               dominance_v=as.numeric(dominance_v),
                                               freqmut_v=as.numeric(freqmut_v),
                                               dom=dominance,
                                               freq=freq,
                                               sel=selplot,
                                               mig=migplot,
                                               revcol = revcol,
                                               log10scale=logscale,
                                               relerror=relerror,
                                               xlabel="2*Ns*", ylabel="2*Nm*")
  # lower recombination rate
  # all simulations
  relstat_sumstat_lowr <- plotheat_rel_nointerpol(dstat=sumstat_lowr$D_fst_mean,
                                                 hdstat=sumstat_lowr$HD_fst_mean,
                                                 tagstat=tagstat,
                                                 zeroval=zeroval,
                                                 breaks=breaks,
                                                 dominance_v=as.numeric(dominance_v),
                                                 freqmut_v=as.numeric(freqmut_v),
                                                 dom=dominance,
                                                 freq=freq,
                                                 sel=selplot,
                                                 mig=migplot,
                                                 revcol = revcol,
                                                 log10scale=logscale,
                                                 relerror=relerror,
                                                 xlabel="2*Ns*", ylabel="2*Nm*")
  # conditional on allele retention
  relstat_sumstat_lowr_cb <- plotheat_rel_nointerpol(dstat=sumstat_lowr_cb$D_fst_mean,
                                                    hdstat=sumstat_lowr_cb$HD_fst_mean,
                                                    tagstat=tagstat,
                                                    zeroval=zeroval,
                                                    breaks=breaks,
                                                    dominance_v=as.numeric(dominance_v),
                                                    freqmut_v=as.numeric(freqmut_v),
                                                    dom=dominance,
                                                    freq=freq,
                                                    sel=selplot,
                                                    mig=migplot,
                                                    revcol = revcol,
                                                    log10scale=logscale,
                                                    relerror=relerror,
                                                    xlabel="2*Ns*", ylabel="2*Nm*")
  
  # merge all the figures into a list
  supfig5list <- c(relstat_sumstat_lowr$figs,
                   relstat_sumstat$figs,
                   relstat_sumstat_lowr_cb$figs,
                   relstat_sumstat_cb$figs)

# Arrange the plots and plot them to file
# for columns corresponding to lowr_q0.1 r_q0.1 lowr_q0.5 r_q0.5
# for rows corresponding to all_h0.01 all_h0.5 cond_h0.01 cond_h0.50
supfig5_v1 <- ggarrange(plotlist=supfig5list[c(1,5,2,6,3,7,4,8,
                                            9,13,10,14,11,15,12,16)], 
                     nrow=4, ncol=2*length(freqmut_v), 
                     common.legend = TRUE, legend="right")
supfig5_v1

# Arrange the plots and plot them to file
# for columns corresponding to all_q0.1 all_q0.5 cond_q0.1 cond_q0.5
# for rows corresponding to lowr_h0.01 r_h0.01 lowr_h0.5 r_h0.50
supfig5_v2 <- ggarrange(plotlist=supfig5list[c(1,2,9,10,
                                               5,6,13,14,
                                               3,4,11,12,
                                               7,8,15,16)], 
                        nrow=4, ncol=2*length(freqmut_v), 
                        common.legend = TRUE, legend="right")
supfig5_v2

# save SVG file that can be opened with Inkscape
svglite(file=paste("./Figures/SupFigS5.svg",sep=""), width=12, height=8)
supfig5_v2
dev.off()

# save PDF file
pdf(file=paste("./Figures/SupFigS5.pdf",sep=""), width=12, height=8)
supfig5_v2
dev.off()

