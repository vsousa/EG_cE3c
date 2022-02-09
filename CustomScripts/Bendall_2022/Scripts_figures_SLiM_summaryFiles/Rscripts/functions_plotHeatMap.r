#' PLOTHEAT_REL_NOINTERPOL
#' Function to plot heatmap of the relative H/D statistic
#' H - hemizygous, D - diploid
#' with a break in the x and y scale to deal with log10(0) values
#' i.e., using log scale for x and y axis with a break to distinguish the log(zero)
#' @param dsatat numeric vector of size n with the values of statistic for Diploid (autosome)
#' @param hdstat numeric vector of size n with the values of statistic for Hemizygous (X-chr or haplodiploid)
#' @param tagstat string with a tag that is the title of the plot
#' @param zeroval numeric value with the value used to replace zero stats (to avoid division by zero)
#' @param breaks  numeric value with the sequence of breaks breaks[1] is the mininum and breaks[length(breaks)] is the maximum value in the heatmap scale. All values above that are pooled.
#' @param dominance_v  numeric vector  with the dominance values considered
#' @param freqmut_v  numeric vector with the initial frequency values considered
#' @param dom numeric vector of size n with the dominance corresponding to each dstat and hdstat values
#' @param freq numeric vector of size n with the initial frequency corresponding to each dstat and hdstat values
#' @param sel numeric vector of size n with the selection coefficient corresponding to each dstat and hdstat values
#' @param mig numeric vector of size n with the migration rate corresponding to each dstat and hdstat values
#' @param revcol boolean (TRUE - the scale is reversed, blue to H/D>1)
#' @param log10scale boolean (TRUE - plot the heatmap scale in log10 scale)
#' @param relerror relative error such that all values between 1*(1-relerror) and 1*(1+relerror) are set to 1
#' @param xlabel string with the label of the x-axis
#' @param ylabel string with the label of the y-axis
plotheat_rel_nointerpol <- function(dstat, hdstat, tagstat, zeroval, 
                                   breaks, dominance_v, freqmut_v, 
                                   dom, freq, sel, mig, revcol, 
                                   log10scale=FALSE,
                                   relerror=0.05,
                                   xlabel="log(2*Ns*)", ylabel="log(2*Nm*)") {
  # load required package
  require(RColorBrewer)
  require(ggtext)
  require(scales)
  require(ggpubr)
  
  # lists to output figures and legend
  fig_df <- list() # save figure
  target <- list() # save title of figure (parameter combinations)
  count <- 1 # counter
  diploid_stat <- dstat # diploid statistic
  hd_stat <- hdstat     # haplodiploid statistic
  
  # to deal with divisions by zero replace zero by a small value 
  diploid_stat[diploid_stat==0] <- zeroval
  hd_stat[hd_stat==0] <- zeroval
  
  # for through the dominance coefficients and initial frequency combinations
  for(h in dominance_v) {
    for(f in 1:length(freqmut_v)) {
      # save the combination of parameters
      target[[count]] <- c(h,freqmut_v[f]) 
      print(target[[count]])
      # get the entries that correspond to HD 
      eval <- dom==h & freq==as.numeric(freqmut_v[f])
      sum(eval)
      
      # relative values as HD_stat/D_stat
      if(!log10scale) {
        tmp_scaled <- (hd_stat[eval])/(diploid_stat[eval])
        # replace NA values by 1
        # tmp_scaled[is.na(tmp_scaled)] <- 1
      } else {
        # relative ratio in log10 scale
        tmp_scaled <- log10(hd_stat[eval])-log10(diploid_stat[eval])
        # replace NA values by 0
        # tmp_scaled[is.na(tmp_scaled)] <- 0
      }
      
      # is set to 1.0
      settoone <- which(tmp_scaled>(1-relerror) & tmp_scaled<(1+relerror))
      tmp_scaled[settoone] <- 1
      
      # create data.frame with combination of selection and migration
      # and the ratio of the sumstat
      df <- data.frame(s=as.factor(sel[eval]), m=as.factor(mig[eval]), stat=tmp_scaled)
      
      # define the color scheme using RColorBrewer pallete "RdBu"
      # check if we use reverse scale (red-less than 1, blue-more than 1)
      # cor not (red-more than 1, blue-less than 1)
      if(!revcol) { 
        mycols <- rev(brewer.pal(n=10,"RdBu"))
        colfunc<-colorRampPalette(mycols)
      } else {
        mycols <- brewer.pal(n=10,"RdBu")
        colfunc<-colorRampPalette(mycols)
      }
      
      # plot the heatmaps 
      # define minimum and maximum values for the stat in the plotted scale
      max_tmp <- breaks[length(breaks)]
      min_tmp <- breaks[1]
      
      # save plot into list
      fig_tmp <-
        ggplot(data=df, aes(s, m, fill=stat))+
        geom_tile() + 
        theme_bw() +
        xlab(xlabel) +
        ylab(ylabel) +
        # labs(title = paste("h", h, "_f", freqmut_v[f], sep="")) +
        theme(text=element_text(size=12), 
              axis.text.x=element_markdown(size=10), 
              axis.text.y=element_markdown(size=10), 
              axis.title.x=element_markdown(size=12), 
              # axis.title.y=element_markdown(size=12, hjust=0), # for Figure 2 uncomment
              axis.title.y=element_markdown(size=12), # for figure 2 comment
              legend.text = element_text(size=10), legend.position="right")
      
      if(log10scale) {
        fig_df[[count]] <- fig_tmp +
          scale_fill_gradientn(colours = c(mycols[1:5],grey(0.9),grey(0.9),mycols[6:10]), values=rescale(breaks),
                               limits=c(min_tmp*1.02,max_tmp*1.02), name=paste(tagstat), na.value = "white",
                               oob=squish, breaks=breaks) #+
        # uncomment to add contour lines  
        # scale_colour_gradientn(colours = c(mycols[1:5],grey(0.9),grey(0.9),mycols[6:10]), values=rescale(breaks),
        #                        limits=c(min_tmp*1.02,max_tmp*1.02), name=paste(tagstat), 
        #                        oob=squish, breaks=breaks)
        # 
      } else {
        fig_df[[count]] <-  fig_tmp +
          scale_fill_gradientn(colours = c(mycols[1:5],grey(0.9),grey(0.9),mycols[6:10]), values=rescale(breaks),
                               limits=c(min_tmp*0.98,max_tmp*1.02), name=paste(tagstat), na.value = "white",
                               oob=squish, breaks=breaks) #+
        # uncomment to add contour lines  
        # scale_colour_gradientn(colours =  c(mycols[1:5],grey(0.9),grey(0.9),mycols[6:10]), values=rescale(breaks),
        #                        limits=c(min_tmp*0.98,max_tmp*1.02), name=paste(tagstat),
        #                        oob=squish, breaks=breaks)
      }
      
      count <- count + 1
    }
  }
  
  return(list(figs=fig_df, target=target))
}


#' PLOTHEAT_DIFF_IQR_NOINTERPOL
#' Function to plot heatmap of the relative (H-D)/IQR(D) statistic
#' that is the difference between Haplodiploid and Diploid 
#' standardized by the inter quantile range (IQR) of the diploid
#' H - hemizygous, D - diploid
#' with a break in the x and y scale to deal with log10(0) values
#' i.e., using log scale for x and y axis with a break to distinguish the log(zero)
#' @param dstat numeric vector of size n with the values of statistic for Diploid (autosome)
#' @param hdstat numeric vector of size n with the values of statistic for Hemizygous (X-chr or haplodiploid)
#' @param dq025 numeric vector of size n with the values of the quantile 0.25 for Diploid
#' @param dq075 numeric vector of size n with the values of the quantile 0.75 for Diploid
#' @param tagstat string with a tag that is the title of the plot
#' @param breaks  numeric value with the sequence of breaks breaks[1] is the mininum and breaks[length(breaks)] is the maximum value in the heatmap scale. All values above that are pooled.
#' @param dominance_v  numeric vector  with the dominance values considered
#' @param freqmut_v  numeric vector with the initial frequency values considered
#' @param dom numeric vector of size n with the dominance corresponding to each dstat and hdstat values
#' @param freq numeric vector of size n with the initial frequency corresponding to each dstat and hdstat values
#' @param sel numeric vector of size n with the selection coefficient corresponding to each dstat and hdstat values
#' @param mig numeric vector of size n with the migration rate corresponding to each dstat and hdstat values
#' @param revcol boolean (TRUE - the scale is reversed, blue to H/D>1)
#' @param log10scale boolean (TRUE - plot the heatmap scale in log10 scale)
#' @param relerror relative error such that all values between -relerror and relerror are set to 0
#' @param xlabel string with the label of the x-axis
#' @param ylabel string with the label of the y-axis
plotheat_diff_iqr_nointerpol <- function(dstat, hdstat, dq025, dq075, tagstat, 
                                    breaks, dominance_v, freqmut_v, 
                                    dom, freq, sel, mig, revcol, 
                                    log10scale=FALSE,
                                    relerror=0.05,
                                    xlabel="log(2*Ns*)", ylabel="log(2*Nm*)") {
  # load required package
  require(RColorBrewer)
  require(ggtext)
  require(scales)
  require(ggpubr)
  
  # lists to output figures and legend
  fig_df <- list() # save figure
  target <- list() # save title of figure (parameter combinations)
  count <- 1 # counter
  diploid_stat <- dstat # diploid statistic
  hd_stat <- hdstat     # haplodiploid statistic
  
  # for through the dominance coefficients and initial frequency combinations
  for(h in dominance_v) {
    for(f in 1:length(freqmut_v)) {
      # save the combination of parameters
      target[[count]] <- c(h,freqmut_v[f]) 
      print(target[[count]])
      # get the entries that correspond to HD 
      eval <- dom==h & freq==as.numeric(freqmut_v[f])
      sum(eval)
      
      # relative values as HD_stat/D_stat
      if(!log10scale) {
        tmp_scaled <- (hd_stat[eval]-diploid_stat[eval])/(dq075[eval]-dq025[eval])
      } else {
        # relative ratio in log10 scale
        tmp_scaled <- log10(hd_stat[eval]-diploid_stat[eval])-log10(dq075[eval]-dq025[eval])
      }
      
      # is set to 0
      # settoone <- which(tmp_scaled > -relerror & tmp_scaled < relerror)
      # tmp_scaled[settoone] <- 0
      # 
      
      # create data.frame with combination of selection and migration
      # and the ratio of the sumstat
      df <- data.frame(s=as.factor(sel[eval]), m=as.factor(mig[eval]), stat=tmp_scaled)
      
      # define the color scheme using RColorBrewer pallete "RdBu"
      # check if we use reverse scale (red-less than 0, blue-more than 0)
      # cor not (red-more than 0, blue-less than 0)
      if(!revcol) { 
        mycols <- rev(brewer.pal(n=10,"RdBu"))
        colfunc<-colorRampPalette(mycols)
      } else {
        mycols <- brewer.pal(n=10,"RdBu")
        colfunc<-colorRampPalette(mycols)
      }
      
      # plot the heatmaps 
      # define minimum and maximum values for the stat in the plotted scale
      max_tmp <- breaks[length(breaks)]
      min_tmp <- breaks[1]
      
      # save plot into list
      fig_tmp <-
        ggplot(data=df, aes(s, m, fill=stat))+
        geom_tile() + 
        theme_bw() +
        xlab(xlabel) +
        ylab(ylabel) +
        theme(text=element_text(size=12), axis.title.y =element_markdown(), axis.title.x =element_markdown(), axis.text.x = element_markdown(), axis.text.y = element_markdown(), legend.position="right")
      
      
      if(log10scale) {
        fig_df[[count]] <- fig_tmp +
          scale_fill_gradientn(colours = c(mycols[1:5],grey(0.9),grey(0.9),mycols[6:10]), values=rescale(breaks),
                               limits=c(min_tmp*1.02,max_tmp*1.02), name=paste(tagstat), na.value = "white",
                               oob=squish, breaks=breaks) #+
        # uncomment to add contour lines  
        # scale_colour_gradientn(colours = c(mycols[1:5],grey(0.9),grey(0.9),mycols[6:10]), values=rescale(breaks),
        #                        limits=c(min_tmp*1.02,max_tmp*1.02), name=paste(tagstat), 
        #                        oob=squish, breaks=breaks)
        # 
      } else {
        fig_df[[count]] <-  fig_tmp +
          scale_fill_gradientn(colours = c(mycols[1:5],grey(0.9),grey(0.9),mycols[6:10]), values=rescale(breaks),
                               limits=c(min_tmp*0.98,max_tmp*1.02), name=paste(tagstat), na.value = "white",
                               oob=squish, breaks=breaks) #+
        # uncomment to add contour lines  
        # scale_colour_gradientn(colours =  c(mycols[1:5],grey(0.9),grey(0.9),mycols[6:10]), values=rescale(breaks),
        #                        limits=c(min_tmp*0.98,max_tmp*1.02), name=paste(tagstat),
        #                        oob=squish, breaks=breaks)
      }
      
      count <- count + 1
    }
  }
  
  return(list(figs=fig_df, target=target))
}


#' PLOTHEAT_REL_MYV_BREAK
#' Function to plot heatmap of the relative H/D statistic
#' H - hemizygous, D - diploid
#' with a break in the x and y scale to deal with log10(0) values
#' i.e., using log scale for x and y axis with a break to distinguish the log(zero)
#' @param dsatat numeric vector of size n with the values of statistic for Diploid (autosome)
#' @param hdstat numeric vector of size n with the values of statistic for Hemizygous (X-chr or haplodiploid)
#' @param tagstat string with a tag that is the title of the plot
#' @param zeroval numeric value with the value used to replace zero stats (to avoid division by zero)
#' @param breaks  numeric value with the sequence of breaks breaks[1] is the mininum and breaks[length(breaks)] is the maximum value in the heatmap scale. All values above that are pooled.
#' @param dominance_v  numeric vector  with the dominance values considered
#' @param freqmut_v  numeric vector with the initial frequency values considered
#' @param dom numeric vector of size n with the dominance corresponding to each dstat and hdstat values
#' @param freq numeric vector of size n with the initial frequency corresponding to each dstat and hdstat values
#' @param sel numeric vector of size n with the selection coefficient corresponding to each dstat and hdstat values
#' @param mig numeric vector of size n with the migration rate corresponding to each dstat and hdstat values
#' @param revcol boolean (TRUE - the scale is reversed, blue to H/D>1)
#' @param log10scale boolean (TRUE - plot the heatmap scale in log10 scale)
#' @param breakzero boolean (TRUE - plot a break to deal with log(0))
#' @param relerror relative error such that all values between 1*(1-relerror) and 1*(1+relerror) are set to 1
#' @param alpha the number of points on each axis of the grid where the interpolation is done
#' @param xlabel string with the label of the x-axis
#' @param ylabel string with the label of the y-axis
#' @param xticks vector with the positions for the x-axis tick marks
#' @param xticklab vector of strings with the labels for the x-axis tick marks
#' @param yticks vector with the positions for the y-axis tick marks
#' @param yticklab vector of strings with the labels for the y-axis tick marks
plotheat_rel_myv_break <- function(dstat, hdstat, tagstat, zeroval, 
                                   breaks, dominance_v, freqmut_v, 
                                   dom, freq, sel, mig, revcol, 
                                   log10scale=FALSE, breakzero=TRUE, 
                                   relerror=0.05, alpha=20,
                                   xlabel="log(2*Ns*)", ylabel="log(2*Nm*)",
                                   xticks=c(0.5,1.0,1.5,2.0),
                                   xticklab=c("*s*=0","1.0","1.5","2.0"),
                                   yticks=c(-0.7,-0.4,-0.2,0.0,0.2,0.4,0.6),
                                   yticklab=c("*m*=0","-0.4","-0.2","0.0","0.2","0.4","0.6")) {
  # load required package
  require(RColorBrewer)
  require(ggtext)
  require(scales)
  require(ggpubr)
  
  # lists to output figures and legend
  fig_df <- list() # save figure
  target <- list() # save title of figure (parameter combinations)
  count <- 1 # counter
  diploid_stat <- dstat # diploid statistic
  hd_stat <- hdstat     # haplodiploid statistic
  
  # to deal with divisions by zero replace zero by a small value 
  diploid_stat[diploid_stat==0] <- zeroval
  hd_stat[hd_stat==0] <- zeroval
  
  # for through the dominance coefficients and initial frequency combinations
  for(h in dominance_v) {
    for(f in 1:length(freqmut_v)) {
      # save the combination of parameters
      target[[count]] <- c(h,freqmut_v[f]) 
      print(target[[count]])
      # get the entries that correspond to HD 
      eval <- dom==h & freq==as.numeric(freqmut_v[f])
      sum(eval)
      
      # relative values as HD_stat/D_stat
      if(!log10scale) {
        tmp_scaled <- (hd_stat[eval])/(diploid_stat[eval])
        # replace NA values by 1
        # tmp_scaled[is.na(tmp_scaled)] <- 1
      } else {
        # relative ratio in log10 scale
        tmp_scaled <- log10(hd_stat[eval])-log10(diploid_stat[eval])
        # replace NA values by 0
        # tmp_scaled[is.na(tmp_scaled)] <- 0
      }
      
      
      # replace values larger than maxval by the maximum value
      # tmp_scaled[tmp_scaled >= maxval] <- maxval
      # # replace values smaller than minval by the minimum value
      # tmp_scaled[tmp_scaled <= minval] <- minval
      
      # define the error to be between 1-relerror and 1+relerror
      # is set to 1.0
      settoone <- which(tmp_scaled>(1-relerror) & tmp_scaled<(1+relerror))
      tmp_scaled[settoone] <- 1
       
      # create data.frame with combination of selection and migration
      # and the ratio of the sumstat
      df <- data.frame(s=sel[eval], m=mig[eval], stat=tmp_scaled)
      
      # To use a simple 2D-barycentric interpolation 
      # compute the range of values in the x-axis (selection)
      # and y-axis (migration)
      minx_int <- min(df$s)
      maxx_int <- max(df$s)
      miny_int <- min(df$m)
      maxy_int <- max(df$m)
      rangex <- abs(minx_int-maxx_int)
      rangey <- abs(miny_int-maxy_int)
      
      # grid of s and m values where stat values
      # will be interpolated
      # define the minimum number of entries in the axis with smaller range
      # alpha <- 20 # given as input to the function
      # define the grid of points for x (x0) and for y (y0)
      # this depends on whether the x axis or the y axis have a larger range
      if(rangex>rangey) {
        scalem <- round(rangex/rangey)
        xo=seq(minx_int, maxx_int, length.out=floor(alpha*scalem))
        yo=seq(miny_int, maxy_int, length.out=alpha)
      } else {
        scalem <- round(rangey/rangex)
        xo=seq(minx_int, maxx_int, length.out=alpha)
        yo=seq(miny_int, maxy_int, length.out=floor(alpha*scalem))
      }
      
      # perform the interpolation using the barycentric option
      # create a data.frame with combination of s and m
      X <- cbind(df$s,df$m) 
      # create a grid with the combination of values of x0 and y0 at which we will interpolate the value of stat
      Xi <- expand.grid(xo,yo) 
      # perform the barycentric interpolation, resulting in a vector
      # with the interpolated values of the statistic for each combination
      # of s and m parameter avlues in xo and yo
      interpstat <- interp.barycentric(X, df$stat, Xi)
      # convert the vector into a matrix (in case of using image() function)
      #tmpint <- matrix(NA, nrow=length(xo), ncol=length(yo))
      #tmpint[cbind(findInterval(Xi[,1],xo),findInterval(Xi[,2],yo))] <- interpstat
      
      # to use ggplot for the plots
      # convert into a data.frame
      dfinterp <- data.frame(s=Xi[,1], m=Xi[,2], stat=interpstat)
      
      # define the color scheme using RColorBrewer pallete "RdBu"
      # check if we use reverse scale (red-less than 1, blue-more than 1)
      # cor not (red-more than 1, blue-less than 1)
      if(!revcol) { 
        mycols <- rev(brewer.pal(n=10,"RdBu"))
        colfunc<-colorRampPalette(mycols)
      } else {
        mycols <- brewer.pal(n=10,"RdBu")
        colfunc<-colorRampPalette(mycols)
      }
      
      # add a break point to deal with the fact that 0 was replaced
      # by a given value to avoid log10 of zero
      # the second row and second column of the grid is replaced by NA
      # assuming that the smallest value in x and y of the grid correspond
      # to the replacement of log(0)
      if(breakzero) {
        dfinterp$stat[dfinterp$s==unique(dfinterp$s)[2]] <- NA
        dfinterp$stat[dfinterp$m==unique(dfinterp$m)[2]] <- NA
      }
      
      # plot the heatmaps 
      # define minimum and maximum values for the stat in the plotted scale
      max_tmp <- breaks[length(breaks)]
      min_tmp <- breaks[1]
      # add the break points at the mid point between the second and third
      auxx <- unique(dfinterp$s)[1:3]
      auxy <- unique(dfinterp$m)[1:3]
      
      # save plot into list
      fig_tmp <-
        ggplot(data=dfinterp, aes(s, m))+
        geom_raster(aes(fill = stat), interpolate = FALSE, hjust = 0.5, vjust = 0.5, na.rm=TRUE, show.legend = TRUE) +
        # uncomment to add title
        # ggtitle(paste(tagstat, paste(paste("h=", h, "f=", freqmut_v[f]), collapse = " "))) +
        # uncomment to add contour lines
        # geom_contour(aes(z = stat, colour=..level..), size=0.8, lineend = "round", linejoin = "mitre", breaks=breaks) +
        xlab(xlabel) +
        ylab(ylabel) +
        scale_x_continuous(breaks = xticks, labels = xticklab, expand = c(0,0)) +
        scale_y_continuous(breaks = yticks, labels = yticklab, expand = c(0,0)) +
        theme_bw() +
        # annotate("segment", size=1, x = (auxx[1]+auxx[2])/2, xend = 0.03+((auxx[1]+auxx[2])/2), y =auxy[1]-0.07, yend = auxy[1]-0.01) +
        # annotate("segment", size=1, x = 0.02+(auxx[1]+auxx[2])/2, xend = 0.05+((auxx[1]+auxx[2])/2), y =auxy[1]-0.07, yend = auxy[1]-0.01) +
        # annotate("segment", size=1, x = auxx[1]-0.05, xend = auxx[1], y = auxy[1]-0.02, yend = auxy[1]+0.02) +
        # annotate("segment", size=1, x = auxx[1]-0.05, xend = auxx[1], y = auxy[1], yend = auxy[1]+0.04) +
        theme(text=element_text(size=12), axis.title.y =element_markdown(), axis.title.x =element_markdown(), axis.text.x = element_markdown(), axis.text.y = element_markdown(), legend.position="right")
        
        
      if(log10scale) {
        fig_df[[count]] <- fig_tmp +
          scale_fill_gradientn(colours = c(mycols[1:5],grey(0.9),grey(0.9),mycols[6:10]), values=rescale(breaks),
                               limits=c(min_tmp*1.02,max_tmp*1.02), name=paste(tagstat), na.value = "white",
                               oob=squish, breaks=breaks) #+
          # uncomment to add contour lines  
          # scale_colour_gradientn(colours = c(mycols[1:5],grey(0.9),grey(0.9),mycols[6:10]), values=rescale(breaks),
          #                        limits=c(min_tmp*1.02,max_tmp*1.02), name=paste(tagstat), 
          #                        oob=squish, breaks=breaks)
          # 
      } else {
        fig_df[[count]] <-  fig_tmp +
          scale_fill_gradientn(colours = c(mycols[1:5],grey(0.9),grey(0.9),mycols[6:10]), values=rescale(breaks),
                               limits=c(min_tmp*0.98,max_tmp*1.02), name=paste(tagstat), na.value = "white",
                               oob=squish, breaks=breaks) #+
          # uncomment to add contour lines  
          # scale_colour_gradientn(colours =  c(mycols[1:5],grey(0.9),grey(0.9),mycols[6:10]), values=rescale(breaks),
          #                        limits=c(min_tmp*0.98,max_tmp*1.02), name=paste(tagstat),
          #                        oob=squish, breaks=breaks)
       }
      
      count <- count + 1
    }
  }
  
  return(list(figs=fig_df, target=target))
}




#' PLOTHEAT_VAR
#' Function to plot heatmap of the values of a statistic
#' @param stat numeric vector of size n with the values of statistic
#' @param tagstat string with a tag that is the title of the plot
#' @param zeroval numeric value with the value used to replace zero stats (to avoid division by zero)
#' @param maxval  numeric value with the maximum value in the heatmap scale. All values above that are pooled.
#' @param minval  numeric value with the minimum value in the heatmap scale. All values below that are pooled.
#' @param dominance_v  numeric vector  with the dominance values considered
#' @param freqmut_v  numeric vector with the initial frequency values considered
#' @param dom numeric vector of size n with the dominance corresponding to each dstat and hdstat values
#' @param freq numeric vector of size n with the initial frequency corresponding to each dstat and hdstat values
#' @param sel numeric vector of size n with the selection coefficient corresponding to each dstat and hdstat values
#' @param mig numeric vector of size n with the migration rate corresponding to each dstat and hdstat values
#' @param revcol boolean (TRUE - the scale is reversed, blue to H/D>1)
#' @param log10scale boolean (TRUE - plot the heatmap scale in log10 scale)
#' @param breakzero boolean (TRUE - plot a break to deal with log(0))
plotheat_var <- function(stat, tagstat, zeroval, maxval, minval, 
                         dominance_v, freqmut_v, dom, freq, 
                         sel, mig, revcol, log10scale=FALSE, breakzero=FALSE) {
  # lists to output
  fig_df <- list()
  target <- list()
  count <- 1 # counter
  
  # loop through the combination of dominance coefficients and initial frequencies
  for(h in dominance_v) {
    for(f in 1:length(freqmut_v)) {
      # save the combination of parameters
      target[[count]] <- c(h,freqmut_v[f]) 
      print(target[[count]])
      # get the entries that correspond to the combination of parameters
      eval <- dom==h & freq==as.numeric(freqmut_v[f])
      sum(eval)

      # transform values if needed 
      if(!log10scale) {
        tmp_scaled <- stat[eval]
      } else {
        tmp_scaled <- log10(stat[eval])
      }
      
      # replace values larger than maxval by the maximum value
      tmp_scaled[tmp_scaled >= maxval] <- maxval
      # replace values smaller than minval by the minimum value
      tmp_scaled[tmp_scaled <= minval] <- minval
      
      # create data.frame with combination of selection and migration
      # and the ratio of the sumstat
      df <- data.frame(s=sel[eval], m=mig[eval], stat=tmp_scaled)
      
      # To use a simple 2D-barycentric interpolation 
      # compute the range of values in the x-axis (selection)
      # and y-axis (migration)
      minx_int <- min(df$s)
      maxx_int <- max(df$s)
      miny_int <- min(df$m)
      maxy_int <- max(df$m)
      rangex <- abs(minx_int-maxx_int)
      rangey <- abs(miny_int-maxy_int)
      
      # grid of s and m values where stat values
      # will be interpolated
      # define the minimum number of entries in the axis with smaller range
      alpha <- 30
      # define the grid of points for x (x0) and for y (y0)
      # this depends on whether the x axis or the y axis have a larger range
      if(rangex>rangey) {
        scalem <- round(rangex/rangey)
        # manual linear interpolation
        xo=seq(minx_int, maxx_int, length.out=floor(alpha*scalem))
        yo=seq(miny_int, maxy_int, length.out=alpha)
      } else {
        scalem <- round(rangey/rangex)
        # manual linear interpolation
        xo=seq(minx_int, maxx_int, length.out=alpha)
        yo=seq(miny_int, maxy_int, length.out=floor(alpha*scalem))
      }
      
      # perform the interpolation using the barycentric option
      # create a data.frame with combination of s and m
      X <- cbind(df$s,df$m) 
      # create a grid with the combination of values of x0 and y0 at which we will interpolate the value of stat
      Xi <- expand.grid(xo,yo) 
      # perform the barycentric interpolation, resulting in a vector
      # with the interpolated values of the statistic for each combination
      # of s and m parameter avlues in xo and yo
      interpstat <- interp.barycentric(X, df$stat, Xi)
      # convert the vector into a matrix (in case of using image() function)
      #tmpint <- matrix(NA, nrow=length(xo), ncol=length(yo))
      #tmpint[cbind(findInterval(Xi[,1],xo),findInterval(Xi[,2],yo))] <- interpstat
      
      # to use ggplot for the plots
      # convert into a data.frame
      dfinterp <- data.frame(s=Xi[,1], m=Xi[,2], stat=interpstat)

      # define the color scheme using RColorBrewer pallete "YlGn"
      if(!revcol) {
        mycols <- rev(brewer.pal(n=9,"YlGn"))
        colfunc<-colorRampPalette(mycols)
      } else {
        mycols <- brewer.pal(n=9,"YlGn")
        colfunc<-colorRampPalette(mycols)
      }
      
      # add a break point to deal with the fact that 0 was replaced
      # by a given value to avoid log10 of zero
      # the second row and second column of the grid is replaced by NA
      # assuming that the smallest value in x and y of the grid correspond
      # to the replacement of log(0)
      if(breakzero) {
        dfinterp$stat[dfinterp$s==unique(dfinterp$s)[2]] <- NA
        dfinterp$stat[dfinterp$m==unique(dfinterp$m)[2]] <- NA
      }
      
      # plot the heatmaps 
      # define minimum and maximum values for the stat in the plotted scale
      max_tmp <- maxval
      min_tmp <- minval
      if(log10scale) {
        fig_df[[count]] <-
          ggplot(data=dfinterp, aes(s, m))+
          geom_raster(aes(fill = stat), interpolate = FALSE, hjust = 0.0, vjust = 0.0) +
          ggtitle(paste(tagstat, paste(paste("h=", h, "f=", freqmut_v[f]), collapse = " "))) +
          geom_contour(aes(z = stat, colour=..level..), size=1, lineend = "round", linejoin = "mitre", breaks=seq(from=min_tmp*1.05,to=max_tmp*1.05, length.out = 16), color=gray(0.6)) +
          # geom_contour(aes(z = stat, colour=..level..), size=1.2, breaks=(min_tmp+max_tmp)/2, color="white") +
          xlab("log10(s)") +
          ylab("log10(m)") +
          #geom_point(data = df, aes(s, m), colour = 'white') +
          #scale_y_reverse() +
          scale_fill_gradientn(colours = mycols, values=rescale(c(min_tmp*1.05,(min_tmp+max_tmp)/2,max_tmp*1.05)),
                               limits=c(min_tmp*1.05,max_tmp*1.05), name=paste(tagstat), na.value = "white" ) +
          scale_colour_gradientn(colours = mycols, values=rescale(c(min_tmp*1.05,(min_tmp+max_tmp)/2,max_tmp*1.05)),
                                 limits=c(min_tmp*1.05,max_tmp*1.05), name=paste(tagstat))
        
      } else {
        fig_df[[count]] <-
          ggplot(data=dfinterp, aes(s, m))+
          geom_raster(aes(fill = stat), interpolate = FALSE, hjust = 0.0, vjust = 0.0) +
          ggtitle(paste(tagstat, paste(paste("h=", h, "f=", freqmut_v[f]), collapse = " "))) +
          geom_contour(aes(z = stat, colour=..level..), size=1, lineend = "round", linejoin = "mitre", breaks=seq(from=min_tmp*0.95,to=max_tmp*1.05, length.out = 16), color=gray(0.6)) +
          # geom_contour(aes(z = stat, colour=..level..), size=1.2, breaks=(min_tmp+max_tmp)/2, color="white") +
          xlab("log10(s)") +
          ylab("log10(m)") +
          #geom_point(data = df, aes(s, m), colour = 'white') +
          #scale_y_reverse() +
          scale_fill_gradientn(colours = mycols, values=rescale(c(min_tmp*0.95,(min_tmp+max_tmp)/2,max_tmp*1.05)),
                               limits=c(min_tmp*0.95,max_tmp*1.05), name=paste(tagstat), na.value = "white") +
          scale_colour_gradientn(colours = mycols, values=rescale(c(min_tmp*0.95,(min_tmp+max_tmp)/2,max_tmp*1.05)),
                                 limits=c(min_tmp*0.95,max_tmp*1.05), name=paste(tagstat))
        
      }
      
      
      count <- count + 1
    }
  }
  
  return(list(figs=fig_df, target=target))
}


#' PLOT_FST_POINTS_v1
#' Function to plot FST as a function of selective coefficient
#' with points with bars corresponding to the interquantile range
#' getting as input a vector with data for A and X and 
#' a vector chrm indicating which entries correspond to A and X
#' @param mstat numeric vector of size n with the values of mean statistic
#' @param q25stat numeric vector of size n with the values of quantile 0.25 statistic
#' @param q75stat numeric vector of size n with the values of quantile 0.75 statistic
#' @param tagstat string with a tag that is the title of the plot
#' @param maxval  numeric value with the maximum FST value 
#' @param minval  numeric value with the minimum FST value (default=0)
#' @param dominance_v  numeric vector  with the dominance values considered
#' @param dom numeric vector of size n with the dominance corresponding to each dstat and hdstat values
#' @param mig numeric vector of size n with the migration rate corresponding to each dstat and hdstat values
#' @param sel numeric vector of size n with the selection coefficient corresponding to each dstat and hdstat values
#' @param chrm numeric vector of size n with the chromosome (A or X) corresponding to each dstat and hdstat values
#' @param case_mig numeric vector of size m with the migration rates that are plotted
#' @param colsvs numeric vector of size 2 with the colors for Haplodiploid (X) and Diploid (A)
#' @param pointsize numeric value with size of point
plot_fst_points_v1 <- function(mstat, q25stat, q75stat, tagstat, maxval, minval=0, 
                         dominance_v, dom, mig, sel, chrm, case_mig, colsvs, pointsize) {
  plotfst_point <- list()
  count <- 1
  for(m in 1:length(case_mig)) {
    for(h in dominance_v) {
      # select the trajectory for D and H that correspond to the selected case
      eval_d <- which(mig==case_mig[m] & dom==h & chrm=="A")
      eval_h <- which(mig==case_mig[m] & dom==h & chrm=="X")
      
      # create a data.frame with the required info
      # diploid
      dfd <- data.frame(chrm=rep("D", length(unique(sel))), 
                        s=as.factor(round(2*Ne*sel[eval_d])), 
                        fst_mean=mstat[eval_d], 
                        fst_25=q25stat[eval_d],
                        fst_75=q75stat[eval_d])
      # hemizygous 
      dfh <- data.frame(chrm=rep("H", length(unique(sel))), 
                        s=as.factor(round(2*Ne*sel[eval_h])), 
                        fst_mean=mstat[eval_h], 
                        fst_25=q25stat[eval_h],
                        fst_75=q75stat[eval_h])
      # merge the two data.frames
      df <- rbind(dfd, dfh)
      
      
      # plot the fst using points
      plotfst_point[[count]] <- ggplot(data=df, aes(x=s, y=fst_mean, group=chrm, colour=chrm)) + 
        geom_pointrange(aes(ymin=fst_25, ymax=fst_75, group=chrm, fill=chrm), position=position_dodge(0.6), show.legend = TRUE, size=pointsize)  +
        theme_classic() +
        scale_color_manual(values=colvs) + 
        scale_y_continuous(limits = c(minval, maxval)) + 
        labs(y="Differentiation (F<sub>ST</sub>)", x = "Selection (2*Ns*)") +
        labs(colour = "") + 
        theme(axis.text.x=element_text(size=10), 
              axis.text.y=element_text(size=10), 
              axis.title.x=element_markdown(size=12), 
              axis.title.y=element_markdown(size=12), 
              legend.text = element_text(size=10))
      count <- count + 1
    }
  }
  # return the list of figures
  plotfst_point
}


#' PLOT_FST_POINTS_v2
#' Function to plot FST as a function of selective coefficient
#' with points with bars corresponding to the interquantile range
#' getting as input different vectors with data for A and X separately 
#' without a vector chrm indicating which entries correspond to A and X
#' as there is a one to one correspondence
#' @param mstatA numeric vector of size n with the values of mean statistic for diploid
#' @param q25statA numeric vector of size n with the values of quantile 0.25 statistic for diploid
#' @param q75statA numeric vector of size n with the values of quantile 0.75 statistic for diploid
#' @param mstatX numeric vector of size n with the values of mean statistic for haplodiploid
#' @param q25statX numeric vector of size n with the values of quantile 0.25 statistic for haplodiploid
#' @param q75statX numeric vector of size n with the values of quantile 0.75 statistic for haplodiploid
#' @param tagstat string with a tag that is the title of the plot
#' @param maxval  numeric value with the maximum FST value 
#' @param minval  numeric value with the minimum FST value (default=0)
#' @param dominance_v  numeric vector  with the dominance values considered
#' @param dom numeric vector of size n with the dominance corresponding to each dstat and hdstat values
#' @param mig numeric vector of size n with the migration rate corresponding to each dstat and hdstat values
#' @param sel numeric vector of size n with the selection coefficient corresponding to each dstat and hdstat values
#' @param case_mig numeric vector of size m with the migration rates that are plotted
#' @param colsvs numeric vector of size 2 with the colors for Haplodiploid (X) and Diploid (A)
#' @param pointsize numeric value with size of point
plot_fst_points_v2 <- function(mstatA, q25statA, q75statA, mstatX, q25statX, q75statX,
                               tagstat, maxval, minval=0, 
                               dominance_v, dom, mig, sel, case_mig, colsvs, pointsize) {
  plotfst_point <- list()
  count <- 1
  for(m in 1:length(case_mig)) {
    for(h in dominance_v) {
      # select the trajectory for D and H that correspond to the selected case
      eval <- which(mig==case_mig[m] & dom==h)
      
      # create a data.frame with the required info
      # diploid
      dfd <- data.frame(chrm=rep("D", length(unique(sel))), 
                        s=as.factor(round(2*Ne*sel[eval])), 
                        fst_mean=mstatA[eval], 
                        fst_25=q25statA[eval],
                        fst_75=q75statA[eval])
      # hemizygous 
      dfh <- data.frame(chrm=rep("H", length(unique(sel))), 
                        s=as.factor(round(2*Ne*sel[eval])), 
                        fst_mean=mstatX[eval], 
                        fst_25=q25statX[eval],
                        fst_75=q75statX[eval])
      # merge the two data.frames
      df <- rbind(dfd, dfh)
      
      
      # plot the fst using points
      plotfst_point[[count]] <- ggplot(data=df, aes(x=s, y=fst_mean, group=chrm, colour=chrm)) + 
        geom_pointrange(aes(ymin=fst_25, ymax=fst_75, group=chrm, fill=chrm), position=position_dodge(0.6), show.legend = TRUE, size=pointsize)  +
        theme_classic() +
        scale_color_manual(values=colvs) + 
        scale_y_continuous(limits = c(minval, maxval)) + 
        labs(y="Differentiation (F<sub>ST</sub>)", x = "Selection (2*Ns*)") +
        labs(colour = "") + 
        theme(axis.text.x=element_text(size=10), 
              axis.text.y=element_text(size=10), 
              axis.title.x=element_markdown(size=12), 
              axis.title.y=element_markdown(size=12), 
              legend.text = element_text(size=10))
      count <- count + 1
    }
  }
  # return the list of figures
  plotfst_point
}


#' PLOT_FST_LINES_v1
#' Function to plot FST as a function of selective coefficient
#' with lines and shaded area corresponding to the interquantile range
#' getting as input a vector with data for A and X and 
#' a vector chrm indicating which entries correspond to A and X
#' @param mstat numeric vector of size n with the values of mean statistic
#' @param q25stat numeric vector of size n with the values of quantile 0.25 statistic
#' @param q75stat numeric vector of size n with the values of quantile 0.75 statistic
#' @param tagstat string with a tag that is the title of the plot
#' @param maxval  numeric value with the maximum value in the heatmap scale. All values above that are pooled.
#' @param minval  numeric value with the minimum value in the heatmap scale. All values below that are pooled.
#' @param dominance_v  numeric vector  with the dominance values considered
#' @param dom numeric vector of size n with the dominance corresponding to each dstat and hdstat values
#' @param mig numeric vector of size n with the migration rate corresponding to each dstat and hdstat values
#' @param sel numeric vector of size n with the selection coefficient corresponding to each dstat and hdstat values
#' @param chrm numeric vector of size n with the chromosome (A or X) corresponding to each dstat and hdstat values
#' @param case_mig numeric vector of size m with the migration rates that are plotted
#' @param colsvs numeric vector of size 2 with the colors for Haplodiploid (X) and Diploid (A)
plot_fst_lines_v1 <- function(mstat, q25stat, q75stat, tagstat, maxval, minval=0, 
                               dominance_v, dom, mig, sel, chrm, case_mig, colsvs) {
  plotfst_lines <- list()
  count <- 1
  for(m in 1:length(case_mig)) {
    for(h in dominance_v) {
      # select the trajectory for D and H that correspond to the selected case
      eval_d <- which(mig==case_mig[m] & dom==h & chrm=="A")
      eval_h <- which(mig==case_mig[m] & dom==h & chrm=="X")
      
      # create a data.frame with the required info
      # diploid
      dfd <- data.frame(chrm=rep("D", length(unique(sel))), 
                        s=2*Ne*sel[eval_d], 
                        fst_mean=mstat[eval_d], 
                        fst_25=q25stat[eval_d],
                        fst_75=q75stat[eval_d])
      # hemizygous 
      dfh <- data.frame(chrm=rep("H", length(unique(sel))), 
                        s=2*Ne*sel[eval_h], 
                        fst_mean=mstat[eval_h], 
                        fst_25=q25stat[eval_h],
                        fst_75=q75stat[eval_h])
      # merge the two data.frames
      df <- rbind(dfd, dfh)
      
      
      # plot the fst using lines
      plotfst_lines[[count]] <- ggplot(data=df, aes(x=s, y=fst_mean, group=chrm, colour=chrm)) + 
        geom_line(size=1)  + 
        geom_ribbon(aes(ymin=fst_25, ymax=fst_75, group=chrm, fill=chrm), alpha=0.2) +
        theme_classic() +
        scale_color_manual(values=colvs) + 
        scale_fill_manual(values=colvs) + 
        scale_y_continuous(limits = c(minval, maxval)) + 
        labs(y="F<sub>ST</sub>", x = "2*Ns*") +
        labs(colour = "") + 
        theme(axis.text.x=element_text(size=10), 
              axis.text.y=element_text(size=10), 
              axis.title.x=element_markdown(size=10), 
              axis.title.y=element_markdown(size=10), 
              legend.text = element_text(size=10))
      
      count <- count + 1
    }
  }
  # return the list of figures
  plotfst_lines
}

#' PLOT_FST_LINES_v2
#' Function to plot FST as a function of selective coefficient
#' with lines and shaded area corresponding to the interquantile range
#' getting as input different vectors with data for A and X separately 
#' without a vector chrm indicating which entries correspond to A and X
#' as there is a one to one correspondence
#' @param mstatA numeric vector of size n with the values of mean statistic for diploid
#' @param q25statA numeric vector of size n with the values of quantile 0.25 statistic for diploid
#' @param q75statA numeric vector of size n with the values of quantile 0.75 statistic for diploid
#' @param mstatX numeric vector of size n with the values of mean statistic for haplodiploid
#' @param q25statX numeric vector of size n with the values of quantile 0.25 statistic for haplodiploid
#' @param q75statX numeric vector of size n with the values of quantile 0.75 statistic for haplodiploid
#' @param tagstat string with a tag that is the title of the plot
#' @param maxval  numeric value with the maximum value in the heatmap scale. All values above that are pooled.
#' @param minval  numeric value with the minimum value in the heatmap scale. All values below that are pooled.
#' @param dominance_v  numeric vector  with the dominance values considered
#' @param dom numeric vector of size n with the dominance corresponding to each dstat and hdstat values
#' @param mig numeric vector of size n with the migration rate corresponding to each dstat and hdstat values
#' @param sel numeric vector of size n with the selection coefficient corresponding to each dstat and hdstat values
#' @param case_mig numeric vector of size m with the migration rates that are plotted
#' @param colsvs numeric vector of size 2 with the colors for Haplodiploid (X) and Diploid (A)
plot_fst_lines_v2 <- function(mstatA, q25statA, q75statA, mstatX, q25statX, q75statX,
                              tagstat, maxval, minval=0, 
                              dominance_v, dom, mig, sel,case_mig, colsvs) {
  plotfst_lines <- list()
  count <- 1
  for(m in 1:length(case_mig)) {
    for(h in dominance_v) {
      # select the trajectory for D and H that correspond to the selected case
      eval <- which(mig==case_mig[m] & dom==h)
      
      # create a data.frame with the required info
      # diploid
      dfd <- data.frame(chrm=rep("D", length(unique(sel))), 
                        s=2*Ne*sel[eval], 
                        fst_mean=mstatA[eval], 
                        fst_25=q25statA[eval],
                        fst_75=q75statA[eval])
      # hemizygous 
      dfh <- data.frame(chrm=rep("H", length(unique(sel))), 
                        s=2*Ne*sel[eval], 
                        fst_mean=mstatX[eval], 
                        fst_25=q25statX[eval],
                        fst_75=q75statX[eval])
      # merge the two data.frames
      df <- rbind(dfd, dfh)
      
      
      # plot the fst using lines
      plotfst_lines[[count]] <- ggplot(data=df, aes(x=s, y=fst_mean, group=chrm, colour=chrm)) + 
        geom_line(size=1)  + 
        geom_ribbon(aes(ymin=fst_25, ymax=fst_75, group=chrm, fill=chrm), alpha=0.2) +
        theme_classic() +
        scale_color_manual(values=colvs) + 
        scale_fill_manual(values=colvs) + 
        scale_y_continuous(limits = c(minval, maxval)) + 
        labs(y="F<sub>ST</sub>", x = "2*Ns*") +
        labs(colour = "") + 
        theme(axis.text.x=element_text(size=10), 
              axis.text.y=element_text(size=10), 
              axis.title.x=element_markdown(size=10), 
              axis.title.y=element_markdown(size=10), 
              legend.text = element_text(size=10))
      
      count <- count + 1
    }
  }
  # return the list of figures
  plotfst_lines
}



