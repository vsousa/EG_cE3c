# load packages
library(stringdist)

#' GET_PIDXYFST_2DSFS
#' computes the estimate of pairwise differences, dxy and FST from the 2D SFS
#' @param sfs2d matrix with the 2D SFS
get_pidxyfst_2dsfs <- function(sfs2d) {
  # get a matrix with the frequencies 
  entries_sfs <-as.matrix(expand.grid(1:nrow(sfs2d),1:ncol(sfs2d)))
  # get the frequencies corresponding to each entry of the sfs
  freq <- entries_sfs-1
  # the first column corresponds to freq of pop2 and 2nd colum to freq of pop1 
  ss <- dim(sfs2d)-1
  # correcting factors to compute expected het
  cf <- ss/(ss-1)
  # transform absolute freq into relative allele frequencies
  for(i in 1:2) {
    freq[,i] <- freq[,i]/ss[i]
  }
  
  # compute number of pairwise differences as the sum of the cf*2*p*(1-p)
  het <- numeric(2)
  for(i in 1:2) {
    het_freq <- cf[i]*2*freq[,i]*(1-freq[,i])
    het[i] <- sum(het_freq*sfs2d[entries_sfs])
  }
  
  # compute dxy based on the frequencies from the SFS
  diff_freq <- (freq[,1]*(1-freq[,2]))+(freq[,2]*(1-freq[,1]))
  dxy <- sum(diff_freq*sfs2d[entries_sfs])
  
  # compute weighted fst based on the frequencies from the SFS
  numerator <- (freq[,1]-freq[,2])^2 - 
    ((freq[,1]*(1-freq[,1]))/(ss[1]-1)) -
    ((freq[,2]*(1-freq[,2]))/(ss[2]-1))
  denominator <-  diff_freq
  fst_weighted <- sum(numerator*sfs2d[entries_sfs])/
    sum(denominator*sfs2d[entries_sfs])
  
  # compute FST for each site
  hw <- (cf[1]*freq[,1]*(1-freq[,1]))+(cf[2]*freq[,2]*(1-freq[,2]))
  hb <- ((freq[,1]*(1-freq[,2]))+((1-freq[,1])*freq[,2]))
  # get the relative probability of each entry
  rel_p <- sfs2d[entries_sfs]/sum(sfs2d)
  # get the mean value defined as sum pi*xi, where xi is 1-(he/hb) and pi is the probability
  mean_fstsite <- sum(rel_p*(1-(hw/hb)), na.rm=TRUE)
  # get the variance as sum(pi*xi^2)-mean
  var_fstsite <- sum(rel_p*(1-(hw/hb))^2, na.rm=TRUE)-(mean_fstsite^2)
  
  list(pi=het, dxy=dxy, fst=fst_weighted, fst_site=mean_fstsite, fst_site_var=var_fstsite)
} 


#' Calculate 2D SFS from haplotype matrix (e.g. ms or scrm)
#'
#' @param haplo numeric matrix with nsites x ngene copies haplotype (0 or 1) - not genotypes
#' @param ss   vector of size 2 with numeric value with sample size for each pop
#'
#' @return matrix with 2D-SFS
#' @export
#'
#' @examples
getsfs2d <- function(haplo, ss) {
  freq1 <- rowSums(haplo[,1:ss[1]])
  freq2 <- rowSums(haplo[,(ss[1]+1):sum(ss)])
  # create a 2D table with the counts of sites
  tmp <- table(freq1,freq2)
  # create empty matrix with the 2D SFS
  sfs2d <- matrix(0, ncol=ss[2]+1, nrow=ss[1]+1)
  # get the index of rows and index of columns 
  # (+1 because the first index corresponds to 0 frequency)
  col_i <- as.numeric(colnames(tmp))+1
  row_i <- as.numeric(rownames(tmp))+1
  sfs2d[row_i,col_i] <- tmp
  sfs2d
}



#' Calculate mean PI for each pop, dxy and pairwise FST values
#'
#' This function computes FST according to Hudson's estimator following Bathia et al (2013).
#'
#' @param haplo numeric matrix with nsites x nind haplotype
#' @param ss   numeric vector of size 2 with sample size for each pop
#'
#' @return vector with pi_pop1, pi_pop2, dxy and pairwise FST between the two populations.
#' @export
#'
#' @examples
getMeanFst_Pi_Dxy <- function(haplo, ss) {
  freq1 <- rowSums(haplo[,1:ss[1]])/ss[1]
  freq2 <- rowSums(haplo[,(ss[1]+1):sum(ss)])/ss[2]
  p1 <- freq1*(1-freq1)
  p2 <- freq2*(1-freq2)
  # FST according to Bathia et al (2013)
  pdiffsquare <- (freq1-freq2)^2
  numerator <- pdiffsquare-(p1/(ss[1]-1))-(p2/(ss[2]-1))
  denominator <- (freq1*(1-freq2))+(freq2*(1-freq1))
  # Pi for each pop using correction factor
  cf <- ss/(ss-1)
  pi1 <- sum(2*cf[1]*p1)
  pi2 <- sum(2*cf[2]*p2)
  # Dxy 
  dxy <- sum(denominator)
  # Get FST
  if(dxy>0) {
    fst <- sum(numerator)/dxy
  } else {
    fst <- 0
  }
  c(pi1,pi2,dxy,fst)
}

#' Calculate Gmin statistic
#'
#' This function computes Gmin according to Geneva et al (2015).
#'
#' @param haplo numeric matrix with nsites x ngene copies haplotype (coded as 0, 1) - not genotypes
#' @param ss   numeric vector of size 2 with sample size of each pop
#'
#' @return Gmin between the two populations.
#' @export
#'
#' @examples
getGmin <- function(haplo, ss) {
  # transform the haplo into a vector of strings
  # each element of haplo_string contains the haplotype of a given individual
  haplo_string <- apply(haplo, 2, function(col) paste(col, collapse="")) 
  # # get the frequency of each haplotype
  # freq_hap <- table(haplo_string)
  # # get the unique haplotypes
  # unique_hap <- names(freq_hap)
  # freq_hap <- as.numeric(freq_hap)
  # # get haplotypes in pop 1
  # index_pop1 <- which(unique_hap %in% haplo_string[1:ss])
  # 
  # # get haplotypes in pop 2
  # index_pop2 <- which(unique_hap %in% haplo_string[(ss+1):(2*ss)])
  # 
  # get list of index of comparison of haplotypes
  #comparisons <- as.matrix(expand.grid(index_pop1,index_pop2))
  #dimnames(comparisons) <- NULL
  #apply(comparisons, 1, function(row) {unique_hap[row]})
  #identical(unique_hap[row[1]],unique_hap[row[2]])
  # dist_seq <- stringdistmatrix(unique_hap[index_pop1], unique_hap[index_pop2], method = "hamming")
  dist_seq <- stringdistmatrix(haplo_string[1:ss[1]], haplo_string[(ss[1]+1):sum(ss)], method = "hamming", nthread = 1)
  min_dxy <- min(dist_seq)
  mean_dxy <- mean(dist_seq)
  c(min_dxy,mean_dxy,min_dxy/mean_dxy)
}


#' Calculate FST value per site
#'
#' This function computes FST according to Hudson's estimator following Bathia (2013).
#'
#' @param haplo numeric matrix with nsites x nind haplotype
#' @param ss   numeric vector of size 2 with sample size of each pop
#'
#' @return FST between the two populations.
#' @export
#'
#' @examples
getFst <- function(haplo, ss) {
  freq1 <- rowSums(haplo[,1:ss[1], drop=FALSE])/ss[1]
  freq2 <- rowSums(haplo[,(ss[1]+1):sum(ss), drop=FALSE])/ss[2]
  p1 <- freq1*(1-freq1)
  p2 <- freq2*(1-freq2)
  pdiffsquare <- (freq1-freq2)^2
  numerator <- pdiffsquare-(p1/(ss[1]-1))-(p2/(ss[2]-1))
  denominator <- (freq1*(1-freq2))+(freq2*(1-freq1))
  res <- numerator/denominator
  # get the ones where denominator are larger than zero
  eval_zero <- denominator==0
  res[eval_zero] <- 0
  res
}


#' Calculate FST value per site
#'
#' This function computes FST according to Hudson's estimator following Hudson et al (1992).
#'
#' @param haplo numeric matrix with nsites x nind haplotype
#' @param ss   numeric vector with sample size for each pop
#'
#' @return FST between the two populations.
#' @export
#'
#' @examples
getFst_hudson <- function(haplo, ss) {
  freq1 <- rowSums(haplo[,1:ss[1], drop=FALSE])/ss[1]
  freq2 <- rowSums(haplo[,(ss[1]+1):sum(ss), drop=FALSE])/ss[2]
  # correction factor
  cf <- ss/(ss-1)
  hw <- ((cf[1]*freq1*(1-freq1))+(cf[2]*freq2*(1-freq2)))
  hb <- ((freq1*(1-freq2))+((1-freq1)*freq2))
  res <- 1-(hw/hb)
  # replace the sites with hb==0 by FST zero, otherwise we get NaN as we have a division by zero
  eval_zero <- hb==0
  res[eval_zero] <- 0 
  res
}

#' Calculate FST given allele freqs
#'
#' This function computes FST according to Hudson's estimator following Hudson et al (1992).
#'
#' @param freq1 numeric value with freq in pop1
#' @param freq2 numeric value with freq in pop2
#'
#' @return FST between the two populations.
#' @export
#'
#' @examples
getFst_freq_hudson <- function(freq1, freq2) {
  if(freq1==freq2) {
    res <- 0  
  } else {
    hw <- ((freq1*(1-freq1))+(freq2*(1-freq2)))
    hb <- ((freq1*(1-freq2))+((1-freq1)*freq2))
    res <- 1-(hw/hb)
  }
  res
}

#' Calculate FST given allele freqs and sample sizes
#'
#' This function computes FST according to Hudson's estimator following Hudson et al (1992).
#'
#' @param freq1 numeric value with freq in pop1
#' @param freq2 numeric value with freq in pop2
#' @param ss numeric vector with sample size of pop1 and pop2
#'
#' @return FST between the two populations.
#' @export
#'
#' @examples
getFst_freq_ss_hudson <- function(freq1, freq2, ss) {
  # correcting factors to compute expected het
  cf <- ss/(ss-1)
  hw <- (cf[1]*freq1*(1-freq1))+(cf[2]*freq2*(1-freq2))
  hb <- ((freq1*(1-freq2))+((1-freq1)*freq2))
  res <- 1-(hw/hb)
  # replace the sites with hb==0 by FST zero, otherwise we get NaN as we have a division by zero
  eval_zero <- hb==0
  res[eval_zero] <- NA 
  res
}

#' GETEXPECTEDNEUTRAL_STATS
#' function that computes the expected stats under a neutral IM model
#' simple approximations for FST
#' based on FST=(Tb-Tw)/Tb
#' where Tb is the coalescent time between 2 lineages sampled in different pops
#' and Tw is the coalescent time for 2 lineages sampled within the same pop
#' using eq. 17 in Zeng and Corcoran (2015) Genetics
#' but here following notation of eq. 22 from Wilkins-Herbots (2008)
#' @param N effective size in number of individuals in SLIM simulation
#' @param sr sex ratio defined as sr=Nm/(Nm+Nf), where Nm is size of males, Nf is size of females
#' @param U neutral mutation rate per site used in SLIM simulations
#' @param L sequence length in bp
#' @param window.size length of windows in bp
#' @param tsplit time of split in generations
#' @param m migration rate per generation (assumed to be symmetrical)
getexpectedneutral_stats <- function(N, sr, U, L, window.size, tsplit, m) {
  # get the Ne (EFFECTIVE SIZE) conditional on number of individuals
  # according to Mendez (2017) TPB "Differences in the effective population sizes of males and females do not require differences in their distribution of offspring number"
  # https://www.sciencedirect.com/science/article/pii/S0040580916300892#fig1
  # The effective size of autosome is:
  # NA=4NmNf/(Nm+Nf)
  # NX=9NmNf/(4Nm+2Nf)
  # The ratio of NX/NA depends on the sex ratio sr=Nm/(Nm+Nf)
  # NX/NA=(9/8)(1/(1+sr))
  ratioNxNa <- (9/8)*(1/(1+sr)) # ratio X/A
  # Effective size given N, sex-ratio and NX/NA
  Ne <- 4*N*sr*(1-sr)*ratioNxNa
  
  # ISOLATION MODEL
  if(m==0) {
    # expected theta for entire chromosome
    theta <- 4*Ne*U*L 
    # expected theta for windows 
    theta_w <- 4*Ne*U*window.size 
    # fst for isolation model, based on expected coalescent times
    expfst <- tsplit/(tsplit+(2*Ne)) 
    # simple approximations for DXY for isolation model
    # dxy for entire chromosome
    expdxy <- (tsplit+(2*Ne))*(U*L*2)
    # dxy per window
    expdxy_w <- (tsplit+(2*Ne))*(U*window.size*2) 
  } else if(m>0) {
    # expected fst and dxy depending for IM model
    # obtained from eq. 17 in Zeng and Corcoran (2015) Genetics
    # but here following eq. 22 from Wilkins-Herbots (2008)
    M <- 4*Ne*m
    tau <- tsplit/(2*Ne)
    D <- (2*M)^2 + 1
    lambda_plus <- ((2*M)+1+sqrt(D))/2
    lambda_minus <- ((2*M)+1-sqrt(D))/2
    elambda_plus <- exp(-lambda_plus*tau)
    elambda_minus <- exp(-lambda_minus*tau)
    lambda_minus_ratio <- 1-(1/lambda_minus)
    lambda_plus_ratio <- 1-(1/lambda_plus)
    firstterm <- (lambda_plus-1)*lambda_minus_ratio*elambda_minus
    secterm <- (1-lambda_minus)*lambda_plus_ratio*elambda_plus
    # time coalescent within for 2 lineages sampled within the same pop
    time_within <- (2+((1/sqrt(D))*(firstterm+secterm)))*(2*Ne)
    # time coalescent between for 2 lineages sampled from different pop
    firstterm <- lambda_plus*lambda_minus_ratio*elambda_minus
    secterm <- lambda_minus*lambda_plus_ratio*elambda_plus
    time_between <- (2+(1/M)+((1/sqrt(D))*(firstterm+secterm)))*(2*Ne)
    
    # expected theta for each population in IM model
    theta <- time_within*(2*U*L) # total chromosome
    theta_w <- time_within*(2*U*window.size) # per window
    # expected DXY in IM model
    expdxy <- time_between*(2*U*L)
    expdxy_w <- time_between*(2*U*window.size)
    # expected FST in IM model
    expfst <- (time_between-time_within)/time_between
    
  }
  
  # output the values
  list(pi=theta, pi_w=theta_w, dxy=expdxy, dxy_w=expdxy_w, fst=expfst)
}


#' T_COL
#' Transparent colors
#' function to make transparent colors from 
#' Mark Gardener 2015, taken from the followin website
#' www.dataanalytics.org.uk
#' @param	color color name
#' @param	percent percentage of transparency
#' @param	name an optional name for the color
t_col <- function(color, percent = 50, name = NULL) {
  
  # Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
  
  ## Save the color
  t.col
  
}
## END

