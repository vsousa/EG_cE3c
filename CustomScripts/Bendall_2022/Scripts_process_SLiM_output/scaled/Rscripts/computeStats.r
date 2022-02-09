# load packages
library(stringdist)

#' Calculate mean PI for each pop, dxy and pairwise FST values
#'
#' This function computes FST according to Hudson's estimator following Bhatia et al (2013). doi: 10.1101/gr.154831.113 
#'
#' @param haplo numeric matrix with nsites x nind haplotype
#' @param ss   numeric value with sample size (assume equal for both pops)
#'
#' @return vector with pi_pop1, pi_pop2, dxy and pairwise FST between the two populations.
#' @export
#'
#' @examples
getMeanFst_Pi_Dxy <- function(haplo, ss) {
  freq1 <- rowSums(haplo[,1:ss])/ss
  freq2 <- rowSums(haplo[,(ss+1):(2*ss)])/ss
  p1 <- freq1*(1-freq1)
  p2 <- freq2*(1-freq2)
  # FST according to Bathia et al (2013)
  pdiffsquare <- (freq1-freq2)^2
  numerator <- pdiffsquare-(p1/(ss-1))-(p2/(ss-1))
  denominator <- (freq1*(1-freq2))+(freq2*(1-freq1))
  # Pi for each pop using correction factor
  cf <- ss/(ss-1)
  pi1 <- sum(2*cf*p1)
  pi2 <- sum(2*cf*p2)
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


#' Calculate FST value per site
#'
#' This function computes FST according to Hudson's estimator following Bathia (2013).
#'
#' @param haplo numeric matrix with nsites x nind haplotype
#' @param ss   numeric value with sample size (assume equal for both pops)
#'
#' @return FST between the two populations.
#' @export
#'
#' @examples
getFst <- function(haplo, ss) {
  freq1 <- rowSums(haplo[,1:ss])/ss
  freq2 <- rowSums(haplo[,(ss+1):(2*ss)])/ss
  p1 <- freq1*(1-freq1)
  p2 <- freq2*(1-freq2)
  pdiffsquare <- (freq1-freq2)^2
  numerator <- pdiffsquare-(p1/(ss-1))-(p2/(ss-1))
  denominator <- (freq1*(1-freq2))+(freq2*(1-freq1))
  res <- numerator/denominator
  res
}


#' Calculate FST value per site
#'
#' This function computes FST according to Hudson's estimator following Hudson et al (1992).
#'
#' @param haplo numeric matrix with nsites x nind haplotype
#' @param ss   numeric value with sample size (assume equal for both pops)
#'
#' @return FST between the two populations.
#' @export
#'
#' @examples
getFst_hudson <- function(haplo, ss) {
  freq1 <- rowSums(haplo[,1:ss, drop=FALSE])/ss
  freq2 <- rowSums(haplo[,(ss+1):(2*ss), drop=FALSE])/ss
  # correction factor
  cf <- ss/(ss-1)
  hw <- ((cf*freq1*(1-freq1))+(cf*freq2*(1-freq2)))
  hb <- ((freq1*(1-freq2))+((1-freq1)*freq2))
  res <- 1-(hw/hb)
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

#' Calculate FST Lasne et al
#'
#' This function computes FST following Lasne et al (2017) Evolution, eq 11.
#'
#' @param freq1 numeric value with freq in pop1
#' @param freq2 numeric value with freq in pop2
#'
#' @return FST between the two populations.
#' @export
#'
#' @examples
getFst_freq_lasne <- function(freq1, freq2) {
  deltafreq <- freq2 - freq1
  if(deltafreq==0) {
    res <- 0
  } else {
    meanfreq <- (freq1+freq2)/2
    res <- deltafreq^2 / (4*meanfreq*(1-meanfreq))  
  }
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
      # but here following eq. 22 from Wilkinson-Herbots (2008)  10.1016/j.tpb.2007.11.001 
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

