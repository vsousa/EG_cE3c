# load packages
library(stringdist)

#' GET_PIDXYFST_2DSFS
#' computes the estimate of pairwise differences, dxy and FST from the 2D SFS
#' @param sfs2d matrix with the 2D SFS
get_pidxyfst_2dsfs_rel <- function(sfs2d) {
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
  # get the relative probability of each entry
  rel_p <- sfs2d[entries_sfs]/sum(sfs2d)
  
  
  # compute number of pairwise differences as the sum of the cf*2*p*(1-p)
  het <- numeric(2)
  for(i in 1:2) {
    het_freq <- cf[i]*2*freq[,i]*(1-freq[,i])
    #het[i] <- sum(het_freq*sfs2d[entries_sfs])
    het[i] <- sum(het_freq*rel_p)
  }
  
  # compute dxy based on the frequencies from the SFS
  diff_freq <- (freq[,1]*(1-freq[,2]))+(freq[,2]*(1-freq[,1]))
  #dxy <- sum(diff_freq*sfs2d[entries_sfs])
  dxy <- sum(diff_freq*rel_p)
  
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
  # remove the monomorphic sites
  rel_p[c(1, length(rel_p))] <- 0
  # get the relative just with polymorphic sites
  rel_p_poly <- rel_p/sum(rel_p)
  # get the mean value defined as sum pi*xi, where xi is 1-(he/hb) and pi is the probability
  #mean_fstsite <- sum(rel_p*(1-(hw/hb)), na.rm=TRUE)
  mean_fstsite <- sum(rel_p_poly*(1-(hw/hb)), na.rm=TRUE)
  # get the variance as sum(pi*xi^2)-mean
  var_fstsite <- sum(rel_p_poly*(1-(hw/hb))^2, na.rm=TRUE)-(mean_fstsite^2)
  
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
  freq1 <- rowSums(haplo[,1:ss[1], drop=FALSE])/ss[1]
  freq2 <- rowSums(haplo[,(ss[1]+1):sum(ss), drop=FALSE])/ss[2]
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
