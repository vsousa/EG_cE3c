#' R2LINKAGE
#' compute the linkage disequilibrium between a pair of sites
#' INPUT
#' @param  genomatrix matrix with 2 columns corresponding to the haplotypes at 2 sites
#'                Each haplotype is a column. Mutations are coded as 0, 1. 
#' @return r2 for a pair of sites
r2linkage <- function(genomatrix) {
  # get the frequency of each SNP
  freq <- colSums(genomatrix)
  # get the frequency of haplotype ab
  freq_ab <- sum(rowSums(genomatrix)==2)
  # get the sample size in the population
  n <- nrow(genomatrix)
  # compute the r2
  r2_simplify <- (prod(freq)-(n*freq_ab))^2/(prod(freq)*(freq[1]-n)*(freq[2]-n))
  # r2_ratio <- ((freq_ab/n)-prod(freq/n))^2/(prod(freq/n)*prod(1-(freq/n)))
  r2_simplify
}

#' GET_MEAN_R2
#' compute the mean r2 for LD of sites within a given window
#' INPUT
#'  @param ms_matrix_pop matrix with 2 columns corresponding to the haplotypes of a pair of sites
#'  @param index_window  vector with the index of the window of each site
#'  @param j number of the index of the window we are considering
#'  RETURN
#'  @return mean r2 for all pairwise comparisons within a window
get_meanr2 <- function(ms_matrix_pop,index_window,j) {
  snps_in_window <- which(index_window == j)
  if(length(snps_in_window)>1) {
    # based on the sites that belong to a given window define each pair of sites
    index_pairs <- combn(snps_in_window,2)
    # get the mean_r2 for each window
    res <- mean(apply(index_pairs, 2, function(col){r2linkage(ms_matrix_pop[,col])}))  
  } else {
    res <- NA
  }
  res
}
