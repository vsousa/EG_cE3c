# PROCESS the output files from SLiM
# the following scripts are used to process the multiple files created for each SLiM simulations
# according to the demographic history of Neodiprion pinetum and N. lecontei
# inferred with fastsimcoal2
# Author: Vitor Sousa
# Last updated: 07/02/2022

# 1. compute statistics for linked selection based on samples of 8 and 12 genomes of 500kb from each pop
# average across the chromosome and window-based analysis
# process ms files: readms_createSummaryFiles_sawflies.r
# requires the script: computeStats_SFS.r
# INPUT: ms files from SLiM runs
# OUTPUT: creates 2 folders with the summary of results across 1000 runs
# ./TAG/sumstat/    : summary statistics for all simulations (irrespective of keeping the beneficial allele)
# ./TAG/sumstat_cb/ : summary statistics conditional on retainining the beneficial a allele in pop1

# 2. process statistics of ms to obtain means and quantiles of pi, dxy and fst
#    obtain means and quantiles of pi and fst along the chromosome
#    and obtain neutral SFS
# process the summary files of ms to create simple tables that can be used for plots later.
# process_ms_summaryFiles_neutralSFS_windows.r
# INPUT: summary files output from step 1 above
# OUTPUT: a) creates 8 files with the summary of pi for pop1, pi for pop2, dxy, fst for cases with all runs ("sumstat") and 
#            cases conditional on keeping allele a in pop1 ("sumstat_cb")
#         b) creates a RDS file with the mean and quantiles 0.05 0.25 0.75 0.95 for 
#            pi_pop1, pi_pop2, dxy, fst across 1000 simulations, for each window of 20Kb
#         c) creates a RDF file with the neutral SFS for diploid case and neutral SFS for haplodiploid case,
#            summing over the 1000 simulations for all the neutral combinations

# 3 - Required file with definition of functions
# computeStats_SFS.r



