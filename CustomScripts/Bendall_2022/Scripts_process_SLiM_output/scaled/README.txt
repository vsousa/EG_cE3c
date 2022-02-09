# PROCESS the output files from SLiM
# the following scripts are used to process the multiple files created for each SLiM simulations
# Author: Vitor Sousa
# Last updated: 07/02/2022

# RSCRIPTS folder
# contains the file with definition of functions to compute summary statistics
# computeStats.r

# 1. compute statistics for linked selection based on samples of 20 genomes of 500kb from each pop
# average across the chromosome and window-based analysis
# process ms files: read_ms_createSummaryFiles.r
# requires the script: ./Rscripts/computeStats.r
# INPUT: ms files from SLiM runs
# OUTPUT: creates 2 folders with the summary of results across 1000 runs
# ./TAG/sumstat/    : summary statistics for all simulations (irrespective of keeping the beneficial allele)
# ./TAG/sumstat_cb/ : summary statistics conditional on retainining the beneficial a allele in pop1

# 2. compute statistics for selected site from simulations of chromosomes with 500kb
# process trajectory files: read_traj_createSummaryFiles.r
# requires the script: ./Rscripts/computeStats.r
# INPUT: trajectory files from SLiM runs
# OUTPUT: creates 2 summary files with the trajectores and a matrix with means across sims


# 3. further process statistics of ms to obtain means and quantiles of pi, dxy and fst
# process the summary files of ms to create simple tables that can be used for plots.
# process_ms_summaryFiles.r
# INPUT: summary files output from step 1 above
# OUTPUT: creates 8 files with the summary of pi for pop1, pi for pop2, dxy, fst for cases with all runs ("sumstat") and 
#         cases conditional on keeping allele a in pop1 ("sumstat_cb")


# 4. further process statistics of ms to obtain means and quantiles of pi and fst along the chromosome
#    on windows base.
#    process_ms_summaryFiles_windowscan.r
# INPUT: summary files output from step 1 above
# OUTPUT: creates a RDS file with the mean and quantiles 0.05 0.25 0.75 0.95 for 
#         pi_pop1, pi_pop2, dxy, fst across 1000 simulations, for each window of 20Kb


# 5. read ms files to compute LD r^2 statistic
# read_ms_compute_LD.r
# INPUT: ms files from SLiM runs
# OUTPUT: creates 1 folders with the mean_r2 at xKb windows for each population across 1000 runs
# ./TAG/sumstat/    : summary statistics for all simulations (irrespective of keeping the beneficial allele)


# 6. process the LD statistics files
# process_ld.r
# INPUT:  summary files output from step 5 above
# OUTPUT: summary table with the mean and quantiles 0.05 0.25 0.75 0.95 for 
#         r2 for pop1 across 1000 simulations, across the 25 windows of 20Kb


