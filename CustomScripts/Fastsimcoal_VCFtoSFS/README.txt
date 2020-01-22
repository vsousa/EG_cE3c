# Vitor Sousa
# 20.01.2020 -----------------

To go from a VCF to SFS you can modify the settings of ProcessVCF.sh.
This will divide your genome into "blocks" of a given length.
For each block, it resamples a given number of individuals
at SNPs without missing data.

# REQUIRED FILES:
# - vcffile
#   	 VCF file after applying your filters. NOTE: For SFS analyses you should
#        have a VCF filtered based on depth of coverage (DP>10), 
#        but also discarding regions with very high coverage and 
#        with excess heterozygosity across all individuals,
#        as those could be due to mapping errors in repetitive or structural variants, 
#        including duplicated regions, discarding sites with low mapping quality.
#        NOTE: you should not have a MAF filter, as SFS is very sensitive
#        to rare variants!
#        
# - indpopinfo 
#        Text file with two columns with individual ID and corresponding POP
#        NOTE: indpopinfo file defines the order of individuals and populations.
#              Thus, this file can have just a sub-set of individuals from the VCF,
#              for which you want to obtain the SFS.
