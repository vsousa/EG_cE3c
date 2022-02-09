# Scripts to run the analysis using SLiM version 3.3
# SLiM version 3.3, built Aug  5 2019 18:29:03
# Author: Vitor Sousa
# Last updated: 07/02/2022

# There is one folder with the template for SLiM input files for
# - scaled: in this case the effective size and recombination rates are adjusted, such that both 
#           the hemizygous (haplodiploid) and diploid chromosomes have the same average 4Nu and 4Nr, 
#           where N is the effective size, u is the mutation rate and r is the recombination rate.

# within that folder you can find the SLiM template input files, ending with *.s
# and the bash scripts used to modify the template files and launch SLiM.

# MODEL
# We assume an isolation with migration model where an ancestral population of size 2N=1500
# which is equivalent to the number of gene copies in a population of X-chromosomes with 
# 1000 individuals (500 females with two copies each, and 500 males with one copy each).
# This ancestral population evolves for Tanc=10,000 generations (i.e. Tanc>4*Ne) to ensure that
# it reaches mutation-drift equilibrium.
# We simulate chromosomes with 500Kb with only neutral mutations.
# The simulation of the ancestral population is done with a particular SLiM input file and corresponding bash scripts.
# At the end of the Tanc generations, the state of the population is saved.

# Then, to simulate two populations diverging, we consider that the ancestral population
# splits into two descending populations of size 2N=1500 that experience a constant and
# symmetric migration rate m=m12=m21, where m12 is the migration rate from pop1 to pop2, and
# m21 is the migration rate from pop2 to pop1.
# The state of the ancestral population is loaded into memory
# and we simulate chromosomes with 500Kb with neutral mutations, and with a site
# under divergent selection at position 250Kb.
# At that site, we follow the trajectory of allele a with initial frequency q0.
# We consider a parallel dominance model, where allele a is favoured in pop1 with 
# selective coefficient 1+s in homozygotes (and male hemizygous), 
# and allele A is favoured in pop2 with selective coefficient 1+s in homozygotes (and male hemizygous).
# Different input template files for SLiM were used for the case of a new mutation with q0=1/(2N) and other values.


