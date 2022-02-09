# SLiM input files and bash scripts to run the simulations
# for different combination of parameters
# Author: Vitor Sousa
# Last updated: 07/02/2022

# This folder contains the following files:
# a) input files for SLiM
# these are templates with TAGs for relevant parameters. 
# Hence, by replacing these TAGs by specific values, 
# these template files can be used for simulations with different parameter combinations
# - sim_createAncestral_template.s : template to simulate the ancestral population
# - sim_divsel3_scaled_new_template.s : template to simulate IM model with two diverging populations with gene flow,
#                                       experiencing divergent selection due to new mutation
# - sim_divsel3_scaled_template.s  : template to simulate IM model with two diverging populations with gene flow,
#                                    experiencing divergent selection due to mutations with an initial frequency FREQ
#                                    mimicking standing genetic variation
#
# b) bash files to run the simulations
# The combination of parameters and the template files used are defined in these files.
# - CreateAncestralPops.sh : script to simulate the ancestral populations
# - SimFromAncestralPops.sh : script to simulate the diverging populations pop1 and pop2 

# To run the simulations follow the steps 
# 1. Create ancestral population
#    define combination of parameters in CreateAncestralPops.sh
./CreateAncestralPops.sh

# 2. Simulate the divergence of populations
#    specify the template file (either sim_divsel3_scaled_new_template.s or sim_divsel3_scaled_template.s)
#    use the same combination of parameters as in CreateAncestralPops.sh
./SimFromAncestralPops.sh