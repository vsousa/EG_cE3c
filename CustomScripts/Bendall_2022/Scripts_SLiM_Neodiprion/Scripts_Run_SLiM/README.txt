# SLiM input files and bash scripts to run the simulations
# for different combination of parameters
# Summary of demographic estimates and re-scaling of parameters to run with SLiM
# Author: Vitor Sousa
# Last updated: 07/02/2022

# In the models we assumed that:
# - N. pinetum corresponds to pop1, where the beneficial mutation allele a (with initial frequency q0) is favoured
# - N. lecontei corresponds to pop2, where the allele A is favoured (with initial frequency 1-q0) 

# This folder contains the following files:
# a) input files for SLiM
# these are templates with TAGs for relevant parameters. 
# Hence, by replacing these TAGs by specific values, 
# these template files can be used for simulations with different parameter combinations
# - sim_createAncestral_template.s : template to simulate the ancestral population
# - sim_divsel3_scaled_template.s : template to simulate IM model with two diverging populations with gene flow,
#                                   experiencing divergent selection due to new mutation
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
#    specify the template file (either sim_divsel3_xaut_newmut_template.s or sim_divsel3_xaut_template.s)
#    use the same combination of parameters as in CreateAncestralPops.sh
./SimFromAncestralPops.sh

# c) File with the demographic estimates and re-scaling of parameters
# FILE: DemographicHistory_ParamEstimates_SLIM_rescaling.xlsx

# d) Script to compute scaling of recombination rate based on 
#    estimates of Neodiprion lecontei crosses
# FILE: haldane_map.r