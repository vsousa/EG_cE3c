# Vitor Sousa
# 20.01.2020 -----------------

To analyse the output of fastsimcoal2 runs, use the r script
- ./Scripts/Analyse_Fsc.r

This R script will call several functions that will collect the runs
of fastsimcoal2. It will:
- output the Maximum Likelihood Estimated parameters, from the best run.
- plot graphically the parameter estimates
- plot the fit of the model to the observed SFS
- re-compute the likelihood and AIC based on an observed SFS just with independent SNPs

The functions are defined in the files:
- ParFileInterpreter_VS.r : functions to visually show the estimates
- utilFscOutput.r : several functions to read, process and plot fastsimcoal2 results
