# Haldane mapping function
# haldane map function
# cm to recombination probability
haldane_cm_to_r <- function(cm) {
  (1-exp(-(2*cm)/100))/2
}


# haldane map function
# prob recombination to cm
haldane_r_to_cm <- function(prob) {
  50*log(1/(1-(2*prob)))
}

# convert recombination cM into probability for SLiM simulations
cmsawflies <- 3.43 #cM/Mb

# prob of recombination per Mb
probr <- haldane_cm_to_r(cmsawflies)
probr
# probability of r per site
probr_site <- probr/1e6
probr_site

# rec to use in slim assuming males do not recombine
# multiply by (3/2)
probr_site_slim <- probr_site*(3/2)

# total probability of r in 10Mb
nsites <- 1e7
probr_chr <- probr_site*nsites

# relationship between mutation in recombination
# assuming a mutrate of 1.75e-8 
# (which would correspond to the same theta=4Nanc*u*L as obtained with fastsimcoal2)
mutsite <- 1.7e-8
probr_site_slim/mutsite

# which was round to 3
(probr_site_slim*(2/3)*1e6)


# rec rate in simulations with recrate=mutrate
recpersite <- 2.5e-7*(2/3) # recombination per site
nsites <- 2.5e5
rec_chr <-recpersite*nsites  # recombination rate across chromosome
# convert to cM
haldane_r_to_cm(rec_chr)

# rec rate in simulations with recrate=0.1*mutrate
recpersite <- 2.5e-8*(2/3) # recombination per site
nsites <- 2.5e5
rec_chr <-recpersite*nsites  # recombination rate across chromosome
# convert to cM
haldane_r_to_cm(rec_chr)

# rec rate in simulations for sawflies
recpersite <- 1.05e-6*(2/3) # recombination per site
nsites <- 2.5e5
rec_chr <-recpersite*nsites  # recombination rate across chromosome

# convert to cM
haldane_r_to_cm(rec_chr)
rec_chr









