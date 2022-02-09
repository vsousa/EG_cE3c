// simulate ancestral population and save its state after ENDGEN generations
initialize()
{
	// DEFINE PARAMETERS
	// mutation rate
	murate = MUTRATE; 
	// sex-ratio
	sr = SEXRATIO;
	// chromosome (A - autosome/diploid, X - x-chromosome/hemizygous)
	XorA = CHR;
	// effective population size (Ne)
	popsize = POPSIZEANC;
	// recombination rate
	recrate = RECRATE;
	
	// SCALE PARAMETERS
	// scale effective size and recombination rate
	// such that the expected diversity and LD patterns are the same
	// Ratio of NX/NA depends on the sex ratio sr=Nm/(Nm+Nf)
	// NX/NA=(9/8)(1/(1+sr))  according to
	// Wright (1969) Evolution and Genetics of Populations. Volume 2 The theory of gene frequencies. Equations 8.10 and 8.12.
	// Mendez (2017) TPB ("Differences in the effective population sizes of males and females do not require differences in their distribution of offspring number")
	// Clemente F, Gautier M, Vitalis R (2018) Inferring sex-specific demographic history from SNP data. PLoS Genet 14(1): e1007191. https:// doi.org/10.1371/journal.pgen.1007191
	if(XorA == "A") {
		// ratio of effective sizes NX/NA
		ratioNxNa = (9/8)*(1/(1+sr));
		// obtain the population size for a diploid population that corresponds to the Ne of the X-chromosome
		popsize = asInteger(round(popsize*ratioNxNa));
		// obtain the recombination rate, assuming that only females recombine at X-chromosome
		recrate = recrate*(2/3);
	}
	
	// DEFINE CONSTANT PARAMETERS
	// Effective population size, after re-scaling
	defineConstant("Npop", popsize);
	// Chromosome
	defineConstant("Chr", CHR);

	// print the effective size and chromosome
	cat(format("Npop=%i\t", Npop));
	cat("Chr=" + Chr);

	// INITIALIZE PARAMETERS
	// set the overall mutation rate
	initializeMutationRate(murate);
	// m1 mutation type: neutral
	initializeMutationType("m1", 0.5, "f", 0.0);
	// do not convert to substitutions (fixed derived alleles are recorded)
        m1.convertToSubstitution = F;
	// avoid mutations of different groups to stack
	m1.mutationStackGroup = -1;
        m1.mutationStackPolicy = "l";

	// g1 genomic element type: uses m1 for all mutations (neutral)
	initializeGenomicElementType("g1", m1, 1.0);

	// initialize chromosome of length 500 kb of neutral mutations
	initializeGenomicElement(g1, 0, 499999);

	// initialize uniform recombination along the chromosome
	initializeRecombinationRate(recrate);
	
        // specify that we are looking at A or X chromosome
        initializeSex(XorA, 1.0);
}

// create a population of Npop individuals with a given sexratio
1
{
        // create a sub-population with sex ration SEXRATIO (NOTE: this is defined in terms of males sr=Nm/(Nm+Nf)
	sim.addSubpop("p1", Npop, SEXRATIO); 
}

// run to generation ENDGEN generation
ENDGEN late()
{
	// save output into file with filename
	filename = TAGFILE;
	// save the final state of the simulation
	sim.outputFull(filename + Chr + "_" + rep + "_ancestralPop.slim");
}


