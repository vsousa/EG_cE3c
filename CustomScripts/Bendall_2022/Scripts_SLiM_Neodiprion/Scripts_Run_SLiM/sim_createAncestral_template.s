// set up a simple neutral simulation
initialize()
{
	// convert variables
	murate = MUTRATE;
	sr = SEXRATIO;
	XorA = CHR;
	popsize = POPSIZEANC;
	recrate = RECRATE;
	
	// Ratio of NX/NA depends on the sex ratio sr=Nm/(Nm+Nf)
	// NX/NA=(9/8)(1/(1+sr)), according to Mendez (2017) TPB ("Differences in the effective population sizes of males and females do not require differences in their distribution of offspring number")
	// and according to Wright (1969) Evolution and Genetics of Populations. Volume 2 The theory of gene frequencies. Equations 8.10 and 8.12.
	if(XorA == "A") {
		ratioNxNa = (9/8)*(1/(1+sr));
		popsize = asInteger(round(popsize*ratioNxNa));
		recrate = recrate*(2/3);
	}
	defineConstant("Npop", popsize);
	defineConstant("Chr", CHR);

	cat(format("Npop=%i\t", Npop));
	cat("Chr=" + Chr);

	// set the overall mutation rate
	initializeMutationRate(murate);
	// m1 mutation type: neutral
	initializeMutationType("m1", 0.5, "f", 0.0);
        m1.convertToSubstitution = T;
	// we do not want mutations of different groups to stack
	m1.mutationStackGroup = -1;
        m1.mutationStackPolicy = "l";

	// g1 genomic element type: uses m1 for all mutations (neutral)
	initializeGenomicElementType("g1", m1, 1.0);

	// chromosome of length SEQSIZE
	initializeGenomicElement(g1, 0, SEQSIZE);

	// uniform recombination along the chromosome
	initializeRecombinationRate(recrate);
	
        // specify that we are looking at X chromosome
        initializeSex(XorA, 1.0);
}

// create a population
1
{
        // create a sub-population with sex ratio SEXRATIO (NOTE: this is defined in terms of males)
	sim.addSubpop("p1", Npop, SEXRATIO); // tag, Ne, sex-ratio
}

// run to generation ENDGEN
ENDGEN late()
{
	filename = TAGFILE;
	// save the final state of the simulation
	sim.outputFull(filename + Chr + "_" + rep + "_ancestralPop.slim");
}


