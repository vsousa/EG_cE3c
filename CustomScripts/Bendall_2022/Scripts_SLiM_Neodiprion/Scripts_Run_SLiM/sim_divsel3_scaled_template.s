// IM model with two populations under divergent selection with migration,
// with different effective sizes and migration rates
// template
initialize()
{
	// DEFINE PARAMETERS
	// mutation rate
	murate = MUTRATE;
	// selective coefficient of mutation under divergent selection
	smutben = SELMUTBEN;
	// selection coefficient of deleterious mutations
	smutdel = SELMUTDEL;
	// sex-ratio defined as sr=Nm/(Nm+Nf) where Nm is the pop size of males and Nf is the pop size of females
	sr = SEXRATIO;
	// chromosome (A - autosome/diploid, X - x-chromosome/hemizygous)
	XorA = CHR;
	// Population size of pop1
	pop1size = POP1SIZE;
	// Population size of pop2
	pop2size = POP2SIZE;
	// Migration rate m=m12
	migrate = MIGRATE;
	// Recombination rate
	recrate = RECRATE;
	
	// SCALE PARAMETERS
	// scale effective size and recombination rate
	// such that the expected diversity and LD patterns are the same
	// Ratio of NX/NA depends on the sex ratio sr=Nm/(Nm+Nf)
	// NX/NA=(9/8)(1/(1+sr)), according to Wright (1969) Evolution and Genetics of Populations. Volume 2 The theory of gene frequencies. Equations 8.10 and 8.12.
	if(XorA == "A") {
		// ratio of effective sizes NX/NA
		ratioNxNa = (9/8)*(1/(1+sr));
		// obtain the population size for a diploid population that corresponds to the Ne of the X-chromosome
		// pop1
		pop1size = asInteger(round(pop1size*ratioNxNa));
		// pop2
		pop2size = asInteger(round(pop2size*ratioNxNa));
		// obtain the recombination rate, assuming that only females recombine at X-chromosome
		recrate = recrate*(2/3);
	}
	
	// DEFINE CONSTANT PARAMETERS
	// these can be used on other parts of the SLiM code
	// Effective population size, after re-scaling
	defineConstant("Npop1", pop1size);
	defineConstant("Npop2", pop2size);
	// Migration rate after re-scaling 
	// Note that the migration rate does not need to be re-scaled as since we change NX and NA the NXm=NAm.
	defineConstant("MigRate", migrate);
	// Selective coefficient of mutation under divergent selection
	defineConstant("SmutBen", smutben);
	// print the constant parameters
	cat(format("Npop1=%i\t", Npop1));
	cat(format("Npop2=%i\t", Npop2));	
	cat(format("MigRate=%6.5f\t", MigRate));
	cat(format("SmutBen=%6.5f\n", SmutBen));

	// INITIALIZE PARAMETERS
	// set the overall mutation rate
	initializeMutationRate(murate);
	
	// m1 mutation type: neutral
	initializeMutationType("m1", 0.5, "f", 0.0);
        m1.convertToSubstitution = T;
	// prevent mutations of different groups to stack
	m1.mutationStackGroup = -1;
        m1.mutationStackPolicy = "l";

	// m2 mutation type: beneficial in pop 1 (allele a at locus under divergent selection)
	initializeMutationType("m2", DOMINANCEBEN, "f", smutben);
        m2.convertToSubstitution = F;
        m2.mutationStackPolicy = "l";
       	m2.mutationStackGroup = -1;
	
	// m3 mutation type: beneficial in pop 2 (allele A at locus under divergent selection)
	// parallel dominance model
	dominanceP2 = 1-DOMINANCEBEN;	
	initializeMutationType("m3", dominanceP2, "f", smutben);
        m3.convertToSubstitution = F;
        m3.mutationStackPolicy = "l";
       	m3.mutationStackGroup = -1;

	// g1 genomic element type: neutral m1 mutations for all sites
	initializeGenomicElementType("g1", m1, 1.0);

	// chromosome of length SEQSIZE in bp
	initializeGenomicElement(g1, 0, SEQSIZE);

	// uniform recombination along the chromosome
	initializeRecombinationRate(recrate);
	
        // specify that we are looking at X chromosome
        initializeSex(XorA, 1.0);
}

// Read the state of the ancestral population simulated for 4N generations
1 late() {
	// Name of file with ancestral slim pop info
	ancfilename_tag = ANCFILE;
	ancfilename = "../FOLDERANC/" + ancfilename_tag + CHR + "_" + rep + "_ancestralPop.slim";
	cat("FileName=" + ancfilename);
	// Read data from file
	sim.readFromPopulationFile(ancfilename);
	// start a newly seeded run by incrementing the previous seed
	setSeed(getSeed());
	sim.generation=1;
	cat("\nLoaded ancestral population file!\n");
}

// Simulate the population split and origin of another population
2 { 
    // population split with a given sex ratio
    sim.addSubpopSplit("p2", Npop2, p1, SEXRATIO); 
    // Set new size for population p1
    p1.setSubpopulationSize(Npop1);

}

// Add divergent selected locus to pop1 and pop2
// with a given initial frequency, assumed to be the same in both pops
3 { 

	// Add a mutation under divergent selection with initial frequency FREQ
	// select randomly a freq from each pop and add the mutation under divergent selection
	// mutation of type m2 corresponds to allele a
	// mutation of type m3 corresponds to allele A
	
	// COMPUTE THE INITIAL FREQUENCY
	// given the different number of gene copies in A and X, need to compute the number
	// of genomes that need to sample to have the same initial frequency in X and A.
	// Get the number of gene copies in haplodiploids
	// Ratio of NX/NA depends on the sex ratio sr=Nm/(Nm+Nf)
	// NX/NA=(9/8)(1/(1+sr)), according to Wright (1969) Evolution and Genetics of Populations. Volume 2 The theory of gene frequencies. Equations 8.10 and 8.12.
	if(CHR == "X") 
	{
		ratioNxNa = (9/8)*(1/(1+SEXRATIO));
		// numgenes is the number of chromosomes in each population
		numgenesPop1 = asInteger(round(2.0*Npop1*ratioNxNa));
		numgenesPop2 = asInteger(round(2.0*Npop2*ratioNxNa));
	}
	if(CHR == "A") 
	{
		// numgenes is the number of chromosomes in each population
		numgenesPop1 = asInteger(round(2.0*Npop1));
		numgenesPop2 = asInteger(round(2.0*Npop2));
		
	}
	// compute the initial number of genomes we need to mutate
	freqThresholdPop1 = asInteger(round(FREQ*numgenesPop1)); // freqThresholdPop1 is the number of genomes we will mutate in Pop1
	freqThresholdPop2 = asInteger(round(FREQ*numgenesPop2)); // freqThresholdPop2 is the number of genomes we will mutate in Pop2

	// POP1
	// pick genomes randomly to which allele a is added, 
	// allele A is added to all remaining genomes 
	// get all the genomes from pop1
	gp1 = p1.genomes;
	aux_p1 = !gp1.isNullGenome; // index of genomes that are not null
	tmp_p1 = c(0:((p1.individualCount*2)-1)); // vector of indices
	tmp_m2_sg1 = sample(tmp_p1[aux_p1], freqThresholdPop1, F); // sample the index of freqThreshold genomes
	diff_sg1=setDifference(tmp_p1[aux_p1],tmp_m2_sg1); // find the genomes without the selected mutation
	// Add m2 mutations to randomly selected genomes of pop1
	sg1_m2 = gp1[tmp_m2_sg1]; // get the genomes that will get a m2 mut in pop1
	mutm2 = sg1_m2.addNewDrawnMutation(m2, POSMUT);  // add m2 mutation at selected position 
	// Add m3 mutations to the other genomes of pop1
	sg1_m3 = gp1[diff_sg1]; // get the genomes that will get a m3 mut in pop1
	mutm3 = sg1_m3.addNewDrawnMutation(m3, POSMUT);  // add m3 mutation at selected position 
	
	// POP2
	// repeat for pop 2
	gp2 = p2.genomes; // get p2 genomes
	aux_p2 = !gp2.isNullGenome; // index of genomes that are not null
	tmp_p2 = c(0:((p2.individualCount*2)-1)); // vector of indices
	tmp_m2_sg2 = sample(tmp_p2[aux_p2], freqThresholdPop2, F); // sample the index of freqThreshold genomes
	diff_sg2=setDifference(tmp_p2[aux_p2],tmp_m2_sg2); // find the genomes without the selected mutation
	// Add m2 mutations to randomly selected genomes of pop2
	sg2_m2 = gp2[tmp_m2_sg2]; // get the genomes that will get a m2 mut in pop2
	sg2_m2.addMutations(mutm2);  // add m2 mutation at selected position 	
	// Add m3 mutations to randomly selected genomes of pop2
	sg2_m3 = gp2[diff_sg2]; // get the genomes that will get a m3 mut in pop2
	sg2_m3.addMutations(mutm3);  // add m3 mutation at selected position 
}


// FITNESS LANDSCAPE UNDER PARALLEL DOMINANCE
// Set the selective coefficient of m2 muts as -s
// All muts of type m2 are beneficial in p1 and neutral in p2
3:ENDGEN fitness(m2,p2) {
   return 1.0;
}
// All muts of type m3 are beneficial in p2 and neutral in p1
3:ENDGEN fitness(m3,p1) {
   return 1.0;
}


// Populations with Ndesc effective size experience migration until ENDGEN generation
3:ENDGEN { 
    // Set new size for population p1
    p1.setSubpopulationSize(Npop1); 
    p2.setSubpopulationSize(Npop2); 
    // set migration rate
    // migration rate is defined as the ratio of the migration rate from Pin to Lec and Vice-versa
    // The fastsimcoal2 estimates backwards in time indicate a ratio of m10/m01=0.0468
    MigRateIntoPop1=MigRate;
    MigRateIntoPop2=MigRate*RATIOM01M10;
    p1.setMigrationRates(c(p2), c(MigRateIntoPop1));     
    p2.setMigrationRates(c(p1), c(MigRateIntoPop2));
    // Every 100 generations output info about the beneficial mutations
    if((sim.generation % 100) == 0) {
	    filenamems = TAGFILEEND;
	    benmuts_m2=sim.mutationsOfType(m2);
	    if(size(benmuts_m2)>0) 
	    {
	    	freqbenp1_m2 = sim.mutationFrequencies(p1, benmuts_m2); // NULL to get frequency from all pops    
	    	freqbenp2_m2 = sim.mutationFrequencies(p2, benmuts_m2); // NULL to get frequency from all pops    
	    }
	    else
	    {
	    	freqbenp1_m2 = 0;
	    	freqbenp2_m2 = 0;
	    }
	    benmuts_m3=sim.mutationsOfType(m3);
	    if(size(benmuts_m3)>0) 
	    {
	    	freqbenp1_m3 = sim.mutationFrequencies(p1, benmuts_m3); // NULL to get frequency from all pops    
	    	freqbenp2_m3 = sim.mutationFrequencies(p2, benmuts_m3); // NULL to get frequency from all pops    
	    }
	    else
	    {
	    	freqbenp1_m3 = 0;
	    	freqbenp2_m3 = 0;
	    }
	    // Write to a file
	    writeFile(filenamems + "_" + rep + ".benmuts", rep + " " + sim.generation + " " + freqbenp1_m2 + " " + freqbenp2_m2 + " " + freqbenp1_m3 + " " + freqbenp2_m3 , append=T);

	    // MEAN FITNESS
	    // Compute the mean fitness of each population
	    fitp1 = mean(p1.cachedFitness(NULL));
	    fitp2 = mean(p2.cachedFitness(NULL));
	    writeFile(filenamems + "_" + rep + ".meanfit", rep + " " + sim.generation + " " + fitp1 + " " + fitp2, append=T);
    }

}

// CHECK fitness landscape
// At generation 5 print the genotypes and corresponding fitness of all individuals
// this is used to check if the fitness landscape is correctly defined
5 {
    // Define output file name tag
    filenamems = TAGFILEEND;
    // Ouput the fitness of individuals of pop1 and pop2
    p1fit=p1.cachedFitness(NULL);
    p2fit=p2.cachedFitness(NULL);

    // Ouput the genotype of each individual
    for(i in 0:(p1.individualCount-1)) 
    {
	// write a file with the following columns for each individual (each individual corresponds to a row)	
	// generation sex genotype_number_a_copies genotype_number_A_copies fitness
	// note that the genotype is coded by the number of copies of a and A alleles:
	// homozygote for a is 2 0, 
	// heterozygote is 1 1,  
	// homozygote for A is 0 2
    	writeFile(filenamems + "_" + rep + ".fitlandscape", sim.generation + " 1" + " " + p1.individuals[i].sex + " " + p1.individuals[i].countOfMutationsOfType(m2) + " " + p1.individuals[i].countOfMutationsOfType(m3) + " " + p1fit[i], append=T);
    }
    
    // Get the genotype and fitness for each individual of pop2    
    for(i in 0:(p2.individualCount-1)) 
    {
    	// write a file with the following columns for each individual (each individual corresponds to a row)	
    	// generation sex genotype_number_a_copies genotype_number_A_copies fitness
    	// note that the genotype is coded by the number of copies of a and A alleles:
    	// homozygote for a is 2 0, 
    	// heterozygote is 1 1,  
	// homozygote for A is 0 2
    	writeFile(filenamems + "_" + rep + ".fitlandscape", sim.generation + " 2" + " " + p2.individuals[i].sex + " " + p2.individuals[i].countOfMutationsOfType(m2)+ " " + p2.individuals[i].countOfMutationsOfType(m3) + " " + p2fit[i], append=T);
    }
}

// sample 20 genomes from each pop at SAVEGEN
SAVEGEN late()
{
	// Output ms file with sample from populations
	filenamems = TAGFILESAVE;

	// get a sample of genomes from p1 are not null (only X chromosome sampled)
	gp1 = p1.genomes;
	gp2 = p2.genomes;

	// get number of individuals and of firt male index (genomes are sorted)
	gp1male = gp1[0:((p1.individualCount*2)-1)]; 
	// sample 8 genomes from pop1
        sg1 = sample(gp1male[!gp1male.isNullGenome], 8, F);        
        // repeat for pop2
	gp2male = gp2[0:((p2.individualCount*2)-1)]; // sample all
        // sample 12 genomes from pop2
        sg2 = sample(gp2male[!gp2male.isNullGenome], 12, F);
        // Print ms output by putting the two genome vectors into a single vector
        sg = c(sg1, sg2); 
        sg.outputMS("ms_" + filenamems + "_" + rep + ".txt", append=F);

}

// run to last generation
ENDGEN late()
{
	// Output ms file with sample from populations
	filenamems = TAGFILEEND;

	// get a sample of genomes from p1 are not null (only X chromosome sampled)
	gp1 = p1.genomes;
	gp2 = p2.genomes;

	// get number of individuals and of first male index (genomes are sorted)
	gp1male = gp1[0:((p1.individualCount*2)-1)]; 
	// sample 20 genomes from pop1
        sg1 = sample(gp1male[!gp1male.isNullGenome], 8, F);        
        // repeat for pop2
	gp2male = gp2[0:((p2.individualCount*2)-1)]; // sample all
        // sample 20 genomes from pop2
        sg2 = sample(gp2male[!gp2male.isNullGenome], 12, F);
        // Print ms output by putting the two genome vectors into a single vector
        sg = c(sg1, sg2); 
        sg.outputMS("ms_" + filenamems + "_" + rep + ".txt", append=F);

}


