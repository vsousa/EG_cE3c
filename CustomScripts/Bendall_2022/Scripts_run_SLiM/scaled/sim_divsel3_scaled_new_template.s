// IM model with two populations under divergent selection with migration
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
	// Population size of ancestral population
	popsize = POPSIZEANC;
	// Population size at time of split into two pops
	nbot = POPSIZEBOT;
	// Population size of descendent pops pop1 and pop2
	ndesc = POPSIZEDESC;
	// Migration rate m=m12=m21
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
		popsize = asInteger(round(popsize*ratioNxNa));
		nbot = 	asInteger(round(nbot*ratioNxNa));
		ndesc = asInteger(round(ndesc*ratioNxNa));
		// obtain the recombination rate, assuming that only females recombine at X-chromosome
		recrate = recrate*(2/3);
	}
	
	// DEFINE CONSTANT PARAMETERS
	// these can be used on other parts of the SLiM code
	// Effective population size, after re-scaling
	defineConstant("Npop", popsize);
	defineConstant("Nbot", nbot);
	defineConstant("Ndesc", ndesc);
	// Migration rate after re-scaling 
	// Note that the migration rate does not need to be re-scaled as since we change NX and NA the NXm=NAm.
	defineConstant("MigRate", migrate);
	// Selective coefficient of mutation under divergent selection
	defineConstant("SmutBen", smutben);
	// print the constant parameters
	cat(format("Npop=%i\t", Npop));
	cat(format("Nbot=%i\t", Nbot));
	cat(format("Ndesc=%i\t", Ndesc));
	cat(format("MigRate=%6.5f\t", MigRate));
	cat(format("SmutBen=%6.5f\n", SmutBen));

	// INITIALIZE PARAMETERS
	// set the overall mutation rate
	initializeMutationRate(murate);
	// m1 mutation type: neutral
	initializeMutationType("m1", 0.5, "f", 0.0);
        m1.convertToSubstitution = F;
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

	// chromosome of length 500 kb of neutral mutations
	initializeGenomicElement(g1, 0, 499999);

	// uniform recombination along the chromosome
	initializeRecombinationRate(recrate);
	
        // specify that we are looking at A or X chromosome
        initializeSex(XorA, 1.0);
}

// initialize population reading the state of the ancestral population 
1 late() {
	// Name of file with ancestral slim pop info
	ancfilename_tag = ANCFILE;
	ancfilename = "../Ancestral_500Kb/" + ancfilename_tag + CHR + "_" + rep + "_ancestralPop.slim";
	cat("FileName=" + ancfilename);
	// Read data from file
	sim.readFromPopulationFile(ancfilename);
	// start a newly seeded run by incrementing the previous seed
	setSeed(getSeed());
	sim.generation=1;
	cat("\nLoaded ancestral population file!\n");
}

// Simulate the fouding of population 2 (split time)
2 { 
    // create population 2 with size Nbot and sex-ratio SEXRATIO
    sim.addSubpopSplit("p2", Nbot, p1, SEXRATIO); 
    // set size for population p1 to Nbot
    p1.setSubpopulationSize(Nbot);

}

// Add divergent selected locus to pop1 and pop2
// with a given initial frequency, assumed to be the same in both pops
3 { 
	// Add a NEW mutation with initial frequency 1/(2N) under divergent selection with initial frequency FREQ
	// select randomly a freq from each pop and add the mutation under divergent selection
	// mutation of type m2 corresponds to allele a
	// mutation of type m3 corresponds to allele A
	
	// POP1
	// pick genomes randomly to which allele a is added, 
	// allele A is added to all remaining genomes 
	// get all the genomes from pop1
	gp1 = p1.genomes;
	aux_p1 = !gp1.isNullGenome; // index of genomes that are not null
	tmp_p1 = c(0:((p1.individualCount*2)-1)); // vector of indices
	tmp_m2_sg1 = sample(tmp_p1[aux_p1], 1, F); // sample the index of the genome where the mutation is added
	diff_sg1=setDifference(tmp_p1[aux_p1],tmp_m2_sg1); // find the genomes without the selected mutation
	// Add m2 mutations to randomly selected genome of pop1
	sg1_m2 = gp1[tmp_m2_sg1]; // get the genomes that will get a m2 mut in pop1
	mutm2 = sg1_m2.addNewDrawnMutation(m2, 249999);  // add m2 mutation at selected position 
	// Add m3 mutations to the other genomes of pop1
	sg1_m3 = gp1[diff_sg1]; // get the genomes that will get a m3 mut in pop1
	mutm3 = sg1_m3.addNewDrawnMutation(m3, 249999);  // add m3 mutation at selected position 
	// POP2
	// repeat for pop 2
	gp2 = p2.genomes; // get p2 genomes
	aux_p2 = !gp2.isNullGenome; // index of genomes that are not null
	tmp_p2 = c(0:((p2.individualCount*2)-1)); // vector of indices
	tmp_m2_sg2 = sample(tmp_p2[aux_p2], 1, F); // sample the index of the genome where the mutation is added
	diff_sg2=setDifference(tmp_p2[aux_p2],tmp_m2_sg2); // find the genomes without the selected mutation
	// Add m2 mutations to randomly selected genomes of pop2
	sg2_m2 = gp2[tmp_m2_sg2]; // get the genomes that will get a m2 mut in pop2
	sg2_m2.addMutations(mutm2);  // add m2 mutation at selected position 	
	// Add m3 mutations to randomly selected genomes of pop2
	sg2_m3 = gp2[diff_sg2]; // get the genomes that will get a m3 mut in pop2
	sg2_m3.addMutations(mutm3);  // add m3 mutation at selected position 
}


/ FITNESS LANDSCAPE UNDER PARALLEL DOMINANCE
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
    p1.setSubpopulationSize(Ndesc); 
    // Set new size for population p2
    p2.setSubpopulationSize(Ndesc); 
    // Set migration rate to migrate (symmetric m=m12=m21)
    // note: this is defined as the immigration rate
    p1.setMigrationRates(c(p2), c(MigRate));
    p2.setMigrationRates(c(p1), c(MigRate));
    
    // Every ten generations output info about the beneficial mutations
    if((sim.generation % 10) == 0) {
	    // ALLELE FREQUENCIES
	    // get the list of alleles a in both pops (mutations of type m2)
	    benmuts_m2=sim.mutationsOfType(m2);
	    // if the allele a exists save their frequency in each population
	    if(size(benmuts_m2)>0) 
	    {
	    	// get the frequency of allele a (mutations of type m2 in each pop)
	    	freqbenp1_m2 = sim.mutationFrequencies(p1, benmuts_m2);
	    	freqbenp2_m2 = sim.mutationFrequencies(p2, benmuts_m2);  
	    }
	    else
	    {
	    	// if there are no mutations of type m2 it means the frequency is zero in both pops
	    	freqbenp1_m2 = 0;
	    	freqbenp2_m2 = 0;
	    }
	    
	    // get the list of alleles A in both pops (mutations of type m3)
	    benmuts_m3=sim.mutationsOfType(m3);
	    // if the allele a exists save their frequency in each population
	    if(size(benmuts_m3)>0) 
	    {
	    	// get the frequency of allele A (mutations of type m3 in each pop)
	    	freqbenp1_m3 = sim.mutationFrequencies(p1, benmuts_m3); 
	    	freqbenp2_m3 = sim.mutationFrequencies(p2, benmuts_m3); 
	    }
	    else
	    {
	   	// if there are no mutations of type m3 it means the frequency is zero in both pops
	    	freqbenp1_m3 = 0;
	    	freqbenp2_m3 = 0;
	    }
	    // Write frequencies to a file 
	    // output file with the following columns
	    // generation freq_a_pop1 freq_a_pop2 freq_A_pop1 freq_A_pop2
	    writeFile(filenamems + "_" + rep + ".benmuts", rep + " " + sim.generation + " " + freqbenp1_m2 + " " + freqbenp2_m2 + " " + freqbenp1_m3 + " " + freqbenp2_m3 , append=T);

	    // MEAN FITNESS
	    // Compute the mean fitness of each population
	    fitp1 = mean(p1.cachedFitness(NULL));
	    fitp2 = mean(p2.cachedFitness(NULL));
	    // Write fitness to a file with the following columns
	    // generation mean_fitness_pop1 mean_fitness_pop2
	    writeFile(filenamems + "_" + rep + ".meanfit", rep + " " + sim.generation + " " + fitp1 + " " + fitp2, append=T);
    }

}

// CHECK fitness landscape
// At generation 5 print the genotypes and corresponding fitness of all individuals
// this is used to check if the fitness landscape is correctly defined
5 {
    // Define output file name tag
    filenamems = TAGFILEEND;
    // Get the fitness of each individual of pop1 and pop2
    p1fit=p1.cachedFitness(NULL);
    p2fit=p2.cachedFitness(NULL);

    // Get the genotype and fitness for each individual of pop1    
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

	// get a sample of genomes from p1 and p2 that are not null (only X chromosome sampled)
	gp1 = p1.genomes;
	gp2 = p2.genomes;

	// get number of individuals and index of firt male index (genomes are sorted)
	// Note that the index can be used to distinguish males and females, but in this case
	// we sample genomes irrespective of males or females since we use index from 0 to ((p1.individualCount*2)-1)
	gp1male = gp1[0:((p1.individualCount*2)-1)]; 
	// sample 20 genomes from pop1 (discarding the null genomes)
        sg1 = sample(gp1male[!gp1male.isNullGenome], 20, F);        
        // repeat for pop2
	gp2male = gp2[0:((p2.individualCount*2)-1)]; // sample all
        // sample 20 genomes from pop2
        sg2 = sample(gp2male[!gp2male.isNullGenome], 20, F);
        // Print ms output with genomes from pop1 first followed by genomes of pop2
        sg = c(sg1, sg2); 
        sg.outputMS("ms_" + filenamems + "_" + rep + ".txt", append=F);
}

// run to last generation
ENDGEN late()
{
	// Output ms file with sample from populations
	filenamems = TAGFILEEND;

	// get a sample of genomes from p1 and p2 that are not null (only X chromosome sampled when looking at hemizygous loci)
	gp1 = p1.genomes;
	gp2 = p2.genomes;

	// get number of individuals and of firt male index (genomes are sorted)
	// Note that the index can be used to distinguish males and females, but in this case
	// we sample genomes irrespective of males or females since we use index from 0 to ((p1.individualCount*2)-1)
	gp1male = gp1[0:((p1.individualCount*2)-1)]; 
	// sample 20 genomes from pop1
        sg1 = sample(gp1male[!gp1male.isNullGenome], 20, F);        
        // repeat for pop2
	gp2male = gp2[0:((p2.individualCount*2)-1)]; // sample all
        // sample 20 genomes from pop2
        sg2 = sample(gp2male[!gp2male.isNullGenome], 20, F);
        // Print ms output with genomes from pop1 first followed by genomes of pop2
        sg = c(sg1, sg2); 
        sg.outputMS("ms_" + filenamems + "_" + rep + ".txt", append=F);
}