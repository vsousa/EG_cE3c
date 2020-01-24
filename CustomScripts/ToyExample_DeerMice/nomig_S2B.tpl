//Parameters for the coalescence simulation program : fastsimcoal2.exe
3 samples to simulate
//Population effective sizes (number of genes) for each pop in the same order as in SFS. In this case all parameters are relative to NPOP_OFFsouth, which is fixed to 100,000.
$N_SOUTH$
$N_ON$
100000
//Samples sizes (min individuals used in ProcessVCF.sh of 1,2,3 diploid individuals corresponds to 2,4,6 gene copies
2
4
6
//Growth rates : negative growth implies population expansion RatioASW logRASW
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
2  historical event
$TDIV1$ 1 0 1  $ResDiv1$ 0 0
$TDIV2$ 2 0 1  $ResDiv2$ 0 0
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 0 OUTEXP
