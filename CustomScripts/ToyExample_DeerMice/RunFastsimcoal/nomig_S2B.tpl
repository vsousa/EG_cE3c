//Parameters for the coalescence simulation program : fastsimcoal2.exe
3 samples to simulate
//Population effective sizes (number of genes)
100000
$N_ON$
$N_SOUTH$
//Samples sizes
6
8
10
//Growth rates : negative growth implies population expansion RatioASW logRASW
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
2  historical event
$TDIV1$ 1 2 1  $ResDiv1$ 0 0
$TDIV2$ 0 2 1  $ResDiv2$ 0 0
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 3.67e-8 OUTEXP
