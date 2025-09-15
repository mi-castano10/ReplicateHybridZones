//Parameters for the coalescence simulation program : simcoal.exe
3 samples to simulate :
//Population effective sizes (number of genes)
NICT
NHYB
NFLA
//Samples sizes and samples age 
10
10
10
//Growth rates: negative growth implies population expansion
GROWI
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//migration matrix 0
0 RECENTMIGRATEIH 0
RECENTRATEMIGHI 0 MIGRATEHF$
0 RATEMIGFH$ 0
//migration matrix 1
0 MIGRATEIH$ 0
RATEMIGHI$ 0 MIGRATEHF$
0 RATEMIGFH$ 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
4 historical event
TADM2 0 0 0 1 0 1 //change migration matrix to 1
TDIV1 1 2 1 1 0 1 //all hybrids (source) goes back into ancestral flam (sink)
TDIV1 0 2 1 1 0 1 //all ict (source) goes back into ancestral flam (sink)
TDIV1 2 2 1 ANCSIZE 0 1 absoluteResize //all flam (source) goes back into ancestral flam (sink)
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0 4.47e-9 OUTEXP //using mutation rate derived from geospiza fortis (1.72e-9 per year) (Nadachowska-Brzyska 2015) and using a generation time of 2.6 years (Bird 2020)
