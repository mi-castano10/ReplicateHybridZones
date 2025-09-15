//Number of populations (demes)
3 
//Population effective sizes (number of genes)
NICT
NHYB
NFLA
//Samples sizes and samples age 
10
10
10
//Growth rates: negative growth implies population expansion
0
0
0
//Number of migration matrices : 0 implies no migration between demes
4
//migration matrix 0
0 RECENTMIGRATEIH 0
RECENTRATEMIGHI 0 RECENTMIGRATEHF
0 RECENTRATEMIGFH 0
//migration matrix 1
0 MIGRATEIH$ 0
RATEMIGHI$ 0 MIGRATEHF$
0 RATEMIGFH$ 0
//migration matrix 2
0 0 0
0 0 MIGRATEHF$
0 RATEMIGFH$ 0
//migration matrix 3
0 0 0
0 0 0
0 0 0
//historical event : time, source, sink, migrants, new deme size, new growth rate, new migration matrix index
4 historical event
TADM3 0 0 0 1 0 1 //no resizing, at TADM3 change to migration matrix 1 (from smaller to larger migration in all demes)
TADM2 0 0 0 1 0 2 //no resizing, at TADM2 change to migration matrix 2 (no migration between hybrids and ict)
TDIV2 1 2 1 RESIZE 0 3 absoluteResize //all hyb (source) goes back into ancestral flam (sink) and TDIV2 = TADM1, gene flow between HF stops change to migration matrix 2
TDIV1 0 2 1 ANCSIZE 0 3 absoluteResize //all ict (source) goes back into ancestral flam (sink)
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0 4.47e-9 OUTEXP //using mutation rate derived from geospiza fortis (1.72e-9 per year) (Nadachowska-Brzyska 2015) and using a generation time of 2.6 years (Bird 2020)
