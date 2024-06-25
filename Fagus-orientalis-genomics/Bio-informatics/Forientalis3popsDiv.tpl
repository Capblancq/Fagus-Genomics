//Number of population samples (demes)
3
//Population effective sizes (number of genes)
N_LC
N_GCW
N_HY
//Samples sizes and samples age
20
20
20
//Growth rates: negative growth implies population expansion
R_LC
R_GCW
R_HY
//Number of migration matrices : 0 implies no migration between demes
3
//Migration matrix 0
0 m_LC_GCW m_LC_HY
m_LC_GCW 0 m_GCW_HY
m_LC_HY m_GCW_HY 0
//Migration matrix 1
0 0 m_anc
0 0 0
m_anc 0 0
//Migration matrix 2
0 0 0
0 0 0
0 0 0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
4 historical event
T1 1 0 1 ResizeTIME1 R_anc 1
T1 1 1 0 0 0 1 // Deme 1 is killed
T2 2 0 1 ResizeTIME2 0 2
T2 2 2 0 0 0 2 // Deme 2 is killed
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 2e-9 OUTEXP
