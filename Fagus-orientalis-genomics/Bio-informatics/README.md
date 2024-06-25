

1. "1-fastqTrimReads.sh", "2-fastqMapping.sh", "3-removingPCRduplicates.sh" and "4-SNPcallingANGSD.sh" describe how we used Trimmomatic, bwa and ANGSD (http://www.popgen.dk/angsd/index.php/ANGSD#Overview) to produce genotype likelihoods from .fastq files

2. "5.1-PCAngsd.sh", "5.2-NGSadmix.sh", "6.1-GenDivSites.sh", "6.2-GenDivRegions.sh", "7.1-LDdecayInput.sh", "7.2-LDwindows.sh" and "8-GeneticLoad.sh" describe how we used genotype likelihoods to produce different estimates of genetic diversity (thetas, genetic load), genetic structure (PCA, admixture proportions and Fst) and linkage disiquilibrium (LD).

3. 
