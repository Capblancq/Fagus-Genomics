### Forest genomics in the Caucasus through the lens of its dominant tree species - Fagus orientalis
---------------------

This repository includes all the scripts used to perform the analyses and produce the figures of the study described in Capblancq et al. 2024, Molecular Ecology (XX):

1. Sequencing.R can be used to draw Figure 1 from summary statistics estimated on the seuqence aligments (.bam files).

2. Diversity_sites.R can be used to draw figure 2 from the output of the subprogram of ANGSD "thetaStat".

3. Structure_Forientalis.R describes how to use PCAngsd and NGSadmix outputs to draw figure 3.

4. Script_stairwayplot.R describes how to use STAIRWAY PLOT outputs to draw figure 4.

5. Diversity_genome.R can be used to draw figure 5 from the output of the subprogram of ANGSD "thetaStat".



3. Genotype_likelihoods_ANGSD.sh describes how we used ANGSD (http://www.popgen.dk/angsd/index.php/ANGSD#Overview) to produce genotype likelihoods from .bam files

4. Downstream_analyses_ANGSD.sh describes how we used genotype likelihoods to produce different estimates of genetic diversity (thetas), genetic structure (PCA and Fst) and linkage disiquilibrium (LD).

5. Annotations.sh describes how we annotated the resulting vcf file using the norway spruce reference genome and annotations.
