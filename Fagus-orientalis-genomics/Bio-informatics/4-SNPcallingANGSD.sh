#!/bin/bash   

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=30G
##SBATCH	--error=/dev/null
##SBATCH	--output=/dev/null

# Required programs
#module load angsd/0.929
#module load htslib

# Input folder
input=~/grant_456/project_data/bam/

# Output folder
mkdir ~/grant_456/project_data/genotypes/
output=~/grant_456/project_data/genotypes/

# Reference genome
ref=~/grant_456/project_data/ReferenceGenome/GCA_907173295.1_Bhaga_Chr_genomic.fna
#/home/users/tcapblan/TOOLS/samtools-1.3.1/samtools faidx ${ref}

#### Genotype Likelihoods ####

# List of path for the .bam files 
#ls ${input}*.final.bam > ${output}FO_bam.list

# Estimating genotype likelihoods for 165 Fagus orientalis samples (only polymorphic sites)
~/TOOLS/angsd/angsd -b ${output}FO_bam.list \
	-ref ${ref} \
	-anc ${ref} \
	-out ${output}Fagus_orientalis_poly \
	-nThreads 10 \
	-uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
	-setMinDepthInd 3 -setMaxDepthInd 50 -minInd 100 -setMaxDepth 10000 \
	-skipTriallelic 1 \
	-SNP_pval 1e-6 \
	-GL 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -doGeno 2 -geno_minDepth 3 -postCutoff 0.9 -doCounts 1 -doVcf 1
#	-doGlf 2

# Estimating genotype likelihoods for 165 Fagus orientalis samples (mono- and poly-morphic sites)
#angsd -b ${output}FO_bam.list \
#        -ref ${ref} \
#        -anc ${ref} \
#        -out ${output}Fagus_orientalis_all \
#        -nThreads 10 \
#        -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
#        -setMinDepthInd 3 -setMaxDepthInd 50 -minInd 100 -setMaxDepth 10000 \
#        -skipTriallelic 1 \
#        -GL 1 \
#        -doMaf 1 -doMajorMinor 1 -doPost 1 \
#        -doDepth 1 -doCounts 1 
#	-doSaf 1	
