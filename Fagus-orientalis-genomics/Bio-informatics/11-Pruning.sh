#!/bin/bash   

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=30G
##SBATCH	--error=/dev/null
##SBATCH	--output=/dev/null

# Required programs
module load vcftools

# Input/Output folder
input=~/grant_456/project_data/genotypes/

# Edit the header of the vcf file because angsd stupidely adds its signature on the vcf version line, which creates problems in many downstream analyses...
zcat ${input}/Fagus_orientalis_poly.vcf.gz | sed 's/##fileformat=VCFv4.2(angsd version)/##fileformat=VCFv4.2/g' | gzip > ${input}/Fagus_orientalis_poly2.vcf.gz

# The vcf file had to be bgzipped and indexed as below before running that script 
vcftools --gzvcf ${input}/Fagus_orientalis_poly2.vcf.gz --thin 1000 --out ${input}/Fagus_orientalis_poly_pruned --recode

# Remove temporary vcf file
rm ${input}/Fagus_orientalis_poly2.vcf.gz
