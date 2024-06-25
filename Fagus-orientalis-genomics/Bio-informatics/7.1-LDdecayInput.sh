#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=30G
#SBATCH --array=0-3%2
##SBATCH        --error=/dev/null
##SBATCH        --output=/dev/null

# Defining the file to process using the Array ID
files=(~/*_bam.list)
file=${files[$SLURM_ARRAY_TASK_ID]}

# Retreive the name of the species from the file name
pop=${file/_bam.list}
name=`basename ${pop}`

# Input/Output folders
input=~/grant_456/scratch/Regions
#mkdir ~/grant_456/scratch/LD
output=~/grant_456/scratch/LD

# Required programs
#module load angsd/0.929
#module load htslib

# Reference genome
#ref=~/grant_456/project_data/ReferenceGenome/GCA_907173295.1_Bhaga_Chr_genomic.fna

# Create beagle file 
#angsd -b ${file} -GL 1 -out ${output}/${name}_poly -ref ${ref} -anc ${ref} -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 3 -setMinDepthInd 3 -setMaxDepthInd 100 -skipTriallelic 1 -snp_pval 1e-6 -doMaf 1 -doMajorMinor 1 -doGlf 2 -doCounts 1 

# Creating a pos file
#zcat ${output}/${name}_poly.beagle.gz | sed '1d' | awk '{split($1,a,"_")}{print a[1]"\t"a[2]}' > ${output}/pos_${name}.txt

# Estimating LD
#COUNT=`cat ${file} | wc -l`
#N_SITES=$((`cat ${output}/pos_${name}.txt | wc -l`))

# Estime LD with ngsLD
#~/TOOLS/ngsLD/ngsLD --geno ${output}/${name}_poly.beagle.gz --probs --n_ind ${COUNT} --n_sites ${N_SITES} --pos ${output}/pos_${name}.txt --out ${output}/${name}.LD --min_maf 0.05 --rnd_sample 0.001 --max_kb_dist 50

# Subsample the output file
cat ${output}/${name}.LD | awk 'BEGIN {srand()} !/^$/ { if (rand() <= .0001) print $0}' > ${output}/${name}_subset.LD
