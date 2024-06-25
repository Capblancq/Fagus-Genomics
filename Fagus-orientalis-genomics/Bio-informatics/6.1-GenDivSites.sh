#!/bin/bash   

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --array=0-32%10
##SBATCH	--error=/dev/null
##SBATCH	--output=/dev/null

# Defining the file to process using the Array ID
lists=(~/grant_456/scratch/Populations/*_bam.list)
list=${lists[$SLURM_ARRAY_TASK_ID]}

### Population genetic diversity estimates ###

# Required programs
module load angsd/0.929
module load htslib

# Output directory
output=~/grant_456/scratch/Populations 

# Reference genome
ref=~/grant_456/project_data/ReferenceGenome/GCA_907173295.1_Bhaga_Chr_genomic.fna

# Name of the population
pop=${list/_bam.list}
name=`basename ${pop}`

# Estimating the SFS for each population for all sites (both mono- and poly-morphic sites)
angsd -b ${list} -anc ${ref} -ref ${ref} -out ${output}/${name} -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 2 -setMinDepthInd 4 -setMaxDepthInd 100 -setMinDepth 50 -setMaxDepth 500 -GL 1  -doCounts 1 -doSaf 1
realSFS ${output}/${name}.saf.idx -P 10 > ${output}/${name}.sfs

# Estimating thetas from the SFS
realSFS saf2theta ${output}/${name}.saf.idx -outname ${output}/${name} -sfs ${output}/${name}.sfs
thetaStat do_stat ${output}/${name}.thetas.idx

# Removing intermediary files
rm ${output}/${name}.saf.idx
rm ${output}/${name}.saf.gz
rm ${output}/${name}.saf.pos.gz
rm ${output}/${name}.thetas.idx
rm ${output}/${name}.thetas.gz
