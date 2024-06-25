#!/bin/bash   

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=16G
##SBATCH	--error=/dev/null
##SBATCH	--output=/dev/null

# Directory with bam files
DIR=~/grant_456/project_data/bam

# Estimates coverage for each window along the genome
for file in ${DIR}/*.final.bam
do
	sample=${file/.final.bam/}
	name=`basename ${sample}`
	
	# Find per base coverage for each window
	~/TOOLS/samtools-1.17/samtools bedcov windows.bed ${file} -Q 20 -d 3 > ${DIR}/coverage_${name}.txt
done

while IFS= read -r line
do
	# Get the information about the window
	CHROM=$( echo ${line} | awk '{print $1}' )
	POS=$( echo ${line} | awk '{print $2}' )
	
	# Number of SNPs on that window
	NB_SNP=`zcat ~/grant_456/project_data/genotypes/Fagus_orientalis_poly.mafs.gz | awk -v chrom=${CHROM} -v pos=${POS} 'FNR>1{if ($1==chrom && $2>=(pos-50000) && $2<=(pos+50000)) print $0}' | wc -l`

	touch NB_snps_windows.txt
	echo ${CHROM} ${POS} ${NB_SNP} >> NB_snps_windows.txt 

done < windows.txt
