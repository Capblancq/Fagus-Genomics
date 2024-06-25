#!/bin/bash   

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=16G
##SBATCH	--error=/dev/null
##SBATCH	--output=/dev/null

#Directory to store the stats on mapping results
stats=~/grant_456/project_data/bam/

#Directory with demultiplexed fastq.gz files
input=~/grant_456/scratch/cleanreads/

#Directory for outputs
mkdir ~/grant_456/scratch/Mapping_bwa/
output=~/grant_456/scratch/Mapping_bwa/

#Reference for mapping (fasta file)
ref=~/grant_456/project_data/ReferenceGenome/GCA_907173295.1_Bhaga_Chr_genomic.fna

#Index the reference 
#/home/users/tcapblan/TOOLS/bwa-0.7.17/bwa index ${ref}

# Aligning individual sequences to the reference
for forward in ${input}9*_1P.fastq.gz
do
	f=${forward/_1P.fastq.gz/}
	name=`basename ${f}`
	echo "@ Aligning $name..."
	reverse=${forward/_1P.fastq.gz/_2P.fastq.gz}
	/home/users/tcapblan/TOOLS/bwa-0.7.17/bwa mem -t 10 -M -a ${ref} ${forward} ${reverse} > ${output}/${name}.sam
	/home/users/tcapblan/TOOLS/samtools-1.3.1/samtools view -bS ${output}/${name}.sam -o ${output}/${name}.bam
	/home/users/tcapblan/TOOLS/samtools-1.3.1/samtools sort ${output}/${name}.bam -o ${output}/${name}.sorted.bam
	rm ${output}/${name}.sam
	rm ${output}/${name}.bam
done

## Stats on bwa alignments
#touch ${stats}/res.aln.reads.out
#touch ${stats}/names.txt
#for INDIV in ${output}/*.sorted.bam
#do
#	f=${INDIV/.bam/}
#	name=`basename ${f}`
#	echo ${name} >> ${stats}/names.txt
#	/home/users/tcapblan/TOOLS/samtools-1.3.1/samtools flagstat ${INDIV} | awk 'NR>=6&&NR<=13 {print $1}' | column -x
#done >> ${stats}/res.aln.reads.out

## Reads mapping quality scores after bwa alignment
#touch ${stats}/reads_mapping_Qscores.txt
#for file in ${output}/*.sorted.bam
#do
#	/home/users/tcapblan/TOOLS/samtools-1.3.1/samtools view ${file} | awk '$5>0{c1++}; $5>29{c29++}; $5>=0{c0++}; $5>19{c19++}; $5>9{c9++} END {print c0 " " c1 " " c9 " " c19 " " c29}'
#done >> ${stats}/reads_mapping_Qscores.txt

## Nucleotide coverage on bwa .bam files
#touch ${stats}/mean_coverage.txt
#for file in ${output}/*.sorted.bam
#do
#	/home/users/tcapblan/TOOLS/samtools-1.3.1/samtools depth ${file} | awk '{sum+=$3} END {print sum/NR}'
#done >> ${stats}/mean_coverage.txt

