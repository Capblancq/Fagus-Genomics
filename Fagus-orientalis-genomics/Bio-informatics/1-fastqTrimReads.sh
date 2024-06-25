#!/bin/bash   

#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=12G
##SBATCH	--error=/dev/null
##SBATCH	--output=/dev/null

DIR=~/grant_456/project_data/fastq
OUT=~/grant_456/scratch

mkdir ${OUT}/cleanreads/

## Trimmomatic 

for f1 in ${DIR}/100_1.fastq.gz  
do 
 	f=${f1/_1.fastq.gz/}
	name=`basename ${f}`
	# Run trimmomatic
	~/TOOLS/jre1.8.0_301/bin/java -classpath ~/TOOLS/Trimmomatic-0.39/trimmomatic-0.39.jar org.usadellab.trimmomatic.TrimmomaticPE \
        -threads 20 \
        -phred33 \
        -basein "$f1" \
        -baseout "${OUT}/cleanreads/${name}.fastq.gz" \
        ILLUMINACLIP:/home/users/tcapblan/TOOLS/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:20 \
        TRAILING:20 \
        SLIDINGWINDOW:6:20 \
        MINLEN:35 
done 

## Statistics on the reads

#touch ${DIR}/NBreads_1P.txt
#touch ${DIR}/NBreads_1U.txt
#touch ${DIR}/NBreads_2U.txt
#touch ${DIR}/NBreads_names.txt
#for forward in ${OUT}/cleanreads/*_1P.fastq.gz
#do
#	R1_only=${forward/_1P.fastq.gz/_1U.fastq.gz}
#	R2_only=${forward/_1P.fastq.gz/_2U.fastq.gz}
#	f=${forward/_1P.fastq.gz/}
#	name=`basename ${f}`
#	echo ${name} >> ${DIR}/NBreads_names.txt
#	echo $(zcat ${forward} | wc -l)/4|bc >> ${DIR}/NBreads_1P.txt
#	echo $(zcat ${R1_only} | wc -l)/4|bc >> ${DIR}/NBreads_1U.txt
#	echo $(zcat ${R2_only} | wc -l)/4|bc >> ${DIR}/NBreads_2U.txt
#done

#R --vanilla <<EOF
#
#	reads_R1 = read.table("NBreads_1P.txt",header = FALSE)
#	reads_R1only = read.table("NBreads_1U.txt",header = FALSE)
#	reads_R2only = read.table("NBreads_2U.txt",header = FALSE)
#	ind_names = read.table("NBreads_names.txt", header = FALSE)
#	seq_res = cbind (ind_names, reads_R1, reads_R1only, reads_R2only)
#	colnames(seq_res) = c("Individuals", "Reads_paired", "Reads_R1_only", "Reads_R2_only")
#	write.table(seq_res, "sequencing_results.txt", row.names = F, quote = F)
#	
#EOF
