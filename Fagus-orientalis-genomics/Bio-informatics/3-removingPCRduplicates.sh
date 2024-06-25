#!/bin/bash   

#SBATCH --ntasks=1
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=20G
##SBATCH	--error=/dev/null
##SBATCH	--output=/dev/null

#Directory to store the final bam files and stats on mapping results
output=~/grant_456/project_data/bam_HZ/

# Removing PCR duplicates
#for file in ${output}/*.sorted.bam
#do
#	f=${file/.sorted.bam/}
#	name=`basename ${f}`
#	/home/users/tcapblan/TOOLS/sambamba-0.6.8 markdup -r -t 2 ${file} ${file}.rmdup.bam --tmpdir=${output}
#	/home/users/tcapblan/TOOLS/samtools-1.3.1/samtools sort ${file}.rmdup.bam -o ${output}${name}.final.bam
#	/home/users/tcapblan/TOOLS/samtools-1.3.1/samtools index ${output}${name}.final.bam
#	rm ${file}
#	rm ${file}.rmdup.bam
#	rm ${file}.rmdup.bam.bai
#done

# Stats after PCR duplicates removal
#rm ${output}/rmdup.res.aln.reads.out
#rm ${output}/rmdup.names.txt
#touch ${output}/rmdup.names.txt
#for INDIV in ${output}/*.final.bam
#do
#	f=${INDIV/.final.bam/}
#	name=`basename ${f}`
#	echo ${name} >> ${output}/rmdup.names.txt
#	/home/users/tcapblan/TOOLS/samtools-1.3.1/samtools flagstat ${INDIV} | awk 'NR>=6&&NR<=13 {print $1}' | column -x
#done >> ${output}/rmdup.res.aln.reads.out

# Reads mapping quality scores
rm ${output}/rmdup_reads_mapping_Qscores.txt
for file in ${output}/*.final.bam
do
	/home/users/tcapblan/TOOLS/samtools-1.3.1/samtools view ${file} | awk '$5>0{c1++}; $5>29{c29++}; $5>=0{c0++}; $5>19{c19++}; $5>9{c9++} END {print c0 " " c1 " " c9 " " c19 " " c29}'
done >> ${output}/rmdup_reads_mapping_Qscores.txt

# Nucleotide coverage
rm ${output}/rmdup_mean_coverage.txt
for file in ${output}/*.final.bam
do
	/home/users/tcapblan/TOOLS/samtools-1.3.1/samtools depth ${file} | awk '{sum+=$3} END {print sum/NR}'
done >> ${output}/rmdup_mean_coverage.txt

