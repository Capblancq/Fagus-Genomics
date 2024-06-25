#!/bin/bash
#SBATCH --job-name=LD_windows
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --array=0-100%50 #5,000 windows
#SBATCH        --error=/dev/null
#SBATCH        --output=/dev/null

##SLURM_ARRAY_TASK_ID=$((${SLURM_ARRAY_TASK_ID}+10000))

# Extract chrom and pos of the target window
CHROM=$( sed -n $((${SLURM_ARRAY_TASK_ID}+1))p windows.txt | awk '{print $1}' )
POS=$( sed -n $((${SLURM_ARRAY_TASK_ID}+1))p windows.txt | awk '{print $2}' )

echo ${CHROM} ${POS}

# Input/Output folders
INPUT=~/grant_456/scratch/LD
OUTPUT=~/grant_456/scratch/LD

# Source R
module load r

for population in GreaterCaucasusWest GreaterCaucasusEast LesserCaucasus Hyrcanian
do

  # Subsetting the beagle files based on the windows
  zcat ${INPUT}/${population}_poly.beagle.gz | awk -v chrom=${CHROM} -v pos=${POS} 'FNR==1{print}{split($1,a,"_")}{if (a[1]==chrom && a[2]>=(pos-50000) && a[2]<=(pos+50000)) print $0}' | gzip > ${OUTPUT}/${population}_window${SLURM_ARRAY_TASK_ID}.beagle.gz

  # Creating a pos file
  zcat ${OUTPUT}/${population}_window${SLURM_ARRAY_TASK_ID}.beagle.gz | sed '1d' | awk '{split($1,a,"_")}{print a[1]"\t"a[2]}' > ${OUTPUT}/pos_${population}_window${SLURM_ARRAY_TASK_ID}.txt

  # Estimating LD
  COUNT=`cat ${population}_bam.list | wc -l`
  N_SITES=$((`cat ${OUTPUT}/pos_${population}_window${SLURM_ARRAY_TASK_ID}.txt | wc -l`))

  # Creating output file
  touch ${OUTPUT}/${population}_LDr2bp.txt
  
  # Retreive mean r2 on the window 
  if (("${N_SITES}" > "200"))
  then
	#Estime LD with ngsLD
	~/TOOLS/ngsLD/ngsLD --geno ${OUTPUT}/${population}_window${SLURM_ARRAY_TASK_ID}.beagle.gz --probs --n_ind ${COUNT} --n_sites ${N_SITES} --pos ${OUTPUT}/pos_${population}_window${SLURM_ARRAY_TASK_ID}.txt --out ${OUTPUT}/${population}_window${SLURM_ARRAY_TASK_ID}.LD

	# Estimate r2 per base
	R2=`awk '{ print $4/$3 }' ${OUTPUT}/${population}_window${SLURM_ARRAY_TASK_ID}.LD | awk '{sum+=$1} END { print sum/NR}'`	

	# Feeding the output file
 	echo ${CHROM} ${POS} ${R2} ${N_SITES} >> ${OUTPUT}/${population}_LDr2bp.txt
  else
	echo ${CHROM} ${POS} "NA" ${N_SITES} >> ${OUTPUT}/${population}_LDr2bp.txt
  fi

  # Removing temporary files
  rm ${OUTPUT}/${population}_window${SLURM_ARRAY_TASK_ID}.beagle.gz
  rm ${OUTPUT}/pos_${population}_window${SLURM_ARRAY_TASK_ID}.txt
  rm ${OUTPUT}/${population}_window${SLURM_ARRAY_TASK_ID}.LD

done
