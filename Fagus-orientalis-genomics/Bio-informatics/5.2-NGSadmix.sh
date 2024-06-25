#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=30G
##SBATCH        --error=/dev/null
##SBATCH        --output=/dev/null

# Create directory to store the results
mkdir ~/pl0278-01/scratch/NGSadmix

# Run NGS admix 5 times for each K and for K ranging from 1 to 6
for K in 5 6 
do
#	~/pl0278-01/scratch/TOOLS/NGSadmix/NGSadmix -likes ~/grant_456/project_data/genotypes/Fagus_orientalis_poly.beagle.gz -K ${K} -o `echo ~/pl0278-01/scratch/NGSadmix/Fagus_orientalis_K${K}_1` -minMaf 0.03 -minInd 120
#	~/pl0278-01/scratch/TOOLS/NGSadmix/NGSadmix -likes ~/grant_456/project_data/genotypes/Fagus_orientalis_poly.beagle.gz -K ${K} -o `echo ~/pl0278-01/scratch/NGSadmix/Fagus_orientalis_K${K}_2` -minMaf 0.03 -minInd 120
#	~/pl0278-01/scratch/TOOLS/NGSadmix/NGSadmix -likes ~/grant_456/project_data/genotypes/Fagus_orientalis_poly.beagle.gz -K ${K} -o `echo ~/pl0278-01/scratch/NGSadmix/Fagus_orientalis_K${K}_3` -minMaf 0.03 -minInd 120
#	~/pl0278-01/scratch/TOOLS/NGSadmix/NGSadmix -likes ~/grant_456/project_data/genotypes/Fagus_orientalis_poly.beagle.gz -K ${K} -o `echo ~/pl0278-01/scratch/NGSadmix/Fagus_orientalis_K${K}_4` -minMaf 0.03 -minInd 120
#	~/pl0278-01/scratch/TOOLS/NGSadmix/NGSadmix -likes ~/grant_456/project_data/genotypes/Fagus_orientalis_poly.beagle.gz -K ${K} -o `echo ~/pl0278-01/scratch/NGSadmix/Fagus_orientalis_K${K}_5` -minMaf 0.03 -minInd 120
  	~/pl0278-01/scratch/TOOLS/NGSadmix/NGSadmix -likes ~/grant_456/project_data/genotypes/Fagus_orientalis_poly.beagle.gz -K ${K} -o `echo ~/pl0278-01/scratch/NGSadmix/Fagus_orientalis_K${K}_6` -minMaf 0.03 -minInd 120
        ~/pl0278-01/scratch/TOOLS/NGSadmix/NGSadmix -likes ~/grant_456/project_data/genotypes/Fagus_orientalis_poly.beagle.gz -K ${K} -o `echo ~/pl0278-01/scratch/NGSadmix/Fagus_orientalis_K${K}_7` -minMaf 0.03 -minInd 120
        ~/pl0278-01/scratch/TOOLS/NGSadmix/NGSadmix -likes ~/grant_456/project_data/genotypes/Fagus_orientalis_poly.beagle.gz -K ${K} -o `echo ~/pl0278-01/scratch/NGSadmix/Fagus_orientalis_K${K}_8` -minMaf 0.03 -minInd 120
        ~/pl0278-01/scratch/TOOLS/NGSadmix/NGSadmix -likes ~/grant_456/project_data/genotypes/Fagus_orientalis_poly.beagle.gz -K ${K} -o `echo ~/pl0278-01/scratch/NGSadmix/Fagus_orientalis_K${K}_9` -minMaf 0.03 -minInd 120
        ~/pl0278-01/scratch/TOOLS/NGSadmix/NGSadmix -likes ~/grant_456/project_data/genotypes/Fagus_orientalis_poly.beagle.gz -K ${K} -o `echo ~/pl0278-01/scratch/NGSadmix/Fagus_orientalis_K${K}_10` -minMaf 0.03 -minInd 120
  	~/pl0278-01/scratch/TOOLS/NGSadmix/NGSadmix -likes ~/grant_456/project_data/genotypes/Fagus_orientalis_poly.beagle.gz -K ${K} -o `echo ~/pl0278-01/scratch/NGSadmix/Fagus_orientalis_K${K}_11` -minMaf 0.03 -minInd 120
        ~/pl0278-01/scratch/TOOLS/NGSadmix/NGSadmix -likes ~/grant_456/project_data/genotypes/Fagus_orientalis_poly.beagle.gz -K ${K} -o `echo ~/pl0278-01/scratch/NGSadmix/Fagus_orientalis_K${K}_12` -minMaf 0.03 -minInd 120
        ~/pl0278-01/scratch/TOOLS/NGSadmix/NGSadmix -likes ~/grant_456/project_data/genotypes/Fagus_orientalis_poly.beagle.gz -K ${K} -o `echo ~/pl0278-01/scratch/NGSadmix/Fagus_orientalis_K${K}_13` -minMaf 0.03 -minInd 120
        ~/pl0278-01/scratch/TOOLS/NGSadmix/NGSadmix -likes ~/grant_456/project_data/genotypes/Fagus_orientalis_poly.beagle.gz -K ${K} -o `echo ~/pl0278-01/scratch/NGSadmix/Fagus_orientalis_K${K}_14` -minMaf 0.03 -minInd 120
        ~/pl0278-01/scratch/TOOLS/NGSadmix/NGSadmix -likes ~/grant_456/project_data/genotypes/Fagus_orientalis_poly.beagle.gz -K ${K} -o `echo ~/pl0278-01/scratch/NGSadmix/Fagus_orientalis_K${K}_15` -minMaf 0.03 -minInd 120
done

# Retreive run likelihoods 
#rm ~/pl0278-01/scratch/NGSadmix/likelihoods.txt
#echo "K" "Likelihood" >> ~/pl0278-01/scratch/NGSadmix/likelihoods.txt
#for log in ~/pl0278-01/scratch/NGSadmix/*.log
#do 
#	f=`basename ${log}`
#	K=`echo ${log} | head -c 28 | tail -c 1`
#	likelihood=`grep -Po 'like=\K[^ ]+' ${log}`
#	echo ${K} ${likelihood} >> ~/pl0278-01/scratch/NGSadmix/likelihoods.txt
#done

