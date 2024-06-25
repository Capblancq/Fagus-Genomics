#!/bin/bash   

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=30G
##SBATCH        --error=/dev/null
##SBATCH        --output=/dev/null


# PCA using genotype likelihoods store in beagle format

/home/users/tcapblan/python/bin/python3.6 ~/TOOLS/pcangsd/pcangsd.py \
	-beagle ~/grant_456/project_data/genotypes/Fagus_orientalis_poly.beagle.gz \
	-o ~/grant_456/project_data/genotypes/Fagus_orientalis_poly \
	-threads 10 \
	-minMaf 0.03 \
	-dosage_save \
	-pi_save \
	-pcadapt \
	-snp_weights \
	-sites_save
	
