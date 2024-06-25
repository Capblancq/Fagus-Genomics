#!/bin/bash   

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=4-00:00:00
#SBATCH --mem-per-cpu=30G
##SBATCH	--error=/dev/null
##SBATCH	--output=/dev/null

# Required programs
module load angsd/0.929
module load htslib
module load python/3.9.2

# Reference genome
ref=~/grant_456/project_data/ReferenceGenome/GCA_907173295.1_Bhaga_Chr_genomic.fna

# Moving to the right directory
cd ~/grant_456/scratch/Regions

# Create dadi compatible SFS using realSFS
realSFS dadi -P 10 LesserCaucasus.saf.idx GreaterCaucasusWest.saf.idx Hyrcanian.saf.idx -sfs LesserCaucasus.sfs -sfs GreaterCaucasusWest.sfs -sfs Hyrcanian.sfs -ref ${ref} -anc ${ref} -nSites 20000000 > LC_GCW_HY.dadi

# Randomly subsample the files that are too big to be treated otherwise...
cat LC_GCW_HY.dadi | awk 'BEGIN  {srand()} !/^$/  { if (rand() <= .1 || FNR==1) print}' > LC_GCW_HY.dadi.sample 

# Transform the SFS created by realSFS into dadi like SFS
perl ~/TOOLS/realsfs2dadi.pl LC_GCW_HY.dadi.sample 35 59 10 > LC_GCW_HY.dadi.sample.sfs
sed -i 's/pop0/LC/g' LC_GCW_HY.dadi.sample.sfs
sed -i 's/pop1/GCW/g' LC_GCW_HY.dadi.sample.sfs
sed -i 's/pop2/HY/g' LC_GCW_HY.dadi.sample.sfs

# Downsample the sfs using dadi
python3.9 - << EOF

import dadi

fs = dadi.Misc.make_data_dict('LC_GCW_HY.dadi.sample.sfs')
fs_LC_GCW_HY = dadi.Spectrum.from_data_dict(fs,pop_ids=['LC','GCW','HY'],projections=[70,118,20],polarized=True)
fs_LC_GCW_HY.to_file('LC_GCW_HY_DSFS.obs')

fs_LC_GCW_HY = dadi.Spectrum.from_data_dict(fs,pop_ids=['LC','GCW','HY'],projections=[20,20,20],polarized=True)
fs_LC_GCW_HY.to_file('LC_GCW_HY_downto20_DSFS.obs')

EOF

# Remove tmp files
rm *.dadi.sample
rm *.dadi.sample.sfs
