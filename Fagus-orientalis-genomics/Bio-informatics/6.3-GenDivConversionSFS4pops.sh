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
module load r

# Reference genome
ref=~/grant_456/project_data/ReferenceGenome/GCA_907173295.1_Bhaga_Chr_genomic.fna

# Moving to the right directory
cd ~/grant_456/scratch/Regions

# Create dadi compatible SFS using realSFS
#realSFS dadi -P 10 LesserCaucasus.saf.idx GreaterCaucasusWest.saf.idx -sfs LesserCaucasus.sfs -sfs GreaterCaucasusWest.sfs -ref ${ref} -anc ${ref} -nSites 20000000 > LC_GCW.dadi
#realSFS dadi -P 10 LesserCaucasus.saf.idx GreaterCaucasusEast.saf.idx -sfs LesserCaucasus.sfs -sfs GreaterCaucasusEast.sfs -ref ${ref} -anc ${ref} -nSites 20000000 > LC_GCE.dadi
realSFS dadi -P 10 LesserCaucasus.saf.idx Hyrcanian.saf.idx -sfs LesserCaucasus.sfs -sfs Hyrcanian.sfs -ref ${ref} -anc ${ref} -nSites 20000000 > LC_HY.dadi
#realSFS dadi -P 10 GreaterCaucasusWest.saf.idx GreaterCaucasusEast.saf.idx -sfs GreaterCaucasusWest.sfs -sfs GreaterCaucasusEast.sfs -ref ${ref} -anc ${ref} -nSites 20000000 > GCW_GCE.dadi
#realSFS dadi -P 10 GreaterCaucasusWest.saf.idx Hyrcanian.saf.idx -sfs GreaterCaucasusWest.sfs -sfs Hyrcanian.sfs -ref ${ref} -anc ${ref} -nSites 20000000 > GCW_HY.dadi
#realSFS dadi -P 10 GreaterCaucasusEast.saf.idx Hyrcanian.saf.idx -sfs GreaterCaucasusEast.sfs -sfs Hyrcanian.sfs -ref ${ref} -anc ${ref} -nSites 20000000 > GCE_HY.dadi

# Randomly subsample the files that are too big to be treated otherwise...
#cat LC_GCW.dadi | awk 'BEGIN  {srand()} !/^$/  { if (rand() <= .1 || FNR==1) print}' > LC_GCW.dadi.sample 
#cat LC_GCE.dadi | awk 'BEGIN  {srand()} !/^$/  { if (rand() <= .1 || FNR==1) print}' > LC_GCE.dadi.sample
cat LC_HY.dadi | awk 'BEGIN  {srand()} !/^$/  { if (rand() <= .1 || FNR==1) print}' > LC_HY.dadi.sample
#cat GCW_GCE.dadi | awk 'BEGIN  {srand()} !/^$/  { if (rand() <= .1 || FNR==1) print}' > GCW_GCE.dadi.sample
#cat GCW_HY.dadi | awk 'BEGIN  {srand()} !/^$/  { if (rand() <= .1 || FNR==1) print}' > GCW_HY.dadi.sample
#cat GCE_HY.dadi | awk 'BEGIN  {srand()} !/^$/  { if (rand() <= .1 || FNR==1) print}' > GCE_HY.dadi.sample

# Transform the SFS created by realSFS into dadi like SFS
#perl ~/TOOLS/realsfs2dadi.pl LC_GCW.dadi.sample 35 59 > LC_GCW.dadi.sample.sfs
#sed -i 's/pop0/LC/g' LC_GCW.dadi.sample.sfs
#sed -i 's/pop1/GCW/g' LC_GCW.dadi.sample.sfs
#perl ~/TOOLS/realsfs2dadi.pl LC_GCE.dadi.sample 35 49 > LC_GCE.dadi.sample.sfs
#sed -i 's/pop0/LC/g' LC_GCE.dadi.sample.sfs
#sed -i 's/pop1/GCE/g' LC_GCE.dadi.sample.sfs
perl ~/TOOLS/realsfs2dadi.pl LC_HY.dadi.sample 35 10 > LC_HY.dadi.sample.sfs
sed -i 's/pop0/LC/g' LC_HY.dadi.sample.sfs
sed -i 's/pop1/HY/g' LC_HY.dadi.sample.sfs
#perl ~/TOOLS/realsfs2dadi.pl GCW_GCE.dadi.sample 59 49 > GCW_GCE.dadi.sample.sfs
#sed -i 's/pop0/GCW/g' GCW_GCE.dadi.sample.sfs
#sed -i 's/pop1/GCE/g' GCW_GCE.dadi.sample.sfs
#perl ~/TOOLS/realsfs2dadi.pl GCW_HY.dadi.sample 59 10 > GCW_HY.dadi.sample.sfs
#sed -i 's/pop0/GCW/g' GCW_HY.dadi.sample.sfs
#sed -i 's/pop1/HY/g' GCW_HY.dadi.sample.sfs
#perl ~/TOOLS/realsfs2dadi.pl GCE_HY.dadi.sample 49 10 > GCE_HY.dadi.sample.sfs
#sed -i 's/pop0/GCE/g' GCE_HY.dadi.sample.sfs
#sed -i 's/pop1/HY/g' GCE_HY.dadi.sample.sfs

# Downsample the sfs using dadi
python3.9 - << EOF

import dadi

#fs = dadi.Misc.make_data_dict('LC_GCW.dadi.sample.sfs')
#fs_LC_GCW = dadi.Spectrum.from_data_dict(fs,pop_ids=['LC','GCW'],projections=[20,20],polarized=True)
#fs_LC_GCW.to_file('LC_GCW_jointDAFpop1_0.obs')

#fs = dadi.Misc.make_data_dict('LC_GCE.dadi.sample.sfs')
#fs_LC_GCE = dadi.Spectrum.from_data_dict(fs,pop_ids=['LC','GCE'],projections=[20,20],polarized=True)
#fs_LC_GCE.to_file('LC_GCE_jointDAFpop2_0.obs')

fs = dadi.Misc.make_data_dict('LC_HY.dadi.sample.sfs')
fs_LC_HY = dadi.Spectrum.from_data_dict(fs,pop_ids=['LC','HY'],projections=[20,20],polarized=True)
fs_LC_HY.to_file('LC_HY_jointDAFpop3_0.obs')

#fs = dadi.Misc.make_data_dict('GCW_GCE.dadi.sample.sfs')
#fs_GCW_GCE = dadi.Spectrum.from_data_dict(fs,pop_ids=['GCW','GCE'],projections=[20,20],polarized=True)
#fs_GCW_GCE.to_file('GCW_GCE_jointDAFpop2_1.obs')

#fs = dadi.Misc.make_data_dict('GCW_HY.dadi.sample.sfs')
#fs_GCW_HY = dadi.Spectrum.from_data_dict(fs,pop_ids=['GCW','HY'],projections=[20,20],polarized=True)
#fs_GCW_HY.to_file('GCW_HY_jointDAFpop3_1.obs')

#fs = dadi.Misc.make_data_dict('GCE_HY.dadi.sample.sfs')
#fs_GCE_HY = dadi.Spectrum.from_data_dict(fs,pop_ids=['GCE','HY'],projections=[20,20],polarized=True)
#fs_GCE_HY.to_file('GCE_HY_jointDAFpop3_2.obs')

EOF

R --vanilla << EOF

N1 <- 10
N2 <- 10

#sfs <- scan(paste("LC_GCW_jointDAFpop1_0.obs", sep=""), quiet=T, skip = 1, nlines = 1)
#tab <- matrix(data = sfs, ncol = (N1*2+1), nrow = (N2*2+1))
#colnames(tab) <- paste(rep("d0", (N1*2+1)), seq(0,(N1*2)), sep = "_")
#row.names(tab) <- paste(rep("d1", (N2*2+1)), seq(0,(N2*2)), sep = "_")
#write.table(tab, "~/grant_456/scratch/Fastsimcoal/LC_GCW_jointDAFpop1_0.obs", row.names = T, col.names = T, sep = "\t", quote = F)

#sfs <- scan(paste("LC_GCE_jointDAFpop2_0.obs", sep=""), quiet=T, skip = 1, nlines = 1)
#tab <- matrix(data = sfs, ncol = (N1*2+1), nrow = (N2*2+1))
#colnames(tab) <- paste(rep("d0", (N1*2+1)), seq(0,(N1*2)), sep = "_")
#row.names(tab) <- paste(rep("d1", (N2*2+1)), seq(0,(N2*2)), sep = "_")
#write.table(tab, "~/grant_456/scratch/Fastsimcoal/LC_GCE_jointDAFpop2_0.obs", row.names = T, col.names = T, sep = "\t", quote = F)

sfs <- scan(paste("LC_HY_jointDAFpop3_0.obs", sep=""), quiet=T, skip = 1, nlines = 1)
tab <- matrix(data = sfs, ncol = (N1*2+1), nrow = (N2*2+1))
colnames(tab) <- paste(rep("d0", (N1*2+1)), seq(0,(N1*2)), sep = "_")
row.names(tab) <- paste(rep("d1", (N2*2+1)), seq(0,(N2*2)), sep = "_")
write.table(tab, "~/grant_456/scratch/Fastsimcoal/LC_HY_jointDAFpop3_0.obs", row.names = T, col.names = T, sep = "\t", quote = F)

#sfs <- scan(paste("GCW_GCE_jointDAFpop2_1.obs", sep=""), quiet=T, skip = 1, nlines = 1)
#tab <- matrix(data = sfs, ncol = (N1*2+1), nrow = (N2*2+1))
#colnames(tab) <- paste(rep("d0", (N1*2+1)), seq(0,(N1*2)), sep = "_")
#row.names(tab) <- paste(rep("d1", (N2*2+1)), seq(0,(N2*2)), sep = "_")
#write.table(tab, "~/grant_456/scratch/Fastsimcoal/GCW_GCE_jointDAFpop2_1.obs", row.names = T, col.names = T, sep = "\t", quote = F)

#sfs <- scan(paste("GCW_HY_jointDAFpop3_1.obs", sep=""), quiet=T, skip = 1, nlines = 1)
#tab <- matrix(data = sfs, ncol = (N1*2+1), nrow = (N2*2+1))
#colnames(tab) <- paste(rep("d0", (N1*2+1)), seq(0,(N1*2)), sep = "_")
#row.names(tab) <- paste(rep("d1", (N2*2+1)), seq(0,(N2*2)), sep = "_")
#write.table(tab, "~/grant_456/scratch/Fastsimcoal/GCW_HY_jointDAFpop3_1.obs", row.names = T, col.names = T, sep = "\t", quote = F)

#sfs <- scan(paste("GCE_HY_jointDAFpop3_2.obs", sep=""), quiet=T, skip = 1, nlines = 1)
#tab <- matrix(data = sfs, ncol = (N1*2+1), nrow = (N2*2+1))
#colnames(tab) <- paste(rep("d0", (N1*2+1)), seq(0,(N1*2)), sep = "_")
#row.names(tab) <- paste(rep("d1", (N2*2+1)), seq(0,(N2*2)), sep = "_")
#write.table(tab, "~/grant_456/scratch/Fastsimcoal/GCE_HY_jointDAFpop3_2.obs", row.names = T, col.names = T, sep = "\t", quote = F)

EOF

# Remove tmp files
rm *.dadi.sample
rm *.dadi.sample.sfs
rm *.dadi
rm *_jointDAFpop*.obs
