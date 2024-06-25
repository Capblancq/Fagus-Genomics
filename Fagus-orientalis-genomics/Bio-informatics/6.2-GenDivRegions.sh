#!/bin/bash   

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=30G
##SBATCH	--error=/dev/null
##SBATCH	--output=/dev/null

### Regional population SFS ###

# Required programs
module load angsd/0.929
module load htslib

# Output directory
output=~/grant_456/scratch/Regions

# Reference genome
ref=~/grant_456/project_data/ReferenceGenome/GCA_907173295.1_Bhaga_Chr_genomic.fna

# Find polymorphic sites for each regional population
#for list in *_bam.list
#do
#	region=${list/_bam.list}
#	name=`basename ${region}`
#	angsd -b ${list} -ref ${ref} -anc ${ref} -out ${output}/${name} -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 5 -setMinDepthInd 4 -setMaxDepthInd 100 -skipTriallelic 1 -GL 1 -doMaf 1 -doMajorMinor 1 -doCounts 1
#	zcat ${output}/${name}.mafs.gz | awk 'BEGIN { OFS = ":" }{ print $1,$2 }' | sed '1d' > ${output}/${name}_sites.txt
#	rm ${output}/${name}.mafs.gz ${output}/${name}.arg
#done

# Finding loci intersection among regions
#comm -12 ${output}/West_sites.txt ${output}/East_sites.txt > ${output}/West_East_sites.txt
#comm -12 ${output}/West_East_sites.txt ${output}/Admixed_sites.txt > ${output}/West_East_Admixed_sites.txt 

#sed 's/:/\t/' ${output}/West_East_Admixed_sites.txt > ${output}/intersect.txt #| sort -b -k1,1
#cut -f1 ${output}/intersect.txt | uniq > ${output}/intersect.chrs
#angsd sites index ${output}/intersect.txt

# Re-estimating the genotype likelihoods with only intersecting sites (monomorphic and polymorphic) + optimization of the SFS
for list in *_bam.list
do
        region=${list/_bam.list}
        name=`basename ${region}`
	# SFS estimaiton
#	angsd -b ${list} -GL 1 -out ${output}/${name} -ref ${ref} -anc ${ref} -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 3 -setMinDepthInd 3 -setMaxDepthInd 100 -skipTriallelic 1 -doSaf 1 -fold 0 -doMaf 1 -doMajorMinor 1 -doCounts 1 #-sites ${output}/intersect.txt -rf ${output}/intersect.chrs
	# EM optimization of the sfs
#	realSFS ${output}/${name}.saf.idx -maxIter 50000 -tole 1e-6 -P 10 > ${output}/${name}.sfs
	# Estimating thetas from the SFS
#	realSFS saf2theta ${output}/${name}.saf.idx -outname ${output}/${name} -sfs ${output}/${name}.sfs
	thetaStat do_stat ${output}/${name}.thetas.idx -win 100000 -step 100000 -type 0
done

# Optimization of the 2D SFS for the four groups
#realSFS ${output}/GreaterCaucasusWest.saf.idx ${output}/GreaterCaucasusEast.saf.idx -P 10 -nSites 50000000 -maxIter 100 > ${output}/GreaterCaucasusWest_GreaterCaucasusEast.ml
# Average the different SFS obtained from each 10mio sites windows above
#awk '{ for(i=1;i<=NF;i++) total[i]+=$i ; } END { for(i=1;i<=NF;i++) printf "%f ",total[i]/NR ;}' ${output}/GreaterCaucasusWest_GreaterCaucasusEast.ml > ${output}/GreaterCaucasusWest_GreaterCaucasusEast.sfs 
# Remove tmp file
#rm ${output}/GreaterCaucasusWest_GreaterCaucasusEast.ml

# Same for all pairs
#realSFS ${output}/GreaterCaucasusWest.saf.idx ${output}/LesserCaucasus.saf.idx -maxIter 100 -tole 1e-6 -P 10 -nSites 50000000 -maxIter 100 > ${output}/GreaterCaucasusWest_LesserCaucasus.ml
#awk '{ for(i=1;i<=NF;i++) total[i]+=$i ; } END { for(i=1;i<=NF;i++) printf "%f ",total[i]/NR ;}' ${output}/GreaterCaucasusWest_LesserCaucasus.ml > ${output}/GreaterCaucasusWest_LesserCaucasus.sfs
#rm ${output}/GreaterCaucasusWest_LesserCaucasus.ml
#realSFS ${output}/GreaterCaucasusWest.saf.idx ${output}/Hyrcanian.saf.idx -maxIter 100 -tole 1e-6 -P 10 -nSites 50000000 -maxIter 100 > ${output}/GreaterCaucasusWest_Hyrcanian.ml
#awk '{ for(i=1;i<=NF;i++) total[i]+=$i ; } END { for(i=1;i<=NF;i++) printf "%f ",total[i]/NR ;}' ${output}/GreaterCaucasusWest_Hyrcanian.ml > ${output}/GreaterCaucasusWest_Hyrcanian.sfs
#rm ${output}/GreaterCaucasusWest_Hyrcanian.ml
#realSFS ${output}/GreaterCaucasusEast.saf.idx ${output}/LesserCaucasus.saf.idx -maxIter 100 -tole 1e-6 -P 4 -nSites 10000000 > ${output}/GreaterCaucasusEast_LesserCaucasus.ml
#awk '{ for(i=1;i<=NF;i++) total[i]+=$i ; } END { for(i=1;i<=NF;i++) printf "%f ",total[i]/NR ;}' ${output}/GreaterCaucasusEast_LesserCaucasus.ml > ${output}/GreaterCaucasusEast_LesserCaucasus.sfs
#rm ${output}/GreaterCaucasusEast_LesserCaucasus.ml
#realSFS ${output}/GreaterCaucasusEast.saf.idx ${output}/Hyrcanian.saf.idx -maxIter 100 -tole 1e-6 -P 4 -nSites 10000000 > ${output}/GreaterCaucasusEast_Hyrcanian.ml
#awk '{ for(i=1;i<=NF;i++) total[i]+=$i ; } END { for(i=1;i<=NF;i++) printf "%f ",total[i]/NR ;}' ${output}/GreaterCaucasusEast_Hyrcanian.ml > ${output}/GreaterCaucasusEast_Hyrcanian.sfs
#rm ${output}/GreaterCaucasusEast_Hyrcanian.ml
#realSFS ${output}/LesserCaucasus.saf.idx ${output}/Hyrcanian.saf.idx -maxIter 100 -tole 1e-6 -P 4 -nSites 10000000 > ${output}/LesserCaucasus_Hyrcanian.ml
#awk '{ for(i=1;i<=NF;i++) total[i]+=$i ; } END { for(i=1;i<=NF;i++) printf "%f ",total[i]/NR ;}' ${output}/LesserCaucasus_Hyrcanian.ml > ${output}/LesserCaucasus_Hyrcanian.sfs
#rm ${output}/LesserCaucasus_Hyrcanian.ml

# Prepare the Fst
#realSFS fst index ${output}/GreaterCaucasusWest.saf.idx ${output}/GreaterCaucasusEast.saf.idx ${output}/Hyrcanian.saf.idx \
#	-sfs ${output}/GreaterCaucasusWest_GreaterCaucasusEast.sfs \
#	-sfs ${output}/GreaterCaucasusWest_Hyrcanian.sfs \
#	-sfs ${output}/GreaterCaucasusEast_Hyrcanian.sfs \
#	-fstout ${output}/GCW_GCE_H_Fst

#realSFS fst index ${output}/GreaterCaucasusWest.saf.idx ${output}/LesserCaucasus.saf.idx ${output}/Hyrcanian.saf.idx \
#        -sfs ${output}/GreaterCaucasusWest_LesserCaucasus.sfs \
#        -sfs ${output}/GreaterCaucasusWest_Hyrcanian.sfs \
#        -sfs ${output}/LesserCaucasus_Hyrcanian.sfs \
#        -fstout ${output}/GCW_LC_H_Fst

#realSFS fst index ${output}/GreaterCaucasusEast.saf.idx ${output}/LesserCaucasus.saf.idx \
#        -sfs ${output}/GreaterCaucasusEast_LesserCaucasus.sfs \
#        -fstout ${output}/GCE_LC_Fst

# Get the estimate along sliding windows
#realSFS fst stats2 ${output}/GCW_GCE_H_Fst.fst.idx -win 10000 -step 10000 -type 0 > ${output}/Fst_sliding_GCW_GCE_H
#realSFS fst stats2 ${output}/GCW_LC_H_Fst.fst.idx -win 10000 -step 10000 -type 0 > ${output}/Fst_sliding_GCW_LC_H
#realSFS fst stats2 ${output}/GCE_LC_Fst.fst.idx -win 10000 -step 10000 -type 0 > ${output}/Fst_sliding_GCE_LC
