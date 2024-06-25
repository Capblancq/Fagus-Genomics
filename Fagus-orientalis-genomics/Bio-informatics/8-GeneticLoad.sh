#!/bin/bash   

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=30G
##SBATCH	--error=/dev/null
##SBATCH	--output=/dev/null

# Required programs
module load java8/jdk1.8.0_40
module load vcftools
module load bcftools

# Input/Output folder
input=~/grant_456/project_data/genotypes/
mkdir ~/grant_456/scratch/Populations
output=~/grant_456/scratch/Populations

## Building a new genome database for snpEff
#nano ~/TOOLS/snpEff/snpEff.config
## Adding 
#    # Fagus sylvatica
#    Fsylvatica.genome : Fagus_sylvatica
#mkdir ~/TOOLS/snpEff/data/Fsylvatica
#cp ~/grant_456/project_data/ReferenceGenome/GCA_907173295.1_Bhaga_Chr_genomic.fna ~/TOOLS/snpEff/data/Fsylvatica/sequences.fa # genome
#cp ~/grant_456/project_data/ReferenceGenome/GCA_907173295.1_Bhaga_Chr_genomic.gff3 ~/TOOLS/snpEff/data/Fsylvatica/genes.gff # annotations

## Build the new dataset/reference for F. sylvatica
#java -jar ~/TOOLS/snpEff/snpEff.jar build -gff3 -v Fsylvatica

# Using a vcf file with all the polymorphic position from ANGSD 
#java -Xmx30G -jar ~/TOOLS/snpEff/snpEff.jar -c ~/TOOLS/snpEff/snpEff.config -v Fsylvatica ${input}/Fagus_orientalis_poly.vcf.gz > ${input}/Fagus_orientalis_poly_annotated.vcf

## The vcf file had to be bgzipped and indexed as below before running that script 
#bcftools view -I ${input}/Fagus_orientalis_poly_annotated.vcf -O z -o ${input}/Fagus_orientalis_poly_annotated.vcf.gz
#bcftools index ${input}/Fagus_orientalis_poly_annotated.vcf.gz

# Split the vcf file based on sample locality
bcftools view -s ind88,ind89,ind90,ind91,ind92 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/AH_01.vcf
#bcftools view -s ind93,ind94,ind95,ind96,ind97 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/AH_02.vcf
#bcftools view -s ind98,ind99,ind100,ind102,ind162 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/AZ_01.vcf
#bcftools view -s ind103,ind104,ind105,ind106,ind163 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/AZ_02.vcf
#bcftools view -s ind108,ind109,ind110,ind112,ind113 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/AZ_03.vcf
#bcftools view -s ind114,ind115,ind116,ind117,ind118 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/AZ_04.vcf
#bcftools view -s ind130,ind131,ind132,ind133,ind134 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/AZ_05.vcf
#bcftools view -s ind125,ind126,ind127,ind128,ind129 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/AZ_06.vcf
#bcftools view -s ind153,ind154,ind155,ind160 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/AZ_07.vcf
#bcftools view -s ind120,ind122,ind123,ind124,ind161 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/AZ_08.vcf
#bcftools view -s ind136,ind137,ind138,ind139,ind140 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/AZ_09.vcf
#bcftools view -s ind141,ind142,ind143,ind144,ind145 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/AZ_10.vcf
#bcftools view -s ind148,ind149,ind150,ind151,ind152 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/AZ_11.vcf
#bcftools view -s ind2,ind3,ind4,ind5,ind6 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/GC_01.vcf
#bcftools view -s ind7,ind16,ind19,ind22,ind30 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/GC_02.vcf
#bcftools view -s ind55,ind60,ind64,ind68,ind73 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/GC_03.vcf
#bcftools view -s ind8,ind9,ind10,ind11,ind12 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/GC_04.vcf
#bcftools view -s ind13,ind14,ind15,ind17,ind18 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/GC_05.vcf
#bcftools view -s ind40,ind41,ind42,ind43,ind44 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/GC_06.vcf
#bcftools view -s ind35,ind36,ind37,ind38,ind39 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/GC_07.vcf
#bcftools view -s ind83,ind84,ind85,ind86,ind87 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/GC_08.vcf
#bcftools view -s ind78,ind79,ind80,ind81,ind82 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/GC_09.vcf
#bcftools view -s ind31,ind32,ind33,ind34,ind159 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/GC_10.vcf
#bcftools view -s ind20,ind21,ind23,ind24,ind156 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/GC_11.vcf
#bcftools view -s ind101,ind107,ind111,ind119,ind121 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/GC_12.vcf
#bcftools view -s ind74,ind75,ind76,ind77,ind158 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/LC_01.vcf
#bcftools view -s ind69,ind70,ind71,ind72,ind157 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/LC_02.vcf
#bcftools view -s ind25,ind26,ind27,ind28,ind29 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/LC_03.vcf
#bcftools view -s ind62,ind63,ind65,ind66,ind67 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/LC_04.vcf
#bcftools view -s ind56,ind57,ind58,ind59,ind61 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/LC_05.vcf
#bcftools view -s ind0,ind1,ind135,ind146,ind147 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/LC_06.vcf
#bcftools view -s ind45,ind46,ind47,ind48,ind49 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/LC_07.vcf
#bcftools view -s ind50,ind51,ind52,ind53,ind54 ${input}/Fagus_orientalis_poly_annotated.vcf.gz > ${output}/LC_08.vcf

# Prepare files with annotations and allele frequencies per site
#for VCF in ${output}/*.vcf
#do 
#  path=${VCF/.vcf/}
#  name=`basename ${path}`
#
#  sed -i 's/##fileformat=VCFv4.2(angsd version)/##fileformat=VCFv4.2/g' ${VCF}
#  vcftools --vcf ${VCF} --freq2 --derived --out ${output}/${name}
#  vcftools --vcf ${VCF} --hardy --out ${output}/${name}

#  cat ${VCF} | grep -v "#" | awk '{print $8}' | awk -F '|' '{print $2}' > ${output}/${name}.tmp
#  cat ${VCF} | grep -v "#" | awk '{print $1"_"$2}' > ${output}/${name}.sites
#  paste ${output}/${name}.sites ${output}/${name}.tmp > ${output}/${name}.annotations
#  cat ${output}/${name}.frq | awk 'FNR>1 {print $1"_"$2,$3,$4,$5}' > ${output}/${name}.freq
#  cat ${output}/${name}.hwe | awk 'FNR>1 {print $1"_"$2,$3}' > ${output}/${name}.gen

#  join ${output}/${name}.freq ${output}/${name}.gen > ${output}/${name}_FreqGen.tmp
#  join ${output}/${name}_FreqGen.tmp ${output}/${name}.annotations > ${output}/${name}_FreqGenAnnot.tmp
  
#  cat ${output}/${name}_FreqGenAnnot.tmp | grep "missense_variant\|synonymous_variant\|stop_gained\|start_lost" > ${output}/${name}_FreqAnnotations.txt 
  
#  rm ${VCF}
#  rm ${output}/${name}.frq
#  rm ${output}/${name}.freq
#  rm ${output}/${name}.hwe
#  rm ${output}/${name}.tmp
#  rm ${output}/${name}.sites
#  rm ${output}/${name}.log
#  rm ${output}/${name}.gen
#  rm ${output}/${name}_FreqGen.tmp
#  rm ${output}/${name}_FreqGenAnnot.tmp
#  rm ${output}/${name}.annotations
#done
