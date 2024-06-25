#!/bin/bash

#SBATCH --job-name=Fastsimcoal2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=20G
##SBATCH        --error=/dev/null
##SBATCH        --output=/dev/null

# Move to the right directory
cd ~/pl0278-01/scratch/Fastsimcoal/RES_3pops

# Find the best run
BEST_RUN=`cat ./run*/Forientalis3popsDiv/Forientalis3popsDiv.bestlhoods | grep -v MaxObsLhood | awk '{print NR,$16}' | sort -k 2 | head -n 1 | awk '{print $1}'`

# Directory
mkdir bootstrap
cd bootstrap

# Retreive the .par file produced for the best run
#cp ../run${BEST_RUN}/Forientalis3popsDiv/Forientalis3popsDiv_maxL.par ./Forientalis3popsDiv_boot.par

# Change a bit the .par file to 
#sed -i 's/^1 0$/200000 0/g' Forientalis3popsDiv_boot.par 
#sed -i 's/^FREQ 1/DNA 100/g' Forientalis3popsDiv_boot.par

# Then generate 100 SFS
#~/TOOLS/fsc27_linux64/fsc27093 -i Forientalis3popsDiv_boot.par -n100 -j -d -s0 -x -I -q -c10

# Running the optimization with each bootstraped SFS
for i in {65..100}
do
	cp ../../Forientalis3popsDiv.tpl ./Forientalis3popsDiv_boot/Forientalis3popsDiv_boot_$i/Forientalis3popsDiv_boot.tpl
	cp ../../Forientalis3popsDiv.est ./Forientalis3popsDiv_boot/Forientalis3popsDiv_boot_$i/Forientalis3popsDiv_boot.est
	cp ../run${BEST_RUN}/Forientalis3popsDiv/Forientalis3popsDiv.pv ./Forientalis3popsDiv_boot/Forientalis3popsDiv_boot_$i/Forientalis3popsDiv_boot.pv
	cd ./Forientalis3popsDiv_boot/Forientalis3popsDiv_boot_$i
	~/TOOLS/fsc27_linux64/fsc27093 -t Forientalis3popsDiv_boot.tpl -e Forientalis3popsDiv_boot.est -n1000000 -d -M -L50 --initValues Forientalis3popsDiv_boot.pv -c10 -q
	cd ../../
done

# Estimating confidence intervals
#cat ./Forientalis3popsDiv_boot/Forientalis3popsDiv_boot_1/Forientalis3popsDiv_boot/Forientalis3popsDiv_boot.bestlhoods | awk 'NR==1' > Parameters_bootstrap.txt
#for i in {1..100}
#do
#	cat ./Forientalis3popsDiv_boot/Forientalis3popsDiv_boot_$i/Forientalis3popsDiv_boot/Forientalis3popsDiv_boot.bestlhoods | awk 'NR==2' >> Parameters_bootstrap.txt
#done

#R --vanilla <<EOF
#
#	TAB <- read.table("Parameters_bootstrap.txt", header=T)
#	values <- read.table("../run${BEST_RUN}/Forientalis3popsDiv/Forientalis3popsDiv.bestlhoods", header=T)
#	n <- nrow(TAB)
#	CI <- list()
#	for (i in 1:10){
#	a <- values[i]
#	s <- sd(TAB[,i])
#	error <- qt(0.975,df=n-1)*s/sqrt(n)
#	left <- a-error
#	right <- a+error
#	CI[[i]] <- cbind(left, a, right)
#	}
#	CI <- lapply(CI, as.numeric)
#	CI <- do.call(cbind,CI)
#	colnames(CI) <- colnames(values)[1:10]
#	write.table(CI, "Parameters_CI.txt", row.names=F, quote=F)
#
#EOF

