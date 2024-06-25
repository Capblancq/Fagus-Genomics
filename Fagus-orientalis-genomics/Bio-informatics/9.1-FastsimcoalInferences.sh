#!/bin/bash

#SBATCH --job-name=Fastsimcoal2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=4-00:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --array=36-50%5
##SBATCH        --error=/dev/null
##SBATCH        --output=/dev/null

# Run ID
i=$SLURM_ARRAY_TASK_ID

# Move to working directory
cd ~/pl0278-01/scratch/Fastsimcoal/

# Create directory for the scenario
mkdir RES_3pops
mkdir ./RES_3pops/run${i}

# Copy and rename input files
cp Forientalis3popsDiv.tpl Forientalis3popsDiv.est ./RES_3pops/run${i}/
cp LC_GCW_HY_downto20_DSFS.obs ./RES_3pops/run${i}/Forientalis3popsDiv_DSFS.obs
#cp Forientalis4popsDiv_simple.tpl Forientalis4popsDiv_simple.est ./RES/run${i}/
#cp LC_GCW_jointDAFpop1_0.obs ./RES/run${i}/Forientalis4popsDiv_jointDAFpop1_0.obs
#cp LC_GCE_jointDAFpop2_0.obs ./RES/run${i}/Forientalis4popsDiv_jointDAFpop2_0.obs
#cp LC_HY_jointDAFpop3_0.obs ./RES/run${i}/Forientalis4popsDiv_jointDAFpop3_0.obs
#cp GCW_GCE_jointDAFpop2_1.obs ./RES/run${i}/Forientalis4popsDiv_jointDAFpop2_1.obs
#cp GCW_HY_jointDAFpop3_1.obs ./RES/run${i}/Forientalis4popsDiv_jointDAFpop3_1.obs
#cp GCE_HY_jointDAFpop3_2.obs ./RES/run${i}/Forientalis4popsDiv_jointDAFpop3_2.obs

# Enter run i directory
cd ./RES_3pops/run${i}

# Run fastsimcoal 3 pops
~/TOOLS/fsc27_linux64/fsc27093 -t Forientalis3popsDiv.tpl -n1000000 -d -e Forientalis3popsDiv.est -M -L60 -q -c12 --multiSFS

# Run fastsimcaol 4 pops
#~/TOOLS/fsc27_linux64/fsc27093 -t Forientalis4popsDiv_simple.tpl -n1000000 -d -e Forientalis4popsDiv_simple.est -M -L60 -q -c10
