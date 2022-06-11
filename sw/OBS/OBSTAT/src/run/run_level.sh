#!/bin/bash


#SBATCH -J level
#SBATCH -o level.out
#SBATCH -p cnall
#SBATCH -t 1:00:00

a=$1
b=$2
cd ../main/
srun --mpi=pmi2 -n $a python consumer-level-nc-update.py $b


