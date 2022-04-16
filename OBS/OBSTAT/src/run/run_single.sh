#!/bin/bash


#SBATCH -J single
#SBATCH -o single.out
#SBATCH -p ib.cluster
#SBATCH --qos=debug
#SBATCH --ntasks-per-node=30
#SBATCH --nodes=1
#SBATCH -t 1:00:00
#SBATCH -w cn002

a=$1
b=$2
cd ../main/
srun -N 1 -n $a python consumer-single.py $b


