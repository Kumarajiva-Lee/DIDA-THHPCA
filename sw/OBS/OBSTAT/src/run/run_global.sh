#!/bin/bash


#SBATCH -J global
#SBATCH -o global.out
#SBATCH -p ib.cluster
#SBATCH --qos=debug
#SBATCH --ntasks-per-node=30
#SBATCH --nodes=1
#SBATCH -t 1:00:00
#SBATCH -w cn003

a=$1
b=$2
cd ../main/
srun -N 1 -n $a python consumer-global-gather.py $b


