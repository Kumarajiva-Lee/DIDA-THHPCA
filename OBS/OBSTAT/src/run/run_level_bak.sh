#!/bin/bash


#SBATCH -J level
#SBATCH -o level.out
#SBATCH -o level.err
#SBATCH -p cnall
#SBATCH -n 8
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1

source /apps/compilers/intel/psxe/2019/u5/compilers_and_libraries/linux/mkl/bin/mklvars.sh intel64
source /apps/compilers/intel/psxe/2019/u5/impi/2019.5.281/intel64/bin/mpivars.sh
source /apps/compilers/intel/psxe/2019/u5/parallel_studio_xe_2019/psxevars.sh intel64
export PATH=/apps/compilers/intel/psxe/2019/u5/parallel_studio_xe_2019/compilers_and_libraries_2019/linux/compiler/include/intel64:/apps/compilers/intel/psxe/2019/u5/parallel_studio_xe_2019/compilers_and_libraries_2019/linux/mpi/intel64/include:$PATH

ulimit -s unlimited
ulimit -c unlimited


cd ../main/
#srun --mpi=pmi2 -n $a python consumer-level-nc-update.py $b
/apps/compilers/intel/psxe/2019/u5/intelpython3/bin/mpiexec.hydra -n $1 /home/xuewei/ssp/3rdParty/anaconda3/bin/python consumer-level-nc-update-corun.py $2 $3 $4

