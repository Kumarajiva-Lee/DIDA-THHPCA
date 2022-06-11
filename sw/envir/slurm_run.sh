#!/bin/bash

work_dir=$1
namelist=$2
nodes=4
nodes=$3

if (( $# > 1 )); then
  #if [[ -d $work_dir ]]; then
    #rm -r $work_dir
  #fi

  #mkdir $work_dir

  if [[ ! -f $namelist ]]; then
    echo '[Error]: namelists.input missing!'
    exit 1
  else
      sh setconf2redis.sh $2
  fi

  atm_group=`grep -i atm_group $namelist |grep -v num| awk '{print $3}' | sed $'s/\'//g'`
  da_group=`grep -i da_group $namelist |grep -v num| awk '{print $3}' | sed $'s/\'//g'`
  atm_group_num=`grep -i atm_group_num $namelist | awk '{print $3}' | sed $'s/\'//g'`
  da_group_num=`grep -i da_group_num $namelist | awk '{print $3}' | sed $'s/\'//g'`
  procs=`expr $atm_group_num \* $atm_group + $da_group_num \* $da_group`
  echo $procs
  cp $namelist $work_dir/
  sed -i "s/#SBATCH -M/#SBATCH -N $nodes/g" dida.sbatch
  sed -i "/procs=/a\procs=$procs" dida.sbatch
  sed -i "/namelist=/a\namelist=$namelist" dida.sbatch
  cp dida.sbatch $work_dir/
  cd $work_dir
  ln -s ../../build/dida.exe dida.exe
  mkdir timing
  sbatch dida.sbatch $procs $namelist

  #srun -N $nodes -n $procs --ntasks-per-socket 14 ./dida.exe $namelist
  cd ..
else
  echo '[Error]: Input missing!'
  echo '[Usage]: ./run.sh work_dir namelist'
  exit 1
fi

