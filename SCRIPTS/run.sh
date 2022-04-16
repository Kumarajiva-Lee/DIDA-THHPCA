#!/bin/bash

work_dir=$1
namelist=$2
atm_num=$3
da_num=$4
let total_num=atm_num+da_num
echo $total_num

if (( $# == 4 )); then
  if [[ -d $work_dir ]]; then
    rm -r $work_dir
  fi
  mkdir $work_dir
  if [[ ! -f "$namelist" ]]; then
    echo '[Error]: namelist.input missing!'
    exit 1
  fi
  cp $namelist $work_dir/namelists.input

  cd $work_dir

  #ln -s /home/export/base/shisuan/swthhpca/online/data/nmsc/share/inifile/gmcore/360x180_0.5h/l32 l32
  #ln -s ../../../../DIDA-SHARE/share/inifile/gmcore/360x180_0.5h/l32 ./l32

  #cp ../../build/dida.exe .
  sed -i "20,20c     atm_group                 = $atm_num" namelists.input
  sed -i "22,22c     atm_stride                = $total_num" namelists.input
  sed -i "24,24c     da_group                  = $da_num" namelists.input
  sed -i "26,26c     da_stride                 = $total_num" namelists.input
  sed -i "27,27c     da_root                   = $atm_num" namelists.input
  
  #bsub -I -q q_test_ss -share_size 15000 -host_stack 1024 -n $3 -cgsp 64 -b ../../build/dida.exe $namelist
  bsub -J $1 -p -o $1.log -q q_share -share_size 15000 -host_stack 1024 -cgsp 64 -b -n $total_num -exclu ../../build/dida.exe namelists.input
  #bsub -I -J $1 -p -o -q q_test_thhpca -share_size 15000 -host_stack 1024 -cgsp 64 -b -n $3 -exclu ../../build/dida.exe $namelist

  cd ..
else
  echo '[Error]: Input missing!'
  echo '[Usage]: ./run.sh work_dir namelist atm_process_number da_process_number'
  exit 1
fi

