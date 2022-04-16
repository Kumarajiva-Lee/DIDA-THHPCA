#!/bin/bash


#######################分层图 ############################################


namepath='./namelist/initial_for_obstat.input'

hashkey_get=`sed '/^   case_name                 = /!d;s/.*=//' $namepath`

hashkey_get=${hashkey_get#*\'}
hashkey_get=${hashkey_get%\'*}

redis_node_get=`sed '/^   redis_node                = /!d;s/.*=//' $namepath`
redis_node_get=${redis_node_get#*\'}
redis_node_get=${redis_node_get%\'*}

draw_level=`sed '/^   draw_level                = /!d;s/.*=//' $namepath`
draw_single=`sed '/^   draw_single               = /!d;s/.*=//' $namepath`
draw_global=`sed '/^   draw_global               = /!d;s/.*=//' $namepath`

echo "draw_level:$draw_level"
echo "draw_single:$draw_single"
echo "draw_global:$draw_global"

echo $hashkey_get
#num_z_int=(`echo $num_z | awk '{print int($0)}'`)

hashkey=$hashkey_get
redis_node=$redis_node_get

echo $hashkey
echo $redis_node_get

mpas_num_lev='mpas_num_lev'
obs_stat_interval_vertical='obs_stat_interval_vertical'


num_z=$(redis-cli -h $redis_node -c hget $hashkey $mpas_num_lev)
num_z_int=(`echo $num_z | awk '{print int($0)}'`)

echo $num_z


delta_z=$(redis-cli -h $redis_node -c hget ${hashkey} ${obs_stat_interval_vertical})
delta_z_int=(`echo $delta_z | awk '{print int($0)}'`)


echo $num_z_int
echo $delta_z_int

length_norm_lev=`expr $num_z_int / $delta_z_int + 1`
echo $length_norm_lev


################################################################################

vari_surf=("ps")

#obs_var_name='obs_var_name'
#vari_get=$(redis-cli -h $redis_node -c hget ${hashkey} ${obs_var_name})

vari_get=`sed '/^   obsstat_var_name          = /!d;s/.*=//' $namepath`
vari_get=${vari_get#*\'}
vari_get=${vari_get%\'*}
vari_norm_surf=(`echo $vari_get | tr ',' ' '`)
length_vari=${#vari_norm_surf[@]}
echo $length_vari

echo ${vari_norm_surf[@]}

length_surf_vari=0
echo $length_surf_vari
for x in ${vari_norm_surf[@]}
do
    [[ ${vari_surf[@]/$x/} != ${vari_surf[@]} ]] && length_surf_vari=`expr $length_surf_vari + 1` || echo $x
done
echo $length_surf_vari


length_norm_vari=`expr $length_vari - $length_surf_vari`



if [ $draw_level == 1 ]
then
   n=`expr $length_norm_lev \* $length_norm_vari + $length_surf_vari`

   echo "$n"

   #mpiexec -np $n python consumer-mpi-multi-real.py

   srun -N 1 -n $n -w cn001 python consumer-level.py
fi


if [ $draw_single == 1 ]
then
   n=$length_vari
   echo "$n"
   srun -N 1 -n $n -w cn002 python consumer-single.py
fi


if [ $draw_global == 1 ]
then
   n=$length_vari
   echo "$n"
   srun -N 1 -n $n -w cn003 python consumer-global.py
fi










