#!/bin/bash


#######################分层图 ############################################


namepath='./namelist/initial_for_obstat.input'


hashkey_get=`grep -i case_name $namepath | awk '{print $3}' | sed $'s/\'//g'`

redis_node_get=`grep -i redis_node $namepath | awk '{print $3}' | sed $'s/\'//g'`

draw_level=`grep -i draw_level $namepath | awk '{print $3}' | sed $'s/\'//g'`
draw_single=`grep -i draw_single $namepath | awk '{print $3}' | sed $'s/\'//g'`
draw_global=`grep -i draw_global $namepath | awk '{print $3}' | sed $'s/\'//g'`

numblock=`grep -i numblock $namepath | awk '{print $3}' | sed $'s/\'//g'`
redis_port=`grep -i redis_port $namepath | awk '{print $3}' | sed $'s/\'//g'`

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
mpas_num_lev_start='mpas_num_lev_start'


num_z=$(redis-cli -h $redis_node -c -p $redis_port hget $hashkey $mpas_num_lev)
num_z_int=(`echo $num_z | awk '{print int($0)}'`)

echo $num_z


delta_z=$(redis-cli -h $redis_node -c -p $redis_port hget ${hashkey} ${obs_stat_interval_vertical})
delta_z_int=(`echo $delta_z | awk '{print int($0)}'`)


echo $num_z_int
echo $delta_z_int

mpas_start=$(redis-cli -h $redis_node -c -p $redis_port hget ${hashkey} ${mpas_num_lev_start})
mpas_start_int=(`echo $mpas_start | awk '{print int($0)}'`)

echo $mpas_start

minus=`expr $num_z_int - $mpas_start_int`
length_norm_lev=`expr $minus / $delta_z_int + 1`
echo $length_norm_lev

m_g=`grep -i mpas_or_gmcore $namepath | awk '{print $3}' | sed $'s/\'//g'`

################################################################################

vari_surf=("ps")


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

if [[ ! -d $1  ]]
then
    mkdir -p $1
fi

pic_path=$(cd $1 && pwd)

echo "$pic_path"

cd ./main/
#cd ./run/

if [ $draw_level == 1 ]
then
   n=`expr $length_norm_lev \* $length_norm_vari + $length_surf_vari`

   echo "$n"

   p=`expr $n \* $numblock`

   if [ $m_g == 0 ]
   then
      #srun -n $p python consumer-level-nc-mpas.py $pic_path
      srun -p cnall --mpi=pmi2 -n $p python consumer-level-nc-mpas.py $pic_path $2 $3
   else
      srun -p cnall --mpi=pmi2 -n $p python consumer-level-nc-update.py $pic_path $2 $3
      #mpiexec -n $p python consumer-level-nc-update.py $pic_path
      #sbatch dida.sbatch $p $pic_path
   fi

fi


if [ $draw_single == 1 ]
then
   n=$length_vari
   echo "$n"
   #srun -N 1 -n $n -w cn003 python consumer-single.py $pic_path
fi


if [ $draw_global == 1 ]
then
   #n=$length_vari
   n=`expr $length_norm_lev \* $length_norm_vari + $length_surf_vari`
   echo "$n"
   #srun -N 1 -n $n python consumer-global-gather.py $pic_path
fi










