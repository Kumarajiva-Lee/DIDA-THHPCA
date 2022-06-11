base_namelist=namelists.input
for sceneid in $( seq 1 8 )
do
#sceneid=1
colid=`expr ${sceneid} - 1`
num_line=`grep -c "" scene.ini`
echo $num_line
target_namelist="$base_namelist$sceneid"
cp $base_namelist $target_namelist

for x in $( seq 1 `expr ${num_line} - 0`  )
#for x in $( seq 1 10  )
do

  line_message=`sed -n "${x}p" scene.ini`
  param_cut=(`echo $line_message | tr '=' ' '`)
  param=${param_cut[0]}

  if [ $param != "start_time" ]; then
    value_get=`grep -i ${param} scene.ini |grep -v mpas_num_lev| awk '{print $3}' `
    value_cut=(`echo $value_get | tr ';' ' '`)
    value=${value_cut[$colid]}
    echo $value


    value_bak=`grep -i $param $target_namelist |grep -v mpas_num_lev| awk '{print $3}'`
    vv=(${value_bak/ //;})
    value_line=`grep -n $param $target_namelist |grep -v mpas_num_lev |awk -F ":" '{print $1}' `
    #echo ${#value_line[@]}
    #echo ${value_line[@]}
    mm=(${value_line/ //,})
    echo ${#mm[@]}
    if (( ${#mm[@]} > 1  )); then
        for y in $( seq 0 `expr ${#mm[@]} - 1` )
        do
          mmid=$((10#${mm[$y]}+0))
          echo $mmid
          sed -i "${mmid}s/${vv[0]}/$value/g" $target_namelist
        done
    else
        sed -i "${value_line}s/$value_bak/$value/g" $target_namelist
    fi


  else


    line_message=`sed -n "${x}p" scene.ini`
    param_cut=(`echo $line_message | tr '=' ' '`)
    param=${param_cut[0]}
    value_get=`grep -i ${param} scene.ini |grep -v da_start_time| awk '{print $3}' `
    value_cut=(`echo $value_get | tr ';' ' '`)
    value=${value_cut[$colid]}
    echo $value


    value_bak=`grep -i $param $target_namelist |grep -v da_start_time| awk '{print $3}'`
    vv=(${value_bak/ //;})
    value_line=`grep -n $param $target_namelist |grep -v da_start_time |awk -F ":" '{print $1}' `
    #echo ${#value_line[@]}
    #echo ${value_line[@]}
    mm=(${value_line/ //,})
    echo ${#mm[@]}
    if (( ${#mm[@]} > 1  )); then
        for y in $( seq 0 `expr ${#mm[@]} - 1` )
        do
          mmid=$((10#${mm[$y]}+0))
          echo $mmid
          sed -i "${mmid}s/${vv[0]}/$value/g" $target_namelist
        done
    else
        sed -i "${value_line}s/$value_bak/$value/g" $target_namelist
    fi

  fi


done

done
