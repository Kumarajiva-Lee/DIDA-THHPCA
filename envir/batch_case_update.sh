dida_dir=/home/xuewei/lustre/yuz/dida2.0/dida-new-fix-0301/dida-v2.0-int
envir_dir=/home/xuewei/lustre/yuz/dida2.0/dida-new-fix-0301/dida-v2.0-int/envir
obstat_dir=${dida_dir}/OBS/OBSTAT/src
create_name=testmpas
mpas_or_gmcore=0
nocompile_or_compile=0


CREATE(){

    if [[ -d $create_name ]]; then
        echo "testwork exists"
    else
        cd ${dida_dir}/SCRIPTS
        ./create_newcase -c $create_name
        cp $envir_dir/slurm_run.sh ${dida_dir}/${create_name}/.
    fi

}


COMPILE(){

    work_dir=${dida_dir}/${create_name}
    cd $work_dir
    sh compile.sh intel

}

NAMELIST(){

    param=(`echo ${case_group[0]} | tr ':' ' '`)
    length_param=${#param[@]}
    echo ${param[@]}
    value=`echo ${case_group[$caseid]} | awk -F "_" '{print $3}'`
    value=(`echo ${value} | tr ':' ' '`)
    echo ${value[@]}
    baseline_namelist_input="namelists.input$sceneid"


    work_dir=${dida_dir}/${create_name}
    cp $envir_dir/define_case/baseline_namelist/$baseline_namelist_input ${work_dir}/namelists.input
    cd $work_dir


    namelist=namelists.input

    case_name_bak=`grep -i case_name $namelist | awk '{print $3}' | sed $'s/\'//g'`
    case_name_line=`grep -n case_name $namelist | awk -F ":" '{print $1}'`
    sed -i "${case_name_line}s/$case_name_bak/$case_name/g" $namelist


    for x in $( seq 0 `expr ${length_param} - 1` )
    do
        if [ ${param[$x]} == 'start_time' ]
        then
                value_bak=`grep -i ${param[$x]} $namelist |grep -v da| awk '{print $3}' | sed $'s/\'//g'`
                value_line=`grep -n ${param[$x]} $namelist |grep -v da| awk -F ":" '{print $1}'`
        else
                value_bak=`grep -i ${param[$x]} $namelist | awk '{print $3}' | sed $'s/\'//g'`
                value_line=`grep -n ${param[$x]} $namelist | awk -F ":" '{print $1}'`
        fi
        echo $value_line
        sed -i "${value_line}s/$value_bak/${value[$x]}/g" $namelist
    done


}



OBSTAT(){

    obstat_ini=initial_for_obstat.input

    #work_dir=${dida_dir}/${create_name}
    cd $case_work_dir

    cp $obstat_dir/namelist/$obstat_ini .
    namelist=$work_dir/namelists.input

    case_name_bak=`grep -i case_name $obstat_ini | awk '{print $3}' | sed $'s/\'//g'`
    case_name_line=`grep -n case_name $obstat_ini | awk -F ":" '{print $1}'`
    sed -i "${case_name_line}s/$case_name_bak/$case_name/g" $obstat_ini


    obsstat_start_time_bak=`grep -i obsstat_start_time $obstat_ini | awk '{print $3}' | sed $'s/\'//g'`
    obsstat_start_time_case=`grep -i da_start_time $namelist | awk '{print $3}' | sed $'s/\'//g'`
    obsstat_start_time_line=`grep -n obsstat_start_time $obstat_ini | awk -F ":" '{print $1}'`
    sed -i "${obsstat_start_time_line}s/$obsstat_start_time_bak/$obsstat_start_time_case/g" $obstat_ini


    obsstat_end_time_bak=`grep -i obsstat_end_time $obstat_ini | awk '{print $3}' | sed $'s/\'//g'`
    echo $obsstat_end_time_bak
    obsstat_end_time_case=`grep -i end_time $namelist | head -1 | awk '{print $3}' | sed $'s/\'//g'`
    echo $obsstat_end_time_case
    obsstat_end_time_line=`grep -n obsstat_end_time $obstat_ini | awk -F ":" '{print $1}'`
    sed -i "${obsstat_end_time_line}s/$obsstat_end_time_bak/$obsstat_end_time_case/g" $obstat_ini


    obsstat_var_name_bak=`grep -i obsstat_var_name $obstat_ini | awk '{print $3}' | sed $'s/\'//g'`
    obsstat_var_name_case=`grep -i da_var_name $namelist | awk '{print $3}' | sed $'s/\'//g'`
    obsstat_var_name_line=`grep -n obsstat_var_name $obstat_ini | awk -F ":" '{print $1}'`
    sed -i "${obsstat_var_name_line}s/$obsstat_var_name_bak/$obsstat_var_name_case/g" $obstat_ini

    mpas_or_gmcore_bak=`grep -i mpas_or_gmcore $obstat_ini | awk '{print $3}' | sed $'s/\'//g'`
    mpas_or_gmcore_line=`grep -n  mpas_or_gmcore $obstat_ini | awk -F ":" '{print $1}'`
    sed -i "${mpas_or_gmcore_line}s/$mpas_or_gmcore_bak/$mpas_or_gmcore/g" $obstat_ini

    redis_address_bak=`grep -i redis_address $obstat_ini | awk '{print $3}' | sed $'s/\'//g'`
    redis_address_case=`grep -i redis_address $namelist | awk '{print $3}' | sed $'s/\'//g'`
    redis_address_line=`grep -n redis_address $obstat_ini | awk -F ":" '{print $1}'`
    sed -i "${redis_address_line}s/$redis_address_bak/$redis_address_case/g" $obstat_ini

    redis_node_case=`grep -i redis_address $namelist | awk '{print $3}' | sed $'s/\'//g' | cut -d : -f 1`
    redis_port_case=`grep -i redis_address $namelist | awk '{print $3}' | sed $'s/\'//g' | cut -d : -f 2 | cut -d , -f 1`

    redis_node_bak=`grep -i redis_node $obstat_ini | awk '{print $3}' | sed $'s/\'//g'`
    redis_node_line=`grep -n redis_node $obstat_ini | awk -F ":" '{print $1}'`
    sed -i "${redis_node_line}s/$redis_node_bak/$redis_node_case/g" $obstat_ini

    redis_port_bak=`grep -i redis_port $obstat_ini | awk '{print $3}' | sed $'s/\'//g'`
    redis_port_line=`grep -n redis_port $obstat_ini | awk -F ":" '{print $1}'`
    sed -i "${redis_port_line}s/$redis_port_bak/$redis_port_case/g" $obstat_ini


    sed -i "s#/share2/uq/ssp/#$case_work_dir/#g" initial_for_obstat.input

    #cp ${dida_dir}/OBS/OBSTAT/src/run.sh .
    cp ${dida_dir}/OBS/OBSTAT/src/run/run_level_bak.sh run_level.sh
    cp ${dida_dir}/OBS/OBSTAT/src/submit.sh .
    #sed -i "s#./main/#$obstat_dir/main/#g" run.sh
    sed -i "s#./run/#./#g" submit.sh
    #sed -i "s#namepath='./namelist/initial_for_obstat.input'#namepath='./initial_for_obstat.input'#g" run.sh
    sed -i "s#namepath='./namelist/initial_for_obstat.input'#namepath='./initial_for_obstat.input'#g" submit.sh
    sed -i "s#../main/#${dida_dir}/OBS/OBSTAT/src/main/#g" run_level.sh
    sed -i "s/SBATCH -J level/SBATCH -J yuz$case_name/g" run_level.sh

}



CREATE

if [ $nocompile_or_compile == 1 ]
then
   COMPILE
fi


cd $envir_dir
num_line=`grep -c "" case_group.ini`

for x in $( seq 1 `expr ${num_line} - 0` )
do
    cd $envir_dir
    line_message=`sed -n "${x}p" case_group.ini`
    group_cut=(`echo $line_message | tr '=' ' '`)
    group_cut=${group_cut[0]}
    echo $group_cut
    flag=";"
    if [[ $group_cut != *$flag* ]]; then
        group=$group_cut
        case_group=`grep -i $group case_group.ini | awk '{print $3}' | sed $'s/\'//g'`
        case_group=(`echo $case_group | tr ';' ' '`)
        length_case=`expr ${#case_group[@]} - 1`
        sceneid=(`echo $group | tr '_' ' '`)
        sceneid=${sceneid[1]}

        for caseid in $( seq 1 ${length_case})
        do
            case_name=${case_group[$caseid]}
	    case_name=(`echo $case_name | tr '_' ' '`)
	    case_name=${case_name[0]}
            echo $case_name
            case_work_dir=${dida_dir}/${create_name}/${case_name}
            if [[ -d $case_work_dir   ]]; then
               rm -r $case_work_dir
            fi
            mkdir $case_work_dir


            NAMELIST
            OBSTAT

	    #cd $dida_dir/$create_name
	    #sh slurm_run.sh ${case_name} namelists.input

	    cd $envir_dir
	    cp /home/xuewei/lustre/share/nmsc_slurm_example/dida.sbatch ${dida_dir}/${create_name}/.

	    sed -i "s/SBATCH -J dida-test/SBATCH -J yuz$case_name/g" ${dida_dir}/${create_name}/dida.sbatch

            ./run-envir-corun.sh $create_name $case_name $dida_dir

        done

    else
        continue
    fi
done


