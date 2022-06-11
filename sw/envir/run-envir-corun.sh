
SUBMIT(){
if [ ! -z $JobCmd -a ! -z $JobFile -a -e $JobCmd  ];then
    sh $JobCmd $case_name $namelist > $JobFile
    if [ $? -eq 0  ];then
        sleep 4
        if [ -e $JobFile  ]; then
            job_num=`tail -1 $JobFile  | cut -d " " -f 4`
            cd $work_dir
            ./submit.sh $picpath $input_path $job_num >& obstat.log
        fi
    fi
fi
}


create_name=$1
case_name=$2
dida_dir=$3
namelist=namelists.input
obstat_ini=initial_for_obstat.input
obstat_dir=$dida_dir/OBS/OBSTAT/src
create_dir=${dida_dir}/${create_name}
work_dir=${dida_dir}/${create_name}/${case_name}
picpath=${work_dir}/pic
input_path=${work_dir}/${obstat_ini}


cd $create_dir

export JobFile=`pwd`/jobid.log
export JobCmd=`pwd`/slurm_run.sh

SUBMIT
iret=$?



