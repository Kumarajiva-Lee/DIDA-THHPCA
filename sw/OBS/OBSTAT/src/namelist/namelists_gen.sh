#!/bin/bash
#temp var
export redis_node='192.168.1.1'

#total
export case_name='case_0_1_1'                       # string  例如 case_0_1_1 试验名称
export redis_address='192.168.1.1:6379,192.168.1.2:6379,192.168.1.3:6379'
export atm_mpas_sceneid=1                           # integer 1/2/3/4 四个mpas初始场起始时间场景
export num_lon=360                              # integer 例如 360    模式水平分辨率（纬向）
export num_lat=180                              # integer 例如 180    模式水平分辨率（经向）
export num_lev=26                               # integer 例如 64 模式垂直分辨率
                                                    # string  例如 node0：6379；node2：6379   redis集群的连接地址
export atm_mpas2atm_time='1 days'                   # string  '1 days'；'3 days'；‘5 days’    从mpas初始场到gmcore试验开始时间
export atm_phase_diff=10800                         # integer 例如 10800（单位：秒）  选取样本相差
export da_endnn_in_seconds=30                       # integer 例如 30（单位：秒） 在线学习步长
export da_start_time='2000,01,01,00,00'                       # string 同化起始时间
export end_time='2000,01,30,00,00'                  #同化结束时间
export da_in_seconds=10800                          # integer 例如 10800（单位：秒）  同化步长
export atm_ensemble_total=1                        #gmcore总集合数
export atm_ensemble_group=1                         #gmcore分组数
export da_ensemble=1                               # integer 16  同化集合数. optional value: 16, 32, 64, 128, 256, 512, 1024;
export da_var_name='u,v,t,ps'                     # String array    'u,v,t,ps,gz'    变量类型名称
export da_asynchronous=0                            # integer 0：不开启；1：开启  异步同化开关
export da_mode=0                                    # integer 0   同化类型: 0-在线学习; 1-LETKF同化; 2-集合代理LETKF; 3-同化代理;
export atm_group=12                                  # integer 例如 4  atm组进程数
export atm_group_num=1                              # integer 例如 10 atm进程组数
export atm_stride=16                                 # integer 例如 6  atm进程跨步
export atm_root=0                                   # integer 例如 0  atm开始进程
export da_group=4                                   # integer 例如 2  da组进程数
export da_group_num=1                               # integer 例如 10 da进程组数
export da_stride=16                                  # integer 例如 6  da进程跨步
export da_root=12                                    # integer 例如 4  da开始进程

#atm
export num_proc_lon=2                            # integer gmcore经向进程划分
export num_proc_lat=6                            # integer gmcore纬向进程划分
export test_case='rh4'                           # string gmcore标准算例名称
export run_days=1                                # 运行时间
export time_scheme='pc2'                         # string 'pc2' 时间差分格式
export dt_in_seconds=90                          # gmcore 积分步长
export history_interval='3 hours'                # string history 文件输出频率
export reduce_factors='20, 20, 20, 15, 15, 15, 12, 12, 10, 10, 10, 10, 6, 6, 4, 4, 2, 2'        # integer 例如 8, 8, 8, 6, 6, 6, 4, 4, 2  简并因子1
export use_div_damp=.false.
export use_vor_damp=.false.
export use_polar_damp=.false.
export try_read=4122
export job_proc_lat=20,30,40,40,30,20               # 新增
export member_num=1

#hybrid_coord
export template='test_l26'
#da
export da_loc_distance=800.0                         # float   例如 800.0d0    水平局地化半径
export da_inflation_factor=1.0                       # float   1.1；1.3；1.5；1.7  膨胀系数，大于1小于2的数
export da_polar_letkf_lat=80.0                       # float   例如 80.0d0 极区纬度（实施集合方差滤波）
export da_block_interval=4                           # integer 例如 4或8   letkf分块大小
export da_block_overlap=2                            # integer 例如 2或4   letkf重叠宽度
export da_halo=3                                     # integer 3 da模块halo区宽度
export da_kernel=5                                   # integer 5   水平卷积核边长
export da_conv_stride=3                              # integer 3   卷积核跳转步长
export da_pc_f=9                                     # integer 9   特征层主成分数量
export da_pc_p1=9                                    # integer 9   预处理层1主成分数量
export da_pc_p2=9                                    # integer 9   预处理层2主成分数量
export da_proc_f=0.001                               # float   0.001   特征层进程比例
export da_proc_p=0.999                               # float   0.999   预处理层进程比例

#obs
export obs_hrz_interval_lon=10                       # integer 例如 10 （0表示单点试验）   水平观测间隔（纬向）
export obs_hrz_interval_lat=10                       # integer 例如 10 （0表示单点试验）   水平观测间隔（经向）
export obs_vtc_interval=1                            # integer 例如 4（0表示单点试验） 垂直观测间隔
export obs_var_name='u,v,t,ps'                    # character 例如'u,v,t,ps,gz'
export obs_stat_interval_lon=1                      # integer 例如 10 计算网格纬向间隔几个gmocre网格点
export obs_stat_interval_lat=1                      # integer 例如 10 计算网格经向间隔几个gmcore网格点
export obs_stat_interval_vertical=4                  # integer 例如 4  计算网格垂向间隔几个gmcore模式层
export mpas_num_lev=26                             # integer 例如60 计算网格垂向层数
export obsstat_start_time='2000,01,05,00,00' 
export obsstat_end_time='2000,01,05,09,00'
export obsstat_var_name='u,v'
export draw_level=1
export draw_single=1
export draw_global=1
export obs_open=0


while getopts :c:r:l:s:h opt
do
    case $opt in
        c)
            case_name=$OPTARG
        ;;
        r)
            num_lon=`echo $OPTARG | cut -d 'x' -f 1`
            num_lat=`echo $OPTARG | cut -d 'x' -f 2`
        ;;
        l)
            num_lev=$OPTARG
        ;;
        s)
            redis_address=$OPTARG
        ;;
        e)
            redis_address=$OPTARG
        ;;
        h)
            echo "-c        case_name, exp: test_0_1_1"
            echo "-s        redis_address, exp: 192.168.1.1"
            echo "-r        reslution, exp: 360x180"
            echo "-l        num_lev, exp: 64"
            echo "-h        help"
            exit 1;
        ;;
        ?)
        echo "Unknown parameter $OPTARG"
        exit 1;;
    esac
done

############################################################
cat > ./namelists.input << EOF
 &namelist_total
   case_name                 = '$case_name'
   redis_address             = '$redis_address'
   atm_mpas_sceneid          = $atm_mpas_sceneid
   num_lon                   = $num_lon
   num_lat                   = $num_lat
   num_lev                   = $num_lev
   atm_mpas2atm_time         = '$atm_mpas2atm_time'
   atm_phase_diff            = $atm_phase_diff
   da_endnn_in_seconds       = $da_endnn_in_seconds
   da_start_time             = $da_start_time
   end_time                  = $end_time
   da_in_seconds             = $da_in_seconds
   atm_ensemble_total        = $atm_ensemble_total
   atm_ensemble_group        = $atm_ensemble_group
   da_ensemble               = $da_ensemble
   da_var_name               = '$da_var_name'
   da_asynchronous           = $da_asynchronous
   da_mode                   = $da_mode
   atm_group                 = $atm_group
   atm_group_num             = $atm_group_num
   atm_stride                = $atm_stride
   atm_root                  = $atm_root
   da_group                  = $da_group
   da_group_num              = $da_group_num
   da_stride                 = $da_stride
   da_root                   = $da_root
/
 &namelist_atm
   num_proc_lon              = $num_proc_lon
   num_proc_lat              = $num_proc_lat
   test_case                 = '$test_case'
   run_days                  = $run_days
   time_scheme               = '$time_scheme'
   dt_in_seconds             = $dt_in_seconds
   history_interval          = '$history_interval'
   reduce_factors            = $reduce_factors
   use_div_damp              = $use_div_damp
   use_vor_damp              = $use_vor_damp
   use_polar_damp            = $use_polar_damp
   try_read                  = $try_read
   job_proc_lat              = $job_proc_lat
   member_num                = $member_num
/
 &hybrid_coord
   template          = 'test_l26'
/
 &namelist_da
   da_loc_distance           = $da_loc_distance
   da_inflation_factor       = $da_inflation_factor
   da_polar_letkf_lat        = $da_polar_letkf_lat
   da_block_interval         = $da_block_interval
   da_block_overlap          = $da_block_overlap
   da_halo                   = $da_halo
   da_kernel                 = $da_kernel
   da_conv_stride            = $da_conv_stride
   da_pc_f                   = $da_pc_f
   da_pc_p1                  = $da_pc_p1
   da_pc_p2                  = $da_pc_p2
   da_proc_f                 = $da_proc_f
   da_proc_p                 = $da_proc_p
/
 &namelist_obs
   obs_hrz_interval_lon      = $obs_hrz_interval_lon
   obs_hrz_interval_lat      = $obs_hrz_interval_lat
   obs_vtc_interval          = $obs_vtc_interval
   obs_var_name              = '$obs_var_name'
   obs_stat_interval_lon     = $obs_stat_interval_lon
   obs_stat_interval_lat     = $obs_stat_interval_lat
   obs_stat_interval_vertical= $obs_stat_interval_vertical
   mpas_num_lev              = $mpas_num_lev
/
EOF
########################
cat > ./initial_for_obstat.input << EOF
 [obstat_ini]
   case_name                 = '$case_name'
   redis_address             = '$redis_address'
   redis_node                = '$redis_node'
   obsstat_start_time        = '$obsstat_start_time'
   obsstat_end_time          = '$obsstat_end_time'
   obsstat_var_name          = '$obsstat_var_name'
   draw_level                = $draw_level
   draw_single               = $draw_single
   draw_global               = $draw_global
   obs_open                  = $obs_open

EOF
#######################
#IFS=$'\n'
#for conf in `grep '=' ./namelists.input`
#do
    #redis-cli -h $redis_node -c hset $case_name \
        #`echo $conf | cut -d '=' -f 1 | awk '{gsub(/^\s+|\s+$/, "");print}'` \
        #`echo $conf | cut -d '=' -f 2 | awk '{gsub(/^\s+|\s+$/, "");print}'` \
         #>/dev/null 2>&1
#done
