module da_config
  implicit none

  !namelist total
   character(30)  :: da_varname
   real*4         :: da_loc_distance
   integer        :: da_inflation_on
   real*4         :: da_inflation_factor
   real*4         :: da_polar_letkf_lat
   integer        :: da_block_interval
   integer        :: da_block_overlap
   integer        :: da_var_num
   integer        :: da_mode
   integer        :: da_ensemble
   integer        :: da_kernel
   integer        :: da_pc_f
   integer        :: da_pc_p1
   integer        :: da_pc_p2
   real*4         :: da_proc_f
   real*4         :: da_proc_p

contains
  subroutine init_config_da()
    use coupler_config
    use redis_mod

    implicit none

    type(c_ptr) :: rds_cc

    rds_cc = RedisConnect(redis_address)

    call RedisHget(rds_cc, case_name, 'da_varname', da_varname)
    call RedisHget(rds_cc, case_name, 'da_loc_distance', da_loc_distance)
    call RedisHget(rds_cc, case_name, 'da_inflation_on', da_inflation_on)
    call RedisHget(rds_cc, case_name, 'da_inflation_factor', da_inflation_factor)
    call RedisHget(rds_cc, case_name, 'da_polar_letkf_lat', da_polar_letkf_lat)
    call RedisHget(rds_cc, case_name, 'da_block_interval', da_block_interval)
    call RedisHget(rds_cc, case_name, 'da_block_overlap', da_block_overlap)
    call RedisHget(rds_cc, case_name, 'da_var_num', da_var_num)
    call RedisHget(rds_cc, case_name, 'da_mode', da_mode)
    call RedisHget(rds_cc, case_name, 'da_ensemble', da_ensemble)
    call RedisHget(rds_cc, case_name, 'da_kernel', da_kernel)
    call RedisHget(rds_cc, case_name, 'da_pc_f', da_pc_f)
    call RedisHget(rds_cc, case_name, 'da_pc_p1', da_pc_p1)
    call RedisHget(rds_cc, case_name, 'da_pc_p2', da_pc_p2)
    call RedisHget(rds_cc, case_name, 'da_proc_f', da_proc_f)
    call RedisHget(rds_cc, case_name, 'da_proc_p', da_proc_p)

  end subroutine
end module da_config 

