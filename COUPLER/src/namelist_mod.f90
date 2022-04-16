module namelist_mod
  use mpi
  implicit none

  !namelist total
  character(256) :: redis_address
  character(30)  :: case_name
  character(30)  :: atm_mpas2atm_time
  character(30)  :: atm_run_time
  character(30)  :: atm2da_spinup_time
  integer :: atm_num_lon
  integer :: atm_num_lat
  integer :: atm_num_lev
  integer :: atm_mpas_sceneid
  integer :: atm_phase_diff
  integer :: atm_ensemble
  integer :: da_in_seconds
  integer :: da_asynchronous
  integer :: da_mode
  integer :: da_ensemble
  integer :: da_endnn_in_seconds
  integer :: atm_group
  integer :: da_group
  integer :: atm_group_num
  integer :: da_group_num
  integer :: atm_stride
  integer :: da_stride
  integer :: atm_root
  integer :: da_root
  integer :: job_proc_lat(100)
  integer :: member_num

  !statements that specify the namelists
  namelist /namelist_total/ &
    redis_address         , &
    case_name             , &
    atm_mpas2atm_time     , &
    atm_run_time          , &
    atm2da_spinup_time    , &
    atm_num_lon           , &
    atm_num_lat           , &
    atm_num_lev           , &
    atm_mpas_sceneid      , &
    atm_phase_diff        , &
    atm_ensemble          , &
    da_in_seconds         , &
    da_asynchronous       , &
    da_mode               , &
    da_ensemble           , &
    da_endnn_in_seconds   , &
    atm_group             , &
    da_group              , &
    atm_group_num         , &
    da_group_num          , &
    atm_stride            , &
    da_stride             , &
    atm_root              , &
    da_root               , &
    job_proc_lat         , &
    member_num

  namelist /namelist_atm/ &
   atm_split_scheme           , &
   atm_dt_in_seconds          , &
   atm_reduce_factors         , &
   atm_damp_orders            , &
   atm_fast_cycles            , &
   atm_pv_scheme              , &
   atm_coarse_polar_lats_exp  , &
   atm_coarse_polar_lats_mul  

  namelist /namelist_da/ &
   da_loc_distance           , &
   da_inflation_factor       , &
   da_polar_letkf_lat        , &
   da_block_interval         , &
   da_block_overlap          , &
   da_var_name               , &
   da_mode                   , &
   da_ensemble               , &
   da_kernel                 , &
   da_pc_f                   , &
   da_pc_p1                  , &
   da_pc_p2                  , &
   da_proc_f                 , &
   da_proc_p                 

  namelist /namelist_obs/ &
   obs_hrz_interval_lon      , &
   obs_hrz_interval_lat      , &
   obs_vtc_interval          , &
   obs_var_name              , &
   obs_stat_interval_lon     , &
   obs_stat_interval_lat     , &
   obs_stat_interval_vertical, &
   mpas_num_lev

end module coupler_config 
