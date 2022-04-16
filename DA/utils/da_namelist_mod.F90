module da_namelist_mod
  implicit none

  !! share
  !character(256) :: atm_mpas_sceneid = 'N/A'
  !character(256) :: case_name = 'N/A'
  !character(256) :: redis_address = 'N/A'

  ! character(2) :: atm_mpas_sceneid = '3'
  ! integer :: num_lon=360
  ! integer :: num_lat=180
  ! integer :: da_in_seconds=10800
  ! integer :: da_start_time_array(5) = (/2019,6,7,0,0/)
  
  !integer atm_num_lev
  !integer atm_ensemble

  !integer :: da_mode = 0
  !integer :: da_asynchronous = 0
  !integer :: da_ensemble = 0
  !integer :: da_in_seconds = 10800

  !integer :: da_x_proc_num = 0
  !integer :: da_y_proc_num = 0

  ! da
  real(8) :: da_loc_distance = 0.0
  real(8) :: da_vert_distance = 0.0
  real(8) :: da_inflation_factor = 1.0
  real(8) :: da_polar_letkf_lat = 90.0
  logical :: da_gc_enable = .false.
  integer :: da_obs_sparse = 1
  integer :: da_obs_dense = 0
  integer :: da_block_interval = 4
  integer :: da_vertical_interval = 32
  integer :: da_block_overlap = 2
  integer :: da_halo = 3
  real(8) :: gamma1 = 0
  real(8) :: gamma2 = 1
  integer :: da_kernel = 5
  integer :: da_conv_stride = 3
  integer :: da_pc_f = 9
  integer :: da_pc_p1 = 9
  integer :: da_pc_p2 = 9
  real(8) :: da_proc_f = 0.001
  real(8) :: da_proc_p = 0.999
  integer :: ps_pos = 1
  integer :: computing_precision = 53
  
  ! da_obs
  logical       :: if_regular
  integer       :: obs_hrz_interval_lon
  integer       :: obs_hrz_interval_lat
  integer       :: obs_vtc_interval
  character(3)  :: obs_var_name(3)
  integer       :: mpas_num_lev_start = 13
  integer       :: mpas_num_lev = 60
  integer       :: onepoint_x=1
  integer       :: onepoint_y=1
  integer       :: onepoint_z=1
  character(30) :: obs_sceneid
  logical       :: if_sparse   ! 是否稀疏化观测
  logical       :: if_aircar
  logical       :: if_aircft
  logical       :: if_satwnd
  logical       :: if_adpupa
  logical       :: if_gpsro
  logical       :: if_sfcshp
  logical       :: if_adpsfc
  logical       :: polar_discard ! 是否丢弃极区观测
  integer       :: adpupa_num_lev = 32
  integer       :: gpsro_num_lev = 32
  integer       :: obs_stat_interval_lon = 10
  integer       :: obs_stat_interval_lat = 10
  integer       :: obs_stat_interval_vertical = 4
  character(256):: randn_filepath= 'N/A'
  
  namelist /namelist_da/ &
    da_loc_distance, da_vert_distance, da_inflation_factor, da_gc_enable, da_obs_dense, &
    da_obs_sparse, da_polar_letkf_lat, da_block_interval, da_vertical_interval, &
    da_block_overlap, da_halo, gamma1, gamma2, da_kernel, &
    da_conv_stride, da_pc_f, da_pc_p1, da_pc_p2, da_proc_f, da_proc_p, computing_precision

  namelist /namelist_obs/ &
    if_regular, obs_hrz_interval_lon, obs_hrz_interval_lat, obs_vtc_interval, &
    obs_var_name, mpas_num_lev_start, mpas_num_lev, onepoint_x, onepoint_y, onepoint_z, obs_sceneid, &
    if_sparse, if_aircar, if_aircft, if_satwnd, if_adpupa, if_gpsro, if_sfcshp, if_adpsfc, &
    polar_discard, adpupa_num_lev, gpsro_num_lev, obs_stat_interval_lon, obs_stat_interval_lat, obs_stat_interval_vertical, randn_filepath

contains

  subroutine parse_namelist(file_path)

    character(*), intent(in) :: file_path

    open(10, file=file_path, status='old')
    rewind(10)
    read(10, nml=namelist_da)
    rewind(10)
    read(10, nml=namelist_obs)
    close(10)

  end subroutine parse_namelist

end module da_namelist_mod
