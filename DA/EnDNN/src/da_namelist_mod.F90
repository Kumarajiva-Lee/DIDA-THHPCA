module da_namelist_mod
  implicit none

  !! share
  !character(256) :: atm_mpas_sceneid = 'N/A'
  !character(256) :: case_name = 'N/A'
  !character(256) :: redis_address = 'N/A'

  !integer atm_num_lon
  !integer atm_num_lat
  !integer atm_num_lev
  !integer atm_ensemble

  !integer :: da_mode = 0
  !integer :: da_asynchronous = 0
  !integer :: da_ensemble = 0
  !integer :: da_in_seconds = 10800

  !integer :: da_x_proc_num = 0
  !integer :: da_y_proc_num = 0

  ! da
  character(20) :: da_var_name = 'N/A'
  real(8) :: da_loc_distance = 0.0
  real(8) :: da_inflation_factor = 1.0
  real(8) :: da_polar_letkf_lat = 90.0
  integer :: da_block_interval = 4
  integer :: da_block_overlap = 2
  integer :: da_halo = 3
  integer :: da_kernel = 5
  integer :: da_conv_stride = 3
  integer :: da_pc_f = 9
  integer :: da_pc_p1 = 9
  integer :: da_pc_p2 = 9
  real(8) :: da_proc_f = 0.001
  real(8) :: da_proc_p = 0.999
  
  ! da_obs
  character(20) :: obs_var_name = 'N/A'
  integer :: obs_hrz_interval_lon = 10
  integer :: obs_hrz_interval_lat = 10
  integer :: obs_vtc_interval = 4
  integer :: obs_stat_interval_lon = 10
  integer :: obs_stat_interval_lat = 10
  integer :: obs_stat_interval_vertical = 4
  integer :: mpas_num_lev = 60
  
  namelist /namelist_da/ &
    da_var_name, da_loc_distance, da_inflation_factor, &
    da_polar_letkf_lat, da_block_interval, &
    da_block_overlap, da_halo, da_kernel, &
    da_conv_stride, da_pc_f, da_pc_p1, da_pc_p2, da_proc_f, da_proc_p 

  namelist /namelist_obs/ &
    obs_hrz_interval_lon, obs_hrz_interval_lat, obs_vtc_interval, obs_var_name, &
    obs_stat_interval_lon, obs_stat_interval_lat, obs_stat_interval_vertical, mpas_num_lev

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
