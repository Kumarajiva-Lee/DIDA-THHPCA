program unit_test
  use redis_mod
  use da_namelist_mod
  use read_obs_mod

  implicit none
  integer lats, late, lons, lone
  integer ntimes
  real(4), allocatable :: lat_grid(:), lon_grid(:)
  type(c_ptr) rc
  type(obs_list), allocatable :: obsinfo(:)
  character(256)       :: namelist_path
  integer :: n

  ! read namelist
  call get_command_argument(1, namelist_path)
  call parse_namelist(namelist_path)

  ! redis connect
  rc = RedisConnect('172.16.81.1:6379,172.16.81.2:6379')
  
  lats = 1
  late = 180
  lons = 1
  lone = 360
  ntimes = 0

  allocate(lat_grid(late-lats+1), lon_grid(lone-lons+1))
  call read_model_grid(rc, lats, late, lons, lone, lat_grid, lon_grid)
  call obs_init
  allocate(obsinfo(obs_type_num))
  n = 0
  if (if_aircar) then
    n = n + 1
    call allocate_obslist(rc, lat_grid(1), lat_grid(late-lats+1), lon_grid(1), lon_grid(lone-lons+1), ntimes, 'aircar', obsinfo(n))
    call read_obs(rc, ntimes, lat_grid(1), lat_grid(late-lats+1), lon_grid(1), lon_grid(lone-lons+1), obsinfo(n))
print *, obsinfo(n)%ph(1,1), obsinfo(n)%yo(:,1,1),  obsinfo(n)%obs_rinv(1,1,1)
  endif
  if (if_aircft) then
    n = n + 1
    call allocate_obslist(rc, lat_grid(1), lat_grid(late-lats+1), lon_grid(1), lon_grid(lone-lons+1), ntimes, 'aircft', obsinfo(n))
    call read_obs(rc, ntimes, lat_grid(1), lat_grid(late-lats+1), lon_grid(1), lon_grid(lone-lons+1), obsinfo(n))
print *, obsinfo(n)%ph(1,1), obsinfo(n)%yo(:,1,1),  obsinfo(n)%obs_rinv(1,1,1)
  endif
  if (if_satwnd) then
    n = n + 1
    call allocate_obslist(rc, lat_grid(1), lat_grid(late-lats+1), lon_grid(1), lon_grid(lone-lons+1), ntimes, 'satwnd', obsinfo(n))
    call read_obs(rc, ntimes, lat_grid(1), lat_grid(late-lats+1), lon_grid(1), lon_grid(lone-lons+1), obsinfo(n))
print *, obsinfo(n)%ph(1,1), obsinfo(n)%yo(:,1,1),  obsinfo(n)%obs_rinv(1,1,1)
  endif
  if (if_adpupa) then
    n = n + 1
    call allocate_obslist(rc, lat_grid(1), lat_grid(late-lats+1), lon_grid(1), lon_grid(lone-lons+1), ntimes, 'adpupa', obsinfo(n))
    call read_obs(rc, ntimes, lat_grid(1), lat_grid(late-lats+1), lon_grid(1), lon_grid(lone-lons+1), obsinfo(n))
print *, obsinfo(n)%ph(1,1), obsinfo(n)%yo(:,1,1),  obsinfo(n)%obs_rinv(1,1,1)
  endif
  if (if_gpsro) then
    n = n + 1
    call allocate_obslist(rc, lat_grid(1), lat_grid(late-lats+1), lon_grid(1), lon_grid(lone-lons+1), ntimes, 'gpsro', obsinfo(n))
    call read_obs(rc, ntimes, lat_grid(1), lat_grid(late-lats+1), lon_grid(1), lon_grid(lone-lons+1), obsinfo(n))
print *, obsinfo(n)%ph(1,1), obsinfo(n)%yo(:,1,1),  obsinfo(n)%obs_rinv(1,1,1)
  endif
  if (if_sfcshp) then
    n = n + 1
    call allocate_obslist(rc, lat_grid(1), lat_grid(late-lats+1), lon_grid(1), lon_grid(lone-lons+1), ntimes, 'sfcshp', obsinfo(n))
    call read_obs(rc, ntimes, lat_grid(1), lat_grid(late-lats+1), lon_grid(1), lon_grid(lone-lons+1), obsinfo(n))
print *, obsinfo(n)%ph(1,1), obsinfo(n)%yo(:,1,1),  obsinfo(n)%obs_rinv(1,1,1)
  endif
  if (if_adpsfc) then
    n = n + 1
    call allocate_obslist(rc, lat_grid(1), lat_grid(late-lats+1), lon_grid(1), lon_grid(lone-lons+1), ntimes, 'adpsfc', obsinfo(n))
    call read_obs(rc, ntimes, lat_grid(1), lat_grid(late-lats+1), lon_grid(1), lon_grid(lone-lons+1), obsinfo(n))
print *, obsinfo(n)%ph(1,1), obsinfo(n)%yo(:,1,1),  obsinfo(n)%obs_rinv(1,1,1)
  endif
      
  call obs_final(obsinfo)
  deallocate(lat_grid, lon_grid, obsinfo)
end program unit_test
