#include "../../utils/da_config.inc"
#define SINGLE_POINT 1

module read_obs_mod

  use coupler_config, da_start_time_array => da_start_time
  use da_namelist_mod
  use redis_mod
  use datetime
  use randn_mod
  use flogger

  implicit none

  type :: obs_list
    character*10         :: obs_type        ! obs type
    integer              :: obs_num         ! number of observation
    integer              :: obs_nz          ! number of obs level
    integer              :: obs_nvar        ! number of obs variable
    character(3),allocatable :: obs_varname(:)  ! variable name of obs; obs_varname(obs_nvar)
    integer, allocatable     :: obs_npos(:)     ! obs start position in background n; obs_npos(obs_nvar)
    real(8), allocatable     :: lat(:), lon(:)  ! latitude and longitude of obs; lat/lon(obs_num)
    real(8), allocatable     :: ph(:,:)         ! pressure of obs; ph(obs_nz, obs_num)
    ayb_type, allocatable    :: dis_yo(:,:,:)   ! obs variable value(obs_nvar, obs_nz, obs_num)
    yb_type, allocatable     :: pro_yb(:,:,:,:) ! obs objection value(ens, obs_nvar, obs_nz, obs_num)
    real(8), allocatable     :: obs_rinv(:,:,:) ! inverse of the observation error; R^-1(obs_nvar, obs_nz, obs_num)
  end type

  integer :: obs_type_num               ! number of obs type
  integer :: ens

  integer :: obs_type_num_max = 7
  character(len = 10), dimension(:), allocatable :: obs_type_list
  integer :: obs_vertex_id(4)

  private ens, obs_type_num_max, obs_vertex_id

contains
  subroutine obs_init(ens_in)
    implicit none
    integer :: n
    integer, intent(in) :: ens_in

    ens = ens_in

    allocate(obs_type_list(obs_type_num_max))

    if (if_regular) then
      obs_type_num = 1
      obs_type_list(1) = 'REGULAR' ! regular distribution
      return
    endif

    n = 0
    if (if_aircar) n = n + 1
    if (if_aircar) obs_type_list(n) = 'AIRCAR'
    if (if_aircft) n = n + 1
    if (if_aircft) obs_type_list(n) = 'AIRCFT'
    if (if_satwnd) n = n + 1
    if (if_satwnd) obs_type_list(n) = 'SATWND'
    if (if_adpupa) n = n + 1
    if (if_adpupa) obs_type_list(n) = 'ADPUPA'
    if (if_gpsro ) n = n + 1
    if (if_gpsro ) obs_type_list(n) = 'GPSRO'
    if (if_sfcshp) n = n + 1
    if (if_sfcshp) obs_type_list(n) = 'SFCSHP'
    if (if_adpsfc) n = n + 1
    if (if_adpsfc) obs_type_list(n) = 'ADPSFC'

    obs_type_num = n
  end subroutine obs_init

  subroutine allocate_obslist(rc, lats_grid, late_grid, lons_grid, lone_grid, ntimes, obstypename, obsinfo)
    type(c_ptr), intent(in)        :: rc
    real(4)    , intent(in)        :: lats_grid, late_grid, lons_grid, lone_grid
    integer, intent(in)            :: ntimes
    character(*), intent(in)       :: obstypename
    type(obs_list), intent(inout)  :: obsinfo
    character(:), allocatable      :: hkey, vkey
    character(12)                  :: date_char
    integer                        :: obs_num
    type(datetime_type)              da_start_time
    type(datetime_type)              curr_time

    obsinfo%obs_type = obstypename

    ! calculate date
    da_start_time = create_datetime(year=da_start_time_array(1),  &
                                   month=da_start_time_array(2), &
                                   day=da_start_time_array(3),   &
                                   hour=da_start_time_array(4),  &
                                   minute=da_start_time_array(5))
    curr_time = da_start_time+create_timedelta(seconds=da_in_seconds*ntimes)

    date_char   = curr_time%format('%Y%m%d%H%M')

    if (trim(adjustl(obstypename)) == 'AIRCAR') then
      obsinfo%obs_nz   = 1
      obsinfo%obs_nvar = 3
      allocate(obsinfo%obs_varname(obsinfo%obs_nvar))
      obsinfo%obs_varname(1) = 'u'
      obsinfo%obs_varname(2) = 'v'
      obsinfo%obs_varname(3) = 'pt'
    elseif (trim(adjustl(obstypename)) == 'AIRCFT') then
      obsinfo%obs_nz   = 1
      obsinfo%obs_nvar = 3
      allocate(obsinfo%obs_varname(obsinfo%obs_nvar))
      obsinfo%obs_varname(1) = 'u'
      obsinfo%obs_varname(2) = 'v'
      obsinfo%obs_varname(3) = 'pt'
    elseif (trim(adjustl(obstypename)) == 'SATWND') then
      obsinfo%obs_nz   = 1
      obsinfo%obs_nvar = 2
      allocate(obsinfo%obs_varname(obsinfo%obs_nvar))
      obsinfo%obs_varname(1) = 'u'
      obsinfo%obs_varname(2) = 'v'
    elseif (trim(adjustl(obstypename)) == 'ADPUPA') then
      obsinfo%obs_nz   = adpupa_num_lev
      obsinfo%obs_nvar = 3
      allocate(obsinfo%obs_varname(obsinfo%obs_nvar))
      obsinfo%obs_varname(1) = 'u'
      obsinfo%obs_varname(2) = 'v'
      obsinfo%obs_varname(3) = 'pt'
    elseif (trim(adjustl(obstypename)) == 'GPSRO') then
      obsinfo%obs_nz   = gpsro_num_lev
      obsinfo%obs_nvar = 1
      allocate(obsinfo%obs_varname(obsinfo%obs_nvar))
      obsinfo%obs_varname(1) = 'pt'
    elseif (trim(adjustl(obstypename)) == 'SFCSHP') then
      obsinfo%obs_nz   = 1
      obsinfo%obs_nvar = 3
      allocate(obsinfo%obs_varname(obsinfo%obs_nvar))
      obsinfo%obs_varname(1) = 'u'
      obsinfo%obs_varname(2) = 'v'
      obsinfo%obs_varname(3) = 'pt'
    elseif (trim(adjustl(obstypename)) == 'ADPSFC') then
      obsinfo%obs_nz   = 1
      obsinfo%obs_nvar = 3
      allocate(obsinfo%obs_varname(obsinfo%obs_nvar))
      obsinfo%obs_varname(1) = 'u'
      obsinfo%obs_varname(2) = 'v'
      obsinfo%obs_varname(3) = 'pt'
    endif


    ! get number of obs in the block
    hkey = "obs:" // trim(adjustl(obs_sceneid)) // ":" // trim(adjustl(obstypename)) //":" // trim(adjustl(date_char))
    call RedisObscount(rc, trim(hkey), lons_grid, lone_grid, lats_grid, late_grid, obs_num)
    if (trim(adjustl(obstypename)) == 'ADPUPA') then
      obsinfo%obs_num = obs_num/adpupa_num_lev
    elseif (trim(adjustl(obstypename)) == 'GPSRO') then
      obsinfo%obs_num = obs_num/gpsro_num_lev
    else
      obsinfo%obs_num = obs_num
    endif

    allocate(obsinfo%obs_npos(obsinfo%obs_nvar))
    allocate(obsinfo%lat(obsinfo%obs_num))
    allocate(obsinfo%lon(obsinfo%obs_num))
    allocate(obsinfo%ph(obsinfo%obs_nz, obsinfo%obs_num))
    allocate(obsinfo%dis_yo(obsinfo%obs_nvar, obsinfo%obs_nz, obsinfo%obs_num))
    allocate(obsinfo%pro_yb(ens, obsinfo%obs_nvar, obsinfo%obs_nz, obsinfo%obs_num))
    allocate(obsinfo%obs_rinv(obsinfo%obs_nvar, obsinfo%obs_nz, obsinfo%obs_num))

  end subroutine allocate_obslist

  subroutine allocate_obslist_regular(rc, lats, late, lons, lone, ntimes, obstypename, obsinfo)
    type(c_ptr), intent(in)        :: rc
    integer    , intent(in)        :: lats, late, lons, lone
    integer, intent(in)            :: ntimes
    character(*), intent(in)       :: obstypename
    type(obs_list), intent(inout)  :: obsinfo
    integer                        :: obs_lon_num, obs_lat_num
    integer                        :: n

    obsinfo%obs_type = obstypename
#if (SINGLE_POINT == 0)
    obsinfo%obs_nz   = (mpas_num_lev - mpas_num_lev_start) / obs_vtc_interval + 1
#elif (SINGLE_POINT == 1)
    obsinfo%obs_nz   = 1
#endif

    obsinfo%obs_nvar = 0
    do n = 1, 4
      if (obs_var_name(n) == 'pt' .or. &
          obs_var_name(n) == 'u'  .or. &
          obs_var_name(n) == 'v') then
        obsinfo%obs_nvar = obsinfo%obs_nvar+1
      endif
    enddo
    allocate(obsinfo%obs_varname(obsinfo%obs_nvar))
    do n = 1, obsinfo%obs_nvar
      obsinfo%obs_varname(n) = obs_var_name(n)
    enddo

#if (SINGLE_POINT == 0)
    obs_vertex_id(1) = obs_hrz_interval_lon-mod(lons-1,obs_hrz_interval_lon)+lons
    obs_vertex_id(2) = -mod(lone-1,obs_hrz_interval_lon)+lone
    if (mod(lons-1,obs_hrz_interval_lon)<1) then
        obs_vertex_id(1) = obs_vertex_id(1) - obs_hrz_interval_lon
    endif

    if (mod(lone-1,obs_hrz_interval_lon)<1) then
        obs_vertex_id(2) = obs_vertex_id(2) - obs_hrz_interval_lon
    endif
    obs_lon_num = (obs_vertex_id(2) - obs_vertex_id(1))/obs_hrz_interval_lon+1

    obs_vertex_id(3) = obs_hrz_interval_lat-mod(lats-mod(num_lat/2, obs_hrz_interval_lat),obs_hrz_interval_lat)+lats
    obs_vertex_id(4) = -mod(late-mod(num_lat/2, obs_hrz_interval_lat),obs_hrz_interval_lat)+late
    if (mod(lats-mod(num_lat/2, obs_hrz_interval_lat),obs_hrz_interval_lat)<1) then
        obs_vertex_id(3) = obs_vertex_id(3) - obs_hrz_interval_lat
    endif
    if (mod(late-mod(num_lat/2, obs_hrz_interval_lat),obs_hrz_interval_lat)<1) then
        obs_vertex_id(4) = obs_vertex_id(4) - obs_hrz_interval_lat
    endif
    obs_lat_num = (obs_vertex_id(4) - obs_vertex_id(3))/obs_hrz_interval_lat+1

    obsinfo%obs_num = obs_lon_num*obs_lat_num
#elif (SINGLE_POINT == 1)
    obs_vertex_id(1) = onepoint_x
    obs_vertex_id(2) = onepoint_x
    obs_vertex_id(3) = onepoint_y
    obs_vertex_id(4) = onepoint_y
    obsinfo%obs_num = 1
#endif

    allocate(obsinfo%obs_npos(obsinfo%obs_nvar))
    allocate(obsinfo%lat(obsinfo%obs_num))
    allocate(obsinfo%lon(obsinfo%obs_num))
    allocate(obsinfo%ph(obsinfo%obs_nz, obsinfo%obs_num))
    allocate(obsinfo%dis_yo(obsinfo%obs_nvar, obsinfo%obs_nz, obsinfo%obs_num))
    allocate(obsinfo%pro_yb(ens, obsinfo%obs_nvar, obsinfo%obs_nz, obsinfo%obs_num))
    allocate(obsinfo%obs_rinv(obsinfo%obs_nvar, obsinfo%obs_nz, obsinfo%obs_num))

  end subroutine allocate_obslist_regular

  subroutine allocate_obslist_by_num(obstypename, obs_num_in, obsinfo)

    integer, intent(in) :: obs_num_in
    character(*), intent(in) :: obstypename
    type(obs_list), intent(inout)  :: obsinfo

    obsinfo%obs_num = obs_num_in
    obsinfo%obs_type = obstypename

    if (trim(adjustl(obstypename)) == 'AIRCAR') then
      obsinfo%obs_nz   = 1
      obsinfo%obs_nvar = 3
      allocate(obsinfo%obs_varname(obsinfo%obs_nvar))
      obsinfo%obs_varname(1) = 'u'
      obsinfo%obs_varname(2) = 'v'
      obsinfo%obs_varname(3) = 'pt'
    elseif (trim(adjustl(obstypename)) == 'AIRCFT') then
      obsinfo%obs_nz   = 1
      obsinfo%obs_nvar = 3
      allocate(obsinfo%obs_varname(obsinfo%obs_nvar))
      obsinfo%obs_varname(1) = 'u'
      obsinfo%obs_varname(2) = 'v'
      obsinfo%obs_varname(3) = 'pt'
    elseif (trim(adjustl(obstypename)) == 'SATWND') then
      obsinfo%obs_nz   = 1
      obsinfo%obs_nvar = 2
      allocate(obsinfo%obs_varname(obsinfo%obs_nvar))
      obsinfo%obs_varname(1) = 'u'
      obsinfo%obs_varname(2) = 'v'
    elseif (trim(adjustl(obstypename)) == 'ADPUPA') then
      obsinfo%obs_nz   = adpupa_num_lev
      obsinfo%obs_nvar = 3
      allocate(obsinfo%obs_varname(obsinfo%obs_nvar))
      obsinfo%obs_varname(1) = 'u'
      obsinfo%obs_varname(2) = 'v'
      obsinfo%obs_varname(3) = 'pt'
    elseif (trim(adjustl(obstypename)) == 'GPSRO') then
      obsinfo%obs_nz   = gpsro_num_lev
      obsinfo%obs_nvar = 1
      allocate(obsinfo%obs_varname(obsinfo%obs_nvar))
      obsinfo%obs_varname(1) = 'pt'
    elseif (trim(adjustl(obstypename)) == 'SFCSHP') then
      obsinfo%obs_nz   = 1
      obsinfo%obs_nvar = 3
      allocate(obsinfo%obs_varname(obsinfo%obs_nvar))
      obsinfo%obs_varname(1) = 'u'
      obsinfo%obs_varname(2) = 'v'
      obsinfo%obs_varname(3) = 'pt'
    elseif (trim(adjustl(obstypename)) == 'ADPSFC') then
      obsinfo%obs_nz   = 1
      obsinfo%obs_nvar = 3
      allocate(obsinfo%obs_varname(obsinfo%obs_nvar))
      obsinfo%obs_varname(1) = 'u'
      obsinfo%obs_varname(2) = 'v'
      obsinfo%obs_varname(3) = 'pt'
    endif

    allocate(obsinfo%lat(obsinfo%obs_num))
    allocate(obsinfo%lon(obsinfo%obs_num))
    allocate(obsinfo%ph(obsinfo%obs_nz, obsinfo%obs_num))
    allocate(obsinfo%dis_yo(obsinfo%obs_nvar, obsinfo%obs_nz, obsinfo%obs_num))
    allocate(obsinfo%pro_yb(ens, obsinfo%obs_nvar, obsinfo%obs_nz, obsinfo%obs_num))
    allocate(obsinfo%obs_rinv(obsinfo%obs_nvar, obsinfo%obs_nz, obsinfo%obs_num))

  end subroutine

  subroutine allocate_obslist_by_num_regular(obs_num_in, obsinfo)

    integer, intent(in) :: obs_num_in
    type(obs_list), intent(inout)  :: obsinfo
    integer :: n

    obsinfo%obs_num = obs_num_in
#if (SINGLE_POINT == 0)
    obsinfo%obs_nz   = (mpas_num_lev - mpas_num_lev_start) / obs_vtc_interval + 1
#elif (SINGLE_POINT == 1)
    obsinfo%obs_nz   = 1
#endif

    obsinfo%obs_nvar = 0
    do n = 1, 4
      if (obs_var_name(n) == 'pt' .or. &
          obs_var_name(n) == 'u'  .or. &
          obs_var_name(n) == 'v') then
        obsinfo%obs_nvar = obsinfo%obs_nvar+1
      endif
    enddo
    allocate(obsinfo%obs_varname(obsinfo%obs_nvar))
    do n = 1, obsinfo%obs_nvar
      obsinfo%obs_varname(n) = obs_var_name(n)
    enddo

    allocate(obsinfo%lat(obsinfo%obs_num))
    allocate(obsinfo%lon(obsinfo%obs_num))
    allocate(obsinfo%ph(obsinfo%obs_nz, obsinfo%obs_num))
    allocate(obsinfo%dis_yo(obsinfo%obs_nvar, obsinfo%obs_nz, obsinfo%obs_num))
    allocate(obsinfo%pro_yb(ens, obsinfo%obs_nvar, obsinfo%obs_nz, obsinfo%obs_num))
    allocate(obsinfo%obs_rinv(obsinfo%obs_nvar, obsinfo%obs_nz, obsinfo%obs_num))

  end subroutine

  subroutine read_obs(rc, ntimes, lats_grid, late_grid, lons_grid, lone_grid, obsinfo)
    implicit none
    type(c_ptr), intent(in)       :: rc                      ! redis connection
    integer, intent(in)           :: ntimes                  ! time of data assimilation cycle
    real(4), intent(in)           :: lats_grid, late_grid    ! latitude of first model grid and last grid in the block
    real(4), intent(in)           :: lons_grid, lone_grid    ! longitude of first model grid and last grid in the block
    type(obs_list), intent(inout) :: obsinfo                 ! obs information
    character(:),allocatable      :: hkey, vkey
    type(datetime_type)              da_start_time
    type(datetime_type)              curr_time
    integer                       :: n, i, j, indx
    character(10)                 :: obstypename
    real(4), allocatable          :: lon_buf(:), lat_buf(:), lev_buf(:), randn_1d(:)
    real(4), allocatable          :: val(:,:), val_buf(:,:)
    character(12)                 :: date_char
    integer                       :: obs_num

    if (obsinfo%obs_num == 0) return

    obstypename = obsinfo%obs_type

    ! calculate date
    da_start_time = create_datetime(year=da_start_time_array(1),  &
                                   month=da_start_time_array(2), &
                                   day=da_start_time_array(3),   &
                                   hour=da_start_time_array(4),  &
                                   minute=da_start_time_array(5))
    curr_time = da_start_time+create_timedelta(seconds=da_in_seconds*ntimes)
    date_char   = curr_time%format('%Y%m%d%H%M')

    ! read variable value from redis
    hkey = "obs:" // trim(adjustl(obs_sceneid)) // ":" // trim(adjustl(obstypename)) //":" // trim(adjustl(date_char))

    do n = 1, obsinfo%obs_nvar
      if (n == 1) then
        vkey = trim(obsinfo%obs_varname(n))
      else
        vkey = trim(vkey) // ":" // trim(obsinfo%obs_varname(n))
      endif
        vkey = trim(vkey) // ":" // trim(obsinfo%obs_varname(n)) // "oe"
    enddo

    if (trim(adjustl(obstypename)) == 'ADPUPA') then
      obs_num = obsinfo%obs_num*adpupa_num_lev
    elseif (trim(adjustl(obstypename)) == 'GPSRO') then
      obs_num = obsinfo%obs_num*gpsro_num_lev
    else
      obs_num = obsinfo%obs_num
    endif

    allocate(lon_buf(obs_num), lat_buf(obs_num), lev_buf(obs_num), &
             val_buf(obsinfo%obs_nvar*2+1, obs_num), randn_1d(obs_num))

    call RedisGetobserveF(rc, trim(hkey), trim(vkey)//':p', lons_grid, lone_grid, lats_grid, late_grid, &
                          lon_buf, lat_buf, lev_buf, val_buf, obsinfo%obs_nvar*2+1)

    ! add perturbation: dis_yo_real + std*N(0,1)
    allocate(val(obsinfo%obs_nvar, obs_num))
    do i = 1, obsinfo%obs_nvar
      randn_1d = randn(obs_num)
      do j = 1, obs_num
        val(i, j) = val_buf(i*2-1, j) + randn_1d(j) * val_buf(i*2, j)
      enddo
    enddo

    ! assign and reshape data
    indx = 0
    if (trim(adjustl(obstypename)) == 'GPSRO' ) then
      do i = 1, obs_num
        if (lev_buf(i) == 1) then
          indx = indx + 1
          obsinfo%lat(indx) = lat_buf(i)
          obsinfo%lon(indx) = lon_buf(i)
        endif
        obsinfo%ph(lev_buf(i), indx) = val_buf(obsinfo%obs_nvar*2+1, i)
        do n = 1, obsinfo%obs_nvar
          obsinfo%dis_yo(n,lev_buf(i), indx) = val(n, i)
          obsinfo%obs_rinv(n,lev_buf(i), indx) = 1.0 / (val_buf(n*2, i)**2)
        enddo
      enddo
    elseif (trim(adjustl(obstypename)) == 'ADPUPA') then
      do i = 1, obs_num
        if (lev_buf(i) == 1) then
          indx = indx + 1
          obsinfo%lat(indx) = lat_buf(i)
          obsinfo%lon(indx) = lon_buf(i)
        endif
        obsinfo%ph(lev_buf(i), indx) = val_buf(obsinfo%obs_nvar*2+1, i)
        do n = 1, obsinfo%obs_nvar
          obsinfo%dis_yo(n,lev_buf(i), indx) = val(n, i)
          obsinfo%obs_rinv(n,lev_buf(i), indx) = 1.0 / (val_buf(n*2, i)**2)
        enddo
      enddo
    else
      do i = 1, obs_num     
        obsinfo%lat(i) = lat_buf(i)
        obsinfo%lon(i) = lon_buf(i)
        obsinfo%ph(1,i)  = val_buf(obsinfo%obs_nvar*2+1, i)
        do n = 1, obsinfo%obs_nvar
          obsinfo%dis_yo(n,1,i) = val(n, i)
          obsinfo%obs_rinv(n,1,i) = 1.0 / (val_buf(n*2, i)**2)
        enddo
      enddo
    endif

    deallocate(lon_buf)
    deallocate(lat_buf)
    deallocate(lev_buf)
    deallocate(val_buf)
    deallocate(randn_1d)
    deallocate(val)
  end subroutine read_obs

  subroutine read_obs_regular(rc, ntimes, lats, late, lons, lone, obsinfo)
    implicit none
    type(c_ptr), intent(in)       :: rc                      ! redis connection
    integer, intent(in)           :: ntimes                  ! time of data assimilation cycle
    integer, intent(in)           :: lats, late              ! latitude of first model grid and last grid in the block
    integer, intent(in)           :: lons, lone              ! longitude of first model grid and last grid in the block
    type(obs_list), intent(inout) :: obsinfo                 ! obs information
    character(:), allocatable     :: hkey, vkey, hash_key_lat, hash_key_lon
    character(5)                  :: num_lon_str, num_lat_str
    integer                       :: latids_start, latids_end, latidn_start, latidn_end
    integer                       :: lonidw_start, lonidw_end, lonide_start, lonide_end
    real(4), allocatable          :: obs_latn_tmp(:), obs_lats_tmp(:)
    real(4), allocatable          :: obs_lone_tmp(:), obs_lonw_tmp(:)
    real(4), allocatable          :: obs_r_tmp(:,:)
    real(4), allocatable          :: yo_nw_tmp_3d(:,:,:,:), yo_sw_tmp_3d(:,:,:,:), &
                                     yo_ne_tmp_3d(:,:,:,:), yo_se_tmp_3d(:,:,:,:), &
                                     yo_tmp_3d(:,:,:,:)
    type(datetime_type)              da_start_time
    type(datetime_type)              curr_time
    integer                       :: n, i, j, k, indx, var
    character(10)                 :: obstypename
    character(12)                 :: date_char
    character(3)                  :: time_str
    real(4), allocatable          :: randn_1d(:), randn_3d(:,:,:)
    integer                       :: obs_lat_num, obs_lon_num

    ! read lon and lat
    write(num_lon_str, "(i5)") num_lon
    write(num_lat_str, "(i5)") num_lat

    hash_key_lat = "lat:" // trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str))
    hash_key_lon = "lon:" // trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str))

    latids_start = obs_vertex_id(3)
    latids_end   = obs_vertex_id(4)
    latidn_start = obs_vertex_id(3)+1
    latidn_end   = obs_vertex_id(4)+1

    lonidw_start = obs_vertex_id(1)
    lonidw_end   = obs_vertex_id(2)
    lonide_start = obs_vertex_id(1)+1
    lonide_end   = obs_vertex_id(2)+1

#if (SINGLE_POINT == 0)
    obs_lat_num  = (latidn_end-latidn_start)/obs_hrz_interval_lat+1
    obs_lon_num  = (lonide_end-lonide_start)/obs_hrz_interval_lon+1
#elif (SINGLE_POINT == 1)
    obs_lat_num  = 1
    obs_lon_num  = 1
#endif
    
    allocate(obs_latn_tmp(obs_lat_num), obs_lats_tmp(obs_lat_num))
    allocate(obs_lonw_tmp(obs_lon_num), obs_lone_tmp(obs_lon_num))
    call RedisHmgetf1d(rc, hash_key_lat, "lat", latidn_start, latidn_end, obs_hrz_interval_lat, 1, obs_latn_tmp)
    call RedisHmgetf1d(rc, hash_key_lat, "lat", latids_start, latids_end, obs_hrz_interval_lat, 1, obs_lats_tmp)
    call RedisHmgetf1d(rc, hash_key_lon, "lon", lonidw_start, lonidw_end, obs_hrz_interval_lon, 1, obs_lonw_tmp)
    call RedisHmgetf1d(rc, hash_key_lon, "lon", lonide_start, lonide_end, obs_hrz_interval_lon, 1, obs_lone_tmp)

    indx = 0
    do i = 1, obs_lon_num
      do j = 1, obs_lat_num
        indx = indx + 1
        obsinfo%lat(indx) = (obs_latn_tmp(j)+obs_lats_tmp(j))/2.0
        obsinfo%lon(indx) = (obs_lonw_tmp(i)+obs_lone_tmp(i))/2.0
      enddo
    enddo

    deallocate(obs_latn_tmp, obs_lats_tmp, obs_lonw_tmp, obs_lone_tmp)

    ! read obs_r
    hkey = "obserror:" // trim(adjustl(atm_mpas_sceneid)) // ":" //&
            trim(adjustl(num_lon_str)) // "x"//trim(adjustl(num_lat_str))
    do n = 1, obsinfo%obs_nvar
      if (n == 1) then
        vkey = trim(obsinfo%obs_varname(n))
      else
        vkey = trim(vkey) // ":" // trim(obsinfo%obs_varname(n))
      endif
    enddo

    allocate(obs_r_tmp(obsinfo%obs_nvar, obsinfo%obs_nz))
#if (SINGLE_POINT == 0)
    call RedisHmgetf1d(rc, hkey, vkey, mpas_num_lev_start, &
                       mpas_num_lev, obs_vtc_interval, obsinfo%obs_nvar, obs_r_tmp)
#elif (SINGLE_POINT == 1)
    call RedisHmgetf1d(rc, hkey, vkey, onepoint_z, onepoint_z, 1, obsinfo%obs_nvar, obs_r_tmp)
#endif

    do n = 1, obsinfo%obs_nvar
      do i = 1, obsinfo%obs_nz
        do j = 1, obsinfo%obs_num
          obsinfo%obs_rinv(n, i, j) = 1.0 / (obs_r_tmp(n, i)**2)
        enddo
      enddo
    enddo
    deallocate(obs_r_tmp)

    ! read yo
    da_start_time = create_datetime(year=da_start_time_array(1),  &
                                   month=da_start_time_array(2), &
                                   day=da_start_time_array(3),   &
                                   hour=da_start_time_array(4),  &
                                   minute=da_start_time_array(5))
    curr_time = da_start_time+create_timedelta(seconds=da_in_seconds*ntimes)
    date_char = curr_time%format('%Y%m%d%H%M')

    hkey = "realfield:" // trim(adjustl(atm_mpas_sceneid)) //":" //&
           trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str)) &
           // ":" //trim(adjustl(date_char))
    vkey = 'p'
    do n = 1, obsinfo%obs_nvar
      vkey = trim(vkey) // ":" // trim(obsinfo%obs_varname(n))
    enddo

    allocate(yo_nw_tmp_3d(obsinfo%obs_nvar+1, obsinfo%obs_nz, obs_lat_num, obs_lon_num), &
             yo_sw_tmp_3d(obsinfo%obs_nvar+1, obsinfo%obs_nz, obs_lat_num, obs_lon_num), &
             yo_ne_tmp_3d(obsinfo%obs_nvar+1, obsinfo%obs_nz, obs_lat_num, obs_lon_num), &
             yo_se_tmp_3d(obsinfo%obs_nvar+1, obsinfo%obs_nz, obs_lat_num, obs_lon_num), & 
             yo_tmp_3d   (obsinfo%obs_nvar+1, obsinfo%obs_nz, obs_lat_num, obs_lon_num))
    allocate(randn_1d(obsinfo%obs_num), randn_3d(obsinfo%obs_nvar, obsinfo%obs_nz, obsinfo%obs_num))
#if (SINGLE_POINT == 0)
    call RedisHmgetf3d(rc, hkey, vkey, obs_vertex_id(1), obs_vertex_id(2), &
         obs_hrz_interval_lon, obs_vertex_id(3)+1, obs_vertex_id(4)+1, obs_hrz_interval_lat, &
         mpas_num_lev_start, mpas_num_lev, obs_vtc_interval, obsinfo%obs_nvar+1, yo_nw_tmp_3d)

    call RedisHmgetf3d(rc, hkey, vkey, obs_vertex_id(1)+1, obs_vertex_id(2)+1, &
         obs_hrz_interval_lon, obs_vertex_id(3)+1,obs_vertex_id(4)+1, obs_hrz_interval_lat, &
         mpas_num_lev_start, mpas_num_lev, obs_vtc_interval, obsinfo%obs_nvar+1, yo_ne_tmp_3d)

    call RedisHmgetf3d(rc, hkey, vkey, obs_vertex_id(1), obs_vertex_id(2), &
         obs_hrz_interval_lon, obs_vertex_id(3), obs_vertex_id(4), obs_hrz_interval_lat, &
         mpas_num_lev_start, mpas_num_lev, obs_vtc_interval, obsinfo%obs_nvar+1, yo_sw_tmp_3d)

    call RedisHmgetf3d(rc, hkey, vkey, obs_vertex_id(1)+1, obs_vertex_id(2)+1,&
         obs_hrz_interval_lat, obs_vertex_id(3), obs_vertex_id(4),obs_hrz_interval_lon, &
         mpas_num_lev_start, mpas_num_lev, obs_vtc_interval, obsinfo%obs_nvar+1, yo_se_tmp_3d)
#elif (SINGLE_POINT == 1)
    call RedisHmgetf3d(rc, hkey, vkey, obs_vertex_id(1), obs_vertex_id(2), &
         obs_hrz_interval_lon, obs_vertex_id(3)+1, obs_vertex_id(4)+1, obs_hrz_interval_lat, &
         onepoint_z, onepoint_z, 1, obsinfo%obs_nvar+1, yo_nw_tmp_3d)

    call RedisHmgetf3d(rc, hkey, vkey, obs_vertex_id(1)+1, obs_vertex_id(2)+1, &
         obs_hrz_interval_lon, obs_vertex_id(3)+1,obs_vertex_id(4)+1, obs_hrz_interval_lat, &
         onepoint_z, onepoint_z, 1, obsinfo%obs_nvar+1, yo_ne_tmp_3d)

    call RedisHmgetf3d(rc, hkey, vkey, obs_vertex_id(1), obs_vertex_id(2), &
         obs_hrz_interval_lon, obs_vertex_id(3), obs_vertex_id(4), obs_hrz_interval_lat, &
         onepoint_z, onepoint_z, 1, obsinfo%obs_nvar+1, yo_sw_tmp_3d)

    call RedisHmgetf3d(rc, hkey, vkey, obs_vertex_id(1)+1, obs_vertex_id(2)+1,&
         obs_hrz_interval_lat, obs_vertex_id(3), obs_vertex_id(4),obs_hrz_interval_lon, &
         onepoint_z, onepoint_z, 1, obsinfo%obs_nvar+1, yo_se_tmp_3d)
#endif
    yo_tmp_3d =(yo_nw_tmp_3d + yo_sw_tmp_3d + yo_ne_tmp_3d + yo_se_tmp_3d)/4.0

    indx = 0
    do i = 1, obs_lon_num
      do j = 1, obs_lat_num
        indx = indx + 1
        do k = 1, obsinfo%obs_nz
          obsinfo%ph(k, indx)  = yo_tmp_3d(1,k,j,i)
        enddo
      enddo
    enddo

    write(time_str,"(i3.3)")ntimes+1 !read gassian perturb
    if (randn_filepath == 'random') then
      do n = 1, obsinfo%obs_nvar
        do k = 1, obsinfo%obs_nz
          randn_1d = randn(obsinfo%obs_num)
          indx = 0
          do i = 1, obs_lon_num
            do j = 1, obs_lat_num
              indx = indx + 1
              obsinfo%dis_yo(n, k, indx) = yo_tmp_3d(n+1,k,j,i) + randn_1d(indx)/sqrt(obsinfo%obs_rinv(n,k,indx))
            end do
          end do
        enddo
      enddo
    else
      open(1, file=trim(randn_filepath) // '/3d_'//trim(time_str) //'.txt')
      read(1, 130)(((randn_3d(var,k,i),var=1,obsinfo%obs_nvar),k=1,obsinfo%obs_nz),i=1,obsinfo%obs_num)
      close(1)
      do n = 1, obsinfo%obs_nvar
        do k = 1, obsinfo%obs_nz
          indx = 0
          do i = 1, obs_lon_num
            do j = 1, obs_lat_num
              indx = indx + 1
              obsinfo%dis_yo(n, k, indx) = yo_tmp_3d(n+1,k,j,i) + randn_3d(n,k,indx)/sqrt(obsinfo%obs_rinv(n,k,indx))
            end do
          end do
        enddo
      enddo
    endif

130     format(20f16.8)
    deallocate(yo_nw_tmp_3d, yo_sw_tmp_3d, yo_ne_tmp_3d, yo_se_tmp_3d, yo_tmp_3d)
    deallocate(randn_1d, randn_3d)
!print *, "lat, lon, rinv, ph, yo ", obsinfo%lat, obsinfo%lon, obsinfo%obs_rinv, obsinfo%ph, obsinfo%dis_yo
  end subroutine read_obs_regular

  subroutine obs_final(obsinfo)
    type(obs_list), intent(inout)   :: obsinfo

    if (allocated(obsinfo%obs_varname)) deallocate(obsinfo%obs_varname)
    if (allocated(obsinfo%obs_npos)) deallocate(obsinfo%obs_npos)
    if (allocated(obsinfo%lat)) deallocate(obsinfo%lat)
    if (allocated(obsinfo%lon)) deallocate(obsinfo%lon)
    if (allocated(obsinfo%ph)) deallocate(obsinfo%ph)
    if (allocated(obsinfo%dis_yo)) deallocate(obsinfo%dis_yo)
    if (allocated(obsinfo%pro_yb)) deallocate(obsinfo%pro_yb)
    if (allocated(obsinfo%obs_rinv)) deallocate(obsinfo%obs_rinv)

  end subroutine obs_final
end module read_obs_mod

