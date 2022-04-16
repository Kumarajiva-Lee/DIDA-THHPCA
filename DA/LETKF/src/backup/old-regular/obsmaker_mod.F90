#include "../../utils/da_config.inc"
#define SINGLE_POINT 0
module obsmaker_mod
! observation make

  use coupler_config, only:num_lon, num_lat, atm_mpas_sceneid, da_in_seconds, case_name, &
                      da_start_time_array=>da_start_time
  use da_namelist_mod, only:obs_var_name, mpas_num_lev_start, mpas_num_lev, obs_hrz_interval_lon, &
                      obs_hrz_interval_lat, obs_vtc_interval, onepoint_x, onepoint_y, onepoint_z
  use redis_mod
  use flogger
  use string
  
  implicit none

  integer    :: obs_num_2d !num of 2d obs var
  integer    :: obs_num_3d !num of 3d obs var
  integer    :: obs_l !num of obs lev
  character(len = 20),dimension(:), allocatable :: obs_var_name_list
!  real(4), allocatable, dimension(:)   :: obs_r !obs error of 2d var

  public obs_num_2d, obs_num_3d, obs_l, obs_var_name_list

  public obs_init
  public obs_get_r
  public obs_get_loc
  public obs_get_yo

contains

! prepare obs_varname_name_list, obs_num_2d, obs_num_3d, obs_l
subroutine obs_init(n_var, n_start, n_end)

    use var_info_mod
    use splitchar
    use coupler_config, only:num_lev, da_var_name

    implicit none

    integer, intent(out) :: n_var
    integer, intent(out) :: n_start(10), n_end(10)
    integer :: obs_var_num, da_var_num !num of total obs var
    character(20), dimension(10):: obs_var_name_list_origin
    character(20), dimension(10):: da_var_name_list_origin
    character(len = 20),dimension(:), allocatable :: da_var_name_list
    integer :: da_num_2d, da_num_3d
    integer :: i, j
! get obs var name

    call stringsplit(obs_var_name, ',', obs_var_name_list_origin, obs_var_num)

    allocate(obs_var_name_list(obs_var_num))
    print*, da_var_name, obs_var_name, obs_var_name_list_origin
    call var_info(obs_var_name_list_origin(1:obs_var_num), obs_var_num, obs_var_name_list, obs_num_2d, obs_num_3d)
    ! get da var name
    print *, obs_num_2d, obs_num_3d
    call stringsplit(da_var_name, ',', da_var_name_list_origin, da_var_num)

    allocate(da_var_name_list(da_var_num))

    call var_info(da_var_name_list_origin(1:da_var_num), da_var_num, da_var_name_list, da_num_2d, da_num_3d)

    print*, 'varnamelist ',da_var_name_list, obs_var_name_list
    n_var = 0
    do i = 1, da_num_2d
        do j = 1, obs_num_2d
            if (trim(da_var_name_list(i)) .eq. trim(obs_var_name_list(j))) then
                n_start(n_var+1) = 1+n_var
                n_end(n_var+1) = 1+n_var
                n_var = n_var+1
            endif
        end do
    end do

    do j = obs_num_2d+1, obs_num_2d+obs_num_3d
        do i = da_num_2d+1, da_num_2d+da_num_3d
            if (trim(da_var_name_list(i)) .eq. trim(obs_var_name_list(j))) then
                n_start(n_var+1) = da_num_2d+1+(i-1-da_num_2d)*num_lev
                n_end(n_var+1) = da_num_2d+(i-da_num_2d)*num_lev
                n_var = n_var+1
            endif
        end do
    end do

!    obs_l = CEILING(real(mpas_num_lev - mpas_num_lev_start)/real(obs_vtc_interval))
    obs_l = floor(real(mpas_num_lev - mpas_num_lev_start)/real(obs_vtc_interval))+1
#if (SINGLE_POINT == 1)
    obs_l = 1
#endif
!    print *,'var list:', da_var_name_list, obs_var_name_list,'n_start', n_start, 'n_end',n_end,'n_var',n_var
end subroutine obs_init

!calculate obs number in a block
subroutine obs_get_num(lats, late, lons, lone, obs_lat_num, obs_lon_num, obs_vertex_id)
    integer, intent(in)     :: lats, late, lons, lone !lonid and latid of block
    integer, intent(out)    :: obs_lat_num, obs_lon_num !obs num in y/x direction in block
    integer, dimension(4), intent(out) :: obs_vertex_id !1:start id of lon; 2:end id of lon 3: start id of lat; 4:end id of lat
    integer :: latids, latide, lonids, lonide !lonid and latid of block

    latids = lats
    lonids = lons
    latide = late
    lonide = lone
    if (late .lt. num_lat) latide = late+1

    if (lone .lt. num_lon) lonide = lone+1
    !num obs in x direction(on the right of gmcore grid)
#if (SINGLE_POINT == 0)
    obs_vertex_id(1) = obs_hrz_interval_lon-mod(lonids-1,obs_hrz_interval_lon)+lonids
    obs_vertex_id(2) = -mod(lonide-1,obs_hrz_interval_lon)+lonide
    if (mod(lonids-1,obs_hrz_interval_lon)<1) then
        obs_vertex_id(1) = obs_vertex_id(1) - obs_hrz_interval_lon
    endif

    if (mod(lonide-1,obs_hrz_interval_lon)<1) then
        obs_vertex_id(2) = obs_vertex_id(2) - obs_hrz_interval_lon
    endif
    obs_lon_num = (obs_vertex_id(2) - obs_vertex_id(1))/obs_hrz_interval_lon+1
#else
    obs_vertex_id(1) = onepoint_x
    obs_vertex_id(2) = onepoint_x
    obs_lon_num = 1
#endif

    !num obs in y direction(symmetrical)
#if (SINGLE_POINT == 0)
    obs_vertex_id(3) = obs_hrz_interval_lat-mod(latids-mod(num_lat/2, obs_hrz_interval_lat),obs_hrz_interval_lat)+latids
    obs_vertex_id(4) = -mod(latide-mod(num_lat/2, obs_hrz_interval_lat),obs_hrz_interval_lat)+latide
    if (mod(latids-mod(num_lat/2, obs_hrz_interval_lat),obs_hrz_interval_lat)<1) then
        obs_vertex_id(3) = obs_vertex_id(3) - obs_hrz_interval_lat
    endif
    if (mod(latide-mod(num_lat/2, obs_hrz_interval_lat),obs_hrz_interval_lat)<1) then
        obs_vertex_id(4) = obs_vertex_id(4) - obs_hrz_interval_lat
    endif
    obs_lat_num = (obs_vertex_id(4) - obs_vertex_id(3))/obs_hrz_interval_lat+1
#else
    obs_vertex_id(3) = onepoint_y
    obs_vertex_id(4) = onepoint_y
    obs_lat_num = 1
#endif
end subroutine obs_get_num

! read obs lon and lat from redis
subroutine obs_get_loc(rc, obs_lon_num, obs_lat_num, obs_vertex_id, obs_lon, obs_lat)

    type(c_ptr), intent(in) :: rc
    integer, dimension(4), intent(in) :: obs_vertex_id
    integer, intent(in)               :: obs_lat_num, obs_lon_num

    real(4), intent(out) :: obs_lat(obs_lat_num)
    real(4), intent(out) :: obs_lon(obs_lon_num)

    integer :: i,j
    character(:),allocatable :: hash_key_lat, hash_key_lon
    character(5) :: num_lon_str, num_lat_str
    integer :: latids_start, latids_end, latidn_start, latidn_end, lonidw_start, lonidw_end, lonide_start, lonide_end
    real(4) :: obs_lone_tmp(obs_lon_num), obs_lonw_tmp(obs_lon_num), obs_latn_tmp(obs_lat_num), obs_lats_tmp(obs_lat_num)
    integer :: obs_latids0, obs_latids, obs_latide, obs_lonids, obs_lonide

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

    call RedisHmgetf1d(rc, hash_key_lat, "lat", latidn_start, latidn_end, obs_hrz_interval_lat, 1, obs_latn_tmp)
    call RedisHmgetf1d(rc, hash_key_lat, "lat", latids_start, latids_end, obs_hrz_interval_lat, 1, obs_lats_tmp)
    call RedisHmgetf1d(rc, hash_key_lon, "lon", lonidw_start, lonidw_end, obs_hrz_interval_lon, 1, obs_lonw_tmp)
    call RedisHmgetf1d(rc, hash_key_lon, "lon", lonide_start, lonide_end, obs_hrz_interval_lon, 1, obs_lone_tmp)

    obs_lat = (obs_latn_tmp+obs_lats_tmp)/2.0
    obs_lon = (obs_lonw_tmp+obs_lone_tmp)/2.0
!---------obs lat lon output------------!
    obs_latids0 = obs_hrz_interval_lat-mod(1-mod(num_lat/2, obs_hrz_interval_lat),obs_hrz_interval_lat)+1
    if (mod(1-mod(num_lat/2, obs_hrz_interval_lat),obs_hrz_interval_lat)<1) then
        obs_latids0 = obs_latids0 - obs_hrz_interval_lat
    endif

    obs_latids = (obs_vertex_id(3)-obs_latids0)/obs_hrz_interval_lat+1
    obs_latide = (obs_vertex_id(4)-obs_latids0)/obs_hrz_interval_lat+1
    obs_lonids = (obs_vertex_id(1)-1)/obs_hrz_interval_lon+1
    obs_lonide = (obs_vertex_id(2)-1)/obs_hrz_interval_lon+1

    hash_key_lat = "obslat:" // trim(adjustl(case_name))// &
                   ":" //trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str))
    hash_key_lon = "obslon:" // trim(adjustl(case_name))// &
                   ":" // trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str))

    call RedisHmsetf1d(rc, hash_key_lat, "lat", obs_latids, obs_latide, 1, 1, obs_lat)
    call RedisHmsetf1d(rc, hash_key_lon, "lon", obs_lonids, obs_lonide, 1, 1, obs_lon)
!---------------------------------------!
end subroutine obs_get_loc

subroutine obs_get_r(rc, obs_r)

    type(c_ptr), intent(in) :: rc
    real(4), intent(out)     :: obs_r(obs_num_2d+obs_num_3d*obs_l)

    character(30)            :: hash_key, var_key
    character(5)            :: num_lon_str, num_lat_str
    integer                 :: i
    real(4)                 :: obs_r_tmp(obs_num_3d, obs_l)

    !print *, 'Read obs err cov'
    write(num_lon_str, "(i5)") num_lon
    write(num_lat_str, "(i5)") num_lat

    hash_key = "obserror:" // trim(adjustl(atm_mpas_sceneid)) // ":" //&
            trim(adjustl(num_lon_str)) // "x"//trim(adjustl(num_lat_str))

!print *, obs_num_2d, obs_num_3d
    if (obs_num_2d .ne. 0) then
        do i = 1, obs_num_2d
            var_key = "0:"//trim(adjustl(obs_var_name_list(i)))
            call RedisHget(rc, hash_key, var_key, obs_r(i))
!            obs_r(i)=1
            obs_r(i)=obs_r(i)**2
        end do
    endif

    if (obs_num_3d .ne. 0) then
        var_key = trim(adjustl(obs_var_name_list(1+obs_num_2d)))
        !print *, obs_var_name_list
        do i=2+obs_num_2d, obs_num_2d+obs_num_3d
            var_key = trim(adjustl(var_key)) // ':' // trim(adjustl(obs_var_name_list(i)))
        end do
#if (SINGLE_POINT == 0)
!        call RedisHmgetf1d(rc, hash_key, var_key, mod(mpas_num_lev-1, obs_vtc_interval)+1, &
!             mpas_num_lev, obs_vtc_interval, obs_num_3d, obs_r_tmp)
        call RedisHmgetf1d(rc, hash_key, var_key, mpas_num_lev_start, &
             mpas_num_lev, obs_vtc_interval, obs_num_3d, obs_r_tmp)
#else
        call RedisHmgetf1d(rc, hash_key, var_key, onepoint_z, onepoint_z, 1, obs_num_3d, obs_r_tmp)
#endif
        do i=1,obs_num_3d
            obs_r(obs_num_2d+1+(i-1)*obs_l:obs_num_2d+i*obs_l) = obs_r_tmp(i,:)**2
           ! obs_r(obs_num_2d+1+(i-1)*obs_l:obs_num_2d+i*obs_l) = 1
        end do
    endif
end subroutine obs_get_r

subroutine obs_get_yo(rc, obs_vertex_id, obs_lat_num, obs_lon_num, ntimes, obs_rinv, obs_p, yo)

    use da_namelist_mod, only:randn_filepath
    use datetime
    use randn_mod

    type(c_ptr), intent(in) :: rc
    integer, dimension(4), intent(in) :: obs_vertex_id
    integer, intent(in)               :: obs_lat_num, obs_lon_num
    integer, intent(in) :: ntimes
    real(4), intent(in) :: obs_rinv(obs_num_2d+obs_num_3d*obs_l)
    !real(4), intent(out) :: obs_p(obs_l, obs_lat_num*obs_lon_num)
    !real(4), intent(out) :: yo(obs_num_2d+obs_num_3d*obs_l, obs_lat_num*obs_lon_num)    
    ayb_type, intent(out) :: obs_p(obs_l, obs_lat_num*obs_lon_num)
    ayb_type, intent(out) :: yo(obs_num_2d+obs_num_3d*obs_l, obs_lat_num*obs_lon_num)

    character(len = 12) :: date_char
    character(:),allocatable :: hash_key, var_key
    integer :: i,k,var,ix,iy,iobs
    integer :: total_obs
    character(5) :: num_lon_str, num_lat_str
    real(4) :: yo_nw_tmp_2d(obs_num_2d, obs_lat_num, obs_lon_num),&
               yo_sw_tmp_2d(obs_num_2d, obs_lat_num, obs_lon_num),&
               yo_ne_tmp_2d(obs_num_2d, obs_lat_num, obs_lon_num),&
               yo_se_tmp_2d(obs_num_2d, obs_lat_num, obs_lon_num)
    real(4) :: yo_nw_tmp_2d_sl(obs_num_2d, obs_lat_num, obs_lon_num)
    real(4) :: yo_nw_tmp_3d((obs_num_3d+1), obs_l, obs_lat_num, obs_lon_num),&
               yo_sw_tmp_3d((obs_num_3d+1), obs_l, obs_lat_num, obs_lon_num),&
               yo_ne_tmp_3d((obs_num_3d+1), obs_l, obs_lat_num, obs_lon_num),&
               yo_se_tmp_3d((obs_num_3d+1), obs_l, obs_lat_num, obs_lon_num)
    real(4) :: yo_tmp_2d(obs_num_2d, obs_lat_num, obs_lon_num)
    real(4) :: yo_tmp_3d(obs_num_3d+1, obs_l, obs_lat_num, obs_lon_num)
    real(4) :: randn_3d(obs_num_3d, obs_l, obs_lat_num*obs_lon_num)
    real(4) :: randn_2d(obs_num_2d, obs_lat_num*obs_lon_num)
    real(4) :: randn_1d(obs_lat_num*obs_lon_num)
    real(4) ::  tmp
    character(3) :: time_str
    type(datetime_type) da_start_time
    type(datetime_type) curr_time
    type(timedelta_type) dt

    integer :: obs_latids0, obs_latids, obs_latide, obs_lonids, obs_lonide
    real(4) :: yo_out_2d(obs_num_2d, obs_lat_num, obs_lon_num)
    real(4) :: yo_out_3d(obs_num_3d, obs_l, obs_lat_num, obs_lon_num)
    !print *, 'Read OBS values'
    total_obs = obs_lon_num*obs_lat_num

    write(num_lon_str, "(i5)") num_lon
    write(num_lat_str, "(i5)") num_lat

    da_start_time = create_datetime(year=da_start_time_array(1),  &
                                   month=da_start_time_array(2), &
                                   day=da_start_time_array(3),   &
                                   hour=da_start_time_array(4),  &
                                   minute=da_start_time_array(5))
    curr_time = da_start_time+create_timedelta(seconds=da_in_seconds*ntimes)
    date_char = curr_time%format('%Y%m%d%H%M')

    hash_key = "realfield:" // trim(adjustl(atm_mpas_sceneid)) //":" //&
               trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str)) &
               // ":" //trim(adjustl(date_char))

    if (obs_num_2d .ne. 0) then
        var_key = trim(adjustl(obs_var_name_list(1)))
        do i=2, obs_num_2d
            var_key = trim(adjustl(var_key)) //':'//trim(adjustl(obs_var_name_list(i)))
        end do

    call log_print("redis_hmget2d start")
    call log_print("da:"//to_str(obs_lon_num*obs_lat_num))

    call RedisHmgetf2d(rc, hash_key, var_key, obs_vertex_id(1), obs_vertex_id(2), &
         obs_hrz_interval_lat, obs_vertex_id(3)+1, obs_vertex_id(4)+1, obs_hrz_interval_lon, &
         obs_num_2d, yo_nw_tmp_2d(1:obs_num_2d, :, :))

    call RedisHmgetf2d(rc, hash_key, var_key, obs_vertex_id(1)+1,obs_vertex_id(2)+1, &
         obs_hrz_interval_lon, obs_vertex_id(3)+1, obs_vertex_id(4)+1, obs_hrz_interval_lat, &
         obs_num_2d, yo_ne_tmp_2d(1:obs_num_2d, :, :))

    call RedisHmgetf2d(rc, hash_key, var_key, obs_vertex_id(1), obs_vertex_id(2), &
         obs_hrz_interval_lon, obs_vertex_id(3), obs_vertex_id(4), obs_hrz_interval_lat, &
         obs_num_2d, yo_sw_tmp_2d(1:obs_num_2d, :, :))

    call RedisHmgetf2d(rc, hash_key, var_key, obs_vertex_id(1)+1, obs_vertex_id(2)+1,&
         obs_hrz_interval_lon, obs_vertex_id(3), obs_vertex_id(4),obs_hrz_interval_lat, &
         obs_num_2d, yo_se_tmp_2d(1:obs_num_2d, :, :))

    endif

    var_key = 'p'
    if (obs_num_3d .ne. 0) then
        do i=obs_num_2d+1, obs_num_2d+obs_num_3d
            var_key = trim(adjustl(var_key)) //':'//trim(adjustl(obs_var_name_list(i)))
        end do
    endif

#if (SINGLE_POINT == 0)
    call RedisHmgetf3d(rc, hash_key, var_key, obs_vertex_id(1), obs_vertex_id(2), &
         obs_hrz_interval_lon, obs_vertex_id(3)+1, obs_vertex_id(4)+1, obs_hrz_interval_lat, &
!         mod(mpas_num_lev-1, obs_vtc_interval)+1, mpas_num_lev, obs_vtc_interval, obs_num_3d+1, yo_nw_tmp_3d)
         mpas_num_lev_start, mpas_num_lev, obs_vtc_interval, obs_num_3d+1, yo_nw_tmp_3d)

    call RedisHmgetf3d(rc, hash_key, var_key, obs_vertex_id(1)+1, obs_vertex_id(2)+1, &
         obs_hrz_interval_lon, obs_vertex_id(3)+1,obs_vertex_id(4)+1, obs_hrz_interval_lat, &
!         mod(mpas_num_lev-1, obs_vtc_interval)+1, mpas_num_lev, obs_vtc_interval, obs_num_3d+1, yo_ne_tmp_3d) 
         mpas_num_lev_start, mpas_num_lev, obs_vtc_interval, obs_num_3d+1, yo_ne_tmp_3d)

    call RedisHmgetf3d(rc, hash_key, var_key, obs_vertex_id(1), obs_vertex_id(2), &
         obs_hrz_interval_lon, obs_vertex_id(3), obs_vertex_id(4), obs_hrz_interval_lat, &
!         mod(mpas_num_lev-1, obs_vtc_interval)+1, mpas_num_lev, obs_vtc_interval, obs_num_3d+1, yo_sw_tmp_3d)
         mpas_num_lev_start, mpas_num_lev, obs_vtc_interval, obs_num_3d+1, yo_sw_tmp_3d)

    call RedisHmgetf3d(rc, hash_key, var_key, obs_vertex_id(1)+1, obs_vertex_id(2)+1,&
         obs_hrz_interval_lat, obs_vertex_id(3), obs_vertex_id(4),obs_hrz_interval_lon, &
!         mod(mpas_num_lev-1, obs_vtc_interval)+1, mpas_num_lev, obs_vtc_interval, obs_num_3d+1, yo_se_tmp_3d)
         mpas_num_lev_start, mpas_num_lev, obs_vtc_interval, obs_num_3d+1, yo_se_tmp_3d)
#else
    call RedisHmgetf3d(rc, hash_key, var_key, obs_vertex_id(1), obs_vertex_id(2), &
         obs_hrz_interval_lon, obs_vertex_id(3)+1, obs_vertex_id(4)+1, obs_hrz_interval_lat, &
         onepoint_z, onepoint_z, 1, obs_num_3d+1, yo_nw_tmp_3d)

    call RedisHmgetf3d(rc, hash_key, var_key, obs_vertex_id(1)+1, obs_vertex_id(2)+1, &
         obs_hrz_interval_lon, obs_vertex_id(3)+1,obs_vertex_id(4)+1, obs_hrz_interval_lat, &
         onepoint_z, onepoint_z, 1, obs_num_3d+1, yo_ne_tmp_3d)

    call RedisHmgetf3d(rc, hash_key, var_key, obs_vertex_id(1), obs_vertex_id(2), &
         obs_hrz_interval_lon, obs_vertex_id(3), obs_vertex_id(4), obs_hrz_interval_lat, &
         onepoint_z, onepoint_z, 1, obs_num_3d+1, yo_sw_tmp_3d)

    call RedisHmgetf3d(rc, hash_key, var_key, obs_vertex_id(1)+1, obs_vertex_id(2)+1,&
         obs_hrz_interval_lat, obs_vertex_id(3), obs_vertex_id(4),obs_hrz_interval_lon, &
         onepoint_z, onepoint_z, 1, obs_num_3d+1, yo_se_tmp_3d)
#endif

    yo_tmp_2d = (yo_nw_tmp_2d + yo_sw_tmp_2d + yo_ne_tmp_2d + yo_se_tmp_2d)/4.0
    yo_tmp_3d =(yo_nw_tmp_3d + yo_sw_tmp_3d + yo_ne_tmp_3d + yo_se_tmp_3d)/4.0

    obs_p  = reshape(yo_tmp_3d(1,:,:,:), (/obs_l, total_obs/))
!gmcore
    write(time_str,"(i3.3)")ntimes+1 !read gassian perturb

    if (obs_num_2d .ne. 0) then

        if (randn_filepath == 'random') then
            do i = 1, obs_num_2d
                randn_1d = randn(obs_lon_num*obs_lat_num)
                iobs = 0
                do iy = 1, obs_lon_num
                    do ix = 1, obs_lat_num
                        iobs = iobs + 1
    !                    yo(i,iobs) = yo_tmp_2d(i, ix, iy)
                        yo(i,iobs) = yo_tmp_2d(i, ix, iy) + randn_1d(iobs)/sqrt(obs_rinv(i))
                    end do
                end do
            enddo
        else
            open(2, file=trim(randn_filepath) // '/2d_'//trim(time_str) //'.txt')
            read(2, 130)((randn_2d(var, i),var=1,obs_num_2d),i=1,obs_lon_num*obs_lat_num)
            close(2)

            do i = 1, obs_num_2d
                iobs = 0
                do iy = 1, obs_lon_num
                    do ix = 1, obs_lat_num
                        iobs = iobs + 1
    !                    yo(i,iobs) = yo_tmp_2d(i, ix, iy)
                        yo(i,iobs) = yo_tmp_2d(i, ix, iy) + randn_2d(i,iobs)/sqrt(obs_rinv(i))
                    end do
                end do
            enddo
        endif
    endif

    if (obs_num_3d .ne. 0) then
 
        if (randn_filepath == 'random') then
            do i = 1, obs_num_3d
                do k = 1, obs_l
                    randn_1d = randn(obs_lon_num*obs_lat_num)
                    iobs = 0
                    do iy = 1, obs_lon_num
                        do ix = 1, obs_lat_num
                            iobs = iobs + 1
    !                        yo(obs_num_2d+(i-1)*obs_l+k,iobs) = yo_tmp_3d(i+1,k,ix,iy)
                            yo(obs_num_2d+(i-1)*obs_l+k,iobs) = yo_tmp_3d(i+1,k,ix,iy) &
                            + randn_1d(iobs)/sqrt(obs_rinv(obs_num_2d+(i-1)*obs_l+k))
                        end do
                    end do
                enddo
            enddo
        else
            open(1, file=trim(randn_filepath) // '/3d_'//trim(time_str) //'.txt')
            read(1, 130)(((randn_3d(var,k,i),var=1,obs_num_3d),k=1,obs_l),i=1,obs_lon_num*obs_lat_num)
            close(1)
            do i = 1, obs_num_3d
                do k = 1, obs_l
                    iobs = 0
                    do iy = 1, obs_lon_num
                        do ix = 1, obs_lat_num
                            iobs = iobs + 1
    !                        yo(obs_num_2d+(i-1)*obs_l+k,iobs) = yo_tmp_3d(i+1,k,ix,iy)
                            yo(obs_num_2d+(i-1)*obs_l+k,iobs) = yo_tmp_3d(i+1,k,ix,iy) &
                            + randn_3d(i,k,iobs)/sqrt(obs_rinv(obs_num_2d+(i-1)*obs_l+k))
                        end do
                    end do
                enddo
            enddo
        endif
    endif
130     format(20f16.8)
!print *, 'OBS read finished!'

!----------obs_output-----------!
!    obs_latids0 = obs_hrz_interval_lat-mod(1-mod(num_lat/2, obs_hrz_interval_lat),obs_hrz_interval_lat)+1
!    if (mod(1-mod(num_lat/2, obs_hrz_interval_lat),obs_hrz_interval_lat)<1) then
!        obs_latids0 = obs_latids0 - obs_hrz_interval_lat
!    endif
!
!    obs_latids = (obs_vertex_id(3)-obs_latids0)/obs_hrz_interval_lat+1
!    obs_latide = (obs_vertex_id(4)-obs_latids0)/obs_hrz_interval_lat+1
!    obs_lonids = (obs_vertex_id(1)-1)/obs_hrz_interval_lon+1
!    obs_lonide = (obs_vertex_id(2)-1)/obs_hrz_interval_lon+1
!
!    hash_key = "obs:" // trim(adjustl(case_name)) // ":"//trim(adjustl(num_lon_str)) &
!               // "x"// trim(adjustl(num_lat_str)) // ":" //trim(adjustl(date_char))
!
!    if (obs_num_2d .ne. 0) then
!
!        var_key = trim(adjustl(obs_var_name_list(1)))
!        do i=2, obs_num_2d
!            var_key = trim(adjustl(var_key)) //':'//trim(adjustl(obs_var_name_list(i)))
!        end do
!     
!        do i=1, obs_num_2d
!            iobs = 0
!            do iy = 1, obs_lon_num
!                do ix = 1, obs_lat_num
!                    iobs = iobs + 1
!                    yo_out_2d(i,ix,iy) = yo(i,iobs) 
!                end do
!            end do
!        end do
!
!        call RedisHmsetf2d(rc, hash_key, var_key, obs_lonids, obs_lonide, 1, &
!             obs_latids, obs_latide, 1, obs_num_2d, yo_out_2d)
!    endif
!    if (obs_num_3d .ne. 0) then
!
!        var_key = trim(adjustl(obs_var_name_list(obs_num_2d+1)))
!        do i=2, obs_num_3d
!            var_key = trim(adjustl(var_key)) //':'//trim(adjustl(obs_var_name_list(obs_num_2d+i)))
!        end do
!
!        do i=1, obs_num_3d
!            do k = 1, obs_l
!                iobs = 0
!                do iy = 1, obs_lon_num
!                    do ix = 1, obs_lat_num
!                        iobs = iobs + 1
!                        yo_out_3d(i,k,ix,iy) = yo(obs_num_2d+(i-1)*obs_l+k,iobs) 
!                    end do
!                end do
!            enddo
!        end do
!
!        call RedisHmsetf3d(rc, hash_key, var_key, obs_lonids, obs_lonide, 1, &
!             obs_latids, obs_latide, 1, 1, obs_l, 1, obs_num_3d, yo_out_3d)
!    end if
!-------------------------------!
end subroutine obs_get_yo

subroutine obs_get_yo_interp(pid,rc, obs_vertex_id, obs_lat_t, obs_lon_t, obs_lat_num, obs_lon_num, ntimes, obs_rinv, obs_p, yo)
    use da_namelist_mod, only:obs_stat_interval_lon, obs_stat_interval_lat, randn_filepath
    use datetime
    use randn_mod

    integer, intent(in) :: pid 
    type(c_ptr), intent(in) :: rc
    integer, dimension(4), intent(in) :: obs_vertex_id
    integer, intent(in)               :: obs_lat_num, obs_lon_num
    real(4), intent(in) :: obs_lat_t(obs_lat_num), obs_lon_t(obs_lon_num)
    integer, intent(in) :: ntimes
    real(4), intent(in) :: obs_rinv(obs_num_2d+obs_num_3d*obs_l)
    ayb_type, intent(out) :: obs_p(obs_l, obs_lat_num*obs_lon_num)
    ayb_type, intent(out) :: yo(obs_num_2d+obs_num_3d*obs_l, obs_lat_num*obs_lon_num)

    character(len = 12) :: date_char
    character(:),allocatable :: hash_key, var_key, hash_key_lat, hash_key_lon
    character(5) :: num_lon_str, num_lat_str
    integer :: lonids, lonide, latids, latide, ix, iy, var, k, iobs, idx, idy
    real(4) :: r1, r2, t1, t2
    real(4) :: randn_1d(obs_lat_num*obs_lon_num)
    real(4),allocatable :: real_lat(:), real_lon(:)
    real(4),allocatable :: yo_tmp_2d(:, :, :), yo_tmp_3d(:,:,:,:)
    real(4) :: yo_tmp(obs_num_2d+(obs_num_3d+1)*obs_l, obs_lat_num, obs_lon_num)

    character(3) :: time_str
    type(datetime_type) da_start_time
    type(datetime_type) curr_time
    type(timedelta_type) dt

    integer :: obs_latids0, obs_latids, obs_latide, obs_lonids, obs_lonide
    real(4) :: yo_out_2d(obs_num_2d, obs_lat_num, obs_lon_num)
    real(4) :: yo_out_3d(obs_num_3d, obs_l, obs_lat_num, obs_lon_num)

    write(num_lon_str, "(i5)") num_lon
    write(num_lat_str, "(i5)") num_lat

    da_start_time = create_datetime(year=da_start_time_array(1),  &
                                   month=da_start_time_array(2), &
                                   day=da_start_time_array(3),   &
                                   hour=da_start_time_array(4),  &
                                   minute=da_start_time_array(5))
    curr_time = da_start_time+create_timedelta(seconds=da_in_seconds*ntimes)
    date_char = curr_time%format('%Y%m%d%H%M')

    hash_key = "realfield:" // trim(adjustl(atm_mpas_sceneid)) //":" //&
               trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str)) &
               // ":" //trim(adjustl(date_char))
    hash_key_lat = "lat:" // trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str))
    hash_key_lon = "lon:" // trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str))

    ! realfield id
    lonids = int(floor((float(obs_vertex_id(1))-1.0)/float(obs_stat_interval_lon))) + 1 
    lonide = int(floor((float(obs_vertex_id(2))-1.0)/float(obs_stat_interval_lon))) + 2 
    latids = int(floor((float(obs_vertex_id(3))-1.0)/float(obs_stat_interval_lat))) + 1
    latide = int(floor((float(obs_vertex_id(4))-1.0)/float(obs_stat_interval_lat))) + 2
!!!!!debug
!if (pid ==0) print*, "lonids, lonide, latids, latide", lonids, lonide, latids, latide
!if (pid ==0) print*, (lonids-1)*obs_stat_interval_lon+1, (lonide-1)*obs_stat_interval_lon+1, (latids-1)*obs_stat_interval_lat, (latide-1)*obs_stat_interval_lat
!!!!!!!!!!

    ! get realfield lat lon
    allocate(real_lon(lonide-lonids+1))
    allocate(real_lat(latide-latids+1))
    if ((lonide-1)*obs_stat_interval_lon+1 > num_lon) then
        call RedisHmgetf1d(rc, hash_key_lon, "lon", (lonids-1)*obs_stat_interval_lon+1, num_lon, obs_stat_interval_lon, 1, real_lon(1:lonide-lonids))
        call RedisHmgetf1d(rc, hash_key_lon, "lon", num_lon, num_lon, 1, 1, real_lon(lonide-lonids+1))
    else
        call RedisHmgetf1d(rc, hash_key_lon, "lon", (lonids-1)*obs_stat_interval_lon+1, (lonide-1)*obs_stat_interval_lon+1, obs_stat_interval_lon, 1, real_lon)
    endif

    if ((latide-1)*obs_stat_interval_lat+1 > num_lat) then
        call RedisHmgetf1d(rc, hash_key_lat, "lat", (latids-1)*obs_stat_interval_lat+1, num_lat, obs_stat_interval_lat, 1, real_lat(1:latide-latids))
        call RedisHmgetf1d(rc, hash_key_lat, "lat", num_lat, num_lat, 1, 1, real_lat(latide-latids+1))
    else
        call RedisHmgetf1d(rc, hash_key_lat, "lat", (latids-1)*obs_stat_interval_lat+1, (latide-1)*obs_stat_interval_lat+1, obs_stat_interval_lat, 1, real_lat)
    endif
!!!!!debug
!if (pid ==0) print *,"real_lon",pid,real_lon,"real_lat",real_lat
!if (pid ==0) print *,"obs_lon", pid,obs_lon_t, "obs_lat", obs_lat_t
!!!!!!!!!!
    ! get realfield
    allocate(yo_tmp_2d(obs_num_2d, latide-latids+1, lonide-lonids+1))
    allocate(yo_tmp_3d(obs_num_3d+1, obs_l, latide-latids+1, lonide-lonids+1))

    if (obs_num_2d .ne. 0) then
        var_key = trim(adjustl(obs_var_name_list(1)))
        do var=2, obs_num_2d
            var_key = trim(adjustl(var_key)) //':'//trim(adjustl(obs_var_name_list(var)))
        end do

        call log_print("redis_hmget2d start")
        call log_print("da:"//to_str(obs_lon_num*obs_lat_num))

        call RedisHmgetf2d(rc, hash_key, var_key, lonids, lonide, 1, latids, latide, 1, obs_num_2d, yo_tmp_2d)
    endif

    var_key = 'p'
    if (obs_num_3d .ne. 0) then
        do var=obs_num_2d+1, obs_num_2d+obs_num_3d
            var_key = trim(adjustl(var_key)) //':'//trim(adjustl(obs_var_name_list(var)))
        end do
    endif
    
    call RedisHmgetf3d(rc, hash_key, var_key, lonids, lonide, 1, latids, latide, 1, &
         mpas_num_lev_start, mpas_num_lev, obs_vtc_interval, obs_num_3d+1, yo_tmp_3d)
!!!!!debug
!if (pid ==0) print*,"yo_tmp_3d(:,:,1,1)", yo_tmp_3d(:,:,1,1)
!if (pid ==0) print*,"yo_tmp_3d(:,:,1,2)", yo_tmp_3d(:,:,1,2)
!if (pid ==0) print*,"yo_tmp_3d(:,:,2,1)", yo_tmp_3d(:,:,2,1)
!if (pid ==0) print*,"yo_tmp_3d(:,:,2,2)", yo_tmp_3d(:,:,2,2)
!!!!!!!!!!
    do ix = 1, obs_lon_num
        do iy = 1, obs_lat_num
            idx = int(floor((float(obs_vertex_id(1)+(ix-1)*obs_hrz_interval_lon)-1.0)/float(obs_stat_interval_lon))) + 1 - lonids + 1
            idy = int(floor((float(obs_vertex_id(3)+(iy-1)*obs_hrz_interval_lat)-1.0)/float(obs_stat_interval_lat))) + 1 - latids + 1

            t1 = (real_lat(idy+1) - obs_lat_t(iy))/(real_lat(idy+1) - real_lat(idy))
            t2 = (real_lon(idx+1) - obs_lon_t(ix))/(real_lon(idx+1) - real_lon(idx))

            if (obs_num_2d .ne. 0) then
                do var=1, obs_num_2d
                    r1 = yo_tmp_2d(var, idy, idx) * t1 + &
                    yo_tmp_2d(var, idy+1, idx) * (1 - t1)
                    r2 = yo_tmp_2d(var, idy, idx+1) * t1 + &
                    yo_tmp_2d(var, idy+1, idx+1) * (1 - t1)
                    yo_tmp(var, iy ,ix) = r1 * t2 + r2 * (1 - t2)
                enddo
            endif

            do var=1,obs_num_3d+1
                do k=1,obs_l
                    r1 = yo_tmp_3d(var, k, idy, idx) * t1 + &
                    yo_tmp_3d(var, k, idy+1, idx) * (1 - t1)
                    r2 = yo_tmp_3d(var, k, idy, idx+1) * t1 + &
                    yo_tmp_3d(var, k, idy+1, idx+1) * (1 - t1)
                    yo_tmp(obs_num_2d+(var-1)*obs_l+k, iy ,ix) = r1 * t2 + r2 * (1 - t2)
                enddo
            enddo
        enddo
    enddo
!!!!!debug
!if (pid ==1) print*,"yo_tmp(:,1,1)",yo_tmp(:,1,1)
!!!!!!!!!!
    obs_p  = reshape(yo_tmp(obs_num_2d+1:obs_num_2d+obs_l,:,:), (/obs_l, obs_lat_num*obs_lon_num/))
!!!!!debug
!if (pid ==1) print*,"obs_p(:,1)",obs_p(:,1)
!!!!!!!!!!
    if (obs_num_2d .ne. 0) then
        do var = 1, obs_num_2d
            randn_1d = randn(obs_lon_num*obs_lat_num)
            iobs = 0
            do iy = 1, obs_lon_num
                do ix = 1, obs_lat_num
                    iobs = iobs + 1
!                    yo(i,iobs) = yo_tmp_2d(i, ix, iy)
                    yo(var,iobs) = yo_tmp(var, ix, iy) !+ randn_1d(iobs)/sqrt(obs_rinv(var))
                end do
            end do
        enddo
    endif

    do var = 1, obs_num_3d
        do k = 1, obs_l
            randn_1d = randn(obs_lon_num*obs_lat_num)
            iobs = 0
            do iy = 1, obs_lon_num
                do ix = 1, obs_lat_num
                    iobs = iobs + 1
                    yo(obs_num_2d+(var-1)*obs_l+k,iobs) = yo_tmp(obs_num_2d+var*obs_l+k,ix,iy) !&
!                    + randn_1d(iobs)/sqrt(obs_rinv(obs_num_2d+(var-1)*obs_l+k))
                end do
            end do
        enddo
    enddo
!!!!!debug
!if (pid ==1) print*,"yo(:,1)", yo(:,1)
!!!!!!!!!
130     format(20f16.8)
call log_print_root("OBS read finished")
do ix = 1, obs_l
    do iy = 1, obs_lat_num*obs_lon_num
        if (isnan(obs_p(ix,iy))) print *,"pid=",pid,ix,iy,"obs_p nan"
    enddo
enddo
do ix=1, obs_num_2d+obs_num_3d*obs_l
    do iy =1,obs_lat_num*obs_lon_num
        if (isnan(yo(ix,iy))) print *,"pid=",pid,ix,iy,"yo nan"
    enddo
enddo
!----------obs_output-----------!
!    obs_latids0 = obs_hrz_interval_lat-mod(1-mod(num_lat/2, obs_hrz_interval_lat),obs_hrz_interval_lat)+1
!    if (mod(1-mod(num_lat/2, obs_hrz_interval_lat),obs_hrz_interval_lat)<1) then
!        obs_latids0 = obs_latids0 - obs_hrz_interval_lat
!    endif
!
!    obs_latids = (obs_vertex_id(3)-obs_latids0)/obs_hrz_interval_lat+1
!    obs_latide = (obs_vertex_id(4)-obs_latids0)/obs_hrz_interval_lat+1
!    obs_lonids = (obs_vertex_id(1)-1)/obs_hrz_interval_lon+1
!    obs_lonide = (obs_vertex_id(2)-1)/obs_hrz_interval_lon+1
!
!    hash_key = "obs:" // trim(adjustl(case_name)) // ":"//trim(adjustl(num_lon_str)) &
!               // "x"// trim(adjustl(num_lat_str)) // ":" //trim(adjustl(date_char))
!
!    if (obs_num_3d .ne. 0) then
!
!        var_key = trim(adjustl(obs_var_name_list(obs_num_2d+1)))
!        do var=2, obs_num_3d
!            var_key = trim(adjustl(var_key)) //':'//trim(adjustl(obs_var_name_list(obs_num_2d+var)))
!        end do
!
!        do var=1, obs_num_3d
!            do k = 1, obs_l
!                iobs = 0
!                do iy = 1, obs_lon_num
!                    do ix = 1, obs_lat_num
!                        iobs = iobs + 1
!                        yo_out_3d(var,k,ix,iy) = yo(obs_num_2d+(var-1)*obs_l+k,iobs) 
!                    end do
!                end do
!            enddo
!        end do
!
!        call RedisHmsetf3d(rc, hash_key, var_key, obs_lonids, obs_lonide, 1, &
!             obs_latids, obs_latide, 1, 1, obs_l, 1, obs_num_3d, yo_out_3d)
!    end if
!-------------------------------!
end subroutine obs_get_yo_interp
end module obsmaker_mod

#undef SINGLE_POINT
