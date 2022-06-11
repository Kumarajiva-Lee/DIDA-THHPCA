! information of mesh grid points of a process for DA

module pro_info
    use coupler_config, only:num_lat, num_lon, atm_mpas_sceneid
    use redis_mod
    implicit none

    ! number of grid points in current process
    integer :: nx, ny
    integer :: global_bx, global_by
    real(4), allocatable :: grid_lat(:), grid_lon(:)
    real(4), allocatable :: grid_lat_global(:), grid_lon_global(:)
    integer :: s_lat, e_lat, s_lon, e_lon
    integer :: x_sec, x_res, y_sec, y_res

    interface latlon_dist
        procedure latlon_dist_r4
        procedure latlon_dist_r8
    end interface

    private nx, ny
    private global_bx, global_by
    private x_sec, x_res, y_sec, y_res

contains
    subroutine pro_info_init(rc, nx_in, ny_in, global_bx_in, global_by_in, &
        x_sec_in, x_res_in, y_sec_in, y_res_in)

        use redis_mod
        type(c_ptr), intent(in) :: rc
        integer, intent(in) :: nx_in, ny_in
        integer, intent(in) :: global_bx_in, global_by_in
        integer, intent(in) :: x_sec_in, x_res_in, y_sec_in, y_res_in

        integer :: i
        character(5) :: i_str, num_lon_str, num_lat_str
        character(256) :: hash_key_lat, hash_key_lon
        integer :: te_lon, te_lat

        !print *, 'read gmcore gird'
        nx = nx_in
        ny = ny_in

        global_bx = global_bx_in
        global_by = global_by_in

        x_sec = x_sec_in
        x_res = x_res_in
        y_sec = y_sec_in
        y_res = y_res_in
        
        call pro_info_get_indx(global_bx, global_by, s_lat, e_lat, s_lon, e_lon)

        write(num_lon_str, "(i5)") num_lon
        write(num_lat_str, "(i5)") num_lat

        if (e_lat .lt. num_lat) then
            te_lat = e_lat+1
        else
            te_lat = e_lat
        endif

        if (e_lon .lt. num_lon) then
            te_lon = e_lon+1
        else
            te_lon = e_lon
        endif

        allocate(grid_lat(1: nx+1))
        allocate(grid_lon(1: ny+1))

        allocate(grid_lat_global(1: num_lat+1))
        allocate(grid_lon_global(1: num_lon+1))

        call read_model_grid(rc, s_lat, te_lat, s_lon, te_lon, grid_lat, grid_lon)

        if (e_lat .eq. num_lat) grid_lat(nx+1) = 90
        if (e_lon .eq. num_lon) grid_lon(ny+1) = 360

        call read_model_grid(rc, 1, num_lat, 1, num_lon, grid_lat_global, grid_lon_global)

        grid_lat_global(num_lat+1) = 90
        grid_lon_global(num_lon+1) = 360

        ! hash_key_lat = "lat:" // trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str))
        ! hash_key_lon = "lon:" // trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str))
 
        ! call RedisHmgetd1d(rc, trim(adjustl(hash_key_lat)), "lat", s_lat, &
        !      te_lat, 1, 1, grid_lat)
        ! call RedisHmgetd1d(rc, trim(adjustl(hash_key_lon)), "lon", s_lon, &
        !      te_lon, 1, 1, grid_lon)

    end subroutine

    ! Given process block position, give start/end lat/lon index
    subroutine pro_info_get_indx(bx, by, tx_s, tx_e, ty_s, ty_e)!上一个模块的s_lat对应这些输出 qsm

        integer, intent(in) :: bx, by
        integer, intent(out):: tx_s, tx_e, ty_s, ty_e

        if(bx - 1 < x_res) then
            tx_s = (bx - 1) * (x_sec + 1) + 1
            tx_e = tx_s + x_sec
        else
            tx_s = (bx - 1) * x_sec + x_res + 1
            tx_e = tx_s + x_sec - 1
        end if
    
        if(by - 1 < y_res) then
            ty_s = (by - 1) * (y_sec + 1) + 1
            ty_e = ty_s + y_sec
        else
            ty_s = (by - 1) * y_sec + y_res + 1
            ty_e = ty_s + y_sec - 1
        end if

        ! tx_s = (bx - 1) * nx + 1
        ! tx_e = bx * nx
        ! ty_s = (by - 1) * ny + 1
        ! ty_e = by * ny

        ! if (tx_e .lt. num_lat) tx_e = tx_e + 1

        ! if (ty_e .lt. num_lon) ty_e = ty_e + 1
    end subroutine

    subroutine read_model_grid(rc, lats, late, lons, lone, lat_grid, lon_grid)
        implicit none
        type(c_ptr), intent(in)       :: rc
        integer, intent(in)           :: lats, late, lons, lone
        real(4), intent(out)          :: lat_grid(late-lats+1)
        real(4), intent(out)          :: lon_grid(lone-lons+1)
        character(5)                  :: num_lon_str, num_lat_str
        character(:),allocatable      :: hkey
    
        ! read model grid
        write(num_lon_str, "(i5)") num_lon
        write(num_lat_str, "(i5)") num_lat
    
        hkey = "lat:" // trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str))
        call RedisHmgetf1d(rc, hkey, "lat", lats, late, 1, 1, lat_grid)
    
        hkey = "lon:" // trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str))
        call RedisHmgetf1d(rc, hkey, "lon", lons, lone, 1, 1, lon_grid)
    
    end subroutine read_model_grid

#define a 6371.393
    real(8) function latlon_dist_r8(x1, y1, x2, y2)

        real(8) :: x1, y1, x2, y2
        real(8) :: pi

        intrinsic atan
        
        pi = atan(1.0) * 4
        latlon_dist_r8 = a*acos(min(cos(x1/180*pi)*cos(x2/180*pi)*cos(y1/180*pi-y2/180*pi)+sin(x1/180*pi)*sin(x2/180*pi),1.0))

        return

    end function

    real(4) function latlon_dist_r4(x1, y1, x2, y2)

        real(4) :: x1, y1, x2, y2
        real(4) :: pi

        intrinsic atan
        
        pi = atan(1.0) * 4
        latlon_dist_r4 = a*acos(min(cos(x1/180*pi)*cos(x2/180*pi)*cos(y1/180*pi-y2/180*pi)+sin(x1/180*pi)*sin(x2/180*pi),1.0))

        return

    end function
#undef a

    ! real(8) function get_lat(grid_x)

    !     integer, intent(in) :: grid_x

    ! end function

    ! real(8) function get_lon(grid_y)

    !     integer, intent(in) :: grid_y
    
    ! end function

    subroutine pro_info_finalize()

        deallocate(grid_lat)
        deallocate(grid_lon)

    end subroutine

    subroutine pro_info_show()
        
        integer :: i, j

        do i = 1, nx
            do j = 1, ny
                write(*, 300) grid_lat(i), grid_lon(j)
            end do
        end do
        300 format(F10.2, " ", F10.2\)
    end subroutine

end module
