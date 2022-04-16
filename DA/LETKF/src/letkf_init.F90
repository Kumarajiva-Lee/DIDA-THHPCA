! letkf initialization
subroutine letkf_init(nx_in, ny_in, n_in, pid_in, &
    global_bx_in, global_by_in, global_nx_in, global_ny_in, &
    comm_letkf_in, letkf_rc_in, pos2id)

    !use namelist_mod
    use da_namelist_mod
    use coupler_config
    use var_info_mod
    use splitchar
    !use mpal

    integer, intent(in) :: nx_in, ny_in, n_in, pid_in
    integer, intent(in) :: global_bx_in, global_by_in
    integer, intent(in) :: global_nx_in, global_ny_in
    integer, intent(in) :: comm_letkf_in
    type(c_ptr), intent(in) :: letkf_rc_in

    external pos2id

    character(20), dimension(10):: da_var_name_list_origin
    integer :: obs_num_2d, obs_num_3d

    integer :: i
    integer :: current_level

    ens = cnn_ensemble
    nx = nx_in
    ny = ny_in
    nz = num_lev
    halo = da_halo
    n = n_in

    inflation_factor = da_inflation_factor
    local_dist = da_loc_distance
    vert_dist = da_vert_distance

    pid = pid_in
    
    global_bx = global_bx_in
    global_by = global_by_in
    global_nx = global_nx_in
    global_ny = global_ny_in
    global_gx = num_lat
    global_gy = num_lon

    gc_filter = da_gc_enable
    ps_indx = ps_pos
    if (polar_discard) then
        obs_polar = da_polar_letkf_lat
    else
        obs_polar = 91.0
    end if

    call stringsplit(da_var_name, ',', da_var_name_list_origin, n_var)
    !allocate(da_var_list (n_var))
    allocate(da_var_start(n_var))
    allocate(da_var_lev  (n_var))
    
    call var_info(da_var_name_list_origin(1:n_var), n_var, da_var_list, obs_num_2d, obs_num_3d)

    current_level = 1
    do i = 1, obs_num_2d
        da_var_start(i) = current_level
        da_var_lev(i) = 1
        current_level = current_level + 1
    end do
    do i = 1, obs_num_3d
        da_var_start(obs_num_2d + i) = current_level
        da_var_lev(obs_num_2d + i) = nz
        current_level = current_level + nz
    end do

    cube_n = da_block_interval
    cube_nz = da_vertical_interval
    overlap = da_block_overlap
    precision = computing_precision

    comm_letkf = comm_letkf_in
    letkf_rc = letkf_rc_in

    call obs_init(ens)

    call comm_info_init(pos2id)

    call mpal_init(precision)

#if (DEBUG_PRINT == 1)
    !print *, "LETKF id", pid, "global_bx", global_bx, "global_by", global_by
    !print *, "LETKF id", pid, "global_nx", global_nx, "global_ny", global_ny
    !print *, "LETKF id", pid, "north_n", north_n, "south_n", south_n
    do i = -north_n, south_n
    !    print*, "LETKF id", pid, "lc_num", i, "west", obs_lc(i)%west_n, "east", obs_lc(i)%east_n
    end do
#endif

end subroutine

subroutine comm_info_init(pos2id)

    external :: pos2id

    integer :: i, j, k

    !integer :: x1, x2, y1, y2
    real(8) :: x1, x2, y1, y2
    integer :: tmp_bx, tmp_by
    integer :: tmp_lat
    integer :: ea_tmp
    real(8) :: lon_dist

    integer :: tx_s, tx_e, ty_s, ty_e

    ! decide processes for communication
    ! get north_n and south_n
    north_n = 0
    tmp_bx = global_bx - 1; tmp_by = global_by
    x1 = grid_lat_global(s_lat); y1 = grid_lon_global(s_lon); y2 = y1
    do while (tmp_bx > 0)
        call pro_info_get_indx(tmp_bx, tmp_by, tx_s, tx_e, ty_s, ty_e)
        x2 = grid_lat_global(tx_e + 1)
        ! if (pid == 1) print *, "haha", s_lat, s_lon, tx_e, s_lon
        ! if (pid == 1) print *, "haha", x1, y1, x2, y2
        ! if (pid == 1) print *, "haha", pid, latlon_dist(x1, y1, x2, y2), local_dist
        if (latlon_dist(x1, y1, x2, y2) > local_dist) exit
        north_n = north_n + 1
        tmp_bx = tmp_bx - 1
    end do
    tmp_bx = global_bx + 1; tmp_by = global_by
    x1 = grid_lat_global(e_lat + 1); y1 = grid_lon_global(s_lon); y2 = y1
    do while (tmp_bx <= global_nx)
        call pro_info_get_indx(tmp_bx, tmp_by, tx_s, tx_e, ty_s, ty_e)
        x2 = grid_lat_global(tx_s)
        if (latlon_dist(x1, y1, x2, y2) > local_dist) exit
        south_n = south_n + 1
        tmp_bx = tmp_bx + 1
    end do

    allocate(obs_lc(-north_n: south_n))

    do i = -north_n, south_n
        tmp_bx = global_bx + i

        ! west
        obs_lc(i)%west_n = 0
        tmp_by = global_by
        do while(obs_lc(i)%west_n <= (global_ny - 1) / 2)
            tmp_by = mod(tmp_by + global_ny - 2, global_ny) + 1
            if (nearest_dist(global_bx, global_by, tmp_bx, tmp_by) > local_dist) exit
            obs_lc(i)%west_n = obs_lc(i)%west_n + 1
        end do

        ! east
        obs_lc(i)%east_n = 0
        tmp_by = global_by
        do while(obs_lc(i)%east_n <= (global_ny - 1) / 2)
            tmp_by = mod(tmp_by, global_ny) + 1
            if (nearest_dist(global_bx, global_by, tmp_bx, tmp_by) > local_dist) exit
            obs_lc(i)%east_n = obs_lc(i)%east_n + 1
        end do

        if (obs_lc(i)%west_n + obs_lc(i)%east_n > global_ny - 1) obs_lc(i)%east_n = obs_lc(i)%east_n - 1
    end do

    ! do i = -north_n, south_n
    !     tmp_bx = global_bx + i
    !     call pro_info_get_indx(global_bx, global_by, tx_s, tx_e, ty_s, ty_e)
    !     if (tmp_bx <= global_nx / 2) then
    !         tmp_lat = tx_s
    !     else
    !         tmp_lat = tx_e + 1
    !     end if

    !     ! west
    !     tmp_by = global_by
    !     call pro_info_get_indx(tmp_bx, tmp_by, tx_s, tx_e, ty_s, ty_e)
    !     x1 = grid_lat_global(tmp_lat); y1 = grid_lon_global(ty_s); x2 = x1
    !     obs_lc(i)%west_n = 0
    !     do while(obs_lc(i)%west_n <= (global_ny - 1) / 2)
    !         tmp_by = mod(tmp_by + global_ny - 2, global_ny) + 1
    !         call pro_info_get_indx(tmp_bx, tmp_by, tx_s, tx_e, ty_s, ty_e)
    !         y2 = grid_lon_global(mod(ty_e, global_gy) + 1)
    !         if (latlon_dist(x1, y1, x2, y2) > local_dist) exit
    !         obs_lc(i)%west_n = obs_lc(i)%west_n + 1
    !         ! if (pid == 0) then
    !         !     print *, "LETKF id", pid, tmp_by, latlon_dist(x1, y1, x2, y2)
    !         !     print *, "x1y1x2y2", tmp_lat, ty_s, tmp_lat, mod(ty_e, global_gy) + 1
    !         ! end if
    !     end do

    !     ! east
    !     tmp_by = global_by
    !     call pro_info_get_indx(tmp_bx, tmp_by, tx_s, tx_e, ty_s, ty_e)
    !     x1 = grid_lat_global(tmp_lat); y1 = grid_lon_global(mod(ty_e, global_gy) + 1); x2 = x1
    !     obs_lc(i)%east_n = 0
    !     do while(obs_lc(i)%east_n <= (global_ny - 1) / 2)
    !         tmp_by = mod(tmp_by, global_ny) + 1
    !         call pro_info_get_indx(tmp_bx, tmp_by, tx_s, tx_e, ty_s, ty_e)
    !         y2 = grid_lon_global(ty_s)
    !         if (latlon_dist(x1, y1, x2, y2) > local_dist) exit
    !         obs_lc(i)%east_n = obs_lc(i)%east_n + 1
    !     end do
    !     if (obs_lc(i)%west_n + obs_lc(i)%east_n > global_ny - 1) obs_lc(i)%east_n = obs_lc(i)%east_n - 1
    ! end do

    do i = -north_n, south_n
        allocate(obs_lc(i)%obs_data(obs_type_num, -obs_lc(i)%west_n: obs_lc(i)%east_n))
        allocate(obs_lc(i)%nid(-obs_lc(i)%west_n: obs_lc(i)%east_n))
        do j = -obs_lc(i)%west_n, obs_lc(i)%east_n
            call pos2id(global_bx + i - 1, mod(global_by + j - 1 + global_ny, global_ny), obs_lc(i)%nid(j))
        end do
    end do

    ! decide processes to comm
    ! get north_n and south_n
    ! x1 = (global_bx - 1) * nx + 1; y1 = (global_by - 1) * ny + 1
    ! x2 = x1 - 1
    ! north_n = 0
    ! do while (x2 > 0 .and. grid_dist(x1, y1, x2, y1) < local_dist)
    !     north_n = north_n + 1
    !     x2 = x2 - nx
    ! end do
    ! x1 = global_bx * nx; y1 = global_by * ny
    ! x2 = x2 + 1
    ! south_n = 0
    ! do while (x2 < global_gx .and. grid_dist(x1, y1, x2, y1) < local_dist)
    !     south_n = south_n + 1
    !     x2 = x2 + nx
    ! end do

    ! allocate(obs_lc(-north_n: south_n))

    ! north
    ! x1 = (global_bx - 1) * nx + 1; y1 = ny
    ! x2 = x1 - 1; y2 = ny + 1
    ! do i = -1, -north_n, -1
    !     obs_lc(i)%east_n = 0
    !     do while (obs_lc(i)%east_n < global_ny / 2 .and. grid_dist(x1, y1, x2, y2) < local_dist)
    !         obs_lc(i)%east_n = 0
    !         y2 = y2 + ny
    !     end do
    !     if (obs_lc(i)%east_n == global_ny / 2) then
    !         obs_lc(i)%west_n = (global_ny - 1) / 2
    !     else
    !         obs_lc(i)%west_n = obs_lc(i)%east_n
    !     end if
    !     x2 = x2 - nx
    ! end do

    ! ! middle
    ! x1 = (global_bx - 1) * nx + 1;
    ! y1 = ny; y2 = ny + 1
    ! if (x1 > global_ny / 2 - ny / 2) then
    !     x1 = global_bx * nx
    ! end if
    ! obs_lc(0)%east_n = 0
    ! do while (obs_lc(0)%east_n < global_ny / 2 .and. grid_dist(x1, y1, x1, y2) < local_dist)
    !     obs_lc(0)%east_n = obs_lc(0)%east_n + 1
    !     y2 = y2 + ny
    ! end do
    ! if (obs_lc(0)%east_n == global_ny / 2) then
    !     obs_lc(0)%west_n = (global_ny - 1) / 2
    ! else
    !     obs_lc(0)%west_n = obs_lc(0)%east_n
    ! end if

    ! ! south
    ! x1 = global_bx * nx; y1 = ny
    ! x2 = x1 + 1; y2 = ny + 1
    ! do i = 1, south_n
    !     obs_lc(i)%east_n = 0
    !     do while (obs_lc(i)%east_n < global_ny / 2 .and. grid_dist(x1, y1, x2, y2) < local_dist)
    !         obs_lc(i)%east_n = 0
    !         y2 = y2 + ny
    !     end do
    !     if (obs_lc(i)%east_n == global_ny / 2) then
    !         obs_lc(i)%west_n = (global_ny - 1) / 2
    !     else
    !         obs_lc(i)%west_n = obs_lc(i)%east_n
    !     end if
    !     x2 = x2 + nx
    ! end do

end subroutine


real(8) function nearest_dist(bx1, by1, bx2, by2)

    integer, intent(in) :: bx1, by1, bx2, by2

    integer :: xs1, ys1, xe1, ye1
    integer :: xs2, ys2, xe2, ye2
    integer :: tx1, ty1, tx2, ty2

    real(4) :: d1, d2
    real(4) :: lat1, lat2

    call pro_info_get_indx(bx1, by1, xs1, xe1, ys1, ye1)
    call pro_info_get_indx(bx2, by2, xs2, xe2, ys2, ye2)

    if (bx1 > bx2) then
        tx1 = xs1
        tx2 = xe2 + 1
    else if (bx1 < bx2) then
        tx1 = xe1 + 1
        tx2 = xs2
    else
        if (bx1 > global_nx / 2) then
            tx1 = xe1 + 1
            tx2 = xe2 + 1
        else
            tx1 = xs1
            tx2 = xs2
        end if
    end if

    lat1 = grid_lat_global(tx1)
    lat2 = grid_lat_global(tx2)

    if (abs(lat1) > obs_polar) lat1 = sign(real(obs_polar, 4), lat1)
    if (abs(lat2) > obs_polar) lat2 = sign(real(obs_polar, 4), lat2)

    if (by1 .eq. by2) then
        ty1 = ys1
        ty2 = ys2
        nearest_dist = latlon_dist(lat1, grid_lon_global(ty1), &
        lat2, grid_lon_global(ty2))
    else
        ty1 = ys1
        ty2 = mod(ye2, global_gy) + 1
        d1 = latlon_dist(lat1, grid_lon_global(ty1), &
        lat2, grid_lon_global(ty2))
        ty1 = mod(ye1, global_gy) + 1
        ty2 = ys2
        d2 = latlon_dist(lat1, grid_lon_global(ty1), &
        lat2, grid_lon_global(ty2))
        nearest_dist = min(d1, d2)
    end if

end function
