! ens, nx, ny, nz, n, obs_dlat_in, obs_dlon_in, : namelist
! halo, pid, bx, by, nx, ny : define
! cube, overlap : algorithm
! pos2id : map function
subroutine letkf_init(ens_in, nx_in, ny_in, nz_in, halo_in, n_in, &
    obs_dlat_in, obs_dlon_in, obs_dlev_in, polar_in, inflation_factor_in, &
    pid_in, global_bx_in, global_by_in, global_nx_in, global_ny_in, &
    cube_n_in, overlap_in, gc_enable_in, neighbour_halo_in, &
    da_local_distance_in, OBS_DENSE_in, pos2id, comm_in, rc_in)

#if (DEBUG == 1)
    use mpi
#endif
    use pro_info, only : s_lat, s_lon, pro_info_get_indx

    integer, intent(in) :: ens_in, nx_in, ny_in, nz_in, halo_in, n_in
    integer, intent(in) :: obs_dlat_in, obs_dlon_in, obs_dlev_in
    real(8), intent(in) :: polar_in, inflation_factor_in
    ! integer, intent(in) :: obs_n_in
    ! integer, intent(in) :: obs_var_start_in(obs_n_in)
    ! integer, intent(in) :: obs_var_lev_in(obs_n_in)
    integer, intent(in) :: pid_in
    integer, intent(in) :: global_bx_in, global_by_in
    integer, intent(in) :: global_nx_in, global_ny_in
    integer, intent(in) :: cube_n_in, overlap_in
    logical, intent(in) :: gc_enable_in, neighbour_halo_in
    real(8), intent(in) :: da_local_distance_in
    integer, intent(in) :: OBS_DENSE_in

    external :: pos2id
    integer, intent(in) :: comm_in
    type(c_ptr), intent(in):: rc_in

    integer :: i, j, k
    integer :: tx_s, tx_e, ty_s, ty_e
    integer :: tx, ty
    integer :: dx(0: 4), dy(0: 4)
    integer :: t_obslats, t_obslate, t_obslons, t_obslone
    integer :: obs_vertex_id_t(4)
    integer :: n_var
    integer :: n_start(10), n_end(10)
    integer, allocatable, dimension(:) :: tmp_obsn
    integer :: tmppid, remote_need, remote_dist1, remote_dist2
    integer :: remote_lat, remote_lon
    integer :: tmp_lat_pos(1: 4), tmp_lon_pos(1: 4)
#if (DEBUG == 1)
    character * (MPI_MAX_PROCESSOR_NAME) processor_name
    integer namelen, ierr
#endif

    ens  = ens_in
    nx   = nx_in
    ny   = ny_in
    nz   = nz_in
    halo = halo_in
    n    = n_in
    obs_dlat = obs_dlat_in
    obs_dlon = obs_dlon_in
    delta_l = obs_dlev_in
    pid = pid_in
    global_bx    = global_bx_in
    global_by    = global_by_in
    global_nx   = global_nx_in
    global_ny   = global_ny_in
    letkf_rc = rc_in
!    ps_pos = 1

    call obs_init(n_var, n_start, n_end)
    ! call obs_get_loc(rc, obs_lon_num, obs_lat_num, obs_vertex_id, obs_lon, obs_lat)
    call obs_get_num(s_lat, e_lat, s_lon, e_lon, &
    obs_list(0)%lat_n, obs_list(0)%lon_n, obs_vertex_id)
    allocate(obs_lat_t(obs_list(0)%lat_n))
    allocate(obs_lon_t(obs_list(0)%lon_n))
    call obs_get_loc(letkf_rc, obs_list(0)%lon_n, obs_list(0)%lat_n, &
    obs_vertex_id, obs_lon_t(:), obs_lat_t(:))
    call obs_get_first(obs_lat_t(1), obs_lon_t(1), obs_1lat, obs_1lon)
    !call obs_get_indx(s_lat, s_lon, obs_lat_t(1), obs_lon_t(1), obs_1gox, obs_1goy)

    polar = polar_in
    inflation_factor = inflation_factor_in

    cube_n  = cube_n_in
    overlap = overlap_in
    gc_enable = gc_enable_in
    neighbour_halo = neighbour_halo_in
    da_local_distance = da_local_distance_in
    OBS_DENSE = OBS_DENSE_in
    
    !obs_l   = obs_l_in
    obs_n   = n_var
    allocate(obs_var_start(obs_n))
    allocate(obs_var_lev(obs_n))

    obs_var_start(1:obs_n) = n_start(1:obs_n)
    obs_var_lev(1:obs_n) = n_end(1:obs_n) - n_start(1:obs_n) + 1

    obs_nz = 0
    do i = 1, obs_n
        if (obs_var_lev(i) == 1) then
            obs_nz = obs_nz + obs_var_lev(i)
        else
            obs_nz = obs_nz + obs_l
        end if
    end do

    comm_letkf = comm_in

    call obs_get_ndist(ndist)

    allocate(local_dist(0: ndist))
    allocate(tmp_obsn(max(global_nx, global_ny)))
    call obs_get_local_dist(ndist, local_dist)

    allocate(rinv(obs_num_2d+obs_num_3d*obs_l))
    call obs_get_r(letkf_rc, rinv)
    do i = 1, obs_num_2d+obs_num_3d*obs_l
        rinv(i) = 1.0/rinv(i)
    end do

    dx(1) = 0; dy(0) = 0
    dx(1) = 0; dy(1) = 1
    dx(2) = 0; dy(2) = -1
    dx(3) = 1; dy(3) = 0
    dx(4) = -1; dy(4) = 0
    
    neighbour(0) = pid
    allocate(obs_list(0)%info(obs_list(0)%lat_n*obs_list(0)%lon_n))
    allocate(obs_list(0)%ph(obs_l, obs_list(0)%lat_n*obs_list(0)%lon_n))
    allocate(obs_list(0)%pro_yb(ens, obs_nz, obs_list(0)%lat_n*obs_list(0)%lon_n))
    allocate(obs_list(0)%dis_yo(obs_nz, obs_list(0)%lat_n*obs_list(0)%lon_n))
    do i = 1, 4
        tx = global_bx + dx(i)
        ty = mod(global_by + dy(i) - 1 + global_ny, global_ny) + 1
        call pos2id(tx - 1, ty - 1, neighbour(i))
        if (neighbour(i) .ne. -1 .and. neighbour(i) .ne. pid) then
            !if (i .ne. 0) call pos2id(tx - 1, ty - 1, neighbour(i))
            call pro_info_get_indx(tx, ty, tx_s, tx_e, ty_s, ty_e)
            call obs_get_num(tx_s, tx_e, ty_s, ty_e, obs_list(i)%lat_n, obs_list(i)%lon_n, obs_vertex_id_t)
            !call obs_get_indx(tx_s, ty_s, t_obslats, t_obslons)
            !call obs_get_indx(tx_e + 1, ty_e + 1, t_obslate, t_obslone)
            !obs_list(i)%lat_n = t_obslate - t_obslats
            !obs_list(i)%lon_n = t_obslone - t_obslons
            !if (i == 0) neighbour(i) = 0
            allocate(obs_list(i)%info(obs_list(i)%lat_n*obs_list(i)%lon_n))
            allocate(obs_list(i)%ph(obs_l, obs_list(i)%lat_n*obs_list(i)%lon_n))
            allocate(obs_list(i)%pro_yb(ens, obs_nz, obs_list(i)%lat_n*obs_list(i)%lon_n))
            allocate(obs_list(i)%dis_yo(obs_nz, obs_list(i)%lat_n*obs_list(i)%lon_n))
        else
            neighbour(i) = -1 
            obs_list(i)%lat_n = 0
            obs_list(i)%lon_n = 0
        end if
    end do

    remote_comm = 0
    rns(:) = 0
    do i = 1, ndist
        do j = 1, 2
            if (obs_list(j)%lon_n < local_dist(i)) remote_comm = 1
        end do
        do j = 3, 4
            if (obs_list(j)%lat_n < local_dist(i)) remote_comm = 1
        end do
    end do

    remote_comm = 0
    rns(:) = 0
    do i = 1, ndist
        do j = 1, 2
            if (neighbour(j) .ne. -1 .and. obs_list(j)%lon_n < local_dist(i)) remote_comm = 1
        end do
        do j = 3, 4
            if (neighbour(j) .ne. -1 .and. obs_list(j)%lat_n < local_dist(i)) remote_comm = 1
        end do
    end do

    allocate(r_neighbour(ndist*2, 4))
    allocate(r_sn(ndist*2, 4))
    allocate(r_en(ndist*2, 4))
    rns(1:4) = 0
    !if (pid == 0) then
    !print *, "begin!!!"
    !print *, "obs_list", obs_list(0)%lon_n, obs_list(1)%lon_n, obs_list(2)%lon_n
    if (remote_comm == 1) then
        do j = 1, 2
            if (neighbour(j) .ne. -1) then
                do i = 1, ndist
                    tx = global_bx
                    ty = mod(global_by + dy(j) + global_ny - 1, global_ny) + 1
                    remote_dist1 = obs_list(0)%lon_n + 1; remote_dist2 = obs_list(0)%lon_n + obs_list(j)%lon_n
                    if (local_dist(i) > obs_list(j)%lon_n) then
                        do while (remote_dist1 - obs_list(0)%lon_n <= local_dist(i))
                            ty = mod(ty + dy(j) + global_ny - 1, global_ny) + 1
                            call pos2id(tx - 1, ty - 1, tmppid)
                            if (tmppid == pid) exit
                            call pro_info_get_indx(tx, ty, tx_s, tx_e, ty_s, ty_e)
                            call obs_get_num(tx_s, tx_e, ty_s, ty_e, &
                            remote_lat, remote_lon, obs_vertex_id_t)
                            remote_dist1 = remote_dist2 + 1
                            remote_dist2 = remote_dist2 + remote_lon
                            ! if (pid == 0) print *, tmppid, remote_dist1, remote_dist2, j, local_dist(i)
                            ! if ((remote_dist1 >= 1 + local_dist(i) .and. remote_dist1 <= obs_list(0)%lon_n + local_dist(i)) .or. &
                            ! remote_dist2 >= 1 + local_dist(i) .and. remote_dist2 <= obs_list(0)%lon_n + local_dist(i)) then
                            if (remote_dist1 <= obs_list(0)%lon_n + local_dist(i) .and. &
                            remote_dist2 >= 1 + local_dist(i)) then
                                !if (pid == 0) print *, "wuhu!", i, local_dist(i)
                                remote_need = 1
                                do k = 1, rns(j)
                                    if (r_neighbour(k, j) == tmppid) remote_need = 0
                                end do
                                if (remote_need == 1) then
                                    rns(j) = rns(j) + 1
                                    r_neighbour(rns(j), j) = tmppid
                                    r_sn(rns(j), j) = remote_dist1
                                    r_en(rns(j), j) = remote_dist2
                                end if
                            end if
                        end do
                    end if
                end do
            end if
        end do
        do j = 3, 4
            if (neighbour(j) .ne. -1) then
                do i = 1, ndist
                    tx = global_bx + dx(j)
                    ty = global_by
                    remote_dist1 = obs_list(0)%lat_n + 1; remote_dist2 = obs_list(0)%lat_n + obs_list(j)%lat_n
                    if (local_dist(i) > obs_list(j)%lat_n) then
                        do while (remote_dist1 - obs_list(0)%lat_n <= local_dist(i))
                            tx = tx + dx(j)
                            call pos2id(tx - 1, ty - 1, tmppid)
                            if (tmppid == pid .or. tmppid == -1) exit
                            call pro_info_get_indx(tx, ty, tx_s, tx_e, ty_s, ty_e)
                            call obs_get_num(tx_s, tx_e, ty_s, ty_e, &
                            remote_lat, remote_lon, obs_vertex_id_t)
                            remote_dist1 = remote_dist2 + 1
                            remote_dist2 = remote_dist2 + remote_lat
                            ! if ((remote_dist1 >= 1 + local_dist(i) .and. remote_dist1 <= obs_list(0)%lat_n + local_dist(i)) .or. &
                            ! remote_dist2 >= 1 + local_dist(i) .and. remote_dist2 <= obs_list(0)%lat_n + local_dist(i)) then
                            if (remote_dist1 <= obs_list(0)%lat_n + local_dist(i) .and. &
                            remote_dist2 >= 1 + local_dist(i)) then
                                remote_need = 1
                                do k = 1, rns(j)
                                    if (r_neighbour(k, j) == tmppid) remote_need = 0
                                end do
                                if (remote_need == 1) then
                                    rns(j) = rns(j) + 1
                                    r_neighbour(rns(j), j) = tmppid
                                    r_sn(rns(j), j) = remote_dist1
                                    r_en(rns(j), j) = remote_dist2
                                end if
                            end if
                        end do
                    end if
                end do
            end if
        end do
    end if
    !end if
    
    deallocate(tmp_obsn)

    !!! remote_obsn : How many observations are sent to remote neighbours        !!!
    remote_obsn = 1
    !!! ------------------------------------------------------------------------ !!!

    obs_list(5)%lat_n = (rns(1) + rns(2) + rns(3) + rns(4)) * remote_obsn
    !if (pid == 0) print *, "rns", rns(1), rns(2), rns(3), rns(4)
    allocate(obs_list(5)%info(obs_list(5)%lat_n))
    allocate(obs_list(5)%ph(obs_l, obs_list(5)%lat_n))
    allocate(obs_list(5)%pro_yb(ens, obs_nz, obs_list(5)%lat_n))
    allocate(obs_list(5)%dis_yo(obs_nz, obs_list(5)%lat_n))

    tmp_lat_pos(1) = 4; tmp_lon_pos(1) = 1
    tmp_lat_pos(2) = 3; tmp_lon_pos(2) = 2
    tmp_lat_pos(3) = 4; tmp_lon_pos(3) = 2
    tmp_lat_pos(4) = 3; tmp_lon_pos(4) = 1
    if (neighbour_halo == .true.) then
        do i = 6, 9
            obs_list(i)%lat_n = min(obs_list(tmp_lat_pos(i-5))%lat_n, OBS_DENSE + 1)
            obs_list(i)%lon_n = min(obs_list(tmp_lon_pos(i-5))%lon_n, OBS_DENSE + 1)
            if (obs_list(i)%lat_n * obs_list(i)%lon_n >= 1) then
                allocate(obs_list(i)%info(obs_list(i)%lat_n*obs_list(i)%lon_n))
                allocate(obs_list(i)%ph(obs_l, obs_list(i)%lat_n*obs_list(i)%lon_n))
                allocate(obs_list(i)%pro_yb(ens, obs_nz, obs_list(i)%lat_n*obs_list(i)%lon_n))
                allocate(obs_list(i)%dis_yo(obs_nz, obs_list(i)%lat_n*obs_list(i)%lon_n))
            end if
        end do
    end if

    allocate(yb_for_cpy(ens, obs_nz, 1-OBS_DENSE-1:obs_list(0)%lat_n+OBS_DENSE+1, &
    1-OBS_DENSE-1:obs_list(0)%lon_n+OBS_DENSE+1))
    allocate(ayb_for_cpy(obs_nz, 1-OBS_DENSE-1:obs_list(0)%lat_n+OBS_DENSE+1, &
    1-OBS_DENSE-1:obs_list(0)%lon_n+OBS_DENSE+1))

#if (DEBUG_PRINT == 1)
    call log_notice("Letkf init done!")
! if (pid <= 3) then
!     write(*, *) pid, "Letkf inited!"
!     write(*, 100) pid, global_bx, global_by, nx, ny, nz, n, cube_n
! 100 format (I2, " At position( ", I3, ", ", I3, " ) block size ", I3, " ", I3, " ", I3, " nmax = ", I3, " cube_n = ", I3)
!     write(*, *) pid, global_nx, global_ny
! !     write(*, 101) pid, obs_dlat, obs_dlon
! ! 101 format (I2, " obsd         ", I3, ",", I3)
! !     write(*, 102) pid, s_lat, s_lon
! ! 102 format (I2, " s_latlon     ", I3, ",", I3)
! !     write(*, 103) pid, obs_1lat, obs_1lon
! ! 103 format (I2, " obs_obs_1l     ", I3, ",", I3)
! !     write(*, 104) pid, obs_l, obs_nz
! ! 104 format (I2, " obs_obs_l ", I3, ", obs_nz", I3)
!     write(*, 105) pid, remote_comm
! 105 format(I2, " remote_comm ", I2)

    !  do i = 0, 4
    !      write(*, 200) pid, i, obs_list(i)%lat_n, obs_list(i)%lon_n, neighbour(i)
    !  end do
    
    !do i = 6, 9
    !    write(*, 202) pid, i, obs_list(i)%lat_n, obs_list(i)%lon_n, neighbour(3), neighbour(4)
    !end do

    ! do i = 1, 4
    !     do j = 1, rns(i)
    !         write(*, 201) pid, i, j, r_neighbour(j, i), r_sn(j, i), r_en(j, i)
    !     end do
    ! end do
    !200 format (I2, " obslist(", I1, ")", " obsn_lat ", I5, " obs_nlon ", I5, ' neighbour ', I5)
    !201 format (I2, " remote ", I3, "  ", I3, " r_neighbour ", I3, " r_sn ", I5, " r_en ", I5)
    !202 format (I2, " obslist(", I1, ")", " obsn_lat ", I5, " obs_nlon ", I5, ' neighbour ', I5, " ", I5)
! end if
#endif

#if (DEBUG == 1)
    call mpi_get_processor_name(processor_name, namelen, ierr)
    call log_notice("Letkf init done on " // processor_name)
#endif

end subroutine

subroutine letkf_finalize()

    integer :: i

    deallocate(obs_var_start)
    deallocate(obs_var_lev)
    deallocate(local_dist)
    deallocate(rinv)
    deallocate(obs_lat_t)
    deallocate(obs_lon_t)
    deallocate(r_neighbour)
    deallocate(r_sn)
    deallocate(r_en)
    do i = 0, 5
        if (obs_list(i)%lat_n >= 1) then
            deallocate(obs_list(i)%info)
            deallocate(obs_list(i)%ph)
            deallocate(obs_list(i)%pro_yb)
            deallocate(obs_list(i)%dis_yo)
        end if
    end do
    if (neighbour_halo == .true.) then
        do i = 6, 9
            if (obs_list(i)%lat_n >= 1) then
                deallocate(obs_list(i)%info)
                deallocate(obs_list(i)%ph)
                deallocate(obs_list(i)%pro_yb)
                deallocate(obs_list(i)%dis_yo)
            end if
        end do
    end if
    deallocate(yb_for_cpy)
    deallocate(ayb_for_cpy)
    
end subroutine
