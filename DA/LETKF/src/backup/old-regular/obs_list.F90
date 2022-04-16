#if (yb_size == 8)
#define YB_MPI_TYPE MPI_REAL8
#endif
#if (yb_size == 4)
#define YB_MPI_TYPE MPI_REAL4
#endif

#if (ayb_size == 8)
#define AYB_MPI_TYPE MPI_REAL8
#endif
#if (ayb_size == 4)
#define AYB_MPI_TYPE MPI_REAL4
#endif

! Generate observation list and communicate with neighbours

subroutine gen_obslist(ntimes, xb, axb)

    use pro_info, only : grid_lat, grid_lon
    use da_namelist_mod, only:obs_stat_interval_lon, obs_stat_interval_lat

    integer,  intent(in) :: ntimes
    xb_type,  intent(in) :: xb(ens, n, 1-halo:nx+halo, 1-halo:ny+halo)
    axb_type, intent(in) :: axb(n, 1-halo:nx+halo, 1-halo:ny+halo)

    integer :: i, j, k, l
    integer :: ens_num, nn
    yb_type :: yb_t(ens)
    

    real(8) :: t_lat, t_lon
    real(8) :: t_obslat, t_obslon
    integer :: tx, ty
    integer :: tl, t_obsid
    integer :: tmp
    integer :: s_obs_gp, e_obs_gp

    real(8) :: sum_y

    if ((obs_stat_interval_lon .ne. 1) .or. (obs_stat_interval_lat .ne. 1)) then
        call obs_get_yo_interp(pid, letkf_rc, obs_vertex_id, obs_lat_t, obs_lon_t, obs_list(0)%lat_n, obs_list(0)%lon_n, ntimes, rinv,&
        obs_list(0)%ph, obs_list(0)%dis_yo(1:obs_nz, 1:total_obs))
    else
        call obs_get_yo(letkf_rc, obs_vertex_id, obs_list(0)%lat_n, obs_list(0)%lon_n, ntimes, rinv,&
        obs_list(0)%ph, obs_list(0)%dis_yo(1:obs_nz, 1:total_obs))
    endif

    t_obsid = 0
    do i = 1, obs_list(0)%lon_n
        do j = 1, obs_list(0)%lat_n

            ! grid point, index start from 1
            tx = obs_1lat + (j - 1) * obs_dlat
            ty = obs_1lon + (i - 1) * obs_dlon

            t_obsid = t_obsid + 1
            ! get obs info and info3d
            obs_list(0)%info(t_obsid)%lat = obs_lat_t(j)
            obs_list(0)%info(t_obsid)%lon = obs_lon_t(i)

            tmp = 0

            do k = 1, obs_n
                if (obs_var_lev(k) .ne. 1) then
                    tl = obs_var_lev(k)-obs_l*delta_l
                    do l = 1, obs_l
                        tl = tl + delta_l
                        tmp = tmp + 1
                        ! notice : this function needs to calculate actual value of observation
                        ! Not disturbation !!!
                        ! 3d not implemented yet
                        call obs_interpolation(ens, 1, xb(1:ens, obs_var_start(k)+tl-1, tx-1:tx, ty-1:ty), &
                        axb(obs_var_start(k)+tl-1, tx-1:tx, ty-1:ty), &
                        obs_list(0)%pro_yb(1:ens, tmp, t_obsid), &
                        grid_lat(tx-1), grid_lon(ty-1), grid_lat(tx), grid_lon(ty), &
                        obs_list(0)%info(t_obsid)%lat, obs_list(0)%info(t_obsid)%lon)
                    end do
                else
                    tmp = tmp + 1
                    call obs_interpolation(ens, 1, xb(1:ens, obs_var_start(k), tx-1:tx, ty-1:ty), &
                    axb(obs_var_start(k), tx-1:tx, ty-1:ty), &
                    obs_list(0)%pro_yb(1:ens, tmp, t_obsid), &
                    grid_lat(tx-1), grid_lon(ty-1), grid_lat(tx), grid_lon(ty), &
                    obs_list(0)%info(t_obsid)%lat, obs_list(0)%info(t_obsid)%lon)
                end if
            end do

        end do
    end do
    total_obs = t_obsid

#if (DEBUG == 1)
    if (tmp .ne. obs_nz) then
        call log_warning("obs interpolation number mismatch!")
    end if
#endif

    do i = 1, total_obs
        !do k = 1, obs_l
        do nn = 1, obs_nz
            sum_y = 0d0
            do ens_num = 1, ens
                sum_y = sum_y + obs_list(0)%pro_yb(ens_num, nn, i)
            end do
            sum_y = sum_y / ens
            obs_list(0)%dis_yo(nn, i) = obs_list(0)%dis_yo(nn, i) - sum_y
            ! convert actual value to disturbation
            do ens_num = 1, ens
                obs_list(0)%pro_yb(ens_num, nn, i) = obs_list(0)%pro_yb(ens_num, nn, i) - sum_y
            end do
        end do
        !end do
    end do
    
end subroutine

subroutine gen_obslist_3d(ntimes, xb, axb)

    use pro_info, only : grid_lat, grid_lon
    use da_namelist_mod, only:obs_stat_interval_lon, obs_stat_interval_lat

    integer,  intent(in) :: ntimes
    xb_type,  intent(in) :: xb(ens, n, 1-halo:nx+halo, 1-halo:ny+halo)
    axb_type, intent(in) :: axb(n, 1-halo:nx+halo, 1-halo:ny+halo)

    integer :: i, j, k, l
    integer :: ens_num, nn
    ayb_type :: xb_ph(ens, nz+1)
    real(8) :: dlat, dlon, t1, t2

    real(8) :: t_lat, t_lon
    real(8) :: t_obslat, t_obslon
    integer :: tx, ty
    integer :: tl, t_obsid
    integer :: tmp
    integer :: s_obs_gp, e_obs_gp

    real(8) :: sum_y

    if ((obs_stat_interval_lon .ne. 1) .or. (obs_stat_interval_lat .ne. 1)) then
        call obs_get_yo_interp(pid,letkf_rc, obs_vertex_id, obs_lat_t, obs_lon_t, obs_list(0)%lat_n, obs_list(0)%lon_n, ntimes, rinv,&
        obs_list(0)%ph, obs_list(0)%dis_yo(1:obs_nz, 1:total_obs))
    else
        call obs_get_yo(letkf_rc, obs_vertex_id, obs_list(0)%lat_n, obs_list(0)%lon_n, ntimes, rinv,&
        obs_list(0)%ph, obs_list(0)%dis_yo(1:obs_nz, 1:total_obs))
    endif

    t_obsid = 0
    do i = 1, obs_list(0)%lon_n
        do j = 1, obs_list(0)%lat_n

            ! grid point, index start from 1
            tx = obs_1lat + (j - 1) * obs_dlat
            ty = obs_1lon + (i - 1) * obs_dlon

            t_obsid = t_obsid + 1
            ! get obs info and info3d
            obs_list(0)%info(t_obsid)%lat = obs_lat_t(j)
            obs_list(0)%info(t_obsid)%lon = obs_lon_t(i)

            !if (pid == 0) print *, "begin ps2ph"

            dlat = grid_lat(tx) - grid_lat(tx-1)
            dlon = grid_lon(ty) - grid_lon(ty-1)
        
            t1 = (dlat - obs_list(0)%info(t_obsid)%lat + grid_lat(tx-1)) / dlat
            t2 = (dlon - obs_list(0)%info(t_obsid)%lon + grid_lon(ty-1)) / dlon

            call ps2ph(xb(1:ens, ps_pos, tx-1:tx, ty-1:ty), axb(ps_pos, tx-1:tx, ty-1:ty), xb_ph, t1, t2)

            call obs_interpolation_3d(ens, n, nz, obs_l, xb(1:ens, 1:n, tx-1:tx, ty-1:ty), &
            axb(1:n, tx-1:tx, ty-1:ty), xb_ph, obs_list(0)%ph(:, t_obsid), &
            obs_list(0)%pro_yb(1:ens, 1:obs_nz, t_obsid), t1, t2)

        end do
    end do

    total_obs = t_obsid

    do i = 1, total_obs
        !do k = 1, obs_l
        do nn = 1, obs_nz
            sum_y = 0d0
            do ens_num = 1, ens
                sum_y = sum_y + obs_list(0)%pro_yb(ens_num, nn, i)
            end do
            sum_y = sum_y / ens
            obs_list(0)%dis_yo(nn, i) = obs_list(0)%dis_yo(nn, i) - sum_y
            ! convert actual value to disturbation
            do ens_num = 1, ens
                obs_list(0)%pro_yb(ens_num, nn, i) = obs_list(0)%pro_yb(ens_num, nn, i) - sum_y
            end do
        end do
        !end do
    end do

#if (DEBUG_PRINT == 1)
    call log_notice("Done interpolation!")
#endif
    
end subroutine

subroutine comm_obslist()

    use mpi

    integer :: i, j, k, l
    integer :: ri, rj
    integer :: ierr
    integer, allocatable, dimension(:) :: request
    integer, allocatable, dimension(:, :) :: status
    integer, allocatable, dimension(:) :: r_request
    integer, allocatable, dimension(:, :) :: r_status
    integer, allocatable, dimension(:) :: halo_request
    integer, allocatable, dimension(:, :) :: halo_status
    integer, allocatable, dimension(:) :: random_map
    integer :: sendrecv_status(mpi_status_size)
    integer :: comm_start, comm_end, comm_max, r_comm_max, random_id, tag_base
    integer :: requests, total_requests, total_r_requests
    integer :: remote_sum
    integer :: t_indx
    integer :: remote_id, remote_recv_id, send_size, recv_size
    integer :: together
    integer :: msize_yb, msize_yo
    real :: random_tmp

    integer :: dx(4), dy(4), map(4), halo_map(6:9), nobs(9)

    integer :: yb_pack_pos, yb_pack_size
    integer :: ayb_pack_pos, ayb_pack_size
    integer :: local_lat, local_lon
    integer :: recv_list_num
    yb_type :: yb_send(obs_nz*ens, (OBS_DENSE+1)*(OBS_DENSE+1))
    ayb_type :: ayb_send(obs_nz, (OBS_DENSE+1)*(OBS_DENSE+1))

    dx(1) = 1; dy(1) = 0
    dx(2) = -1; dy(2) = 0
    dx(3) = 0; dy(3) = 1
    dx(4) = 0; dy(4) = -1
    map(1) = 2; map(2) = 1; map(3) = 4; map(4) = 3
    halo_map(6) = 1; halo_map(7) = 2; halo_map(8) = 2; halo_map(9) = 1

    total_requests = 0
    total_r_requests = 0

    !!! together : How many n_max(s) of observations are sent as a whole message !!!
    together = 256
    !!! ------------------------------------------------------------------------ !!!
    !!! msize_yb & msize_yo : message size of yb and yo                          !!!
    msize_yb = 4
    msize_yo = msize_yb * ens
    !!! ------------------------------------------------------------------------ !!!

    do j = 1, 4
        nobs(j) = obs_list(j)%lat_n * obs_list(j)%lon_n
        total_requests = total_requests + together * 4
        total_r_requests = total_r_requests + (rns(j) * remote_obsn) * 4
    end do
    do j = 6, 9
        nobs(j) = obs_list(j)%lat_n * obs_list(j)%lon_n
    end do

#if (DEBUG_BARRIER == 1)
    call mpi_barrier(comm_letkf, ierr)
#endif

#if (DEBUG == 1)
    call GPTLstart("Neighbour communication")
    call GPTLstart("Neighbour update-halo")
#endif

    do i = 1, 2
        send_size = obs_nz * obs_list(0)%lat_n * obs_list(0)%lon_n
        recv_size = obs_nz * obs_list(map(i))%lat_n * obs_list(map(i))%lon_n
        if (neighbour(i) .ne. -1 .and. neighbour(map(i)) .ne. -1) then
            call mpi_sendrecv(obs_list(0)%pro_yb(1, 1, 1), ens*send_size, YB_MPI_TYPE, neighbour(i), pid, &
            obs_list(map(i))%pro_yb(1, 1, 1), ens*recv_size, YB_MPI_TYPE, neighbour(map(i)), neighbour(map(i)), &
            comm_letkf, sendrecv_status, ierr)

            call mpi_sendrecv(obs_list(0)%dis_yo(1, 1), send_size, AYB_MPI_TYPE, neighbour(i), pid, &
            obs_list(map(i))%dis_yo(1, 1), recv_size, AYB_MPI_TYPE, neighbour(map(i)), neighbour(map(i)), &
            comm_letkf, sendrecv_status, ierr)
        end if
    end do

    do i = 3, 4
        send_size = obs_nz * obs_list(0)%lat_n * obs_list(0)%lon_n
        recv_size = obs_nz * obs_list(map(i))%lat_n * obs_list(map(i))%lon_n
        if (neighbour(i) == -1 .and. neighbour(map(i)) .ne. -1) then
            call mpi_recv(obs_list(map(i))%pro_yb(1, 1, 1), ens*recv_size, YB_MPI_TYPE, &
            neighbour(map(i)), neighbour(map(i)), comm_letkf, sendrecv_status, ierr)

            call mpi_recv(obs_list(map(i))%dis_yo(1, 1), recv_size, AYB_MPI_TYPE, &
            neighbour(map(i)), neighbour(map(i)), comm_letkf, sendrecv_status, ierr)
        end if
        if (neighbour(map(i)) == -1 .and. neighbour(i) .ne. -1) then
            call mpi_send(obs_list(0)%pro_yb(1, 1, 1), ens*send_size, YB_MPI_TYPE, &
            neighbour(i), pid, comm_letkf, ierr)

            call mpi_send(obs_list(0)%dis_yo(1, 1), send_size, AYB_MPI_TYPE, &
            neighbour(i), pid, comm_letkf, ierr)
        end if
        if (neighbour(i) .ne. -1 .and. neighbour(map(i)) .ne. -1) then
            call mpi_sendrecv(obs_list(0)%pro_yb(1, 1, 1), ens*send_size, YB_MPI_TYPE, neighbour(i), pid, &
            obs_list(map(i))%pro_yb(1, 1, 1), ens*recv_size, YB_MPI_TYPE, neighbour(map(i)), neighbour(map(i)), &
            comm_letkf, sendrecv_status, ierr)

            call mpi_sendrecv(obs_list(0)%dis_yo(1, 1), send_size, AYB_MPI_TYPE, neighbour(i), pid, &
            obs_list(map(i))%dis_yo(1, 1), recv_size, AYB_MPI_TYPE, neighbour(map(i)), neighbour(map(i)), &
            comm_letkf, sendrecv_status, ierr)
        end if
    end do

#if (DEBUG == 1)
    call GPTLstop("Neighbour update-halo")
#endif

#if (DEBUG_BARRIER == 1)
    call mpi_barrier(comm_letkf, ierr)
#endif

#if (DEBUG == 1)
    call GPTLstart("Neighbour update-obs_dense")
#endif

    if (neighbour_halo == .true.) then
        ! Give halo to north process
        do i = 1, 2
            if (i == 1) recv_list_num = 9
            if (i == 2) recv_list_num = 7
            remote_id = neighbour(4)
            remote_recv_id = neighbour(3)
            local_lat = min(obs_list(i)%lat_n, OBS_DENSE+1)
            local_lon = min(obs_list(i)%lon_n, OBS_DENSE+1)
            yb_pack_size = obs_nz*ens*local_lat*local_lon
            ayb_pack_size = obs_nz*local_lat*local_lon
            if (remote_id .ne. -1 .and. obs_list(i)%lat_n * obs_list(i)%lon_n >= 1) then
                yb_pack_pos = 0
                ayb_pack_pos = 0
                do j = 1, local_lon
                    do k = 1, local_lat
                        if (i == 1) then
                            t_indx = (j - 1) * obs_list(i)%lat_n + k
                        end if
                        if (i == 2) then
                            t_indx = (obs_list(i)%lon_n - j) * obs_list(i)%lat_n + k
                        end if
                        call mpi_pack(obs_list(i)%pro_yb(:, :, t_indx), ens*obs_nz, YB_MPI_TYPE, yb_send, &
                        yb_pack_size*yb_size, yb_pack_pos, comm_letkf, ierr)
                        call mpi_pack(obs_list(i)%dis_yo(:, t_indx), obs_nz, AYB_MPI_TYPE, ayb_send, &
                        ayb_pack_size*ayb_size, ayb_pack_pos, comm_letkf, ierr)
                    end do
                end do
            else
                remote_id = MPI_PROC_NULL
            end if
            if (remote_recv_id == -1) remote_recv_id = MPI_PROC_NULL
            call mpi_sendrecv(yb_send, yb_pack_size, YB_MPI_TYPE, remote_id, pid, &
            obs_list(recv_list_num)%pro_yb(1, 1, 1), nobs(recv_list_num)*ens*obs_nz, YB_MPI_TYPE, &
            remote_recv_id, remote_recv_id, comm_letkf, sendrecv_status, ierr)
            call mpi_sendrecv(ayb_send, ayb_pack_size, AYB_MPI_TYPE, remote_id, pid, &
            obs_list(recv_list_num)%dis_yo(1, 1), nobs(recv_list_num)*obs_nz, AYB_MPI_TYPE, &
            remote_recv_id, remote_recv_id, comm_letkf, sendrecv_status, ierr)
        end do
        ! give halo to south process
        do i = 1, 2
            if (i == 1) recv_list_num = 6
            if (i == 2) recv_list_num = 8
            remote_id = neighbour(3)
            remote_recv_id = neighbour(4)
            local_lat = min(obs_list(i)%lat_n, OBS_DENSE+1)
            local_lon = min(obs_list(i)%lon_n, OBS_DENSE+1)
            yb_pack_size = obs_nz*ens*local_lat*local_lon
            ayb_pack_size = obs_nz*local_lat*local_lon
            if (remote_id .ne. -1 .and. obs_list(i)%lat_n * obs_list(i)%lon_n >= 1) then
                yb_pack_pos = 0
                ayb_pack_pos = 0
                do j = 1, local_lon
                    do k = 1, local_lat
                        if (i == 1) then
                            t_indx = (j - 1) * obs_list(i)%lat_n + obs_list(i)%lat_n - k + 1
                        end if
                        if (i == 2) then
                            t_indx = (obs_list(i)%lon_n - j) * obs_list(i)%lat_n + obs_list(i)%lat_n - k + 1
                        end if
                        call mpi_pack(obs_list(i)%pro_yb(:, :, t_indx), ens*obs_nz, YB_MPI_TYPE, yb_send, &
                        yb_pack_size*yb_size, yb_pack_pos, comm_letkf, ierr)
                        call mpi_pack(obs_list(i)%dis_yo(:, t_indx), obs_nz, AYB_MPI_TYPE, ayb_send, &
                        ayb_pack_size*ayb_size, ayb_pack_pos, comm_letkf, ierr)
                    end do
                end do
            else
                remote_id = MPI_PROC_NULL
            end if
            if (remote_recv_id == -1) remote_recv_id = MPI_PROC_NULL
            call mpi_sendrecv(yb_send, yb_pack_size, YB_MPI_TYPE, remote_id, pid, &
            obs_list(recv_list_num)%pro_yb(1, 1, 1), nobs(recv_list_num)*ens*obs_nz, YB_MPI_TYPE, &
            remote_recv_id, remote_recv_id, comm_letkf, sendrecv_status, ierr)
            call mpi_sendrecv(ayb_send, ayb_pack_size, AYB_MPI_TYPE, remote_id, pid, &
            obs_list(recv_list_num)%dis_yo(1, 1), nobs(recv_list_num)*obs_nz, AYB_MPI_TYPE, &
            remote_recv_id, remote_recv_id, comm_letkf, sendrecv_status, ierr)
        end do
    end if

#if (DEBUG == 1)
    call GPTLstop("Neighbour update-obs_dense")
    call GPTLstop("Neighbour communication")
#endif

#if (DEBUG_BARRIER == 1)
    call mpi_barrier(comm_letkf, ierr)
#endif

    allocate(r_request(total_r_requests))
    allocate(r_status(mpi_status_size, total_r_requests))
    allocate(random_map(total_obs))

#if (DEBUG == 1)
    call GPTLstart("Remote observation")
#endif

    requests = 0
    remote_sum = 0
    do i = 1, 4
        do j = 1, rns(i)
            do k = 1, remote_obsn
                send_size = ens * obs_nz
                requests = requests + 1
                remote_sum = remote_sum + 1
                call mpi_irecv(obs_list(5)%pro_yb(:, :, remote_sum), send_size, YB_MPI_TYPE, r_neighbour(j, i), 2*k, comm_letkf, r_request(requests), ierr)
                
                send_size = obs_nz
                requests = requests + 1
                call mpi_irecv(obs_list(5)%dis_yo(:, remote_sum), send_size, AYB_MPI_TYPE, r_neighbour(j, i), 2*k+1, comm_letkf, r_request(requests), ierr)
            end do
        end do
        do j = 1, rns(map(i))
            random_map(:) = 0
            do k = 1, remote_obsn
                call random_number(random_tmp)
                random_id = floor(random_tmp * total_obs) + 1
                do while (random_map(random_id) .ne. 0)
                    random_id = mod(random_id, total_obs) + 1
                end do
                random_map(random_id) = 1
                send_size = ens * obs_nz
                requests = requests + 1
                call mpi_isend(obs_list(0)%pro_yb(:, :, random_id), send_size, YB_MPI_TYPE, r_neighbour(j, map(i)), 2*k, comm_letkf, r_request(requests), ierr)

                send_size = obs_nz
                requests = requests + 1
                call mpi_isend(obs_list(0)%dis_yo(:, random_id), send_size, AYB_MPI_TYPE, r_neighbour(j, map(i)), 2*k+1, comm_letkf, r_request(requests), ierr)
            end do
        end do
    end do

    if (requests > 0) call mpi_waitall(requests, r_request(1:requests), r_status(:, 1:requests), ierr)

    deallocate(r_request)
    deallocate(r_status)
    deallocate(random_map)

#if (DEBUG == 1)
    call GPTLstop("Remote observation")
#endif

#if (DEBUG_BARRIER == 1)
    call mpi_barrier(comm_letkf, ierr)
#endif

#if (DEBUG_PRINT == 1)
    call log_notice("Obslist communication done!")
#endif

end subroutine

#define obs_list_cpy(direct) \
yb_for_cpy(:, :, j, i) = obs_list(direct)%pro_yb(:, :, tindx);\
ayb_for_cpy(:, j, i) = obs_list(direct)%dis_yo(:, tindx)

#define get_tmp_latlon(direct) \
tmp_lat = obs_list(direct)%lat_n; tmp_lon = obs_list(direct)%lon_n

subroutine gen_obs_for_cpy()

    integer :: i, j
    integer :: tmp_lat, tmp_lon
    integer :: mylat, mylon
    integer :: tindx

    yb_for_cpy(:, :, :, :) = 0
    ayb_for_cpy(:, :, :) = 0

    mylat = obs_list(0)%lat_n
    mylon = obs_list(0)%lon_n

    ! if (pid == 0) then
    !     do i = 0, 9
    !         print *, "obs_list", i, obs_list(i)%lat_n, obs_list(i)%lon_n
    !     end do
    ! end if

    get_tmp_latlon(8)
    do i = 1 - tmp_lon, 0
        do j = 1 - tmp_lat, 0
            tindx = (i + tmp_lon - 1) * tmp_lat + j + tmp_lat
            obs_list_cpy(8)
        end do
    end do
    get_tmp_latlon(2)
    do i = 1 - min(OBS_DENSE + 1, tmp_lon), 0
        do j = 1, tmp_lat
            tindx = (i + tmp_lon - 1) * tmp_lat + j
            obs_list_cpy(2)
        end do
    end do
    get_tmp_latlon(7)
    do i = 1 - tmp_lon, 0
        do j = mylat + 1, mylat + tmp_lat
            tindx = (i + tmp_lon - 1) * tmp_lat + j - mylat
            obs_list_cpy(7)
        end do
    end do

    get_tmp_latlon(4)
    do i = 1, tmp_lon
        do j = 1 - min(OBS_DENSE + 1, tmp_lat), 0
            tindx = (i - 1) * tmp_lat + j + tmp_lat
            obs_list_cpy(4)
        end do
    end do
    get_tmp_latlon(0)
    do i = 1, tmp_lon
        do j = 1, tmp_lat
            tindx = (i - 1) * tmp_lat + j
            obs_list_cpy(0)
        end do
    end do
    get_tmp_latlon(3)
    do i = 1, tmp_lon
        do j = mylat + 1, mylat + min(OBS_DENSE + 1, tmp_lat)
            tindx = (i - 1) * tmp_lat + j - mylat
            obs_list_cpy(3)
        end do
    end do

    get_tmp_latlon(6)
    do i = mylon + 1, mylon + tmp_lon
        do j = 1 - tmp_lat, 0
            tindx = (i - mylon - 1) * tmp_lat + j + tmp_lat
            obs_list_cpy(6)
        end do
    end do
    get_tmp_latlon(1)
    do i = mylon + 1, mylon + min(OBS_DENSE + 1, tmp_lon)
        do j = 1, tmp_lat
            tindx = (i - mylon - 1) * tmp_lat + j
            obs_list_cpy(1)
        end do
    end do
    get_tmp_latlon(9)
    do i = mylon + 1, mylon + tmp_lon
        do j = mylat + 1, mylat + tmp_lat
            tindx = (i - mylon - 1) * tmp_lat + j - mylat
            obs_list_cpy(9)
        end do
    end do
    
end subroutine

#undef obs_list_cpy
#undef get_tmp_latlon

#undef YB_MPI_TYPE
#undef AYB_MPI_TYPE
