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

! generate observation structure
subroutine gen_obs_struct(ntimes, xb, axb)

    integer,  intent(in) :: ntimes
    xb_type,  intent(in) :: xb(ens, n, 1-halo:nx+halo, 1-halo:ny+halo)
    axb_type, intent(in) :: axb(n, 1-halo:nx+halo, 1-halo:ny+halo)

    integer :: i, j, k
    integer :: nv, nl, ens_num

    real(8) :: t_lat, t_lon

    integer :: tx, ty
    real(8) :: ax, ay

    ayb_type :: xb_ph(ens, nz+1)
    ayb_type :: sum_y

    integer :: obs_vid(n_var)

    integer :: lats, late, lons, lone

#define self obs_lc(0)%obs_data(i, 0)

    do i = 1, obs_type_num
        if (if_regular) then
            call pro_info_get_indx(global_bx, global_by, lats, late, lons, lone)
            call allocate_obslist_regular(letkf_rc, lats, late, lons, lone, &
                                          ntimes, obs_type_list(i), self)
        else
            call allocate_obslist(letkf_rc, grid_lat(1), grid_lat(nx+1), &
            grid_lon(1), grid_lon(ny+1), ntimes, obs_type_list(i), self)
        endif

        ! print *, "LETKFhaha", pid, self%obs_num
        ! print *, "LETKFhaha", pid, grid_lat(1), grid_lat(nx+1), &
        ! grid_lon(1), grid_lon(ny+1)

        if (if_regular) then
            call read_obs_regular(letkf_rc, ntimes, lats, late, lons, lone, &
                                  self)
!print *, "qsm", lats, late, lons, lone, self%obs_type, self%obs_num, self%obs_nz, self%obs_nvar,&
!         self%obs_varname, self%lat(1), self%lon(1), self%ph(:,1), "yo", self%dis_yo(:,:,1), "rinv", self%obs_rinv(:,:,1)
        else
            call read_obs(letkf_rc, ntimes, grid_lat(1), grid_lat(nx+1), &
                          grid_lon(1), grid_lon(ny+1), self)
        endif
        

#if (DEBUG == 1)
        if (pid == 0 .and. isnan(self%dis_yo(1, 1, 1))) then
            print *, i, j, nv
            call log_error("nan redis")
        end if
#endif

        !if (pid == 0) print *, "LETKFhaha", self%ph(1,1), self%dis_yo(:,1,1),  self%obs_rinv(1,1,1)
    end do

    do i = 1, obs_type_num

        ! find obs_vid
        obs_vid = 0
        do k = 1, n_var
            do j = 1, self%obs_nvar
                if (trim(da_var_list(k)) .eq. trim(self%obs_varname(j))) then
                    obs_vid(j) = k
                end if
            !if (pid == 0) print *, "LETKF haha", pid, da_var_list(k), self%obs_varname(j)
            end do
        end do
#if (DEBUG == 1)
        do j = 1, self%obs_nvar
            if (obs_vid(j) == 0) then
                call log_error("observation variable name mismatching!")
            end if
        end do
#endif

        do j = 1, self%obs_num

            ! interpolation
            t_lat = self%lat(j)
            t_lon = self%lon(j)

            !if (pid == 0) print *, t_lat, t_lon

            call find_obs_grid(t_lat, t_lon, tx, ty, ax, ay)

            call ps2ph(xb(1:ens, ps_indx, tx-1:tx, ty-1:ty), axb(ps_indx, tx-1:tx, ty-1:ty), xb_ph, ax, ay)

            ! if (pid == 0) then
            !     print *, "txty", tx, ty
            !     print *, xb(1, 1, tx-1:tx, ty-1:ty)
            ! end if

            call obs_interpolation_3d_2(xb(1:ens, 1:n, tx-1:tx, ty-1:ty), &
            axb(1:n, tx-1:tx, ty-1:ty), xb_ph, self, obs_vid, j, ax, ay)

#if (DEBUG == 1)
            if (pid == 0 .and. isnan(self%pro_yb(1, 1, 1, j))) then
                print *, t_lat, t_lon
                print *, tx, ty, ps_indx, ax, ay
                print *, "ph", self%ph(1, j)
                print *, "ph", xb_ph(1, :)
                print *, "vid", obs_vid(1: self%obs_nvar)
                print *, i, j
                call log_error("interpolation nan")
            end if
#endif

        end do

    end do

    do i = 1, obs_type_num
        do j = 1, self%obs_num 
            do nl = 1, self%obs_nz
                do nv = 1, self%obs_nvar
                    sum_y = 0d0
                    do ens_num = 1, ens
                        sum_y = sum_y + self%pro_yb(ens_num, nv, nl, j)
                    end do
                    sum_y = sum_y / ens
                    self%dis_yo(nv, nl, j) = self%dis_yo(nv, nl, j) - sum_y
                    do ens_num = 1, ens
                        self%pro_yb(ens_num, nv, nl, j) = self%pro_yb(ens_num, nv, nl, j) - sum_y
                    end do
#if (DEBUG == 1)
                    if (pid == 0 .and. isnan(self%dis_yo(nv, nl, j))) then
                        print *, i, j, nv
                        call log_error("nan obs")
                    end if
#endif
                end do
            end do
        end do
    end do

#undef self

end subroutine

! find four grid points corresponding to obs node
! tx/ty: grid index
! ax/ay: interpolation 
subroutine find_obs_grid(lat, lon, tx, ty, ax, ay)

    use pro_info, only : grid_lat, grid_lon

    real(8), intent(in) :: lat, lon
    integer, intent(out) :: tx, ty
    real(8), intent(out) :: ax, ay

    real(8) :: dlat, dlon

    integer :: i, j

    tx = 0
    do i = 1, nx + 1
        if (grid_lat(i) >= lat) then
            tx = i
            exit
        end if
    end do
#if (DEBUG == 1)
    if (tx == 0 .or. grid_lat(tx) < lat) then
        print *, pid, "lat", tx, grid_lat(tx), lat
        print *, pid, "grid_lat", grid_lat(1: nx+1)
        call log_error("Outside interpolation!")
    endif
#endif

    ty = 0
    do j = 1, ny + 1
        if (grid_lon(j) >= lon) then
            ty = j
            exit
        end if
    end do
#if (DEBUG == 1)
    if (ty == 0 .or. grid_lon(ty) < lon) then
        print *, pid, "lon", ty, grid_lon(ty), lon
        print *, pid, "grid_lon", grid_lon(1: ny+1)
        call log_error("Outside interpolation!")
    endif
#endif

    if (tx > 1) then
        dlat = grid_lat(tx) - grid_lat(tx-1)
        ax = (dlat - lat + grid_lat(tx-1)) / dlat
    else
        ax = 0d0
    end if

    if (ty > 1) then
        dlon = grid_lon(ty) - grid_lon(ty-1)
        ay = (dlon - lon + grid_lon(ty-1)) / dlon
    else
        ay = 0d0
    end if

end subroutine

! observation communication
subroutine comm_obs_data()

    integer :: ti, tj, tk
    integer :: tx, ty
    integer :: i, j, k

    integer :: total_requests
    integer :: neighbour_id

    integer :: send_size
    integer :: recv_size  

    integer :: ierr
    integer :: requests
    integer, allocatable, dimension(:) :: request
    integer, allocatable, dimension(:, :) :: status

    ! 2(send/recv) * 5(lat + lon + pro_yb + dis_yo + obs_rinv) * obs_type_num
    total_requests = 0
    do i = -north_n, south_n
        do j = -obs_lc(i)%west_n, obs_lc(i)%east_n
            total_requests = total_requests + 2 * 5 * obs_type_num
        end do
    end do

    allocate(request(total_requests))
    allocate(status(mpi_status_size, total_requests))

    requests = 0
    do i = -north_n, south_n
        do j = -obs_lc(i)%west_n, obs_lc(i)%east_n
            if (i .ne. 0 .or. j .ne. 0) then

                ! tx = global_bx + i
                ! ty = mod(global_by + j - 1, global_ny) + 1
                !call pos2id(tx - 1, ty - 1, neighbour_id)

                neighbour_id = obs_lc(i)%nid(j)

                do k = 1, obs_type_num

                    requests = requests + 1

                    call mpi_irecv(obs_lc(i)%obs_data(k, j)%obs_num, 1, &
                    MPI_INTEGER, neighbour_id, k, comm_letkf, &
                    request(requests), ierr)

                    requests = requests + 1

                    call mpi_isend(obs_lc(0)%obs_data(k, 0)%obs_num, 1, &
                    MPI_INTEGER, neighbour_id, k, comm_letkf, &
                    request(requests), ierr)

                    !print *, "LETKF id:", pid, trim(obs_type_list(k)), obs_lc(0)%obs_data(k, 0)%obs_num

                end do
            end if
        end do
    end do

    call mpi_waitall(requests, request(1:requests), status(:, 1:requests), ierr)

    !call log_notice("obs num comm done!")

    ! allocate
    do i = -north_n, south_n
        do j = -obs_lc(i)%west_n, obs_lc(i)%east_n
            do k = 1, obs_type_num
                if (i .ne. 0 .or. j .ne. 0) then
                    if (if_regular) then
                        call allocate_obslist_by_num_regular(obs_lc(i)%obs_data(k, j)%obs_num, obs_lc(i)%obs_data(k, j))
                    else
                        call allocate_obslist_by_num(obs_type_list(k), obs_lc(i)%obs_data(k, j)%obs_num, obs_lc(i)%obs_data(k, j))
                        !print *, "LETKF id:", pid, i, j, k, trim(obs_type_list(k)), obs_lc(i)%obs_data(k, j)%obs_num
                    end if
                end if
            end do
        end do
    end do

#define dest obs_lc(i)%obs_data(k, j)
#define self obs_lc(0)%obs_data(k, 0)
    do k = 1, obs_type_num

    requests = 0
    do i = -north_n, south_n
        do j = -obs_lc(i)%west_n, obs_lc(i)%east_n
            if (i .ne. 0 .or. j .ne. 0) then

                ! tx = global_bx + i
                ! ty = mod(global_by + j - 1, global_ny) + 1
                ! call pos2id(tx - 1, ty - 1, neighbour_id)

                neighbour_id = obs_lc(i)%nid(j)

                !print *, "LETKF id", pid, "lc", i, "ea", j, neighbour_id

                !do k = 1, obs_type_num

                    requests = requests + 1
                    recv_size = dest%obs_num
                    call mpi_irecv(dest%lat(1), recv_size, &
                    MPI_REAL8, neighbour_id, k * 5 + 1, comm_letkf, request(requests), ierr)

                    requests = requests + 1
                    recv_size = dest%obs_num
                    call mpi_irecv(dest%lon(1), recv_size, &
                    MPI_REAL8, neighbour_id, k * 5 + 2, comm_letkf, request(requests), ierr)

                    requests = requests + 1
                    recv_size = dest%obs_num * dest%obs_nz * dest%obs_nvar
                    call mpi_irecv(dest%dis_yo(1, 1, 1), recv_size, &
                    AYB_MPI_TYPE, neighbour_id, k * 5 + 3, comm_letkf, request(requests), ierr)

                    requests = requests + 1
                    recv_size = dest%obs_num * dest%obs_nz * dest%obs_nvar * ens
                    call mpi_irecv(dest%pro_yb(1, 1, 1, 1), recv_size, &
                    YB_MPI_TYPE, neighbour_id, k * 5 + 4, comm_letkf, request(requests), ierr)

                    requests = requests + 1
                    recv_size = dest%obs_num * dest%obs_nz * dest%obs_nvar
                    call mpi_irecv(dest%obs_rinv(1, 1, 1), recv_size, &
                    MPI_REAL8, neighbour_id, k * 5 + 5, comm_letkf, request(requests), ierr)

                !end do
            end if
        end do
    end do

    do i = south_n, -north_n, -1
        do j = obs_lc(i)%east_n, -obs_lc(i)%west_n, -1

            if (i .ne. 0 .or. j .ne. 0) then

                neighbour_id = obs_lc(i)%nid(j)

                !do k = 1, obs_type_num

                    requests = requests + 1
                    send_size = self%obs_num
                    call mpi_isend(self%lat(1), send_size, &
                    MPI_REAL8, neighbour_id, k * 5 + 1, comm_letkf, request(requests), ierr)

                    requests = requests + 1
                    send_size = self%obs_num
                    call mpi_isend(self%lon(1), send_size, &
                    MPI_REAL8, neighbour_id, k * 5 + 2, comm_letkf, request(requests), ierr)

                    requests = requests + 1
                    send_size = self%obs_num * self%obs_nz * self%obs_nvar
                    call mpi_isend(self%dis_yo(1, 1, 1), send_size, &
                    AYB_MPI_TYPE, neighbour_id, k * 5 + 3, comm_letkf, request(requests), ierr)

                    requests = requests + 1
                    send_size = self%obs_num * self%obs_nz * self%obs_nvar * ens
                    call mpi_isend(self%pro_yb(1, 1, 1, 1), send_size, &
                    YB_MPI_TYPE, neighbour_id, k * 5 + 4, comm_letkf, request(requests), ierr)
                    
                    requests = requests + 1
                    send_size = self%obs_num * self%obs_nz * self%obs_nvar
                    call mpi_isend(self%obs_rinv(1, 1, 1), send_size, &
                    MPI_REAL8, neighbour_id, k * 5 + 5, comm_letkf, request(requests), ierr)

                !end do

            end if

        end do
    end do

#undef self
#undef dest

    call mpi_waitall(requests, request(1:requests), status(:, 1:requests), ierr)

    end do

end subroutine
