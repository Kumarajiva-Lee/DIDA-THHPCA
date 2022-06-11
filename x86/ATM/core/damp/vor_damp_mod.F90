module vor_damp_mod
  use mpi
  use flogger
  use string
  use const_mod
  use process_mod
  use namelist_mod
  use parallel_mod
  use parallel_types_mod
  use block_mod
  use tridiag_mod
  use operators_mod
  use pa_mod

  implicit none

  private

  public vor_damp_init
  public vor_damp_final
  public vor_damp_run

  real(r8), allocatable :: c_lon(:,:)
  real(r8), allocatable :: c_lat(:,:)

  real(r8) , allocatable :: rhs_all(:) , v_all(:)  ! send to solver , one member
  real(r8) , allocatable :: rhs_tran(:,:) , v_tran(:,:) ! mpi all member

  logical, allocatable :: use_implicit_solver(:)
  real(r8), parameter :: beta = 0.5_r8
  type(tridiag_solver_type), allocatable :: zonal_solver(:,:)

contains

  subroutine vor_damp_init(blocks)

    type(block_type), intent(in) :: blocks(:)

    integer j, j0, jr, k
    integer iblk, js, je
    real(r8) a, b
    logical :: trag_flag = .false.


    if (.not. use_vor_damp) return

    call vor_damp_final()

    ! Only do vorticity damping in reduced regions.
    ! First, find the interface when reduce starts.
    j0 = 0
    do j = global_mesh%full_lat_ibeg, global_mesh%full_lat_iend
      if (global_mesh%full_lat(j) < 0 .and. -global_mesh%full_lat_deg(j) >= vor_damp_lat0) j0 = j
    end do

    allocate(c_lon(global_mesh%num_full_lat,global_mesh%num_full_lev))
    allocate(c_lat(global_mesh%num_half_lat,global_mesh%num_full_lev))

    

    select case (vor_damp_order)
    case (2)
      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          jr = merge(j - global_mesh%full_lat_ibeg_no_pole + 1, global_mesh%full_lat_iend_no_pole - j + 1, global_mesh%full_lat(j) < 0)
          if (j0 == 0) then
            c_lon(j,k) = vor_damp_coef2 * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
          else
            c_lon(j,k) = vor_damp_coef2 * &
              exp(jr**2 * log(1.0e-10) / j0**2) * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
          end if
        end do
      end do

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
          jr = merge(j - global_mesh%half_lat_ibeg + 1, global_mesh%half_lat_iend - j + 1, global_mesh%half_lat(j) < 0)
          if (j0 == 0) then
            c_lat(j,k) = vor_damp_coef2 * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
          else
            c_lat(j,k) = vor_damp_coef2 * &
              exp(jr**2 * log(1.0e-10) / j0**2) * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
          end if
        end do
      end do
    case default
      call log_error('Unsupported vor_damp_order ' // trim(to_str(vor_damp_order)) // '!')
    end select

    ! Initialize cyclic tridiagonal solvers on each zonal circles if need implicit integration.
    js =  1e8
    je = -1e8
    do iblk = 1, size(proc%blocks)
      js = min(proc%blocks(iblk)%mesh%half_lat_ibeg_no_pole, js)
      je = max(proc%blocks(iblk)%mesh%half_lat_iend_no_pole, je)
    end do
    allocate(use_implicit_solver(js:je))
    use_implicit_solver = .false.
    allocate(zonal_solver(js:je,global_mesh%num_full_lev))
    do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
      do j = js, je
        if (abs(global_mesh%half_lat_deg(j)) > vor_damp_imp_lat0) then
          use_implicit_solver(j) = .true.
          trag_flag = .true.
          if (k > 1) then
            if (c_lat(j,k) == c_lat(j,k-1)) then
              call zonal_solver(j,k)%clone(zonal_solver(j,k-1))
              cycle
            end if
          end if
          b = -c_lat(j,k) * (1 - beta) * dt_in_seconds / global_mesh%le_lat(j)**2
          a = 2 * (-b) + 1
          if (zonal_tridiag_solver == 'mkl') then
            call zonal_solver(j,k)%init_sym_const(global_mesh%num_full_lon, a, b, zonal_tridiag_solver)
            ! write(*,*) proc%id, j, k, zonal_solver(j,k)%solver
          else
            call zonal_solver(j,k)%init_sym_const(blocks(1)%mesh%num_full_lon, a, b, zonal_tridiag_solver)
            ! write(*,*) proc%id, j, k, zonal_solver(j,k)%solver
          endif
        end if
      end do
    end do


    if (zonal_tridiag_solver == 'mkl') then
      if (trag_flag) then
        allocate(rhs_all(global_mesh%num_full_lon))
        allocate(v_all(global_mesh%num_full_lon))
        allocate(rhs_tran(member_num ,global_mesh%num_full_lon)) 
        allocate(v_tran(member_num , global_mesh%num_full_lon))
      end if
    endif

  end subroutine vor_damp_init

  subroutine vor_damp_final()

    if (allocated(c_lon)) deallocate(c_lon)
    if (allocated(c_lat)) deallocate(c_lat)

    if (allocated(rhs_all))  deallocate(rhs_all)
    if (allocated(rhs_tran)) deallocate(rhs_tran)
    if (allocated(v_all))    deallocate(v_all)
    if (allocated(v_tran))   deallocate(v_tran)

    if (allocated(use_implicit_solver)) deallocate(use_implicit_solver)
    if (allocated(zonal_solver)) deallocate(zonal_solver)

  end subroutine vor_damp_final

  subroutine vor_damp_run(block, dt, state)

    type(block_type), intent(in) :: block
    real(8), intent(in) :: dt
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    real(r8) rhs(member_num , block%mesh%full_lon_ibeg:block%mesh%full_lon_iend)
    integer status(MPI_STATUS_SIZE), ierr
    integer i, j, k, im
    integer proc_length
    integer tmp_dims(2) , tmp_coords(2)

    associate (mesh => block%mesh, &
      vor  => state%vor , &
      u    => state%u   , &
      v    => state%v)

    !!!!!!!!    
    call calc_vor(block, state, u, v)
#ifdef Detail_Time
    call Get_Start_Time(tran_time_start)
#endif

    call state%async(async_vor)%wait()

#ifdef Detail_Time
    call Get_End_Time(tran_time_end)
    tran_time = tran_time + tran_time_end - tran_time_start
#endif

    if (partition_type  == 'regular') then
      tmp_dims = proc%cart_dims
      tmp_coords = proc%cart_coords
    else
      tmp_dims = proc%irr_comm%cart_dims
      tmp_coords = proc%irr_comm%cart_coords
    end if

    proc_length = global_mesh%num_half_lon / tmp_dims(1)



    select case (vor_damp_order)
    case (2)
#ifdef Detail_Time
      call Get_Start_Time(tran_time_start)
#endif
      call state%async(async_u)%wait()
#ifdef Detail_Time
      call Get_End_Time(tran_time_end)
      tran_time = tran_time + tran_time_end - tran_time_start
#endif
#ifdef Detail_Time
      call Get_Start_Time(cal_time_start)
#endif
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            u(:,i,j,k) = u(:,i,j,k) - dt * c_lon(j,k) * ( &
              vor(:,i,j,k) - vor(:,i,j-1,k)) / mesh%le_lon(j)
          end do
        end do
      end do
#ifdef Detail_Time
      call Get_End_Time(cal_time_end)
      cal_time = cal_time + cal_time_end - cal_time_start
#endif
#ifdef Detail_Time
      call Get_Start_Time(tran_time_start)
#endif
      call fill_halo_member(block, state%u, full_lon=.false., full_lat=.true., full_lev=.true., async=state%async(async_u))
#ifdef Detail_Time
      call Get_End_Time(tran_time_end)
      tran_time = tran_time + tran_time_end - tran_time_start
#endif
      ! call fill_halo_member(block, state%u, full_lon=.false., full_lat=.true., full_lev=.true.)

#ifdef Detail_Time
      call Get_Start_Time(tran_time_start)
#endif
      call state%async(async_v)%wait()
      if (use_implicit_solver(mesh%half_lat_ibeg_no_pole) == .true. .or. &
          use_implicit_solver(mesh%half_lat_iend_no_pole) == .true.) call state%async(async_u)%wait() ! 如果是隐式耗散需要等待用新的u

#ifdef Detail_Time
      call Get_End_Time(tran_time_end)
      tran_time = tran_time + tran_time_end - tran_time_start
#endif

#ifdef Detail_Time
            call Get_Start_Time(cal_time_start)
#endif
      
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg, mesh%half_lat_iend
          if (use_implicit_solver(j)) then
            ! Set right hand side.

            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              rhs(:,i) = v(:,i,j,k) + &
                         c_lat(j,k) * beta * dt / mesh%le_lat(j) * (       &
                         vor(:,i,j,k) - vor(:,i-1,j,k)                 &
                         ) +                                                     &
                       c_lat(j,k) * (1 - beta) * dt / mesh%le_lat(j) * ( &
                         u(:,i-1,j+1,k) - u(:,i-1,j,k) -               &
                         u(:,i  ,j+1,k) + u(:,i  ,j,k)                 &
                       ) / mesh%de_lat(j)
            end do
            if (zonal_tridiag_solver == 'mkl') then
              if (tmp_dims(1) == 1) then
                do im = 1 , member_num
                  call zonal_solver(j,k)%solve(rhs(:,:), v(:, mesh%full_lon_ibeg:mesh%full_lon_iend,j,k) , im )
                end do
              else
                !all to 1
                if (tmp_coords(1) == 0) then 

                  rhs_tran(:,1:proc_length) = rhs

                  do i = 2 , tmp_dims(1) 
                    call MPI_RECV(rhs_tran(: , (i - 1) * proc_length + 1 : i  * proc_length ) , proc_length * member_num , MPI_DOUBLE , proc%id + tmp_dims(2) * (i-1)  , 100 , proc%comm , status,ierr)
                  end do

                  do im = 1 , member_num
                    !rhs_all = rhs_tran(im , :)
                    call zonal_solver(j,k)%solve(rhs_tran, v_tran , im)
                    !v_tran(im , :) = v_all
                  end do

                  do i = 2 , tmp_dims(1) 
                    call MPI_SEND(v_tran(: ,(i - 1) * proc_length + 1 : i  * proc_length) , proc_length * member_num , MPI_DOUBLE , proc%id + tmp_dims(2) * (i-1) , 100 , proc%comm , status,ierr)
                  end do
                  v(:,mesh%full_lon_ibeg:mesh%full_lon_iend,j,k) = v_tran(: ,1:proc_length)
                else
                  call MPI_SEND(rhs , size(rhs) , MPI_DOUBLE , proc%id - tmp_dims(2) * tmp_coords(1) , 100 , proc%comm , ierr)
                  call MPI_RECV(v(:,mesh%full_lon_ibeg:mesh%full_lon_iend,j,k) , proc_length * member_num , MPI_DOUBLE , proc%id - tmp_dims(2) * tmp_coords(1) , 100 , proc%comm , status,ierr)
                end if
              end if
            else if (zonal_tridiag_solver == 'spk' ) then
              !do im = 1 , member_num
                call zonal_solver(j,k)%solve(rhs(:,:), v(: , mesh%full_lon_ibeg:mesh%full_lon_iend,j,k) )
              !end do
            else
              if (is_root_proc()) call log_error('Wrong  zonal_tridiag_solver input ' // trim(zonal_tridiag_solver) // '!')
            end if
          else
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              v(:,i,j,k) = v(:,i,j,k) + dt * c_lat(j,k) * ( &
                vor(:,i,j,k) - vor(:,i-1,j,k)) / mesh%le_lat(j)
            end do
          end if
        end do
      end do
    
#ifdef Detail_Time
      call Get_End_Time(cal_time_end)
      cal_time = cal_time + cal_time_end - cal_time_start
#endif


#ifdef Detail_Time
      call Get_Start_Time(tran_time_start)
#endif
      call fill_halo_member(block, state%v, full_lon=.true., full_lat=.false., full_lev=.true., async=state%async(async_v))
#ifdef Detail_Time
      call Get_End_Time(tran_time_end)
      tran_time = tran_time + tran_time_end - tran_time_start
#endif
      ! call fill_halo_member(block, state%v, full_lon=.true., full_lat=.false., full_lev=.true.)
   
    case (4)
    end select
    end associate

  end subroutine vor_damp_run

end module vor_damp_mod
