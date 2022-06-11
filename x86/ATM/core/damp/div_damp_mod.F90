module div_damp_mod
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

  public div_damp_init
  public div_damp_final
  public div_damp_run

  real(r8), allocatable :: c_lon(:,:)
  real(r8), allocatable :: c_lat(:,:)
  real(r8) , allocatable :: rhs_all(:) , u_all(:)  ! send to solver , one member
  real(r8) , allocatable :: rhs_tran(:,:) , u_tran(:,:) ! mpi all member

  logical, allocatable :: use_implicit_solver(:)
  real(r8), parameter :: beta = 0.5_r8
  type(tridiag_solver_type), allocatable :: zonal_solver(:,:)

contains

  subroutine div_damp_init(blocks)

    type(block_type), intent(in) :: blocks(:)

    integer j, k, r, jr, j0
    integer iblk, js, je
    real(r8) a, b
    logical :: trag_flag = .false.

    if (.not. use_div_damp) return

    call div_damp_final()

    j0 = 0
    do j = global_mesh%full_lat_ibeg, global_mesh%full_lat_iend
      if (global_mesh%full_lat(j) < 0 .and. -global_mesh%full_lat_deg(j) >= div_damp_lat0) j0 = j
    end do

    allocate(c_lon(global_mesh%num_full_lat,global_mesh%num_full_lev))
    allocate(c_lat(global_mesh%num_half_lat,global_mesh%num_full_lev))


    select case (div_damp_order)
    case (2)
      r = 1

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          jr = merge(j - global_mesh%full_lat_ibeg_no_pole + 1, global_mesh%full_lat_iend_no_pole - j + 1, global_mesh%full_lat(j) < 0)
          if (baroclinic) then
            if (j0 == 0) then
              c_lon(j,k) = div_damp_coef2 * &
                (1.0_r8 + div_damp_top * exp((k-1)**2 * log(0.2_r8) / (div_damp_k0-1)**2)) * &
                global_mesh%full_cos_lat(j)**r * &
                radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
            else
              c_lon(j,k) = div_damp_coef2 * &
                (1.0_r8 + div_damp_top * exp((k-1)**2 * log(0.2_r8) / (div_damp_k0-1)**2)) * &
                (global_mesh%full_cos_lat(j)**r + div_damp_pole * exp(jr**2 * log(0.01_r8) / j0**2)) * &
                radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
            end if
          else
            c_lon(j,k) = div_damp_coef2 * &
              global_mesh%full_cos_lat(j)**r * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
          end if
        end do
      end do

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg, global_mesh%half_lat_iend
          jr = merge(j - global_mesh%half_lat_ibeg + 1, global_mesh%half_lat_iend - j + 1, global_mesh%half_lat(j) < 0)
          if (baroclinic) then
            if (j0 == 0) then
              c_lat(j,k) = div_damp_coef2 * &
                (1.0_r8 + div_damp_top * exp((k-1)**2 * log(0.2_r8) / (div_damp_k0-1)**2)) * &
                global_mesh%half_cos_lat(j)**r * &
                radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
            else
              c_lat(j,k) = div_damp_coef2 * &
                (1.0_r8 + div_damp_top * exp((k-1)**2 * log(0.2_r8) / (div_damp_k0-1)**2)) * &
                (global_mesh%half_cos_lat(j)**r + div_damp_pole * exp(jr**2 * log(0.01_r8) / j0**2)) * &
                radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
            end if
          else
            c_lat(j,k) = div_damp_coef2 * &
              global_mesh%half_cos_lat(j)**r * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
          end if
        end do
      end do
    case default
      call log_error('Unsupported div_damp_order ' // trim(to_str(div_damp_order)) // '!')
    end select

    ! Initialize cyclic tridiagonal solvers on each zonal circles if need implicit integration.

    js =  1e8
    je = -1e8
    do iblk = 1, size(proc%blocks)
      js = min(proc%blocks(iblk)%mesh%full_lat_ibeg_no_pole, js)
      je = max(proc%blocks(iblk)%mesh%full_lat_iend_no_pole, je)
    end do

    allocate(use_implicit_solver(js:je))

    use_implicit_solver = .false.
    allocate(zonal_solver(js:je,global_mesh%num_full_lev))
    do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
      do j = js, je
        if (abs(global_mesh%full_lat_deg(j)) > div_damp_imp_lat0) then
          use_implicit_solver(j) = .true.
          trag_flag = .true.
          if (k > 1) then
            if (c_lon(j,k) == c_lon(j,k-1)) then
              call zonal_solver(j,k)%clone(zonal_solver(j,k-1))
              cycle
            end if
          end if
          b = -c_lon(j,k) * (1 - beta) * dt_in_seconds / global_mesh%de_lon(j)**2
          a = 2 * (-b) + 1
          if (zonal_tridiag_solver == 'mkl') then
            call zonal_solver(j,k)%init_sym_const(global_mesh%num_half_lon, a, b,zonal_tridiag_solver)
          else
            call zonal_solver(j,k)%init_sym_const(blocks(1)%mesh%num_half_lon, a, b, zonal_tridiag_solver)
          endif
        end if
      end do
    end do

    if (zonal_tridiag_solver == 'mkl') then
      if (trag_flag) then
        allocate(rhs_all(global_mesh%num_half_lon))
        allocate(u_all(global_mesh%num_half_lon))
        allocate(rhs_tran(member_num ,global_mesh%num_half_lon))
        allocate(u_tran(member_num , global_mesh%num_half_lon))
      end if
    endif

  end subroutine div_damp_init

  subroutine div_damp_final()

    if (allocated(c_lon)) deallocate(c_lon)
    if (allocated(c_lat)) deallocate(c_lat)

    if (allocated(rhs_all)) deallocate(rhs_all)
    if (allocated(u_all)) deallocate(u_all)
    if (allocated(rhs_tran)) deallocate(rhs_tran)
    if (allocated(u_tran)) deallocate(u_tran)

    if (allocated(use_implicit_solver)) deallocate(use_implicit_solver)
    if (allocated(zonal_solver)) deallocate(zonal_solver)

  end subroutine div_damp_final

  subroutine div_damp_run(block, dt, state)

    type(block_type), intent(in) :: block
    real(8), intent(in) :: dt
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    real(r8) rhs(member_num , block%mesh%half_lon_ibeg:block%mesh%half_lon_iend)
    integer i, j, k , im
    integer proc_length
    integer status(MPI_STATUS_SIZE), ierr
    integer tmp_dims(2) , tmp_coords(2)

    mesh => state%mesh

    if (partition_type  == 'regular') then
      tmp_dims = proc%cart_dims
      tmp_coords = proc%cart_coords
    else
      tmp_dims = proc%irr_comm%cart_dims
      tmp_coords = proc%irr_comm%cart_coords
    end if

    proc_length = global_mesh%num_half_lon / tmp_dims(1)

    
    !maybe useless
    call calc_div(block, state)

#ifdef Detail_Time
    call Get_Start_Time(tran_time_start)
#endif

    call state%async(async_div)%wait()

#ifdef Detail_Time
    call Get_End_Time(tran_time_end)
    tran_time = tran_time + tran_time_end - tran_time_start
#endif

    select case (div_damp_order)
    case (2)

#ifdef Detail_Time
      call Get_Start_Time(tran_time_start)
#endif

      call state%async(async_v)%wait() ! 对v耗散

#ifdef Detail_Time
      call Get_End_Time(tran_time_end)
      tran_time = tran_time + tran_time_end - tran_time_start
      call Get_Start_Time(cal_time_start)
#endif

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg, mesh%half_lat_iend   
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            state%v(:,i,j,k) = state%v(:,i,j,k) + dt * c_lat(j,k) * ( &         
                               state%div(:,i,j+1,k) - state%div(:,i,j,k)) / mesh%de_lat(j)
          end do
        end do
      end do

#ifdef Detail_Time
      call Get_End_Time(cal_time_end)
      cal_time = cal_time + cal_time_end - cal_time_start
      call Get_Start_Time(tran_time_start)
#endif

      call fill_halo_member(block, state%v, full_lon=.true., full_lat=.false., full_lev=.true., async=state%async(async_v))

#ifdef Detail_Time
      call Get_End_Time(tran_time_end)
      tran_time = tran_time + tran_time_end - tran_time_start
      call Get_Start_Time(tran_time_start)
#endif
      
      call state%async(async_u)%wait() ! 对u耗散
      if (use_implicit_solver(mesh%full_lat_ibeg_no_pole) == .true. .or. &
          use_implicit_solver(mesh%full_lat_iend_no_pole) == .true.) call state%async(async_v)%wait() ! 如果是隐式耗散需要等待用新的v

#ifdef Detail_Time
      call Get_End_Time(tran_time_end)
      tran_time = tran_time + tran_time_end - tran_time_start
      call Get_Start_Time(cal_time_start)
#endif

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          if (use_implicit_solver(j)) then
            ! Set right hand side.
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              rhs(:,i) = state%u(:,i,j,k) + &
                         c_lon(j,k) * beta * dt / mesh%de_lon(j) * (       &
                         state%div(:,i+1,j,k) - state%div(:,i,j,k)                 &
                       ) +                                                     &
                       c_lon(j,k) * (1 - beta) * dt / mesh%de_lon(j) * ( &
                         state%v(:,i+1,j,k) - state%v(:,i+1,j-1,k) -               &
                         state%v(:,i  ,j,k) + state%v(:,i  ,j-1,k)                 &
                       ) / mesh%le_lon(j)
            end do


            if (zonal_tridiag_solver == 'mkl') then
              if (tmp_dims(1) == 1) then
                do im = 1 , member_num
                  call zonal_solver(j,k)%solve(rhs(:,:), state%u(:,mesh%half_lon_ibeg:mesh%half_lon_iend,j,k) , im)
                end do
              else  
                !all to 1 
                if (tmp_coords(1) == 0) then 
                  rhs_tran(:,1:proc_length) = rhs
                  
                  do i = 2 , tmp_dims(1) 
                    call MPI_RECV(rhs_tran(: , (i - 1) * proc_length + 1 : i  * proc_length ) , proc_length * member_num , MPI_DOUBLE , proc%id + tmp_dims(2) * (i-1)  , 100 , proc%comm , status,ierr)
                  end do

                  do im = 1 , member_num
                    !rhs_all(:) = rhs_tran(im , :)
                    call zonal_solver(j,k)%solve(rhs_tran, u_tran , im)
                    !u_tran(im , :) = u_all(:)
                  end do

                  do i = 2 , tmp_dims(1) 
                    call MPI_SEND(u_tran(: ,(i - 1) * proc_length + 1 : i  * proc_length) , proc_length * member_num , MPI_DOUBLE , proc%id + tmp_dims(2) * (i-1) , 100 , proc%comm , status,ierr)
                  end do

                  state%u(:,1:proc_length,j,k) = u_tran(: ,1:proc_length)
                else
                  call MPI_SEND(rhs(:,:) , size(rhs , 1) * size(rhs , 2) , MPI_DOUBLE , proc%id - tmp_dims(2) * tmp_coords(1) , 100 , proc%comm , ierr)
                  call MPI_RECV(state%u(:,mesh%half_lon_ibeg:mesh%half_lon_iend,j,k) , proc_length * member_num , MPI_DOUBLE , proc%id - tmp_dims(2) * tmp_coords(1) , 100 , proc%comm , status,ierr)
                end if
              end if
            else if (zonal_tridiag_solver == 'spk' ) then
              !do im = 1 , member_num
                call zonal_solver(j,k)%solve(rhs(:,:), state%u(:,mesh%half_lon_ibeg:mesh%half_lon_iend,j,k))
              !end do
            else
              if (is_root_proc()) call log_error('Wrong  zonal_tridiag_solver input ' // trim(zonal_tridiag_solver) // '!')
            end if
            !#####################################################
          else
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              state%u(:,i,j,k) = state%u(:,i,j,k) + dt * c_lon(j,k) * ( &
                state%div(:,i+1,j,k) - state%div(:,i,j,k)) / mesh%de_lon(j)
            end do
          end if
        end do
      end do

#ifdef Detail_Time
      call Get_End_Time(cal_time_end)
      cal_time = cal_time + cal_time_end - cal_time_start
      call Get_Start_Time(tran_time_start)
#endif

      call fill_halo_member(block, state%u, full_lon=.false., full_lat=.true., full_lev=.true., async=state%async(async_u))

#ifdef Detail_Time
      call Get_End_Time(tran_time_end)
      tran_time = tran_time + tran_time_end - tran_time_start
#endif

    case (4)
    end select

  end subroutine div_damp_run

end module div_damp_mod

