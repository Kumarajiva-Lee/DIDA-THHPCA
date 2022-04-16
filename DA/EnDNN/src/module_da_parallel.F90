#include "../../utils/da_config.inc"

!#define HALF_PRECISION 0

module da_parallel
  use mpi
  !use coupler_mod
  use letkf_mod, only:letkf_init, letkf_run
#ifdef HALF_PRECISION
  use rp_emulator
#endif
  !use toy_model
  use pro_info
  use string
  use splitchar
  use var_info_mod
  !use obsmaker_mod
  use da_vert_coord_mod
  use coupler_config, ens => cnn_ensemble
  use da_namelist_mod
  ! use coupler_config, only:vert_coord_scheme, vert_coord_template
  ! use da_namelist_mod, only:obs_hrz_interval_lon, &
  ! obs_hrz_interval_lat, obs_vtc_interval, &
  ! da_polar_letkf_lat, da_inflation_factor, parse_namelist, &
  ! da_halo, da_loc_distance, da_gc_enable, da_obs_dense, &
  ! da_block_interval, da_block_overlap

  implicit none
  integer :: rank, lrank, lsize
  integer :: comm_cart
  integer :: dims(2), coords(2)
  logical :: periods(2)

  integer :: x_sec, y_sec, x_res, y_res
  integer :: xst, xed, yst, yed
  integer :: x_size, y_size
  integer :: my_test
  
  integer :: n_max
  integer, allocatable, dimension(:) :: obs_var_lev
  integer, allocatable, dimension(:) :: obs_var_start
  
  xb_type, allocatable, dimension(:,:,:,:) :: xb
  axb_type, allocatable, dimension(:,:,:) :: axb
  xb_type, allocatable, dimension(:,:,:,:) :: xa
  axb_type, allocatable, dimension(:,:,:) :: axa

  type(c_ptr) :: da_rc

  contains

  ! get kernel & stride number, calculate layers
  ! mpi_comm_split, mct_init, change the numbers of mct modules
  subroutine init_da_parallel(namelist_path, MPI_DA_GROUP) 
    character(*), intent(in) :: namelist_path
    integer, intent(in) :: MPI_DA_GROUP
    integer :: ierrs
    integer :: nx, ny
    integer :: i
    integer :: n
    integer :: n_var
    integer, allocatable, dimension(:) :: var_lev
    character(len = 20), dimension(10) :: var_name
    character(len = 20), allocatable, dimension(:) :: var_name_ordered
    integer :: obs_num_2d, obs_num_3d
  
    dims(1:2) = 0

    !wait snode info by ssp 
    !call wait_snode_info('DA')

    call parse_namelist(namelist_path)
    call MPI_Comm_rank(MPI_DA_GROUP, rank, ierrs)
    call MPI_Comm_size(MPI_DA_GROUP, lsize, ierrs)
    call MPI_Dims_create(lsize, 2, dims,ierrs)
    call MPI_Cart_create(MPI_DA_GROUP, 2, dims, periods, 0, comm_cart, ierrs)
    call MPI_Comm_rank(comm_cart, lrank, ierrs)
    call MPI_Cart_coords(comm_cart, lrank, 2, coords, ierrs)

    da_rc=RedisConnect(trim(redis_address))
    !da_rc=RedisConnect_Balanced(redis_address)


    ! print *, "da info:", lsize, rank, dims
    x_sec = num_lon / dims(1)
    x_res = mod(num_lon, dims(1))
    y_sec = num_lat / dims(2)
    y_res = mod(num_lat, dims(2))
    x_size = xed - xst
    y_size = yed - yst

    if(coords(1) < x_res) then
        xst = coords(1) * (x_sec + 1) + 1
        xed = xst + x_sec
    else
        xst = coords(1) * x_sec + x_res + 1
        xed = xst + x_sec - 1
    end if

    if(coords(2) < y_res) then
        yst = coords(2) * (y_sec + 1) + 1
        yed = yst + y_sec 
    else
        yst = coords(2) * y_sec + y_res + 1
        yed = yst + y_sec - 1
    end if

    !print *, "da:", lrank, xst, xed, yst, yed
    nx = xed-xst+1
    ny = yed-yst+1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !!!!!!
    ! vars parsing
    n_max = 0
    call stringsplit(da_var_name,',',var_name,n_var)
    allocate(var_name_ordered(n_var))
    allocate(var_lev(n_var))
    call var_info(var_name, n_var, var_name_ordered, obs_num_2d, obs_num_3d)
    n_max = 0
    do i = 1, obs_num_2d
      var_lev(i) = 1
      n_max = n_max + var_lev(i)
    end do
    do i = obs_num_2d + 1, n_var
      var_lev(i) = num_lev
      n_max = n_max + var_lev(i)
    end do
    deallocate(var_lev)
    deallocate(var_name_ordered)
    ! !!!!!!!!!!!!!
    ! !!!!!!!!!!!!!
    ! call obs_init(ens)
    !print *, "pro_info_init"
    call pro_info_init(da_rc, ny, nx, coords(2)+1, coords(1)+1, y_sec, y_res, x_sec, x_res)
    ! call pro_info_show()
    !print *, "pro_info_init end"

    call da_vert_coord_init(num_lev, vert_coord_scheme, vert_coord_template)

    
    ! call letkf_init(ens, ny, nx, num_lev, da_halo, n_max, &
    ! obs_hrz_interval_lat, obs_hrz_interval_lon, obs_vtc_interval, &
    ! da_polar_letkf_lat, da_inflation_factor, &
    ! rank, coords(2)+1, coords(1)+1, dims(2), dims(1), &
    ! da_block_interval, da_block_overlap, da_gc_enable, .true., &
    ! da_loc_distance, da_obs_dense, pos2id, MPI_DA_GROUP, da_rc)
    !print *, "letkf init start"
    call letkf_init(ny, nx, n_max, rank, coords(2)+1, coords(1)+1, dims(2), dims(1), &
    MPI_DA_GROUP, da_rc, pos2id)
 
    !print *,'lekf init sucess'   
    allocate(xb(ens, n*num_lev, 1-da_halo:ny+da_halo, 1-da_halo:nx+da_halo))
    allocate(axb(n*num_lev, 1-da_halo:ny+da_halo, 1-da_halo:nx+da_halo))
    allocate(xa(ens, n*num_lev, 1-da_halo:ny+da_halo, 1-da_halo:nx+da_halo))
    allocate(axa(n*num_lev, 1-da_halo:ny+da_halo, 1-da_halo:nx+da_halo))

    !call pos2id(0, 0, my_test)

    ! call toy_model_init(ens, n, ny, nx, num_lev, halo, &
    ! yst, xst, num_lat, num_lon, xb, axb)
    ! call pos2id(1, -1, my_test)
    ! print *, 'Before run my_test', rank, xst, xed, yst, yed

    !call letkf_run(xb, axb, xa, axa)

  end subroutine init_da_parallel

  subroutine run_da_parallel(n_max, nx, ny, num_times, xb_in, xa_out, spread_b, spread_a)

    integer, intent(in) :: n_max
    integer, intent(in) :: nx, ny
    integer, intent(in) :: num_times
    axb_type, intent(in) :: xb_in(ens, n_max, 1-da_halo:nx+da_halo, 1-da_halo:ny+da_halo)
    axb_type, intent(out) :: xa_out(ens, n_max, 1-da_halo:nx+da_halo, 1-da_halo:ny+da_halo)
    real(8), intent(out) :: spread_b(n_max, 1:nx, 1:ny)
    real(8), intent(out) :: spread_a(n_max, 1:nx, 1:ny)

#if (HALF_PRECISION == 1)
    type(rpe_var) :: tmp
#endif  
    integer :: i, j, k, ens_num

    xb_type :: xb(ens, n_max, 1-da_halo:nx+da_halo, 1-da_halo:ny+da_halo)
    axb_type:: axb(n_max, 1-da_halo:nx+da_halo, 1-da_halo:ny+da_halo)
    xb_type :: xa(ens, n_max, 1-da_halo:nx+da_halo, 1-da_halo:ny+da_halo)
    axb_type:: axa(n_max, 1-da_halo:nx+da_halo, 1-da_halo:ny+da_halo)

    ! dnn2letkf
    do i = 1-da_halo, ny+da_halo
      do j = 1-da_halo, nx+da_halo
        do k = 1, n_max
          do ens_num = 1, ens
            if (isnan(xb_in(ens_num, k, j, i))) then
              print *, rank, "bg nan!!!", k, j, i
              stop
            endif
          enddo
          axb(k, j, i) = sum(xb_in(1:ens, k, j, i)) / real(ens)
          xb(1:ens, k, j, i) = xb_in(1:ens, k, j, i) - axb(k, j, i)
#if (HALF_PRECISION == 1)
          do ens_num = 1, ens
            tmp%sbits = 10
            tmp%val = xb(ens_num, k, j, i)
            call apply_truncation(tmp)
            xb(ens_num, k, j, i) = tmp%val
          end do
#endif  
!          if (rank == 0 .and. i == 10 .and. j == 10 .and. k == 10) then
!            print *, i, j, xb(:, k, j, i)
!          end if
        end do

      end do
    end do
    
    call letkf_run(num_times, xb, axb, xa, axa)
    ! xa(:,:,:,:)=xb(:,:,:,:)
    ! axa(:,:,:)=axb(:,:,:)
    ! print*, "Only transform"

    ! letkf2dnn
    do i = 1-da_halo, ny+da_halo
      do j = 1-da_halo, nx+da_halo
        do k = 1, n_max
          if (i < 1 .or. i > ny .or. j < 1 .or. j > nx) then
            xa_out(1:ens, k, j, i) = 0
          else if (abs(grid_lat(j)) <= da_polar_letkf_lat) then
!            xa_out(1:ens, k, j, i) = xa(1:ens, k, j, i) + axa(k, j, i)
            do ens_num = 1, ens
!              xa_out(ens_num, k, j, i) = axb(k, j, i) + gamma1*xb(ens_num, k, j, i) + &
!                                         gamma2*(xa(ens_num, k, j, i) + axa(k, j, i) -&
!                                         axb(k, j, i))
              ! xa_out(ens_num, k, j, i) = axa(k, j, i) + gamma1*xb(ens_num, k, j, i) &
              !                          + gamma2*xa(ens_num, k, j, i) ! Kotsuki et al.(2017) fixed RTPP
              xa_out(ens_num, k, j, i) = gamma1*(axb(k, j, i) + xb(ens_num, k, j, i)) &
                                       + gamma2*(axa(k, j, i) + xa(ens_num, k, j, i)) ! hzd RTPP
            enddo

            spread_b(k, j, i) = 0
            spread_a(k, j, i) = 0
            do ens_num = 1, ens
              spread_b(k, j, i) = spread_b(k, j, i) + xb(ens_num, k, j, i)&
                                  * xb(ens_num, k, j, i)
              spread_a(k, j, i) = spread_a(k, j, i) + (xa_out(ens_num, k, j, i) - axa(k, j, i))&
                                  * (xa_out(ens_num, k, j, i) - axa(k, j, i))
            enddo
            spread_b(k, j, i) = sqrt(spread_b(k, j, i)/real(ens))
            spread_a(k, j, i) = sqrt(spread_a(k, j, i)/real(ens))
          else
            xa_out(1:ens, k, j, i) = xb(1:ens, k, j, i) + axb(k, j, i)
            spread_b(k, j, i) = 0
            spread_a(k, j, i) = 0
            do ens_num = 1, ens
              spread_b(k, j, i) = spread_b(k, j, i) + xb(ens_num, k, j, i)&
                                  * xb(ens_num, k, j, i)
              spread_a(k, j, i) = spread_a(k, j, i) + xb(ens_num, k, j, i)&
                                  * xb(ens_num, k, j, i)
            enddo
            spread_b(k, j, i) = sqrt(spread_b(k, j, i)/real(ens))
            spread_a(k, j, i) = sqrt(spread_a(k, j, i)/real(ens))
          end if
        end do
      end do
    end do

  end subroutine
!   subroutine run_da_parallel(n_max, nx, ny, num_times, xb_in, xa_out)

!     integer, intent(in) :: n_max
!     integer, intent(in) :: nx, ny
!     integer, intent(in) :: num_times
!     axb_type, intent(in) :: xb_in(ens, n_max, 1-da_halo:nx+da_halo, 1-da_halo:ny+da_halo)
!     axb_type, intent(out) :: xa_out(ens, n_max, 1-da_halo:nx+da_halo, 1-da_halo:ny+da_halo)

! #if (HALF_PRECISION == 1)
!     type(rpe_var) :: tmp
! #endif  
!     integer :: i, j, k, ens_num

!     xb_type :: xb(ens, n_max, 1-da_halo:nx+da_halo, 1-da_halo:ny+da_halo)
!     axb_type:: axb(n_max, 1-da_halo:nx+da_halo, 1-da_halo:ny+da_halo)
!     xb_type :: xa(ens, n_max, 1-da_halo:nx+da_halo, 1-da_halo:ny+da_halo)
!     axb_type:: axa(n_max, 1-da_halo:nx+da_halo, 1-da_halo:ny+da_halo)

!     ! dnn2letkf
!     do i = 1-da_halo, ny+da_halo
!       do j = 1-da_halo, nx+da_halo
!         do k = 1, n_max
!           do ens_num = 1, ens
!             if (isnan(xb_in(ens_num, k, j, i))) then
!               print *, rank, "bg nan!!!", k, j, i
!               stop
!             endif
!           enddo
!           axb(k, j, i) = sum(xb_in(1:ens, k, j, i)) / real(ens)
!           xb(1:ens, k, j, i) = xb_in(1:ens, k, j, i) - axb(k, j, i)
! #if (HALF_PRECISION == 1)
!           do ens_num = 1, ens
!             tmp%sbits = 10
!             tmp%val = xb(ens_num, k, j, i)
!             call apply_truncation(tmp)
!             xb(ens_num, k, j, i) = tmp%val
!           end do
! #endif  
! !          if (rank == 0 .and. i == 10 .and. j == 10 .and. k == 10) then
! !            print *, i, j, xb(:, k, j, i)
! !          end if
!         end do

!       end do
!     end do
    
!     call letkf_run(num_times, xb, axb, xa, axa)
!     ! xa(:,:,:,:)=xb(:,:,:,:)
!     ! axa(:,:,:)=axb(:,:,:)
!     ! print*, "Only transform"

!     ! letkf2dnn
!     do i = 1-da_halo, ny+da_halo
!       do j = 1-da_halo, nx+da_halo
!         do k = 1, n_max
!           if (i < 1 .or. i > ny .or. j < 1 .or. j > nx) then
!             xa_out(1:ens, k, j, i) = 0
!           else if (abs(grid_lat(j)) <= da_polar_letkf_lat) then
!             xa_out(1:ens, k, j, i) = xa(1:ens, k, j, i) + axa(k, j, i)
!           else
!             xa_out(1:ens, k, j, i) = xb(1:ens, k, j, i) + axb(k, j, i)
!           end if
!         end do
!         ! do k = 1, ens
!         ! end do
!         ! if (rank == 0 .and. mod(i, 10) == 0 .and. mod(j, 10) == 0) then
!         !   print *, i, j, xa(1, 1, j, i)
!         ! end if
!         ! 666 format(10.2F)
!       end do
!     end do

!   end subroutine

  ! function pos2id
#include "letkf_utils.F90"

end module da_parallel

#undef HALF_PRECISION
