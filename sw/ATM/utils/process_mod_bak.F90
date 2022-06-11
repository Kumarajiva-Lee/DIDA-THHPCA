#define Only_Gmcore

module process_mod

  use mpi
  use flogger
  use string
  use namelist_mod
  use mesh_mod
  use block_mod

  implicit none

  private

  public process_init
  public process_create_blocks
  public process_stop
  public process_final
  public proc
  public is_root_proc
  public zonal_circle_type

  integer, public, parameter :: decomp_1d_lat = 1
  integer, public, parameter :: decomp_2d_simple = 2
  integer, public, parameter :: decomp_2d_irregular = 3

  integer, public, parameter :: decomp_reduce_south_region   = 1
  integer, public, parameter :: decomp_reduce_south_boundary = 2
  integer, public, parameter :: decomp_reduce_north_region   = 3
  integer, public, parameter :: decomp_reduce_north_boundary = 4
  integer, public, parameter :: decomp_normal_region         = 5
  integer, public, parameter :: decomp_normal_south_boundary = 6
  integer, public, parameter :: decomp_normal_north_boundary = 7

  type process_neighbor_type
    integer :: id       = MPI_PROC_NULL
    integer :: cart_id  = MPI_PROC_NULL
    integer :: orient   = 0
    integer :: lon_ibeg = inf_i4
    integer :: lon_iend = inf_i4
    integer :: lat_ibeg = inf_i4
    integer :: lat_iend = inf_i4
    integer :: lon_ibeg_r = inf_i4
    integer :: lon_ibeg_s = inf_i4
    integer :: lon_iend_r = inf_i4
    integer :: lon_iend_s = inf_i4
    logical :: unequal  = .false.
  contains
    procedure :: init => process_neighbor_init
  end type process_neighbor_type

  type zonal_circle_type
    integer :: group = MPI_GROUP_NULL
    integer :: comm  = MPI_COMM_NULL
    integer :: np    = 0
    integer :: id    = MPI_PROC_NULL
    integer, allocatable :: recv_type_r8(:,:) ! 0: one level, 1: full_lev, 2: half_lev
  contains
    procedure :: init => zonal_circle_init
    final :: zonal_circle_final
  end type zonal_circle_type

  type irr_comm_type
    integer :: group = MPI_GROUP_NULL
    integer :: comm  = MPI_COMM_NULL
    integer :: cart_comm = MPI_COMM_NULL
    integer :: np    = 0
    integer :: id    = MPI_PROC_NULL
    integer :: cart_id = MPI_PROC_NULL
    integer :: part_id
    integer :: cart_dims(2)   = 0
    integer :: cart_coords(2) = 0
    integer, allocatable :: recv_type_r8(:,:) ! 0: one level, 1: full_lev, 2: half_lev

  end type irr_comm_type

  type process_type
    integer :: comm           = MPI_COMM_NULL
    integer :: cart_comm      = MPI_COMM_NULL
    integer :: group          = MPI_GROUP_NULL
    integer :: cart_group     = MPI_GROUP_NULL
    integer :: cart_dims(2)   = 0
    integer :: cart_coords(2) = 0
    integer :: id             = MPI_PROC_NULL          ! MPI process ID
    integer :: cart_id        = MPI_PROC_NULL          ! MPI process ID in cart_comm
    integer idom                                       ! Nest domain index (root domain is 1)
    integer np
    integer num_lon
    integer num_lat
    integer lon_ibeg
    integer lon_iend
    integer lat_ibeg
    integer lat_iend

    ! integer :: irr_comm           = MPI_COMM_NULL
    ! integer :: irr_group          = MPI_GROUP_NULL
    ! integer :: irr_id             = MPI_PROC_NULL  
    ! integer :: irr_cart_dims(2)   = MPI_COMM_NULL
    ! integer :: irr_cart_coords(2) = MPI_COMM_NULL

    logical NeedReduce
    integer member_num

    integer ngb_num

    type(irr_comm_type) irr_comm
    type(zonal_circle_type) zonal_circle
    type(process_neighbor_type), allocatable :: ngb(:) ! Neighbor processes
    type(block_type), allocatable :: blocks(:)

    integer decomp_type
    integer decomp_loc
  end type process_type

  type(process_type) proc

contains

  subroutine process_init(comm)

    integer, intent(in), optional :: comm

    integer ierr

    if (present(comm)) then
      proc%comm = comm
    else
      call MPI_INIT(ierr)
      proc%comm = MPI_COMM_WORLD
    end if
    call MPI_COMM_GROUP(proc%comm, proc%group, ierr)
    call MPI_COMM_SIZE(proc%comm, proc%np, ierr)
    call MPI_COMM_RANK(proc%comm, proc%id, ierr)

    proc%NeedReduce = 0
    proc%member_num = member_num
    
    if (partition_type == 'regular') then
      call setup_mpi_simple()
      call decompose_domains()
      call setup_zonal_comm_for_reduce()
    else if (partition_type == 'irregular') then
      call setup_mpi_irregular()
      call decompose_domains_irregular()
      call setup_zonal_comm_for_reduce()
    else 
      if (is_root_proc()) call log_error('Wrong partition_type !')
    end if


  end subroutine process_init

  subroutine process_stop(code)

    integer, intent(in) :: code

    integer ierr

    call MPI_ABORT(proc%comm, code, ierr)

  end subroutine process_stop

  subroutine process_final()

    integer ierr

    if (allocated(proc%ngb   )) deallocate(proc%ngb   )
    if (allocated(proc%blocks)) deallocate(proc%blocks)
    if (proc%group       /= MPI_GROUP_NULL) call MPI_GROUP_FREE(proc%group      , ierr)
    if (proc%cart_group  /= MPI_GROUP_NULL) call MPI_GROUP_FREE(proc%cart_group , ierr)

#ifdef Only_Gmcore
    call MPI_FINALIZE(ierr)
#endif

  end subroutine process_final

  pure logical function is_root_proc()

    is_root_proc = proc%id == 0

  end function is_root_proc

  subroutine setup_mpi_simple()

    integer ierr, np, tmp_comm, i
    logical periods(2)

    proc%decomp_type = decomp_1d_lat
    proc%decomp_loc  = decomp_normal_region

    if (num_proc_lon(1) * num_proc_lat(1) == proc%np) then
      ! Check if process topology in namelist is compatible with MPI runtime.
      np = 0
      do i = 1, nest_max_dom
        np = np + num_proc_lon(i) * num_proc_lat(i)
      end do
      if (proc%np /= np .and. is_root_proc()) then
        call log_notice('Namelist num_proc_lon and num_proc_lat are not compatible with MPI runtime. Reset to MPI runtime.')
        num_proc_lat(1) = proc%np
      end if
      ! Set the process topology into proc object.
      np = 0
      do i = 1, nest_max_dom
        np = np + num_proc_lon(i) * num_proc_lat(i)
        if (proc%id + 1 <= np) then
          proc%cart_dims(1) = num_proc_lon(i)
          proc%cart_dims(2) = num_proc_lat(i)
          proc%idom = i
          exit
        end if
      end do
    else
      proc%cart_dims = [1, proc%np]
      proc%idom = 1
    end if
    periods = [.true.,.false.]

    ! Set MPI process topology.
    call MPI_COMM_SPLIT(proc%comm, proc%idom, proc%id, tmp_comm, ierr)
    call MPI_CART_CREATE(tmp_comm, 2, proc%cart_dims, periods, .true., proc%cart_comm, ierr)
    call MPI_COMM_GROUP(proc%cart_comm, proc%cart_group, ierr)
    call MPI_COMM_FREE(tmp_comm, ierr)
    call MPI_COMM_RANK(proc%cart_comm, proc%cart_id, ierr)
    call MPI_CART_COORDS(proc%cart_comm, proc%cart_id, 2, proc%cart_coords, ierr)

  end subroutine setup_mpi_simple

  subroutine decompose_domains()

    integer ierr, tmp_id(1), i, j
    integer hw

    ! Set neighborhood of the process.
    if (allocated(proc%ngb)) deallocate(proc%ngb)
    select case (proc%decomp_loc)
    case (decomp_normal_region)
      allocate(proc%ngb(4))
      call MPI_CART_SHIFT(proc%cart_comm, 0, 1, proc%ngb(west )%cart_id, proc%ngb(east )%cart_id, ierr)
      call MPI_CART_SHIFT(proc%cart_comm, 1, 1, proc%ngb(south)%cart_id, proc%ngb(north)%cart_id, ierr)
    end select

    ! Translate Cartesian ID of neighbors to global ID.
    do i = 1, size(proc%ngb)
      if (proc%ngb(i)%id == MPI_PROC_NULL) then
        call MPI_GROUP_TRANSLATE_RANKS(proc%cart_group, 1, [proc%ngb(i)%cart_id], proc%group, tmp_id, ierr)
        proc%ngb(i)%id = tmp_id(1)
      end if
    end do

    proc%ngb_num = 4

    ! Set initial values for num_lon, num_lat, lon_ibeg, lat_ibeg.
    proc%num_lon = global_mesh%num_full_lon
    select case (proc%decomp_loc)
    case (decomp_normal_region)
#ifdef V_POLE
      proc%num_lat = global_mesh%num_half_lat
#else
      proc%num_lat = global_mesh%num_full_lat
#endif
    end select

    call round_robin(proc%cart_dims(1), proc%cart_coords(1), proc%num_lon, proc%lon_ibeg, proc%lon_iend)
    call round_robin(proc%cart_dims(2), proc%cart_coords(2), proc%num_lat, proc%lat_ibeg, proc%lat_iend)

    hw = global_mesh%lon_halo_width

    ! TODO: Support 2D decomposition.
    call proc%ngb(west )%init(west , lat_ibeg=proc%lat_ibeg, lat_iend=proc%lat_iend)
    call proc%ngb(east )%init(east , lat_ibeg=proc%lat_ibeg, lat_iend=proc%lat_iend)
    call proc%ngb(south)%init(south, lon_ibeg=proc%lon_ibeg-hw, lon_iend=proc%lon_iend+hw)
    call proc%ngb(north)%init(north, lon_ibeg=proc%lon_ibeg-hw, lon_iend=proc%lon_iend+hw)

  end subroutine decompose_domains

  subroutine setup_mpi_irregular()

    integer ierr, np, tmp_comm, i , irr_iter
    integer tot_proc , tot_lat 
    integer, allocatable :: irr_group_rank(:)
    integer, allocatable :: irr_group_arr(:)
    integer, allocatable :: irr_comm_arr(:)
    logical periods(2)

    proc%decomp_type = decomp_2d_irregular
    proc%decomp_loc  = decomp_normal_region
    proc%idom = 1

    if (num_proc_lon(1) * num_proc_lat(1) /= proc%np) then
      if (is_root_proc()) call log_error('Wrong ATM total proc !')
    end if

    if (irr_part == 0) then
      if (is_root_proc()) call log_error('Missing irr_part !')
    end if

    tot_proc = 0
    tot_lat  = 0
    do i = 1 , irr_part
      tot_lat  = tot_lat + irr_num_lat(i)
      tot_proc = tot_proc + irr_num_proc_lon(i) * irr_num_proc_lat(i)
    end do

    if (tot_lat /= num_lat .or. tot_proc /= num_proc_lon(1) * num_proc_lat(1)) then
      if (is_root_proc()) call log_error('Wrong irregular input !')
    end if

    allocate(irr_group_arr(irr_part))
    allocate(irr_comm_arr(irr_part))

    tot_proc = 0
    
    do irr_iter = 1 , irr_part
      allocate(irr_group_rank(irr_num_proc_lon(irr_iter) * irr_num_proc_lat(irr_iter)))

      do i = 1 , irr_num_proc_lon(irr_iter) * irr_num_proc_lat(irr_iter)
        irr_group_rank(i) = tot_proc + i - 1
      end do

      call MPI_GROUP_INCL(proc%group , size(irr_group_rank), irr_group_rank, irr_group_arr(irr_iter) , ierr)
      call MPI_COMM_CREATE(proc%comm , irr_group_arr(irr_iter) , irr_comm_arr(irr_iter) , ierr)

      do i = 1 , irr_num_proc_lon(irr_iter) * irr_num_proc_lat(irr_iter)
        if (proc%id == tot_proc + i - 1) then
          proc%irr_comm%group   = irr_group_arr(irr_iter)
          proc%irr_comm%comm    = irr_comm_arr(irr_iter)
          proc%irr_comm%part_id = irr_iter
        end if
      end do

      tot_proc = tot_proc + irr_num_proc_lon(irr_iter) * irr_num_proc_lat(irr_iter)
      deallocate(irr_group_rank)
    end do

    call MPI_COMM_SIZE(proc%irr_comm%comm , proc%irr_comm%np , ierr)
    call MPI_COMM_RANK(proc%irr_comm%comm,  proc%irr_comm%id, ierr)

    deallocate(irr_group_arr)
    deallocate(irr_comm_arr)

    periods = [.true.,.false.]
    do irr_iter = 1 , irr_part
      if (proc%irr_comm%part_id == irr_iter) then
        proc%irr_comm%cart_dims(1)=irr_num_proc_lon(irr_iter)
        proc%irr_comm%cart_dims(2)=irr_num_proc_lat(irr_iter)

        call MPI_CART_CREATE(proc%irr_comm%comm, 2, proc%irr_comm%cart_dims, periods, .true., proc%irr_comm%cart_comm, ierr)
        call MPI_COMM_RANK(proc%irr_comm%cart_comm, proc%irr_comm%cart_id, ierr)
        call MPI_CART_COORDS(proc%irr_comm%cart_comm, proc%irr_comm%cart_id, 2, proc%irr_comm%cart_coords, ierr)
      end if
    end do

    !write(*,*) proc%id , proc%irr_comm%id , proc%irr_comm%np , proc%irr_comm%cart_id , proc%irr_comm%cart_coords(1) , proc%irr_comm%cart_coords(2)

  end subroutine setup_mpi_irregular

  
  subroutine decompose_domains_irregular()

    integer ierr, tmp_id(1), i, j, p
    integer hw
    integer tmp_num_lat(0:irr_part) , tmp_num_proc(0:irr_part)
    integer tmp_part, proc_beg, proc_end, tmp_num, tmp_beg, tmp_end, tmp_beg_r, tmp_end_r, tmp_beg_s, tmp_end_s 

    tmp_num_lat(0:1) = 0
    tmp_num_proc(0:1) = 0
    do i = 2 , irr_part
      tmp_num_lat(i) = sum(irr_num_lat(1:i-1))
      tmp_num_proc(i) = tmp_num_proc(i-1) + irr_num_proc_lon(i-1) * irr_num_proc_lat(i-1)
    end do


    ! Set neighborhood of the process.
    if (allocated(proc%ngb)) deallocate(proc%ngb)
    allocate(proc%ngb(10))
    proc%ngb_num = 4
    call MPI_CART_SHIFT(proc%irr_comm%cart_comm, 0, 1, proc%ngb(west )%cart_id, proc%ngb(east )%cart_id, ierr)
    call MPI_CART_SHIFT(proc%irr_comm%cart_comm, 1, 1, proc%ngb(south)%cart_id, proc%ngb(north)%cart_id, ierr)
    proc%ngb(west )%orient = west
    proc%ngb(east )%orient = east
    proc%ngb(south)%orient = south
    proc%ngb(north)%orient = north

    ! Translate Cartesian ID of neighbors to global ID.
    do i = 1, proc%ngb_num
      if (proc%ngb(i)%id == MPI_PROC_NULL) then
        call MPI_GROUP_TRANSLATE_RANKS(proc%irr_comm%group, 1, [proc%ngb(i)%cart_id], proc%group, tmp_id, ierr)
        proc%ngb(i)%id = tmp_id(1)
      end if
    end do

    !write(*,*) proc%id , proc%irr_comm%id , proc%irr_comm%np , proc%ngb(west )%id, proc%ngb(east )%id , proc%ngb(south)%id, proc%ngb(north)%id


    ! Set initial values for num_lon, num_lat, lon_ibeg, lat_ibeg.
    proc%num_lon = global_mesh%num_full_lon
    proc%num_lat = irr_num_lat(proc%irr_comm%part_id)

    call round_robin(proc%irr_comm%cart_dims(1), proc%irr_comm%cart_coords(1), proc%num_lon, proc%lon_ibeg, proc%lon_iend)
    call round_robin(proc%irr_comm%cart_dims(2), proc%irr_comm%cart_coords(2), proc%num_lat, proc%lat_ibeg, proc%lat_iend)
    proc%lat_ibeg = proc%lat_ibeg + tmp_num_lat(proc%irr_comm%part_id)
    proc%lat_iend = proc%lat_iend + tmp_num_lat(proc%irr_comm%part_id)

    hw = global_mesh%lon_halo_width

    ! Set muti ngb & init
    do i = 1, proc%ngb_num

      if (proc%ngb(i)%id == -1) then
        tmp_part = proc%irr_comm%part_id

        if (proc%ngb(i)%orient == south) then
          if (proc%lat_ibeg == 1) then 
            call proc%ngb(i)%init(south, lon_ibeg=proc%lon_ibeg-hw, lon_iend=proc%lon_iend+hw) 
          else
            proc_beg = tmp_num_proc(tmp_part - 1) + irr_num_proc_lat(tmp_part - 1) - 1
            proc_end = tmp_num_proc(tmp_part) - 1
            p = 0

            if (irr_num_proc_lon(tmp_part) < irr_num_proc_lon(tmp_part - 1)) then
              !bigger than ngb
              do j = proc_beg , proc_end , irr_num_proc_lat(tmp_part - 1)
                tmp_num = num_lon
                call round_robin(irr_num_proc_lon(tmp_part - 1), p, tmp_num , tmp_beg, tmp_end)
                if (proc%lon_ibeg <= tmp_beg .and.  tmp_end <= proc%lon_iend) then

                  
                  
                  if (proc%lon_ibeg == tmp_beg) then
                    tmp_beg_r = tmp_beg - hw
                  else 
                    tmp_beg_r = tmp_beg
                  end if

                  if (proc%lon_iend == tmp_end) then
                    tmp_end_r = tmp_end + hw
                  else
                    tmp_end_r = tmp_end
                  end if  

                  if (proc%ngb(i)%id == -1) then 
                    proc%ngb(i)%id = j 
                    call proc%ngb(i)%init(south, is_unequal = .true. , lon_ibeg_r = tmp_beg_r,    lon_iend_r = tmp_end_r, &
                                                                       lon_ibeg_s = tmp_beg - hw, lon_iend_s = tmp_end + hw)

                  else !muti ngb
                    proc%ngb_num = proc%ngb_num + 1
                    proc%ngb(proc%ngb_num)%id = j
                    call proc%ngb(proc%ngb_num)%init(south, is_unequal = .true. , lon_ibeg_r = tmp_beg_r,    lon_iend_r = tmp_end_r, &
                                                                                  lon_ibeg_s = tmp_beg - hw, lon_iend_s = tmp_end + hw)
                  end if
                end if

                p = p + 1
              end do
            else
              !smaller than ngb
              do j = proc_beg , proc_end , irr_num_proc_lat(tmp_part - 1)
                tmp_num = num_lon
                call round_robin(irr_num_proc_lon(tmp_part - 1), p, tmp_num , tmp_beg, tmp_end)
                if (tmp_beg <= proc%lon_ibeg .and. proc%lon_iend <= tmp_end) then

                  if(tmp_beg == proc%lon_ibeg) then
                    tmp_beg_s = proc%lon_ibeg - hw
                  else
                    tmp_beg_s = proc%lon_ibeg
                  end if

                  if (tmp_end == proc%lon_iend) then
                    tmp_end_s = proc%lon_iend + hw
                  else
                    tmp_end_s = proc%lon_iend
                  end if

                  proc%ngb(i)%id = j
                  call proc%ngb(i)%init(south, is_unequal = .true. , lon_ibeg_r = proc%lon_ibeg - hw, lon_iend_r = proc%lon_iend + hw, &
                                                                     lon_ibeg_s = tmp_beg_s         , lon_iend_s = tmp_end_s              )
                end if

                p = p + 1
              end do
            end if
          end if
        end if

        if (proc%ngb(i)%orient == north) then
          if (proc%lat_iend == num_lat) then 
            call proc%ngb(i)%init(north, lon_ibeg=proc%lon_ibeg-hw, lon_iend=proc%lon_iend+hw) 
          else
            proc_beg = tmp_num_proc(tmp_part + 1)
            proc_end = proc_beg + irr_num_proc_lon(tmp_part + 1) * irr_num_proc_lat(tmp_part + 1)
            p = 0

            if (irr_num_proc_lon(tmp_part) < irr_num_proc_lon(tmp_part + 1)) then
              !bigger than ngb
              do j = proc_beg , proc_end , irr_num_proc_lat(tmp_part + 1)
                tmp_num = num_lon
                call round_robin(irr_num_proc_lon(tmp_part + 1), p, tmp_num , tmp_beg, tmp_end)
                if (proc%lon_ibeg <= tmp_beg .and.  tmp_end <= proc%lon_iend) then

                  if (proc%lon_ibeg == tmp_beg) then
                    tmp_beg_r = tmp_beg - hw
                  else 
                    tmp_beg_r = tmp_beg
                  end if

                  if (proc%lon_iend == tmp_end) then
                    tmp_end_r = tmp_end + hw
                  else
                    tmp_end_r = tmp_end
                  end if

                  if (proc%ngb(i)%id == -1) then 
                    proc%ngb(i)%id = j 
                    call proc%ngb(i)%init(north, is_unequal = .true. , lon_ibeg_r = tmp_beg_r,    lon_iend_r = tmp_end_r, &
                                                                       lon_ibeg_s = tmp_beg - hw, lon_iend_s = tmp_end + hw)
                  else !muti ngb
                    proc%ngb_num = proc%ngb_num + 1
                    proc%ngb(proc%ngb_num)%id = j
                    call proc%ngb(proc%ngb_num)%init(north, is_unequal = .true. , lon_ibeg_r = tmp_beg_r,    lon_iend_r = tmp_end_r, &
                                                                                  lon_ibeg_s = tmp_beg - hw, lon_iend_s = tmp_end + hw)
                  end if
                end if
                
                p = p + 1
              end do
            else
              !smaller than ngb
              do j = proc_beg , proc_end , irr_num_proc_lat(tmp_part - 1)
                tmp_num = num_lon
                call round_robin(irr_num_proc_lon(tmp_part + 1), p, tmp_num , tmp_beg, tmp_end)
                if (tmp_beg <= proc%lon_ibeg .and. proc%lon_iend <= tmp_end) then

                  if(tmp_beg == proc%lon_ibeg) then
                    tmp_beg_s = proc%lon_ibeg - hw
                  else
                    tmp_beg_s = proc%lon_ibeg
                  end if

                  if (tmp_end == proc%lon_iend) then
                    tmp_end_s = proc%lon_iend + hw
                  else
                    tmp_end_s = proc%lon_iend
                  end if

                  proc%ngb(i)%id = j
                  call proc%ngb(i)%init(north, is_unequal = .true. , lon_ibeg_r = proc%lon_ibeg - hw, lon_iend_r = proc%lon_iend + hw, &
                                                                     lon_ibeg_s = tmp_beg_s         , lon_iend_s = tmp_end_s             )
                end if

                p = p + 1
              end do
            end if
          end if
        end if
      else
        if (proc%ngb(i)%orient == south .or. proc%ngb(i)%orient == north) then
          call proc%ngb(i)%init(proc%ngb(i)%orient, lon_ibeg=proc%lon_ibeg-hw, lon_iend=proc%lon_iend+hw)
        else
          call proc%ngb(i)%init(proc%ngb(i)%orient, lat_ibeg=proc%lat_ibeg, lat_iend=proc%lat_iend)
        end if
      end if  
    end do

    ! write(*,*) proc%id , proc%irr_comm%id , proc%irr_comm%np , proc%lon_ibeg, proc%lon_iend , proc%lat_ibeg, proc%lat_iend



    ! hw = global_mesh%lon_halo_width

    ! ! TODO: Support 2D decomposition.
    ! call proc%ngb(west )%init(west , lat_ibeg=proc%lat_ibeg, lat_iend=proc%lat_iend)
    ! call proc%ngb(east )%init(east , lat_ibeg=proc%lat_ibeg, lat_iend=proc%lat_iend)
    ! call proc%ngb(south)%init(south, lon_ibeg=proc%lon_ibeg-hw, lon_iend=proc%lon_iend+hw)
    ! call proc%ngb(north)%init(north, lon_ibeg=proc%lon_ibeg-hw, lon_iend=proc%lon_iend+hw)

  end subroutine decompose_domains_irregular

  subroutine setup_zonal_comm_for_reduce()

    ! Create zonal communicator for reduce algorithm.
    if (proc%idom == 1) then ! Only root domain has reduce region.
      call proc%zonal_circle%init()
    end if

  end subroutine setup_zonal_comm_for_reduce

  subroutine process_create_blocks(para_member_num)
    
    integer, intent(in) :: para_member_num
    integer i, j, dtype

    if (.not. allocated(proc%blocks)) allocate(proc%blocks(1))

    call proc%blocks(1)%init(proc%id, para_member_num , global_mesh%lon_halo_width, global_mesh%lat_halo_width, &
                             proc%lon_ibeg, proc%lon_iend, proc%lat_ibeg, proc%lat_iend)

    

    select case (r8)
    case (4)
      dtype = MPI_REAL
    case (8)
      dtype = MPI_DOUBLE
    case (16)
      dtype = MPI_REAL16
    case default
      call log_error('Unsupported parameter r8!')
    end select

    ! Setup halos (only normal halos for the time being).
    allocate(proc%blocks(1)%halo(proc%ngb_num))
    do i = 1, proc%ngb_num
      select case (proc%ngb(i)%orient)
      case (west, east)
        call proc%blocks(1)%halo(i)%init(proc%blocks(1)%mesh, proc%ngb(i)%orient, dtype,               &
                                         host_id=proc%id, ngb_proc_id=proc%ngb(i)%id, member_num = proc%member_num , &
                                         lat_ibeg=proc%ngb(i)%lat_ibeg, lat_iend=proc%ngb(i)%lat_iend)
      case (south, north)
        if (.not. proc%ngb(i)%unequal) then
          call proc%blocks(1)%halo(i)%init(proc%blocks(1)%mesh, proc%ngb(i)%orient, dtype,               &
                                          host_id=proc%id, ngb_proc_id=proc%ngb(i)%id, member_num = proc%member_num , &
                                          lon_ibeg=proc%ngb(i)%lon_ibeg, lon_iend=proc%ngb(i)%lon_iend)
        else
          call proc%blocks(1)%halo(i)%init(proc%blocks(1)%mesh, proc%ngb(i)%orient, dtype,               &
                                          host_id=proc%id, ngb_proc_id=proc%ngb(i)%id, member_num = proc%member_num , &
                                          is_unequal = .true. , lon_ibeg_r=proc%ngb(i)%lon_ibeg_r, lon_iend_r=proc%ngb(i)%lon_iend_r , &
                                                                lon_ibeg_s=proc%ngb(i)%lon_ibeg_s, lon_iend_s=proc%ngb(i)%lon_iend_s)
        end if
      end select
    end do
    ! Initialize async objects.
    do i = 1, size(proc%blocks(1)%state)
      do j = 1, size(proc%blocks(1)%state(i)%async)
        call proc%blocks(1)%state(i)%async(j)%init(proc%ngb_num)
      end do
    end do

  end subroutine process_create_blocks

    subroutine process_neighbor_init(this, orient, lon_ibeg, lon_iend, lat_ibeg, lat_iend, is_unequal, lon_ibeg_r, lon_iend_r, lon_ibeg_s, lon_iend_s)

    class(process_neighbor_type), intent(inout) :: this
    integer, intent(in) :: orient
    integer, intent(in), optional :: lon_ibeg
    integer, intent(in), optional :: lon_iend
    integer, intent(in), optional :: lat_ibeg
    integer, intent(in), optional :: lat_iend
    logical, intent(in), optional :: is_unequal
    integer, intent(in), optional :: lon_ibeg_r
    integer, intent(in), optional :: lon_iend_r
    integer, intent(in), optional :: lon_ibeg_s
    integer, intent(in), optional :: lon_iend_s


    this%orient = orient

    if (present(is_unequal)) this%unequal = is_unequal

    if (.not. this%unequal) then
      select case (orient)
      case (west, east)
        this%lat_ibeg = lat_ibeg
        this%lat_iend = lat_iend
      case (south, north)
        this%lon_ibeg   = lon_ibeg
        this%lon_iend   = lon_iend
        this%lon_ibeg_r = this%lon_ibeg
        this%lon_ibeg_s = this%lon_ibeg
        this%lon_iend_r = this%lon_iend
        this%lon_iend_s = this%lon_iend
      end select
    else
      select case (orient)
      case (west, east)
        this%lat_ibeg = lat_ibeg
        this%lat_iend = lat_iend
      case (south, north)
        this%lon_ibeg_r = lon_ibeg_r
        this%lon_iend_r = lon_iend_r
        this%lon_ibeg_s = lon_ibeg_s
        this%lon_iend_s = lon_iend_s
      end select
    end if


  end subroutine process_neighbor_init

  subroutine zonal_circle_init(this)

    class(zonal_circle_type), intent(inout) :: this

    integer ierr, i, num_lon, ibeg, iend
    integer, allocatable :: zonal_proc_id(:)
    
    if (partition_type == 'regular') then
      allocate(zonal_proc_id(proc%cart_dims(1)))
      do i = 1, proc%cart_dims(1)
        call MPI_CART_RANK(proc%cart_comm, [i-1,proc%cart_coords(2)], zonal_proc_id(i), ierr)
      end do
      call MPI_GROUP_INCL(proc%cart_group, size(zonal_proc_id), zonal_proc_id, this%group, ierr)
      call MPI_COMM_CREATE_GROUP(proc%cart_comm, this%group, sum(zonal_proc_id), this%comm, ierr)
      call MPI_COMM_SIZE(this%comm, this%np, ierr)
      call MPI_COMM_RANK(this%comm, this%id, ierr)
      deallocate(zonal_proc_id)
    else if (partition_type == 'irregular') then
      allocate(zonal_proc_id(proc%irr_comm%cart_dims(1)))
      do i = 1, proc%irr_comm%cart_dims(1)
        call MPI_CART_RANK(proc%irr_comm%cart_comm, [i-1,proc%irr_comm%cart_coords(2)], zonal_proc_id(i), ierr)
      end do
      call MPI_GROUP_INCL(proc%irr_comm%group, size(zonal_proc_id), zonal_proc_id, this%group, ierr)
      call MPI_COMM_CREATE_GROUP(proc%irr_comm%cart_comm, this%group, sum(zonal_proc_id), this%comm, ierr)
      call MPI_COMM_SIZE(this%comm, this%np, ierr)
      call MPI_COMM_RANK(this%comm, this%id, ierr)
      deallocate(zonal_proc_id)
    end if

    !write(*,*) proc%id , this%id , this%np 

    if (this%id == 0) then
      allocate(this%recv_type_r8(this%np,0:2))
      do i = 1, this%np
        num_lon = global_mesh%num_full_lon
        call round_robin(this%np, i - 1, num_lon, ibeg, iend)
        call MPI_TYPE_CREATE_SUBARRAY(2, [member_num , global_mesh%num_full_lon], &
                                         [member_num ,             num_lon], &
                                         [0,ibeg-1], MPI_ORDER_FORTRAN, MPI_DOUBLE, &
                                         this%recv_type_r8(i,0), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_r8(i,0), ierr)
        call MPI_TYPE_CREATE_SUBARRAY(3, [member_num , global_mesh%num_full_lon,global_mesh%num_full_lev], &
                                         [member_num ,                  num_lon,global_mesh%num_full_lev], &
                                         [0,ibeg-1,0], MPI_ORDER_FORTRAN, MPI_DOUBLE, &
                                         this%recv_type_r8(i,1), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_r8(i,1), ierr)
        call MPI_TYPE_CREATE_SUBARRAY(3, [member_num , global_mesh%num_full_lon,global_mesh%num_half_lev], &
                                         [member_num ,                  num_lon,global_mesh%num_half_lev], &
                                         [0,ibeg-1,0], MPI_ORDER_FORTRAN, MPI_DOUBLE, &
                                         this%recv_type_r8(i,2), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_r8(i,2), ierr)
      end do
    end if

  end subroutine zonal_circle_init

  subroutine zonal_circle_final(this)

    type(zonal_circle_type), intent(inout) :: this

    integer i, k, ierr

    if (allocated(this%recv_type_r8)) then
      do k = 0, 2
        do i = 1, this%np
          call MPI_TYPE_FREE(this%recv_type_r8(i,k), ierr)
        end do
        deallocate(this%recv_type_r8)
      end do
    end if

    if (this%group /= MPI_GROUP_NULL) call MPI_GROUP_FREE(this%group, ierr)

  end subroutine zonal_circle_final

  subroutine round_robin(dim, coord, num, ibeg, iend)

    integer, intent(in) :: dim
    integer, intent(in) :: coord
    integer, intent(inout) :: num
    integer, intent(out) :: ibeg ! Start from 1.
    integer, intent(out) :: iend ! Start from 1.

    integer res_num, tmp_num, i

    res_num = mod(num, dim)
    ibeg = 1
    do i = 0, coord - 1
      if (res_num /= 0 .and. i < res_num) then
        tmp_num = num / dim + 1
      else
        tmp_num = num / dim
      end if
      ibeg = ibeg + tmp_num
    end do
    if (res_num /= 0 .and. coord < res_num) then
      num = num / dim + 1
    else
      num = num / dim
    end if
    iend = ibeg + num - 1

  end subroutine round_robin

end module process_mod
