!#define Only_Gmcore

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
  public process_init_zonal
  public process_create_blocks_stage1
  public process_create_blocks_stage2
  public process_stop
  public process_final
  public proc
  public is_root_proc
  public zonal_circle_type

  integer, public, parameter :: decomp_1d_lat = 1
  integer, public, parameter :: decomp_2d_simple = 2
  integer, public, parameter :: decomp_2d_irregular = 3
  integer, public, parameter :: decomp_normal_region  = 5

  type process_neighbor_type
    integer :: id       = MPI_PROC_NULL
    integer :: cart_id  = MPI_PROC_NULL
    integer :: orient   = 0                         !east west south north
    integer :: ngb_type = 0                         !send recv
    integer :: tag      = 0                         !message tag 
    integer :: lon_ibeg = inf_i4
    integer :: lon_iend = inf_i4
    integer :: lat_ibeg = inf_i4
    integer :: lat_iend = inf_i4
  contains
    procedure :: init => process_neighbor_init
  end type process_neighbor_type

  type zonal_circle_type
    integer :: group = MPI_GROUP_NULL
    integer :: comm  = MPI_COMM_NULL
    integer :: np    = 0
    integer :: id    = MPI_PROC_NULL
    integer :: west_ngb_id    = MPI_PROC_NULL
    integer :: east_ngb_id    = MPI_PROC_NULL
    integer, allocatable :: recv_type_r4(:,:) ! 0: one level, 1: full_lev, 2: half_lev
    integer, allocatable :: recv_type_r8(:,:) ! 0: one level, 1: full_lev, 2: half_lev
  contains
    procedure :: init => zonal_circle_init
    final :: zonal_circle_final
  end type zonal_circle_type

  type irr_comm_type
    integer :: group = MPI_GROUP_NULL
    integer :: comm  = MPI_COMM_NULL
    integer :: cart_group = MPI_GROUP_NULL
    integer :: cart_comm = MPI_COMM_NULL
    integer :: np    = 0
    integer :: id    = MPI_PROC_NULL
    integer :: cart_id = MPI_PROC_NULL
    integer :: part_id
    integer :: pid_ibeg = -1
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
    integer west_pids(2), east_pids(2)                 ! 给reduced_chunk用
    integer np
    integer num_lon
    integer num_lat
    integer lon_ibeg
    integer lon_iend
    integer lat_ibeg
    integer lat_iend
    integer comp_lat_ibeg
    integer comp_lat_iend
    integer lon_halo_width , lon_halo_width_north , lon_halo_width_south
    integer lat_halo_width , lat_halo_width_north , lat_halo_width_south

    integer member_value
    
    integer ngb_num
    integer ngb_send_num , ngb_recv_num
    integer ngb_north_send_beg , ngb_south_send_beg
    integer ngb_north_recv_beg , ngb_south_recv_beg

    type(irr_comm_type) irr_comm
    type(zonal_circle_type) zonal_circle
    type(process_neighbor_type), allocatable :: ngb(:) ! Neighbor processes send
    type(process_neighbor_type), allocatable :: ngb_send(:) ! Neighbor processes send
    type(process_neighbor_type), allocatable :: ngb_recv(:) ! Neighbor processes recv
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

    proc%member_value = member_num
    
    call setup_mpi_all_ngb()
    call decompose_domains_all_ngb()
    !call setup_zonal_comm_for_reduce()
 
  end subroutine process_init

  subroutine process_init_zonal(comm)

    integer, intent(in), optional :: comm

    call setup_zonal_comm_for_reduce()

  end subroutine process_init_zonal

  subroutine process_stop(code)

    integer, intent(in) :: code

    integer ierr

    call MPI_ABORT(proc%comm, code, ierr)

  end subroutine process_stop

  subroutine process_final()

    integer i, ierr

    do i = 1, size(proc%blocks)
      call proc%blocks(i)%clear()
    end do

    
    if (allocated(proc%blocks)) deallocate(proc%blocks)
    if (allocated(proc%ngb   )) deallocate(proc%ngb   )
    if (proc%group       /= MPI_GROUP_NULL) call MPI_GROUP_FREE(proc%group      , ierr)
    if (proc%cart_group  /= MPI_GROUP_NULL) call MPI_GROUP_FREE(proc%cart_group , ierr)

#ifdef Only_Gmcore
    call MPI_FINALIZE(ierr)
#endif

  end subroutine process_final

  pure logical function is_root_proc()

    is_root_proc = proc%id == 0

  end function is_root_proc

  subroutine setup_mpi_all_ngb()

    integer ierr, np, tmp_comm, i , irr_iter, near
    integer tot_proc , tot_lat 
    integer, allocatable :: irr_group_rank(:)
    integer, allocatable :: irr_group_arr(:)
    integer, allocatable :: irr_comm_arr(:)
    logical periods(2)

    proc%decomp_type = decomp_2d_irregular
    proc%decomp_loc  = decomp_normal_region
    proc%idom = 1

    if (num_proc_total /= proc%np) then
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

    if (tot_lat /= num_lat .or. tot_proc /= num_proc_total) then
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
          proc%irr_comm%pid_ibeg = tot_proc ! 这个irr_comm内起始的进程号
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
        call MPI_COMM_GROUP(proc%irr_comm%cart_comm,proc%irr_comm%cart_group , ierr)
        call MPI_COMM_RANK(proc%irr_comm%cart_comm, proc%irr_comm%cart_id, ierr)
        call MPI_CART_COORDS(proc%irr_comm%cart_comm, proc%irr_comm%cart_id, 2, proc%irr_comm%cart_coords, ierr)
      end if
    end do

    

    do near = 1, size(proc%west_pids)
      proc%west_pids(near) = proc%id - near * irr_num_proc_lat(proc%irr_comm%part_id)
      if (proc%west_pids(near) < proc%irr_comm%pid_ibeg) then
        proc%west_pids(near) = proc%west_pids(near) + irr_num_proc_lat(proc%irr_comm%part_id) * irr_num_proc_lon(proc%irr_comm%part_id)
      end if

      proc%east_pids(near) = proc%id + near * irr_num_proc_lat(proc%irr_comm%part_id)
      if (proc%east_pids(near) >= proc%irr_comm%pid_ibeg + proc%irr_comm%np) then
        proc%east_pids(near) = proc%east_pids(near) - irr_num_proc_lat(proc%irr_comm%part_id) * irr_num_proc_lon(proc%irr_comm%part_id)
      end if
    end do



  end subroutine setup_mpi_all_ngb

  subroutine decompose_domains_all_ngb()

    integer ierr, tmp_id(1), i, j, p, near
    integer hw
    integer tmp_num_lat(0:irr_part) , tmp_num_proc(0:irr_part)
    integer tmp_part, proc_beg, proc_end, tmp_num, tmp_lon_beg, tmp_lon_end , tmp_lat_beg, tmp_lat_end
    integer tmp_beg_r, tmp_end_r, tmp_beg_s, tmp_end_s 

    tmp_num_lat(0:1) = 0
    tmp_num_proc(0:1) = 0
    do i = 2 , irr_part
      tmp_num_lat(i) = sum(irr_num_lat(1:i-1))
      tmp_num_proc(i) = tmp_num_proc(i-1) + irr_num_proc_lon(i-1) * irr_num_proc_lat(i-1)
    end do


    !Set neighborhood of the process.
    !Todo:allocate ngb more precise
    if (allocated(proc%ngb))      deallocate(proc%ngb)
    if (allocated(proc%ngb_send)) deallocate(proc%ngb_send)
    if (allocated(proc%ngb_recv)) deallocate(proc%ngb_recv)
    allocate(proc%ngb     (-2:4))
    allocate(proc%ngb_send(-2:30))
    allocate(proc%ngb_recv(-2:30))

    proc%ngb_send_num = 2
    proc%ngb_recv_num = 2

    call MPI_CART_SHIFT(proc%irr_comm%cart_comm, 0, 1, proc%ngb(west )%cart_id, proc%ngb(east )%cart_id, ierr)
    call MPI_CART_SHIFT(proc%irr_comm%cart_comm, 1, 1, proc%ngb(south)%cart_id, proc%ngb(north)%cart_id, ierr)

    ! Translate Cartesian ID of neighbors to global ID.
    do i = 1, 2
      if (proc%ngb(i)%id == MPI_PROC_NULL) then
        call MPI_GROUP_TRANSLATE_RANKS(proc%irr_comm%cart_group, 1, [proc%ngb(i)%cart_id], proc%group, tmp_id, ierr)
        proc%ngb_send(i)%id = tmp_id(1)
        proc%ngb_recv(i)%id = tmp_id(1)
      end if
    end do

    do i = 3, 4
      if (proc%ngb(i)%id == MPI_PROC_NULL) then
        call MPI_GROUP_TRANSLATE_RANKS(proc%irr_comm%cart_group, 1, [proc%ngb(i)%cart_id], proc%group, tmp_id, ierr)
        proc%ngb(i)%id = tmp_id(1)
      end if
    end do

    proc%ngb_send(east)%orient = east
    proc%ngb_send(west)%orient = west
    proc%ngb_recv(east)%orient = east
    proc%ngb_recv(west)%orient = west

    ! Set initial values for num_lon, num_lat, lon_ibeg, lat_ibeg.
    proc%num_lon = num_lon
    proc%num_lat = irr_num_lat(proc%irr_comm%part_id)

    call round_robin(proc%irr_comm%cart_dims(1), proc%irr_comm%cart_coords(1), proc%num_lon, proc%lon_ibeg, proc%lon_iend)
    call round_robin(proc%irr_comm%cart_dims(2), proc%irr_comm%cart_coords(2), proc%num_lat, proc%lat_ibeg, proc%lat_iend)
    proc%lat_ibeg = proc%lat_ibeg + tmp_num_lat(proc%irr_comm%part_id)
    proc%lat_iend = proc%lat_iend + tmp_num_lat(proc%irr_comm%part_id)

    call process_create_blocks_stage1

    call calc_lon_halo_width(proc%lon_halo_width_south , proc%lon_halo_width_north , proc%lon_halo_width )
    proc%lat_halo_width = global_mesh%lat_halo_width

    proc%lat_halo_width_south = proc%lat_halo_width
    proc%lat_halo_width_north = proc%lat_halo_width

    !write(*,*) , 'id = ' ,  proc%id , proc%lon_halo_width , proc%lat_halo_width


    tmp_part = proc%irr_comm%part_id

    !Set west Ngb
      call proc%ngb_send(west)%init(west, send_type, tag = proc%ngb_send(west)%id + east,                                               &
                              lon_ibeg = proc%lon_ibeg                       , lon_iend =proc%lon_ibeg + proc%lon_halo_width - 1, &
                              lat_ibeg = proc%lat_ibeg                       , lat_iend =proc%lat_iend                            )
      call proc%ngb_recv(west)%init(west, recv_type, tag = proc%id + west,                                                              &
                              lon_ibeg = proc%lon_ibeg - proc%lon_halo_width , lon_iend =proc%lon_ibeg - 1 ,                      &
                              lat_ibeg = proc%lat_ibeg                       , lat_iend =proc%lat_iend                            )  
    
    !Set east Ngb
      call proc%ngb_send(east)%init(east, send_type, tag = proc%ngb_send(east)%id + west,                                               &
                              lon_ibeg = proc%lon_iend - proc%lon_halo_width + 1 , lon_iend =proc%lon_iend,                       &
                              lat_ibeg = proc%lat_ibeg                           , lat_iend =proc%lat_iend                        )
      call proc%ngb_recv(east)%init(east, recv_type, tag = proc%id + east ,                                                             &
                              lon_ibeg = proc%lon_iend + 1                       , lon_iend =proc%lon_iend + proc%lon_halo_width, &
                              lat_ibeg = proc%lat_ibeg                           , lat_iend =proc%lat_iend                        )                               

    !Set North Ngb
      proc%ngb_north_send_beg = proc%ngb_send_num + 1
      proc%ngb_north_recv_beg = proc%ngb_recv_num + 1

      if (proc%ngb(north)%id == -1) then !irr boundary
        
        if (proc%lat_iend == num_lat) then
          !north_pole,nothing to do
        else
          tmp_num = irr_num_lat(tmp_part + 1)
          call round_robin(irr_num_proc_lat(tmp_part + 1),0,tmp_num , tmp_lat_beg , tmp_lat_end)

          proc_beg = tmp_num_proc(tmp_part + 1)
          proc_end = proc_beg + irr_num_proc_lon(tmp_part + 1) * irr_num_proc_lat(tmp_part + 1) - 1
          p = 0
          !all proc in ngb_irr_part
          do j = proc_beg , proc_end  , irr_num_proc_lat(tmp_part + 1)
            tmp_num = num_lon
            call round_robin(irr_num_proc_lon(tmp_part + 1), p, tmp_num , tmp_lon_beg, tmp_lon_end)
            call calc_proc_area_intersection(j,north,tmp_lon_beg, tmp_lon_end, tmp_lat_beg, tmp_lat_end)
            p = p + 1
          end do
        endif
      else!inside irr part
        !north west
        proc%ngb_send_num = proc%ngb_send_num + 1
        proc%ngb_recv_num = proc%ngb_recv_num + 1
        proc%ngb_send(proc%ngb_send_num)%id = proc%ngb_send(west)%id + 1
        proc%ngb_recv(proc%ngb_recv_num)%id = proc%ngb_recv(west)%id + 1

        call proc%ngb_send(proc%ngb_send_num)%init(north , send_type, tag = proc%ngb_send(proc%ngb_send_num)%id + east,                       &
                      lon_ibeg = proc%lon_ibeg ,                               lon_iend = proc%lon_ibeg + proc%lon_halo_width_north - 1,&
                      lat_ibeg = proc%lat_iend - proc%lat_halo_width + 1,      lat_iend = proc%lat_iend )
        call proc%ngb_recv(proc%ngb_recv_num)%init(north , recv_type, tag = proc%id + west,                                         &
                      lon_ibeg = proc%lon_ibeg - min(proc%lon_halo_width,proc%lon_halo_width_north), lon_iend = proc%lon_ibeg - 1,                        &
                      lat_ibeg = proc%lat_iend + 1 ,                           lat_iend = proc%lat_iend + proc%lat_halo_width )

        !north
        proc%ngb_send_num = proc%ngb_send_num + 1
        proc%ngb_recv_num = proc%ngb_recv_num + 1
        proc%ngb_send(proc%ngb_send_num)%id = proc%id + 1
        proc%ngb_recv(proc%ngb_recv_num)%id = proc%id + 1

        call proc%ngb_send(proc%ngb_send_num)%init(north , send_type, tag = proc%ngb_send(proc%ngb_send_num)%id + south, &
                      lon_ibeg = proc%lon_ibeg                            ,    lon_iend = proc%lon_iend ,&
                      lat_ibeg = proc%lat_iend - proc%lat_halo_width + 1,      lat_iend = proc%lat_iend  )
        call proc%ngb_recv(proc%ngb_recv_num)%init(north , recv_type, tag = proc%id + north ,                            &
                      lon_ibeg = proc%lon_ibeg                            ,    lon_iend = proc%lon_iend , &
                      lat_ibeg = proc%lat_iend + 1 ,                           lat_iend = proc%lat_iend + proc%lat_halo_width)

        !north east
        proc%ngb_send_num = proc%ngb_send_num + 1
        proc%ngb_recv_num = proc%ngb_recv_num + 1
        proc%ngb_send(proc%ngb_send_num)%id = proc%ngb_send(east)%id + 1
        proc%ngb_recv(proc%ngb_recv_num)%id = proc%ngb_recv(east)%id + 1

        call proc%ngb_send(proc%ngb_send_num)%init(north, send_type , tag = proc%ngb_send(proc%ngb_send_num)%id + west,   &
                      lon_ibeg = proc%lon_iend - proc%lon_halo_width_north + 1, lon_iend = proc%lon_iend,&
                      lat_ibeg = proc%lat_iend - proc%lat_halo_width + 1      , lat_iend = proc%lat_iend )
        call proc%ngb_recv(proc%ngb_recv_num)%init(north, recv_type , tag = proc%id + east,                               &
                      lon_ibeg = proc%lon_iend + 1                            , lon_iend = proc%lon_iend + min(proc%lon_halo_width,proc%lon_halo_width_north),&
                      lat_ibeg = proc%lat_iend + 1 ,                            lat_iend = proc%lat_iend + proc%lat_halo_width )
      
      end if
    !north end

    !Set South Ngb
      proc%ngb_south_send_beg = proc%ngb_send_num + 1
      proc%ngb_south_recv_beg = proc%ngb_recv_num + 1

      if (proc%ngb(south)%id == -1) then !irr boundary
        if (proc%lat_ibeg == 1) then
          !north_pole,nothing to do
        else
          tmp_num = irr_num_lat(tmp_part - 1)
          call round_robin(irr_num_proc_lat(tmp_part - 1),0,tmp_num , tmp_lat_beg , tmp_lat_end)
          
          proc_beg = tmp_num_proc(tmp_part - 1) + irr_num_proc_lat(tmp_part - 1) - 1
          proc_end = tmp_num_proc(tmp_part) - 1
          p = 0

          do j = proc_beg , proc_end , irr_num_proc_lat(tmp_part - 1)
            tmp_num = num_lon
            call round_robin(irr_num_proc_lon(tmp_part - 1), p, tmp_num , tmp_lon_beg, tmp_lon_end)
            call calc_proc_area_intersection(j,south,tmp_lon_beg, tmp_lon_end, tmp_lat_beg, tmp_lat_end)
            p = p + 1
          end do
        endif 
      else !inside irr part
        !south west
        proc%ngb_send_num = proc%ngb_send_num + 1
        proc%ngb_recv_num = proc%ngb_recv_num + 1
        proc%ngb_send(proc%ngb_send_num)%id = proc%ngb_send(west)%id - 1
        proc%ngb_recv(proc%ngb_recv_num)%id = proc%ngb_recv(west)%id - 1  

        call proc%ngb_send(proc%ngb_send_num)%init(south , send_type , tag = proc%ngb_send(proc%ngb_send_num)%id + east,                   &
                      lon_ibeg = proc%lon_ibeg                            , lon_iend = proc%lon_ibeg + proc%lon_halo_width_south - 1,&
                      lat_ibeg = proc%lat_ibeg                            , lat_iend = proc%lat_ibeg + proc%lat_halo_width - 1 )
        call proc%ngb_recv(proc%ngb_recv_num)%init(south , recv_type , tag = proc%id + west,             &
                      lon_ibeg = proc%lon_ibeg - min(proc%lon_halo_width,proc%lon_halo_width_south), lon_iend = proc%lon_ibeg - 1,&
                      lat_ibeg = proc%lat_ibeg - proc%lat_halo_width, lat_iend = proc%lat_ibeg - 1 )
      
        !south
        proc%ngb_send_num = proc%ngb_send_num + 1
        proc%ngb_recv_num = proc%ngb_recv_num + 1
        proc%ngb_send(proc%ngb_send_num)%id = proc%id - 1
        proc%ngb_recv(proc%ngb_recv_num)%id = proc%id - 1

        call proc%ngb_send(proc%ngb_send_num)%init(south , send_type , tag = proc%ngb_send(proc%ngb_send_num)%id + north,                &
                      lon_ibeg = proc%lon_ibeg                            , lon_iend = proc%lon_iend,         &
                      lat_ibeg = proc%lat_ibeg                            , lat_iend = proc%lat_ibeg + proc%lat_halo_width - 1)
        call proc%ngb_recv(proc%ngb_recv_num)%init(south , recv_type , tag = proc%id + south,                   &
                      lon_ibeg = proc%lon_ibeg                            , lon_iend = proc%lon_iend,     &
                      lat_ibeg = proc%lat_ibeg - proc%lat_halo_width, lat_iend = proc%lat_ibeg - 1  )

        !south east
        proc%ngb_send_num = proc%ngb_send_num + 1
        proc%ngb_recv_num = proc%ngb_recv_num + 1
        proc%ngb_send(proc%ngb_send_num)%id = proc%ngb_send(east)%id - 1
        proc%ngb_recv(proc%ngb_recv_num)%id = proc%ngb_recv(east)%id - 1

        call proc%ngb_send(proc%ngb_send_num)%init(south , send_type , tag = proc%ngb_send(proc%ngb_send_num)%id + west,              &
                      lon_ibeg = proc%lon_iend - proc%lon_halo_width_south + 1 , lon_iend = proc%lon_iend,   &
                      lat_ibeg = proc%lat_ibeg                            ,      lat_iend = proc%lat_ibeg + proc%lat_halo_width - 1)
        call proc%ngb_recv(proc%ngb_recv_num)%init(south , recv_type , tag = proc%id + east,                                           &
                      lon_ibeg = proc%lon_iend + 1 ,                             lon_iend = proc%lon_iend + min(proc%lon_halo_width,proc%lon_halo_width_south), &
                      lat_ibeg = proc%lat_ibeg - proc%lat_halo_width,            lat_iend = proc%lat_ibeg - 1)
      endif
    !south end
    
    ! do i = 1 , proc%ngb_send_num
    !   if (proc%ngb_send(i)%orient == south) then
    !     write(*,*) proc%id ,'south s', proc%ngb_send(i)%id ,proc%ngb_send(i)%tag , proc%ngb_send(i)%lon_ibeg , proc%ngb_send(i)%lon_iend , proc%ngb_send(i)%lat_ibeg , proc%ngb_send(i)%lat_iend
    !   else if (proc%ngb_send(i)%orient == north) then
    !     write(*,*) proc%id ,'north s', proc%ngb_send(i)%id ,proc%ngb_send(i)%tag ,  proc%ngb_send(i)%lon_ibeg , proc%ngb_send(i)%lon_iend , proc%ngb_send(i)%lat_ibeg , proc%ngb_send(i)%lat_iend
    !   else if (proc%ngb_send(i)%orient == west) then
    !     write(*,*) proc%id ,'west s' , proc%ngb_send(i)%id ,proc%ngb_send(i)%tag ,  proc%ngb_send(i)%lon_ibeg , proc%ngb_send(i)%lon_iend , proc%ngb_send(i)%lat_ibeg , proc%ngb_send(i)%lat_iend
    !   else
    !     write(*,*) proc%id ,'east s' , proc%ngb_send(i)%id ,proc%ngb_send(i)%tag ,  proc%ngb_send(i)%lon_ibeg , proc%ngb_send(i)%lon_iend , proc%ngb_send(i)%lat_ibeg , proc%ngb_send(i)%lat_iend
    !   end if
    ! end do

    ! do i = 1 , proc%ngb_recv_num
    !   if (proc%ngb_recv(i)%orient == south) then
    !     write(*,*) proc%id ,'south r', proc%ngb_recv(i)%id ,proc%ngb_recv(i)%tag ,  proc%ngb_recv(i)%lon_ibeg , proc%ngb_recv(i)%lon_iend , proc%ngb_recv(i)%lat_ibeg , proc%ngb_recv(i)%lat_iend
    !   else if (proc%ngb_recv(i)%orient == north) then
    !     write(*,*) proc%id ,'north r', proc%ngb_recv(i)%id ,proc%ngb_recv(i)%tag ,  proc%ngb_recv(i)%lon_ibeg , proc%ngb_recv(i)%lon_iend , proc%ngb_recv(i)%lat_ibeg , proc%ngb_recv(i)%lat_iend
    !   else if (proc%ngb_recv(i)%orient == west) then
    !     write(*,*) proc%id ,'west r' , proc%ngb_recv(i)%id ,proc%ngb_recv(i)%tag ,  proc%ngb_recv(i)%lon_ibeg , proc%ngb_recv(i)%lon_iend , proc%ngb_recv(i)%lat_ibeg , proc%ngb_recv(i)%lat_iend
    !   else
    !     write(*,*) proc%id ,'east r' , proc%ngb_recv(i)%id ,proc%ngb_recv(i)%tag ,  proc%ngb_recv(i)%lon_ibeg , proc%ngb_recv(i)%lon_iend , proc%ngb_recv(i)%lat_ibeg , proc%ngb_recv(i)%lat_iend
    !   end if
    ! end do

    ! write(*,*) 'proc' , proc%id ,  proc%lon_ibeg , proc%lon_iend , proc%lat_ibeg , proc%lat_iend  
    
  end subroutine decompose_domains_all_ngb

  subroutine setup_zonal_comm_for_reduce()

    ! Create zonal communicator for reduce algorithm.
    if (proc%idom == 1) then ! Only root domain has reduce region.
      call proc%zonal_circle%init()
    end if

  end subroutine setup_zonal_comm_for_reduce

  subroutine process_create_blocks_stage1()
    
    integer i, j, iv, dtype
    logical flag
    integer max_hw

    if (.not. allocated(proc%blocks)) allocate(proc%blocks(1))

    call proc%blocks(1)%init_stage_1(proc%id, global_mesh%lon_halo_width, global_mesh%lat_halo_width, &
                                     proc%lon_ibeg, proc%lon_iend, proc%lat_ibeg, proc%lat_iend)

    ! Each process calculate lon_halo_width from its filter%ngrid_lat(:).
    max_hw = 2
    do j = proc%blocks(1)%mesh%half_lat_ibeg, proc%blocks(1)%mesh%half_lat_iend
      max_hw = max(max_hw, ( max(proc%blocks(1)%filter%ngrid_lat(j),proc%blocks(1)%filter%ngrid_lon(j)) - 1) / 2)
    end do
    proc%lon_halo_width = max(max_hw, global_mesh%lon_halo_width)
    call global_mesh%reinit(proc%lon_halo_width)

  end subroutine process_create_blocks_stage1

  subroutine process_create_blocks_stage2()
    
    integer i, j, iv, dtype
    logical flag
     

    call proc%blocks(1)%init_stage_2(global_mesh%lon_halo_width)

    ! call proc%blocks(1)%init(proc%id, global_mesh%lon_halo_width, &
    !                          global_mesh%lat_halo_width_north,global_mesh%lat_halo_width_south, &
    !                          proc%regional_north , proc%regional_south , &
    !                          proc%lon_ibeg, proc%lon_iend, proc%lat_ibeg, proc%lat_iend)

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
    allocate(proc%blocks(1)%halo_send(proc%ngb_send_num))
    allocate(proc%blocks(1)%halo_Recv(proc%ngb_recv_num))

    do i = 1, proc%ngb_send_num

      call proc%blocks(1)%halo_send(i)%init(proc%blocks(1)%mesh, proc%ngb_send(i)%orient,dtype, &
                                            host_id=proc%id, ngb_proc_id=proc%ngb_send(i)%id, &
                                            tag = proc%ngb_send(i)%tag, &
                                            lon_ibeg = proc%lon_ibeg - proc%lon_halo_width, &
                                            lat_ibeg = proc%lat_ibeg - proc%lat_halo_width, &
                                            halo_lon_ibeg = proc%ngb_send(i)%lon_ibeg, &
                                            halo_lon_iend = proc%ngb_send(i)%lon_iend, &
                                            halo_lat_ibeg = proc%ngb_send(i)%lat_ibeg, &
                                            halo_lat_iend = proc%ngb_send(i)%lat_iend  )
    end do

    do i = 1, proc%ngb_recv_num
      
      call proc%blocks(1)%halo_recv(i)%init(proc%blocks(1)%mesh, proc%ngb_recv(i)%orient,dtype, &
                                            host_id=proc%id, ngb_proc_id=proc%ngb_recv(i)%id, &
                                            tag = proc%ngb_recv(i)%tag, &
                                            lon_ibeg = proc%lon_ibeg - proc%lon_halo_width, &
                                            lat_ibeg = proc%lat_ibeg - proc%lat_halo_width, &
                                            halo_lon_ibeg = proc%ngb_recv(i)%lon_ibeg, &
                                            halo_lon_iend = proc%ngb_recv(i)%lon_iend, &
                                            halo_lat_ibeg = proc%ngb_recv(i)%lat_ibeg, &
                                            halo_lat_iend = proc%ngb_recv(i)%lat_iend  )
    end do

    ! Initialize async objects.
    do iv = 1 , vector_num
      do i = 1, size(proc%blocks(1)%state,1)
        do j = 1, size(proc%blocks(1)%state(i,iv)%async)
          call proc%blocks(1)%state(i,iv)%async(j)%init(proc%ngb_send_num , proc%ngb_recv_num)
        end do
      end do
    end do

  end subroutine process_create_blocks_stage2

  subroutine process_neighbor_init(this, orient, ngb_type, tag, lon_ibeg, lon_iend, lat_ibeg, lat_iend)

    class(process_neighbor_type), intent(inout) :: this
    integer, intent(in) :: orient
    integer, intent(in), optional :: ngb_type
    integer, intent(in), optional :: tag
    integer, intent(in), optional :: lon_ibeg
    integer, intent(in), optional :: lon_iend
    integer, intent(in), optional :: lat_ibeg
    integer, intent(in), optional :: lat_iend

    this%orient      = orient
    this%ngb_type    = ngb_type
    this%tag         = tag

    this%lon_ibeg = lon_ibeg
    this%lon_iend = lon_iend
    this%lat_ibeg = lat_ibeg
    this%lat_iend = lat_iend

  end subroutine process_neighbor_init

  subroutine zonal_circle_init(this)

    class(zonal_circle_type), intent(inout) :: this

    integer ierr, i, num_lon, ibeg, iend
    integer west_cart_id, east_cart_id, tmp_id(1)
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

      ! Get IDs of the west and east neighbors in zonal circle comm.
      call MPI_CART_SHIFT(proc%cart_comm, 0, 1, west_cart_id, east_cart_id, ierr)
      call MPI_GROUP_TRANSLATE_RANKS(proc%cart_group, 1, [west_cart_id], this%group, tmp_id, ierr); this%west_ngb_id = tmp_id(1)
      call MPI_GROUP_TRANSLATE_RANKS(proc%cart_group, 1, [east_cart_id], this%group, tmp_id, ierr); this%east_ngb_id = tmp_id(1)
    else if (partition_type == 'irregular') then
      allocate(zonal_proc_id(proc%irr_comm%cart_dims(1)))
      do i = 1, proc%irr_comm%cart_dims(1)
        call MPI_CART_RANK(proc%irr_comm%cart_comm, [i-1,proc%irr_comm%cart_coords(2)], zonal_proc_id(i), ierr)
      end do
      call MPI_GROUP_INCL(proc%irr_comm%cart_group, size(zonal_proc_id), zonal_proc_id, this%group, ierr)
      call MPI_COMM_CREATE_GROUP(proc%irr_comm%cart_comm, this%group, sum(zonal_proc_id), this%comm, ierr)
      call MPI_COMM_SIZE(this%comm, this%np, ierr)
      call MPI_COMM_RANK(this%comm, this%id, ierr)
      deallocate(zonal_proc_id)

      ! Get IDs of the west and east neighbors in zonal circle comm.
      call MPI_CART_SHIFT(proc%irr_comm%cart_comm, 0, 1, west_cart_id, east_cart_id, ierr)
      call MPI_GROUP_TRANSLATE_RANKS(proc%irr_comm%cart_group, 1, [west_cart_id], this%group, tmp_id, ierr); this%west_ngb_id = tmp_id(1)
      call MPI_GROUP_TRANSLATE_RANKS(proc%irr_comm%cart_group, 1, [east_cart_id], this%group, tmp_id, ierr); this%east_ngb_id = tmp_id(1)
    end if

    

    if (this%id == 0) then
      ! Single precision
      allocate(this%recv_type_r4(this%np,0:2))
      do i = 1, this%np
        num_lon = global_mesh%num_full_lon
        call round_robin(this%np, i - 1, num_lon, ibeg, iend)
        call MPI_TYPE_CREATE_SUBARRAY(2, [member_num , global_mesh%num_full_lon], &
                                         [member_num ,             num_lon], &
                                         [0,ibeg-1], MPI_ORDER_FORTRAN, MPI_REAL, &
                                         this%recv_type_r4(i,0), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_r4(i,0), ierr)
        call MPI_TYPE_CREATE_SUBARRAY(3, [member_num , global_mesh%num_full_lon,global_mesh%num_full_lev], &
                                         [member_num ,                  num_lon,global_mesh%num_full_lev], &
                                         [0,ibeg-1,0], MPI_ORDER_FORTRAN, MPI_REAL, &
                                         this%recv_type_r4(i,1), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_r4(i,1), ierr)
        call MPI_TYPE_CREATE_SUBARRAY(3, [member_num , global_mesh%num_full_lon,global_mesh%num_half_lev], &
                                         [member_num ,                  num_lon,global_mesh%num_half_lev], &
                                         [0,ibeg-1,0], MPI_ORDER_FORTRAN, MPI_REAL, &
                                         this%recv_type_r4(i,2), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_r4(i,2), ierr)
      end do
      ! Double precision
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

    if (allocated(this%recv_type_r4)) then
      do k = 0, 2
        do i = 1, this%np
          call MPI_TYPE_FREE(this%recv_type_r4(i,k), ierr)
        end do
        deallocate(this%recv_type_r4)
      end do
    end if

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

  subroutine calc_lon_halo_width(lon_halo_width_south , lon_halo_width_north , lon_halo_width)

    integer, intent(inout) :: lon_halo_width_south
    integer, intent(inout) :: lon_halo_width_north
    integer, intent(inout) :: lon_halo_width


    integer res_num, tmp_num, i,j,ngb_lat_ibeg, ngb_lat_iend, tmp_num_lat 
    integer max_hw 


    !south_halo_width 
    ! if (proc%lat_ibeg == 1) then
    !   lon_halo_width_south = (proc%blocks(1)%filter%ngrid_lat(1) - 1) / 2 
    ! else
    !   if (proc%irr_comm%cart_coords(2) == 0) then
    !     tmp_num_lat = irr_num_lat(proc%irr_comm%part_id - 1)
    !     call round_robin(irr_num_proc_lat(proc%irr_comm%part_id - 1), irr_num_proc_lat(proc%irr_comm%part_id - 1) - 1, tmp_num_lat, ngb_lat_ibeg, ngb_lat_iend)
    !     if (proc%irr_comm%part_id > 2) then
    !       ngb_lat_ibeg = ngb_lat_ibeg + sum(irr_num_lat(1:proc%irr_comm%part_id - 2))
    !       ngb_lat_iend = ngb_lat_iend + sum(irr_num_lat(1:proc%irr_comm%part_id - 2))
    !     end if
    !   else! in same irr part
    !     tmp_num_lat = irr_num_lat(proc%irr_comm%part_id)
    !     call round_robin(irr_num_proc_lat(proc%irr_comm%part_id), proc%irr_comm%cart_coords(2) - 1, tmp_num_lat, ngb_lat_ibeg, ngb_lat_iend)
    !     if (proc%irr_comm%part_id > 1) then
    !       ngb_lat_ibeg = ngb_lat_ibeg + sum(irr_num_lat(1:proc%irr_comm%part_id - 1))
    !       ngb_lat_iend = ngb_lat_iend + sum(irr_num_lat(1:proc%irr_comm%part_id - 1))
    !     end if
    !   end if
      
    !   max_hw = 2
    !   do j = ngb_lat_ibeg,ngb_lat_iend
    !     max_hw = max(max_hw, (proc%blocks(1)%filter%ngrid_lat(j) - 1) / 2)
    !   end do
    lon_halo_width_south = 2

    ! end if


    ! !north_halo_width
    ! if (proc%lat_iend == num_lat) then
    !   lon_halo_width_south = (proc%blocks(1)%filter%ngrid_lat(num_lat) - 1) / 2 
    ! else
    !   if (proc%irr_comm%cart_coords(2) ==  proc%irr_comm%cart_dims(2) - 1) then
    !     tmp_num_lat = irr_num_lat(proc%irr_comm%part_id + 1)
    !     call round_robin(irr_num_proc_lat(proc%irr_comm%part_id + 1), 0 , tmp_num_lat, ngb_lat_ibeg, ngb_lat_iend)
    !     ngb_lat_ibeg = ngb_lat_ibeg + sum(irr_num_lat(1:proc%irr_comm%part_id))
    !     ngb_lat_iend = ngb_lat_iend + sum(irr_num_lat(1:proc%irr_comm%part_id))
    !   else! in same irr part
    !     tmp_num_lat = irr_num_lat(proc%irr_comm%part_id)
    !     call round_robin(irr_num_proc_lat(proc%irr_comm%part_id), proc%irr_comm%cart_coords(2) + 1, tmp_num_lat, ngb_lat_ibeg, ngb_lat_iend)
    !     if (proc%irr_comm%part_id > 1) then
    !       ngb_lat_ibeg = ngb_lat_ibeg + sum(irr_num_lat(1:proc%irr_comm%part_id - 1))
    !       ngb_lat_iend = ngb_lat_iend + sum(irr_num_lat(1:proc%irr_comm%part_id - 1))
    !     end if
    !   end if

    !   max_hw = 2
    !   do j = ngb_lat_ibeg,ngb_lat_iend
    !     max_hw = max(max_hw, (proc%blocks(1)%filter%ngrid_lat(j) - 1) / 2)
    !   end do
      lon_halo_width_north = 2

    ! end if

    !my_halo_width
    ! max_hw = 2
    ! do j = proc%lat_ibeg, proc%lat_iend
    !   max_hw = max(max_hw, (max(proc%blocks(1)%filter%ngrid_lat(j),proc%blocks(1)%filter%ngrid_lon(j)) - 1) / 2)
    ! end do
    ! lon_halo_width = max(max_hw, global_mesh%lon_halo_width)
    
    ! if (lon_halo_width_south < lon_halo_width)  lon_halo_width_south = lon_halo_width
    ! if (lon_halo_width_north < lon_halo_width)  lon_halo_width_north = lon_halo_width

    ! lon_halo_width = max(lon_halo_width_south , lon_halo_width_north)

  end subroutine calc_lon_halo_width

  subroutine calc_proc_area_intersection(ngb_id, orient, ngb_lon_ibeg, ngb_lon_iend, ngb_lat_ibeg, ngb_lat_iend)
    integer, intent(in) :: ngb_id
    integer, intent(in) :: orient
    integer, intent(in) :: ngb_lon_ibeg
    integer, intent(in) :: ngb_lon_iend
    integer, intent(in) :: ngb_lat_ibeg
    integer, intent(in) :: ngb_lat_iend

    integer halo_lon_ibeg , halo_lon_iend , halo_lon_ibeg_normal , halo_lon_iend_normal
    integer buf_lon_ibeg(3), buf_lon_iend(3), buf_tag(3), buf_lat_ibeg, buf_lat_iend
    integer buf_num
    integer i

    if (orient == north) then !north
      
      !send
        buf_num = 0
        buf_lat_ibeg  = proc%lat_iend - proc%lat_halo_width + 1
        buf_lat_iend  = proc%lat_iend
        halo_lon_ibeg_normal = -1
        halo_lon_iend_normal = -1
    
        if (ngb_lon_ibeg == 1) then !包含西循环边界
          halo_lon_ibeg = num_lon - proc%lon_halo_width_north + 1
          halo_lon_iend = num_lon
          call calc_proc_buf_halo(normal, ngb_id + west, buf_tag, buf_lon_ibeg, buf_lon_iend, buf_num, halo_lon_ibeg, halo_lon_iend, proc%lon_ibeg, proc%lon_iend)
          halo_lon_ibeg_normal = ngb_lon_ibeg
        end if
         
        if(ngb_lon_iend == num_lon) then!包含东循环边界
          halo_lon_ibeg = 1
          halo_lon_iend = proc%lon_halo_width_north
          call calc_proc_buf_halo(normal, ngb_id + east, buf_tag, buf_lon_ibeg, buf_lon_iend, buf_num, halo_lon_ibeg, halo_lon_iend, proc%lon_ibeg, proc%lon_iend)
          halo_lon_iend_normal = ngb_lon_iend
        end if

        !常规部分
        if (halo_lon_ibeg_normal == -1) halo_lon_ibeg_normal = ngb_lon_ibeg - proc%lon_halo_width_north
        if (halo_lon_iend_normal == -1) halo_lon_iend_normal = ngb_lon_iend + proc%lon_halo_width_north
        call calc_proc_buf_halo(normal, ngb_id + south, buf_tag, buf_lon_ibeg, buf_lon_iend, buf_num, halo_lon_ibeg_normal, halo_lon_iend_normal, proc%lon_ibeg, proc%lon_iend) 

        do i = 1 , buf_num
          proc%ngb_send_num = proc%ngb_send_num + 1
          proc%ngb_send(proc%ngb_send_num)%id = ngb_id
          call proc%ngb_send(proc%ngb_send_num)%init(north, send_type, tag = buf_tag(i), &
                    lon_ibeg = buf_lon_ibeg(i) , lon_iend = buf_lon_iend(i) , &
                    lat_ibeg = buf_lat_ibeg    , lat_iend = buf_lat_iend      )
        end do
       
      !recv
        buf_num = 0
        buf_lat_ibeg = proc%lat_iend + 1
        buf_lat_iend = proc%lat_iend + proc%lat_halo_width_north
        halo_lon_ibeg_normal = -1
        halo_lon_iend_normal = -1

        if (proc%lon_ibeg == 1) then !包含西循环边界
          halo_lon_ibeg = num_lon - min(proc%lon_halo_width,proc%lon_halo_width_north) + 1
          halo_lon_iend = num_lon
          call calc_proc_buf_halo(east_cycle, proc%id + west, buf_tag, buf_lon_ibeg, buf_lon_iend, buf_num, ngb_lon_ibeg, ngb_lon_iend, halo_lon_ibeg, halo_lon_iend)
          halo_lon_ibeg_normal = proc%lon_ibeg
        end if

        if (proc%lon_iend == num_lon) then !包含东循环边界
          halo_lon_ibeg = 1
          halo_lon_iend = min(proc%lon_halo_width,proc%lon_halo_width_north)
          call calc_proc_buf_halo(west_cycle, proc%id + east, buf_tag, buf_lon_ibeg, buf_lon_iend, buf_num, ngb_lon_ibeg, ngb_lon_iend, halo_lon_ibeg, halo_lon_iend)
          halo_lon_iend_normal = proc%lon_iend
        end if
         
        !常规部分
        if (halo_lon_ibeg_normal == -1) halo_lon_ibeg_normal = proc%lon_ibeg - min(proc%lon_halo_width,proc%lon_halo_width_north)
        if (halo_lon_iend_normal == -1) halo_lon_iend_normal = proc%lon_iend + min(proc%lon_halo_width,proc%lon_halo_width_north)
        call calc_proc_buf_halo(normal, proc%id + north, buf_tag, buf_lon_ibeg, buf_lon_iend, buf_num, ngb_lon_ibeg, ngb_lon_iend, halo_lon_ibeg_normal, halo_lon_iend_normal)

        do i = 1 , buf_num
          proc%ngb_recv_num = proc%ngb_recv_num + 1
          proc%ngb_recv(proc%ngb_recv_num)%id = ngb_id
          call proc%ngb_recv(proc%ngb_recv_num)%init(north, recv_type, tag = buf_tag(i), &
                    lon_ibeg = buf_lon_ibeg(i) , lon_iend = buf_lon_iend(i) , &
                    lat_ibeg = buf_lat_ibeg    , lat_iend = buf_lat_iend      )
        end do
        
    end if !north end

    if (orient == south) then!south 

      !send 
        buf_num = 0
        buf_lat_ibeg = proc%lat_ibeg
        buf_lat_iend = proc%lat_ibeg + proc%lat_halo_width - 1
        halo_lon_ibeg_normal = -1
        halo_lon_iend_normal = -1
      
        if (ngb_lon_ibeg == 1) then !包含西循环边界
          halo_lon_ibeg = num_lon - proc%lon_halo_width_south + 1
          halo_lon_iend = num_lon
          call calc_proc_buf_halo(normal,ngb_id + west, buf_tag, buf_lon_ibeg, buf_lon_iend, buf_num, halo_lon_ibeg, halo_lon_iend, proc%lon_ibeg, proc%lon_iend)
          halo_lon_ibeg_normal = ngb_lon_ibeg
        end if
       
        if (ngb_lon_iend == num_lon) then!包含东循环边界
          halo_lon_ibeg = 1
          halo_lon_iend = proc%lon_halo_width_south
          call calc_proc_buf_halo(normal,ngb_id + east, buf_tag, buf_lon_ibeg, buf_lon_iend, buf_num, halo_lon_ibeg, halo_lon_iend, proc%lon_ibeg, proc%lon_iend)
          halo_lon_iend_normal = ngb_lon_iend
        end if
        
        !常规部分
        if (halo_lon_ibeg_normal == -1) halo_lon_ibeg_normal = ngb_lon_ibeg - proc%lon_halo_width_south
        if (halo_lon_iend_normal == -1) halo_lon_iend_normal = ngb_lon_iend + proc%lon_halo_width_south  
        call calc_proc_buf_halo(normal,ngb_id + north, buf_tag, buf_lon_ibeg, buf_lon_iend, buf_num, halo_lon_ibeg_normal, halo_lon_iend_normal, proc%lon_ibeg, proc%lon_iend)

        do i = 1 , buf_num
          proc%ngb_send_num = proc%ngb_send_num + 1
          proc%ngb_send(proc%ngb_send_num)%id = ngb_id
          call proc%ngb_send(proc%ngb_send_num)%init(south, send_type, tag = buf_tag(i), &
                    lon_ibeg = buf_lon_ibeg(i) , lon_iend = buf_lon_iend(i) , &
                    lat_ibeg = buf_lat_ibeg    , lat_iend = buf_lat_iend      )
        end do

      !recv
        buf_num = 0
        buf_lat_ibeg = proc%lat_ibeg - proc%lat_halo_width
        buf_lat_iend = proc%lat_ibeg - 1
        halo_lon_ibeg_normal = -1
        halo_lon_iend_normal = -1

        if (proc%lon_ibeg == 1) then !包含西循环边界
          halo_lon_ibeg = num_lon - min(proc%lon_halo_width,proc%lon_halo_width_south) + 1
          halo_lon_iend = num_lon
          call calc_proc_buf_halo(east_cycle,proc%id + west, buf_tag, buf_lon_ibeg, buf_lon_iend, buf_num, ngb_lon_ibeg, ngb_lon_iend, halo_lon_ibeg, halo_lon_iend)
          halo_lon_ibeg_normal = proc%lon_ibeg
        end if
         
        if (proc%lon_iend == num_lon) then!包含东循环边界
          halo_lon_ibeg = 1
          halo_lon_iend = min(proc%lon_halo_width,proc%lon_halo_width_south)
          call calc_proc_buf_halo(west_cycle,proc%id + east, buf_tag, buf_lon_ibeg, buf_lon_iend, buf_num, ngb_lon_ibeg, ngb_lon_iend, halo_lon_ibeg, halo_lon_iend)
          halo_lon_iend_normal = proc%lon_iend
        end if

        !常规部分
        if (halo_lon_ibeg_normal == -1) halo_lon_ibeg_normal = proc%lon_ibeg - min(proc%lon_halo_width,proc%lon_halo_width_south)
        if (halo_lon_iend_normal == -1) halo_lon_iend_normal = proc%lon_iend + min(proc%lon_halo_width,proc%lon_halo_width_south)
        call calc_proc_buf_halo(normal,proc%id + south, buf_tag, buf_lon_ibeg, buf_lon_iend, buf_num, ngb_lon_ibeg, ngb_lon_iend, halo_lon_ibeg_normal, halo_lon_iend_normal)

        do i = 1 , buf_num
          proc%ngb_recv_num = proc%ngb_recv_num + 1
          proc%ngb_recv(proc%ngb_recv_num)%id = ngb_id
          call proc%ngb_recv(proc%ngb_recv_num)%init(south, recv_type, tag = buf_tag(i), &
                    lon_ibeg = buf_lon_ibeg(i) , lon_iend = buf_lon_iend(i) , &
                    lat_ibeg = buf_lat_ibeg    , lat_iend = buf_lat_iend      )
        end do

      end if!south end
 
  end subroutine calc_proc_area_intersection

  subroutine calc_proc_buf_halo(halo_type, halo_tag, buf_tag, buf_lon_ibeg, buf_lon_iend, buf_num, ngb_lon_ibeg, ngb_lon_iend, proc_lon_ibeg, proc_lon_iend)
    integer, intent(in)    :: halo_type
    integer, intent(in)    :: halo_tag
    integer, intent(inout) :: buf_tag(:)
    integer, intent(inout) :: buf_lon_ibeg(:)
    integer, intent(inout) :: buf_lon_iend(:)
    integer, intent(inout) :: buf_num
    integer, intent(in) :: ngb_lon_ibeg
    integer, intent(in) :: ngb_lon_iend
    integer, intent(in) :: proc_lon_ibeg
    integer, intent(in) :: proc_lon_iend

    if (.not. (ngb_lon_iend < proc_lon_ibeg .or. proc_lon_iend < ngb_lon_ibeg)) then
      buf_num = buf_num + 1
      buf_lon_ibeg(buf_num) = max(ngb_lon_ibeg , proc_lon_ibeg)
      buf_lon_iend(buf_num) = min(ngb_lon_iend , proc_lon_iend)
      buf_tag(buf_num)      = halo_tag

      if (halo_type == east_cycle) then
        buf_lon_ibeg(buf_num) = buf_lon_ibeg(buf_num) - num_lon
        buf_lon_iend(buf_num) = buf_lon_iend(buf_num) - num_lon
      else if (halo_type == west_cycle) then
        buf_lon_ibeg(buf_num) = buf_lon_ibeg(buf_num) + num_lon
        buf_lon_iend(buf_num) = buf_lon_iend(buf_num) + num_lon
      end if
    end if 
  
  end subroutine calc_proc_buf_halo


end module process_mod
