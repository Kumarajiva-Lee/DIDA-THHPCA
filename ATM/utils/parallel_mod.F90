module parallel_mod

  use mpi
  use flogger
  use const_mod
  use mesh_mod
  use block_mod
  use process_mod
  use parallel_types_mod
  use namelist_mod

  implicit none

  private

  public proc
  public process_init
  public process_stop
  public process_final
  public is_root_proc
  public async_type

  public fill_zonal_halo
  public fill_zonal_halo_member
  ! public fill_reduced_chunk_halo_member
  public fill_halo
  public fill_halo_member
  public zero_halo
  public zero_halo_member
  public zonal_sum
  public zonal_sum_member
  public zonal_sum_ensure_order
  public global_sum
  ! public overlay_inner_halo
  ! public overlay_inner_halo_member
  ! public overlay_inner_halo_chunkwise_member
  public barrier
  public gather_zonal_array
  public scatter_zonal_array

  interface fill_zonal_halo
    module procedure fill_zonal_halo_1d_r8
    module procedure fill_zonal_halo_1d_r8_async
    module procedure fill_zonal_halo_2d_r8
    module procedure fill_zonal_halo_2d_r8_async
  end interface fill_zonal_halo

  interface fill_zonal_halo_member
    module procedure fill_zonal_halo_1d_r4_member
    module procedure fill_zonal_halo_1d_r8_member
    module procedure fill_zonal_halo_1d_r8_async_member
    module procedure fill_zonal_halo_2d_r4_member
    module procedure fill_zonal_halo_2d_r8_member
    module procedure fill_zonal_halo_2d_r8_async_member
  end interface fill_zonal_halo_member

  ! interface fill_reduced_chunk_halo_member
  !   module procedure fill_reduced_chunk_halo_2d_r4_member
  !   module procedure fill_reduced_chunk_halo_2d_r8_member
  !   module procedure fill_reduced_chunk_halo_3d_r4_member
  !   module procedure fill_reduced_chunk_halo_3d_r8_member    
  ! end interface fill_reduced_chunk_halo_member

  interface fill_halo
    module procedure fill_halo_2d_r8
    module procedure fill_halo_2d_r8_async
    module procedure fill_halo_3d_r8
    module procedure fill_halo_3d_r8_async
  end interface fill_halo

  interface fill_halo_member
    module procedure fill_halo_2d_r4_member
    module procedure fill_halo_2d_r8_member
    module procedure fill_halo_2d_r8_async_member
    module procedure fill_halo_3d_r4_member
    module procedure fill_halo_3d_r8_member
    module procedure fill_halo_3d_r8_async_member
  end interface fill_halo_member

  interface zero_halo
    module procedure zero_halo_1d_r8
  end interface zero_halo

  interface zero_halo_member
    module procedure zero_halo_1d_r4_member
    module procedure zero_halo_1d_r8_member
  end interface zero_halo_member

  ! interface overlay_inner_halo_member
  !   module procedure overlay_inner_halo_r4_member
  !   module procedure overlay_inner_halo_r8_member
  ! end interface overlay_inner_halo_member

  ! interface overlay_inner_halo_chunkwise_member
  !   module procedure overlay_inner_halo_chunkwise_r4_member
  !   module procedure overlay_inner_halo_chunkwise_r8_member
  ! end interface overlay_inner_halo_chunkwise_member

  interface zonal_sum_old
    module procedure zonal_sum_old_0d_r8
    module procedure zonal_sum_old_1d_r8
    module procedure zonal_sum_old_2d_r8
    end interface zonal_sum_old

  interface zonal_sum
    module procedure zonal_sum_0d_r8
    module procedure zonal_sum_1d_r8
  end interface zonal_sum

  interface zonal_sum_member
    module procedure zonal_sum_0d_r4_member
    module procedure zonal_sum_0d_r8_member
    module procedure zonal_sum_1d_r4_member
    module procedure zonal_sum_1d_r8_member
  end interface zonal_sum_member

  interface zonal_sum_ensure_order
    module procedure zonal_sum_ensure_order_1d_r8
    module procedure zonal_sum_ensure_order_2d_r8
  end interface zonal_sum_ensure_order

  interface global_sum
    module procedure global_sum_0d_r4
    module procedure global_sum_0d_r8
    module procedure global_sum_1d_r4
    module procedure global_sum_1d_r8
  end interface global_sum

  interface gather_zonal_array  ! need add member
    module procedure gather_zonal_array_1d_r8
    module procedure gather_zonal_array_2d_r8
  end interface gather_zonal_array

  interface scatter_zonal_array ! need add member 
    module procedure scatter_zonal_array_1d_r8
    module procedure scatter_zonal_array_2d_r8
  end interface scatter_zonal_array

  interface gather_zonal_array_member
    module procedure gather_zonal_array_1d_r4_member
    module procedure gather_zonal_array_1d_r8_member
    module procedure gather_zonal_array_2d_r4_member
    module procedure gather_zonal_array_2d_r8_member
  end interface gather_zonal_array_member

  interface scatter_zonal_array_member 
    module procedure scatter_zonal_array_1d_r4_member
    module procedure scatter_zonal_array_1d_r8_member
    module procedure scatter_zonal_array_2d_r4_member
    module procedure scatter_zonal_array_2d_r8_member
  end interface scatter_zonal_array_member

contains

  subroutine fill_zonal_halo_1d_r8(block, halo_width, array, west_halo, east_halo)

    type(block_type), intent(in) :: block
    integer, intent(in) :: halo_width
    real(8), intent(inout)  :: array(:)
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer status(MPI_STATUS_SIZE), ierr
    integer i1, i2, i3, i4

    if (merge(west_halo, .true., present(west_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! | i3  | i4  |     |     |     |     | i1  | i2  |     |  n  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !    |                                   |
      !    1                              n - 2 w + 1
      i1 = size(array) - 2 * halo_width + 1
      i2 = i1 + halo_width - 1
      i3 = 1
      i4 = i3 + halo_width - 1
      call MPI_SENDRECV(array(i1:i2), halo_width, block%halo(east)%dtype, block%halo(east)%proc_id, 1, &
                        array(i3:i4), halo_width, block%halo(west)%dtype, block%halo(west)%proc_id, 1, &
                        proc%comm, status, ierr)
    end if
    if (merge(east_halo, .true., present(east_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! |     |     | i1  | i2  |     |     |     |     | i3  | i4  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !                |                                   |
      !              1 + w                             n - w + 1
      i1 = 1 + halo_width
      i2 = i1 + halo_width - 1
      i3 = size(array) - halo_width + 1
      i4 = i3 + halo_width - 1
      call MPI_SENDRECV(array(i1:i2), halo_width, block%halo(west)%dtype, block%halo(west)%proc_id, 2, &
                        array(i3:i4), halo_width, block%halo(east)%dtype, block%halo(east)%proc_id, 2, &
                        proc%comm, status, ierr)
    end if

  end subroutine fill_zonal_halo_1d_r8

  subroutine fill_zonal_halo_1d_r8_async(block, halo_width, array, async, west_halo, east_halo)

    type(block_type), intent(in) :: block
    integer, intent(in) :: halo_width
    real(8), intent(inout)  :: array(:)
    type(async_type), intent(inout) :: async
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer ierr
    integer i1, i2, i3, i4

    call fill_zonal_halo_1d_r8(block, halo_width, array, west_halo, east_halo)
    return

    if (merge(west_halo, .true., present(west_halo))) then
      if (async%send_req(west) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(west), MPI_STATUS_IGNORE, ierr)
      if (async%recv_req(west) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(west), MPI_STATUS_IGNORE, ierr)
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! | i3  | i4  |     |     |     |     | i1  | i2  |     |  n  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !    |                                   |
      !    1                              n - 2 w + 1
      i1 = size(array) - 2 * halo_width + 1
      i2 = i1 + halo_width - 1
      i3 = 1
      i4 = i3 + halo_width - 1
      call MPI_ISEND(array(i1:i2), halo_width, block%halo(east)%dtype, block%halo(east)%proc_id, 1, proc%comm, async%send_req(west), ierr)
      call MPI_IRECV(array(i3:i4), halo_width, block%halo(west)%dtype, block%halo(west)%proc_id, 1, proc%comm, async%recv_req(west), ierr)
    end if
    if (merge(east_halo, .true., present(east_halo))) then
      if (async%send_req(east) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(east), MPI_STATUS_IGNORE, ierr)
      if (async%recv_req(east) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(east), MPI_STATUS_IGNORE, ierr)
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! |     |     | i1  | i2  |     |     |     |     | i3  | i4  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !                |                                   |
      !              1 + w                             n - w + 1
      i1 = 1 + halo_width
      i2 = i1 + halo_width - 1
      i3 = size(array) - halo_width + 1
      i4 = i3 + halo_width - 1
      call MPI_ISEND(array(i1:i2), halo_width, block%halo(west)%dtype, block%halo(west)%proc_id, 2, proc%comm, async%send_req(east), ierr)
      call MPI_IRECV(array(i3:i4), halo_width, block%halo(east)%dtype, block%halo(east)%proc_id, 2, proc%comm, async%recv_req(east), ierr)
    end if

  end subroutine fill_zonal_halo_1d_r8_async

  subroutine fill_zonal_halo_2d_r8(block, halo_width, array, west_halo, east_halo)

    type(block_type), intent(in) :: block
    integer, intent(in) :: halo_width
    real(8), intent(inout)  :: array(:,:)
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer status(MPI_STATUS_SIZE), ierr
    integer nx, nz, i1, i2, i3, i4

    nz = size(array, 1)
    nx = size(array, 2)

    if (merge(west_halo, .true., present(west_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! | i3  | i4  |     |     |     |     | i1  | i2  |     |  n  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !    |                                   |
      !    1                             nx - 2 w + 1
      i1 = nx - 2 * halo_width + 1
      i2 = i1 + halo_width - 1
      i3 = 1
      i4 = i3 + halo_width - 1
      call MPI_SENDRECV(array(:,i1:i2), halo_width * nz, block%halo(east)%dtype, block%halo(east)%proc_id, 1, &
                        array(:,i3:i4), halo_width * nz, block%halo(west)%dtype, block%halo(west)%proc_id, 1, &
                        proc%comm, status, ierr)
    end if
    if (merge(east_halo, .true., present(east_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! |     |     | i1  | i2  |     |     |     |     | i3  | i4  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !                |                                   |
      !              1 + w                            nx - w + 1
      i1 = 1 + halo_width
      i2 = i1 + halo_width - 1
      i3 = nx - halo_width + 1
      i4 = i3 + halo_width - 1
      call MPI_SENDRECV(array(:,i1:i2), halo_width * nz, block%halo(west)%dtype, block%halo(west)%proc_id, 2, &
                        array(:,i3:i4), halo_width * nz, block%halo(east)%dtype, block%halo(east)%proc_id, 2, &
                        proc%comm, status, ierr)
    end if

  end subroutine fill_zonal_halo_2d_r8

  subroutine fill_zonal_halo_2d_r8_async(block, halo_width, array, async, west_halo, east_halo)

    type(block_type), intent(in) :: block
    integer, intent(in) :: halo_width
    real(8), intent(inout)  :: array(:,:)
    type(async_type), intent(inout) :: async
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer ierr
    integer nx, nz, i1, i2, i3, i4

    call fill_zonal_halo_2d_r8(block, halo_width, array, west_halo, east_halo)
    return

    nz = size(array, 1)
    nx = size(array, 2)

    if (merge(west_halo, .true., present(west_halo))) then
      if (async%send_req(west) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(west), MPI_STATUS_IGNORE, ierr)
      if (async%recv_req(west) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(west), MPI_STATUS_IGNORE, ierr)
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! | i3  | i4  |     |     |     |     | i1  | i2  |     |  n  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !    |                                   |
      !    1                             nx - 2 w + 1
      i1 = nx - 2 * halo_width + 1
      i2 = i1 + halo_width - 1
      i3 = 1
      i4 = i3 + halo_width - 1
      call MPI_ISEND(array(:,i1:i2), halo_width * nz, block%halo(east)%dtype, block%halo(east)%proc_id, 1, proc%comm, async%send_req(west), ierr)
      call MPI_IRECV(array(:,i3:i4), halo_width * nz, block%halo(west)%dtype, block%halo(west)%proc_id, 1, proc%comm, async%recv_req(west), ierr)
    end if
    if (merge(east_halo, .true., present(east_halo))) then
      if (async%send_req(east) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(east), MPI_STATUS_IGNORE, ierr)
      if (async%recv_req(east) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(east), MPI_STATUS_IGNORE, ierr)
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! |     |     | i1  | i2  |     |     |     |     | i3  | i4  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !                |                                   |
      !              1 + w                            nx - w + 1
      i1 = 1 + halo_width
      i2 = i1 + halo_width - 1
      i3 = nx - halo_width + 1
      i4 = i3 + halo_width - 1
      call MPI_ISEND(array(:,i1:i2), halo_width * nz, block%halo(west)%dtype, block%halo(west)%proc_id, 2, proc%comm, async%send_req(east), ierr)
      call MPI_IRECV(array(:,i3:i4), halo_width * nz, block%halo(east)%dtype, block%halo(east)%proc_id, 2, proc%comm, async%recv_req(east), ierr)
    end if

  end subroutine fill_zonal_halo_2d_r8_async

  subroutine fill_zonal_halo_1d_r4_member(block, halo_width, array, west_halo, east_halo)

    type(block_type), intent(in) :: block
    integer, intent(in) :: halo_width
    real(4), intent(inout)  :: array(:,:)
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer status(MPI_STATUS_SIZE), ierr
    integer i1, i2, i3, i4
    integer nm

    if (merge(west_halo, .true., present(west_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! | i3  | i4  |     |     |     |     | i1  | i2  |     |  n  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !    |                                   |
      !    1                              n - 2 w + 1
      
      nm = size(array,1)

      i1 = size(array,2) - 2 * halo_width + 1
      i2 = i1 + halo_width - 1
      i3 = 1
      i4 = i3 + halo_width - 1
      call MPI_SENDRECV(array(:,i1:i2), halo_width * nm, block%halo(east)%dtype, block%halo(east)%proc_id, 1, &
                        array(:,i3:i4), halo_width * nm, block%halo(west)%dtype, block%halo(west)%proc_id, 1, &
                        proc%comm, status, ierr)
    end if
    if (merge(east_halo, .true., present(east_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! |     |     | i1  | i2  |     |     |     |     | i3  | i4  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !                |                                   |
      !              1 + w                             n - w + 1
      i1 = 1 + halo_width
      i2 = i1 + halo_width - 1
      i3 = size(array,2) - halo_width + 1
      i4 = i3 + halo_width - 1
      call MPI_SENDRECV(array(:,i1:i2), halo_width * nm, block%halo(west)%dtype, block%halo(west)%proc_id, 2, &
                        array(:,i3:i4), halo_width * nm, block%halo(east)%dtype, block%halo(east)%proc_id, 2, &
                        proc%comm, status, ierr)
    end if

  end subroutine fill_zonal_halo_1d_r4_member

  subroutine fill_zonal_halo_1d_r8_member(block, halo_width, array, west_halo, east_halo)

    type(block_type), intent(in) :: block
    integer, intent(in) :: halo_width
    real(8), intent(inout)  :: array(:,:)
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer status(MPI_STATUS_SIZE), ierr
    integer i1, i2, i3, i4
    integer nm

    if (merge(west_halo, .true., present(west_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! | i3  | i4  |     |     |     |     | i1  | i2  |     |  n  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !    |                                   |
      !    1                              n - 2 w + 1
      

      nm = size(array,1)

      i1 = size(array,2) - 2 * halo_width + 1
      i2 = i1 + halo_width - 1
      i3 = 1
      i4 = i3 + halo_width - 1
      call MPI_SENDRECV(array(:,i1:i2), halo_width * nm, block%halo_send(east)%dtype, block%halo_send(east)%proc_id, 1, &
                        array(:,i3:i4), halo_width * nm, block%halo_recv(west)%dtype, block%halo_recv(west)%proc_id, 1, &
                        proc%comm, status, ierr)
    end if
    if (merge(east_halo, .true., present(east_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! |     |     | i1  | i2  |     |     |     |     | i3  | i4  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !                |                                   |
      !              1 + w                             n - w + 1
      i1 = 1 + halo_width
      i2 = i1 + halo_width - 1
      i3 = size(array,2) - halo_width + 1
      i4 = i3 + halo_width - 1
      call MPI_SENDRECV(array(:,i1:i2), halo_width * nm, block%halo_send(west)%dtype, block%halo_send(west)%proc_id, 2, &
                        array(:,i3:i4), halo_width * nm, block%halo_recv(east)%dtype, block%halo_recv(east)%proc_id, 2, &
                        proc%comm, status, ierr)
    end if

  end subroutine fill_zonal_halo_1d_r8_member

  subroutine fill_zonal_halo_1d_r8_async_member(block, halo_width, array, async, west_halo, east_halo)

    type(block_type), intent(in) :: block
    integer, intent(in) :: halo_width
    real(8), intent(inout)  :: array(:,:)
    type(async_type), intent(inout) :: async
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer ierr
    integer i1, i2, i3, i4
    integer nm

    call fill_zonal_halo_1d_r8_member(block, halo_width, array, west_halo, east_halo)
    return

    if (merge(west_halo, .true., present(west_halo))) then
      if (async%send_req(west) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(west), MPI_STATUS_IGNORE, ierr)
      if (async%recv_req(west) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(west), MPI_STATUS_IGNORE, ierr)
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! | i3  | i4  |     |     |     |     | i1  | i2  |     |  n  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !    |                                   |
      !    1                              n - 2 w + 1

      nm = size(array , 1)

      i1 = size(array, 2) - 2 * halo_width + 1
      i2 = i1 + halo_width - 1
      i3 = 1
      i4 = i3 + halo_width - 1
      call MPI_ISEND(array(:,i1:i2), halo_width * nm, block%halo(east)%dtype, block%halo(east)%proc_id, 1, proc%comm, async%send_req(west), ierr)
      call MPI_IRECV(array(:,i3:i4), halo_width * nm, block%halo(west)%dtype, block%halo(west)%proc_id, 1, proc%comm, async%recv_req(west), ierr)
    end if
    if (merge(east_halo, .true., present(east_halo))) then
      if (async%send_req(east) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(east), MPI_STATUS_IGNORE, ierr)
      if (async%recv_req(east) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(east), MPI_STATUS_IGNORE, ierr)
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! |     |     | i1  | i2  |     |     |     |     | i3  | i4  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !                |                                   |
      !              1 + w                             n - w + 1
      i1 = 1 + halo_width
      i2 = i1 + halo_width - 1
      i3 = size(array,2) - halo_width + 1
      i4 = i3 + halo_width - 1
      call MPI_ISEND(array(:,i1:i2), halo_width * nm, block%halo(west)%dtype, block%halo(west)%proc_id, 2, proc%comm, async%send_req(east), ierr)
      call MPI_IRECV(array(:,i3:i4), halo_width * nm, block%halo(east)%dtype, block%halo(east)%proc_id, 2, proc%comm, async%recv_req(east), ierr)
    end if

  end subroutine fill_zonal_halo_1d_r8_async_member

  subroutine fill_zonal_halo_2d_r4_member(block, halo_width, array, west_halo, east_halo)

    type(block_type), intent(in) :: block
    integer, intent(in) :: halo_width
    real(4), intent(inout)  :: array(:,:,:)
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer status(MPI_STATUS_SIZE), ierr
    integer nx, nz, nm, i1, i2, i3, i4
    
    nm = size(array, 1)
    nz = size(array, 2)
    nx = size(array, 3)

    if (merge(west_halo, .true., present(west_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! | i3  | i4  |     |     |     |     | i1  | i2  |     |  n  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !    |                                   |
      !    1                             nx - 2 w + 1
      i1 = nx - 2 * halo_width + 1
      i2 = i1 + halo_width - 1
      i3 = 1
      i4 = i3 + halo_width - 1
      call MPI_SENDRECV(array(:,:,i1:i2), halo_width * nz * nm, block%halo(east)%dtype, block%halo(east)%proc_id, 1, &
                        array(:,:,i3:i4), halo_width * nz * nm, block%halo(west)%dtype, block%halo(west)%proc_id, 1, &
                        proc%comm, status, ierr)
    end if
    if (merge(east_halo, .true., present(east_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! |     |     | i1  | i2  |     |     |     |     | i3  | i4  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !                |                                   |
      !              1 + w                            nx - w + 1
      i1 = 1 + halo_width
      i2 = i1 + halo_width - 1
      i3 = nx - halo_width + 1
      i4 = i3 + halo_width - 1
      call MPI_SENDRECV(array(:,:,i1:i2), halo_width * nz * nm, block%halo(west)%dtype, block%halo(west)%proc_id, 2, &
                        array(:,:,i3:i4), halo_width * nz * nm, block%halo(east)%dtype, block%halo(east)%proc_id, 2, &
                        proc%comm, status, ierr)
    end if

  end subroutine fill_zonal_halo_2d_r4_member

  subroutine fill_zonal_halo_2d_r8_member(block, halo_width, array, west_halo, east_halo)

    type(block_type), intent(in) :: block
    integer, intent(in) :: halo_width
    real(8), intent(inout)  :: array(:,:,:)
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer status(MPI_STATUS_SIZE), ierr
    integer nx, nz, nm, i1, i2, i3, i4
    
    nm = size(array, 1)
    nz = size(array, 2)
    nx = size(array, 3)

    if (merge(west_halo, .true., present(west_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! | i3  | i4  |     |     |     |     | i1  | i2  |     |  n  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !    |                                   |
      !    1                             nx - 2 w + 1
      i1 = nx - 2 * halo_width + 1
      i2 = i1 + halo_width - 1
      i3 = 1
      i4 = i3 + halo_width - 1
      call MPI_SENDRECV(array(:,:,i1:i2), halo_width * nz * nm, block%halo_send(east)%dtype, block%halo_send(east)%proc_id, 1, &
                        array(:,:,i3:i4), halo_width * nz * nm, block%halo_recv(west)%dtype, block%halo_recv(west)%proc_id, 1, &
                        proc%comm, status, ierr)
    end if
    if (merge(east_halo, .true., present(east_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! |     |     | i1  | i2  |     |     |     |     | i3  | i4  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !                |                                   |
      !              1 + w                            nx - w + 1
      i1 = 1 + halo_width
      i2 = i1 + halo_width - 1
      i3 = nx - halo_width + 1
      i4 = i3 + halo_width - 1
      call MPI_SENDRECV(array(:,:,i1:i2), halo_width * nz * nm, block%halo_send(west)%dtype, block%halo_send(west)%proc_id, 2, &
                        array(:,:,i3:i4), halo_width * nz * nm, block%halo_recv(east)%dtype, block%halo_recv(east)%proc_id, 2, &
                        proc%comm, status, ierr)
    end if

  end subroutine fill_zonal_halo_2d_r8_member

  subroutine fill_zonal_halo_2d_r8_async_member(block, halo_width, array, async, west_halo, east_halo)

    type(block_type), intent(in) :: block
    integer, intent(in) :: halo_width
    real(8), intent(inout)  :: array(:,:,:)
    type(async_type), intent(inout) :: async
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer ierr
    integer nx, nz, nm ,i1, i2, i3, i4

    call fill_zonal_halo_2d_r8_member(block, halo_width, array, west_halo, east_halo)
    return

    nm = size(array, 1)
    nz = size(array, 2)
    nx = size(array, 3)

    if (merge(west_halo, .true., present(west_halo))) then
      if (async%send_req(west) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(west), MPI_STATUS_IGNORE, ierr)
      if (async%recv_req(west) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(west), MPI_STATUS_IGNORE, ierr)
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! | i3  | i4  |     |     |     |     | i1  | i2  |     |  n  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !    |                                   |
      !    1                             nx - 2 w + 1
      i1 = nx - 2 * halo_width + 1
      i2 = i1 + halo_width - 1
      i3 = 1
      i4 = i3 + halo_width - 1
      call MPI_ISEND(array(:,:,i1:i2), halo_width * nz * nm, block%halo(east)%dtype, block%halo(east)%proc_id, 1, proc%comm, async%send_req(west), ierr)
      call MPI_IRECV(array(:,:,i3:i4), halo_width * nz * nm, block%halo(west)%dtype, block%halo(west)%proc_id, 1, proc%comm, async%recv_req(west), ierr)
    end if
    if (merge(east_halo, .true., present(east_halo))) then
      if (async%send_req(east) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(east), MPI_STATUS_IGNORE, ierr)
      if (async%recv_req(east) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(east), MPI_STATUS_IGNORE, ierr)
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! |     |     | i1  | i2  |     |     |     |     | i3  | i4  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !                |                                   |
      !              1 + w                            nx - w + 1
      i1 = 1 + halo_width
      i2 = i1 + halo_width - 1
      i3 = nx - halo_width + 1
      i4 = i3 + halo_width - 1
      call MPI_ISEND(array(:,:,i1:i2), halo_width * nz * nm, block%halo(west)%dtype, block%halo(west)%proc_id, 2, proc%comm, async%send_req(east), ierr)
      call MPI_IRECV(array(:,:,i3:i4), halo_width * nz * nm, block%halo(east)%dtype, block%halo(east)%proc_id, 2, proc%comm, async%recv_req(east), ierr)
    end if

  end subroutine fill_zonal_halo_2d_r8_async_member

  ! ! fill_reduced_chunk 的2d版本专门为gzs通信设计，因此只有一种MPI_subarray的类型
  ! subroutine fill_reduced_chunk_halo_2d_r4_member(block, reduced_chunk, array, west_halo, east_halo)

  !   type(block_type), intent(in) :: block
  !   type(reduced_chunk_type), intent(in) :: reduced_chunk
  !   real(4), intent(inout) :: array(:,:,:)  ! (im, i, j)
  !   logical, intent(in), optional :: west_halo
  !   logical, intent(in), optional :: east_halo

  !   integer status(MPI_STATUS_SIZE), ierr, near

  !   if (merge(west_halo, .true., present(west_halo))) then  ! 这个变量（的后续计算）需要西侧的halo，因此给东边发，从西边收
  !     do near = 1, reduced_chunk%num_ngbProc
  !       call MPI_SENDRECV(array, 1, reduced_chunk%east_send_type_2d(near), proc%east_pids(near), 100, &
  !                         array, 1, reduced_chunk%west_recv_type_2d(near), proc%west_pids(near), 100, proc%comm, status, ierr)
  !     end do
  !   end if

  !   if (merge(east_halo, .true., present(east_halo))) then   ! 这个变量（的后续计算）需要东侧的halo，因此给西边发，从东边收
  !     do near = 1, reduced_chunk%num_ngbProc
  !       call MPI_SENDRECV(array, 1, reduced_chunk%west_send_type_2d(near), proc%west_pids(near), 101, &
  !                         array, 1, reduced_chunk%east_recv_type_2d(near), proc%east_pids(near), 101, proc%comm, status, ierr)
  !     end do
  !   end if

  ! end subroutine fill_reduced_chunk_halo_2d_r4_member

  ! subroutine fill_reduced_chunk_halo_2d_r8_member(block, reduced_chunk, array, west_halo, east_halo)

  !   type(block_type), intent(in) :: block
  !   type(reduced_chunk_type), intent(in) :: reduced_chunk
  !   real(8), intent(inout) :: array(:,:,:)  ! (im, i, j)
  !   logical, intent(in), optional :: west_halo
  !   logical, intent(in), optional :: east_halo

  !   integer status(MPI_STATUS_SIZE), ierr, near

  !   if (merge(west_halo, .true., present(west_halo))) then  ! 这个变量（的后续计算）需要西侧的halo，因此给东边发，从西边收
  !     do near = 1, reduced_chunk%num_ngbProc
  !       call MPI_SENDRECV(array, 1, reduced_chunk%east_send_type_2d(near), proc%east_pids(near), 100, &
  !                         array, 1, reduced_chunk%west_recv_type_2d(near), proc%west_pids(near), 100, proc%comm, status, ierr)
  !     end do
  !   end if

  !   if (merge(east_halo, .true., present(east_halo))) then   ! 这个变量（的后续计算）需要东侧的halo，因此给西边发，从东边收
  !     do near = 1, reduced_chunk%num_ngbProc
  !       call MPI_SENDRECV(array, 1, reduced_chunk%west_send_type_2d(near), proc%west_pids(near), 101, &
  !                         array, 1, reduced_chunk%east_recv_type_2d(near), proc%east_pids(near), 101, proc%comm, status, ierr)
  !     end do
  !   end if

  ! end subroutine fill_reduced_chunk_halo_2d_r8_member

  ! ! fill_reduced_chunk 的3d版本支持多种MPI_subarray的类型
  ! subroutine fill_reduced_chunk_halo_3d_r4_member(block, reduced_chunk, array, rela_lat_size, full_lev, west_halo, east_halo)

  !   type(block_type), intent(in) :: block
  !   type(reduced_chunk_type), intent(in) :: reduced_chunk
  !   real(4), intent(inout) :: array(:,:,:,:)  ! (im, i, j, k)
  !   integer, intent(in) :: rela_lat_size  ! 调用这个通信函数的上层，决定要传的数组在lat维相对于reduced_chunk%num_full_lat的相对大小
  !   logical, intent(in), optional :: full_lev
  !   logical, intent(in), optional :: west_halo
  !   logical, intent(in), optional :: east_halo

  !   integer status(MPI_STATUS_SIZE), k, ierr, near

  !   k = merge(1, 2, merge(full_lev, .true., present(full_lev))) ! 1: full, 2: half

  !   if (merge(west_halo, .true., present(west_halo))) then
  !     do near = 1, reduced_chunk%num_ngbProc
  !       call MPI_SENDRECV(array, 1, reduced_chunk%east_send_type_3d(rela_lat_size, k, near), proc%east_pids(near), 100, &
  !                         array, 1, reduced_chunk%west_recv_type_3d(rela_lat_size, k, near), proc%west_pids(near), 100, proc%comm, status, ierr)
  !     end do
  !   end if

  !   if (merge(east_halo, .true., present(east_halo))) then
  !     do near = 1, reduced_chunk%num_ngbProc
  !       call MPI_SENDRECV(array, 1, reduced_chunk%west_send_type_3d(rela_lat_size, k, near), proc%west_pids(near), 101, &
  !                         array, 1, reduced_chunk%east_recv_type_3d(rela_lat_size, k, near), proc%east_pids(near), 101, proc%comm, status, ierr)
  !     end do
  !   end if

  ! end subroutine fill_reduced_chunk_halo_3d_r4_member

  ! subroutine fill_reduced_chunk_halo_3d_r8_member(block, reduced_chunk, array, rela_lat_size, full_lev, west_halo, east_halo)

  !   type(block_type), intent(in) :: block
  !   type(reduced_chunk_type), intent(in) :: reduced_chunk
  !   real(8), intent(inout) :: array(:,:,:,:)  ! (im, i, j, k)
  !   integer, intent(in) :: rela_lat_size  ! 调用这个通信函数的上层，决定要传的数组在lat维相对于reduced_chunk%num_full_lat的相对大小
  !   logical, intent(in), optional :: full_lev
  !   logical, intent(in), optional :: west_halo
  !   logical, intent(in), optional :: east_halo

  !   integer status(MPI_STATUS_SIZE), k, ierr, near

  !   k = merge(1, 2, merge(full_lev, .true., present(full_lev))) ! 1: full, 2: half

  !   if (merge(west_halo, .true., present(west_halo))) then
  !     do near = 1, reduced_chunk%num_ngbProc
  !       call MPI_SENDRECV(array, 1, reduced_chunk%east_send_type_3d(rela_lat_size, k, near), proc%east_pids(near), 100, &
  !                         array, 1, reduced_chunk%west_recv_type_3d(rela_lat_size, k, near), proc%west_pids(near), 100, proc%comm, status, ierr)
  !     end do
  !   end if

  !   if (merge(east_halo, .true., present(east_halo))) then
  !     do near = 1, reduced_chunk%num_ngbProc
  !       call MPI_SENDRECV(array, 1, reduced_chunk%west_send_type_3d(rela_lat_size, k, near), proc%west_pids(near), 101, &
  !                         array, 1, reduced_chunk%east_recv_type_3d(rela_lat_size, k, near), proc%east_pids(near), 101, proc%comm, status, ierr)
  !     end do
  !   end if

  ! end subroutine fill_reduced_chunk_halo_3d_r8_member




  subroutine fill_halo_2d_r8(block, array, full_lon, full_lat, west_halo, east_halo, south_halo, north_halo)

    type(block_type), intent(in) :: block
    real(8), intent(inout) :: array(:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    integer status(MPI_STATUS_SIZE), i, j, ierr

    i = merge(1, 2, full_lon)
    j = merge(1, 2, full_lat)

    if (merge(west_halo, .true., present(west_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(east)%send_type_2d(i,j), block%halo(east)%proc_id, 3, &
                        array, 1, block%halo(west)%recv_type_2d(i,j), block%halo(west)%proc_id, 3, &
                        proc%comm, status, ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(west)%send_type_2d(i,j), block%halo(west)%proc_id, 7, &
                        array, 1, block%halo(east)%recv_type_2d(i,j), block%halo(east)%proc_id, 7, &
                        proc%comm, status, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(north)%send_type_2d(i,j), block%halo(north)%proc_id, 11, &
                        array, 1, block%halo(south)%recv_type_2d(i,j), block%halo(south)%proc_id, 11, &
                        proc%comm, status, ierr)
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(south)%send_type_2d(i,j), block%halo(south)%proc_id, 15, &
                        array, 1, block%halo(north)%recv_type_2d(i,j), block%halo(north)%proc_id, 15, &
                        proc%comm, status, ierr)
    end if

  end subroutine fill_halo_2d_r8

  subroutine fill_halo_2d_r8_async(block, array, full_lon, full_lat, async, west_halo, east_halo, south_halo, north_halo)

    type(block_type), intent(in) :: block
    real(8), intent(inout) :: array(:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    type(async_type), intent(inout) :: async
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    integer i, j, ierr

    call fill_halo_2d_r8(block, array, full_lon, full_lat, west_halo, east_halo, south_halo, north_halo)
    return

    i = merge(1, 2, full_lon)
    j = merge(1, 2, full_lat)

    if (merge(west_halo, .true., present(west_halo))) then
      if (async%send_req(west) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(west), MPI_STATUS_IGNORE, ierr)
      if (async%recv_req(west) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(west), MPI_STATUS_IGNORE, ierr)
      call MPI_ISEND(array, 1, block%halo(east)%send_type_2d(i,j), block%halo(east)%proc_id, 3, proc%comm, async%send_req(west), ierr)
      call MPI_IRECV(array, 1, block%halo(west)%recv_type_2d(i,j), block%halo(west)%proc_id, 3, proc%comm, async%recv_req(west), ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      if (async%send_req(east) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(east), MPI_STATUS_IGNORE, ierr)
      if (async%recv_req(east) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(east), MPI_STATUS_IGNORE, ierr)
      call MPI_ISEND(array, 1, block%halo(west)%send_type_2d(i,j), block%halo(west)%proc_id, 7, proc%comm, async%send_req(east), ierr)
      call MPI_IRECV(array, 1, block%halo(east)%recv_type_2d(i,j), block%halo(east)%proc_id, 7, proc%comm, async%recv_req(east), ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      if (async%send_req(south) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(south), MPI_STATUS_IGNORE, ierr)
      if (async%recv_req(south) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(south), MPI_STATUS_IGNORE, ierr)
      call MPI_ISEND(array, 1, block%halo(north)%send_type_2d(i,j), block%halo(north)%proc_id, 11, proc%comm, async%send_req(south), ierr)
      call MPI_IRECV(array, 1, block%halo(south)%recv_type_2d(i,j), block%halo(south)%proc_id, 11, proc%comm, async%recv_req(south), ierr)
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      if (async%send_req(north) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(north), MPI_STATUS_IGNORE, ierr)
      if (async%recv_req(north) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(north), MPI_STATUS_IGNORE, ierr)
      call MPI_ISEND(array, 1, block%halo(south)%send_type_2d(i,j), block%halo(south)%proc_id, 15, proc%comm, async%send_req(north), ierr)
      call MPI_IRECV(array, 1, block%halo(north)%recv_type_2d(i,j), block%halo(north)%proc_id, 15, proc%comm, async%recv_req(north), ierr)
    end if

  end subroutine fill_halo_2d_r8_async

  subroutine fill_halo_3d_r8(block, array, full_lon, full_lat, full_lev, west_halo, east_halo, south_halo, north_halo)

    type(block_type), intent(in) :: block
    real(8), intent(inout) :: array(:,:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    logical, intent(in), optional :: full_lev
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    integer status(MPI_STATUS_SIZE), i, j, k, ierr

    i = merge(1, 2, full_lon)
    j = merge(1, 2, full_lat)
    k = merge(1, 2, merge(full_lev, .true., present(full_lev)))

    if (merge(west_halo, .true., present(west_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(east)%send_type_3d(i,j,k), block%halo(east)%proc_id, 3, &
                        array, 1, block%halo(west)%recv_type_3d(i,j,k), block%halo(west)%proc_id, 3, &
                        proc%comm, status, ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(west)%send_type_3d(i,j,k), block%halo(west)%proc_id, 7, &
                        array, 1, block%halo(east)%recv_type_3d(i,j,k), block%halo(east)%proc_id, 7, &
                        proc%comm, status, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(north)%send_type_3d(i,j,k), block%halo(north)%proc_id, 11, &
                        array, 1, block%halo(south)%recv_type_3d(i,j,k), block%halo(south)%proc_id, 11, &
                        proc%comm, status, ierr)
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(south)%send_type_3d(i,j,k), block%halo(south)%proc_id, 15, &
                        array, 1, block%halo(north)%recv_type_3d(i,j,k), block%halo(north)%proc_id, 15, &
                        proc%comm, status, ierr)
    end if

  end subroutine fill_halo_3d_r8

  subroutine fill_halo_3d_r8_async(block, array, full_lon, full_lat, async, full_lev, west_halo, east_halo, south_halo, north_halo)

    type(block_type), intent(in) :: block
    real(8), intent(inout) :: array(:,:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    type(async_type), intent(inout) :: async
    logical, intent(in), optional :: full_lev
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    integer i, j, k, ierr

    call fill_halo_3d_r8(block, array, full_lon, full_lat, full_lev, west_halo, east_halo, south_halo, north_halo)
    return

    i = merge(1, 2, full_lon)
    j = merge(1, 2, full_lat)
    k = merge(1, 2, merge(full_lev, .true., present(full_lev)))

    if (merge(west_halo, .true., present(west_halo))) then
      if (async%send_req(west) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(west), MPI_STATUS_IGNORE, ierr)
      if (async%recv_req(west) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(west), MPI_STATUS_IGNORE, ierr)
      call MPI_ISEND(array, 1, block%halo(east)%send_type_3d(i,j,k), block%halo(east)%proc_id, 3, proc%comm, async%send_req(west), ierr)
      call MPI_IRECV(array, 1, block%halo(west)%recv_type_3d(i,j,k), block%halo(west)%proc_id, 3, proc%comm, async%recv_req(west), ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      if (async%send_req(east) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(east), MPI_STATUS_IGNORE, ierr)
      if (async%recv_req(east) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(east), MPI_STATUS_IGNORE, ierr)
      call MPI_ISEND(array, 1, block%halo(west)%send_type_3d(i,j,k), block%halo(west)%proc_id, 7, proc%comm, async%send_req(east), ierr)
      call MPI_IRECV(array, 1, block%halo(east)%recv_type_3d(i,j,k), block%halo(east)%proc_id, 7, proc%comm, async%recv_req(east), ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      if (async%send_req(south) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(south), MPI_STATUS_IGNORE, ierr)
      if (async%recv_req(south) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(south), MPI_STATUS_IGNORE, ierr)
      call MPI_ISEND(array, 1, block%halo(north)%send_type_3d(i,j,k), block%halo(north)%proc_id, 11, proc%comm, async%send_req(south), ierr)
      call MPI_IRECV(array, 1, block%halo(south)%recv_type_3d(i,j,k), block%halo(south)%proc_id, 11, proc%comm, async%recv_req(south), ierr)
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      if (async%send_req(north) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(north), MPI_STATUS_IGNORE, ierr)
      if (async%recv_req(north) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(north), MPI_STATUS_IGNORE, ierr)
      call MPI_ISEND(array, 1, block%halo(south)%send_type_3d(i,j,k), block%halo(south)%proc_id, 15, proc%comm, async%send_req(north), ierr)
      call MPI_IRECV(array, 1, block%halo(north)%recv_type_3d(i,j,k), block%halo(north)%proc_id, 15, proc%comm, async%recv_req(north), ierr)
    end if

  end subroutine fill_halo_3d_r8_async

  subroutine fill_halo_2d_r4_member(block, array, full_lon, full_lat, west_halo, east_halo, south_halo, north_halo)

    type(block_type), intent(in) :: block
    real(4), intent(inout) :: array(:,:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    integer status(MPI_STATUS_SIZE), i, j, p, ierr

    i = merge(1, 2, full_lon)
    j = merge(1, 2, full_lat)

    if (merge(west_halo, .true., present(west_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(east)%send_type_2d(i,j), block%halo(east)%proc_id, 3, &
                        array, 1, block%halo(west)%recv_type_2d(i,j), block%halo(west)%proc_id, 3, &
                        proc%comm, status, ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(west)%send_type_2d(i,j), block%halo(west)%proc_id, 7, &
                        array, 1, block%halo(east)%recv_type_2d(i,j), block%halo(east)%proc_id, 7, &
                        proc%comm, status, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(north)%send_type_2d(i,j), block%halo(north)%proc_id, 11, &
                        array, 1, block%halo(south)%recv_type_2d(i,j), block%halo(south)%proc_id, 11, &
                        proc%comm, status, ierr)
      do p = 5 , proc%ngb_num
        if (proc%ngb(p)%orient == south) then
          call MPI_RECV(array , 1 , block%halo(p)%recv_type_2d(i,j), block%halo(p)%proc_id, 11, proc%comm, status, ierr)
        end if
        if (proc%ngb(p)%orient == north) then
          call MPI_SEND(array , 1 , block%halo(p)%send_type_2d(i,j), block%halo(p)%proc_id, 11, proc%comm, ierr)
        end if
      end do
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(south)%send_type_2d(i,j), block%halo(south)%proc_id, 15, &
                        array, 1, block%halo(north)%recv_type_2d(i,j), block%halo(north)%proc_id, 15, &
                        proc%comm, status, ierr)
      do p = 5 , proc%ngb_num
        if (proc%ngb(p)%orient == north) then
          call MPI_RECV(array , 1 , block%halo(p)%recv_type_2d(i,j), block%halo(p)%proc_id, 15, proc%comm, status, ierr)
        end if
        if (proc%ngb(p)%orient == south) then
          call MPI_SEND(array , 1 , block%halo(p)%send_type_2d(i,j), block%halo(p)%proc_id, 15, proc%comm, ierr)
        end if
      end do
    end if

  end subroutine fill_halo_2d_r4_member

  subroutine fill_halo_2d_r8_member(block, array, full_lon, full_lat, west_halo, east_halo, south_halo, north_halo, south_ngb, north_ngb)

    type(block_type), intent(in) :: block
    real(8), intent(inout) :: array(:,:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo
    logical, intent(in), optional :: south_ngb
    logical, intent(in), optional :: north_ngb

    integer status(MPI_STATUS_SIZE), i, j, p, ierr
    logical south_ngb_in , north_ngb_in

    !xxx_ngb_in = true , need to interact with xxx_ngb
    south_ngb_in = merge(south_ngb , .true. , present(south_ngb))
    north_ngb_in = merge(north_ngb , .true. , present(north_ngb))

    i = merge(1, 2, full_lon)
    j = merge(1, 2, full_lat)

    if (merge(west_halo, .true., present(west_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(east)%send_type_2d(i,j), block%halo(east)%proc_id, 3, &
                        array, 1, block%halo(west)%recv_type_2d(i,j), block%halo(west)%proc_id, 3, &
                        proc%comm, status, ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(west)%send_type_2d(i,j), block%halo(west)%proc_id, 7, &
                        array, 1, block%halo(east)%recv_type_2d(i,j), block%halo(east)%proc_id, 7, &
                        proc%comm, status, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      if (south_ngb_in .and. north_ngb_in) then
        call MPI_SENDRECV(array, 1, block%halo(north)%send_type_2d(i,j), block%halo(north)%proc_id, 11, &
                          array, 1, block%halo(south)%recv_type_2d(i,j), block%halo(south)%proc_id, 11, &
                          proc%comm, status, ierr)
      else
        if (.not. (south_ngb_in .or. north_ngb_in)) then
          !nothing to do
        else if (.not. south_ngb_in) then
          call MPI_SEND(array, 1, block%halo(north)%send_type_2d(i,j), block%halo(north)%proc_id, 11,proc%comm, ierr)
        else if (.not. north_ngb_in) then
          call MPI_RECV(array, 1, block%halo(south)%recv_type_2d(i,j), block%halo(south)%proc_id, 11,proc%comm, status, ierr)
        end if
      end if


      do p = 5 , proc%ngb_num
        if (proc%ngb(p)%orient == south) then
          if (south_ngb_in) call MPI_RECV(array , 1 , block%halo(p)%recv_type_2d(i,j), block%halo(p)%proc_id, 11, proc%comm, status, ierr)
        end if
        if (proc%ngb(p)%orient == north) then
          if (north_ngb_in) call MPI_SEND(array , 1 , block%halo(p)%send_type_2d(i,j), block%halo(p)%proc_id, 11, proc%comm, ierr)
        end if
      end do
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      if (south_ngb_in .and. north_ngb_in) then
        call MPI_SENDRECV(array, 1, block%halo(south)%send_type_2d(i,j), block%halo(south)%proc_id, 15, &
                          array, 1, block%halo(north)%recv_type_2d(i,j), block%halo(north)%proc_id, 15, &
                          proc%comm, status, ierr)
      else
        if (.not. (south_ngb_in .or. north_ngb_in)) then
          !nothing to do
        else if (.not. south_ngb_in) then
          call MPI_RECV(array, 1, block%halo(north)%recv_type_2d(i,j), block%halo(north)%proc_id, 11,proc%comm, status, ierr)
        else if (.not. north_ngb_in) then
          call MPI_SEND(array, 1, block%halo(south)%send_type_2d(i,j), block%halo(south)%proc_id, 11,proc%comm, ierr)
        end if
      end if
      do p = 5 , proc%ngb_num
        if (proc%ngb(p)%orient == north) then
          if (north_ngb_in) call MPI_RECV(array , 1 , block%halo(p)%recv_type_2d(i,j), block%halo(p)%proc_id, 15, proc%comm, status, ierr)
        end if
        if (proc%ngb(p)%orient == south) then
          if (south_ngb_in) call MPI_SEND(array , 1 , block%halo(p)%send_type_2d(i,j), block%halo(p)%proc_id, 15, proc%comm, ierr)
        end if
      end do
    end if

  end subroutine fill_halo_2d_r8_member

  subroutine fill_halo_2d_r8_async_member(block, array, full_lon, full_lat, async, west_halo, east_halo, south_halo, north_halo, south_ngb, north_ngb)

    type(block_type), intent(in) :: block
    real(8), intent(inout) :: array(:,:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    type(async_type), intent(inout) :: async
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo
    logical, intent(in), optional :: south_ngb
    logical, intent(in), optional :: north_ngb

    integer i, j, ierr, p, status(MPI_STATUS_SIZE)
    logical south_ngb_in , north_ngb_in

    if (blocked_comm) then
      call fill_halo_2d_r8_member(block, array, full_lon, full_lat, west_halo, east_halo, south_halo, north_halo, south_ngb, north_ngb)
      return
    end if

    !xxx_ngb_in = true , need to interact with xxx_ngb
    south_ngb_in = merge(south_ngb , .true. , present(south_ngb))
    north_ngb_in = merge(north_ngb , .true. , present(north_ngb))

    i = merge(1, 2, full_lon)
    j = merge(1, 2, full_lat)
    
    if (merge(west_halo, .true., present(west_halo))) then
      if (async%send_req(east) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(east), MPI_STATUS_IGNORE, ierr)
      if (async%recv_req(west) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(west), MPI_STATUS_IGNORE, ierr)

      call MPI_ISEND(array, 1, block%halo_send(east)%mpi_type_2d(i,j) , block%halo_send(east)%proc_id , block%halo_send(east)%tag , &
                     proc%comm, async%send_req(east) , ierr)
      call MPI_IRECV(array, 1, block%halo_recv(west)%mpi_type_2d(i,j) , block%halo_recv(west)%proc_id , block%halo_recv(west)%tag , &
                     proc%comm, async%recv_req(west) , ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      if (async%send_req(west) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(west), MPI_STATUS_IGNORE, ierr)
      if (async%recv_req(east) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(east), MPI_STATUS_IGNORE, ierr)

      call MPI_ISEND(array, 1, block%halo_send(west)%mpi_type_2d(i,j) , block%halo_send(west)%proc_id , block%halo_send(west)%tag , &
                     proc%comm, async%send_req(west), ierr)
      call MPI_IRECV(array, 1, block%halo_recv(east)%mpi_type_2d(i,j) , block%halo_recv(east)%proc_id , block%halo_recv(east)%tag , &
                     proc%comm, async%recv_req(east), ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      do p = proc%ngb_north_send_beg , proc%ngb_south_send_beg - 1
        if (async%send_req(p) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(p), MPI_STATUS_IGNORE, ierr)
        call MPI_ISEND(array, 1, block%halo_send(p)%mpi_type_2d(i,j) , block%halo_send(p)%proc_id , block%halo_send(p)%tag , &
                       proc%comm, async%send_req(p), ierr)
      end do

      do p = proc%ngb_south_recv_beg , proc%ngb_recv_num
        if (async%recv_req(p) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(p), MPI_STATUS_IGNORE, ierr)
        call MPI_IRECV(array, 1, block%halo_recv(p)%mpi_type_2d(i,j) , block%halo_recv(p)%proc_id , block%halo_recv(p)%tag , &
                       proc%comm, async%recv_req(p), ierr)
      end do
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      do p = proc%ngb_south_send_beg , proc%ngb_send_num
        if (async%send_req(p) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(p), MPI_STATUS_IGNORE, ierr)
        call MPI_ISEND(array, 1, block%halo_send(p)%mpi_type_2d(i,j) , block%halo_send(p)%proc_id , block%halo_send(p)%tag , &
                       proc%comm, async%send_req(p), ierr)
      end do

      do p = proc%ngb_north_recv_beg , proc%ngb_south_recv_beg - 1
        if (async%recv_req(p) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(p), MPI_STATUS_IGNORE, ierr)
        call MPI_IRECV(array, 1, block%halo_recv(p)%mpi_type_2d(i,j) , block%halo_recv(p)%proc_id , block%halo_recv(p)%tag , &
                       proc%comm, async%recv_req(p), ierr)
      end do
    end if


    ! if (merge(west_halo, .true., present(west_halo))) then
    !   call MPI_SENDRECV(array, 1, block%halo(east)%send_type_2d(i,j), block%halo(east)%proc_id, 3, &
    !                     array, 1, block%halo(west)%recv_type_2d(i,j), block%halo(west)%proc_id, 3, &
    !                     proc%comm, status, ierr)
    ! end if

    ! if (merge(east_halo, .true., present(east_halo))) then
    !   call MPI_SENDRECV(array, 1, block%halo(west)%send_type_2d(i,j), block%halo(west)%proc_id, 7, &
    !                     array, 1, block%halo(east)%recv_type_2d(i,j), block%halo(east)%proc_id, 7, &
    !                     proc%comm, status, ierr)
    ! end if

    ! if (merge(south_halo, .true., present(south_halo))) then
    !   if (south_ngb_in) then
    !     if (async%recv_req(south) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(south), MPI_STATUS_IGNORE, ierr)
    !     call MPI_IRECV(array, 1, block%halo(south)%recv_type_2d(i,j), block%halo(south)%proc_id, 11, proc%comm, async%recv_req(south), ierr)
    !   end if
    !   if (north_ngb_in) then
    !     if (async%send_req(north) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(north), MPI_STATUS_IGNORE, ierr)
    !     call MPI_ISEND(array, 1, block%halo(north)%send_type_2d(i,j), block%halo(north)%proc_id, 11, proc%comm, async%send_req(north), ierr)
    !   end if
    !   ! 多邻居
    !   do p = 5 , proc%ngb_num
    !     if (proc%ngb(p)%orient == south) then
    !       if (south_ngb_in) then
    !         if (async%recv_req(p) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(p), MPI_STATUS_IGNORE, ierr)
    !         call MPI_IRECV(array , 1 , block%halo(p)%recv_type_2d(i,j), block%halo(p)%proc_id, 11, proc%comm, async%recv_req(p), ierr)
    !       end if
    !     end if
    !     if (proc%ngb(p)%orient == north) then
    !       if (north_ngb_in) then
    !         if (async%send_req(p) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(p), MPI_STATUS_IGNORE, ierr)
    !         call MPI_ISEND(array , 1 , block%halo(p)%send_type_2d(i,j), block%halo(p)%proc_id, 11, proc%comm, async%send_req(p), ierr)
    !       end if
    !     end if
    !   end do
    ! end if

    ! if (merge(north_halo, .true., present(north_halo))) then
    !   if (north_ngb_in) then
    !     if (async%recv_req(north) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(north), MPI_STATUS_IGNORE, ierr)
    !     call MPI_IRECV(array, 1, block%halo(north)%recv_type_2d(i,j), block%halo(north)%proc_id, 15, proc%comm, async%recv_req(north), ierr)
    !   end if
    !   if (south_ngb_in) then
    !     if (async%send_req(south) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(south), MPI_STATUS_IGNORE, ierr)
    !     call MPI_ISEND(array, 1, block%halo(south)%send_type_2d(i,j), block%halo(south)%proc_id, 15, proc%comm, async%send_req(south), ierr)
    !   end if
    !   ! 多邻居
    !   do p = 5 , proc%ngb_num
    !     if (proc%ngb(p)%orient == north) then
    !       if (north_ngb_in) then
    !         if (async%recv_req(p) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(p), MPI_STATUS_IGNORE, ierr)
    !         call MPI_IRECV(array , 1 , block%halo(p)%recv_type_2d(i,j), block%halo(p)%proc_id, 15, proc%comm, async%recv_req(p), ierr)
    !       end if
    !     end if
    !     if (proc%ngb(p)%orient == south) then
    !       if (south_ngb_in) then
    !         if (async%send_req(p) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(p), MPI_STATUS_IGNORE, ierr)
    !         call MPI_ISEND(array , 1 , block%halo(p)%send_type_2d(i,j), block%halo(p)%proc_id, 15, proc%comm, async%send_req(p), ierr)
    !       end if
    !     end if
    !   end do
    ! end if

  end subroutine fill_halo_2d_r8_async_member

  subroutine fill_halo_3d_r4_member(block, array, full_lon, full_lat, full_lev, west_halo, east_halo, south_halo, north_halo)

    type(block_type), intent(in) :: block
    real(4), intent(inout) :: array(:,:,:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    logical, intent(in), optional :: full_lev
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    integer status(MPI_STATUS_SIZE), i, j, k, p, ierr

    i = merge(1, 2, full_lon)
    j = merge(1, 2, full_lat)
    k = merge(1, 2, merge(full_lev, .true., present(full_lev)))

    if (merge(west_halo, .true., present(west_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(east)%send_type_3d(i,j,k), block%halo(east)%proc_id, 3, &
                        array, 1, block%halo(west)%recv_type_3d(i,j,k), block%halo(west)%proc_id, 3, &
                        proc%comm, status, ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(west)%send_type_3d(i,j,k), block%halo(west)%proc_id, 7, &
                        array, 1, block%halo(east)%recv_type_3d(i,j,k), block%halo(east)%proc_id, 7, &
                        proc%comm, status, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(north)%send_type_3d(i,j,k), block%halo(north)%proc_id, 11, &
                        array, 1, block%halo(south)%recv_type_3d(i,j,k), block%halo(south)%proc_id, 11, &
                        proc%comm, status, ierr)
      do p = 5 , proc%ngb_num
        if (proc%ngb(p)%orient == south) then
          call MPI_RECV(array , 1 , block%halo(p)%recv_type_3d(i,j,k), block%halo(p)%proc_id, 11, proc%comm, status, ierr)
        end if
        if (proc%ngb(p)%orient == north) then
          call MPI_SEND(array , 1 , block%halo(p)%send_type_3d(i,j,k), block%halo(p)%proc_id, 11, proc%comm, ierr)
        end if
      end do
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(south)%send_type_3d(i,j,k), block%halo(south)%proc_id, 15, &
                        array, 1, block%halo(north)%recv_type_3d(i,j,k), block%halo(north)%proc_id, 15, &
                        proc%comm, status, ierr)
      do p = 5 , proc%ngb_num

        if (proc%ngb(p)%orient == north) then
          call MPI_RECV(array , 1 , block%halo(p)%recv_type_3d(i,j,k), block%halo(p)%proc_id, 15, proc%comm, status, ierr)
        end if
        if (proc%ngb(p)%orient == south) then
          call MPI_SEND(array , 1 , block%halo(p)%send_type_3d(i,j,k), block%halo(p)%proc_id, 15, proc%comm, ierr)
        end if
      end do
    end if

  end subroutine fill_halo_3d_r4_member

  subroutine fill_halo_3d_r8_member(block, array, full_lon, full_lat, full_lev, west_halo, east_halo, south_halo, north_halo, south_ngb, north_ngb)

    type(block_type), intent(in) :: block
    real(8), intent(inout) :: array(:,:,:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    logical, intent(in), optional :: full_lev
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo
    logical, intent(in), optional :: south_ngb
    logical, intent(in), optional :: north_ngb

    integer status(MPI_STATUS_SIZE), i, j, k, p, ierr
    logical south_ngb_in , north_ngb_in

    !xxx_ngb_in = true , need to interact with xxx_ngb
    south_ngb_in = merge(south_ngb , .true. , present(south_ngb))
    north_ngb_in = merge(north_ngb , .true. , present(north_ngb))

    i = merge(1, 2, full_lon)
    j = merge(1, 2, full_lat)
    k = merge(1, 2, merge(full_lev, .true., present(full_lev)))

    if (merge(west_halo, .true., present(west_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(east)%send_type_3d(i,j,k), block%halo(east)%proc_id, 3, &
                        array, 1, block%halo(west)%recv_type_3d(i,j,k), block%halo(west)%proc_id, 3, &
                        proc%comm, status, ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(west)%send_type_3d(i,j,k), block%halo(west)%proc_id, 7, &
                        array, 1, block%halo(east)%recv_type_3d(i,j,k), block%halo(east)%proc_id, 7, &
                        proc%comm, status, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      if (south_ngb_in .and. north_ngb_in) then
        call MPI_SENDRECV(array, 1, block%halo(north)%send_type_3d(i,j,k), block%halo(north)%proc_id, 11, &
                          array, 1, block%halo(south)%recv_type_3d(i,j,k), block%halo(south)%proc_id, 11, &
                          proc%comm, status, ierr)
      else
        if (.not. (south_ngb_in .or. north_ngb_in)) then
          !nothing to do
        else if (south_ngb_in) then
          call MPI_RECV(array, 1, block%halo(south)%recv_type_3d(i,j,k), block%halo(south)%proc_id, 11, proc%comm, status, ierr)
        else if (north_ngb_in) then
          call MPI_SEND(array, 1, block%halo(north)%send_type_3d(i,j,k), block%halo(north)%proc_id, 11, proc%comm, ierr)
        end if
      end if

      do p = 5 , proc%ngb_num
        if (proc%ngb(p)%orient == south .and. south_ngb_in) then
          call MPI_RECV(array , 1 , block%halo(p)%recv_type_3d(i,j,k), block%halo(p)%proc_id, 11, proc%comm, status, ierr)
        end if
        if (proc%ngb(p)%orient == north .and. north_ngb_in) then
          call MPI_SEND(array , 1 , block%halo(p)%send_type_3d(i,j,k), block%halo(p)%proc_id, 11, proc%comm, ierr)
        end if
      end do
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      if (south_ngb_in .and. north_ngb_in) then
        call MPI_SENDRECV(array, 1, block%halo(south)%send_type_3d(i,j,k), block%halo(south)%proc_id, 15, &
                          array, 1, block%halo(north)%recv_type_3d(i,j,k), block%halo(north)%proc_id, 15, &
                          proc%comm, status, ierr)
      else
        if (.not. (south_ngb_in .or. north_ngb_in)) then
          !nothing to do
        else if (south_ngb_in) then
          call MPI_SEND(array, 1, block%halo(south)%send_type_3d(i,j,k), block%halo(south)%proc_id, 15,proc%comm, ierr)
        else if (north_ngb_in) then
          call MPI_RECV(array, 1, block%halo(north)%recv_type_3d(i,j,k), block%halo(north)%proc_id, 15,proc%comm, status, ierr)
        end if
      end if
      do p = 5 , proc%ngb_num

        if (proc%ngb(p)%orient == north .and. north_ngb_in) then
          call MPI_RECV(array , 1 , block%halo(p)%recv_type_3d(i,j,k), block%halo(p)%proc_id, 15, proc%comm, status, ierr)
        end if
        if (proc%ngb(p)%orient == south .and. south_ngb_in) then
          call MPI_SEND(array , 1 , block%halo(p)%send_type_3d(i,j,k), block%halo(p)%proc_id, 15, proc%comm, ierr)
        end if
      end do
    end if

  end subroutine fill_halo_3d_r8_member

  subroutine fill_halo_3d_r8_async_member(block, array, full_lon, full_lat, async, full_lev, west_halo, east_halo, south_halo, north_halo, south_ngb, north_ngb)

    type(block_type), intent(in) :: block
    real(8), intent(inout) :: array(:,:,:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    type(async_type), intent(inout) :: async
    logical, intent(in), optional :: full_lev
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo
    logical, intent(in), optional :: south_ngb
    logical, intent(in), optional :: north_ngb


    integer status(MPI_STATUS_SIZE), i, j, k, ierr, p
    logical south_ngb_in , north_ngb_in

    if (blocked_comm) then
      call fill_halo_3d_r8_member(block, array, full_lon, full_lat, full_lev, west_halo, east_halo, south_halo, north_halo, south_ngb, north_ngb)
      return
    end if

    ! !xxx_ngb_in = true , need to interact with xxx_ngb
    ! south_ngb_in = merge(south_ngb , .true. , present(south_ngb))
    ! north_ngb_in = merge(north_ngb , .true. , present(north_ngb))

    i = merge(1, 2, full_lon)
    j = merge(1, 2, full_lat)
    k = merge(1, 2, merge(full_lev, .true., present(full_lev)))

    if (merge(west_halo, .true., present(west_halo))) then
      if (async%send_req(east) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(east), MPI_STATUS_IGNORE, ierr)
      if (async%recv_req(west) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(west), MPI_STATUS_IGNORE, ierr)

      call MPI_ISEND(array, 1, block%halo_send(east)%mpi_type_3d(i,j,k) , block%halo_send(east)%proc_id , block%halo_send(east)%tag , &
                     proc%comm, async%send_req(east), ierr)
      call MPI_IRECV(array, 1, block%halo_recv(west)%mpi_type_3d(i,j,k) , block%halo_recv(west)%proc_id , block%halo_recv(west)%tag , &
                     proc%comm, async%recv_req(west), ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      if (async%send_req(west) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(west), MPI_STATUS_IGNORE, ierr)
      if (async%recv_req(east) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(east), MPI_STATUS_IGNORE, ierr)

      call MPI_ISEND(array, 1, block%halo_send(west)%mpi_type_3d(i,j,k) , block%halo_send(west)%proc_id , block%halo_send(west)%tag , &
                     proc%comm, async%send_req(west), ierr)
      call MPI_IRECV(array, 1, block%halo_recv(east)%mpi_type_3d(i,j,k) , block%halo_recv(east)%proc_id , block%halo_recv(east)%tag , &
                     proc%comm, async%recv_req(east), ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      do p = proc%ngb_north_send_beg , proc%ngb_south_send_beg - 1
        if (async%send_req(p) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(p), MPI_STATUS_IGNORE, ierr)
        call MPI_ISEND(array, 1, block%halo_send(p)%mpi_type_3d(i,j,k) , block%halo_send(p)%proc_id , block%halo_send(p)%tag , &
                       proc%comm, async%send_req(p), ierr)
      end do

      do p = proc%ngb_south_recv_beg , proc%ngb_recv_num
        if (async%recv_req(p) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(p), MPI_STATUS_IGNORE, ierr)
        call MPI_IRECV(array, 1, block%halo_recv(p)%mpi_type_3d(i,j,k) , block%halo_recv(p)%proc_id , block%halo_recv(p)%tag , &
                       proc%comm, async%recv_req(p), ierr)
      end do
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      do p = proc%ngb_south_send_beg , proc%ngb_send_num
        if (async%send_req(p) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(p), MPI_STATUS_IGNORE, ierr)
        call MPI_ISEND(array, 1, block%halo_send(p)%mpi_type_3d(i,j,k) , block%halo_send(p)%proc_id , block%halo_send(p)%tag , &
                       proc%comm, async%send_req(p), ierr)
      end do

      do p = proc%ngb_north_recv_beg , proc%ngb_south_recv_beg - 1
        if (async%recv_req(p) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(p), MPI_STATUS_IGNORE, ierr)
        call MPI_IRECV(array, 1, block%halo_recv(p)%mpi_type_3d(i,j,k) , block%halo_recv(p)%proc_id , block%halo_recv(p)%tag , &
                       proc%comm, async%recv_req(p), ierr)
      end do
    end if

    !==========================================================================================

    ! if (merge(west_halo, .true., present(west_halo))) then
    !   call MPI_SENDRECV(array, 1, block%halo(east)%send_type_3d(i,j,k), block%halo(east)%proc_id, 3, &
    !                     array, 1, block%halo(west)%recv_type_3d(i,j,k), block%halo(west)%proc_id, 3, &
    !                     proc%comm, status, ierr)
    ! end if

    ! if (merge(east_halo, .true., present(east_halo))) then
    !   call MPI_SENDRECV(array, 1, block%halo(west)%send_type_3d(i,j,k), block%halo(west)%proc_id, 7, &
    !                     array, 1, block%halo(east)%recv_type_3d(i,j,k), block%halo(east)%proc_id, 7, &
    !                     proc%comm, status, ierr)
    ! end if

    ! if (merge(south_halo, .true., present(south_halo))) then
    !   if (south_ngb_in) then
    !     if (async%recv_req(south) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(south), MPI_STATUS_IGNORE, ierr)
    !     call MPI_IRECV(array, 1, block%halo(south)%recv_type_3d(i,j,k), block%halo(south)%proc_id, 11, proc%comm, async%recv_req(south), ierr)
    !   end if

    !   if (north_ngb_in) then
    !     if (async%send_req(north) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(north), MPI_STATUS_IGNORE, ierr)
    !     call MPI_ISEND(array, 1, block%halo(north)%send_type_3d(i,j,k), block%halo(north)%proc_id, 11, proc%comm, async%send_req(north), ierr)
    !   end if
    !   ! 多邻居
    !   do p = 5 , proc%ngb_num
    !     if (proc%ngb(p)%orient == south .and. south_ngb_in) then
    !       if (async%recv_req(p) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(p), MPI_STATUS_IGNORE, ierr)
    !       call MPI_IRECV(array , 1 , block%halo(p)%recv_type_3d(i,j,k), block%halo(p)%proc_id, 11, proc%comm, async%recv_req(p), ierr)
    !     end if
    !     if (proc%ngb(p)%orient == north .and. north_ngb_in) then
    !       if (async%send_req(p) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(p), MPI_STATUS_IGNORE, ierr)
    !       call MPI_ISEND(array , 1 , block%halo(p)%send_type_3d(i,j,k), block%halo(p)%proc_id, 11, proc%comm, async%send_req(p), ierr)
    !     end if
    !   end do
    ! end if

    ! if (merge(north_halo, .true., present(north_halo))) then
    !   if (south_ngb_in) then
    !     if (async%send_req(south) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(south), MPI_STATUS_IGNORE, ierr)
    !     call MPI_ISEND(array, 1, block%halo(south)%send_type_3d(i,j,k), block%halo(south)%proc_id, 15, proc%comm, async%send_req(south), ierr)
    !   end if

    !   if (north_ngb_in) then
    !     if (async%recv_req(north) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(north), MPI_STATUS_IGNORE, ierr)
    !     call MPI_IRECV(array, 1, block%halo(north)%recv_type_3d(i,j,k), block%halo(north)%proc_id, 15, proc%comm, async%recv_req(north), ierr)
    !   end if

    !     ! 多邻居
    !   do p = 5 , proc%ngb_num
    !     if (proc%ngb(p)%orient == north .and. north_ngb_in) then
    !       if (async%recv_req(p) /= MPI_REQUEST_NULL) call MPI_WAIT(async%recv_req(p), MPI_STATUS_IGNORE, ierr)
    !       call MPI_IRECV(array , 1 , block%halo(p)%recv_type_3d(i,j,k), block%halo(p)%proc_id, 15, proc%comm, async%recv_req(p), ierr)
    !     end if
    !     if (proc%ngb(p)%orient == south .and. south_ngb_in) then
    !       if (async%send_req(p) /= MPI_REQUEST_NULL) call MPI_WAIT(async%send_req(p), MPI_STATUS_IGNORE, ierr)
    !       call MPI_ISEND(array , 1 , block%halo(p)%send_type_3d(i,j,k), block%halo(p)%proc_id, 15, proc%comm, async%send_req(p), ierr)
    !     end if
    !   end do
    ! end if

  end subroutine fill_halo_3d_r8_async_member

  subroutine zero_halo_1d_r8(block, array, west_halo, east_halo)

    type(block_type), intent(in) :: block
    real(8), intent(inout) :: array(:)
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer i1, i2

    if (merge(west_halo, .false., present(west_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! | i1  | i2  |     |     |     |     |     |     |     |     |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !    |     |
      !    1  1 + w - 1
      i1 = 1
      i2 = 1 + block%mesh%lon_halo_width - 1
      array(i1:i2) = 0.0d0
    end if

    if (merge(east_halo, .false., present(east_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! |     |     |     |     |     |     |     |     | i1  | i2  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !                                                    |     |
      !                                                n - w + 1 n
      i1 = size(array) - block%mesh%lon_halo_width + 1
      i2 = size(array)
      array(i1:i2) = 0.0d0
    end if

  end subroutine zero_halo_1d_r8

  subroutine zero_halo_1d_r4_member(block, array, west_halo, east_halo)

    type(block_type), intent(in) :: block
    real(4), intent(inout) :: array(:,:)
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer i1, i2

    if (merge(west_halo, .false., present(west_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! | i1  | i2  |     |     |     |     |     |     |     |     |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !    |     |
      !    1  1 + w - 1
      i1 = 1
      i2 = 1 + block%mesh%lon_halo_width - 1
      array(:,i1:i2) = 0.0d0
    end if

    if (merge(east_halo, .false., present(east_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! |     |     |     |     |     |     |     |     | i1  | i2  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !                                                    |     |
      !                                                n - w + 1 n
      i1 = size(array,2) - block%mesh%lon_halo_width + 1
      i2 = size(array,2)
      array(:,i1:i2) = 0.0d0
    end if

  end subroutine zero_halo_1d_r4_member

  subroutine zero_halo_1d_r8_member(block, array, west_halo, east_halo)

    type(block_type), intent(in) :: block
    real(8), intent(inout) :: array(:,:)
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer i1, i2

    if (merge(west_halo, .false., present(west_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! | i1  | i2  |     |     |     |     |     |     |     |     |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !    |     |
      !    1  1 + w - 1
      i1 = 1
      i2 = 1 + block%mesh%lon_halo_width - 1
      array(:,i1:i2) = 0.0d0
    end if

    if (merge(east_halo, .false., present(east_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! |     |     |     |     |     |     |     |     | i1  | i2  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !                                                    |     |
      !                                                n - w + 1 n
      i1 = size(array,2) - block%mesh%lon_halo_width + 1
      i2 = size(array,2)
      array(:,i1:i2) = 0.0d0
    end if

  end subroutine zero_halo_1d_r8_member

  ! subroutine overlay_inner_halo(block, array, west_halo, east_halo)

  !   type(block_type), intent(in) :: block
  !   real(8), intent(inout) :: array(:)
  !   logical, intent(in), optional :: west_halo
  !   logical, intent(in), optional :: east_halo

  !   integer i1, i2, i3, i4
  !   integer status(MPI_STATUS_SIZE), ierr
  !   real(8) buffer(block%mesh%lon_halo_width)

  !   if (merge(west_halo, .false., present(west_halo))) then
  !     !   west halo |           |                       | east_halo
  !     !  ___________|___________|_______________________|___________
  !     ! |     |     |     |     |     |     |     |     |     |     |
  !     ! |     |     | i3  | i4  |     |     |     |     | i1  | i2  |
  !     ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
  !     !                |                                   |
  !     !              1 + w                             n - w + 1
  !     i1 = size(array) - block%mesh%lon_halo_width + 1
  !     i2 = i1 + block%mesh%lon_halo_width - 1
  !     i3 =  1 + block%mesh%lon_halo_width
  !     i4 = i3 + block%mesh%lon_halo_width - 1
  !     call MPI_SENDRECV(array(i1:i2), block%mesh%lon_halo_width, block%halo(east)%dtype, block%halo(east)%proc_id, 19, &
  !                       buffer      , block%mesh%lon_halo_width, block%halo(west)%dtype, block%halo(west)%proc_id, 19, &
  !                       proc%comm, status, ierr)
  !     array(i3:i4) = array(i3:i4) + buffer
  !   end if

  ! end subroutine overlay_inner_halo

  ! subroutine overlay_inner_halo_r4_member(block, array, west_halo, east_halo, rf)

  !   type(block_type), intent(in) :: block
  !   real(4), intent(inout) :: array(:,:)
  !   logical, intent(in), optional :: west_halo
  !   logical, intent(in), optional :: east_halo
  !   integer, intent(in), optional :: rf ! reduce_factor

  !   integer i1, i2, i3, i4
  !   integer status(MPI_STATUS_SIZE), ierr
  !   real(4) buffer(member_num , block%mesh%lon_halo_width)
  !   integer nm

  !   if (merge(west_halo, .false., present(west_halo))) then
  !     !   west halo |           |                       | east_halo
  !     !  ___________|___________|_______________________|___________
  !     ! |     |     |     |     |     |     |     |     |     |     |
  !     ! |     |     | i3  | i4  |     |     |     |     | i1  | i2  |
  !     ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
  !     !                |                                   |
  !     !              1 + w                             n - w + 1
      
  !     nm = size(array,1) 

  !     i1 = size(array,2) - block%mesh%lon_halo_width + 1
  !     i2 = i1 + merge(rf, block%mesh%lon_halo_width, present(rf)) - 2
  !     ! i2 = i1 + block%mesh%lon_halo_width - 1
  !     i3 =  1 + block%mesh%lon_halo_width
  !     ! i4 = i3 + block%mesh%lon_halo_width - 1
  !     i4 = i3 + merge(rf, block%mesh%lon_halo_width, present(rf)) - 2
  !     call MPI_SENDRECV(array(:,i1:i2), block%mesh%lon_halo_width * nm, block%halo_send(east)%dtype, block%halo_send(east)%proc_id, 19, &
  !                       buffer        , block%mesh%lon_halo_width * nm, block%halo_recv(west)%dtype, block%halo_recv(west)%proc_id, 19, &
  !                       proc%comm, status, ierr)
  !     array(:,i3:i4) = array(:,i3:i4) + buffer
  !   end if

  ! end subroutine overlay_inner_halo_r4_member

  ! subroutine overlay_inner_halo_r8_member(block, array, west_halo, east_halo, rf)

  !   type(block_type), intent(in) :: block
  !   real(8), intent(inout) :: array(:,:)
  !   logical, intent(in), optional :: west_halo
  !   logical, intent(in), optional :: east_halo
  !   integer, intent(in), optional :: rf ! reduce_factor

  !   integer i1, i2, i3, i4
  !   integer status(MPI_STATUS_SIZE), ierr
  !   real(8) buffer(member_num , block%mesh%lon_halo_width)
  !   integer nm

  !   if (merge(west_halo, .false., present(west_halo))) then
  !     !   west halo |           |                       | east_halo
  !     !  ___________|___________|_______________________|___________
  !     ! |     |     |     |     |     |     |     |     |     |     |
  !     ! |     |     | i3  | i4  |     |     |     |     | i1  | i2  |
  !     ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
  !     !                |                                   |
  !     !              1 + w                             n - w + 1
      
  !     nm = size(array,1) 

  !     i1 = size(array,2) - block%mesh%lon_halo_width + 1
  !     i2 = i1 + merge(rf, block%mesh%lon_halo_width, present(rf)) - 2
  !     ! i2 = i1 + block%mesh%lon_halo_width - 1
  !     i3 =  1 + block%mesh%lon_halo_width
  !     ! i4 = i3 + block%mesh%lon_halo_width - 1
  !     i4 = i3 + merge(rf, block%mesh%lon_halo_width, present(rf)) - 2
  !     call MPI_SENDRECV(array(:,i1:i2), block%mesh%lon_halo_width * nm, block%halo_send(east)%dtype, block%halo_send(east)%proc_id, 19, &
  !                       buffer        , block%mesh%lon_halo_width * nm, block%halo_recv(west)%dtype, block%halo_recv(west)%proc_id, 19, &
  !                       proc%comm, status, ierr)
  !     array(:,i3:i4) = array(:,i3:i4) + buffer
  !   end if

  ! end subroutine overlay_inner_halo_r8_member

  ! subroutine overlay_inner_halo_chunkwise_r4_member(block, reduced_chunk, array, full_lat, full_lev)

  !   type(block_type), intent(in) :: block
  !   type(reduced_chunk_type), intent(in) :: reduced_chunk
  !   real(4), intent(inout) :: array(:,:,:,:)
  !   logical, intent(in) :: full_lat
  !   logical, intent(in) :: full_lev

  !   integer j, k, ibeg, iend, jbeg, jend, kbeg, kend, status(MPI_STATUS_SIZE), ierr
    
  !   j = merge(1, 2, full_lat)
  !   k = merge(1, 2, full_lev)

  !   call MPI_SENDRECV(array                          , 1, reduced_chunk%inner_halo_send_type(j,k), block%halo_send(east)%proc_id, 19, &
  !                     reduced_chunk%inner_halo_buffer, 1, reduced_chunk%inner_halo_recv_type(j,k), block%halo_recv(west)%proc_id, 19, &
  !                     proc%comm, status, ierr)
    
  !   ibeg = 1 + block%mesh%lon_halo_width
  !   iend = ibeg + reduced_chunk%reduce_factor - 2
  !   jbeg = merge(reduced_chunk%full_lat_ibeg, reduced_chunk%half_lat_ibeg, full_lat) &
  !         - merge(   block%mesh%full_lat_lb  ,    block%mesh%half_lat_lb  , full_lat) + 1  ! 这里是Fortran的从1开始，需要+1
  !   jend = jbeg + merge(reduced_chunk%num_full_lat, reduced_chunk%num_half_lat, full_lat) - 1
  !   kbeg = 1
  !   kend = merge(reduced_chunk%num_full_lev, reduced_chunk%num_half_lev, full_lev)

  !   array(:, ibeg:iend, jbeg:jend, :) = array(:, ibeg:iend, jbeg:jend, :)&
  !           + reduced_chunk%inner_halo_buffer(:, :, 1:(jend-jbeg+1), 1:(kend-kbeg+1) )

  ! end subroutine overlay_inner_halo_chunkwise_r4_member

  ! subroutine overlay_inner_halo_chunkwise_r8_member(block, reduced_chunk, array, full_lat, full_lev, debug)

  !   type(block_type), intent(in) :: block
  !   type(reduced_chunk_type), intent(inout) :: reduced_chunk
  !   real(8), intent(inout) :: array(:,:,:,:)
  !   logical, intent(in) :: full_lat
  !   logical, intent(in) :: full_lev
  !   logical, intent(in), optional :: debug

  !   integer j, k, ibeg, iend, jbeg, jend, kbeg, kend, status(MPI_STATUS_SIZE), ierr
    
  !   j = merge(1, 2, full_lat)
  !   k = merge(1, 2, full_lev)

  !   ! write(*,*), "oly", proc%id, &
  !   !             size(reduced_chunk%inner_halo_buffer, 1), size(reduced_chunk%inner_halo_buffer, 2), &
  !   !             size(reduced_chunk%inner_halo_buffer, 3), size(reduced_chunk%inner_halo_buffer, 4)
    
  !   call MPI_SENDRECV(array                          , 1, reduced_chunk%inner_halo_send_type(j,k), block%halo_send(east)%proc_id, 19, &
  !                     reduced_chunk%inner_halo_buffer, 1, reduced_chunk%inner_halo_recv_type(j,k), block%halo_recv(west)%proc_id, 19, &
  !                     proc%comm, status, ierr)

  !   ! 要注意这里并没有关于array的dimension的信息（各维度的大小、起始下标等）
    
  !   ibeg = 1 + block%mesh%lon_halo_width
  !   iend = ibeg + reduced_chunk%reduce_factor - 2
  !   jbeg = merge(reduced_chunk%full_lat_ibeg, reduced_chunk%half_lat_ibeg, full_lat) &
  !        - merge(   block%mesh%full_lat_lb  ,    block%mesh%half_lat_lb  , full_lat) + 1  ! 这里是Fortran的从1开始，需要+1
  !   jend = jbeg + merge(reduced_chunk%num_full_lat, reduced_chunk%num_half_lat, full_lat) - 1
  !   kbeg = 1
  !   kend = merge(reduced_chunk%num_full_lev, reduced_chunk%num_half_lev, full_lev)

  !   ! if (merge(debug, .false., present(debug))) then
  !   !   write(*,*), proc%id, ibeg, iend, jbeg, jend, kbeg, kend
  !   ! end if

  !   array(:, ibeg:iend, jbeg:jend, :) = array(:, ibeg:iend, jbeg:jend, :)&
  !           + reduced_chunk%inner_halo_buffer(:, :, 1:(jend-jbeg+1), 1:(kend-kbeg+1) )

  ! end subroutine overlay_inner_halo_chunkwise_r8_member

  subroutine zonal_sum_0d_r8(zonal_circle, work, value)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: work(:)
    real(8), intent(out) :: value

#ifdef ENSURE_ORDER
    real(8) allvalue(global_mesh%num_full_lon)
#endif
    integer ierr

#ifdef ENSURE_ORDER
    if (zonal_circle%np == 1) then
      value = sum(work)
    else
      call gather_zonal_array(zonal_circle, work, allvalue)
      if (zonal_circle%id == 0) value = sum(allvalue)
      call MPI_BCAST(value, 1, MPI_DOUBLE, 0, zonal_circle%comm, ierr)
    end if
#else
    call MPI_ALLREDUCE(sum(work), value, 1, MPI_DOUBLE, MPI_SUM, zonal_circle%comm, ierr)
#endif

  end subroutine zonal_sum_0d_r8

  subroutine zonal_sum_1d_r8(zonal_circle, work, value)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: work(:,:)
    real(8), intent(out) :: value(:)

#ifdef ENSURE_ORDER
    real(8) allvalue(global_mesh%num_full_lon,size(value))
#endif
    integer ierr

#ifdef ENSURE_ORDER
    if (zonal_circle%np == 1) then
      value = sum(work, dim=1)
    else
      call gather_zonal_array(zonal_circle, work, allvalue)
      if (zonal_circle%id == 0) value = sum(allvalue, dim=1)
      call MPI_BCAST(value, size(value), MPI_DOUBLE, 0, zonal_circle%comm, ierr)
    end if
#else
    call MPI_ALLREDUCE(sum(work, dim=1), value, size(value), MPI_DOUBLE, MPI_SUM, zonal_circle%comm, ierr)
#endif

  end subroutine zonal_sum_1d_r8

  subroutine zonal_sum_0d_r4_member(zonal_circle, work, value)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(4), intent(in) :: work(:,:)
    real(4), intent(out) :: value(:)

#ifdef ENSURE_ORDER
    real(4) allvalue(member_num , global_mesh%num_full_lon)
#endif
    integer ierr

#ifdef ENSURE_ORDER
    if (zonal_circle%np == 1) then
      value(:) = sum(work,2)
    else
      call gather_zonal_array_member(zonal_circle, work, allvalue)
      if (zonal_circle%id == 0) value(:) = sum(allvalue,2)
      call MPI_BCAST(value, size(value), MPI_REAL, 0, zonal_circle%comm, ierr)
    end if
#else
    call MPI_ALLREDUCE(sum(work,2), value, size(value), MPI_REAL, MPI_SUM, zonal_circle%comm, ierr)
#endif

  end subroutine zonal_sum_0d_r4_member

  subroutine zonal_sum_0d_r8_member(zonal_circle, work, value)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: work(:,:)
    real(8), intent(out) :: value(:)

#ifdef ENSURE_ORDER
    real(8) allvalue(member_num , global_mesh%num_full_lon)
#endif
    integer ierr

#ifdef ENSURE_ORDER
    if (zonal_circle%np == 1) then
      value(:) = sum(work,2)
    else
      call gather_zonal_array_member(zonal_circle, work, allvalue)
      if (zonal_circle%id == 0) value(:) = sum(allvalue,2)
      call MPI_BCAST(value, size(value), MPI_DOUBLE, 0, zonal_circle%comm, ierr)
    end if
#else
    call MPI_ALLREDUCE(sum(work,2), value, size(value), MPI_DOUBLE, MPI_SUM, zonal_circle%comm, ierr)
#endif

  end subroutine zonal_sum_0d_r8_member

  subroutine zonal_sum_1d_r4_member(zonal_circle, work, value)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(4), intent(in) :: work(:,:,:)
    real(4), intent(out) :: value(:,:)

#ifdef ENSURE_ORDER
    real(4) allvalue(size(value,1),global_mesh%num_full_lon,size(value,2))
#endif
    integer ierr

#ifdef ENSURE_ORDER
    if (zonal_circle%np == 1) then
      value = sum(work, dim=2)
    else
      call gather_zonal_array_member(zonal_circle, work, allvalue)
      if (zonal_circle%id == 0) value = sum(allvalue, dim=2)
      call MPI_BCAST(value, size(value), MPI_REAL, 0, zonal_circle%comm, ierr)
    end if
#else
    call MPI_ALLREDUCE(sum(work, dim=2), value, size(value), MPI_REAL, MPI_SUM, zonal_circle%comm, ierr)
#endif

  end subroutine zonal_sum_1d_r4_member

  subroutine zonal_sum_1d_r8_member(zonal_circle, work, value)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: work(:,:,:)
    real(8), intent(out) :: value(:,:)

#ifdef ENSURE_ORDER
    real(8) allvalue(size(value,1),global_mesh%num_full_lon,size(value,2))
#endif
    integer ierr

#ifdef ENSURE_ORDER
    if (zonal_circle%np == 1) then
      value = sum(work, dim=2)
    else
      call gather_zonal_array_member(zonal_circle, work, allvalue)
      if (zonal_circle%id == 0) value = sum(allvalue, dim=2)
      call MPI_BCAST(value, size(value), MPI_DOUBLE, 0, zonal_circle%comm, ierr)
    end if
#else
    call MPI_ALLREDUCE(sum(work, dim=2), value, size(value), MPI_DOUBLE, MPI_SUM, zonal_circle%comm, ierr)
#endif

  end subroutine zonal_sum_1d_r8_member

  subroutine zonal_sum_old_0d_r8(comm, value)

    integer, intent(in) :: comm
    real(8), intent(inout) :: value

    integer ierr
    real(8) res

    call MPI_ALLREDUCE(value, res, 1, MPI_DOUBLE, MPI_SUM, comm, ierr)
    value = res

  end subroutine zonal_sum_old_0d_r8

  subroutine zonal_sum_old_1d_r8(comm, value)

    integer, intent(in) :: comm
    real(8), intent(inout) :: value(:)

    integer ierr
    real(8) res(size(value))

    call MPI_ALLREDUCE(value, res, size(value), MPI_DOUBLE, MPI_SUM, comm, ierr)
    value = res

  end subroutine zonal_sum_old_1d_r8

  subroutine zonal_sum_old_2d_r8(comm, value)

    integer, intent(in) :: comm
    real(8), intent(inout) :: value(:,:)

    integer ierr
    real(8) res(size(value,1) , size(value,2))

    call MPI_ALLREDUCE(value, res, size(value), MPI_DOUBLE, MPI_SUM, comm, ierr)
    value = res

  end subroutine zonal_sum_old_2d_r8


  subroutine zonal_sum_ensure_order_1d_r8(comm, value , sum)

    integer, intent(in) :: comm
    real(8), intent(in) :: value(:)
    real(8), intent(out):: sum 

    integer status(MPI_STATUS_SIZE), ierr
    integer numproc , myid
    integer i 
    real(8), allocatable :: allvalue(:)  


    call MPI_COMM_SIZE(comm,numproc,ierr)
    call MPI_COMM_RANK(comm , myid , ierr)
    sum = 0.0_r8

    if (numproc == 1) then
      do i = 1 , size(value)
        sum = sum + value(i)
      end do
    else
   
      if (proc%cart_coords(1) == 0) then  !root proc
        allocate(allvalue(global_mesh%num_full_lon))
        allvalue(1:size(value)) = value(:)
        do i = 1 , numproc - 1
          call MPI_RECV(allvalue( i * size(value) + 1 : (i+1)  * size(value) ) , size(value) , MPI_DOUBLE , i  , 0 , comm , status,ierr)
        end do
        do i = 1 , global_mesh%num_full_lon
          sum = sum + allvalue(i)
        end do
        do i = 1 , numproc - 1
          call MPI_SEND(sum , 1 , MPI_DOUBLE , i , 0 , comm , status,ierr)
        end do
        deallocate(allvalue)
      else !other proc
        call MPI_SEND(value(:) , size(value) , MPI_DOUBLE , 0 , 0 , comm , ierr)
        call MPI_RECV(sum , 1 , MPI_DOUBLE , 0 , 0 , comm , status,ierr)
      end if  

    end if

  end subroutine zonal_sum_ensure_order_1d_r8

  subroutine zonal_sum_ensure_order_2d_r8(comm, value , sum , pos_name)

    integer, intent(in) :: comm
    real(8), intent(in) :: value(:,:)
    real(8), intent(out):: sum(:) 
    character(*), optional , intent(in) :: pos_name

    integer status(MPI_STATUS_SIZE), ierr
    integer numproc , myid
    integer i , nm , nx
    real(8), allocatable :: allvalue(:,:)  


    call MPI_COMM_SIZE(comm,numproc,ierr)
    call MPI_COMM_RANK(comm,myid,ierr)
    sum = 0.0_r8
    nm  = size(value,1)
    nx  = size(value,2)


    if (numproc == 1) then
      do i = 1 , nx
        sum = sum + value(:,i)
      end do
    else

      
      if (proc%cart_coords(1) == 0) then  !root proc
        allocate(allvalue(member_num , global_mesh%num_full_lon))
        allvalue(:,1:nx) = value
        do i = 1 , numproc - 1
          call MPI_RECV(allvalue( : , i * nx + 1 : (i+1)  * nx ) , nx * nm , MPI_DOUBLE , i  , 0 , comm , status,ierr)
        end do
        do i = 1 , global_mesh%num_full_lon
          sum = sum + allvalue(:,i)
        end do
        do i = 1 , numproc - 1
          call MPI_SEND(sum , nm , MPI_DOUBLE , i , 0 , comm , status,ierr)
        end do
        deallocate(allvalue)
        
      else !other proc
        call MPI_SEND(value(:,:) , nx * nm , MPI_DOUBLE , 0 , 0 , comm , ierr)
        call MPI_RECV(sum , nm , MPI_DOUBLE , 0 , 0 , comm , status,ierr)
      end if  



    end if

  end subroutine zonal_sum_ensure_order_2d_r8

  subroutine global_sum_0d_r4(comm, value)

    integer, intent(in) :: comm
    real(4), intent(inout) :: value

    integer ierr
    real(4) res

    call MPI_ALLREDUCE(value, res, 1, MPI_REAL, MPI_SUM, comm, ierr)
    value = res

  end subroutine global_sum_0d_r4

  subroutine global_sum_0d_r8(comm, value)

    integer, intent(in) :: comm
    real(8), intent(inout) :: value

    integer ierr
    real(8) res

    call MPI_ALLREDUCE(value, res, 1, MPI_DOUBLE, MPI_SUM, comm, ierr)
    value = res

  end subroutine global_sum_0d_r8

  subroutine global_sum_1d_r4(comm, value)

    integer, intent(in) :: comm
    real(4), intent(inout) :: value(:)

    integer ierr
    real(4) res(size(value))

    call MPI_ALLREDUCE(value, res, size(value), MPI_REAL, MPI_SUM, comm, ierr)
    value = res

  end subroutine global_sum_1d_r4

  subroutine global_sum_1d_r8(comm, value)

    integer, intent(in) :: comm
    real(8), intent(inout) :: value(:)

    integer ierr
    real(8) res(size(value))

    call MPI_ALLREDUCE(value, res, size(value), MPI_DOUBLE, MPI_SUM, comm, ierr)
    value = res

  end subroutine global_sum_1d_r8
  

  subroutine barrier()

    integer ierr

    call MPI_BARRIER(proc%comm, ierr)

  end subroutine barrier

  subroutine gather_zonal_array_1d_r8(zonal_circle, local_array, array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: local_array(:)
    real(8), intent(out) :: array(:)

    integer ierr, i, status

    if (zonal_circle%id == 0) then
      array(1:size(local_array)) = local_array
      do i = 2, zonal_circle%np
        call MPI_RECV(array, 1, zonal_circle%recv_type_r8(i,0), i - 1, 30, zonal_circle%comm, status, ierr)
      end do
    else
      call MPI_SEND(local_array, size(local_array), MPI_DOUBLE, 0, 30, zonal_circle%comm, ierr)
    end if

  end subroutine gather_zonal_array_1d_r8

  subroutine gather_zonal_array_2d_r8(zonal_circle, local_array, array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: local_array(:,:)
    real(8), intent(inout) :: array(:,:)

    integer ierr, i, k, status

    if (zonal_circle%id == 0) then
      k = merge(1, 2, size(local_array, 2) == global_mesh%num_full_lev)
      array(1:size(local_array, 1),:) = local_array
      do i = 2, zonal_circle%np
        call MPI_RECV(array, 1, zonal_circle%recv_type_r8(i,k), i - 1, 31, zonal_circle%comm, status, ierr)
      end do
    else
      call MPI_SEND(local_array, size(local_array), MPI_DOUBLE, 0, 31, zonal_circle%comm, ierr)
    end if

  end subroutine gather_zonal_array_2d_r8

  subroutine scatter_zonal_array_1d_r8(zonal_circle, array, local_array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: array(:)
    real(8), intent(out) :: local_array(:)

    integer ierr, i, status

    if (zonal_circle%id == 0) then
      local_array = array(1:size(local_array))
      do i = 2, zonal_circle%np
        call MPI_SEND(array, 1, zonal_circle%recv_type_r8(i,0), i - 1, 32, zonal_circle%comm, ierr)
      end do
    else
      call MPI_RECV(local_array, size(local_array), MPI_DOUBLE, 0, 32, zonal_circle%comm, status, ierr)
    end if

  end subroutine scatter_zonal_array_1d_r8

  subroutine scatter_zonal_array_2d_r8(zonal_circle, array, local_array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: array(:,:)
    real(8), intent(out) :: local_array(:,:)

    integer ierr, i, k, status

    if (zonal_circle%id == 0) then
      k = merge(1, 2, size(local_array, 2) == global_mesh%num_full_lev)
      local_array = array(1:size(local_array, 1),:)
      do i = 2, zonal_circle%np
        call MPI_SEND(array, 1, zonal_circle%recv_type_r8(i,k), i - 1, 33, zonal_circle%comm, ierr)
      end do
    else
      call MPI_RECV(local_array, size(local_array), MPI_DOUBLE, 0, 33, zonal_circle%comm, status, ierr)
    end if

  end subroutine scatter_zonal_array_2d_r8

  subroutine gather_zonal_array_1d_r4_member(zonal_circle, local_array, array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(4), intent(in) :: local_array(:,:)
    real(4), intent(out) :: array(:,:)

    integer ierr, i, status(MPI_STATUS_SIZE)

    if (zonal_circle%id == 0) then
      array(:,1:size(local_array,2)) = local_array
      do i = 2, zonal_circle%np
        call MPI_RECV(array, 1, zonal_circle%recv_type_r4(i,0), i - 1, 30, zonal_circle%comm, status, ierr)
      end do
    else
      call MPI_SEND(local_array, size(local_array) , MPI_REAL, 0, 30, zonal_circle%comm, ierr)
    end if

  end subroutine gather_zonal_array_1d_r4_member

  subroutine gather_zonal_array_1d_r8_member(zonal_circle, local_array, array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: local_array(:,:)
    real(8), intent(out) :: array(:,:)

    integer ierr, i, status(MPI_STATUS_SIZE)

    if (zonal_circle%id == 0) then
      array(:,1:size(local_array,2)) = local_array
      do i = 2, zonal_circle%np
        call MPI_RECV(array, 1, zonal_circle%recv_type_r8(i,0), i - 1, 30, zonal_circle%comm, status, ierr)
      end do
    else
      call MPI_SEND(local_array, size(local_array) , MPI_DOUBLE, 0, 30, zonal_circle%comm, ierr)
    end if

  end subroutine gather_zonal_array_1d_r8_member

  subroutine gather_zonal_array_2d_r4_member(zonal_circle, local_array, array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(4), intent(in) :: local_array(:,:,:)
    real(4), intent(inout) :: array(:,:,:)

    integer ierr, i, k, status(MPI_STATUS_SIZE)

    if (zonal_circle%id == 0) then
      k = merge(1, 2, size(local_array, 2) == global_mesh%num_full_lev)
      array(:,1:size(local_array, 2),:) = local_array
      do i = 2, zonal_circle%np
        call MPI_RECV(array, 1, zonal_circle%recv_type_r4(i,k), i - 1, 31, zonal_circle%comm, status, ierr)
      end do
    else
      call MPI_SEND(local_array, size(local_array), MPI_REAL, 0, 31, zonal_circle%comm, ierr)
    end if

  end subroutine gather_zonal_array_2d_r4_member

  subroutine gather_zonal_array_2d_r8_member(zonal_circle, local_array, array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: local_array(:,:,:)
    real(8), intent(inout) :: array(:,:,:)

    integer ierr, i, k, status(MPI_STATUS_SIZE)

    if (zonal_circle%id == 0) then
      k = merge(1, 2, size(local_array, 2) == global_mesh%num_full_lev)
      array(:,1:size(local_array, 2),:) = local_array
      do i = 2, zonal_circle%np
        call MPI_RECV(array, 1, zonal_circle%recv_type_r8(i,k), i - 1, 31, zonal_circle%comm, status, ierr)
      end do
    else
      call MPI_SEND(local_array, size(local_array), MPI_DOUBLE, 0, 31, zonal_circle%comm, ierr)
    end if

  end subroutine gather_zonal_array_2d_r8_member

  subroutine scatter_zonal_array_1d_r4_member(zonal_circle, array, local_array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(4), intent(in) :: array(:,:)
    real(4), intent(out) :: local_array(:,:)

    integer ierr, i, status(MPI_STATUS_SIZE)

    if (zonal_circle%id == 0) then
      local_array = array(:,1:size(local_array))
      do i = 2, zonal_circle%np
        call MPI_SEND(array, 1, zonal_circle%recv_type_r4(i,0), i - 1, 32, zonal_circle%comm, ierr)
      end do
    else
      call MPI_RECV(local_array, size(local_array), MPI_REAL, 0, 32, zonal_circle%comm, status, ierr)
    end if

  end subroutine scatter_zonal_array_1d_r4_member

  subroutine scatter_zonal_array_1d_r8_member(zonal_circle, array, local_array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: array(:,:)
    real(8), intent(out) :: local_array(:,:)

    integer ierr, i, status(MPI_STATUS_SIZE)

    if (zonal_circle%id == 0) then
      local_array = array(:,1:size(local_array))
      do i = 2, zonal_circle%np
        call MPI_SEND(array, 1, zonal_circle%recv_type_r8(i,0), i - 1, 32, zonal_circle%comm, ierr)
      end do
    else
      call MPI_RECV(local_array, size(local_array), MPI_DOUBLE, 0, 32, zonal_circle%comm, status, ierr)
    end if

  end subroutine scatter_zonal_array_1d_r8_member

  subroutine scatter_zonal_array_2d_r4_member(zonal_circle, array, local_array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(4), intent(in) :: array(:,:,:)
    real(4), intent(out) :: local_array(:,:,:)

    integer ierr, i, k, status(MPI_STATUS_SIZE)

    if (zonal_circle%id == 0) then
      k = merge(1, 2, size(local_array, 2) == global_mesh%num_full_lev)
      local_array = array(:,1:size(local_array, 1),:)
      do i = 2, zonal_circle%np
        call MPI_SEND(array, 1, zonal_circle%recv_type_r4(i,k), i - 1, 33, zonal_circle%comm, ierr)
      end do
    else
      call MPI_RECV(local_array, size(local_array), MPI_REAL, 0, 33, zonal_circle%comm, status, ierr)
    end if

  end subroutine scatter_zonal_array_2d_r4_member

  subroutine scatter_zonal_array_2d_r8_member(zonal_circle, array, local_array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: array(:,:,:)
    real(8), intent(out) :: local_array(:,:,:)

    integer ierr, i, k, status(MPI_STATUS_SIZE)

    if (zonal_circle%id == 0) then
      k = merge(1, 2, size(local_array, 2) == global_mesh%num_full_lev)
      local_array = array(:,1:size(local_array, 1),:)
      do i = 2, zonal_circle%np
        call MPI_SEND(array, 1, zonal_circle%recv_type_r8(i,k), i - 1, 33, zonal_circle%comm, ierr)
      end do
    else
      call MPI_RECV(local_array, size(local_array), MPI_DOUBLE, 0, 33, zonal_circle%comm, status, ierr)
    end if

  end subroutine scatter_zonal_array_2d_r8_member

end module parallel_mod
