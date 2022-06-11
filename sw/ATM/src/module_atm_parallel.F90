module atm_parallel
  use mpi
  use coupler_config

  implicit none
  integer :: rank, lrank, lsize
  integer :: comm_cart
  integer :: dims(2), coords(2)
  logical :: periods(2)

  integer :: x_sec, y_sec, x_res, y_res
  integer :: xst, xed, yst, yed
  integer :: x_size, y_size

  contains

  subroutine init_atm_parallel(MPI_ATM_GROUP)

    integer, intent(in) :: MPI_ATM_GROUP
    integer :: ierrs
    dims(1:2) = 0

    call MPI_Comm_rank(MPI_ATM_GROUP, rank, ierrs)
    call MPI_Comm_size(MPI_ATM_GROUP, lsize, ierrs)
    call MPI_Dims_create(lsize, 2, dims, ierrs)
    call MPI_Cart_create(MPI_ATM_GROUP, 2, dims, periods, 0, comm_cart, ierrs)
    call MPI_Comm_rank(comm_cart, lrank, ierrs)
    call MPI_Cart_coords(comm_cart, lrank, 2, coords)

    x_sec = atm_num_lon / dims(1)
    x_res = mod(atm_num_lon, dims(1))
    y_sec = atm_num_lat / dims(2)
    y_res = mod(atm_num_lat, dims(2))
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
  end subroutine init_atm_parallel
end module atm_parallel
