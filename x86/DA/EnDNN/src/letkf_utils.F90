
subroutine pos2id(pi, pj, id)
    integer, intent(in) :: pi, pj
    integer, intent(out):: id
    integer :: myrank, ierrs

    integer :: my_dims(2), my_coords(2)
    logical :: my_periods(2)

    call MPI_Cart_get(comm_cart, 2, my_dims, my_periods, my_coords, ierrs)
    !print *, my_dims(1), my_dims(2)
    if (pi < 0 .or. pi >= my_dims(2)) then
      id = -1
    else
      id = mod(pj + my_dims(1), my_dims(1)) * my_dims(2) + pi
    end if
  
end subroutine