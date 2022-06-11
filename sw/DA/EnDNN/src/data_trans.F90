! nx : latitude, ny : longitude
subroutine struct_atm2endnn(ens, n, nx, ny, &
    ens_max, n_max, ens_start, n_start, xb, xb2)
    integer, intent(in) :: ens, n, nx, ny
    integer, intent(in) :: ens_max, n_max, ens_start, n_start
    axb_type, intent(in) :: xb(ens, ny, nx, n)
    axb_type, intent(inout) :: xb2(ens_max, n_max, nx, ny)
    integer :: ens_end
    integer :: nn, i, j, k, ens_num
    ens_end = ens_start + ens - 1
    do nn = n_start, n_start + n - 1
        do i = 1, nx
            do j = 1, ny
                xb2(ens_start:ens_end, nn, i, j) = xb(1:ens, j, i, nn - n_start + 1)
            end do
        end do
    end do
end subroutine

! nx : latitude, ny : longitude
subroutine struct_endnn2atm(ens, n, nx, ny, &
    ens_max, n_max, ens_start, n_start, xb, xb2)
    integer, intent(in) :: ens, n, nx, ny
    integer, intent(in) :: ens_max, n_max, ens_start, n_start
    axb_type, intent(in) :: xb(ens_max, n_max, nx, ny)
    axb_type, intent(out) :: xb2(ens, ny, nx, n)

    integer :: ens_end
    integer :: nn, i, j, ens_num
    ens_end = ens_start + ens - 1
    do j = 1, ny
        do i = 1, nx
            do nn = 1, n
                xb2(1:ens, j, i, nn) = xb(ens_start:ens_end, nn + n_start - 1, i, j)
            end do
        end do
    end do
end subroutine

! nx : latitude, ny : longitude
subroutine struct_atm2endnn_a(n, nx, ny, n_max, n_start, axb, axb2)
    integer, intent(in) :: n, nx, ny
    integer, intent(in) :: n_max, n_start
    axb_type, intent(in) :: axb(ny, nx, n)
    axb_type, intent(inout) :: axb2(n_max, nx, ny)
    integer :: nn, i, j
    do nn = n_start, n_start + n - 1
        do i = 1, nx
            do j = 1, ny
                axb2(nn, i, j) = axb(j, i, nn - n_start + 1)
            end do
        end do
    end do
end subroutine

! nx : latitude, ny : longitude
subroutine struct_endnn2atm_a(n, nx, ny, n_max, n_start, axb, axb2)
    integer, intent(in) :: n, nx, ny
    integer, intent(in) :: n_max, n_start
    axb_type, intent(in) :: axb(n_max, nx, ny)
    axb_type, intent(out) :: axb2(ny, nx, n)
    integer :: nn, i, j
    do j = 1, ny
        do i = 1, nx
            do nn = 1, n
                axb2(j, i, nn) = axb(nn + n_start - 1, i, j)
            end do
        end do
    end do
end subroutine


