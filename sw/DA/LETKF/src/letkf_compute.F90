! Calculate the data assimilation in current process

! ens        : ensemble number
! nx, ny     : block size in a process, nx, ny horizonal
! halo       : size of halo region
! n          : number of prognostic variables
! cube_n     : letkf cube size
! overlap    : overlap of letkf cubes
! total_obs  : total observation number of a single process
! obs_list   : a structure for observation data

subroutine letkf_compute_nooverlap(xb, axb, xa, axa)

    ! background
    xb_type,  intent(in) :: xb(ens, n, 1-halo:nx+halo, 1-halo:ny+halo)
    axb_type, intent(in) :: axb(n, 1-halo:nx+halo, 1-halo:ny+halo)
    ! analysis
    xb_type,  intent(out) :: xa(ens, n, 1-halo:nx+halo, 1-halo:ny+halo)
    axb_type, intent(out) :: axa(n, 1-halo:nx+halo, 1-halo:ny+halo)

    ! yb       : objected observation variation
    ! ayb      : average objected observation
    ! dis      : y_o - ayb
    ! rloc     : localization function
    ! rdiaginv : observation error
    yb_type  :: yb(ens, MAX_OBS)
    ayb_type :: ayb(MAX_OBS)
    ayb_type :: dis(MAX_OBS)
    yb_type  :: rloc(MAX_OBS)
    real(4)  :: rdiaginv(MAX_OBS)

    ! result of letkf_core
    yb_type :: trans(ens, ens)
    yb_type  :: overline_w(ens)

    ! number of observations of current cube
    integer :: nobsl

    ! number of observations of current horizonal position
    integer :: nobsl_h

    ! for error calc
    logical :: calc_again

    integer :: i, j
    ! cube start index and end index
    integer :: si, sj, ei, ej

#if (DEBUG == 1)
    integer :: debugi, debugj
    real(8) :: debug_sum
    integer :: pflag
#endif

    do i = 1, ny, cube_n

        si = i
        ei = i + cube_n - 1
        if (si < 1)  si = 1
        if (ei > ny) ei = ny

        do j = 1, nx, cube_n

            sj = j
            ej = j + cube_n - 1
            if (sj < 1)  sj = 1
            if (ej > nx) ej = nx    

            !sk = 1
            !ek = nz
#if (DEBUG == 1)
            call GPTLstart("Generate Yb matrix")
#endif
            call gen_yb_dis(sj, si, ej, ei, nobsl, yb, dis, rdiaginv, rloc)

#if (DEBUG == 1)
            call GPTLstop("Generate Yb matrix")
#endif

#if (DEBUG == 1)
            call GPTLstart("Lektf computation core")
#endif
            pflag = 0
            !if ((pid == 0 .or. pid == 2) .and. i == 1 .and. j == 1) then
            if (pid == 69 .and. i == 6 .and. j == 2) then
                !pflag = 1
                !print *, "id", pid, "ij", i, j, "nobsl:", nobsl

                open(1, file="data.out", form="FORMATTED")
                write(1, *) nobsl
                write(1, *) ens
                write(1, *) yb
                write(1, *) dis
                write(1, *) rdiaginv
                write(1, *) rloc
                !write(*, *) "yb"
                !do debugi = 1, nobsl
                !    do debugj = 1, ens
                !        write(*, 668) yb(debugj, debugi)
                !    end do
                !    write(*, *) " "
                !end do
                !write (*, *) "dis"
                !do debugi = 1, nobsl
                !    write(*, 668) dis(debugi)
                !end do
                !write(*, *) "rdiaginv"
                !do debugi = 1, nobsl
                !    write(*, 668) rdiaginv(debugi)
                !end do
                !write(*, *) "rloc"
                !do debugi = 1, nobsl
                !    write(*, 668) rloc(debugi)
                !end do
                !write(*, *) " "
668 FORMAT(F8.4,$)
                !print *, "yb", yb
                !print *, "rdiaginv", rdiaginv
                !print *, "rloc", rloc
                !pflag = 1
            end if

            !if (pid == 46 .and. i == 1 .and. j == 1) then
            !    print *, "yb", yb
            !    print *, "rdiaginv", rdiaginv
            !    print *, "dis", dis
            !    print *, "rloc", rloc
            !end if

            if (precision .ge. 52) then
                call letkf_core(ens, nobsl, yb(1: ens, 1: nobsl), rdiaginv(1: nobsl), &
                dis(1: nobsl), rloc(1: nobsl), trans, overline_w, inflation_factor, pflag)
            else
                call letkf_core_half(ens, nobsl, yb(1: ens, 1: nobsl), rdiaginv(1: nobsl), &
                dis(1: nobsl), rloc(1: nobsl), trans, overline_w, inflation_factor, pflag)
            end if

#if (DEBUG == 1)
            calc_again = .false.
            do debugi = 1, ens
                do debugj = 1, ens
                    if (isnan(trans(debugj, debugi))) then
                        print *, "pid ", pid, i, j, "nan!!!"
                        !calc_again = .true.
                        !goto 20
                        call log_error("trans nan")
                    end if
                end do
                if (isnan(overline_w(debugi))) then
                    print *, "pid", pid, i, j, "nan!!!"
                    !calc_again = .true.
                    !goto 20
                    call log_error("overline_w nan ")
                end if
                !print *, "trans", trans
            end do
20          if (calc_again) then
                call letkf_core(ens, nobsl, yb(1: ens, 1: nobsl), rdiaginv(1: nobsl), &
                dis(1: nobsl), rloc(1: nobsl), trans, overline_w, inflation_factor, pflag)
            end if
#endif

#if (DEBUG == 1)
            call GPTLstop("Lektf computation core")
#endif

#if (DEBUG == 1)
            call GPTLstart("Convert xb to xa")
#endif
            call xb_2_xa(ej-sj+1, ei-si+1, &
            xb(:, :, sj:ej, si:ei), &
            axb(:, sj:ej, si:ei), &
            xa(:, :, sj:ej, si:ei), &
            axa(:, sj:ej, si:ei), &
            trans, overline_w)

#if (DEBUG == 1)
            call GPTLstop("Convert xb to xa")
#endif

            !end do

        end do
    end do
end subroutine

subroutine letkf_compute(xb, axb, xa, axa)

    ! background
    xb_type,  intent(in) :: xb(ens, n, 1-halo:nx+halo, 1-halo:ny+halo)
    axb_type, intent(in) :: axb(n, 1-halo:nx+halo, 1-halo:ny+halo)
    ! analysis
    xb_type,  intent(out) :: xa(ens, n, 1-halo:nx+halo, 1-halo:ny+halo)
    axb_type, intent(out) :: axa(n, 1-halo:nx+halo, 1-halo:ny+halo)

    ! halo region for overlapped cubes
    integer :: halo_s

    ! yb       : objected observation variation
    ! ayb      : average objected observation
    ! dis      : y_o - ayb
    ! rloc     : localization function
    ! rdiaginv : observation error
    yb_type  :: yb(ens, MAX_OBS)
    ayb_type :: ayb(MAX_OBS)
    ayb_type :: dis(MAX_OBS)
    yb_type  :: rloc(MAX_OBS)
    real(4)  :: rdiaginv(MAX_OBS)

    ! result of letkf_core
    yb_type :: trans(ens, ens)
    yb_type :: overline_w(ens)

    ! number of observations of current cube
    integer :: nobsl

    ! number of observations of current horizonal position
    integer :: nobsl_h

    integer :: i, j, k

    ! loop for a single cube
    integer :: ii, jj
    ! cube start index and end index
    integer :: si, sj, ei, ej

    ! handle overlap of cubes
    yb_type :: cube_count(1-overlap/2:nx+overlap/2, 1-overlap/2:ny+overlap/2)

    halo_s = overlap / 2

    xa(:, :, :, :) = 0
    axa(:, :, :) = 0
    cube_count(:, :) = 0

    do i = 1, ny, cube_n-overlap

        si = i
        ei = i + cube_n - 1
        if (si < 1)  si = 1
        if (ei > ny) ei = ny

        do j = 1, nx, cube_n-overlap

            sj = j
            ej = j + cube_n - 1
            if (sj < 1)  sj = 1
            if (ej > nx) ej = nx

            ! sk = 1
            ! ek = nz
            
            call gen_yb_dis(sj, si, ej, ei, nobsl, yb, dis, rdiaginv, rloc)

            call letkf_core(ens, nobsl, yb(1: ens, 1: nobsl), rdiaginv(1: nobsl), &
            dis(1: nobsl), rloc(1: nobsl), trans, overline_w, inflation_factor)

            ! call xb_2_xa(ej-sj+1, ei-si+1, &
            ! xb(:, :, sj:ej, si:ei), &
            ! axb(:, sj:ej, si:ei), &
            ! xa(:, :, sj:ej, si:ei), &
            ! axa(:, sj:ej, si:ei), &
            ! trans, overline_w)

            call xb_add2_xa(ej-sj+1, ei-si+1, ej, ei, &
            xb(:, :, sj:ej, si:ei), &
            axb(:, sj:ej, si:ei), &
            xa(:, :, sj:ej, si:ei), &
            axa(:, sj:ej, si:ei), &
            trans, overline_w)

            do ii = si, ei
                do jj = sj, ej
                    cube_count(jj, ii) = cube_count(jj, ii) + 1
                end do
            end do
        end do
    end do

    do i = 1, ny
        do j = 1, nx
            xa(:, :, j, i) = xa(:, :, j, i) / cube_count(j, i)
            axa(:, j, i) = axa(:, j, i) / cube_count(j, i)
            !if (pid == 0 .and. mod(j, 10) == 0 .and. mod(i, 10) == 0) print *, pid, axa(1, j, i), j, i
        end do
    end do
end subroutine

! use trans and overline_w to get analysis from background
subroutine xb_2_xa(dy, dx, xb, axb, xa, axa, trans, overline_w)

    !cube size
    integer, intent(in) :: dy, dx

    ! background
    xb_type,  intent(in) :: xb(ens, n, 1:dy, 1:dx)
    axb_type, intent(in) :: axb(n, 1:dy, 1:dx)
    ! analysis
    xb_type,  intent(out) :: xa(ens, n, 1:dy, 1:dx)
    axb_type, intent(out) :: axa(n, 1:dy, 1:dx)

    yb_type, intent(in) :: trans(ens, ens)
    yb_type,  intent(in) :: overline_w(ens)

    integer :: i, j, k, m, var, ens_num, debugi

    pointer(xb_ptr, xb_matrix)
    pointer(xa_ptr, xa_matrix)

    xb_type :: xb_matrix(ens, n*dy*dx)
    xb_type :: xa_matrix(ens, n*dy*dx)

    if (kind(xb) == kind(trans)) then

        xb_ptr = loc(xb(1, 1, 1, 1))
        xa_ptr = loc(xa(1, 1, 1, 1))

        m = dx * dy * n
        
        if (yb_size == 8) then
            call dgemm('t', 'n', ens, m, ens, 1.d0, trans, ens, &
            xb_matrix, ens, 0.d0, xa_matrix, ens)
        else
            call sgemm('t', 'n', ens, m, ens, 1.e0, trans, ens, &
            xb_matrix, ens, 0.e0, xa_matrix, ens)
        end if
        do i = 1, dx
            do j = 1, dy
                do var = 1, n
                    xa(:, var, j, i) = xa(:, var, j, i) + axb(var, j, i)
                end do
            end do
        end do
    else
        do i = 1, dx
            do j = 1, dy
                do var = 1, n
                    do ens_num = 1, ens
                        xa(ens_num, var, j, i) = axb(var, j, i) + dot_product(xb(1:ens, var, j, i), trans(1:ens, ens_num))
                    end do
                end do
            end do
        end do
    end if
    do i = 1, dx
        do j = 1, dy
            do var = 1, n
                axa(var, j, i) = 0
                do ens_num = 1, ens
                    axa(var, j, i) = axa(var, j, i) + xa(ens_num, var, j, i)
                end do
                axa(var, j, i) = axa(var, j, i) / ens
                do ens_num = 1, ens
                    xa(ens_num, var, j, i) = xa(ens_num, var, j, i) - axa(var, j, i)
                end do
            end do
        end do
    end do

end subroutine

! use trans and overline_w to get analysis from background
subroutine xb_add2_xa(dy, dx, sj, si, xb, axb, xa, axa, trans, overline_w)

    !cube size
    integer, intent(in) :: dy, dx

    !cube position
    integer, intent(in) :: si, sj

    ! background
    xb_type,  intent(in) :: xb(ens, n, 1:dy, 1:dx)
    axb_type, intent(in) :: axb(n, 1:dy, 1:dx)
    ! analysis
    xb_type,  intent(inout) :: xa(ens, n, 1:dy, 1:dx)
    axb_type, intent(inout) :: axa(n, 1:dy, 1:dx)

    xb_type  :: dist(ens)
    axb_type :: mean
    ! integer :: halo_s
    ! yb_type :: halo_factor

    yb_type, intent(in) :: trans(ens, ens)
    yb_type,  intent(in) :: overline_w(ens)

    integer :: i, j, var, ens_num

    !halo_s = overlap / 2

    do i = 1, dx
        do j = 1, dy
            ! halo_factor = 1
            ! if ((i <= overlap .and. si .ne. 1 - halo_s) .or. (dx - i + 1 <= overlap .and. si + dx - 1 .ne. nx + halo_s)) &
            ! halo_factor = halo_factor * 2
            ! if ((j <= overlap .and. sj .ne. 1 - halo_s) .or. (dy - j + 1 <= overlap .and. sj + dy - 1 .ne. ny + halo_s)) &
            ! halo_factor = halo_factor * 2
            do var = 1, n
                mean = 0
                do ens_num = 1, ens
                    dist(ens_num) = axb(var, j, i) + dot_product(xb(1:ens, var, j, i), trans(1:ens, ens_num))
                    ! dist(ens_num) = dist(ens_num) / halo_factor
                    mean = mean + dist(ens_num)
                end do
                mean = mean / ens
                axa(var, j, i) = axa(var, j, i) + mean
                do ens_num = 1, ens
                    xa(ens_num, var, j, i) = xa(ens_num, var, j, i) + dist(ens_num) - mean
                end do
            end do
        end do
    end do

end subroutine

! With vertical localization

subroutine letkf_compute_nooverlap_vertical(xb, axb, xa, axa)

    ! background
    xb_type,  intent(in) :: xb(ens, n, 1-halo:nx+halo, 1-halo:ny+halo)
    axb_type, intent(in) :: axb(n, 1-halo:nx+halo, 1-halo:ny+halo)
    ! analysis
    xb_type,  intent(out) :: xa(ens, n, 1-halo:nx+halo, 1-halo:ny+halo)
    axb_type, intent(out) :: axa(n, 1-halo:nx+halo, 1-halo:ny+halo)

    ! yb       : objected observation variation
    ! ayb      : average objected observation
    ! dis      : y_o - ayb
    ! rloc     : localization function
    ! rdiaginv : observation error
    yb_type  :: yb(ens, MAX_OBS)
    ayb_type :: ayb(MAX_OBS)
    ayb_type :: dis(MAX_OBS)
    yb_type  :: rloc_no_vertical(MAX_OBS)
    yb_type  :: rloc(MAX_OBS)
    real(4)  :: rdiaginv(MAX_OBS)

    ayb_type :: alnp(nz)
    ayb_type :: obs_lnp(MAX_OBS)

    ! result of letkf_core
    yb_type :: trans(ens, ens)
    yb_type  :: overline_w(ens)

    !yb_type, allocatable, dimension(:, :) :: yb_tmp

    ! number of observations of current cube
    integer :: nobsl

    ! number of observations of current horizonal position
    integer :: nobsl_h

    integer :: i, j, k
    ! cube start index and end index
    integer :: si, sj, ei, ej, sk, ek

    integer :: var_iter, k_iter, ens_iter

    integer :: center_i, center_j, center_k

    integer :: dx, dy
    integer :: ni, nj, nk, ens_num, var_pos

    integer :: pflag

    do i = 1, ny, cube_n

        si = i
        ei = i + cube_n - 1
        if (si < 1)  si = 1
        if (ei > ny) ei = ny

        center_i = (si + ei) / 2

        do j = 1, nx, cube_n

            sj = j
            ej = j + cube_n - 1
            if (sj < 1)  sj = 1
            if (ej > nx) ej = nx

            center_j = (sj + ej) / 2

            call gen_yb_dis_vertical(sj, si, ej, ei, nobsl, yb, dis, rdiaginv, rloc_no_vertical, obs_lnp)

            call ps2alnp(axb(ps_indx, center_j, center_i), alnp(1: nz))

            pflag = 0
            if ((pid == 0 .or. pid == 2) .and. i == 1 .and. j == 1) then
                !pflag = 1
                print *, "vertical id", pid, "ij", i, j, "nobsl", nobsl
            end if

            !allocate(yb_tmp(ens, nobsl))

            do k = 1, nz, cube_nz

                sk = k
                ek = k + cube_nz - 1
                if (sk < 1)  sk = 1
                if (ek > nz) ek = nz

                center_k = (sk + ek) / 2

                call vertical_localization(rloc_no_vertical, rloc, nobsl, alnp(center_k), obs_lnp(1: nobsl))

                if (pid == 0 .and. i == 1 .and. j == 1) then
                    !print *, "no vertical", rloc_no_vertical(1: nobsl)
                    !print *, "after vertical", rloc(1: nobsl)
                end if

                !yb_tmp(1: ens, 1: nobsl) = yb(1: ens, 1: nobsl)

                if (precision .ge. 52) then
                    call letkf_core(ens, nobsl, yb(1: ens, 1: nobsl), rdiaginv(1: nobsl), &
                    dis(1: nobsl), rloc(1: nobsl), trans, overline_w, inflation_factor, pflag)
                else
                    !call letkf_core_mpal(ens, nobsl, yb(1: ens, 1: nobsl), rdiaginv(1: nobsl), &
                    !dis(1: nobsl), rloc(1: nobsl), trans, overline_w, inflation_factor, pflag)
                end if

                do var_iter = 1, n_var
                    if (da_var_lev(var_iter) == 1) then
                        if (ek .eq. nz) then
                            
                            var_pos = da_var_start(var_iter)

                            do ni = si, ei
                                do nj = sj, ej
                                    do ens_num = 1, ens
                                        xa(ens_num, var_pos, nj, ni) = axb(var_pos, nj, ni) + dot_product(xb(1:ens, var_pos, nj, ni), trans(1:ens, ens_num))
                                    end do
                                end do
                            end do

                        end if
                    else

                        do ni = si, ei
                            do nj = sj, ej
                                do nk = da_var_start(var_iter) + sk - 1, &
                                    da_var_start(var_iter) + ek - 1
                                    do ens_num = 1, ens
                                        xa(ens_num, nk, nj, ni) = axb(nk, nj, ni) + dot_product(xb(1:ens, nk, nj, ni), trans(1:ens, ens_num))
                                    end do
                                end do
                            end do
                        end do

                    end if
                end do

            end do

            !deallocate(yb_tmp)

        end do
    end do

    do i = 1, ny
        do j = 1, nx
            do var_iter = 1, n
                axa(var_iter, j, i) = 0
                do ens_num = 1, ens
                    axa(var_iter, j, i) = axa(var_iter, j, i) + xa(ens_num, var_iter, j, i)
                end do
                axa(var_iter, j, i) = axa(var_iter, j, i) / ens
                do ens_num = 1, ens
                    xa(ens_num, var_iter, j, i) = xa(ens_num, var_iter, j, i) - axa(var_iter, j, i)
                end do
            end do
        end do
    end do

end subroutine
