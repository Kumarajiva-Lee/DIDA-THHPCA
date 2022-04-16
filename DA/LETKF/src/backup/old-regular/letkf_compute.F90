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
            ! if (pid == 0 .and. mod(i, 10) == 0 .and. mod(j, 10) == 0) then
            !     write (*,*) pid, 'nobsl', nobsl
            !     do debugi = 1, nobsl
            !         do debugj = 1, ens
            !             write(*,668) yb(debugj, debugi)
            !         end do
            !         write(*, *) "dis,rinv,rloc ", dis(debugi), rdiaginv(debugi), rloc(debugi)
            !     end do
            ! end if
#endif

#if (DEBUG == 1)
            call GPTLstart("Lektf computation core")
#endif
            pflag = 0
            if (pid == 1 .and. i == 1 .and. j == 1) pflag = 1

            call letkf_core2(ens, nobsl, yb(1: ens, 1: nobsl), rdiaginv(1: nobsl), &
            dis(1: nobsl), rloc(1: nobsl), trans, overline_w, inflation_factor, pflag, 0.4d0)

#if (DEBUG == 1)
            call GPTLstop("Lektf computation core")
            ! if (pid == 0) then
            ! !     write(*,*),'overline_w',overline_w
            !     write (*,*) pid, i, j, 'trans: '
            !     do debugi = 1, ens
            ! !         debug_sum = 0
            !         do debugj = 1, ens
            !             write(*,668) trans(debugj, debugi)-overline_w(debugj)
            ! !             debug_sum = debug_sum + trans(debugj, debugi)
            !         end do
            ! !         write(*, 668) debug_sum
            !         write(*, *) " "
            !     end do
            ! end if
            ! 668 FORMAT(F10.2\)
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
            ! if (pflag) then
            !     write (*,*) pid, 'xa'
            !     do debugi = 1, ens
            !         write(*,101) xa(debugi, 2, j, i), xb(debugi, 2, j, i)
            !         write(*, *) " "
            !     end do
            !     write (*,*) pid, 'axa'
            !     write(*,101) axa(2, j, i)  ,axb(2, j, i)
            !     write(*, *) " "   
            ! end if
101 FORMAT(F10.2\)

            !if (pflag) call log_error("end")

            ! finding nan
            do debugi = 1, n
                if (isnan(axa(debugi, j, i))) then
                    print *, pid, i, j, "error"
                    call log_error("nan")
                end if
                do debugj = 1, ens
                    if (isnan(xa(debugj, debugi, j, i))) then
                        print *, pid, i, j, "error"
                        call log_error("nan")
                    end if
                end do
            end do
!             if (grid_lat(j) >= 29d0 .and. grid_lat(j) <= 31d0 .and. grid_lon(i) >= 59d0 .and. grid_lon(i) <= 61d0) then
!                 write (*,*) pid, 'xa'
!                 debug_sum = 0
!                 do debugi = 1, ens
!                     write(*,101) xa(debugi, 1, j, i)
!                     debug_sum = debug_sum + xa(debugi, 1, j, i)**2
!                 end do
!                 write(*,101) axa(1, j, i)
!                 write(*,101) debug_sum / ens
!                 write(*, *) " "
!                 write (*,*) pid, 'xb'
!                 debug_sum = 0
!                 do debugi = 1, ens
!                     write(*,101) xb(debugi, 1, j, i)
!                     debug_sum = debug_sum + xb(debugi, 1, j, i)**2
!                 end do
!                 write(*,101) axb(1, j, i)  
!                 write(*,101) debug_sum / ens
!                 write(*, *) " "
!             end if
! 101 FORMAT(F10.2\)
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

!    do i = 1, dx
!        do j = 1, dy
!            do var = 1, n
!!                axa(var, j, i) = 0
!                axa(var, j, i) = axb(var, j, i) + dot_product(xb(1:ens, var, j, i), overline_w)
!                do ens_num = 1, ens
!                    xa(ens_num, var, j, i) = axa(var, j, i) + dot_product(xb(1:ens, var, j, i), trans(1:ens, ens_num))
!!                    xa(ens_num, var, j, i) = axb(var, j, i) + dot_product(xb(1:ens, var, j, i), trans(1:ens, ens_num))
!!                    axa(var, j, i) = axa(var, j, i) + xa(ens_num, var, j, i)
!                end do
!                do ens_num = 1, ens
!                    xa(ens_num, var, j, i) = xa(ens_num, var, j, i) - axa(var, j, i)
!                end do
!            end do
!        end do
!    end do
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
