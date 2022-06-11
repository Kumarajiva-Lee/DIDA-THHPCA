! observation make

! subroutine obs_get_ndist(ndist)
!     use da_namelist_mod, only:da_obs_sparse
!     integer, intent(out) :: ndist
!     if (da_obs_sparse == 1) then
!         ndist = floor(2.0*da_local_distance/obs_dlon)
!     elseif (da_obs_sparse > 1) then
!         ndist = ceiling(log(2.0*da_local_distance/obs_dlon/da_obs_sparse) / log(2.0))+2.0
!     endif
! !   print *, 'ndist=',ndist
! end subroutine

! subroutine obs_get_local_dist(ndist, local_dist)
!     use da_namelist_mod, only:da_obs_sparse

!     integer, intent(in) :: ndist
!     integer, intent(out):: local_dist(0: ndist)
    
!     integer :: i
!     local_dist(0) = 0
!     do i = 1, ndist
!         if (da_obs_sparse == 1) then
!             local_dist(i) = i
!         elseif (da_obs_sparse > 1) then
!             local_dist(i) = 2.0 ** (i - 2.0) * da_obs_sparse
!             if (i == ndist) local_dist(i) = floor(2.0*da_local_distance/obs_dlon)
!         endif
!         if (i==1) local_dist(i) = i
!     end do
! !   print *,'local_dist',local_dist
! end subroutine

! subroutine obs_get_first(obs_lat_val, obs_lon_val, obs_1lat, obs_1lon)

!     real(4), intent(in) :: obs_lat_val, obs_lon_val
!     integer, intent(out) :: obs_1lat, obs_1lon

!     integer :: i

!     obs_1lat = 0
!     obs_1lon = 0

!     do i = 1, nx
!         if (obs_lat_val < grid_lat(i)) then
!             obs_1lat = i
!             exit
!         end if
!     end do
!     do i = 1, ny
!         if (obs_lon_val < grid_lon(i)) then
!             obs_1lon = i
!             exit
!         end if
!     end do

! #if (DEBUG == 1)
!     !print *, pid, "obslatlon", obs_lat_val, obs_lon_val
!     !print *, pid, "grid_latlon", grid_lat(nx), grid_lon(ny)
!     if (obs_1lon == 0 .or. obs_1lat == 0) then
!         call log_warning("obs get first warning!")
!     end if
! #endif
! end subroutine

! ! giving index grid point(start from 1), calculate the nearest observation
! subroutine obs_get_indx(lat_indx, lon_indx, obs_lat, obs_lon, obs_lat_indx, obs_lon_indx)

!     integer, intent(in) :: lat_indx, lon_indx
!     real, intent(in) :: obs_lat(obs_list(0)%lat_n), obs_lon(obs_list(0)%lon_n)
!     integer, intent(out):: obs_lat_indx, obs_lon_indx

!     integer :: i

!     obs_lat_indx = obs_list(0)%lat_n + 1
!     do i = 1, obs_list(0)%lat_n
!         if (obs_lat(i) > grid_lat(lat_indx)) then
!             obs_lat_indx = i
!             exit
!         end if
!     end do
!     !if (obs_lat(i) == 0) grid_lat = 1
!     obs_lon_indx = obs_list(0)%lon_n + 1
!     do i = 1,  obs_list(0)%lon_n
!         if (obs_lon(i) > grid_lon(lon_indx)) then
!             obs_lon_indx = i
!             exit
!         end if
!     end do
!     !if (obs_lon(i) == 0) grid_lon = 1
!     ! obs_lat_indx = minloc(abs(grid_lat(lat_indx) - obs_lat))
!     ! obs_lon_indx = minloc(abs(grid_lon(lon_indx) - obs_lon))
!     ! if (grid_lat(lat_indx) - obs_lat(obs_lat_indx) < 0) obs_lat_indx = obs_lat_indx - 1
!     ! if (grid_lon(lon_indx) - obs_lon(obs_lon_indx) > 0) obs_lon_indx = obs_lon_indx + 1

! end subroutine

! !!! Notice : obs_interpolation is different from obs_interpolation_3d
! !!!          obs_interpolation does it on a single level
! !!!          while obs_interpolation_3d does it on a whole column

! subroutine obs_interpolation(ens, n, xb, axb, yb, lat1, lon1, &
!     lat3, lon3, lat2, lon2)
!     integer, intent(in) :: ens, n
!     xb_type, intent(in) :: xb(1:ens, 1:n, 2, 2)
!     axb_type,intent(in) :: axb(1:n, 2, 2)
!     yb_type, intent(out):: yb(1:ens, 1:n, 1)
!     real(8), intent(in) :: lat1, lon1, lat3, lon3
!     real(8), intent(in) :: lat2, lon2
    
!     real(8) :: dlat, dlon

!     integer :: i, j
!     real(8) :: t1, t2
!     yb_type :: r1(ens, n), r2(ens, n)

!     dlat = lat3 - lat1
!     dlon = lon3 - lon1

!     t1 = (dlat - lat2 + lat1) / dlat
!     t2 = (dlon - lon2 + lon1) / dlon

!     do i = 1, n
!         do j = 1, ens
!             r1(j, i) = (xb(j, i, 1, 1) + axb(i, 1, 1)) * t1 + &
!             (xb(j, i, 2, 1) + axb(i, 2, 1)) * (1 - t1)
!             r2(j, i) = (xb(j, i, 1, 2) + axb(i, 1, 2)) * t1 + &
!             (xb(j, i, 2, 2) + axb(i, 2, 2)) * (1 - t1)
!             yb(j, i, 1) = r1(j, i) * t2 + r2(j, i) * (1 - t2)
!         end do
!     end do

! end subroutine

! subroutine obs_interpolation_3d(ens, n, nz, obs_l, xb, axb, ph, obs_ph, yb, &
!     t1, t2)
!     integer, intent(in) :: ens, n, nz, obs_l
!     xb_type, intent(in) :: xb(1:ens, 1:n, 2, 2)
!     axb_type,intent(in) :: axb(1:n, 2, 2)
!     ayb_type, intent(in):: ph(ens, 1:nz+1)
!     ayb_type, intent(in):: obs_ph(1:obs_l)
!     yb_type, intent(out):: yb(1:ens, 1:obs_nz, 1)
!     real(8), intent(in) :: t1, t2

!     yb_type :: r1, r2, r3, r4

!     integer :: i, j, k, l

!     integer :: ph_pos(ens, obs_l)

!     real(8) :: dlat, dlon

!     integer :: tmp, tl
!     real(8) :: t3

!     ! dlat = lat3 - lat1
!     ! dlon = lon3 - lon1

!     ! t1 = (dlat - lat2 + lat1) / dlat
!     ! t2 = (dlon - lon2 + lon1) / dlon

!     !if (pid == 0) print *, "obs_ph", obs_ph

!     do k = 1, ens
!         j = 1
!         do i = 1, obs_l
!             do while (j <= nz .and. obs_ph(i) - ph(k, j) > 0.00001d0)
!                 j = j + 1
!             end do
!             !if (pid == 0) print *, "obs_ph", obs_ph(i), j
!             ph_pos(k, i) = j - 1
!             !if (ph_pos(k, i) == nz) ph_pos(k, i) = nz - 1
! #if (DEBUG == 1)
!             ! if (j > nz) then
!             !     call log_notice("Outside interpolation!")
!             !     stop
!             ! end if
! #endif
!         end do
!     end do

!     tmp = 0
!     do k = 1, obs_n
!         if (obs_var_lev(k) .ne. 1) then
!             do l = 1, obs_l
!                 tmp = tmp + 1
!                 do i = 1, ens
!                     if (ph_pos(i, l) >= 1 .and. ph_pos(i, l) < nz) then
!                             tl = obs_var_start(k) + ph_pos(i, l) - 1
!                             t3 = (log(ph(i, ph_pos(i, l)+1)) - log(obs_ph(l))) / (log(ph(i, ph_pos(i, l)+1)) - log(ph(i, ph_pos(i, l))))
!                             r1 = (xb(i, tl, 1, 1) + axb(tl, 1, 1)) * t1 + &
!                             (xb(i, tl, 2, 1) + axb(tl, 2, 1)) * (1 - t1)
!                             r2 = (xb(i, tl, 1, 2) + axb(tl, 1, 2)) * t1 + &
!                             (xb(i, tl, 2, 2) + axb(tl, 2, 2)) * (1 - t1)
!                             r3 = (xb(i, tl+1, 1, 1) + axb(tl+1, 1, 1)) * t1 + &
!                             (xb(i, tl+1, 2, 1) + axb(tl+1, 2, 1)) * (1 - t1)
!                             r4 = (xb(i, tl+1, 1, 2) + axb(tl+1, 1, 2)) * t1 + &
!                             (xb(i, tl+1, 2, 2) + axb(tl+1, 2, 2)) * (1 - t1)
!                             r1 = r1 * t2 + r2 * (1 - t2)
!                             r3 = r3 * t2 + r4 * (1 - t2)
!                             yb(i, tmp, 1) = r1 * t3 + r3 * (1 - t3)
!                         !if (pid == 0) print *, "ph_pos", ph_pos(i, l), &
!                         !"l", l, "tl", tl, "yb", yb(:, tmp, 1)
!                     else
!                         if (ph_pos(i, l) == 0) then
!                             tl = obs_var_start(k)
!                         else
!                             tl = obs_var_start(k) + nz - 1
!                         end if
!                         r1 = (xb(i, tl, 1, 1) + axb(tl, 1, 1)) * t1 + &
!                         (xb(i, tl, 2, 1) + axb(tl, 2, 1)) * (1 - t1)
!                         r2 = (xb(i, tl, 1, 2) + axb(tl, 1, 2)) * t1 + &
!                         (xb(i, tl, 2, 2) + axb(tl, 2, 2)) * (1 - t1)
!                         yb(i, tmp, 1) = r1 * t2 + r2 * (1 - t2)
!                     end if
!                 end do
!             end do
!         else
!             tl = obs_var_start(k)
!             tmp = tmp + 1
!             do i = 1, ens
!                 r1 = (xb(i, tl, 1, 1) + axb(tl, 1, 1)) * t1 + &
!                 (xb(i, tl, 2, 1) + axb(tl, 2, 1)) * (1 - t1)
!                 r2 = (xb(i, tl, 1, 2) + axb(tl, 1, 2)) * t1 + &
!                 (xb(i, tl, 2, 2) + axb(tl, 2, 2)) * (1 - t1)
!                 yb(i, tmp, 1) = r1 * t2 + r2 * (1 - t2)
!             end do
!         end if
!     end do

! end subroutine

subroutine obs_interpolation_3d_2(xb, axb, ph, obs_d, obs_vid, obs_i, t1, t2)

    xb_type, intent(in) :: xb(1:ens, 1:n, 2, 2)
    axb_type,intent(in) :: axb(1:n, 2, 2)

    ayb_type, intent(in):: ph(ens, 1:nz+1)

    type(obs_list), intent(inout) :: obs_d
    integer, intent(in) :: obs_vid(n_var)
    integer, intent(in) :: obs_i
    real(8), intent(in) :: t1, t2

    yb_type :: r1, r2, r3, r4
    real(8) :: t3

    integer :: i, j, k, l, tl

    integer :: ph_pos(ens, obs_d%obs_nz)

    do k = 1, ens
        j = 1
        do i = 1, obs_d%obs_nz
            do while (j <= nz .and. obs_d%ph(i, obs_i) > ph(k, j))
                j = j + 1
            end do
            ph_pos(k, i) = j - 1
#if (DEBUG == 1)
            if (j > nz) then
                ! if (pid == 0) then
                    ! print *, "LETKF haha", obs_d%ph(i, obs_i)
                    ! do j = 1, nz
                    !     print 123, ph(k, j)
                    ! end do
                    ! print *, " "
                ! end if
                !call log_error("Outside interpolation!")
            end if
#endif
        end do
    end do

    123 FORMAT(F20.3, $) 

    !tmp = 0
    do l = 1, obs_d%obs_nz
        do k = 1, obs_d%obs_nvar
            if (da_var_lev(obs_vid(k)) > 1) then
                do i = 1, ens
                    if (ph_pos(i, l) >= 1 .and. ph_pos(i, l) < nz) then
                        tl = da_var_start(obs_vid(k)) + ph_pos(i, l) - 1
                        t3 = (log(ph(i, ph_pos(i, l)+1)) - log(obs_d%ph(l, obs_i)) ) &
                        / (log(ph(i, ph_pos(i, l)+1)) - log(ph(i, ph_pos(i, l))))
                        r1 = (xb(i, tl, 1, 1) + axb(tl, 1, 1)) * t1 + &
                        (xb(i, tl, 2, 1) + axb(tl, 2, 1)) * (1 - t1)
                        r2 = (xb(i, tl, 1, 2) + axb(tl, 1, 2)) * t1 + &
                        (xb(i, tl, 2, 2) + axb(tl, 2, 2)) * (1 - t1)
                        r3 = (xb(i, tl+1, 1, 1) + axb(tl+1, 1, 1)) * t1 + &
                        (xb(i, tl+1, 2, 1) + axb(tl+1, 2, 1)) * (1 - t1)
                        r4 = (xb(i, tl+1, 1, 2) + axb(tl+1, 1, 2)) * t1 + &
                        (xb(i, tl+1, 2, 2) + axb(tl+1, 2, 2)) * (1 - t1)
                        r1 = r1 * t2 + r2 * (1 - t2)
                        r3 = r3 * t2 + r4 * (1 - t2)
                        obs_d%pro_yb(i, k, l, obs_i) = r1 * t3 + r3 * (1 - t3)
                        ! if (pid == 0 .and. obs_d%pro_yb(i, k, l, obs_i) == 0.0) then
                        !     print *, "LETKF id: ", pid, obs_vid(k), tl
                        !     print *, t1, t2, t3
                        !     print *, xb(i, :, 1: 2, 1: 2)
                        !     call log_error("0 obs")
                        ! end if
                    else
                        if (ph_pos(i, l) == 0) then
                            tl = da_var_start(obs_vid(k))
                        else
                            tl = da_var_start(obs_vid(k)) + nz - 1
                        end if
                        r1 = (xb(i, tl, 1, 1) + axb(tl, 1, 1)) * t1 + &
                        (xb(i, tl, 2, 1) + axb(tl, 2, 1)) * (1 - t1)
                        r2 = (xb(i, tl, 1, 2) + axb(tl, 1, 2)) * t1 + &
                        (xb(i, tl, 2, 2) + axb(tl, 2, 2)) * (1 - t1)
                        obs_d%pro_yb(i, k, l, obs_i) = r1 * t2 + r2 * (1 - t2)
                        ! if (pid == 0 .and. obs_d%pro_yb(i, k, l, obs_i) == 0.0) then
                        !     print *, "LETKF id: ", pid, obs_vid(k), tl, "pos", ph_pos(i, l), "n", n
                        !     print *, xb(1, 1, 1:2, 1:2)
                        !     print *, xb(1, 2:33, 1:2, 1:2)
                        !     print *, xb(1, 33:55, 1:2, 1:2)
                        !     print *, xb(1, 56:97, 1:2, 1:2)
                        !     call log_error("0 obs2")
                        ! end if
                    end if
                end do
            else
                tl = da_var_start(obs_vid(k))
                !tmp = tmp + 1
                do i = 1, ens
                    r1 = (xb(i, tl, 1, 1) + axb(tl, 1, 1)) * t1 + &
                    (xb(i, tl, 2, 1) + axb(tl, 2, 1)) * (1 - t1)
                    r2 = (xb(i, tl, 1, 2) + axb(tl, 1, 2)) * t1 + &
                    (xb(i, tl, 2, 2) + axb(tl, 2, 2)) * (1 - t1)
                    obs_d%pro_yb(i, k, l, obs_i) = r1 * t2 + r2 * (1 - t2)

                    ! if (pid == 0 .and. obs_d%pro_yb(i, k, l, obs_i) == 0.0) then
                    !     print *, "LETKF id: ", pid, obs_vid(k), tl
                    !     call log_error("0 obs3")
                    ! end if
                end do
            end if
        end do
    end do

end subroutine


subroutine ps2ph_3d_interp(ps, aps, ph, t1, t2)

    use da_vert_coord_mod

    xb_type, intent(in) :: ps(ens, 1, 2, 2)
    axb_type,intent(in) :: aps(1, 2, 2)
    ayb_type, intent(out) :: ph(ens, nz+1)
    real(8), intent(in) :: t1, t2

    axb_type :: ps_tmp(ens)
    axb_type :: r1, r2
    axb_type :: ph_lev(ens, nz+1)

    integer :: i, j

    do i = 1, ens
        r1 = (ps(i, 1, 1, 1) + aps(1, 1, 1))* t1 + &
        (ps(i, 1, 2, 1) + aps(1, 2, 1)) * (1 - t1)
        r2 = (ps(i, 1, 1, 2) + aps(1, 1, 2))* t1 + &
        (ps(i, 1, 2, 2) + aps(1, 2, 2)) * (1 - t1)
        ps_tmp(i) = r1 * t2 + r2 * (1 - t2)
    end do
    !if (pid == 0) print *, pid, "haha, begin calc_ph_lev"
    do i = 1, nz + 1
        do j = 1, ens
            ph_lev(j, i) = da_vert_coord_calc_ph_lev(i, ps_tmp(j))
        end do
    end do
    !if (pid == 0) print *, pid, "haha, done calc_ph_lev"
    do i = 1, nz
        do j = 1, ens
            ph(j, i) = 0.5d0 * (ph_lev(j, i) + ph_lev(j, i+1))
        end do
    end do
    ph(:, nz+1) = ph(:, nz)

end subroutine

! subroutine ps2lnp(ps, aps, lnp)

!     use da_vert_coord_mod

!     xb_type, intent(in) :: ps(ens, nz)
!     axb_type,intent(in) :: aps(nz)
!     ayb_type, intent(out) :: lnp(ens, nz)

!     axb_type :: ps_tmp(ens)

!     axb_type :: ph_lev(ens, nz+1)

!     integer :: i, j

!     do i = 1, nz
!         do j = 1, ens
!             ps_tmp(i) = ps(j, i) + aps(i)
!         end do
!     end do

!     do i = 1, nz + 1
!         do j = 1, ens
!             ph_lev(j, i) = da_vert_coord_calc_ph_lev(i, ps_tmp(j))
!         end do
!     end do
!     do i = 1, nz
!         do j = 1, ens
!             lnp(j, i) = 0.5d0 * (ph_lev(j, i) + ph_lev(j, i+1))
!             lnp(j, i) = log(lnp(j, i))
!         end do
!     end do
!     !lnp(:, nz+1) = lnp(:, nz)

! end subroutine

subroutine ps2alnp(aps, alnp)

    use da_vert_coord_mod

    axb_type, intent(in) :: aps
    ayb_type, intent(out) :: alnp(nz)

    axb_type :: ph_lev(nz+1)

    integer :: i

    do i = 1, nz
        ph_lev(i) = da_vert_coord_calc_ph_lev(i, aps)
    end do
    do i = 1, nz
        alnp(i) = 0.5d0 * (ph_lev(i) + ph_lev(i+1))
        alnp(i) = log(alnp(i))
    end do

end subroutine
