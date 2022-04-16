! Select observations for a single cube and generate yb and dis

! select_obs : select observations in obs_list

!!! assume one cube covers the whole k axis !!!

! j : lat index in process block
! i : lon index in process block

!             1  
!          ---------
!       .  |  *  . |   .
!       .  |  .  * |   .
!          ---------
!                2

!       *              
!       2  ---------
!          |       |   
!          |       |
!          ---------  1
!                     *

#define obs_add2(nobsl, direct, tindx, gc, yb, dis, rdiaginv, rloc) \
do t2 = 1, obs_nz; \
    nobsl = nobsl + 1; \
    yb(:, nobsl) = obs_list(direct)%pro_yb(:, t2, tindx); \
    dis(nobsl) = obs_list(direct)%dis_yo(t2, tindx); \
    rdiaginv(nobsl) = rinv(t2); \
    rloc(nobsl) = gc; \
end do

subroutine gen_yb_dis(sj, si, ej, ei, nobsl, yb, dis, rdiaginv, rloc)

    use pro_info, only : s_lat, e_lat, s_lon, e_lon, grid_lat, grid_lon
    use coupler_config, only : num_lat

    integer, intent(in) :: sj, si, ej, ei
    integer, intent(out) :: nobsl
    yb_type, intent(out) :: yb(ens, MAX_OBS)
    ayb_type,intent(out) :: dis(MAX_OBS)
    real(4), intent(out) :: rdiaginv(MAX_OBS)
    yb_type, intent(out) :: rloc(MAX_OBS)

    integer :: iter, tindx, direct
    ! real(8) :: cube_lat_s, cube_lat_e
    ! real(8) :: cube_lon_s, cube_lon_e
    integer :: cube_gslat, cube_gelat
    integer :: cube_gslon, cube_gelon
    integer :: lat1, lon1, lat2, lon2
    integer :: obs_indx
    real(8) :: lat1_v, lon1_v, lat2_v, lon2_v
    real(8) :: d1, d2
    real(8) :: gc, obs_dist
    integer :: t1, t2, t3
    integer :: temp, tlat, tlon, gtlat, gtlon
    integer :: center_j, center_i, center_lat, center_lon
    integer :: tmp_lat, tmp_lon, tmpsj, tmpsi, tmpej, tmpei
    integer :: i, j, k
    integer :: remote_id

    real(8) :: rci,rj,rk

    cube_gslat = sj
    cube_gslon = si
    cube_gelat = ej
    cube_gelon = ei
    center_j = (sj + ej) / 2
    center_i = (si + ei) / 2
    ! cube_lat_s = s_lat + (j - 1) * dlat
    ! cube_lon_s = s_lon + (i - 1) * dlon
    ! cube_lat_e = s_lat + (j + cube_n - 1) * dlat
    ! cube_lon_e = s_lon + (i + cube_n - 1) * dlon

    ! call obs_get_indx(cube_gslat, cube_gslon, lat1, lon1)
    call obs_get_indx(cube_gslat, cube_gslon, obs_lat_t, obs_lon_t, lat1, lon1)
    lat1 = lat1! - obs_1gox + 1
    lon1 = lon1! - obs_1goy + 1
    ! lat1 = int((cube_lat_s - obs_1lat) / obs_dlat) + 1 + 1
    ! lon1 = int((cube_lon_s - obs_1lon) / obs_dlon) + 1 + 1
    if (lat1 <= 0) lat1 = 1
    if (lat1 >= obs_list(0)%lat_n) lat1 = obs_list(0)%lat_n
    
    if (lon1 <= 0) lon1 = 1
    if (lon1 >= obs_list(0)%lon_n) lon1 = obs_list(0)%lon_n

    !call obs_get_indx(cube_gelat + 1, cube_gelon + 1, lat2, lon2)
    call obs_get_indx(cube_gelat, cube_gelon, obs_lat_t, obs_lon_t, lat2, lon2)
#if (DEBUG == 1)
!     if (pid == 0 .and. si ==1 .and. sj == 90) then
! !     if (pid == 0 .and. mod(si,40) == 0 .and. mod(sj,40) == 0) then
!         write(*, *) "cube", cube_gslat, cube_gelat, cube_gslon, cube_gelon
!         write(*, *) "obs", lat1, lon1, lat2, lon2
!         write(*, *) " "
!     endif
#endif
    lat2 = lat2! - 1! - obs_1gox + 1
    lon2 = lon2! - 1! - obs_1goy + 1
    ! lat2 = int((cube_lat_e - obs_1lat) / obs_dlat) + 1
    ! lon2 = int((cube_lon_e - obs_1lon) / obs_dlon) + 1
    ! if (lat2 <= 0) lat2 = 1
    ! if (lat2 >= obs_list(0)%lat_n) lat2 = obs_list(0)%lat_n

    ! if (lon2 <= 0) lon2 = 1
    ! if (lon2 >= obs_list(0)%lon_n) lon2 = obs_list(0)%lon_n

    ! lat1_v = obs_1lat + (lat1 - 1) * obs_dlat
    ! lon1_v = obs_1lon + (lon1 - 1) * obs_dlon
    ! lat2_v = obs_1lat + (lat2 - 1) * obs_dlat
    ! lon2_v = obs_1lon + (lon2 - 1) * obs_dlon

    ! there're only one observation in a cube
!    if (lat2 == lat1) then

        ! obs_indx = (lon1 - 1) * obs_list(0)%lat_n + lat1
        ! lat1_v = obs_list(0)%info(obs_indx)%lat

        ! d1 = lat1_v - grid_lat(sj-1)
        ! d2 = grid_lat(ej) - lat1_v
        ! if ((lat1 > 1 .and. d1 > d2) .or. lat1 == obs_list(0)%lat_n) then
        !     lat1 = lat1 - 1
        ! else
        !     lat2 = lat2 + 1
        ! end if
    ! there're no observation in a cube
!    else if (lat2 < lat1) then
        ! temp = lat1
        ! lat1 = lat2
        ! lat2 = temp
!    end if

    ! there're only one observation in a cube
!    if (lon2 == lon1) then

        ! obs_indx = (lon1 - 1) * obs_list(0)%lat_n + lat1
        ! lon1_v = obs_list(0)%info(obs_indx)%lon

        ! d1 = lon1_v - grid_lon(si-1)
        ! d2 = grid_lon(ei) - lon1_v
        ! ! d1 = lon1_v - cube_lon_s
        ! ! d2 = cube_lon_e - lon1_v
        ! if ((lon1 > 1 .and. d1 > d2) .or. lon1 == obs_list(0)%lon_n) then
        !     lon1 = lon1 - 1
        ! else
        !     lon2 = lon2 + 1
        ! end if
    ! there're no observation in a cube
!    else if (lon2 < lon1) then
        ! temp = lon1
        ! lon1 = lon2
        ! lon2 = temp
!    end if

    center_lat = (lat2 + lat1) / 2
    center_lon = (lon2 + lon1) / 2
    lat1 = center_lat; lat2 = center_lat
    lon1 = center_lon; lon2 = center_lon

    if (lat1 <= 0) lat1 = 1
    if (lat1 >= obs_list(0)%lat_n) lat1 = obs_list(0)%lat_n
    if (lon1 <= 0) lon1 = 1
    if (lon1 >= obs_list(0)%lon_n) lon1 = obs_list(0)%lon_n
    if (lat2 <= 0) lat2 = 1
    if (lat2 >= obs_list(0)%lat_n) lat2 = obs_list(0)%lat_n
    if (lon2 <= 0) lon2 = 1
    if (lon2 >= obs_list(0)%lon_n) lon2 = obs_list(0)%lon_n

    nobsl = 0
    ! select in yb_for_cpy
    if (OBS_DENSE >= 0) then
        do i = lon1 - OBS_DENSE, lon2 + OBS_DENSE
            do j = lat1 - OBS_DENSE, lat2 + OBS_DENSE
                tlat = j
                tlon = i
                gtlat = obs_1lat + (tlat - 1) * obs_dlat
                gtlon = obs_1lon + (tlon - 1) * obs_dlon
                tmp_lat = abs(center_j - gtlat)
                tmp_lon = abs(center_i - gtlon)
#if (SINGLE_POINT == 1)
                tmp_lon = min(tmp_lon, ny - gtlon + center_i, ny - center_i + gtlon)
#endif
                obs_dist = sqrt(tmp_lat**2*1d0 + (tmp_lon*cos((grid_lat(center_j)+(gtlat - center_j)/2*180/num_lat)/180d0*pi))**2*1d0)
                if (gc_enable == .true.) then
                    gc = schr2(obs_dist * 1d0 / da_local_distance)
                else
                    gc = schr(obs_dist * 1d0 / da_local_distance)
                end if

                do k = 1, obs_nz
                    nobsl = nobsl + 1
                    yb(:, nobsl) = yb_for_cpy(:, k, j, i)
                    dis(nobsl) = ayb_for_cpy(k, j, i)
                    rdiaginv(nobsl) = rinv(k)
                    rloc(nobsl) = gc
                end do
            end do
        end do
    end if
!     if (OBS_DENSE >= 0) then
!         do i = lon1 - OBS_DENSE, lon2 + OBS_DENSE
!             do j = lat1 - OBS_DENSE, lat2 + OBS_DENSE
!                 tlat = j
!                 tlon = i
!                 gtlat = obs_1lat + (tlat - 1) * obs_dlat
!                 gtlon = obs_1lon + (tlon - 1) * obs_dlon
!                 tmp_lat = abs(center_j - gtlat)
!                 tmp_lon = abs(center_i - gtlon)
!#if (SINGLE_POINT == 1)
!                 tmp_lon = min(tmp_lon, ny - gtlon + center_i, ny - center_i + gtlon)
!#endif
!                 if ((i >= 1 - obs_list(2)%lon_n .and. i <= obs_list(0)%lon_n + obs_list(1)%lon_n .and. &
!                     j >= 1 .and. j <= obs_list(0)%lat_n) .or. &
!                     (j >= 1 - obs_list(4)%lat_n .and. j <= obs_list(0)%lat_n + obs_list(3)%lat_n .and. &
!                     i >= 1 .and. i <= obs_list(0)%lon_n)) then
!                     !if (i .ne. lon1 .and. i .ne. lon2 .and. j .ne. lat1 .and. j .ne. lat2) then
!#if (DEBUG == 1)
! !     if (pid == 0 .and. si ==20 .and. sj == 11) then
! !         write(*, *) "dense", OBS_DENSE
! !         write(*, *) "cube", lat1, lon1, lat2, lon2
! !         write(*, *) "gtlat, gtlon, center_j, center_i ", gtlat, gtlon, center_j, center_i
! !         write(*, *) " "
! !     endif
!#endif
!                     direct = 0
!                     if (tlon > obs_list(0)%lon_n) then
!                         tlon = tlon - obs_list(0)%lon_n
!                         direct = 1
!                     endif
!                     if (tlon <= 0) then
!                         tlon = tlon + obs_list(2)%lon_n
!                         direct = 2
!                     end if
!                     if (tlat > obs_list(0)%lat_n) then
!                         tlat = tlat - obs_list(0)%lat_n
!                         direct = 3
!                     end if
!                     if (tlat <= 0) then
!                         tlat = tlat + obs_list(4)%lat_n
!                         direct = 4
!                     end if
!                     tindx = tlat + (tlon - 1) * obs_list(direct)%lat_n
!                     obs_dist = sqrt(tmp_lat**2*1d0 + (tmp_lon*cos((grid_lat(center_j)+(gtlat - center_j)/2*180/num_lat)/180d0*pi))**2*1d0)
!                     if (gc_enable == .true.) then
!                         gc = schr2(obs_dist * 1d0 / da_local_distance)
!                     else
!                         gc = schr(obs_dist * 1d0 / da_local_distance)
!                     end if
!                     obs_add2(nobsl, direct, tindx, gc, yb, dis, rdiaginv, rloc)
!                     if (ayb_for_cpy(1, j, i) .ne. obs_list(direct)%dis_yo(1, tindx)) then
!                         print *, pid, "error at", direct, j, i
!                     end if
!                 else if (neighbour_halo == .true.) then
!                     if (i < 1 .and. j < 1) then
!                         direct = 8
!                         tlat = tlat + obs_list(8)%lat_n
!                         tlon = tlon + obs_list(8)%lon_n
!                     end if
!                     if (i > obs_list(0)%lon_n .and. j < 1) then
!                         direct = 6
!                         tlat = tlat + obs_list(6)%lat_n
!                         tlon = tlon - obs_list(0)%lon_n
!                     end if
!                     if (i < 1 .and. j > obs_list(0)%lat_n) then
!                         direct = 7
!                         tlat = tlat - obs_list(0)%lat_n
!                         tlon = tlon + obs_list(7)%lon_n
!                     end if
!                     if (i > obs_list(0)%lon_n .and. j > obs_list(0)%lat_n) then
!                         direct = 9
!                         tlat = tlat - obs_list(0)%lat_n
!                         tlon = tlon - obs_list(0)%lon_n
!                     end if
!                     if (tlat >= 1 .and. tlat <= obs_list(direct)%lat_n .and. &
!                     tlon >= 1 .and. tlon <= obs_list(direct)%lon_n) then
!                         !if (pid == 0 .and. ei == ny .and. ej == nx) print *, "chosen hiahia!!!!!", j, i, direct, tlat, tlon
!                         !if (pid == 5 .and. si == 1 .and. sj == 1) print *, "chosen hiahia!!!!!", j, i, direct, tlat, tlon
!                         tindx = tlat + (tlon - 1) * obs_list(direct)%lat_n
!                         obs_dist = sqrt(tmp_lat**2*1d0 + (tmp_lon*cos((grid_lat(center_j)+(gtlat - center_j)/2*180/num_lat)/180d0*pi))**2*1d0)
!                         if (gc_enable == .true.) then
!                             gc = schr2(obs_dist * 1d0 / da_local_distance)
!                         else
!                             gc = schr(obs_dist * 1d0 / da_local_distance)
!                         end if
!                         obs_add2(nobsl, direct, tindx, gc, yb, dis, rdiaginv, rloc)                            
!                         if (ayb_for_cpy(1, j, i) .ne. obs_list(direct)%dis_yo(1, tindx)) then
!                             print *, pid, "error at", direct, j, i
!                         end if
!                     end if
!                 end if
!             end do
!         end do
!     end if

    do iter = 0, ndist
        if (local_dist(iter) > OBS_DENSE) then
            ! east, lat1, lon2
            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            tlat = lat1
            tlon = lon2 + local_dist(iter)
            gtlat = obs_1lat + (tlat - 1) * obs_dlat
            gtlon = obs_1lon + (tlon - 1) * obs_dlon
            tmp_lat = abs(center_j - gtlat)
            tmp_lon = abs(center_i - gtlon)
            obs_dist = sqrt(tmp_lat**2*1d0 + (tmp_lon*cos((grid_lat(center_j)+(gtlat - center_j)/2*180/num_lat)/180d0*pi))**2*1d0)
#if (DEBUG == 1)
!     if (pid == 0 .and. si ==20 .and. sj == 11) then
!         write(*, *) "gtlat, gtlon, center_j, center_i ", gtlat, gtlon, center_j, center_i
!         write(*, *) " "
!     endif
#endif
            if (gc_enable == .true.) then
                gc = schr2(obs_dist*1d0 /da_local_distance)
            else
                gc = schr(obs_dist*1d0 /da_local_distance)
            end if
            if (tlon <= obs_list(0)%lon_n + obs_list(1)%lon_n) then
                direct = 0
                if (tlon > obs_list(0)%lon_n) then
                    tlon = tlon - obs_list(0)%lon_n
                    direct = 1
                endif
                tindx = tlat + (tlon - 1) * obs_list(direct)%lat_n
                if (gc_enable == .true.) then
                    gc = schr2(obs_dist*1d0 /da_local_distance)
                else
                    gc = schr(obs_dist*1d0 /da_local_distance)
                end if
                obs_add2(nobsl, direct, tindx, gc, yb, dis, rdiaginv, rloc)
            else if (remote_comm == 1) then
                do i = 1, rns(1)
                    if (center_lon + local_dist(iter) >= r_sn(i, 1) .and. center_lon + local_dist(iter) <= r_en(i, 1)) then
                        direct = 5
                        do j = 1, remote_obsn
                            tindx = (i - 1) * remote_obsn + j
                            !if (pid == 0) print *, "hahaha", pid, tindx, obs_list(5)%lat_n
                            obs_add2(nobsl, direct, tindx, gc, yb, dis, rdiaginv, rloc)
                        end do
                    end if
                end do
            end if
            ! west, lat2, lon1
            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            tlat = lat2
            tlon = lon1 - local_dist(iter)
            gtlat = obs_1lat + (tlat - 1) * obs_dlat
            gtlon = obs_1lon + (tlon - 1) * obs_dlon
            tmp_lat = abs(center_j - gtlat)
            tmp_lon = abs(center_i - gtlon)
            obs_dist = sqrt(tmp_lat**2*1d0 + (tmp_lon*cos((grid_lat(center_j)+(gtlat - center_j)/2*180/num_lat)/180d0*pi))**2*1d0)
#if (DEBUG == 1)
!     if (pid == 0 .and. si ==20 .and. sj == 11) then
!         write(*, *) "gtlat, gtlon, center_j, center_i ", gtlat, gtlon, center_j, center_i
!         write(*, *) " "
!     endif
#endif
            if (gc_enable == .true.) then
                gc = schr2(obs_dist*1d0 /da_local_distance)
            else
                gc = schr(obs_dist*1d0 /da_local_distance)
            end if
            if (tlon > -obs_list(2)%lon_n) then
                direct = 0
                if (tlon <= 0) then
                    tlon = tlon + obs_list(2)%lon_n
                    direct = 2
                end if
                if (gc_enable == .true.) then
                    gc = schr2(obs_dist*1d0 /da_local_distance)
                else
                    gc = schr(obs_dist*1d0 /da_local_distance)
                end if
                tindx = tlat + (tlon - 1) * obs_list(direct)%lat_n
                obs_add2(nobsl, direct, tindx, gc, yb, dis, rdiaginv, rloc)
            else if (remote_comm == 1) then
                do i = 1, rns(2)
                    if (obs_list(0)%lon_n - center_lon + 1 + local_dist(iter) >= r_sn(i, 2) .and. &
                    obs_list(0)%lon_n - center_lon + 1 + local_dist(iter) <= r_en(i, 2)) then
                        direct = 5
                        do j = 1, remote_obsn
                            tindx = (rns(1) + i - 1) * remote_obsn + j
                            obs_add2(nobsl, direct, tindx, gc, yb, dis, rdiaginv, rloc)
                        end do
                    end if
                end do
            end if
            ! south, lat2, lon2
            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            tlat = lat2 + local_dist(iter)
            tlon = lon2
            gtlat = obs_1lat + (tlat - 1) * obs_dlat
            gtlon = obs_1lon + (tlon - 1) * obs_dlon
            tmp_lat = abs(center_j - gtlat)
            tmp_lon = abs(center_i - gtlon)
            obs_dist = sqrt(tmp_lat**2*1d0 + (tmp_lon*cos((grid_lat(center_j)+(gtlat - center_j)/2*180/num_lat)/180d0*pi))**2*1d0)
#if (DEBUG == 1)
!     if (pid == 0 .and. si ==20 .and. sj == 11) then
!         write(*, *) "gtlat, gtlon, center_j, center_i ", gtlat, gtlon, center_j, center_i
!         write(*, *) " "
!     endif
#endif
            if (gc_enable == .true.) then
                gc = schr2(obs_dist*1d0 /da_local_distance)
            else
                gc = schr(obs_dist*1d0 /da_local_distance)
            end if
            if (tlat <= obs_list(0)%lat_n + obs_list(3)%lat_n) then
                direct = 0
                if (tlat > obs_list(0)%lat_n) then
                    tlat = tlat - obs_list(0)%lat_n
                    direct = 3
                end if
                if (gc_enable == .true.) then
                    gc = schr2(obs_dist*1d0 /da_local_distance)
                else
                    gc = schr(obs_dist*1d0 /da_local_distance)
                end if
                tindx = tlat + (tlon - 1) * obs_list(direct)%lat_n
                obs_add2(nobsl, direct, tindx, gc, yb, dis, rdiaginv, rloc)
            else if (remote_comm == 1) then
                do i = 1, rns(3)
                    if (center_lat + local_dist(iter) >= r_sn(i, 3) .and. &
                    center_lat + local_dist(iter) <= r_en(i, 3)) then
                        direct = 5
                        do j = 1, remote_obsn
                            tindx = (rns(1) + rns(2) + i - 1) * remote_obsn + j
                            obs_add2(nobsl, direct, tindx, gc, yb, dis, rdiaginv, rloc)
                        end do
                    end if
                end do
            end if
            ! north, lat1, lon1
            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            tlat = lat1 - local_dist(iter)
            tlon = lon1
            gtlat = obs_1lat + (tlat - 1) * obs_dlat
            gtlon = obs_1lon + (tlon - 1) * obs_dlon
            tmp_lat = abs(center_j - gtlat)
            tmp_lon = abs(center_i - gtlon)
            obs_dist = sqrt(tmp_lat**2*1d0 + (tmp_lon*cos((grid_lat(center_j)+(gtlat - center_j)/2*180/num_lat)/180d0*pi))**2*1d0)
#if (DEBUG == 1)
!     if (pid == 0 .and. si ==20 .and. sj == 11) then
!         write(*, *) "gtlat, gtlon, center_j, center_i ", gtlat, gtlon, center_j, center_i
!         write(*, *) " "
!     endif
#endif
            if (gc_enable == .true.) then
                gc = schr2(obs_dist*1d0 /da_local_distance)
            else
                gc = schr(obs_dist*1d0 /da_local_distance)
            end if
            if (tlat > -obs_list(4)%lat_n) then
                direct = 0
                if (tlat <= 0) then
                    tlat = tlat + obs_list(4)%lat_n
                    direct = 4
                end if
                if (gc_enable == .true.) then
                    gc = schr2(obs_dist*1d0 /da_local_distance)
                else
                    gc = schr(obs_dist*1d0 /da_local_distance)
                end if
                tindx = tlat + (tlon - 1) * obs_list(direct)%lat_n
                obs_add2(nobsl, direct, tindx, gc, yb, dis, rdiaginv, rloc)
            else if (remote_comm == 1) then
                do i = 1, rns(4)
                    if (obs_list(0)%lat_n - center_lat + 1 + local_dist(iter) >= r_sn(i, 4) .and. &
                    obs_list(0)%lat_n - center_lat + 1 + local_dist(iter) <= r_en(i, 4)) then
                        direct = 5
                        do j = 1, remote_obsn
                            tindx = (rns(1) + rns(2) + rns(3) + i - 1) * remote_obsn + j
                            obs_add2(nobsl, direct, tindx, gc, yb, dis, rdiaginv, rloc)
                        end do
                    end if
                end do
            end if
        end if
    end do

    ! select add the observations in obs_list(5)
    ! if (remote_comm .and. obs_list(5)%lat_n >= 0) then
    !     do i = 1, 4
    !         if (i == 4) then
    !             tmpsj = sj; tmpsi = si
    !             tmpej = ej; tmpei = ei
    !         end if
    !         k = 1
    !         do j = 1, dist
    !             if (tmp )
    !             do while (r_sn(k) .or. r_en(k))

    !             end do
    !         end do
    !     end do
    !     ! do iter = 1, obs_list(5)%lat_n
    !     !     tindx = iter
    !     !     direct = 5
    !     !     call obs_add(nobsl, direct, tindx, gc, yb, dis, rdiaginv, rloc)
    !     ! end do
    ! end if

end subroutine

subroutine obs_add(nobsl, direct, tindx, gc, yb, dis, rdiaginv, rloc)
    integer, intent(inout) :: nobsl
    integer, intent(in) :: direct
    integer, intent(in) :: tindx
    real(8), intent(in)  :: gc
    yb_type, intent(inout) :: yb(ens, MAX_OBS)
    ayb_type,intent(inout) :: dis(MAX_OBS)
    real(4), intent(inout) :: rdiaginv(MAX_OBS)
    yb_type, intent(inout) :: rloc(MAX_OBS)

    integer :: t1, t2

    !do t1 = 1, obs_l
        do t2 = 1, obs_nz
            nobsl = nobsl + 1
            yb(:, nobsl) = obs_list(direct)%pro_yb(:, t2, tindx)
            dis(nobsl) = obs_list(direct)%dis_yo(t2, tindx)
            rdiaginv(nobsl) = rinv(t2)
            rloc(nobsl) = gc
            ! if (isnan(dis(nobsl))) then
            !     print *, nobsl, "error!"
            ! end if
        end do
    !end do
#if (DEBUG == 1)
    ! do t2 = 1, ens
    !     if (isnan(yb(t2, nobsl))) then
    !         print *, pid, direct, nobsl, "error1"
    !     end if
    ! end do
#endif

end subroutine

real(8) function schr2(ri)

    implicit none

    real(8) :: rci,ri,rj,rk

    if (ri.gt.2.0) then
        rci=0.0
    else if (ri.gt.1.0) then
        rj=ri/12.0-0.5
        rj=rj*ri+0.625
        rj=rj*ri+1.666667
        rci=(rj*ri-5.0)*ri+4.0-0.6666667/ri
    else if (ri.gt.0.0) then
        rk=ri*ri
        rj=(-0.25*ri+0.5)*ri+0.625
        rci=(rj*ri-1.666667)*rk+1.0
    else
        rci=1.0
    endif

    schr2=rci

    return
end

real(8) function schr(ri)

    implicit none

    real(8) :: rci,ri,rj,rk

    if (ri.gt.2.0) then
        rci=0.0
    else
        rci=1.0
    endif

    schr=rci

    return
end