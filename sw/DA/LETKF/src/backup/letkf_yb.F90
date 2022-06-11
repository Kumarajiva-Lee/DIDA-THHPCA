! generate yb & ayb matrix

subroutine gen_yb_dis(sj, si, ej, ei, nobsl, yb, dis, rdiaginv, rloc)

    integer, intent(in) :: sj, si, ej, ei
    integer, intent(out) :: nobsl

    yb_type, intent(out) :: yb(ens, MAX_OBS)
    ayb_type,intent(out) :: dis(MAX_OBS)
    real(4), intent(out) :: rdiaginv(MAX_OBS)
    yb_type, intent(out) :: rloc(MAX_OBS)

    real(8) :: center_lat, center_lon
    real(8) :: obs_dist
    real(8) :: gc

    integer :: i, j, k, l, m, n

    center_lat = grid_lat((sj + ej) / 2)
    center_lon = grid_lon((si + ei) / 2)

    nobsl = 0

    do i = -north_n, south_n

        do j = -obs_lc(i)%west_n, obs_lc(i)%east_n

    ! do i = 0, 0
    !     do j = 0, 0

            do k = 1, obs_type_num

#define current obs_lc(i)%obs_data(k, j)

                do l = 1, current%obs_num

                    if (abs(current%lat(l)) <= obs_polar) then

                        obs_dist = latlon_dist(current%lat(l), current%lon(l),&
                        center_lat, center_lon)

                        if (gc_filter == .true.) then
                            gc = schr2(obs_dist * 1d0 / local_dist)
                        else
                            gc = schr(obs_dist * 1d0 / local_dist)
                        end if

                        if (gc > 0) then
                            ! if (pid == 0 .and. nobsl < 100) then
                            !     print *, obs_dist, i, j, k, l, obs_dist, gc
                            ! end if
                            do m = 1, current%obs_nz
                                do n = 1, current%obs_nvar
                                    nobsl = nobsl + 1

                                    yb(:, nobsl) = current%pro_yb(:, n, m, l)
                                    dis(nobsl) = current%dis_yo(n, m, l)
                                    rdiaginv(nobsl) = current%obs_rinv(n, m, l)
                                    rloc(nobsl) = gc

                                    ! if (pid == 0 .and. current%pro_yb(1, n, m, l) == 0) then
                                    !     print *, "obs wrong", "lc", i, "ea", j, "num", l, current%dis_yo(n, m, l)
                                    !     print *, "obs wrong", "nvar", m, "nz", n
                                    !     print *, current%pro_yb(:, :, :, l)
                                    !     call log_error("0 obs!")
                                    ! end if

                                    if (nobsl >= MAX_OBS) then
#if (DEBUG == 1)
                                        !call log_notice("number of observations out of limit!")
#endif
                                        goto 10
                                    end if

                                end do
                            end do
                        end if

                    end if
                end do

#undef current

            end do

        end do

    end do

10  return!generate yb end

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