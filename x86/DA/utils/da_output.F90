#include "da_config.inc"
module da_output_mod

    use coupler_config, only:case_name, atm_mpas_sceneid, num_lon, num_lat, num_lev
    use da_namelist_mod, only:lonstep=>obs_stat_interval_lon, latstep=>obs_stat_interval_lat, &
                         zstep=>obs_stat_interval_vertical, mpas_num_lev_start, mpas_num_lev, ps_pos
    use redis_mod
    use flogger
    use string
    use fiona
    use da_vert_coord_mod

    implicit none

contains

    subroutine da_output_redis(pid, MPI_DA_GROUP,rc, latids, latide, lonids, lonide, &
               da_2d_num, da_3d_num, da_var_name_list, date_char, axb, axa, spread_b, spread_a)
        integer, intent(in) :: pid
        integer, intent(in) :: MPI_DA_GROUP
        type(c_ptr), intent(in) :: rc !redis connect
        integer, intent(in) :: latids, latide, lonids, lonide
!进程起止经纬度id号
        integer, intent(in) :: da_2d_num, da_3d_num
        character(len = 20), dimension(da_2d_num+da_3d_num), intent(in) :: da_var_name_list ! var sort of axb and axa
        character(12), intent(in) :: date_char ! time %y%m%d%hh%mm
        axb_type, intent(in) :: axb(da_2d_num+da_3d_num*num_lev, latide-latids+1, lonide-lonids+1) ! mean bg field
        axb_type, intent(in) :: axa(da_2d_num+da_3d_num*num_lev, latide-latids+1, lonide-lonids+1) ! mean an field
        real(8), intent(in) :: spread_b(da_2d_num+da_3d_num*num_lev, latide-latids+1, lonide-lonids+1)
        real(8), intent(in) :: spread_a(da_2d_num+da_3d_num*num_lev, latide-latids+1, lonide-lonids+1)

        axb_type :: axb_tmp(da_2d_num+da_3d_num*num_lev)
        axb_type :: axa_tmp(da_2d_num+da_3d_num*num_lev)

        integer :: i,j,k,m,n,var,tl,tk,ik
        real(8) :: t3
        character(5) :: num_lon_str, num_lat_str
        axb_type, allocatable, dimension(:,:,:) :: axa_4d, axb_4d
        real(4), allocatable, dimension(:,:,:) :: ph_ref
        axb_type :: ph_lev_axb(num_lev+1), ph_lev_axa(num_lev+1), ph_axb(num_lev), ph_axa(num_lev)
        integer :: ph_pos_axb((mpas_num_lev - mpas_num_lev_start)/zstep+1), &
                   ph_pos_axa((mpas_num_lev - mpas_num_lev_start)/zstep+1)
        character(256) :: hash_key, hash_key_amb, hash_key_axb, hash_key_spread_b, hash_key_spread_a, var_key
        integer :: latids_out, latide_out, lonids_out, lonide_out
        integer :: latids_out_loc, latide_out_loc, lonids_out_loc, lonide_out_loc

!!!!!!!debug
! if ( pid==339 ) print *, "latids, latide, lonids, lonide", latids, latide, lonids, lonide      
!!!!!!!!!!!!
        do i = 1, da_2d_num+da_3d_num*num_lev
            do j = 1, latide-latids+1
                do k = 1, lonide-lonids+1
                     if (isnan(axb(i,j,k))) print *,"pid=", pid, 'bg nan'
                     if (isnan(axa(i,j,k))) print *,"pid=", pid, 'an nan'
                enddo
            enddo
        enddo
        lonids_out = int(ceiling((float(lonids)-1.0)/float(lonstep))) + 1
        lonids_out_loc = (lonids_out-1)*lonstep+1-lonids+1
        if (lonide == num_lon .and. (mod((float(num_lon)-1.0),float(lonstep)) .ne. 0)) then
            lonide_out = int(floor((float(lonide)-1.0)/float(lonstep))) + 2
            lonide_out_loc = num_lon-lonids+1
        else
            lonide_out = int(floor((float(lonide)-1.0)/float(lonstep))) + 1
            lonide_out_loc = (lonide_out-1)*lonstep+1-lonids+1
        endif

        latids_out = int(ceiling((float(latids)-1.0)/float(latstep))) + 1
        latids_out_loc = (latids_out-1)*latstep+1-latids+1
        if (latide ==num_lat .and. (mod((float(num_lat)-1.0),float(latstep)) .ne. 0)) then
            latide_out = int(floor((float(latide)-1.0)/float(latstep))) + 2
            latide_out_loc = num_lat-latids+1
        else
            latide_out = int(floor((float(latide)-1.0)/float(latstep))) + 1
            latide_out_loc = (latide_out-1)*latstep+1-latids+1
        endif

        write(num_lon_str, "(i5)") num_lon
        write(num_lat_str, "(i5)") num_lat
        
        hash_key_amb = "ambfield:" //trim(adjustl(case_name)) //":" // &
                       trim(adjustl(atm_mpas_sceneid)) // ':' //&
                       trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str)) // &
                       ":" //trim(adjustl(date_char))
        hash_key_axb = "axbfield:" //trim(adjustl(case_name)) //":" // &
                       trim(adjustl(atm_mpas_sceneid)) // ':' //&
                       trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str)) // &
                       ":" //trim(adjustl(date_char))
        hash_key_spread_b = "spreadb:" //trim(adjustl(case_name)) //":" // &
                       trim(adjustl(atm_mpas_sceneid)) // ':' //&
                       trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str)) // &
                       ":" //trim(adjustl(date_char))
        hash_key_spread_a = "spreada:" //trim(adjustl(case_name)) //":" // &
                       trim(adjustl(atm_mpas_sceneid)) // ':' //&
                       trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str)) // &
                       ":" //trim(adjustl(date_char))

        var_key = trim(adjustl(da_var_name_list(1)))
        do i=2, da_2d_num+da_3d_num
            var_key = trim(adjustl(var_key)) // ':' // trim(adjustl(da_var_name_list(i)))
        end do

        call log_print_root("varlist: " // var_key)
        call log_print_root("da_2d_num: " // to_str(da_2d_num))
        call log_print_root("da_3d_num: " // to_str(da_3d_num))
        call log_print_root("date_char: " // date_char)
        call log_print_root("lonstep: " // to_str(lonstep))
        call log_print_root("latstep: " // to_str(latstep))
        call log_print_root("ztep: " //    to_str(zstep))

        allocate(axa_4d(da_2d_num+da_3d_num*mpas_num_lev, latide_out-latids_out+1, lonide_out-lonids_out+1))
        allocate(axb_4d(da_2d_num+da_3d_num*mpas_num_lev, latide_out-latids_out+1, lonide_out-lonids_out+1))
        allocate(ph_ref((mpas_num_lev - mpas_num_lev_start)/zstep+1, latide_out-latids_out+1, lonide_out-lonids_out+1))

        axa_4d = 0
        axb_4d = 0
        if (ps_pos == 0) then            
            axa_4d = axa
            axb_4d = axb
        else
            if (mod(num_lon-1, lonstep)==0 .and. mod(num_lat-1, latstep)==0) then
                axa_4d(1:da_2d_num, :, :) = axa(1:da_2d_num, &
                latids_out_loc:latide_out_loc:latstep, &
                lonids_out_loc:lonide_out_loc:lonstep)
                axb_4d(1:da_2d_num, :, :) = axb(1:da_2d_num, &
                latids_out_loc:latide_out_loc:latstep, &
                lonids_out_loc:lonide_out_loc:lonstep)
            else
                if ((latide /= num_lat) .and. (lonide /= num_lon)) then
                    axa_4d(1:da_2d_num, :, :) = axa(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    axb_4d(1:da_2d_num, :, :) = axb(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                elseif ((latide /= num_lat) .and. (lonide == num_lon)) then
                    axa_4d(1:da_2d_num, :, 1:lonide_out-lonids_out) = axa(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    axa_4d(1:da_2d_num, :, lonide_out-lonids_out+1) = axa(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, lonide_out_loc)

                    axb_4d(1:da_2d_num, :, 1:lonide_out-lonids_out) = axb(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    axb_4d(1:da_2d_num, :, lonide_out-lonids_out+1) = axb(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, lonide_out_loc)
                elseif ((latide == num_lat) .and. (lonide /= num_lon)) then
                    axa_4d(1:da_2d_num, 1:latide_out-latids_out, :) = axa(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    axa_4d(1:da_2d_num, latide_out-latids_out+1, :) = axa(1:da_2d_num, &
                    latide_out_loc, lonids_out_loc:lonide_out_loc:lonstep)

                    axb_4d(1:da_2d_num, 1:latide_out-latids_out, :) = axb(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    axb_4d(1:da_2d_num, latide_out-latids_out+1, :) = axb(1:da_2d_num, &
                    latide_out_loc, lonids_out_loc:lonide_out_loc:lonstep)
                elseif ((latide == num_lat) .and. (lonide == num_lon)) then
                    axa_4d(1:da_2d_num, 1:latide_out-latids_out, 1:lonide_out-lonids_out)&
                    = axa(1:da_2d_num, latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    axa_4d(1:da_2d_num, latide_out-latids_out+1, 1:lonide_out-lonids_out)&
                    = axa(1:da_2d_num, latide_out_loc, lonids_out_loc:lonide_out_loc:lonstep)
                    axa_4d(1:da_2d_num, 1:latide_out-latids_out, lonide_out-lonids_out+1)&
                    = axa(1:da_2d_num, latids_out_loc:latide_out_loc:latstep, lonide_out_loc)
                    axa_4d(1:da_2d_num, latide_out-latids_out+1, lonide_out-lonids_out+1)&
                    = axa(1:da_2d_num, latide_out_loc, lonide_out_loc)

                    axb_4d(1:da_2d_num, 1:latide_out-latids_out, 1:lonide_out-lonids_out)&
                    = axb(1:da_2d_num, latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    axb_4d(1:da_2d_num, latide_out-latids_out+1, 1:lonide_out-lonids_out)&
                    = axb(1:da_2d_num, latide_out_loc, lonids_out_loc:lonide_out_loc:lonstep)
                    axb_4d(1:da_2d_num, 1:latide_out-latids_out, lonide_out-lonids_out+1)&
                    = axb(1:da_2d_num, latids_out_loc:latide_out_loc:latstep, lonide_out_loc)
                    axb_4d(1:da_2d_num, latide_out-latids_out+1, lonide_out-lonids_out+1)&
                    = axb(1:da_2d_num, latide_out_loc, lonide_out_loc)
                endif
            endif
            
!            print*, "ps_pos",ps_pos, "pid",pid
            hash_key = "realfield:" // trim(adjustl(atm_mpas_sceneid)) //":" //&
               trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str)) &
               // ":" //trim(adjustl(date_char))
            if ((lonstep .ne. 1) .and. (latstep .ne. 1)) then
                call RedisHmgetf3d(rc, hash_key, 'p', lonids_out, lonide_out,&
                    1, latids_out, latide_out, 1, &
                    mpas_num_lev_start, mpas_num_lev, zstep, 1, ph_ref)
            elseif((lonstep == 1) .and. (latstep == 1)) then
                call RedisHmgetf3d(rc, hash_key, 'p', lonids, lonide, &
                     lonstep, latids, latide, latstep, &
                     mpas_num_lev_start, mpas_num_lev, zstep, 1, ph_ref)
            else
                call log_error("obs_stat_interval_lon and &
                     obs_stat_interval_lat are wrong" ,__FILE__, __LINE__)
            endif

            do i = 1, latide_out-latids_out+1
                do j = 1, lonide_out-lonids_out+1
                    ! calculate gmcore ph
                    if ((latids_out_loc+(i-1)*latstep+latids-1)<num_lat .and. &
                        (lonids_out_loc+(j-1)*lonstep+lonids-1)<num_lon) then
                        axb_tmp = axb(:, latids_out_loc+(i-1)*latstep, lonids_out_loc+(j-1)*lonstep)
                        axa_tmp = axa(:, latids_out_loc+(i-1)*latstep, lonids_out_loc+(j-1)*lonstep)
                    elseif ((latids_out_loc+(i-1)*latstep+latids-1)>=num_lat .and. &
                        (lonids_out_loc+(j-1)*lonstep+lonids-1)<num_lon) then
                        axb_tmp = axb(:, latide_out_loc, lonids_out_loc+(j-1)*lonstep)
                        axa_tmp = axa(:, latide_out_loc, lonids_out_loc+(j-1)*lonstep)
                    elseif ((latids_out_loc+(i-1)*latstep+latids-1)<num_lat .and. &
                        (lonids_out_loc+(j-1)*lonstep+lonids-1)>=num_lon) then
                        axb_tmp = axb(:, latids_out_loc+(i-1)*latstep, lonide_out_loc)
                        axa_tmp = axa(:, latids_out_loc+(i-1)*latstep, lonide_out_loc)
                    elseif ((latids_out_loc+(i-1)*latstep+latids-1)>=num_lat .and. &
                        (lonids_out_loc+(j-1)*lonstep+lonids-1)>=num_lon) then
                        axb_tmp = axb(:, latide_out_loc, lonide_out_loc)
                        axa_tmp = axa(:, latide_out_loc, lonide_out_loc)
                    endif

                    do k = 1, num_lev+1
                        ph_lev_axb(k) = da_vert_coord_calc_ph_lev(k, axb_tmp(ps_pos))
                        ph_lev_axa(k) = da_vert_coord_calc_ph_lev(k, axa_tmp(ps_pos))
                    enddo

                    do k = 1, num_lev
                        ph_axb(k) = 0.5d0 * (ph_lev_axb(k) + ph_lev_axb(k+1))
                        ph_axa(k) = 0.5d0 * (ph_lev_axa(k) + ph_lev_axa(k+1))
                    enddo
!!!!!!!!debug
!    if (pid == 339 .and. i == 4 .and. j == 96) then
!        print *, "ph_axb", ph_axb, "ph_axa", ph_axa
!        print *, "ph_ref", ph_ref(:, i, j)
!    endif
!!!!!!!!!!!!
                    ! interp
!                    m = 1
!                    n = 1
                    do k = 1, (mpas_num_lev - mpas_num_lev_start)/zstep+1
                        do ik = 1, num_lev
                            if (ph_ref(k, i, j) - ph_axb(ik) > 0.00001d0) m = ik+1
                        end do
!                        do while (m <= num_lev .and. ph_ref(k, i, j) - ph_axb(m) > 0.00001d0)
!                            m = m + 1
!                            if (m > num_lev) exit
!                        end do
                        ph_pos_axb(k) = m - 1
                        do ik = 1, num_lev
                            if (ph_ref(k, i, j) - ph_axa(ik) > 0.00001d0) n = ik+1
                        end do
!                        do while (n <= num_lev .and. ph_ref(k, i, j) - ph_axa(n) > 0.00001d0)
!                            n = n + 1
!                            if (n > num_lev) exit
!                        end do
                        ph_pos_axa(k) = n - 1
                    enddo
!!!!!!!!debug
!    if (pid == 339 .and. i == 4 .and. j == 96) then
!        print *, "ph_pos_axb", ph_pos_axb, "ph_pos_axa", ph_pos_axa
!    endif
!!!!!!!!!!!!
                    do k = 1, (mpas_num_lev - mpas_num_lev_start)/zstep+1
                        if (ph_pos_axb(k) >= 1 .and. ph_pos_axb(k) < num_lev) then
                            do var = 1, da_3d_num
                                tl = da_2d_num + num_lev*(var-1)+ph_pos_axb(k)
                                tk = da_2d_num + mpas_num_lev*(var-1)+mpas_num_lev_start+(k-1)*zstep
                                t3 = (log(ph_axb(ph_pos_axb(k)+1)) - log(ph_ref(k, i, j))) / &
                                     (log(ph_axb(ph_pos_axb(k)+1)) - log(ph_axb(ph_pos_axb(k))))
                                axb_4d(tk, i, j) = axb_tmp(tl) * t3 + axb_tmp(tl+1) * (1 - t3)
                            enddo
                        elseif (ph_pos_axb(k) == 0) then
                            do var = 1, da_3d_num
                                tl = da_2d_num + num_lev*(var-1)+1
                                tk = da_2d_num + mpas_num_lev*(var-1)+mpas_num_lev_start+(k-1)*zstep
                                axb_4d(tk, i, j) = axb_tmp(tl)
                            enddo
                        elseif (ph_pos_axb(k) >= num_lev) then
                            do var = 1, da_3d_num
                                tl = da_2d_num + num_lev*var
                                tk = da_2d_num + mpas_num_lev*(var-1)+mpas_num_lev_start+(k-1)*zstep
                                axb_4d(tk, i, j) = axb_tmp(tl)
                            enddo
                        endif

                        if (ph_pos_axa(k) >= 1 .and. ph_pos_axa(k) < num_lev) then
                            do var = 1, da_3d_num
                                tl = da_2d_num + num_lev*(var-1)+ph_pos_axa(k)
                                tk = da_2d_num + mpas_num_lev*(var-1)+mpas_num_lev_start+(k-1)*zstep
                                t3 = (log(ph_axa(ph_pos_axa(k)+1)) - log(ph_ref(k, i, j))) / &
                                     (log(ph_axa(ph_pos_axa(k)+1)) - log(ph_axa(ph_pos_axa(k))))
                                axa_4d(tk, i, j) = axa_tmp(tl) * t3 + axa_tmp(tl+1) * (1 - t3)
                            enddo
                        elseif (ph_pos_axa(k) == 0) then
                            do var = 1, da_3d_num
                                tl = da_2d_num + num_lev*(var-1)+1
                                tk = da_2d_num + mpas_num_lev*(var-1)+mpas_num_lev_start+(k-1)*zstep
                                axa_4d(tk, i, j) = axa_tmp(tl)
                            enddo
                        elseif (ph_pos_axa(k) >= num_lev) then
                            do var = 1, da_3d_num
                                tl = da_2d_num + num_lev*var
                                tk = da_2d_num + mpas_num_lev*(var-1)+mpas_num_lev_start+(k-1)*zstep
                                axa_4d(tk, i, j) = axa_tmp(tl)
                            enddo
                        endif
                    enddo
!!!!!!!!!!debug
!    if (pid == 339 .and. i == 4 .and. j == 96) then
!        print *, "tl,tk,t3",tl,tk,t3
!        print *, "axb_4d", axb_4d(:,i,j), "axa_4d", axa_4d(:,i,j)
!        print *, "axb", axb(:,i,j), "axa", axa(:,i,j)
!    endif
!!!!!!!!!!!!!
                enddo
            enddo
        endif

        if ((lonstep .ne. 1) .and. (latstep .ne. 1)) then
            call redis_da_outputd(rc, hash_key_amb, var_key, lonids_out, lonide_out, &
                1, latids_out, latide_out, 1, mpas_num_lev_start, mpas_num_lev, zstep, &
                da_2d_num, da_3d_num, axa_4d-axb_4d)
            call redis_da_outputd(rc, hash_key_axb, var_key, lonids_out, lonide_out, &
                1, latids_out, latide_out, 1, mpas_num_lev_start, mpas_num_lev, zstep, &
                da_2d_num, da_3d_num, axb_4d)
            call redis_da_outputd(rc, hash_key_spread_b, var_key, lonids_out, lonide_out, &
                1, latids_out, latide_out, 1, mpas_num_lev_start, mpas_num_lev, zstep, &
                da_2d_num, da_3d_num, spread_b)
            call redis_da_outputd(rc, hash_key_spread_a, var_key, lonids_out, lonide_out, &
                1, latids_out, latide_out, 1, mpas_num_lev_start, mpas_num_lev, zstep, &
                da_2d_num, da_3d_num, spread_a)

        elseif ((lonstep == 1) .and. (latstep == 1)) then
            call redis_da_outputd(rc, hash_key_amb, var_key, lonids, lonide, &
                lonstep, latids, latide, latstep, mpas_num_lev_start, &
                mpas_num_lev, zstep, da_2d_num, da_3d_num, axa_4d-axb_4d)
            call redis_da_outputd(rc, hash_key_axb, var_key, lonids, lonide, &
                lonstep, latids, latide, latstep, mpas_num_lev_start, &
                mpas_num_lev, zstep, da_2d_num, da_3d_num, axb_4d)
            call redis_da_outputd(rc, hash_key_spread_b, var_key, lonids, lonide, &
                lonstep, latids, latide, latstep, mpas_num_lev_start, &
                mpas_num_lev, zstep, da_2d_num, da_3d_num, spread_b)
            call redis_da_outputd(rc, hash_key_spread_a, var_key, lonids, lonide, &
                lonstep, latids, latide, latstep, mpas_num_lev_start, &
                mpas_num_lev, zstep, da_2d_num, da_3d_num, spread_a)
        else
            call log_error("obs_stat_interval_lon and &
                 obs_stat_interval_lat are wrong" ,__FILE__, __LINE__)
        endif

        deallocate(axa_4d)
        deallocate(axb_4d)
        deallocate(ph_ref)
  end subroutine da_output_redis

    subroutine da_output_nc(pid, MPI_DA_GROUP, rc, latids, latide, lonids, lonide, &
            da_2d_num, da_3d_num, da_var_name_list, date_char, axb, axa, spread_b, spread_a)

        type(c_ptr), intent(in) :: rc !redis connect
        integer, intent(in) :: pid
        integer, intent(in) :: latids, latide, lonids, lonide !进程起止经纬度id号
        integer, intent(in) :: da_2d_num, da_3d_num
        character(len = 20), dimension(da_2d_num+da_3d_num), intent(in) :: da_var_name_list ! var sort of axb and axa
        character(12), intent(in) :: date_char ! time %y%m%d%hh%mm
        axb_type, intent(in) :: axb(da_2d_num+da_3d_num*num_lev, latide-latids+1, lonide-lonids+1) ! mean bg field
        axb_type, intent(in) :: axa(da_2d_num+da_3d_num*num_lev, latide-latids+1, lonide-lonids+1) ! mean an field
        real(8), intent(in) :: spread_b(da_2d_num+da_3d_num*num_lev, latide-latids+1, lonide-lonids+1)
        real(8), intent(in) :: spread_a(da_2d_num+da_3d_num*num_lev, latide-latids+1, lonide-lonids+1)

        axb_type :: axb_tmp(da_2d_num+da_3d_num*num_lev)
        axb_type :: axa_tmp(da_2d_num+da_3d_num*num_lev)
        real(8)  :: spread_b_tmp(da_2d_num+da_3d_num*num_lev)
        real(8)  :: spread_a_tmp(da_2d_num+da_3d_num*num_lev)
        axb_type, allocatable, dimension(:,:,:) :: tmp_3d
        axb_type, allocatable, dimension(:,:) :: tmp_2d

        integer :: i, j, k, ii, jj, kk, m, n, var, tl, tk, ik
        real(8) :: t3
        axb_type, allocatable, dimension(:,:,:) :: axa_4d, axb_4d
        real(8), allocatable, dimension(:,:,:) :: spread_b_4d, spread_a_4d
        real(4), allocatable, dimension(:,:,:) :: ph_ref
!        axb_type :: axa_4d(da_2d_num+da_3d_num*mpas_num_lev, latide-latids+1, lonide-lonids+1)
!        axb_type :: axb_4d(da_2d_num+da_3d_num*mpas_num_lev, latide-latids+1, lonide-lonids+1)
!        real(4) :: ph_ref((mpas_num_lev - mpas_num_lev_start)/zstep+1, latide-latids+1, lonide-lonids+1)
        axb_type :: ph_lev_axb(num_lev+1), ph_lev_axa(num_lev+1), ph_axb(num_lev), ph_axa(num_lev)
        integer :: ph_pos_axb((mpas_num_lev - mpas_num_lev_start)/zstep+1), &
                   ph_pos_axa((mpas_num_lev - mpas_num_lev_start)/zstep+1)
        character(5) :: num_lon_str, num_lat_str
        character(256) :: hash_key, hash_key_amb, hash_key_axb, var_key

        ! netcdf var
        integer, intent(in) :: MPI_DA_GROUP
        character*20 fid0,fid1
        integer start1(3), count1(3)
        integer start2(2), count2(2)
        integer start3(1), count3(1)
        integer, allocatable, dimension(:) :: levdeg!, londeg, latdeg,
        integer :: num_lon_out, num_lat_out, num_lev_out
        integer :: latids_out, latide_out, lonids_out, lonide_out
        integer :: latids_out_loc, latide_out_loc, lonids_out_loc, lonide_out_loc

        do i = 1, da_2d_num+da_3d_num*num_lev
            do j = 1, latide-latids+1
                do k = 1, lonide-lonids+1
                     if (isnan(axb(i,j,k))) print *,"pid=", pid, 'bg nan'
                     if (isnan(axa(i,j,k))) print *,"pid=", pid, 'an nan'
                enddo
            enddo
        enddo
        lonids_out = int(ceiling((float(lonids)-1.0)/float(lonstep))) + 1
        lonids_out_loc = (lonids_out-1)*lonstep+1-lonids+1
        if (lonide == num_lon .and. (mod((float(num_lon)-1.0),float(lonstep)) .ne. 0)) then
            lonide_out = int(floor((float(lonide)-1.0)/float(lonstep))) + 2
            lonide_out_loc = num_lon-lonids+1
        else
            lonide_out = int(floor((float(lonide)-1.0)/float(lonstep))) + 1
            lonide_out_loc = (lonide_out-1)*lonstep+1-lonids+1
        endif

        latids_out = int(ceiling((float(latids)-1.0)/float(latstep))) + 1
        latids_out_loc = (latids_out-1)*latstep+1-latids+1
        if (latide ==num_lat .and. (mod((float(num_lat)-1.0),float(latstep)) .ne. 0)) then
            latide_out = int(floor((float(latide)-1.0)/float(latstep))) + 2
            latide_out_loc = num_lat-latids+1
        else
            latide_out = int(floor((float(latide)-1.0)/float(latstep))) + 1
            latide_out_loc = (latide_out-1)*latstep+1-latids+1
        endif

        if (mod(float(num_lon)-1.0, float(lonstep)) == 0) then
            num_lon_out=floor((float(num_lon)-1.0)/float(lonstep))+1
        else
            num_lon_out=floor((float(num_lon)-1.0)/float(lonstep))+2
        endif

        if (mod(float(num_lat)-1.0, float(latstep)) == 0) then
            num_lat_out=floor((float(num_lat)-1.0)/float(latstep))+1
        else
            num_lat_out=floor((float(num_lat)-1.0)/float(latstep))+2
        endif
!!!!debug
!if (pid==101) print *, pid,"lonids, lonide, latids, latide", lonids, lonide, latids, latide,"lonids_out, lonide_out, latids_out, latide_out, lonids_out_loc, lonide_out_loc, latids_out_loc, latide_out_loc, ", &
!lonids_out, lonide_out, latids_out, latide_out, lonids_out_loc, lonide_out_loc, latids_out_loc, latide_out_loc
!!!!!!!!

        num_lev_out=(mpas_num_lev-mpas_num_lev_start)/zstep+1

        write(num_lon_str, "(i5)") num_lon
        write(num_lat_str, "(i5)") num_lat
        fid0 ='amb'
        fid1 ='axb'

        hash_key_amb = "ambfield:" //trim(adjustl(case_name)) //":" // &
                       trim(adjustl(atm_mpas_sceneid)) // ':' //&
                       trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str)) // &
                       ":" //trim(adjustl(date_char)) // ".nc"
        hash_key_axb = "axbfield:" //trim(adjustl(case_name)) //":" // &
                       trim(adjustl(atm_mpas_sceneid)) // ':' //&
                       trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str)) // &
                       ":" //trim(adjustl(date_char)) // ".nc"

        call log_print_root("da_2d_num: " // to_str(da_2d_num))
        call log_print_root("da_3d_num: " // to_str(da_3d_num))
        call log_print_root("date_char: " // date_char)
        call log_print_root("lonstep: " // to_str(lonstep))
        call log_print_root("latstep: " // to_str(latstep))
        call log_print_root("ztep: " //    to_str(zstep))

        call fiona_init()
        call fiona_create_dataset(fid0, file_path=trim(hash_key_amb), mpi_comm=MPI_DA_GROUP)
        call fiona_add_att(fid0, 'num_lon', num_lon)
        call fiona_add_att(fid0, 'num_lat', num_lat)
        call fiona_add_att(fid0, 'num_lev', num_lev)
        call fiona_add_att(fid0, 'num_lev_stride', zstep)
        call fiona_add_att(fid0, 'output_num_lev', num_lev_out)
        call fiona_add_att(fid0, 'da_2d_num', da_2d_num)
        call fiona_add_att(fid0, 'da_3d_num', da_3d_num)
        call fiona_add_dim(fid0, 'lon' , size=num_lon_out, add_var=.true., decomp=.true.)
        call fiona_add_dim(fid0, 'lat' , size=num_lat_out, add_var=.true., decomp=.true.)
        call fiona_add_dim(fid0, 'lev' , size=num_lev_out, add_var=.true., decomp=.false.)

        call fiona_create_dataset(fid1, file_path=trim(hash_key_axb), mpi_comm=MPI_DA_GROUP)
        call fiona_add_att(fid1, 'num_lon', num_lon)
        call fiona_add_att(fid1, 'num_lat', num_lat)
        call fiona_add_att(fid1, 'num_lev', num_lev)
        call fiona_add_att(fid1, 'num_lev_stride', zstep)
        call fiona_add_att(fid1, 'output_num_lev', num_lev_out)
        call fiona_add_att(fid1, 'da_2d_num', da_2d_num)
        call fiona_add_att(fid1, 'da_3d_num', da_3d_num)
        call fiona_add_dim(fid1, 'lon' , size=num_lon_out, add_var=.false., decomp=.true.)
        call fiona_add_dim(fid1, 'lat' , size=num_lat_out, add_var=.false., decomp=.true.)
        call fiona_add_dim(fid1, 'lev' , size=num_lev_out, add_var=.true., decomp=.false.)

        do i=1, da_2d_num
            call fiona_add_var(fid0, trim(adjustl(da_var_name_list(i))), &
                 long_name=trim(adjustl(da_var_name_list(i))), &
                 units='', dim_names=['lon', 'lat'], data_type='r8')
            call fiona_add_var(fid1, trim(adjustl(da_var_name_list(i))), &
                 long_name=trim(adjustl(da_var_name_list(i))), &
                 units='', dim_names=['lon', 'lat'], data_type='r8')
            call fiona_add_var(fid1, 'spread_b_'//trim(adjustl(da_var_name_list(i))), &
                 long_name='spread of '//trim(adjustl(da_var_name_list(i)))//' in background field', &
                 units='', dim_names=['lon', 'lat'], data_type='r8')
            call fiona_add_var(fid1, 'spread_a_'//trim(adjustl(da_var_name_list(i))), &
                 long_name='spread of '//trim(adjustl(da_var_name_list(i)))//' in analysis field', &
                 units='', dim_names=['lon', 'lat'], data_type='r8')
        end do

        do i=1, da_3d_num
            call fiona_add_var(fid0, trim(adjustl(da_var_name_list(da_2d_num+i))), &
                 long_name=trim(adjustl(da_var_name_list(da_2d_num+i))), &
                 units='', dim_names=['lon', 'lat', 'lev'], data_type='r8')
            call fiona_add_var(fid1, trim(adjustl(da_var_name_list(da_2d_num+i))), &
                 long_name=trim(adjustl(da_var_name_list(da_2d_num+i))), &
                 units='', dim_names=['lon', 'lat', 'lev'], data_type='r8')
            call fiona_add_var(fid1, 'spread_b_'//trim(adjustl(da_var_name_list(da_2d_num+i))), &
                 long_name='spread of '//trim(adjustl(da_var_name_list(da_2d_num+i)))//' in background field', &
                 units='', dim_names=['lon', 'lat', 'lev'], data_type='r8')
            call fiona_add_var(fid1, 'spread_a_'//trim(adjustl(da_var_name_list(da_2d_num+i))), &
                 long_name='spread of '//trim(adjustl(da_var_name_list(da_2d_num+i)))//' in analysis field', &
                 units='', dim_names=['lon', 'lat', 'lev'], data_type='r8')
            
        end do

        call log_print_root("da_output nc init")
        call fiona_start_output(fid0)
        call fiona_start_output(fid1)

        allocate(tmp_3d(lonide_out-lonids_out+1, latide_out-latids_out+1, num_lev_out ))
        allocate(tmp_2d(lonide_out-lonids_out+1, latide_out-latids_out+1))
        allocate(axa_4d(da_2d_num+da_3d_num*mpas_num_lev, latide_out-latids_out+1, lonide_out-lonids_out+1))
        allocate(axb_4d(da_2d_num+da_3d_num*mpas_num_lev, latide_out-latids_out+1, lonide_out-lonids_out+1))
        allocate(spread_a_4d(da_2d_num+da_3d_num*mpas_num_lev, latide_out-latids_out+1, lonide_out-lonids_out+1))
        allocate(spread_b_4d(da_2d_num+da_3d_num*mpas_num_lev, latide_out-latids_out+1, lonide_out-lonids_out+1))
        allocate(ph_ref((mpas_num_lev - mpas_num_lev_start)/zstep+1, latide_out-latids_out+1, lonide_out-lonids_out+1))
!        allocate(londeg(lonids:lonide))
!        allocate(latdeg(latids:latide))
        allocate(levdeg(1:num_lev_out))
!        do i = lonids, lonide
!            londeg(i) = i
!        end do
!
!        do i = latids, latide
!            latdeg(i) = i
!        end do

        k = 1
        do i = mpas_num_lev_start, mpas_num_lev, zstep
            levdeg(k) = i
            k = k + 1
        end do

        call log_print_root("output_num_lev:" // to_str(num_lev_out))
!        start=[lonids, latids,1]
!        count=[lonide-lonids+1, latide-latids+1,num_lev_out]
!        call fiona_output(fid0, 'lonid', londeg(lonids:lonide), start=start, count=count)
!        call fiona_output(fid1, 'lonid', londeg(lonids:lonide), start=start, count=count)
!
!        start3=[latids]
!        count3=[latide-latids+1]
!        call fiona_output(fid0, 'lat', latdeg(latids:latide), start=start3, count=count3)
!        call fiona_output(fid1, 'lat', latdeg(latids:latide), start=start3, count=count3)

        call fiona_output(fid0, 'lev', levdeg(1:num_lev_out))
        call fiona_output(fid1, 'lev', levdeg(1:num_lev_out))

        ! interp
        axa_4d = 0
        axb_4d = 0
        spread_a_4d = 0
        spread_b_4d = 0
        if (ps_pos == 0) then            
            axa_4d = axa
            axb_4d = axb
            spread_a_4d = spread_a
            spread_b_4d = spread_b
        else
            if (mod(num_lon-1, lonstep)==0 .and. mod(num_lat-1, latstep)==0) then
                axa_4d(1:da_2d_num, :, :) = axa(1:da_2d_num, &
                latids_out_loc:latide_out_loc:latstep, &
                lonids_out_loc:lonide_out_loc:lonstep)
                axb_4d(1:da_2d_num, :, :) = axb(1:da_2d_num, &
                latids_out_loc:latide_out_loc:latstep, &
                lonids_out_loc:lonide_out_loc:lonstep)
                spread_a_4d(1:da_2d_num, :, :) = spread_a(1:da_2d_num, &
                latids_out_loc:latide_out_loc:latstep, &
                lonids_out_loc:lonide_out_loc:lonstep)
                spread_b_4d(1:da_2d_num, :, :) = spread_b(1:da_2d_num, &
                latids_out_loc:latide_out_loc:latstep, &
                lonids_out_loc:lonide_out_loc:lonstep)
            else
                if ((latide /= num_lat) .and. (lonide /= num_lon)) then
                    axa_4d(1:da_2d_num, :, :) = axa(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    axb_4d(1:da_2d_num, :, :) = axb(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    spread_a_4d(1:da_2d_num, :, :) = spread_a(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    spread_b_4d(1:da_2d_num, :, :) = spread_b(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                elseif ((latide /= num_lat) .and. (lonide == num_lon)) then
                    axa_4d(1:da_2d_num, :, 1:lonide_out-lonids_out) = axa(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    axa_4d(1:da_2d_num, :, lonide_out-lonids_out+1) = axa(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, lonide_out_loc)

                    axb_4d(1:da_2d_num, :, 1:lonide_out-lonids_out) = axb(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    axb_4d(1:da_2d_num, :, lonide_out-lonids_out+1) = axb(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, lonide_out_loc)
                    
                    spread_a_4d(1:da_2d_num, :, 1:lonide_out-lonids_out) = spread_a(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    spread_a_4d(1:da_2d_num, :, lonide_out-lonids_out+1) = spread_a(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, lonide_out_loc)

                    spread_b_4d(1:da_2d_num, :, 1:lonide_out-lonids_out) = spread_b(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    spread_b_4d(1:da_2d_num, :, lonide_out-lonids_out+1) = spread_b(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, lonide_out_loc)
                elseif ((latide == num_lat) .and. (lonide /= num_lon)) then
                    axa_4d(1:da_2d_num, 1:latide_out-latids_out, :) = axa(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    axa_4d(1:da_2d_num, latide_out-latids_out+1, :) = axa(1:da_2d_num, &
                    latide_out_loc, lonids_out_loc:lonide_out_loc:lonstep)

                    axb_4d(1:da_2d_num, 1:latide_out-latids_out, :) = axb(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    axb_4d(1:da_2d_num, latide_out-latids_out+1, :) = axb(1:da_2d_num, &
                    latide_out_loc, lonids_out_loc:lonide_out_loc:lonstep)

                    spread_a_4d(1:da_2d_num, 1:latide_out-latids_out, :) = spread_a(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    spread_a_4d(1:da_2d_num, latide_out-latids_out+1, :) = spread_a(1:da_2d_num, &
                    latide_out_loc, lonids_out_loc:lonide_out_loc:lonstep)

                    spread_b_4d(1:da_2d_num, 1:latide_out-latids_out, :) = spread_b(1:da_2d_num, &
                    latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    spread_b_4d(1:da_2d_num, latide_out-latids_out+1, :) = spread_b(1:da_2d_num, &
                    latide_out_loc, lonids_out_loc:lonide_out_loc:lonstep)
                elseif ((latide == num_lat) .and. (lonide == num_lon)) then
                    axa_4d(1:da_2d_num, 1:latide_out-latids_out, 1:lonide_out-lonids_out)&
                    = axa(1:da_2d_num, latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    axa_4d(1:da_2d_num, latide_out-latids_out+1, 1:lonide_out-lonids_out)&
                    = axa(1:da_2d_num, latide_out_loc, lonids_out_loc:lonide_out_loc:lonstep)
                    axa_4d(1:da_2d_num, 1:latide_out-latids_out, lonide_out-lonids_out+1)&
                    = axa(1:da_2d_num, latids_out_loc:latide_out_loc:latstep, lonide_out_loc)
                    axa_4d(1:da_2d_num, latide_out-latids_out+1, lonide_out-lonids_out+1)&
                    = axa(1:da_2d_num, latide_out_loc, lonide_out_loc)

                    axb_4d(1:da_2d_num, 1:latide_out-latids_out, 1:lonide_out-lonids_out)&
                    = axb(1:da_2d_num, latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    axb_4d(1:da_2d_num, latide_out-latids_out+1, 1:lonide_out-lonids_out)&
                    = axb(1:da_2d_num, latide_out_loc, lonids_out_loc:lonide_out_loc:lonstep)
                    axb_4d(1:da_2d_num, 1:latide_out-latids_out, lonide_out-lonids_out+1)&
                    = axb(1:da_2d_num, latids_out_loc:latide_out_loc:latstep, lonide_out_loc)
                    axb_4d(1:da_2d_num, latide_out-latids_out+1, lonide_out-lonids_out+1)&
                    = axb(1:da_2d_num, latide_out_loc, lonide_out_loc)

                    spread_a_4d(1:da_2d_num, 1:latide_out-latids_out, 1:lonide_out-lonids_out)&
                    = spread_a(1:da_2d_num, latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    spread_a_4d(1:da_2d_num, latide_out-latids_out+1, 1:lonide_out-lonids_out)&
                    = spread_a(1:da_2d_num, latide_out_loc, lonids_out_loc:lonide_out_loc:lonstep)
                    spread_a_4d(1:da_2d_num, 1:latide_out-latids_out, lonide_out-lonids_out+1)&
                    = spread_a(1:da_2d_num, latids_out_loc:latide_out_loc:latstep, lonide_out_loc)
                    spread_a_4d(1:da_2d_num, latide_out-latids_out+1, lonide_out-lonids_out+1)&
                    = spread_a(1:da_2d_num, latide_out_loc, lonide_out_loc)

                    spread_b_4d(1:da_2d_num, 1:latide_out-latids_out, 1:lonide_out-lonids_out)&
                    = spread_b(1:da_2d_num, latids_out_loc:latide_out_loc:latstep, &
                    lonids_out_loc:lonide_out_loc:lonstep)
                    spread_b_4d(1:da_2d_num, latide_out-latids_out+1,1:lonide_out-lonids_out)&
                    = spread_b(1:da_2d_num, latide_out_loc, lonids_out_loc:lonide_out_loc:lonstep)
                    spread_b_4d(1:da_2d_num, 1:latide_out-latids_out, lonide_out-lonids_out+1)&
                    = spread_b(1:da_2d_num, latids_out_loc:latide_out_loc:latstep, lonide_out_loc)
                    spread_b_4d(1:da_2d_num, latide_out-latids_out+1, lonide_out-lonids_out+1)&
                    = spread_b(1:da_2d_num, latide_out_loc, lonide_out_loc)
                endif
            endif

            hash_key = "realfield:" // trim(adjustl(atm_mpas_sceneid)) //":" //&
               trim(adjustl(num_lon_str)) // "x"// trim(adjustl(num_lat_str)) &
               // ":" //trim(adjustl(date_char))
            if ((lonstep .ne. 1) .and. (latstep .ne. 1)) then
                call RedisHmgetf3d(rc, hash_key, 'p', lonids_out, lonide_out, 1, latids_out, latide_out, 1, &
                     mpas_num_lev_start, mpas_num_lev, zstep, 1, ph_ref)
            elseif((lonstep == 1) .and. (latstep == 1)) then
                call RedisHmgetf3d(rc, hash_key, 'p', lonids, lonide, lonstep, latids, latide, latstep, &
                     mpas_num_lev_start, mpas_num_lev, zstep, 1, ph_ref)
            else
                call log_error("obs_stat_interval_lon and obs_stat_interval_lat are wrong" ,__FILE__, __LINE__)
            endif

            do i = 1, latide_out-latids_out+1
                do j = 1, lonide_out-lonids_out+1
                    ! calculate gmcore ph
                    if ((latids_out_loc+(i-1)*latstep+latids-1)<num_lat .and. &
                        (lonids_out_loc+(j-1)*lonstep+lonids-1)<num_lon) then
                        axb_tmp = axb(:, latids_out_loc+(i-1)*latstep, lonids_out_loc+(j-1)*lonstep)
                        axa_tmp = axa(:, latids_out_loc+(i-1)*latstep, lonids_out_loc+(j-1)*lonstep)
                        spread_b_tmp = spread_b(:, latids_out_loc+(i-1)*latstep, lonids_out_loc+(j-1)*lonstep)
                        spread_a_tmp = spread_a(:, latids_out_loc+(i-1)*latstep, lonids_out_loc+(j-1)*lonstep)
                    elseif ((latids_out_loc+(i-1)*latstep+latids-1)>=num_lat .and. &
                        (lonids_out_loc+(j-1)*lonstep+lonids-1)<num_lon) then
                        axb_tmp = axb(:, latide_out_loc, lonids_out_loc+(j-1)*lonstep)
                        axa_tmp = axa(:, latide_out_loc, lonids_out_loc+(j-1)*lonstep)
                        spread_b_tmp = spread_b(:, latide_out_loc, lonids_out_loc+(j-1)*lonstep)
                        spread_a_tmp = spread_a(:, latide_out_loc, lonids_out_loc+(j-1)*lonstep)
                    elseif ((latids_out_loc+(i-1)*latstep+latids-1)<num_lat .and. &
                        (lonids_out_loc+(j-1)*lonstep+lonids-1)>=num_lon) then
                        axb_tmp = axb(:, latids_out_loc+(i-1)*latstep, lonide_out_loc)
                        axa_tmp = axa(:, latids_out_loc+(i-1)*latstep, lonide_out_loc)
                        spread_b_tmp = spread_b(:, latids_out_loc+(i-1)*latstep, lonide_out_loc)
                        spread_a_tmp = spread_a(:, latids_out_loc+(i-1)*latstep, lonide_out_loc)
                    elseif ((latids_out_loc+(i-1)*latstep+latids-1)>=num_lat .and. &
                        (lonids_out_loc+(j-1)*lonstep+lonids-1)>=num_lon) then
                        axb_tmp = axb(:, latide_out_loc, lonide_out_loc)
                        axa_tmp = axa(:, latide_out_loc, lonide_out_loc)
                        spread_b_tmp = spread_b(:, latide_out_loc, lonide_out_loc)
                        spread_a_tmp = spread_a(:, latide_out_loc, lonide_out_loc)
                    endif

                    do k = 1, num_lev+1
                        ph_lev_axb(k) = da_vert_coord_calc_ph_lev(k, axb_tmp(ps_pos))
                        ph_lev_axa(k) = da_vert_coord_calc_ph_lev(k, axa_tmp(ps_pos))
                    enddo

                    do k = 1, num_lev
                        ph_axb(k) = 0.5d0 * (ph_lev_axb(k) + ph_lev_axb(k+1))
                        ph_axa(k) = 0.5d0 * (ph_lev_axa(k) + ph_lev_axa(k+1))
                    enddo
!!!!!!!!debug
!    if (pid == 101 .and. i == 20 .and. j == 95) then
!        print *, "ph_axb", ph_axb, "ph_axa", ph_axa
!        print *, "ph_ref", ph_ref(:, i, j)
!    endif
!!!!!!!!!!!!
                    ! interp
!                    m = 1
!                    n = 1
                    do k = 1, (mpas_num_lev - mpas_num_lev_start)/zstep+1
                        do ik = 1, num_lev
                            if (ph_ref(k, i, j) - ph_axb(ik) > 0.00001d0) m = ik+1
                        end do
!                        do while (m <= num_lev .and. ph_ref(k, i, j) - ph_axb(m) > 0.00001d0)
!                            m = m + 1
!                        end do
                        ph_pos_axb(k) = m - 1
                        do ik = 1, num_lev
                            if ( ph_ref(k, i, j) - ph_axa(ik) > 0.00001d0) n = ik+1
                        end do
!                        do while (n <= num_lev .and. ph_ref(k, i, j) - ph_axa(n) > 0.00001d0)
!                            n = n + 1
!                        end do
                        ph_pos_axa(k) = n - 1
                    enddo
!!!!!!!!debug
!     if (pid == 101 .and. i == 20 .and. j == 95) then
!        print *, "ph_pos_axb", ph_pos_axb, "ph_pos_axa", ph_pos_axa
!    endif
!!!!!!!!!!!!
                    do k = 1, (mpas_num_lev - mpas_num_lev_start)/zstep+1
                        if (ph_pos_axb(k) >= 1 .and. ph_pos_axb(k) < num_lev) then
                            do var = 1, da_3d_num
                                tl = da_2d_num + num_lev*(var-1)+ph_pos_axb(k)
                                tk = da_2d_num + mpas_num_lev*(var-1)+mpas_num_lev_start+(k-1)*zstep
                                t3 = (log(ph_axb(ph_pos_axb(k)+1)) - log(ph_ref(k, i, j))) / &
                                     (log(ph_axb(ph_pos_axb(k)+1)) - log(ph_axb(ph_pos_axb(k))))
                                axb_4d(tk, i, j) = axb_tmp(tl) * t3 + axb_tmp(tl+1) * (1 - t3)
                                spread_b_4d(tk, i, j) = spread_b_tmp(tl) * t3 + spread_b_tmp(tl+1) * (1 - t3)
                            enddo
                        elseif (ph_pos_axb(k) == 0) then
                            do var = 1, da_3d_num
                                tl = da_2d_num + num_lev*(var-1)+1
                                tk = da_2d_num + mpas_num_lev*(var-1)+mpas_num_lev_start+(k-1)*zstep
                                axb_4d(tk, i, j) = axb_tmp(tl)
                                spread_b_4d(tk, i, j) = spread_b_tmp(tl)
                            enddo
                        elseif (ph_pos_axb(k) >= num_lev) then
                            do var = 1, da_3d_num
                                tl = da_2d_num + num_lev*var
                                tk = da_2d_num + mpas_num_lev*(var-1)+mpas_num_lev_start+(k-1)*zstep
                                axb_4d(tk, i, j) = axb_tmp(tl)
                                spread_b_4d(tk, i, j) = spread_b_tmp(tl)
                            enddo
                        endif

                        if (ph_pos_axa(k) >= 1 .and. ph_pos_axa(k) < num_lev) then
                            do var = 1, da_3d_num
                                tl = da_2d_num + num_lev*(var-1)+ph_pos_axa(k)
                                tk = da_2d_num + mpas_num_lev*(var-1)+mpas_num_lev_start+(k-1)*zstep
                                t3 = (log(ph_axa(ph_pos_axa(k)+1)) - log(ph_ref(k, i, j))) / &
                                     (log(ph_axa(ph_pos_axa(k)+1)) - log(ph_axa(ph_pos_axa(k))))
                                axa_4d(tk, i, j) = axa_tmp(tl) * t3 + axa_tmp(tl+1) * (1 - t3)
                                spread_a_4d(tk, i, j) = spread_a_tmp(tl) * t3 + spread_a_tmp(tl+1) * (1 - t3)
                            enddo
                        elseif (ph_pos_axa(k) == 0) then
                            do var = 1, da_3d_num
                                tl = da_2d_num + num_lev*(var-1)+1
                                tk = da_2d_num + mpas_num_lev*(var-1)+mpas_num_lev_start+(k-1)*zstep
                                axa_4d(tk, i, j) = axa_tmp(tl)
                                spread_a_4d(tk, i, j) = spread_a_tmp(tl)
                            enddo
                        elseif (ph_pos_axa(k) >= num_lev) then
                            do var = 1, da_3d_num
                                tl = da_2d_num + num_lev*var
                                tk = da_2d_num + mpas_num_lev*(var-1)+mpas_num_lev_start+(k-1)*zstep
                                axa_4d(tk, i, j) = axa_tmp(tl)
                                spread_a_4d(tk, i, j) = spread_a_tmp(tl)
                            enddo
                        endif
                    enddo
!!!!!!!!!!debug
!    if (pid == 101 .and. i == 20 .and. j == 95) then
!        print *, "axb_4d", axb_4d(:,i,j), "axa_4d", axa_4d(:,i,j)
!        print *, "axb",i,j, axb(:,i,j), "axa", axa(:,i,j)
!    endif
!!!!!!!!!!!!!
                enddo
            enddo
        endif

        ! write into nc
        do i=1, da_2d_num
            start2=[lonids_out,latids_out]
            count2=[lonide_out-lonids_out+1,latide_out-latids_out+1]
            do ii = 1,lonide_out-lonids_out+1
                do jj = 1,latide_out-latids_out+1
                    tmp_2d(ii,jj) = axa_4d(i, jj, ii) - axb_4d(i, jj, ii)
                end do
            end do
            call fiona_output(fid0, trim(adjustl(da_var_name_list(i))), tmp_2d, start=start2, count=count2)

            do ii = 1,lonide_out-lonids_out+1
                do jj = 1,latide_out-latids_out+1
                    tmp_2d(ii,jj) = axb_4d(i, jj, ii)
                end do
            end do
            call fiona_output(fid1, trim(adjustl(da_var_name_list(i))), tmp_2d, start=start2, count=count2)

            do ii = 1,lonide_out-lonids_out+1
                do jj = 1,latide_out-latids_out+1
                    tmp_2d(ii,jj) = spread_b_4d(i, jj, ii)
                end do
            end do
            call fiona_output(fid1, 'spread_b_'//trim(adjustl(da_var_name_list(i))), tmp_2d, start=start2, count=count2)

            do ii = 1,lonide_out-lonids_out+1
                do jj = 1,latide_out-latids_out+1
                    tmp_2d(ii,jj) = spread_a_4d(i, jj, ii)
                end do
            end do
            call fiona_output(fid1, 'spread_a_'//trim(adjustl(da_var_name_list(i))), tmp_2d, start=start2, count=count2)
        end do

        do i=1, da_3d_num
            start1=[lonids_out,latids_out,1]
            count1=[lonide_out-lonids_out+1,latide_out-latids_out+1,num_lev_out]
            do ii = 1,lonide_out-lonids_out+1
                do jj = 1,latide_out-latids_out+1
                    k=1
                    do kk = mpas_num_lev_start, mpas_num_lev, zstep
                        tmp_3d(ii,jj,k) = axa_4d(da_2d_num + (i - 1)*mpas_num_lev + kk, jj, ii) &
                                        - axb_4d(da_2d_num + (i - 1)*mpas_num_lev + kk, jj, ii)
                        k = k + 1
                    end do
                end do
            end do
            call fiona_output(fid0, trim(adjustl(da_var_name_list(da_2d_num+i))), tmp_3d, start=start1, count=count1)

            do ii = 1,lonide_out-lonids_out+1
                do jj = 1,latide_out-latids_out+1
                    k=1
                    do kk = mpas_num_lev_start, mpas_num_lev, zstep
                        tmp_3d(ii,jj,k) = axb_4d(da_2d_num + (i - 1)*mpas_num_lev + kk, jj, ii)
                        k = k + 1

                    end do
                end do
            end do
            call fiona_output(fid1, trim(adjustl(da_var_name_list(da_2d_num+i))), tmp_3d, start=start1, count=count1)

            do ii = 1,lonide_out-lonids_out+1
                do jj = 1,latide_out-latids_out+1
                    k=1
                    do kk = mpas_num_lev_start, mpas_num_lev, zstep
                        tmp_3d(ii,jj,k) = spread_b_4d(da_2d_num + (i - 1)*mpas_num_lev + kk, jj, ii)
                        k = k + 1

                    end do
                end do
            end do
            call fiona_output(fid1, 'spread_b_'//trim(adjustl(da_var_name_list(da_2d_num+i))), tmp_3d, &
                              start=start1, count=count1)

            do ii = 1,lonide_out-lonids_out+1
                do jj = 1,latide_out-latids_out+1
                    k=1
                    do kk = mpas_num_lev_start, mpas_num_lev, zstep
                        tmp_3d(ii,jj,k) = spread_a_4d(da_2d_num + (i - 1)*mpas_num_lev + kk, jj, ii)
                        k = k + 1

                    end do
                end do
            end do
            call fiona_output(fid1, 'spread_a_'//trim(adjustl(da_var_name_list(da_2d_num+i))), tmp_3d, &
                              start=start1, count=count1)
        end do

        call fiona_end_output(fid0)
        call fiona_end_output(fid1)
        deallocate(levdeg)
        deallocate(tmp_3d)
        deallocate(tmp_2d)
        deallocate(axa_4d)
        deallocate(axb_4d)
        deallocate(spread_a_4d)
        deallocate(spread_b_4d)
        deallocate(ph_ref)
        call log_print_root("da netcdf output done!")

  end subroutine da_output_nc

end module da_output_mod

