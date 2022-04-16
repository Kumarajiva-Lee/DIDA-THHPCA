module dcmip31_test_mod

    use flogger
    use namelist_mod
    use const_mod
    use parallel_mod
    use parallel_types_mod
    use block_mod
    use vert_coord_mod
    use formula_mod
    use operators_mod
    use diag_state_mod
    use member_mod
  
    implicit none
  
    private
  
    public dcmip31_test_set_params
    public dcmip31_test_set_ic
  
    real(r8), parameter :: X    = 125.0_r8
    real(r8), parameter :: u0   = 20.0_r8      ! m s-1
    real(r8), parameter :: teq  = 300.0_r8     ! K
    real(r8), parameter :: peq  = 1.0e5_r8     ! Pa
    real(r8), parameter :: d    = 5000.0_r8    ! m
    real(r8), parameter :: lonc = 2 * pi / 3
    real(r8), parameter :: latc = 0.0_r8
    real(r8), parameter :: dpt  = 1.0_r8       ! K
    real(r8), parameter :: lz   = 20000.0_r8   ! m
    real(r8), parameter :: N    = 0.01_r8      ! s-1
    real(r8), parameter :: N2   = N**2
    real(r8)            :: t0
    real(r8), parameter :: p0   = 1.0e5_r8     ! Pa
  
  contains
  
    subroutine dcmip31_test_set_params()
  
      omega = 0.0_r8
      radius = radius / X
  
    end subroutine dcmip31_test_set_params
  
    subroutine dcmip31_test_set_ic(block)
  
      type(block_type), intent(inout), target :: block
  
      integer i, j, k, im
      real(r8) ts, cos_2lat, r, local_z, local_ztop
  
      associate (mesh   => block%mesh           , &
                 state  => block%state(1,ivector)       , &
                 u      => block%state(1,ivector)%u     , &
                 v      => block%state(1,ivector)%v     , &
                 w      => block%state(1,ivector)%w_lev , &
                 phs    => block%state(1,ivector)%phs   , &
                 ph_lev => block%state(1,ivector)%ph_lev, &
                 ph     => block%state(1,ivector)%ph    , &
                 pt     => block%state(1,ivector)%pt    , &
                 t      => block%state(1,ivector)%t     , &
                 gz_lev => block%state(1,ivector)%gz_lev, &
                 gz     => block%state(1,ivector)%gz    , &
                 gzs    => block%static(ivector)%gzs)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              u(:,i,j,k) = u0 * mesh%full_cos_lat(j)
            end do
          end do
        end do
        call fill_halo_member(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)
  
        v = 0.0_r8
        gzs = 0.0_r8
  
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          cos_2lat = cos(2 * mesh%full_lat(j))
          ts = t0 + (teq - t0) * exp(-u0 * N2 / (4 * g**2) * (u0 + 2 * omega * radius) * (cos_2lat - 1))
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            phs(:,i,j) = peq * exp(u0 / (4 * t0 * Rd) * (u0 + 2 * omega * radius) * (cos_2lat - 1)) * &
                       (ts / teq)**(1 / Rd_o_cp)
          end do
        end do
        call fill_halo_member(block, phs, full_lon=.true., full_lat=.true.,async=state%async(async_phs))
        call state%async(async_phs)%wait()
  
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              do im = 1 , member_num
                ph_lev(im,i,j,k) = vert_coord_calc_ph_lev(k, phs(im,i,j))
              end do
            end do
          end do
        end do
        call fill_halo_member(block, ph_lev, full_lon=.true., full_lat=.true., full_lev=.false.,async=state%async(async_ph_lev))
        call state%async(async_ph_lev)%wait()
  
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              ph(:,i,j,k) = 0.5d0 * (ph_lev(:,i,j,k) + ph_lev(:,i,j,k+1))
            end do
          end do
        end do
        call fill_halo_member(block, ph, full_lon=.true., full_lat=.true., full_lev=.true.,async=state%async(async_ph))
        call state%async(async_ph)%wait()
  
        if (nonhydrostatic) then
          w = 0.0_r8
          do k = mesh%half_lev_ibeg, mesh%half_lev_iend
            do j = mesh%full_lat_ibeg, mesh%full_lat_iend
              cos_2lat = cos(2 * mesh%full_lat(j))
              ts = t0 + (teq - t0) * exp(-u0 * N2 / (4 * g**2) * (u0 + 2 * omega * radius) * (cos_2lat - 1))
              do i = mesh%full_lon_ibeg, mesh%full_lon_iend
                gz_lev(:,i,j,k) = - g**2 / N2 * log(ts / t0 * ((ph_lev(:,i,j,k) / phs(:,i,j))**Rd_o_cp - 1) + 1)
              end do
            end do
          end do
          call fill_halo_member(block, gz_lev, full_lon=.true., full_lat=.true., full_lev=.false.,async=state%async(async_gz_lev))
          call state%async(async_gz_lev)%wait()
        end if
  
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            cos_2lat = cos(2 * mesh%full_lat(j))
            ts = t0 + (teq - t0) * exp(-u0 * N2 / (4 * g**2) * (u0 + 2 * omega * radius) * (cos_2lat - 1))
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              pt(:,i,j,k) = ts * (p0 / phs(:,i,j))**Rd_o_cp / (ts / t0 * ((ph(:,i,j,k) / phs(:,i,j))**Rd_o_cp - 1) + 1)
            end do
          end do
        end do
  
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              do im = 1 , member_num
                ! Perturbation
                local_z = 0.5_r8 * (gz_lev(im,i,j,k+1) + gz_lev(im,i,j,k)) / g
                local_ztop = gz_lev(im,i,j,mesh%half_lev_ibeg) / g 
                r = radius * acos(sin(latc) * mesh%full_sin_lat(j) + cos(latc) * mesh%full_cos_lat(j) * cos(mesh%full_lon(i) - lonc))
                pt(im,i,j,k) = pt(im,i,j,k) + + dpt * d**2 / (d**2 + r**2) * sin(pi * local_z / local_ztop)
              end do
            end do
          end do
        end do
        call fill_halo_member(block, pt, full_lon=.true., full_lat=.true., full_lev=.true.,async=state%async(async_pt))
        call state%async(async_pt)%wait()
      end associate
  
    end subroutine dcmip31_test_set_ic
  
  end module dcmip31_test_mod