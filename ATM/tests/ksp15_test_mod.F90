module ksp15_test_mod

    use namelist_mod
    use const_mod
    use parallel_mod
    use parallel_types_mod
    use block_mod
    use vert_coord_mod
    use formula_mod
    use member_mod
  
    implicit none
  
    private
  
    public ksp15_test_set_params
    public ksp15_01_test_set_ic
    public ksp15_02_test_set_ic
  
    real(r8), parameter :: teq  = 300.0_r8     ! K
    real(r8), parameter :: peq  = 1.0e5_r8     ! Pa
    real(r8), parameter :: ueq  = 20.0_r8      ! m s-1
    real(r8), parameter :: X    = 166.7_r8
    real(r8)            :: h0
    real(r8), parameter :: d0   = 5000.0_r8    ! m
    real(r8), parameter :: xi0  = 4000.0_r8    ! m
    real(r8), parameter :: lonc = pi
    real(r8), parameter :: latc = 0.0_r8
  
  contains
  
    subroutine ksp15_test_set_params()
  
      omega = 0.0_r8
      radius = radius / X
  
    end subroutine ksp15_test_set_params
  
    subroutine ksp15_01_test_set_ic(block)
  
      type(block_type), intent(inout), target :: block
  
      real(r8) dlon, r0
      integer i, j, k, im

      h0 = 25.0_r8
  
      associate (mesh   => block%mesh           , &
                 state  => block%state(1,ivector)       , &
                 u      => block%state(1,ivector)%u     , &
                 v      => block%state(1,ivector)%v     , &
                 t      => block%state(1,ivector)%t     , &
                 pt     => block%state(1,ivector)%pt    , &
                 gzs    => block%static(ivector) %gzs   , &
                 phs    => block%state(1,ivector)%phs   , &
                 ph     => block%state(1,ivector)%ph    , &
                 ph_lev => block%state(1,ivector)%ph_lev, &
                 gz_lev => block%state(1,ivector)%gz_lev)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              u(:,i,j,k) = ueq * mesh%full_cos_lat(j)
            end do
          end do
        end do
        call fill_halo_member(block, u, full_lon=.false., full_lat=.true., full_lev=.true.,async=state%async(async_u))
  
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            dlon = abs(mesh%full_lon(i) - lonc)
            dlon = min(pi2 - dlon, dlon)
            r0 = radius * dlon
            gzs(:,i,j) = g * h0 * exp(-r0**2 / d0**2) * cos(pi * r0 / xi0)**2 * mesh%full_cos_lat(j)
          end do
        end do
        call fill_halo_member(block, gzs, full_lon=.true., full_lat=.true.,async=state%async(async_gzs))
  
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            phs(:,i,j) = peq * exp(-0.5_r8 * ueq**2 / Rd / teq * mesh%full_sin_lat(j)**2 - gzs(:,i,j) / Rd / teq)
          end do
        end do
        call fill_halo_member(block, phs, full_lon=.true., full_lat=.true.,async=state%async(async_phs))
  
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
  
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              ph(:,i,j,k) = 0.5d0 * (ph_lev(:,i,j,k) + ph_lev(:,i,j,k+1))
            end do
          end do
        end do
        call fill_halo_member(block, ph, full_lon=.true., full_lat=.true., full_lev=.true.,async=state%async(async_ph))
  
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              do im = 1 , member_num
                t(im,i,j,k) = teq
                pt(im,i,j,k) = potential_temperature(t(im,i,j,k), ph(im,i,j,k))
              end do
            end do
          end do
        end do
        call fill_halo_member(block, t , full_lon=.true., full_lat=.true., full_lev=.true.,async=state%async(async_t))
        call fill_halo_member(block, pt, full_lon=.true., full_lat=.true., full_lev=.true.,async=state%async(async_pt))
  
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              gz_lev(:,i,j,k) = Rd * teq * log(peq / ph_lev(:,i,j,k)) - 0.5_r8 * ueq**2 * mesh%full_sin_lat(j)**2
            end do
          end do
        end do
        call fill_halo_member(block, gz_lev, full_lon=.true., full_lat=.true., full_lev=.false.,async=state%async(async_gz_lev))

        call state%async(async_u)%wait()
        call state%async(async_gzs)%wait()
        call state%async(async_phs)%wait()
        call state%async(async_ph_lev)%wait()
        call state%async(async_ph)%wait()
        call state%async(async_t)%wait()
        call state%async(async_pt)%wait()
        call state%async(async_gz_lev)%wait()
        
      end associate
  
    end subroutine ksp15_01_test_set_ic

    subroutine ksp15_02_test_set_ic(block)

      type(block_type), intent(inout), target :: block
  
      real(r8) r
      integer i, j, k, im
  
      h0 = 250.0_r8
  
      associate (mesh   => block%mesh           , &
                 state  => block%state(1,ivector)       , &
                 u      => block%state(1,ivector)%u     , &
                 v      => block%state(1,ivector)%v     , &
                 t      => block%state(1,ivector)%t     , &
                 pt     => block%state(1,ivector)%pt    , &
                 gzs    => block%static(ivector) %gzs   , &
                 phs    => block%state(1,ivector)%phs   , &
                 ph     => block%state(1,ivector)%ph    , &
                 ph_lev => block%state(1,ivector)%ph_lev, &
                 gz_lev => block%state(1,ivector)%gz_lev)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              u(:,i,j,k) = ueq * mesh%full_cos_lat(j)
            end do
          end do
        end do
        call fill_halo_member(block, u, full_lon=.false., full_lat=.true., full_lev=.true.,async=state%async(async_u))
  
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            r = radius * acos(sin(latc) * mesh%full_sin_lat(j) + cos(latc) * mesh%full_cos_lat(j) * cos(mesh%full_lon(i) - lonc))
            gzs(:,i,j) = g * h0 * exp(-r**2 / d0**2) * cos(pi * r / xi0)**2
          end do
        end do
        call fill_halo_member(block, gzs, full_lon=.true., full_lat=.true.,async=state%async(async_gzs))
  
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            phs(:,i,j) = peq * exp(-0.5_r8 * ueq**2 / Rd / teq * mesh%full_sin_lat(j)**2 - gzs(:,i,j) / Rd / teq)
          end do
        end do
        call fill_halo_member(block, phs, full_lon=.true., full_lat=.true.,async=state%async(async_phs))
  
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
  
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              ph(:,i,j,k) = 0.5d0 * (ph_lev(:,i,j,k) + ph_lev(:,i,j,k+1))
            end do
          end do
        end do
        call fill_halo_member(block, ph, full_lon=.true., full_lat=.true., full_lev=.true.,async=state%async(async_ph))
  
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              do im = 1 , member_num
                t(im,i,j,k) = teq
                pt(im,i,j,k) = potential_temperature(t(im,i,j,k), ph(im,i,j,k))
              end do
            end do
          end do
        end do
        call fill_halo_member(block, t , full_lon=.true., full_lat=.true., full_lev=.true.,async=state%async(async_t))
        call fill_halo_member(block, pt, full_lon=.true., full_lat=.true., full_lev=.true.,async=state%async(async_pt))
  
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              gz_lev(:,i,j,k) = Rd * teq * log(peq / ph_lev(:,i,j,k)) - 0.5_r8 * ueq**2 * mesh%full_sin_lat(j)**2
            end do
          end do
        end do
        call fill_halo_member(block, gz_lev, full_lon=.true., full_lat=.true., full_lev=.false.,async=state%async(async_gz_lev))
        
        call state%async(async_u)%wait()
        call state%async(async_gzs)%wait()
        call state%async(async_phs)%wait()
        call state%async(async_ph_lev)%wait()
        call state%async(async_ph)%wait()
        call state%async(async_t)%wait()
        call state%async(async_pt)%wait()
        call state%async(async_gz_lev)%wait()

      end associate
  
    end subroutine ksp15_02_test_set_ic
  
  end module ksp15_test_mod