! module mars_cold_run_mod

!     use flogger
!     use namelist_mod, only: topo_file
!     use const_mod
!     use parallel_mod
!     use block_mod
!     use vert_coord_mod
!     use formula_mod
!     use topo_mod
!     use member_mod
  
!     implicit none
  
!     private
  
!     public mars_cold_run_set_ic
  
!     real(r8), parameter :: t0   = 170.0_r8   ! K
!     real(r8), parameter :: ps0  = 6.10_r8    ! Pa
!     real(r8), parameter :: ptop = 0.008_r8   ! Pa
  
!   contains
  
!     subroutine mars_cold_run_set_ic(block)
  
!       type(block_type), intent(inout), target :: block
  
!       integer i, j, k, im
  
!       call topo_read(topo_file)
!       call topo_regrid(block)
  
!       associate (mesh   => block%mesh           , &
!                  state => block%state(1,1)       ,&
!                  u      => block%state(1)%u     , &
!                  v      => block%state(1)%v     , &
!                  t      => block%state(1)%t     , &
!                  pt     => block%state(1)%pt    , &
!                  ph     => block%state(1)%ph    , &
!                  ph_lev => block%state(1)%ph_lev, &
!                  phs    => block%state(1)%phs   , &
!                  gzs    => block%static%gzs)
!       u = 0
!       v = 0
!       t = t0
!       do j = mesh%full_lat_ibeg, mesh%full_lat_iend
!         do i = mesh%full_lon_ibeg, mesh%full_lon_iend
!           phs(:,i,j) = ps0 * exp(gzs(:,i,j) / (g * Rd * t0)) - ptop
!         end do
!       end do
!       call fill_halo_member(block, phs, full_lon=.true., full_lat=.true.,async=state%async(async_phs))
!       call state%async(async_phs)%wait()
  
!       do k = mesh%half_lev_ibeg, mesh%half_lev_iend
!         do j = mesh%full_lat_ibeg, mesh%full_lat_iend
!           do i = mesh%full_lon_ibeg, mesh%full_lon_iend
!             do im = 1 , member_num
!                 ph_lev(i,j,k) = vert_coord_calc_ph_lev(k, phs(i,j))
!             end do
!           end do
!         end do
!       end do
!       call fill_halo(block, ph_lev, full_lon=.true., full_lat=.true., full_lev=.false.,async=state%async(async_ph_lev))
!       call state%async(async_ph_lev)%wait()
  
!       do k = mesh%full_lev_ibeg, mesh%full_lev_iend
!         do j = mesh%full_lat_ibeg, mesh%full_lat_iend
!           do i = mesh%full_lon_ibeg, mesh%full_lon_iend
!             ph(:,i,j,k) = 0.5d0 * (ph_lev(:,i,j,k) + ph_lev(:,i,j,k+1))
!           end do
!         end do
!       end do
!       call fill_halo(block, ph, full_lon=.true., full_lat=.true., full_lev=.true.,async=state%async(async_ph_lev))
!       call state%async(async_ph_lev)%wait()
  
!       do k = mesh%full_lev_ibeg, mesh%full_lev_iend
!         do j = mesh%full_lat_ibeg, mesh%full_lat_iend
!           do i = mesh%full_lon_ibeg, mesh%full_lon_iend
!             do im = 1 , member_num
!                 pt(im,i,j,k) = potential_temperature(t(im,i,j,k), ph(im,i,j,k))
!             end do
!           end do
!         end do
!       end do
!       call fill_halo(block, pt, full_lon=.true., full_lat=.true., full_lev=.true.,async=state%async(async_pt))
!       call state%async(async_pt)%wait()
!       end associate
  
!     end subroutine mars_cold_run_set_ic
  
!   end module mars_cold_run_mod
  