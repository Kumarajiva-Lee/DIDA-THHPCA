module cross_pole_flow_test_mod
  
  use flogger
  use const_mod
  use parallel_mod
  use parallel_types_mod
  use block_mod
  use state_mod
  use member_mod

  implicit none

  private

  public cross_pole_flow_test_set_ic

  real, parameter :: v0 = 20 ! m s-1
  real, parameter :: gz0 = 5.7684e4 ! m2 s-2

contains
  
  subroutine cross_pole_flow_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    real(r8) cos_lat, sin_lat, cos_lon, sin_lon
    integer i, j
    type(mesh_type), pointer :: mesh
    type(state_type), pointer :: state

    mesh => block%mesh
    state => block%state(1,ivector)

    block%static(ivector)%gzs(:,:,:) = 0.0

    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        sin_lon = mesh%half_sin_lon(i)
        block%state(1,ivector)%u(:,i,j,1) = -v0 * sin_lon * sin_lat * (4.0 * cos_lat**2 - 1.0)
      end do
    end do
    call fill_halo_member(block, block%state(1,ivector)%u, full_lon=.false., full_lat=.true.,async=state%async(async_u))

    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      sin_lat = mesh%half_sin_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        cos_lon = mesh%full_cos_lon(i)
        block%state(1,ivector)%v(:,i,j,1) = v0 * sin_lat**2 * cos_lon 
      end do 
    end do 
    call fill_halo_member(block, block%state(1,ivector)%v, full_lon=.true., full_lat=.false.,async=state%async(async_v))

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        sin_lon = mesh%full_sin_lon(i)
        block%state(1,ivector)%gz(:,i,j,1) = gz0 + 2 * radius * omega * v0 * sin_lat**3 * cos_lat * sin_lon 
      end do 
    end do 
    call fill_halo_member(block, block%state(1,ivector)%gz, full_lon=.true., full_lat=.true.,async=state%async(async_gz))

    call state%async(async_u)%wait()
    call state%async(async_v)%wait()
    call state%async(async_gz)%wait()

  end subroutine cross_pole_flow_test_set_ic

end module cross_pole_flow_test_mod
