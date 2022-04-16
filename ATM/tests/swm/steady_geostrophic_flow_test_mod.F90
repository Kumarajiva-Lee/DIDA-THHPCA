module steady_geostrophic_flow_test_mod

  use flogger
  use const_mod
  use parallel_mod
  use parallel_types_mod
  use block_mod
  use state_mod
  use member_mod

  implicit none

  private

  public steady_geostrophic_flow_test_set_ic

  real(r8), parameter :: alpha = 0.0_r8
  real(r8)            :: u0
  real(r8), parameter :: gz0 = 2.94e4_r8 ! m2 s-2

contains

  subroutine steady_geostrophic_flow_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    real(r8) cos_lat, sin_lat, cos_lon, sin_lon, cos_alpha, sin_alpha
    integer i, j
    type(mesh_type), pointer :: mesh
    type(state_type), pointer :: state

    u0 = 2.0_r8 * pi * radius / (12.0_r8 * 86400.0_r8)

    mesh => block%mesh
    state => block%state(1,ivector)

    cos_alpha = cos(alpha)
    sin_alpha = sin(alpha)
    block%static(ivector)%gzs(:,:,:) = 0.0

    do j = block%state(1,ivector)%mesh%full_lat_ibeg, block%state(1,ivector)%mesh%full_lat_iend
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = block%state(1,ivector)%mesh%half_lon_ibeg, block%state(1,ivector)%mesh%half_lon_iend
        cos_lon = mesh%half_cos_lon(i)
        block%state(1,ivector)%u(:,i,j,1) = u0 * (cos_lat * cos_alpha + cos_lon * sin_lat * sin_alpha)
      end do
    end do
    call fill_halo_member(block, block%state(1,ivector)%u, full_lon=.false., full_lat=.true.,async=state%async(async_u))

    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        sin_lon = mesh%full_cos_lon(i)
        block%state(1,ivector)%v(:,i,j,1) = - u0 * sin_lon * sin_alpha
      end do
    end do
    call fill_halo_member(block, block%state(1,ivector)%v, full_lon=.true., full_lat=.false.,async=state%async(async_v))

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        cos_lon = mesh%full_cos_lon(i)
        block%state(1,ivector)%gz(:,i,j,1) = gz0 - (radius * omega * u0 + u0**2 * 0.5) * (sin_lat * cos_alpha - cos_lon * cos_lat * sin_alpha)**2
      end do
    end do
    call fill_halo_member(block, block%state(1,ivector)%gz, full_lon=.true., full_lat=.true.,async=state%async(async_gz))

    call state%async(async_u)%wait()
    call state%async(async_v)%wait()
    call state%async(async_gz)%wait()

  end subroutine steady_geostrophic_flow_test_set_ic

end module steady_geostrophic_flow_test_mod
