module pgf_lin97_mod

  use flogger
  use const_mod
  use namelist_mod
  use parallel_mod
  use parallel_types_mod
  use block_mod
  use pa_mod

  implicit none

contains

  subroutine pgf_lin97_prepare(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

  end subroutine pgf_lin97_prepare

  subroutine pgf_lin97_run(block, state, tend)

    type(block_type), intent(inout), target :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend

    real(r8) dph1(member_num), dph2(member_num), dgz1(member_num), dgz2(member_num), dp1(member_num), dp2(member_num), dpdph(member_num)
    integer i, j, k, move, rf

    !                    o
    !                   /|
    !                  / |
    !                 /  |
    !                /   |
    !   o-----------/------------o
    !   |          /|            |
    !   |         / |            |
    !   |        /  |            |
    !   |       /   |            |
    !   |      o    |            |
    !   o------|    -------------o
    !          |   /
    !          |  /
    !          | /
    !          |/
    !          o

    ! write(*,*) "pgf_lon"

#ifdef Detail_Time
  call Get_Time_Pa(stencil_time_start)
#endif

    associate (mesh          => block%mesh         , & ! in
               ph_lev        => state%ph_lev       , & ! in
               gz_lev        => state%gz_lev       , & ! in
               p_lev         => state%p_lev        , & ! in    ! For nonhydrostatic
               rhod_lon      => state%rhod_lon     , & ! in    !
               rhod_lat      => state%rhod_lat     , & ! in    !
               p_lev_lon     => state%p_lev_lon    , & ! in    !
               p_lev_lat     => state%p_lev_lat    , & ! in    !
               m_lon         => state%m_lon        , & ! in    !
               m_lat         => state%m_lat        , & ! in    !
               pgf_lon       => tend%pgf_lon       , & ! out
               pgf_lat       => tend%pgf_lat)          ! out
      if (hydrostatic) then

        call state%async(async_gz_lev)%wait()
        call state%async(async_ph_lev)%wait()

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          !
          !   4             3
          ! i,j,k        i+1,j,k
          !   o-------------o
          !   |             |
          !   |             |
          !   |    i,j,k    |
          !   |             |
          !   |             |
          !   o-------------o
          ! i,j,k+1      i+1,j,k+1  --> east
          !   1             2
          !
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              dph1 = ph_lev(:,i+1,j,k+1)**Rd_o_cp - ph_lev(:,i  ,j,k  )**Rd_o_cp ! 2 - 4
              dph2 = ph_lev(:,i  ,j,k+1)**Rd_o_cp - ph_lev(:,i+1,j,k  )**Rd_o_cp ! 1 - 3
              dgz1 = gz_lev(:,i  ,j,k+1)          - gz_lev(:,i+1,j,k  )          ! 1 - 3
              dgz2 = gz_lev(:,i  ,j,k  )          - gz_lev(:,i+1,j,k+1)          ! 4 - 2
              pgf_lon(:,i,j,k) = -(dph1 * dgz1 + dph2 * dgz2) / mesh%de_lon(j) / (dph1 + dph2)
            end do
          end do
          !
          !   4             3
          ! i,j,k        i,j+1,k
          !   o-------------o
          !   |             |
          !   |             |
          !   |    i,j,k    |
          !   |             |
          !   |             |
          !   o-------------o
          ! i,j,k+1      i,j+1,k+1  --> north
          !   1             2
          !
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              dph1 = ph_lev(:,i,j+1,k+1)**Rd_o_cp - ph_lev(:,i,j  ,k  )**Rd_o_cp ! 2 - 4
              dph2 = ph_lev(:,i,j  ,k+1)**Rd_o_cp - ph_lev(:,i,j+1,k  )**Rd_o_cp ! 1 - 3
              dgz1 = gz_lev(:,i,j  ,k+1)          - gz_lev(:,i,j+1,k  )          ! 1 - 3
              dgz2 = gz_lev(:,i,j  ,k  )          - gz_lev(:,i,j+1,k+1)          ! 4 - 2
              pgf_lat(:,i,j,k) = -(dph1 * dgz1 + dph2 * dgz2) / mesh%de_lat(j) / (dph1 + dph2)
            end do
          end do
        end do
      else if (nonhydrostatic) then
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              dph1 = ph_lev(:,i+1,j,k+1)**Rd_o_cp - ph_lev(:,i  ,j,k  )**Rd_o_cp ! 2 - 4
              dph2 = ph_lev(:,i  ,j,k+1)**Rd_o_cp - ph_lev(:,i+1,j,k  )**Rd_o_cp ! 1 - 3
              dgz1 = gz_lev(:,i  ,j,k+1)          - gz_lev(:,i+1,j,k  )          ! 1 - 3
              dgz2 = gz_lev(:,i  ,j,k  )          - gz_lev(:,i+1,j,k+1)          ! 4 - 2
              dp1  = p_lev (:,i  ,j,k+1)          - p_lev (:,i+1,j,k  )          ! 1 - 3
              dp2  = p_lev (:,i  ,j,k  )          - p_lev (:,i+1,j,k+1)          ! 4 - 2
              dpdph = (p_lev_lon(:,i,j,k+1) - p_lev_lon(:,i,j,k)) / m_lon(:,i,j,k)
              pgf_lon(:,i,j,k) = -(                             &
                (dph1 * dgz1 + dph2 * dgz2) * dpdph +         &
                (dph1 * dp1  + dph2 * dp2 ) / rhod_lon(:,i,j,k) &
              ) / mesh%de_lon(j) / (dph1 + dph2)
            end do
          end do
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              dph1 = ph_lev(:,i,j+1,k+1)**Rd_o_cp - ph_lev(:,i,j  ,k  )**Rd_o_cp ! 2 - 4
              dph2 = ph_lev(:,i,j  ,k+1)**Rd_o_cp - ph_lev(:,i,j+1,k  )**Rd_o_cp ! 1 - 3
              dgz1 = gz_lev(:,i,j  ,k+1)          - gz_lev(:,i,j+1,k  )          ! 1 - 3
              dgz2 = gz_lev(:,i,j  ,k  )          - gz_lev(:,i,j+1,k+1)          ! 4 - 2
              dp1  = p_lev (:,i,j  ,k+1)          - p_lev (:,i,j+1,k  )          ! 1 - 3
              dp2  = p_lev (:,i,j  ,k  )          - p_lev (:,i,j+1,k+1)          ! 4 - 2
              dpdph = (p_lev_lat(:,i,j,k+1) - p_lev_lat(:,i,j,k)) / m_lat(:,i,j,k)
              pgf_lat(:,i,j,k) = -(                             &
                (dph1 * dgz1 + dph2 * dgz2) * dpdph +         &
                (dph1 * dp1  + dph2 * dp2 ) / rhod_lat(:,i,j,k) &
              ) / mesh%de_lat(j) / (dph1 + dph2)
            end do
          end do
        end do
      end if
    end associate

#ifdef Detail_Time
  call Get_Time_Pa(stencil_time_end)
  stencil_time(12) = stencil_time(12) + stencil_time_end - stencil_time_start
#endif

  end subroutine pgf_lin97_run

end module pgf_lin97_mod