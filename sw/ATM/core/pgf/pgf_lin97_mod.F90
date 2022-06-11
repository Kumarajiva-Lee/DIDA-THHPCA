module pgf_lin97_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use process_mod
  use parallel_mod
  use parallel_types_mod
  use block_mod
  use pa_mod

  implicit none

  
  external :: slave_pgf_run , slave_pgf_run_north , slave_pgf_run_south

  type var_pgf
    integer :: v0,v1,v2,v3,v4,v5,v6,v7
    real*8  :: v8
    integer*8 :: va,vb,vc,vd
  end type var_pgf

contains

  subroutine pgf_lin97_prepare(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

  end subroutine pgf_lin97_prepare


  subroutine pgf_lin97_run(block, state, tend)

    type(block_type), intent(inout), target :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend

    type(var_pgf) :: vp
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

#ifdef Detail_Time
        call Get_Start_Time(tran_time_start)
#endif
        call state%async(async_gz_lev)%wait()
        call state%async(async_ph_lev)%wait()
#ifdef Detail_Time
        call Get_End_Time(tran_time_end)
        tran_time = tran_time + tran_time_end - tran_time_start
#endif

#ifdef Detail_Time
        call Get_Start_Time(cal_time_start)
#endif

        vp%v0 = mesh%full_lon_ibeg
        vp%v1 = mesh%full_lon_iend
        vp%v2 = mesh%full_lat_ibeg
        vp%v3 = mesh%full_lat_iend
        vp%v4 = mesh%full_lev_ibeg
        vp%v5 = mesh%full_lev_iend
        vp%v6 = proc%lon_halo_width
        vp%v7 = 2  
        vp%v8 = Rd_o_cp
        
        vp%va = loc(ph_lev)
        vp%vb = loc(gz_lev)
        vp%vc = loc(pgf_lon)
        vp%vd = loc(pgf_lat)
        
        
        if (mesh%has_north_pole()) then
          call athread_spawn(slave_pgf_run_north, vp)
          call athread_join()
        else if (mesh%has_south_pole()) then
          call athread_spawn(slave_pgf_run_south, vp)
          call athread_join()
        else
          call athread_spawn(slave_pgf_run, vp)
          call athread_join()
        end if
          


#ifdef Detail_Time
        call Get_End_Time(cal_time_end)
        cal_time = cal_time + cal_time_end - cal_time_start
        lin_time = lin_time + cal_time_end - cal_time_start
#endif
      end if
    end associate

  end subroutine pgf_lin97_run

end module pgf_lin97_mod