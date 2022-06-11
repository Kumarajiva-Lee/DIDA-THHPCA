module damp_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use parallel_mod
  use parallel_types_mod
  use block_mod
  use operators_mod
  use div_damp_mod
  use vor_damp_mod
  use smag_damp_mod
  use laplace_damp_mod
  use filter_mod
  use pa_mod
  use member_mod

  implicit none

  private

  public damp_init
  public damp_run
  public damp_final


contains

  subroutine damp_init(blocks)

    type(block_type), intent(in) :: blocks(:)

    call vor_damp_init(blocks)
    call div_damp_init(blocks)
    call laplace_damp_init()

  end subroutine damp_init

  subroutine damp_final()

    call vor_damp_final()
    call div_damp_final()
    call laplace_damp_final()

  end subroutine damp_final

  subroutine damp_run(dt, new, blocks)

    real(8), intent(in) :: dt
    integer, intent(in) :: new
    type(block_type), intent(inout) :: blocks(:)

    integer iblk,cyc

    do iblk = 1, size(blocks)
      if (use_div_damp) then
        do cyc = 1 , div_damp_cycles
          call div_damp_run(blocks(iblk), dt, blocks(iblk)%state(new,ivector))
        end do
      end if
      if (use_vor_damp) then
        do cyc = 1, div_damp_cycles
          call vor_damp_run(blocks(iblk), dt, blocks(iblk)%state(new,ivector))
        end do
      end if
      if (use_smag_damp) then
        call smag_damp_run(blocks(iblk), dt, blocks(iblk)%tend(new,ivector), blocks(iblk)%state(new,ivector))
      end if 
      if (use_div_damp .or. use_vor_damp) then
        associate (block => blocks(iblk)                , &
          state => blocks(iblk)%state(new,ivector)      , &
          u     => blocks(iblk)%state(new,ivector)%u    , &
          v     => blocks(iblk)%state(new,ivector)%v    , &
          u_f   => blocks(iblk)%state(new,ivector)%u_f  , &
          v_f   => blocks(iblk)%state(new,ivector)%v_f)
#ifdef Detail_Time
        call Get_Start_Time(tran_time_start)
#endif
        
        call state%async(async_u)%wait()

#ifdef Detail_Time
        call Get_End_Time(tran_time_end)
        tran_time = tran_time + tran_time_end - tran_time_start
#endif
      
        call filter_on_lon_edge(block, u, u_f)

#ifdef Detail_Time
        call Get_Start_Time(tran_time_start)
#endif

        call fill_halo_member(block, u_f, full_lon=.false., full_lat=.true., full_lev=.true.,async=state%async(async_u_f))
        call state%async(async_u_f)%wait()
        call state%async(async_v)%wait()

#ifdef Detail_Time
        call Get_End_Time(tran_time_end)
        tran_time = tran_time + tran_time_end - tran_time_start
#endif

        call filter_on_lat_edge(block, v, v_f)

#ifdef Detail_Time
        call Get_Start_Time(tran_time_start)
#endif
      
        
        call fill_halo_member(block, v_f, full_lon=.true., full_lat=.false., full_lev=.true.,async=state%async(async_v_f))
        call state%async(async_v_f)%wait()
        
#ifdef Detail_Time
        call Get_End_Time(tran_time_end)
        tran_time = tran_time + tran_time_end - tran_time_start
#endif
        end associate
      end if
    end do

  end subroutine damp_run

end module damp_mod
