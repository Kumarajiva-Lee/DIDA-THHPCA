module gmcore_mod
 
  use flogger
  use string
  use const_mod
  use namelist_mod
  use parallel_mod
  use parallel_types_mod
  use process_mod
  use datetime
  use time_mod, old => old_time_idx, new => new_time_idx
  use history_mod
  use restart_mod
  use initial_mod
  use prepare_mod
  use tend_mod
  use block_mod
  use diag_state_mod
  use vert_coord_mod
  use time_schemes_mod
  use operators_mod
  use interp_mod
  use debug_mod
  use pgf_mod
  use damp_mod
  use diag_state_mod
  use test_forcing_mod
  use pa_mod
  use member_mod
  

  use mountain_zonal_flow_test_mod
  use rossby_haurwitz_wave_test_mod
  use jet_zonal_flow_test_mod
  use steady_geostrophic_flow_test_mod
  use cross_pole_flow_test_mod
  use shallow_water_waves_test_mod
  use vortex_erosion_test_mod
  use steady_state_test_mod
  use rossby_haurwitz_wave_3d_test_mod
  use mountain_wave_test_mod
  use baroclinic_wave_test_mod
  use held_suarez_test_mod
  use steady_state_pgf_test_mod
  use ksp15_test_mod
  use dcmip31_test_mod
  !use mars_cold_run_mod


  use gmcore_coupler, only: init_gmcore_coupler,gmcore2da_send_coupling,gmcore2da_recv_coupling,clean_gmcore_coupler

  implicit none

  private

  public atm_init
  public atm_run
  public atm_final

  interface
    subroutine splitter_interface(dt, blocks)
      import block_type
      real(8), intent(in) :: dt
      type(block_type), intent(inout) :: blocks(:)
    end subroutine splitter_interface

    subroutine set_ic_interface(block)
      import block_type
      type(block_type), intent(inout), target :: block
    end subroutine set_ic_interface
  end interface

  procedure(splitter_interface), pointer :: splitter
  procedure(set_ic_interface), pointer :: set_ic
  procedure(space_operators_interface), pointer :: operators

contains

  subroutine atm_init(namelist_path , comm , group_id , proc_lon_ibeg , proc_lon_iend , proc_lat_ibeg , proc_lat_iend)

    character(*), intent(in) :: namelist_path
    integer , intent(in) :: comm
    integer , intent(in) , optional :: group_id
    integer , intent(inout) , optional :: proc_lon_ibeg
    integer , intent(inout) , optional :: proc_lon_iend
    integer , intent(inout) , optional :: proc_lat_ibeg
    integer , intent(inout) , optional :: proc_lat_iend

    character(10) time_value, time_units
    real(r8) seconds
    real(8) :: new_second
    integer :: iblk
    integer :: im, iv, now_ensemble
    integer origin_time_step 
    type(datetime_type) origin_start_time
    type(datetime_type) origin_end_time

    ivector = 1
    call parse_namelist(namelist_path)
    call const_init(planet)
    call log_init()
    call global_mesh%init_global(num_lon, num_lat, num_lev, &
                                 lon_halo_width=2, &
                                 lat_halo_width=2)
    call process_init(comm)

    if (is_root_proc()) write(*,*) 'Finish Init'
    
    call process_init_zonal()
    if (is_root_proc()) call print_namelist()
    call vert_coord_init(num_lev, namelist_path)
    call process_create_blocks_stage2
    if (present(group_id)) then
      call pa_init(proc%comm,group_id);
    else
      call pa_init(proc%comm)
    endif
    call time_init()
    !call diag_state_init(proc%blocks)
    if (present(group_id)) then
      call history_init(group_id)
    else
      call history_init()
    endif
    call restart_init()
    call time_scheme_init()
    call pgf_init()
    call interp_init()
    call damp_init(proc%blocks) 
#ifndef Only_Gmcore
    call init_gmcore_coupler(proc%comm, proc%blocks(1)%state(1,1))
#endif

    if (present(proc_lon_ibeg)) proc_lon_ibeg = proc%lon_ibeg
    if (present(proc_lon_iend)) proc_lon_iend = proc%lon_iend
    if (present(proc_lat_ibeg)) proc_lat_ibeg = proc%lat_ibeg
    if (present(proc_lat_iend)) proc_lat_iend = proc%lat_iend

    operators => space_operators

    if (initial_file == 'N/A' .and. .not. restart) then
      select case (test_case)
      case ('pgf_test')
        call steady_state_pgf_test_set_params()
      case ('ksp15_01', 'ksp15_02')
        call ksp15_test_set_params()
      case ('dcmip31')
        call dcmip31_test_set_params()
      end select
    end if

    do iv = 1 , vector_num
      ivector = iv
      if (initial_file /= 'N/A') then
        if (present(group_id)) then
          call initial_read(initial_file, group_id)
        else
          call initial_read(initial_file)
        endif
      else if (topo_file /= 'N/A' .and. bkg_file /= 'N/A') then
        if (present(group_id)) then
          call prepare_run(group_id_in = group_id)
        else
          call prepare_run()
        end if
      else if (restart) then
        if (present(group_id)) then
          call restart_read(group_id)
        else
          call restart_read()
        end if
      else
        select case (test_case)
        case ('steady_state')
          set_ic => steady_state_test_set_ic
        case ('rossby_haurwitz_wave')
          set_ic => rossby_haurwitz_wave_3d_test_set_ic
        case ('mountain_wave')
          set_ic => mountain_wave_test_set_ic
        case ('baroclinic_wave')
          set_ic => baroclinic_wave_test_set_ic
        case ('held_suarez')
          set_ic => held_suarez_test_set_ic
        case ('pgf_test')
          set_ic => steady_state_pgf_test_set_ic
        case ('ksp15_01')
          set_ic => ksp15_01_test_set_ic
        case ('ksp15_02')
          set_ic => ksp15_02_test_set_ic
        case ('dcmip31')
          set_ic => dcmip31_test_set_ic
        case ('mars_cold_run')
          !set_ic => mars_cold_run_set_ic
        case ('mountain_zonal_flow')
          set_ic => mountain_zonal_flow_test_set_ic
        case ('rossby_haurwitz_wave_swm')
          set_ic => rossby_haurwitz_wave_test_set_ic
        case ('jet_zonal_flow')
          set_ic => jet_zonal_flow_test_set_ic
        case ('steady_geostrophic_flow')
          set_ic => steady_geostrophic_flow_test_set_ic
        case ('cross_pole_flow')
          set_ic => cross_pole_flow_test_set_ic
        case ('shallow_water_waves')
          set_ic => shallow_water_waves_test_set_ic
        case ('vortex_erosion')
          !set_ic => vortex_erosion_test_set_ic
        case default
          call log_error('Unknown test case ' // trim(test_case) // '!')
        end select

        do iblk = 1, size(proc%blocks)
          call set_ic(proc%blocks(iblk))
        end do

      end if
    end do

    time_value = split_string(print_interval, ' ', 1)
    time_units = split_string(print_interval, ' ', 2)
    read(time_value, *) seconds
    select case (time_units)
    case ('days')
      seconds = seconds * 86400
    case ('hours')
      seconds = seconds * 3600
    case ('minutes')
      seconds = seconds * 60
    case ('seconds')
      seconds = seconds
    case default
      call log_error('Invalid print interval ' // trim(print_interval) // '!')
    end select

    call time_add_alert('print', seconds=seconds, first_alert=.true.)

#ifndef Only_Gmcore
    new_second = da_in_seconds
    call time_add_alert('with_letkf', seconds=new_second, first_alert=.false.)
#endif

    if (use_create_ensemble) then 
      !fix ivector


      ! origin_start_time = start_time
      ! origin_end_time   = end_time
      ! origin_time_step  = time_step

      ! seconds = 60 * 60 * 1.0
      ! call time_add_alert('spin', seconds =  seconds, first_alert=.false.)
      
      ! seconds = 60 * create_ensemble_interval * 1.0
      ! call time_add_alert('new_assem' , seconds = seconds, first_alert=.false.)

      ! if (is_root_proc()) call log_print('Begin Spin' )

      ! call operators_prepare(proc%blocks, old, dt_in_seconds)
      ! call diagnose(proc%blocks, old , 1)
      ! !call output(proc%blocks, old) 

      ! do while (.not. time_is_alerted('spin'))
      !   call time_integrate(dt_in_seconds, proc%blocks)
      !   call time_advance(dt_in_seconds)
      !   call operators_prepare(proc%blocks, old, dt_in_seconds)
      !   call diagnose(proc%blocks, old, 1)
      !   if (is_root_proc() .and. time_is_alerted('print')) call log_print_diag(curr_time%isoformat())
      !   ! call output(proc%blocks, old)
      ! end do

      ! if (is_root_proc()) call log_print('Finish Spin' )

      ! now_ensemble = 1

      ! if (is_root_proc()) call log_print('Begin Create ensembles')
  
      ! do while (.true.)
      !   call time_integrate(dt_in_seconds, proc%blocks)
      !   call time_advance(dt_in_seconds)
      !   if (time_is_alerted('new_assem')) then
      !      if (is_root_proc()) call log_print('Create ensembles ' // to_str(now_ensemble))
      !      call move_new_ensemble(proc%blocks(1) , now_ensemble)
      !      now_ensemble = now_ensemble + 1
      !      if (now_ensemble > member_num) exit
      !   end if
      !   call operators_prepare(proc%blocks, old, dt_in_seconds)
      !   call diagnose(proc%blocks, old , 1)
      !   if (is_root_proc() .and. time_is_alerted('print')) call log_print_diag(curr_time%isoformat())
      !   !call output(proc%blocks, old)
      ! end do
      
      ! !tmp = proc%blocks(1)%ensemble_position
      ! proc%blocks(1)%state(old)%u = proc%blocks(1)%state(0)%u
      ! proc%blocks(1)%state(old)%v = proc%blocks(1)%state(0)%v
      ! proc%blocks(1)%state(old)%pt = proc%blocks(1)%state(0)%pt
      ! proc%blocks(1)%state(old)%phs = proc%blocks(1)%state(0)%phs
      ! proc%blocks(1)%state(old)%gz = proc%blocks(1)%state(0)%gz
      ! proc%blocks(1)%static%gzs = proc%blocks(1)%store_static%gzs
  
      ! if (is_root_proc()) call log_print('Finish Create ensembles')
  
      ! call time_reset_start_time(origin_start_time)
      ! call time_reset_end_time(origin_end_time)
      ! time_step = origin_time_step

    end if

  end subroutine atm_init

  subroutine atm_run(comm , group_id)

    integer, intent(in) :: comm
    integer , intent(in) , optional :: group_id
    
    integer iter_member, iv
    integer flag , now_ensemble
    integer tmp 
    integer iblk, itime
    logical new_file 



    if (is_root_proc()) call log_notice("Begin atm_run")
    new_file = .true.

    call Get_Time_Pa(total_time_start)

    do iblk = 1, size(proc%blocks)
      do iv = 1 , vector_num
        associate (block => proc%blocks(iblk)     , &
                  mesh  => proc%blocks(iblk)%mesh, &
                  state => proc%blocks(iblk)%state(old,iv))
        if (baroclinic) then 
          call prepare_static(block,iv)
          ! Ensure bottom gz_lev is the same as gzs.
          do itime = lbound(block%state, 1), ubound(block%state, 1)
            block%state(itime,iv)%gz_lev(:,:,:,global_mesh%half_lev_iend) = block%static(iv)%gzs
          end do
        end if
        block%state(old,iv)%u_f = block%state(old,iv)%u
        block%state(old,iv)%v_f = block%state(old,iv)%v
        if (baroclinic) then
          block%state(old,iv)%pt_f = block%state(old,iv)%pt
          block%state(old,iv)%phs_f = block%state(old,iv)%phs
        else
          block%state(old,iv)%gz_f = block%state(old,iv)%gz
        end if
        end associate
      end do
    end do




    do iv = 1 , vector_num
      ivector = iv
      call operators_prepare(proc%blocks, old, dt_in_seconds)
      ! if (nonhydrostatic) call nh_prepare(proc%blocks)
      call diagnose(proc%blocks, old , 1)
      if (present(group_id)) then
        call output(old,group_id,new_file_in = .true.)
      else
        call output(old,new_file_in = .true.)
      endif
    end do

    do while (.not. time_is_finished())
      
      do iv = 1 , vector_num
        ivector = iv
        call time_integrate(dt_in_seconds, proc%blocks)
        if (is_root_proc() .and. time_is_alerted('print')) call log_print_diag(curr_time%isoformat())
#ifndef Only_Gmcore
        !fix ivector
        if (time_is_alerted('with_letkf')) then
          call gmcore2da_send_coupling(comm, proc%blocks(1), proc%blocks(1)%state(new,:))
          call gmcore2da_recv_coupling(comm, proc%blocks(1), proc%blocks(1)%state(new,:))
        end if
#endif
        call time_advance(dt_in_seconds)
        call operators_prepare(proc%blocks, old, dt_in_seconds)
        call diagnose(proc%blocks, old , 1)
        if (present(group_id)) then
          call output(old,group_id)
        else
          call output(old)
        endif
      end do




      

    end do

    !call operators_prepare(proc%blocks, old, dt_in_seconds)
    !call diagnose(proc%blocks, old , 1)

    call Get_Time_Pa(total_time_end)

    if (is_root_proc()) call log_notice("Finish atm_run")
    if (is_root_proc()) call log_notice('Atm_run time cost ' // to_str(total_time_end - total_time_start, 5) // ' seconds.')
    

  end subroutine atm_run

  subroutine atm_final()

#ifndef Only_Gmcore
    call clean_gmcore_coupler()
#endif
    call interp_final()
    call damp_final()
    !call diag_state_final()
    call history_final()
#ifdef Detail_Time
    call pa_final()
#endif
    call process_final()

    if (is_root_proc()) call log_notice("Finish atm_final")

  end subroutine atm_final

  subroutine move_new_ensemble(block , now_ensemble)
    type(block_type), intent(inout) :: block
    integer , intent(in) :: now_ensemble


    ! block%state(    block%ensemble_position)%u(now_ensemble , : , : , :) = block%state(new)%u(1,:,:,:)
    ! block%state(    block%ensemble_position)%v(now_ensemble , : , : , :) = block%state(new)%v(1,:,:,:)
    ! block%store_static%gzs(now_ensemble , : , :) = block%static%gzs(1 , : , :)
    ! block%state(    block%ensemble_position)%gz(now_ensemble , : , : , :) = block%state(new)%gz(1 , : , : , :)
    ! block%state(    block%ensemble_position)%pt(now_ensemble , : , : , :) = block%state(new)%pt(1 , : , : , :) 
    ! block%state(    block%ensemble_position)%phs(now_ensemble , : , : )   = block%state(new)%phs(1 , : , :)
        
  end subroutine move_new_ensemble

  subroutine output(itime , group_id , new_file_in)

    integer, intent(in) :: itime
    integer , intent(in) , optional :: group_id
    logical , intent(in) , optional :: new_file_in

    real(8), save :: time1 = 0, time2
    integer i, j, k, iblk
    real*8 io_time_start , io_time_end
    logical new_file

    if (time_is_alerted('history_write')) then
      if (time_step == 0) time1 = mpi_wtime()
      time2 = mpi_wtime()
      if (time_step /= 0) then
        if (is_root_proc()) call log_notice('Time cost ' // to_str(time2 - time1, 5) // ' seconds.')
        time1 = time2
      end if

      if (present(new_file_in)) then
        new_file = new_file_in
      else
        new_file = .false.
      end if

#ifdef Detail_io_Time
    io_time_start = mpi_wtime()
#endif

      if (output_nc) then
        if (present(group_id)) then
          call history_write_state(proc%blocks, itime,group_id,new_file = new_file)
          if (output_debug) call history_write_debug(proc%blocks, itime ,group_id,new_file = new_file)
        else
          call history_write_state(proc%blocks, itime,new_file = new_file)
          if (output_debug) call history_write_debug(proc%blocks, itime,new_file = new_file)
        endif
      end if

#ifdef Detail_io_Time
    io_time_end = mpi_wtime()
    write((23333 + myid),*) , io_time_end - io_time_start
#endif

    end if
    if (time_is_alerted('restart_write')) then
      if (present(group_id)) then
        call restart_write(itime , group_id)
      else
        call restart_write(itime)
      end if
    end if

  end subroutine output

  subroutine diagnose(blocks, itime , member_output)

    type(block_type), intent(inout), target :: blocks(:)
    integer, intent(in) :: itime
    integer, intent(in) :: member_output

    type(mesh_type), pointer :: mesh
    type(state_type), pointer :: state
    type(static_type), pointer :: static
    integer i, j, k, iblk
    real(r8) tm(member_num), te(member_num), tav(member_num), tpe(member_num)

    tm = 0.0_r8
    te = 0.0_r8
    tav = 0.0_r8
    tpe = 0.0_r8
    do iblk = 1, size(blocks)
      mesh => blocks(iblk)%mesh
      state => blocks(iblk)%state(itime,ivector)
      static => blocks(iblk)%static(ivector)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            tm(:) = tm(:) + state%m(:,i,j,k) * mesh%area_cell(j)
          end do
        end do
      end do

      call state%async(async_mf_lon_n)%wait()
      call state%async(async_u)%wait()
      call state%async(async_mf_lat_n)%wait()
      call state%async(async_v)%wait()

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            te(:) = te(:) + state%mf_lon_n(:,i,j,k) * 0.5_r8 * state%u(:,i,j,k) * mesh%area_lon(j) * 2
          end do
        end do
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            te(:) = te(:) + state%mf_lat_n(:,i,j,k) * 0.5_r8 * state%v(:,i,j,k) * mesh%area_lat(j) * 2
          end do
        end do
      end do
      if (hydrostatic) then
        call state%async(async_t)%wait()
        call state%async(async_phs)%wait()

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              te(:) = te(:) + state%m(:,i,j,k) * cp * state%t(:,i,j,k) * mesh%area_cell(j)
            end do
          end do
        end do
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            te(:) = te(:) + static%gzs(:,i,j) * state%phs(:,i,j) * mesh%area_cell(j)
          end do
        end do
      else if (nonhydrostatic) then
      else
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            te(:) = te(:) + (state%m(:,i,j,1)**2 * g * 0.5_r8 + state%m(:,i,j,1) * static%gzs(:,i,j)) * mesh%area_cell(j)
          end do
        end do
      end if

      call state%async(async_pv)%wait()

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg, mesh%half_lat_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            tav(:) = tav(:) + state%pv(:,i,j,k) * mesh%area_vtx(j)
          end do
        end do
      end do

      call state%async(async_m_lon)%wait()
      call state%async(async_m_lat)%wait()
      call state%async(async_pv_lon)%wait()
      call state%async(async_pv_lat)%wait()

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            tpe(:) = tpe(:) + state%m_lon(:,i,j,k) * state%pv_lon(:,i,j,k)**2 * 0.5_r8 * mesh%area_lon(j)
          end do
        end do
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            tpe(:) = tpe(:) + state%m_lat(:,i,j,k) * state%pv_lat(:,i,j,k)**2 * 0.5_r8 * mesh%area_lat(j)
          end do
        end do
      end do
    end do
    call global_sum(proc%comm, tm)
    call global_sum(proc%comm, te)
    call global_sum(proc%comm, tav)
    call global_sum(proc%comm, tpe)

    do iblk = 1, size(blocks)
      blocks(iblk)%state(itime,ivector)%tm  = tm
      blocks(iblk)%state(itime,ivector)%te  = te
      blocks(iblk)%state(itime,ivector)%tav = tav
      blocks(iblk)%state(itime,ivector)%tpe = tpe
      !call diag_state(iblk)%run(proc%blocks(iblk)%state(itime))
    end do

    call log_add_diag('tm' , tm(member_output) )
    call log_add_diag('te' , te(member_output) )
    call log_add_diag('tpe', tpe(member_output) )

  end subroutine diagnose

  subroutine space_operators(block, old_state, star_state, new_state, tend1, tend2, dt, pass)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout)    :: old_state
    type(state_type), intent(inout) :: star_state
    type(state_type), intent(inout) :: new_state
    type(tend_type ), intent(inout) :: tend1
    type(tend_type ), intent(in   ) :: tend2
    real(8), intent(in) :: dt
    integer, intent(in) :: pass

    type(mesh_type), pointer :: mesh
    integer i, j, k


    mesh => star_state%mesh

    call tend1%reset_flags()

    select case (pass)
    case (all_pass)
      call operators_prepare(block, star_state, dt, pass)
      if (hydrostatic) then
        call calc_dmfdlon_dmfdlat  (block, star_state, tend1, dt)
        call calc_dphs             (block, star_state, tend1, dt)
        call calc_wedphdlev_lev    (block, star_state, tend1, dt)
        call calc_wedudlev_wedvdlev(block, star_state, tend1, dt)
        call calc_dptfdlon_dptfdlat(block, star_state, tend1, dt)
        call calc_dptfdlev         (block, star_state, tend1, dt)
        call calc_qhu_qhv          (block, star_state, tend1, dt)
        call calc_dkedlon_dkedlat  (block, star_state, tend1, dt)
        call pgf_run               (block, star_state, tend1)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend1%du(:,i,j,k) =   tend1%qhv(:,i,j,k) - tend1%pgf_lon(:,i,j,k) - tend1%dkedlon(:,i,j,k) - tend1%wedudlev(:,i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dv(:,i,j,k) = - tend1%qhu(:,i,j,k) - tend1%pgf_lat(:,i,j,k) - tend1%dkedlat(:,i,j,k) - tend1%wedvdlev(:,i,j,k)

            end do
          end do

          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dpt(:,i,j,k) = - tend1%dptfdlon(:,i,j,k) - tend1%dptfdlat(:,i,j,k) - tend1%dptfdlev(:,i,j,k)
            end do
          end do
        end do

        tend1%update_u   = .true.
        tend1%update_v   = .true.
        tend1%update_phs = .true.
        tend1%update_pt  = .true.
      else if (nonhydrostatic) then
        call calc_dmfdlon_dmfdlat  (block, star_state, tend1, dt)
        call calc_dphs             (block, star_state, tend1, dt)
        call calc_wedphdlev_lev    (block, star_state, tend1, dt)
        call calc_dptfdlon_dptfdlat(block, star_state, tend1, dt)
        call calc_dptfdlev         (block, star_state, tend1, dt)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dpt(:,i,j,k) = - tend1%dptfdlon(:,i,j,k) - tend1%dptfdlat(:,i,j,k) - tend1%dptfdlev(:,i,j,k)
            end do
          end do
        end do

        tend1%update_phs = .true.
        tend1%update_pt  = .true.
        call update_state(block, tend1, old_state, new_state, dt)

        !call nh_solve(block, tend1, old_state, star_state, new_state, dt)

        call calc_qhu_qhv          (block, star_state, tend1, dt)
        call calc_dkedlon_dkedlat  (block, star_state, tend1, dt)
        call calc_wedudlev_wedvdlev(block, star_state, tend1, dt)
        call pgf_run               (block, new_state, tend1)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend1%du(:,i,j,k) =   tend1%qhv(:,i,j,k) - tend1%pgf_lon(:,i,j,k) - tend1%dkedlon(:,i,j,k) - tend1%wedudlev(:,i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dv(:,i,j,k) = - tend1%qhu(:,i,j,k) - tend1%pgf_lat(:,i,j,k) - tend1%dkedlat(:,i,j,k) - tend1%wedvdlev(:,i,j,k)
            end do
          end do
        end do

        tend1%update_u   = .true.
        tend1%update_v   = .true.

      else
        call calc_dmfdlon_dmfdlat(block, star_state, tend1, dt)
        call calc_qhu_qhv        (block, star_state, tend1, dt)
        call calc_dkedlon_dkedlat(block, star_state, tend1, dt)
        call pgf_run             (block, star_state, tend1    )

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend1%du(:,i,j,k) =   tend1%qhv(:,i,j,k) - tend1%pgf_lon(:,i,j,k) - tend1%dkedlon(:,i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dv(:,i,j,k) = - tend1%qhu(:,i,j,k) - tend1%pgf_lat(:,i,j,k) - tend1%dkedlat(:,i,j,k)
            end do
          end do

          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dgz(:,i,j,k) = - (tend1%dmfdlon(:,i,j,k) + tend1%dmfdlat(:,i,j,k)) * g
            end do
          end do
        end do

        tend1%update_u  = .true.
        tend1%update_v  = .true.
        tend1%update_gz = .true.
      end if
    case (forward_pass)
      call operators_prepare(block, star_state, dt, pass)
      if (hydrostatic) then
        call calc_dmfdlon_dmfdlat  (block, star_state, tend1, dt)
        call calc_dphs             (block, star_state, tend1, dt)
        call calc_wedphdlev_lev    (block, star_state, tend1, dt)
        call calc_wedudlev_wedvdlev(block, star_state, tend1, dt)
        call calc_dptfdlon_dptfdlat(block, star_state, tend1, dt)
        call calc_dptfdlev         (block, star_state, tend1, dt)
        call calc_qhu_qhv          (block, star_state, tend1, dt)
        call calc_dkedlon_dkedlat  (block, star_state, tend1, dt)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend1%du(:,i,j,k) =   tend1%qhv(:,i,j,k) - tend1%dkedlon(:,i,j,k) - tend1%wedudlev(:,i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dv(:,i,j,k) = - tend1%qhu(:,i,j,k) - tend1%dkedlat(:,i,j,k) - tend1%wedvdlev(:,i,j,k)
            end do
          end do

          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dpt(:,i,j,k) = - tend1%dptfdlon(:,i,j,k) - tend1%dptfdlat(:,i,j,k) - tend1%dptfdlev(:,i,j,k)
            end do
          end do
        end do

        tend1%update_phs = .true.
        tend1%update_pt  = .true.
      else if (nonhydrostatic) then

      else
        call calc_dmfdlon_dmfdlat(block, star_state, tend1, dt)
        call calc_qhu_qhv        (block, star_state, tend1, dt)
        call calc_dkedlon_dkedlat(block, star_state, tend1, dt)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend1%du(:,i,j,k) =   tend1%qhv(:,i,j,k) - tend1%dkedlon(:,i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dv(:,i,j,k) = - tend1%qhu(:,i,j,k) - tend1%dkedlat(:,i,j,k)
            end do
          end do

          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dgz(:,i,j,k) = - (tend1%dmfdlon(:,i,j,k) + tend1%dmfdlat(:,i,j,k)) * g
            end do
          end do
        end do

        tend1%update_gz = .true.
      end if
    case (backward_pass)
      call operators_prepare(block, new_state, dt, pass)
      if (hydrostatic) then
        call pgf_run(block, new_state, tend1)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend1%du(:,i,j,k) = tend2%du(:,i,j,k) - tend1%pgf_lon(:,i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dv(:,i,j,k) = tend2%dv(:,i,j,k) - tend1%pgf_lat(:,i,j,k)
            end do
          end do
        end do

        tend1%update_u   = .true.
        tend1%update_v   = .true.
      else if (nonhydrostatic) then

      else
        call pgf_run(block, new_state, tend1)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend1%du(:,i,j,k) = tend2%du(:,i,j,k) - tend1%pgf_lon(:,i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dv(:,i,j,k) = tend2%dv(:,i,j,k) - tend1%pgf_lat(:,i,j,k)
            end do
          end do
        end do

        tend1%update_u  = .true.
        tend1%update_v  = .true.
      end if
    end select

    ! call debug_check_space_operators(block, state, tend)

  end subroutine space_operators

  subroutine time_integrate(dt, blocks)

    real(8), intent(in) :: dt
    type(block_type), intent(inout) :: blocks(:)

    integer iblk

    do iblk = 1, size(blocks)

      call time_integrator(operators, blocks(iblk), old, new, dt)
      !call test_forcing_run(blocks(iblk), dt, blocks(iblk)%static(ivector), blocks(iblk)%state(new,ivector))

    end do
    call damp_run(dt, new, blocks)

  end subroutine time_integrate

end module gmcore_mod
