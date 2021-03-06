module time_schemes_mod

  use flogger
  use const_mod
  use namelist_mod
  use tend_mod
  use block_mod
  use operators_mod
  use parallel_mod
  use parallel_types_mod
  use filter_mod
  use pa_mod
  use time_mod
  use nh_mod
  use member_mod

  implicit none

  private

  public time_scheme_init
  public time_integrator
  public update_state
  public space_operators_interface

  interface
    subroutine space_operators_interface(block, old_state, star_state, new_state, tend1, tend2, dt, pass)
      import block_type, state_type, tend_type
      type(block_type), intent(inout) :: block
      type(state_type), intent(inout) :: old_state
      type(state_type), intent(inout) :: star_state
      type(state_type), intent(inout) :: new_state
      type(tend_type ), intent(inout) :: tend1
      type(tend_type ), intent(in   ) :: tend2
      real(8), intent(in) :: dt
      integer, intent(in) :: pass
    end subroutine space_operators_interface

    subroutine step_interface(space_operators, block, old_state, star_state, new_state, tend1, tend2, dt)
      import space_operators_interface, block_type, state_type, tend_type
      procedure(space_operators_interface), intent(in), pointer :: space_operators
      type(block_type), intent(inout) :: block
      type(state_type), intent(inout) :: old_state
      type(state_type), intent(inout) :: star_state
      type(state_type), intent(inout) :: new_state
      type(tend_type ), intent(inout) :: tend1
      type(tend_type ), intent(inout) :: tend2
      real(8), intent(in) :: dt
    end subroutine step_interface

    subroutine time_integrator_interface(space_operators, block, old, new, dt)
      import block_type, tend_type, state_type, space_operators_interface
      procedure(space_operators_interface), intent(in), pointer :: space_operators
      type(block_type), intent(inout) :: block
      integer, intent(in) :: old
      integer, intent(in) :: new
      real(8), intent(in) :: dt
    end subroutine time_integrator_interface
  end interface

  procedure(step_interface), pointer :: step
  procedure(time_integrator_interface), pointer :: time_integrator

contains

  subroutine time_scheme_init()

    select case (time_scheme)
    case ('euler')
      time_integrator => euler
    case ('pc2')
      time_integrator => predict_correct
    case ('wrfrk3')
      time_integrator => wrf_runge_kutta_3rd
    case default
      time_integrator => predict_correct
    end select

    step => step_forward_backward

  end subroutine time_scheme_init

  subroutine step_all(space_operators, block, old_state, star_state, new_state, tend1, tend2, dt)

    procedure(space_operators_interface), intent(in), pointer :: space_operators
    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: old_state
    type(state_type), intent(inout) :: star_state
    type(state_type), intent(inout) :: new_state
    type(tend_type ), intent(inout) :: tend1
    type(tend_type ), intent(inout) :: tend2
    real(8), intent(in) :: dt
    call space_operators(block, old_state, star_state, new_state, tend1, tend2, dt, all_pass)
    call update_state(block, tend1, old_state, new_state, dt)

  end subroutine step_all

  subroutine step_forward_backward(space_operators, block, old_state, star_state, new_state, tend1, tend2, dt)

    procedure(space_operators_interface), intent(in), pointer :: space_operators
    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: old_state
    type(state_type), intent(inout) :: star_state
    type(state_type), intent(inout) :: new_state
    type(tend_type ), intent(inout) :: tend1
    type(tend_type ), intent(inout) :: tend2
    real(8), intent(in) :: dt
    
    call space_operators(block, old_state, star_state, new_state, tend1, tend2, dt, forward_pass)
    call update_state(block, tend1, old_state, new_state, dt)
    call space_operators(block, old_state, star_state, new_state, tend2, tend1, dt, backward_pass)
    call update_state(block, tend2, old_state, new_state, dt)

  end subroutine step_forward_backward

  subroutine update_state(block, tend, old_state, new_state, dt)

    type(block_type), intent(inout) :: block
    type(tend_type), intent(in) :: tend
    type(state_type), intent(inout) :: old_state
    type(state_type), intent(inout) :: new_state
    real(8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, k

    associate (mesh => block%mesh)

    if (baroclinic) then
      if (tend%update_phs) then
        call old_state%async(async_phs)%wait()

        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            new_state%phs(:,i,j) = old_state%phs(:,i,j) + dt * tend%dphs(:,i,j)
          end do
        end do

        call fill_halo_member(block, new_state%phs, full_lon=.true., full_lat=.true., async=new_state%async(async_phs))
        call new_state%async(async_phs)%wait()
        call filter_on_cell(block, new_state%phs, new_state%phs_f)
        call fill_halo_member(block, new_state%phs_f, full_lon=.true., full_lat=.true., async=new_state%async(async_phs_f))
        call new_state%async(async_phs_f)%wait()

        call diag_ph(block, new_state)
        call diag_m (block, new_state)
        if (nonhydrostatic) then
          !call diag_m_lev(block, new_state)
        end if
      else if (tend%copy_phs) then
        call new_state%async(async_phs)%wait()
        call new_state%async(async_phs_f)%wait()
        call new_state%async(async_ph)%wait()
        call new_state%async(async_ph_lev)%wait()
        new_state%phs    = old_state%phs
        new_state%phs_f  = old_state%phs_f
        new_state%ph_lev = old_state%ph_lev
        new_state%ph     = old_state%ph
        new_state%m      = old_state%m
      end if

      if (tend%update_pt) then
        if (.not. tend%update_phs .and. .not. tend%copy_phs .and. is_root_proc()) call log_error('Mass is not updated or copied!')

        call old_state%async(async_pt)%wait()

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              new_state%pt(:,i,j,k) = (old_state%pt(:,i,j,k) * old_state%m(:,i,j,k) + dt * tend%dpt(:,i,j,k)) / new_state%m(:,i,j,k)
            end do
          end do
        end do

        call fill_halo_member(block, new_state%pt, full_lon=.true., full_lat=.true., full_lev=.true., async=new_state%async(async_pt))
        call new_state%async(async_pt)%wait()
        ! NOTE: When we use FFSL to transport pt, it may be not needed to filter pt.
        call filter_on_cell(block, new_state%pt, new_state%pt_f)
        call fill_halo_member(block, new_state%pt_f, full_lon=.true., full_lat=.true., full_lev=.true., async=new_state%async(async_pt_f))
        call new_state%async(async_pt_f)%wait()

      else if (tend%copy_pt) then
        call new_state%async(async_pt)%wait()
        call new_state%async(async_pt_f)%wait()
        new_state%pt = old_state%pt
        new_state%pt_f = old_state%pt_f
      end if
    else
      if (tend%update_gz) then

        call old_state%async(async_gz)%wait()

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              new_state%gz(:,i,j,k) = old_state%gz(:,i,j,k) + dt * tend%dgz(:,i,j,k)
            end do
          end do
        end do

        call fill_halo_member(block, new_state%gz, full_lon=.true., full_lat=.true., async=new_state%async(async_gz))
        call filter_on_cell(block, new_state%gz, new_state%gz_f)
        call fill_halo_member(block, new_state%gz_f, full_lon=.true., full_lat=.true., async=new_state%async(async_gz_f))
        call new_state%async(async_gz_f)%wait()

      else if (tend%copy_gz) then
        call new_state%async(async_gz)%wait()
        call new_state%async(async_gz_f)%wait()
        new_state%gz = old_state%gz
        new_state%gz_f = old_state%gz_f
      end if
    end if

    if (tend%update_u .and. tend%update_v) then
      call old_state%async(async_u)%wait()
      call old_state%async(async_v)%wait()
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            new_state%u(:,i,j,k) = old_state%u(:,i,j,k) + dt * tend%du(:,i,j,k)
          end do
        end do
      end do
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            new_state%v(:,i,j,k) = old_state%v(:,i,j,k) + dt * tend%dv(:,i,j,k)
          end do
        end do
      end do   

      call fill_halo_member(block, new_state%u, full_lon=.false., full_lat=.true., full_lev=.true., async=new_state%async(async_u))
      call fill_halo_member(block, new_state%v, full_lon=.true., full_lat=.false., full_lev=.true., async=new_state%async(async_v))
      call new_state%async(async_u)%wait()
      call new_state%async(async_v)%wait()
      if (mesh%has_south_pole()) then
        j = mesh%half_lat_ibeg
        new_state%v(:,:,j,:) = 0.4_r8 * new_state%v(:,:,j,:) + 0.6_r8 * new_state%v(:,:,j+1,:)
      end if
      if (mesh%has_north_pole()) then
        j = mesh%half_lat_iend
        new_state%v(:,:,j,:) = 0.4_r8 * new_state%v(:,:,j,:) + 0.6_r8 * new_state%v(:,:,j-1,:)
      end if
      
      call filter_on_lon_edge(block, new_state%u, new_state%u_f)
      call filter_on_lat_edge(block, new_state%v, new_state%v_f)
      call fill_halo_member(block, new_state%u_f, full_lon=.false., full_lat=.true., full_lev=.true., async=new_state%async(async_u_f))
      call fill_halo_member(block, new_state%v_f, full_lon=.true., full_lat=.false., full_lev=.true., async=new_state%async(async_v_f))

      call new_state%async(async_u_f)%wait()
      call new_state%async(async_v_f)%wait()
    end if

    if (time_is_alerted('filter_reset')) then
      if (tend%update_phs) new_state%phs = new_state%phs_f
      if (tend%update_pt ) new_state%pt  = new_state%pt_f
      if (tend%update_u  ) new_state%u   = new_state%u_f
      if (tend%update_v  ) new_state%v   = new_state%v_f
    end if
    end associate

  end subroutine update_state

  subroutine predict_correct(space_operators, block, old, new, dt)

    procedure(space_operators_interface), intent(in), pointer :: space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(8), intent(in) :: dt
    
    associate (state => block%state, tend => block%tend)
      call step(space_operators, block, state(old,ivector), state(old,ivector), state(new,ivector), tend(old,ivector), tend(new,ivector), dt / 2.0_r8)
      call step(space_operators, block, state(old,ivector), state(new,ivector), state(3  ,ivector), tend(old,ivector), tend(new,ivector), dt / 2.0_r8)
      call step(space_operators, block, state(old,ivector), state(3  ,ivector), state(new,ivector), tend(old,ivector), tend(new,ivector), dt         )
    end associate

  end subroutine predict_correct

  subroutine wrf_runge_kutta_3rd(space_operators, block, old, new, dt)

    procedure(space_operators_interface), intent(in), pointer :: space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(8), intent(in) :: dt

    associate (state => block%state, tend => block%tend)
      call step(space_operators, block, state(old,ivector), state(old,ivector), state(new,ivector), tend(old,ivector), tend(new,ivector), dt / 3.0_r8)
      call step(space_operators, block, state(old,ivector), state(new,ivector), state(3  ,ivector), tend(old,ivector), tend(new,ivector), dt / 2.0_r8)
      call step(space_operators, block, state(old,ivector), state(3  ,ivector), state(new,ivector), tend(old,ivector), tend(new,ivector), dt         )
    end associate
 
  end subroutine wrf_runge_kutta_3rd

  
  subroutine euler(space_operators, block, old, new, dt)

    procedure(space_operators_interface), intent(in), pointer :: space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(8), intent(in) :: dt

    associate (state => block%state, tend => block%tend)
      call step(space_operators, block, state(old,ivector), state(old,ivector), state(new,ivector), tend(old,ivector), tend(new,ivector), dt)
    end associate

  end subroutine euler

end module time_schemes_mod
