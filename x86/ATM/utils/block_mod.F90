module block_mod

  use mpi
  use flogger
  use namelist_mod
  use mesh_mod
  use state_mod
  use static_mod
  use tend_mod
  use filter_types_mod
  use halo_mod
  use allocator_mod

  implicit none

  private

  public block_type
  public global_mesh
  public mesh_type
  public state_type
  public static_type
  public tend_type

  !      |         |         |
  !  ____|_________|_________|_____
  !      |                   |
  !      |                   |
  !      |                   |
  !      |                   |
  !      |       BLOCK       |
  !      |                   |
  !      |                   |
  !      |                   |
  !  ____|___________________|_____
  !      |                   |
  !      |                   |

  type block_type
    integer id
    type(mesh_type) mesh
    type(state_type), allocatable :: state(:,:)
    type(static_type),allocatable :: static(:)
    type(tend_type), allocatable :: tend(:,:)
    type(halo_type), allocatable :: halo(:)
    type(halo_type), allocatable :: halo_send(:)
    type(halo_type), allocatable :: halo_recv(:)
    type(filter_type) filter
    integer member_value
    integer ensemble_position
  contains
    procedure :: init_stage_1 => block_init_stage_1
    procedure :: init_stage_2 => block_init_stage_2
    procedure :: clear        => block_clear
    final :: block_final
  end type block_type

contains

  subroutine block_init_stage_1(this, id, lon_halo_width, lat_halo_width, lon_ibeg, lon_iend, lat_ibeg, lat_iend)

    class(block_type), intent(inout) :: this
    integer, intent(in) :: id
    integer, intent(in) :: lon_halo_width
    integer, intent(in) :: lat_halo_width
    integer, intent(in) :: lon_ibeg
    integer, intent(in) :: lon_iend
    integer, intent(in) :: lat_ibeg
    integer, intent(in) :: lat_iend

    this%id = id

    call this%mesh%init_from_parent(global_mesh, this%id, lon_halo_width, lat_halo_width, lon_ibeg, lon_iend, lat_ibeg, lat_iend)
    call this%filter%init(this%mesh)

  end subroutine block_init_stage_1

  subroutine block_init_stage_2(this, lon_halo_width)

    class(block_type), intent(inout) :: this
    integer, intent(in) :: lon_halo_width


    integer i,iv

    call this%mesh%reinit(lon_halo_width)
    

    if (.not. allocated(this%state)) then
      select case (trim(time_scheme))
      case ('euler')
        allocate(this%state(2,vector_num))
        allocate(this%tend (2,vector_num))
      case ('pc2', 'wrfrk3')
        if (use_create_ensemble) then
          allocate(this%state(0:3,vector_num))
          allocate(this%tend (0:3,vector_num))
          this%ensemble_position = 0
        else
          allocate(this%state(3,vector_num))
          allocate(this%tend (3,vector_num))
        end if
      case ('rk3', 'ssprk3')
        allocate(this%state(4,vector_num))
        allocate(this%tend (4,vector_num))
      case ('rk4')
        allocate(this%state(5,vector_num))
        allocate(this%tend (5,vector_num))
      case ('N/A')
        allocate(this%state(1,vector_num))
      case default
        if (this%id == 0) call log_error('Unknown time scheme ' // trim(time_scheme))
      end select
      allocate(this%static(vector_num))
      do iv = 1 , vector_num
        do i = lbound(this%state, 1), ubound(this%state, 1)
          call this%state(i,iv)%init(this%mesh)
        end do
        do i = lbound(this%tend, 1), ubound(this%tend, 1)
          call this%tend(i,iv)%init(this%mesh)
        end do
        call this%static(iv)%init(this%mesh)
      end do
    end if

  end subroutine block_init_stage_2

  subroutine block_clear(this)

    class(block_type), intent(inout) :: this

    integer i,iv
    
    do iv = 1 , vector_num
      do i = 1, size(this%state,1)
        call this%state(i,iv)%clear()
      end do
      do i = 1, size(this%tend,1)
        call this%tend(i,iv)%clear()
      end do
    end do

    do i = 1, size(this%halo)
      call this%halo(i)%clear()
    end do

    do i = 1, size(this%halo_send)
      call this%halo_send(i)%clear()
    end do

    do i = 1, size(this%halo_recv)
      call this%halo_recv(i)%clear()
    end do




    if (allocated(this%halo )) deallocate(this%halo )
    if (allocated(this%state)) deallocate(this%state)
    if (allocated(this%tend )) deallocate(this%tend )
    if (allocated(this%halo )) deallocate(this%halo )

  end subroutine block_clear

  subroutine block_final(this)

    type(block_type), intent(inout) :: this

    call this%clear()

  end subroutine block_final

end module block_mod

