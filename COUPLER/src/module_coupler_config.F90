module coupler_config
  use mpi
  implicit none

  type coupler_info
    integer :: ATM_DA_GROUP            = MPI_COMM_NULL
    integer :: globalrank            = MPI_PROC_NULL !MPI rank in MPI_COMM_WORLD
    integer :: globalsize                              !MPI_COMM_WORLD size
    integer :: total_np                                !total  process num
    integer :: atm_np                                  !total  process num
    integer :: da_np                                   !total  process num
    integer :: ncomps                                  !MCT component model  num
    integer :: atm_comps                                  !MCT component model  num
    integer :: da_comps                                  !MCT component model  num
  end type coupler_info

  type super_node
    integer :: snode_comm    = MPI_COMM_NULL      !super node comm
    integer :: snode_id                !super node id
    integer :: snode_rank            = MPI_PROC_NULL     !super node id
    integer :: snode_size              !super node id
    integer s_num_lon
    integer s_num_lat
    integer s_lon_ibeg
    integer s_lon_iend
    integer s_lat_ibeg
    integer s_lat_iend
    integer :: snodeinfo(4)
  end type super_node

  type atm_group_member
      integer :: atm_subgid          !atm group member num
      integer :: atm_subcomm         !atm group member comm
      integer :: atm_snodesubcomm    !atm group member comm in on super node
      integer :: snode_id            !super node tag
      integer :: atm_subrank            = MPI_PROC_NULL     !mpi_rank in atm group member comm 
      integer :: atm_subsize                                !mpi_size in atm group member comm
      integer :: atm_snodesubrank            = MPI_PROC_NULL     !mpi_rank in super node
      integer :: atm_snodesubsize            = MPI_PROC_NULL     !mpi_size in one super node
  end type atm_group_member

  !coupler constant
  integer :: DEFAULT_Id = 0
  integer :: CNN_Id = 3
  integer :: DA_Id = 2
  integer :: ATM_Id = 1
  integer :: DEFAULT_Color = 9998
  integer :: ATM_Color = 9999
  integer :: DA_Color = 10000
  integer :: CNN_Color = 10001

  integer, dimension(:), allocatable  :: atm_subid
  type(super_node) :: snode
  type(coupler_info) :: cinfo
  type(atm_group_member)   :: atm_groupmember !atm member number
  integer :: ncomps                                  !MCT component model  num
  integer :: myglobalrank, globalsize

  !namelist total
  character(30)  :: case_name
  character(256) :: redis_address
  character(30) :: atm_mpas_sceneid
  integer :: num_lon
  integer :: num_lat
  integer :: num_lev

  character(30)   :: vert_coord_scheme    = 'hybrid'
  character(30)   :: vert_coord_template  = 'N/A'

  character(30)  :: atm_mpas2atm_time
  integer :: atm_phase_diff
  integer :: da_endnn_in_seconds
  integer :: da_start_time(5)
  integer :: end_time(5)
  integer :: da_in_seconds
  integer :: cnn_ensemble
  integer :: atm_ensemble_total
  integer :: atm_ensemble_group
  integer :: da_ensemble
  character(20) :: da_var_name
  integer :: da_asynchronous
  integer :: da_mode
  character(200) :: cnn_model_path
  character(200) :: cnn_mean_path
  character(200) :: cnn_std_path
  integer :: cnn_group
  integer :: cnn_halo
  integer :: atm_group
  integer :: atm_group_num
  integer :: atm_stride
  integer :: atm_root
  integer :: da_group
  integer :: da_group_num
  integer :: da_stride
  integer :: da_root

  !statements that specify the namelists
  namelist /namelist_total/ &
    case_name             , &
    redis_address         , &
    atm_mpas_sceneid      , &
    num_lon               , &
    num_lat               , &
    num_lev               , &

    vert_coord_scheme     , &
    vert_coord_template   , &

    atm_mpas2atm_time     , &
    atm_phase_diff        , &
    da_endnn_in_seconds   , &
    da_start_time         , &
    end_time              , &
    da_in_seconds         , &
    cnn_ensemble          , &
    atm_ensemble_total    , &
    atm_ensemble_group    , &
    da_ensemble           , &
    da_var_name           , &
    da_asynchronous       , &
    da_mode               , &
    cnn_model_path        , &
    cnn_mean_path         , &
    cnn_std_path          , &
    cnn_group             , &
    cnn_halo              , &
    atm_group             , &
    atm_group_num         , &
    atm_stride            , &
    atm_root              , &
    da_group              , &
    da_group_num          , &
    da_stride             , &
    da_root               

end module coupler_config 
