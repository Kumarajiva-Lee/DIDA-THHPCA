program dida

  use mpi
  use coupler_config
  use coupler_mod
  use flogger

  use gmcore_mod

  use da_mod
  use da_parallel
  use da_coupler

  implicit none

  integer :: i,ierr
  character(256) namelist_path_dida , namelist_path_ATM

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, globalsize, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myglobalrank, ierr)
  call GPTLinitialize()
  call GPTLstart("total")

#ifdef RSL_LOG
  call rsl_error_dup1(myglobalrank)
#endif

  call log_print_master("Start DIDA Program!")
  call get_command_argument(1, namelist_path_dida)
  call log_print_master("read namelist:"//trim(namelist_path_dida))
  call init_coupler_config(namelist_path_dida)
  call log_print_master("init_coupler_config done!")
  call init_coupler()
  call log_print_master("init_coupler done!")

#if defined ATM_COUPLING
  call GPTLstart("atm total time")
  do i = 0, atm_ensemble_group - 1
    if (atm_gcolor(myglobalrank) .eq. i) then
      call log_print_root("ATM MODEL BEGIN!")
      call atm_init(namelist_path_dida , atm_groupmember%atm_subcomm,group_id = i)
      call atm_run(atm_groupmember%atm_subcomm,group_id = i)
      call atm_final()
      call log_print_root("ATM MODEL END!")
    end if
  end do
  call GPTLstop("atm total time")
#endif

#if defined DA_COUPLING
  call GPTLstart("da total time")
  if (atm_gcolor(cinfo%globalrank) .eq. DA_Color) then
    call log_print_root("DA MODEL BEGIN!")
    call init_da_parallel(namelist_path_dida, atm_groupmember%atm_subcomm)
    ! print *,'init_da_parallel sucess'
    call init_da_coupler(atm_groupmember%atm_subcomm)
    call dafatm_coupling(atm_groupmember%atm_subcomm)
    call log_print_root("DA MODEL END!")
   ! call clean_da_coupler()
  end if
  call GPTLstop("da total time")
#endif

  if (atm_gcolor(cinfo%globalrank) .eq. DEFAULT_Color) then
    call log_print_root("SOME NODES ARE FREE IN DEFAULE COMM!")
    call Default_driver(atm_groupmember%atm_subcomm)
    call clean_default_coupler()
  end if

  call GPTLstop("total")
  !call GPTLpr(myglobalrank)
  if(myglobalrank == 0) call GPTLpr(myglobalrank)
  !call GPTLpr_summary()
  call GPTLfinalize()

  call log_print_master("DIDA End!")
  call MPI_Finalize(ierr)

END program
