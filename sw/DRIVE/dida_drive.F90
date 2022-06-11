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
  character(256) hname
  integer hlen

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, globalsize, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myglobalrank, ierr)
  call MPI_GET_PROCESSOR_NAME(hname, hlen, ierr)
  call GPTLinitialize()
  call GPTLstart("total")

#ifdef RSL_LOG
  call rsl_error_dup1(myglobalrank)
#endif

  call log_print_master("Start DIDA Program!")
  !call log_notice("HOST_NAME: " // trim(hname))
  call get_command_argument(1, namelist_path_dida)
  call log_print_master("read namelist:"//trim(namelist_path_dida))
  call init_coupler_config(namelist_path_dida)
  call log_print_master("init_coupler_config done!")
  call init_coupler()
  call log_print_master("init_coupler done!")

#if defined ATM_COUPLING
  do i = 0, atm_ensemble_group - 1
    if (atm_gcolor(myglobalrank) .eq. i) then
      call GPTLstart("atm total time")
      call log_print_root("ATM MODEL BEGIN!")
      call atm_init(namelist_path_dida , atm_groupmember%atm_subcomm , group_id = i)
      call atm_run(atm_groupmember%atm_subcomm,group_id = i)
      call atm_final()
      call log_print_root("ATM MODEL END!")
      call GPTLstop("atm total time")
    end if
  end do
#endif

#if defined DA_COUPLING
  if (atm_gcolor(cinfo%globalrank) .eq. DA_Color) then
    call GPTLstart("da total time")
    call log_print_root("DA MODEL BEGIN!")
    call init_da_parallel(namelist_path_dida, atm_groupmember%atm_subcomm)
    !! print *,'init_da_parallel sucess'
    call init_da_coupler(atm_groupmember%atm_subcomm)
    call dafatm_coupling(atm_groupmember%atm_subcomm)
    call log_print_root("DA MODEL END!")
    call GPTLstop("da total time")
   ! call clean_da_coupler()
  end if
#endif

  !if (atm_gcolor(cinfo%globalrank) .eq. CNN_Color) then
  !  call GPTLstart("cnn total time")
  !  call log_print_root("CNN MODEL BEGIN!")
  !  call cnn_init(atm_groupmember%atm_subcomm)
  !  call cnn_run(atm_groupmember%atm_subcomm)
  !  call cnn_final()
  !  call log_print_root("CNN MODEL END!")
  !  call GPTLstop("cnn total time")
  !end if

  if (atm_gcolor(cinfo%globalrank) .eq. DEFAULT_Color) then
    call log_print_root("SOME NODES ARE FREE IN DEFAULE COMM!")
    call Default_driver(atm_groupmember%atm_subcomm)
    call clean_default_coupler()
  end if

  call GPTLstop("total")
  call GPTLpr(myglobalrank)
  !call GPTLpr_summary()
  call GPTLfinalize()

  call log_print_master("DIDA End!")
  call MPI_Finalize(ierr)

END program
