module coupler_log

  use mpi
  use coupler_config

  implicit none

  integer :: fid=6
  character*8, mytime, mydate
  integer  :: myrank
  integer :: ierr
  private

#ifdef LOG_MINFO
  public log_init_file
  public log_close_file
#endif
  public log_notice_root
  public log_print_namelist
  public log_print_couplerinfo
  public log_warning_root
  public log_error_root
  public log_notice
  public log_warning
  public log_error
  public log_print
  public log_print_root

contains

  subroutine log_print_namelist()
    use string
    implicit none

    call log_print_root("=========== total namelists ========")
    call log_print_root("case_name          :" // case_name)
    call log_print_root("redis_address      :" // (redis_address))
    call log_print_root("atm_mpas_sceneid   :" // (atm_mpas_sceneid   ))
    call log_print_root("num_lon            :" // to_str(num_lon            ))
    call log_print_root("num_lat            :" // to_str(num_lat            ))
    call log_print_root("num_lev            :" // to_str(num_lev            ))
    call log_print_root("atm_mpas2atm_time  :" // (atm_mpas2atm_time  ))
    call log_print_root("atm_phase_diff     :" // to_str(atm_phase_diff     ))
    call log_print_root("da_endnn_in_seconds:" // to_str(da_endnn_in_seconds))
    call log_print_root("da_start_time      :" // to_str(da_start_time      ))
    call log_print_root("end_time           :" // to_str(end_time           ))
    call log_print_root("da_in_seconds      :" // to_str(da_in_seconds      ))
    call log_print_root("atm_ensemble_total :" // to_str(atm_ensemble_total ))
    call log_print_root("atm_ensemble_group :" // to_str(atm_ensemble_group ))
    call log_print_root("da_ensemble        :" // to_str(da_ensemble        ))
    call log_print_root("da_var_name        :" // (da_var_name        ))
    call log_print_root("da_asynchronous    :" // to_str(da_asynchronous    ))
    call log_print_root("da_mode            :" // to_str(da_mode            ))
    call log_print_root("atm_group          :" // to_str(atm_group          ))
    call log_print_root("atm_group_num      :" // to_str(atm_group_num      ))
    call log_print_root("atm_stride         :" // to_str(atm_stride         ))
    call log_print_root("atm_root           :" // to_str(atm_root           ))
    call log_print_root("da_group           :" // to_str(da_group           ))
    call log_print_root("da_group_num       :" // to_str(da_group_num       ))
    call log_print_root("da_stride          :" // to_str(da_stride          ))
    call log_print_root("da_root            :" // to_str(da_root            ))
  end subroutine

  subroutine log_print_couplerinfo()
    use string
    implicit none
    call log_print_root("=========== COUPLER INFO ========")
    call log_print_root("Total PROCS        :" // to_str(cinfo%globalsize))
    call log_print_root("ATM TOTAL PROCS    :" // to_str(cinfo%atm_np))
    call log_print_root("DA TOTAL PROCS     :" // to_str(cinfo%da_np))
    call log_print_root("MCT ncomps         :" // to_str(cinfo%ncomps))
    call log_print_root("ATM comps          :" // to_str(cinfo%atm_comps))
    call log_print_root("DA comps           :" // to_str(cinfo%da_comps))

  end subroutine

#ifdef LOG_MINFO
  subroutine log_init_file()
    character(128) :: filename

    call time(mytime)
    call date_and_time(mydate)

    if(atm_groupmember%atm_subgid .eq. (DA_Color + 1)) then
        filename = 'da_log_' // trim(mydate) //"-"//trim(mytime) // '.log'
        fid=60
        open(UNIT = fid , FILE = trim(filename))
    else
        filename = 'atm_log_' // trim(mydate) //"-"//trim(mytime) // '.log'
        fid=70
        open(UNIT = fid , FILE = trim(filename))
    endif

  end subroutine

  subroutine log_close_file()
    character(128) :: filename

    call time(mytime)
    call date_and_time(mydate)

    if(atm_groupmember%atm_subgid .eq. (DA_Color + 1)) then
        close(UNIT = fid)
    else
        close(UNIT = fid)
    endif

  end subroutine
#endif

  subroutine log_notice(message, file, line)

    character(*), intent(in) :: message
    character(*), intent(in), optional :: file
    integer, intent(in), optional :: line

    call time(mytime)
    call date_and_time(mydate)

    if (present(file) .and. present(line)) then
#ifdef LOG_MINFO
        write(fid, '("[", A7, "]","[", A17, "]",  A, ":", I0, A)') 'Notice', &
              mydate // "/"//mytime, trim(file), line, "[M: COUPLER"// "] " // trim(message)
#else
        write(6, '("[", A, "]: ", A, ":", I0, ": ", A)') 'Notice', trim(file), line, trim(message)
#endif
    else
#ifdef LOG_MINFO
        write(fid, '("[", A7, "]","[", A17, "]", A)') 'Notice', &
              mydate // "/"//mytime, "[M: COUPLER"// "] " // trim(message)
#else
        write(6, '("[", A, "]: ", A)') 'Notice', trim(message)
#endif
    end if

  end subroutine log_notice

  subroutine log_warning(message, file, line)

    character(*), intent(in) :: message
    character(*), intent(in), optional :: file
    integer, intent(in), optional :: line

    call time(mytime)
    call date_and_time(mydate)

    if (present(file) .and. present(line)) then
#ifdef LOG_MINFO
        write(fid, '("[", A7, "]","[", A17, "]",  A, ":", I0, A)') 'Warning', &
              mydate // "/"//mytime, trim(file), line, "[M: COUPLER"// "] " // trim(message)
#else
        write(6, '("[", A, "]: ", A, ":", I0, ": ", A)') 'Warning', trim(file), line, trim(message)
#endif
    else
#ifdef LOG_MINFO
        write(fid, '("[", A7, "]","[", A17, "]", A)') 'Warning', &
              mydate // "/"//mytime, "[M: COUPLER"// "] " // trim(message)
#else
        write(6, '("[", A, "]: ", A)') 'Warning', trim(message)
#endif
    end if

  end subroutine log_warning

  subroutine log_error(message, file, line)

    character(*), intent(in) :: message
    character(*), intent(in), optional :: file
    integer, intent(in), optional :: line

    call time(mytime)
    call date_and_time(mydate)

    if (present(file) .and. present(line)) then
#ifdef LOG_MINFO
          write(fid, '("[", A7, "]: ","[", A17, "]: ",  A, ":", I0, ": ", A)') 'Error', &
              mydate // "/"//mytime, trim(file), line, "[M: COUPLER"// "] " // trim(message)
#else
        write(6, '("[", A, "]: ", A, ":", I0, ": ", A)') 'Error', trim(file), line, trim(message)
#endif
    else
#ifdef LOG_MINFO
        write(fid, '("[", A7, "]","[", A17, "]", A)') 'Error', &
              mydate // "/"//mytime, "[M: COUPLER"// "] " // trim(message)
#else
        write(6, '("[", A, "]: ", A)') 'Error', trim(message)
#endif
    end if
    stop 1

  end subroutine log_error

  subroutine log_print(message)

    character(*), intent(in) :: message

    call time(mytime)
    call date_and_time(mydate)

#ifdef LOG_MINFO
          write(fid, '("[", A7, "]","[", A17, "]", A)') '===> ', &
              mydate // "/"//mytime, "[M: COUPLER"// "] " // trim(message)
#else
      write(6, '(A)') '==> ' // trim(message)
#endif

  end subroutine log_print

  subroutine log_print_root(message)

    character(*), intent(in) :: message

    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

    if(myrank == 0) then
      call log_print(message)
    end if
  end subroutine

  subroutine log_print_master(message)

    character(*), intent(in) :: message
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

    if(myrank == 0) then
      call log_print(message)
    end if
  end subroutine

  subroutine log_notice_root(message, file, line)

    character(*), intent(in) :: message
    character(*), intent(in), optional :: file
    integer, intent(in), optional :: line

    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

    if(myrank == 0) then
      call log_notice(message, file, line)
    end if
  end subroutine

  subroutine log_warning_root(message, file, line)

    character(*), intent(in) :: message
    character(*), intent(in), optional :: file
    integer, intent(in), optional :: line

    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

    if(myrank == 0) then
      call log_warning(message, file, line)
    end if
  end subroutine

  subroutine log_error_root(message, file, line)

    character(*), intent(in) :: message
    character(*), intent(in), optional :: file
    integer, intent(in), optional :: line

    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

    if(myrank == 0) then
      call log_error(message, file, line)
    end if
  end subroutine

end module coupler_log
