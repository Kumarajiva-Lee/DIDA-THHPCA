module log_mod

  use hash_table_mod
  use string
  use face
  use mpi
#ifdef INTEL
  use ifport
#endif
#ifdef LOG_MINFO
  use coupler_config
#endif

  implicit none

  integer :: fid=6
  character*8, mytime, mydate
  integer  :: myrank
  integer :: ierr
  private

  public log_init
#ifdef LOG_MINFO
  public log_init_file
  public log_close_file
#endif
  public log_add_diag
  public log_notice_root
  public log_warning_root
  public log_error_root
  public log_notice
  public log_warning
  public log_error
  public log_print_diag
  public log_print
  public log_print_root
  public log_print_master

  type(hash_table_type) diags

contains

  subroutine log_init()

    diags = hash_table()

  end subroutine log_init

#ifdef LOG_MINFO
  subroutine log_init_file()
    character(128) :: filename

    !call time(mytime)
    mytime = ' '
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

    !call time(mytime)
    mytime = ' '
    call date_and_time(mydate)

    if(atm_groupmember%atm_subgid .eq. (DA_Color + 1)) then
        close(UNIT = fid)
    else
        close(UNIT = fid)
    endif

  end subroutine
#endif

  subroutine log_add_diag(name, value)

    character(*), intent(in) :: name
    class(*), intent(in) :: value

    select type (value)
    type is (integer)
    type is (real(4))
    type is (real(8))
    class default
      call log_error('Unsupported diagnostic value type!')
    end select

    call diags%insert(name, value)

  end subroutine log_add_diag

  subroutine log_notice(message, file, line)

    character(*), intent(in) :: message
    character(*), intent(in), optional :: file
    integer, intent(in), optional :: line

    logical is_redirected
    is_redirected = .not. isatty(6)

    !call time(mytime)
    mytime = ' '
    call date_and_time(mydate)

    if (present(file) .and. present(line)) then
      if (is_redirected) then
#ifdef LOG_MINFO
        if(atm_groupmember%atm_subgid .eq. (DA_Color + 1)) then
          write(fid, '("[", A7, "]","[", A17, "]",  A, ":", I0, A)') 'Notice', mydate // "/"//mytime, trim(file), line, &
              "[M: DA"// "][RANK: " // to_str(atm_groupmember%atm_subrank) // &
              "/" // to_str(atm_groupmember%atm_subsize) //"] "// trim(message)
        else
          write(fid, '("[", A7, "]","[", A17, "]",  A, ":", I0, A)') 'Notice', mydate // "/"//mytime, trim(file), line, &
              "[M: ATMGROUP"// " " //to_str(atm_groupmember%atm_subgid) // &
              "/" // to_str(atm_ensemble_group) //")"// "]" // &
              "[RANK: " // to_str(atm_groupmember%atm_subrank) // &
              "/" // to_str(atm_groupmember%atm_subsize) //"] "//trim(message)
        end if
#else
        write(6, '("[", A, "]: ", A, ":", I0, ": ", A)') 'Notice', trim(file), line, trim(message)
#endif
      else
        write(6, '("[", A, "]: ", A, ":", I0, ": ", A)') colorize('Notice', color_fg='green'), trim(file), line, trim(message)
      end if
    else
      if (is_redirected) then
#ifdef LOG_MINFO
        if(atm_groupmember%atm_subgid .eq. (DA_Color + 1)) then
          write(fid, '("[", A7, "]","[", A17, "]",  A)') 'Notice', mydate // "/"//mytime,&
              "[M: DA"// "][RANK: " // to_str(atm_groupmember%atm_subrank) // &
                          "/" // to_str(atm_groupmember%atm_subsize) //"] "// trim(message)
        else
          write(fid, '("[", A7, "]","[", A17, "]",  A)') 'Notice', mydate // "/"//mytime, &
              "[M: ATMGROUP"// " " //to_str(atm_groupmember%atm_subgid) // &
                                          "/" // to_str(atm_ensemble_group) //")"// "]" // &
                          "[RANK: " // to_str(atm_groupmember%atm_subrank) // &
                          "/" // to_str(atm_groupmember%atm_subsize) //"] "// trim(message)
        end if
#else
        write(6, '("[", A, "]: ", A)') 'Notice', trim(message)
#endif
      else
        write(6, '("[", A, "]: ", A)') colorize('Notice', color_fg='green'), trim(message)
      end if
    end if

  end subroutine log_notice

  subroutine log_warning(message, file, line)

    character(*), intent(in) :: message
    character(*), intent(in), optional :: file
    integer, intent(in), optional :: line

    logical is_redirected
    is_redirected = .not. isatty(6)

    !call time(mytime)
    mytime = ' '
    call date_and_time(mydate)

    if (present(file) .and. present(line)) then
      if (is_redirected) then
#ifdef LOG_MINFO
        if(atm_groupmember%atm_subgid .eq. (DA_Color + 1)) then
          write(fid, '("[", A7, "]","[", A17, "]",  A, ":", I0, A)') 'Warning', mydate // "/"//mytime, trim(file), line, &
              "[M: DA"// "][RANK: " // to_str(atm_groupmember%atm_subrank) // &
              "/" // to_str(atm_groupmember%atm_subsize) //"] "// trim(message)
        else
          write(fid, '("[", A7, "]","[", A17, "]",  A, ":", I0, A)') 'Warning', mydate // "/"//mytime, trim(file), line, &
              "[M: ATMGROUP"// " " //to_str(atm_groupmember%atm_subgid) // &
              "/" // to_str(atm_ensemble_group) //")"// "]" // &
              "[RANK: " // to_str(atm_groupmember%atm_subrank) // &
              "/" // to_str(atm_groupmember%atm_subsize) //"] "//trim(message)
        end if
#else
        write(6, '("[", A, "]: ", A, ":", I0, ": ", A)') 'Warning', trim(file), line, trim(message)
#endif
      else
        write(6, '("[", A, "]: ", A, ":", I0, ": ", A)') colorize('Warning', color_fg='yellow'), trim(file), line, trim(message)
      end if
    else
      if (is_redirected) then
#ifdef LOG_MINFO
        if(atm_groupmember%atm_subgid .eq. (DA_Color + 1)) then
          write(fid, '("[", A7, "]","[", A17, "]",  A)') 'Warning', mydate // "/"//mytime, &
              "[M: DA"// "][RANK: " // to_str(atm_groupmember%atm_subrank) // &
                          "/" // to_str(atm_groupmember%atm_subsize) //"] "// trim(message)
        else
          write(fid, '("[", A7, "]","[", A17, "]",  A)') 'Warning', mydate // "/"//mytime, &
              "[M: ATMGROUP"// " " //to_str(atm_groupmember%atm_subgid) // &
                                          "/" // to_str(atm_ensemble_group) //")"// "]" // &
                          "[RANK: " // to_str(atm_groupmember%atm_subrank) // &
                          "/" // to_str(atm_groupmember%atm_subsize) //"] "// trim(message)
        end if
#else
        write(6, '("[", A, "]: ", A)') 'Warning', trim(message)
#endif
      else
        write(6, '("[", A, "]: ", A)') colorize('Warning', color_fg='yellow'), trim(message)
      end if
    end if

  end subroutine log_warning

  subroutine log_error(message, file, line)

    character(*), intent(in) :: message
    character(*), intent(in), optional :: file
    integer, intent(in), optional :: line
    integer :: errorcode, ierrs

    logical is_redirected
    is_redirected = .not. isatty(6)

    !call time(mytime)
    mytime = ' '
    call date_and_time(mydate)

    if (present(file) .and. present(line)) then
      if (is_redirected) then
#ifdef LOG_MINFO
        if(atm_groupmember%atm_subgid .eq. (DA_Color + 1)) then
          write(fid, '("[", A7, "]","[", A17, "]",  A, ":", I0, A)') 'Error', mydate // "/"//mytime, trim(file), line, &
              "[M: DA"// "][RANK: " // to_str(atm_groupmember%atm_subrank) // &
              "/" // to_str(atm_groupmember%atm_subsize) //"] "// trim(message)
        else
          write(fid, '("[", A7, "]","[", A17, "]",  A, ":", I0, A)') 'Error', mydate // "/"//mytime, trim(file), line, &
              "[M: ATMGROUP"// " " //to_str(atm_groupmember%atm_subgid) // &
              "/" // to_str(atm_ensemble_group) //")"// "]" // &
              "[RANK: " // to_str(atm_groupmember%atm_subrank) // &
              "/" // to_str(atm_groupmember%atm_subsize) //"] "//trim(message)
        end if
#else
        write(6, '("[", A, "]: ", A, ":", I0, ": ", A)') 'Error', trim(file), line, trim(message)
#endif
      else
        write(6, '("[", A, "]: ", A, ":", I0, ": ", A)') colorize('Error', color_fg='red'), trim(file), line, trim(message)
      end if
    else
      if (is_redirected) then
#ifdef LOG_MINFO
        if(atm_groupmember%atm_subgid .eq. (DA_Color + 1)) then
          write(fid, '("[", A7, "]","[", A17, "]",  A, ":", I0, A)') 'Error', mydate // "/"//mytime, &
              "[M: DA"// "][RANK: " // to_str(atm_groupmember%atm_subrank) // &
              "/" // to_str(atm_groupmember%atm_subsize) //"] "// trim(message)
        else
          write(fid, '("[", A7, "]","[", A17, "]",  A)') 'Error', mydate // "/"//mytime, &
              "[M: ATMGROUP"// " " //to_str(atm_groupmember%atm_subgid) // &
              "/" // to_str(atm_ensemble_group) //")"// "]" // &
              "[RANK: " // to_str(atm_groupmember%atm_subrank) // &
              "/" // to_str(atm_groupmember%atm_subsize) //"] "//trim(message)
        end if
#else
        write(6, '("[", A, "]: ", A)') 'Error', trim(message)
#endif
      else
        write(6, '("[", A, "]: ", A)') colorize('Error', color_fg='red'), trim(message)
      end if
    end if
    call MPI_Abort(MPI_COMM_WORLD, errorcode, ierrs)

  end subroutine log_error

  subroutine log_print_diag(prefix)

    character(*), intent(in) :: prefix

    type(hash_table_iterator_type) iter

    logical is_redirected
    is_redirected = .not. isatty(6)

    if (is_redirected) then
      write(6, '(A)', advance='no') '==> ' // trim(prefix)
    else
      write(6, '(A)', advance='no') colorize('==> ', color_fg='blue') // trim(prefix)
    end if

    iter = hash_table_iterator(diags)
    do while (.not. iter%ended())
      select type (value => iter%value)
      type is (integer)
        write(6, '(X, A)', advance='no') trim(to_str(value))
      type is (real(4))
        write(6, '(X, A)', advance='no') trim(to_str(value, 14, 20))
      type is (real(8))
        write(6, '(X, A)', advance='no') trim(to_str(value, 14, 20))
      class default
        write(6, '(X, A)', advance='no') iter%key
      end select
      call iter%next()
    end do
    write(6, *)

  end subroutine log_print_diag

  subroutine log_print(message)

    character(*), intent(in) :: message

    logical is_redirected
    is_redirected = .not. isatty(6)

    !call time(mytime)
    mytime = ' '
    call date_and_time(mydate)

    if (is_redirected) then
#ifdef LOG_MINFO
      if(atm_groupmember%atm_subgid .eq. (DA_Color + 1)) then
        write(fid, '("[", A7, "]","[", A17, "]", A)') "===> ", mydate // "/"//mytime, &
           "[M: DA"// "][RANK: " // to_str(atm_groupmember%atm_subrank) // &
              "/" // to_str(atm_groupmember%atm_subsize) //"] "// trim(message)
      else
        write(fid, '("[", A7, "]","[", A17, "]",  A)') '===> ', mydate // "/"//mytime,  &
              "[M: ATMGROUP"// " " //to_str(atm_groupmember%atm_subgid) // &
              "/" // to_str(atm_ensemble_group) //")"// "]" // &
              "[RANK: " // to_str(atm_groupmember%atm_subrank) // &
              "/" // to_str(atm_groupmember%atm_subsize) //"] "//trim(message)
      end if
#else
      write(6, '(A)') '==> ' // trim(message)
#endif
    else
      write(6, '(A)') colorize('==> ', color_fg='blue') // trim(message)
    end if

  end subroutine log_print

  subroutine log_print_root(message)

    character(*), intent(in) :: message

    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

#ifdef LOG_MINFO
    if(atm_groupmember%atm_subrank == 0) then
#else
    if(myrank == 0) then
#endif
      call log_print(message)
    end if
  end subroutine

  subroutine log_print_master(message)

    character(*), intent(in) :: message

    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

    !call time(mytime)
    mytime = ' '
    call date_and_time(mydate)

    if(myrank == 0) then
#ifdef LOG_MINFO
          write(fid, '("[", A7, "]","[", A17, "]", A)') '===> ', &
              mydate // "/"//mytime, "[M: COUPLER] " // trim(message)
#else
      write(6, '(A)') '==> ' // trim(message)
#endif
    end if
  end subroutine

  subroutine log_notice_root(message, file, line)

    character(*), intent(in) :: message
    character(*), intent(in), optional :: file
    integer, intent(in), optional :: line

    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

#ifdef LOG_MINFO
    if(atm_groupmember%atm_subrank == 0) then
#else
    if(myrank == 0) then
#endif
      call log_notice(message, file, line)
    end if
  end subroutine

  subroutine log_warning_root(message, file, line)

    character(*), intent(in) :: message
    character(*), intent(in), optional :: file
    integer, intent(in), optional :: line

    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

#ifdef LOG_MINFO
    if(atm_groupmember%atm_subrank == 0) then
#else
    if(myrank == 0) then
#endif
      call log_warning(message, file, line)
    end if
  end subroutine

  subroutine log_error_root(message, file, line)

    character(*), intent(in) :: message
    character(*), intent(in), optional :: file
    integer, intent(in), optional :: line

    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

#ifdef LOG_MINFO
    if(atm_groupmember%atm_subrank == 0) then
#else
    if(myrank == 0) then
#endif
      call log_error(message, file, line)
    end if
  end subroutine

end module log_mod
