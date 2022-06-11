module pa_mod
  use mpi

  implicit none



  integer myid;
  character*50 file_name;             !output file name
  character*50 detail_file_name       !optional output file name
  character*50 stencil_time_file      !stencil_time output file
  logical pa_operator_prepare_flag    !whether do time counting
  logical pa_space_operator_flag
  logical pa_reduce_flag
  logical pa_flag

  real*8    sub_time_start , sub_time_end , cal_time_start , cal_time_end 
  real*8    tran_time_start , tran_time_end ,step_time_start , step_time_end
  real*8    sub_time , cal_time , tran_time
  real*8    time_start , time_end 
  real*8    detail_time_start , detail_time_end
  real*8    wait_time_start , wait_time_end , wait_time
  real*8    total_time_start , total_time_end
  real*8    stencil_time(20)
  real*8    stencil_time_start,stencil_time_end,stencil_total



contains

subroutine pa_init(comm , group_id_in)
  
  integer, intent(in), optional :: comm
  integer, intent(in), optional :: group_id_in

  integer ierr
  integer group_id
  integer proc_np , proc_id
  character*10 date
  character*10 time
  character*10 zone
  integer date_time(8)
  character*10 s_id 
  character*50 command_cd , command_mk
  logical istatus1,istatus2
  integer i


  if (present(group_id_in)) then
    group_id = group_id_in
  else
    group_id = 0
  endif
  
  call MPI_COMM_SIZE(comm, proc_np, ierr)
  call MPI_COMM_RANK(comm, proc_id, ierr)
  myid = proc_np * group_id + proc_id

#ifdef Detail_Time
  write(s_id,"(i6.6)") myid
  file_name(1:13)  = 'mpi_timeline_'
  file_name(14:19) = s_id
  file_name(20:23) = '.txt'
  open(unit=(10000+myid),POSITION='APPEND',file=file_name)

  stencil_time_file(1:13) = 'stencil_time_'
  stencil_time_file(14:19) = s_id
  stencil_time_file(20:23) = '.txt'
  open(unit=(1000000+myid),POSITION='APPEND',file=stencil_time_file)

  time_start = mpi_wtime()

#endif

#ifdef Detail_Calc_Time
  detail_file_name(1:17)  = 'detail_calc_time_'
  detail_file_name(18:23) = s_id
  detail_file_name(24:27) = '.txt'
  open(unit=(100000+myid),POSITION='APPEND',file=detail_file_name)
#endif 

  wait_time = 0
  do i = 1 , 20 
    stencil_time(i) = 0
  end do

  

  ! file_name(1:21)  = 'mpi_stateupdate_time_'
  ! file_name(22:25) = s_id
  ! file_name(26:29) = '.txt'
  ! open(unit=(20000+myid),POSITION='APPEND',file=file_name)

end subroutine pa_init

subroutine Get_Time_Pa(get_time)
  real*8 , intent(inout) :: get_time

  get_time = mpi_wtime();
end subroutine Get_Time_Pa

subroutine Get_Start_Time(get_time)
  real*8 , intent(inout) :: get_time

  get_time = mpi_wtime();
end subroutine Get_Start_Time

subroutine Get_End_Time(get_time)
  real*8 , intent(inout) :: get_time

  get_time = mpi_wtime();
end subroutine Get_End_Time

subroutine Get_Time_Init()
  cal_time = 0;
  tran_time = 0;
end subroutine Get_Time_Init

subroutine pa_final()

  integer i

  call Get_Time_Pa(time_end)
  write((10000+myid),*) , 97 , time_end - time_start , 'calc'
  close((10000+myid))

  do i = 1 , 12 
    stencil_total = stencil_total + stencil_time(i)
  end do
  write((1000000+myid),*) , stencil_total
  do i = 1 , 12 
    write((1000000+myid),*) , i, stencil_time(i) , stencil_time(i) / stencil_total
  end do
   close((1000000+myid))


#ifdef Detail_Calc_Time
  close((100000+myid))
#endif

end subroutine pa_final

end module pa_mod
