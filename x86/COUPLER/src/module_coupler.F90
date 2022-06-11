module coupler_mod

  use mpi
  use coupler_config
  use string

  implicit none

  character*16, dimension(:), allocatable :: hostnames
  character*16                            :: hostname
  integer, dimension(:), allocatable      :: color       !顶层atm全部资源/da分组
  integer, dimension(:), allocatable      :: atm_gcolor  !atm分组的染色
  integer, dimension(:), allocatable      :: snode_color !同一个超节点的染色
  integer, dimension(:), allocatable      :: snode_atmgcolor  !同一个超节点中的gmcore的不同分组的染色
  integer :: namelen
  integer :: local_color
  integer :: new_comm

  contains
  subroutine init_coupler_config(file_path)
    use coupler_log
    implicit none

    character(*), intent(in) :: file_path
    open(10, file=file_path, status='old')
    read(10, nml=namelist_total)
    close(10)

    call log_print_namelist()
  end subroutine

  subroutine init_coupler()
    use coupler_log

    implicit none

    integer :: i, j, n, ierr, errcode, ens
    integer :: atm_ensemble_num  !atm ensemble num in one atm ensemble group
    integer :: range(3,1)

    range(1,1) = 0
    range(2,1) = atm_stride - 1 !! only one stride!
    range(3,1) = 1
    call MPI_COMM_GROUP(MPI_COMM_WORLD, global_group, ierr)
    call MPI_GROUP_RANGE_INCL(global_group, 1, range, local_group, ierr)
    call MPI_COMM_CREATE_GROUP(MPI_COMM_WORLD, local_group, 0, dida_comm, ierr)
  
    call MPI_COMM_SIZE(dida_comm, cinfo%globalsize, ierr)
    call MPI_COMM_RANK(dida_comm, cinfo%globalrank, ierr)

    allocate(hostnames(0:cinfo%globalsize-1))
    allocate(color(0:cinfo%globalsize-1))
    allocate(atm_gcolor(0:cinfo%globalsize-1))
    allocate(snode_color(0:cinfo%globalsize-1))
    allocate(snode_atmgcolor(0:cinfo%globalsize-1))

    snode_color=DEFAULT_Color
    color=DEFAULT_Color
    atm_gcolor=DEFAULT_Color
    snode_color=DEFAULT_Color
    snode_atmgcolor=DEFAULT_Color

    call MPI_GET_PROCESSOR_NAME(hostname, namelen, ierr)
    call MPI_ALLGATHER(hostname, 16, MPI_CHARACTER, hostnames, 16, &
        MPI_CHARACTER,  dida_comm, ierr)

    if(mod(atm_ensemble_total, atm_ensemble_group) /= 0) then
        call log_error("atm_ensemble_total can not be divided by &
            atm_ensemble_num,please check namelist file")
    else
        atm_ensemble_num = atm_ensemble_total / atm_ensemble_group
    end if

    cinfo%da_comps = 1

    if ((atm_group + da_group) .lt. atm_stride) then
        cinfo%ncomps = atm_ensemble_group + cinfo%da_comps + 1
        cinfo%atm_comps = atm_ensemble_group + 1
    else
        cinfo%atm_comps = atm_ensemble_group
        cinfo%ncomps = cinfo%atm_comps + cinfo%da_comps
    endif

    ncomps = cinfo%ncomps
    cinfo%atm_np = atm_group * atm_group_num
    cinfo%da_np = da_group * da_group_num
    cinfo%total_np = cinfo%globalsize

    allocate(atm_subid(1:atm_ensemble_group))
    do i = 1, atm_ensemble_group
       atm_subid(i) = i
    end do

    DA_Id = atm_ensemble_group + 1
    DEFAULT_Id = DA_Id + 1

    if((cinfo%atm_np + cinfo%da_np) .lt. cinfo%globalsize) then
        call log_warning_root("atm procs add da procs < total proc size! &
        maybe some node is free!")
    end if

    if((cinfo%atm_np + cinfo%da_np) .gt. cinfo%globalsize) then
        call log_error_root("atm procs add da procs > total proc size!")
        call MPI_Abort(dida_comm, errcode, ierr)
    end if

    color = DEFAULT_Color
    snode_color = DEFAULT_Color
    atm_gcolor = DEFAULT_Color
    snode_atmgcolor = DEFAULT_Color

    do i = 0, atm_group_num - 1
      do j = 0, atm_group - 1
        color(atm_stride * i + atm_root + j) = ATM_Color
        snode_color(atm_stride * i + atm_root + j) = i
        atm_gcolor(atm_stride * i + atm_root + j) =  &
            (j/(atm_group/atm_ensemble_group))
        snode_atmgcolor(atm_stride * i + atm_root + j) = &
            (i)*atm_group+(j/(atm_group/atm_ensemble_group))
      end do
    end do

    do i = 0, da_group_num - 1
      do j = 0, da_group - 1
        color(da_stride * i + da_root + j) = DA_Color
        snode_color(da_stride * i + da_root + j) = i
        atm_gcolor(da_stride * i + da_root + j) = DA_Color
        snode_atmgcolor(da_stride * i + da_root + j) = DA_Color
      end do
    end do

    snode%snode_id = snode_color(cinfo%globalrank)
    atm_groupmember%atm_subgid = atm_gcolor(cinfo%globalrank) + 1
    atm_groupmember%snode_id = snode_color(cinfo%globalrank)

    call MPI_COMM_SPLIT(dida_comm, snode_color(cinfo%globalrank) ,&
                        0, snode%snode_comm, ierr)
    call MPI_COMM_SIZE(snode%snode_comm, snode%snode_size, ierr)
    call MPI_COMM_RANK(snode%snode_comm, snode%snode_rank, ierr)

    call MPI_COMM_SPLIT(dida_comm, snode_atmgcolor(cinfo%globalrank) , &
                        0, atm_groupmember%atm_snodesubcomm, ierr)
    call MPI_COMM_SIZE(atm_groupmember%atm_snodesubcomm, &
                       atm_groupmember%atm_snodesubsize, ierr)
    call MPI_COMM_RANK(atm_groupmember%atm_snodesubcomm, &
                       atm_groupmember%atm_snodesubrank, ierr)

    call MPI_COMM_SPLIT(dida_comm, atm_gcolor(cinfo%globalrank) , &
                        0, atm_groupmember%atm_subcomm, ierr)
    call MPI_COMM_SIZE(atm_groupmember%atm_subcomm, &
                       atm_groupmember%atm_subsize, ierr)
    call MPI_COMM_RANK(atm_groupmember%atm_subcomm, &
                       atm_groupmember%atm_subrank, ierr)

    call MPI_COMM_SPLIT(dida_comm, color(cinfo%globalrank), &
                        0, new_comm, ierr)

    !if(atm_gcolor(cinfo%globalrank) == 0) then
    !  call init_snode_info()
    !end if

    !call MPI_Bcast(snode%s_lon_ibeg, 1, MPI_INT, 0,snode%snode_comm, ierr)
    !call MPI_Bcast(snode%s_lon_iend, 1, MPI_INT, 0,snode%snode_comm, ierr)
    !call MPI_Bcast(snode%s_lat_ibeg, 1, MPI_INT, 0,snode%snode_comm, ierr)
    !call MPI_Bcast(snode%s_lat_iend, 1, MPI_INT, 0,snode%snode_comm, ierr)

    call log_print_couplerinfo()
    end subroutine

  subroutine init_snode_info(group_id, xst, xed, yst, yed)
    use coupler_config

    implicit none
    integer, intent(in) :: xst
    integer, intent(in) :: xed
    integer, intent(in) :: yst
    integer, intent(in) :: yed
    integer, intent(in) :: group_id
    integer :: ierrs
    integer :: local_rank
    integer :: i, j, k
    
    !if (group_id == 0 .and. snode%snode_rank == 0) then
    if (group_id == 0) then
         
        call MPI_Reduce(xst,snode%snodeinfo(1),1,MPI_INT,MPI_MIN,0, &
                        atm_groupmember%atm_snodesubcomm, ierrs)
        call MPI_Reduce(xed,snode%snodeinfo(2),1,MPI_INT,MPI_MAX,0, &
                        atm_groupmember%atm_snodesubcomm, ierrs)
        call MPI_Reduce(yst,snode%snodeinfo(3),1,MPI_INT,MPI_MIN,0, &
                        atm_groupmember%atm_snodesubcomm, ierrs)
        call MPI_Reduce(yed,snode%snodeinfo(4),1,MPI_INT,MPI_MAX,0, &
                        atm_groupmember%atm_snodesubcomm, ierrs)

        !call MPI_Bcast(snode%snodeinfo, 4, MPI_INT, 0, atm_groupmember%atm_snodesubcomm, ierrs)
    end if


    if (snode%snode_rank == 0) then
        do i = 1, snode%snode_size - 1
            call MPI_Send(snode%snodeinfo, 4, MPI_INT, i, 0, snode%snode_comm, ierrs)
        end do

        snode%s_lon_ibeg = snode%snodeinfo(1)
        snode%s_lon_iend = snode%snodeinfo(2)
        snode%s_lat_ibeg = snode%snodeinfo(3)
        snode%s_lat_iend = snode%snodeinfo(4)
    else
        call wait_snode_info("ATM")
    end if

  end subroutine

  subroutine wait_snode_info(model_tag)
    use coupler_config
    use coupler_log

    implicit none
    character(*), intent(in) :: model_tag
    integer ierrs
    integer status(MPI_STATUS_SIZE)
    character(256) :: hname 
    integer hlen

    call MPI_Recv(snode%snodeinfo, 4, MPI_INT, 0, 0, snode%snode_comm, status, ierrs)

    snode%s_lon_ibeg = snode%snodeinfo(1)
    snode%s_lon_iend = snode%snodeinfo(2)
    snode%s_lat_ibeg = snode%snodeinfo(3)
    snode%s_lat_iend = snode%snodeinfo(4)

!#ifdef DEBUG_OUT
    call MPI_GET_PROCESSOR_NAME(hname, hlen, ierrs)
    call log_notice(trim(model_tag) // '-SINFO: sid-host['// to_str(snode%snode_id) // &
          ' ' // trim(hname) // & 
          '] its-ite['// to_str(snode%s_lon_ibeg) // ' ' // to_str(snode%s_lon_iend) // &
          '] jts-jte['// to_str(snode%s_lat_ibeg) // ' ' // to_str(snode%s_lat_iend) // ']')
!#endif
  end subroutine
  subroutine DEFAULT_DRIVER(MPI_DEFAULT_GROUP)

    use m_MCTWorld,only: MCTWorld_init => init

    implicit none
    integer,intent(in) :: MPI_DEFAULT_GROUP

    call MCTWorld_init(ncomps, dida_comm, MPI_DEFAULT_GROUP, DEFAULT_Id)

  end subroutine DEFAULT_DRIVER

  subroutine clean_default_coupler()

    use m_MCTWorld,only: MCTWorld_clean => clean

    implicit none

    call MCTWorld_clean()

    deallocate(color)
    deallocate(atm_gcolor)
    deallocate(snode_atmgcolor)
    deallocate(snode_color)
    deallocate(hostnames)
    deallocate(atm_subid)

  end subroutine

end module coupler_mod

