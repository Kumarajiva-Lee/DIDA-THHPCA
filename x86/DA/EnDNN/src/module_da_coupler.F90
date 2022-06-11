#include "../../../mp_config/mct_config.inc"

#define INDEX(i, j, k, N, L) ((i-1)*(N)*(L)+(j-1)*(L)+(k-1))

module da_coupler

  use m_MCTWorld,only: MCTWorld_init => init
  use m_MCTWorld,only: MCTWorld_clean => clean
  use m_GlobalSegMap,only: GlobalSegMap
  use m_GlobalSegMap,only: GlobalSegMap_init => init
  use m_GlobalSegMap,only: GlobalSegMap_lsize => lsize
  use m_GlobalSegMap,only: GlobalSegMap_clean => clean
  use m_GlobalSegMap,only: GlobalSegMap_Ordpnts => OrderedPoints
  use m_Rearranger,only: Rearranger
  use m_Rearranger,only: Rearranger_init => init
  use m_Rearranger,only: Rearranger_clean => clean
  use m_Rearranger,only: Rearranger_rearrange => rearrange
  use m_AttrVect,only    : AttrVect
  use m_AttrVect,only    : AttrVect_init => init
  use m_AttrVect,only    : AttrVect_zero => zero
  use m_AttrVect,only    : AttrVect_clean => clean
  use m_AttrVect,only    : AttrVect_indxR => indexRA
  use m_AttrVect,only    : AttrVect_importRAttr => importRAttr
  use m_AttrVect,only    : AttrVect_exportRAttr => exportRAttr
  use m_Router,only: Router
  use m_Router,only: Router_init => init
  use m_Router,only: Router_clean => clean
  use m_Transfer,only : MCT_Send => send
  use m_Transfer,only : MCT_Recv => recv

  use mpi
  use da_parallel
  use coupler_config
  use pro_info
  use redis_mod
  use da_output_mod
  use da_namelist_mod, only:obs_stat_interval_lon, obs_stat_interval_lat, &
  obs_stat_interval_vertical, da_halo
  use string
  use splitchar
  use var_info_mod
  use flogger
  !use process_mod

  implicit none

  type(GlobalSegMap)   :: DA_GSMap, DA_N_GSMap            ! MCT defined type
  type(AttrVect) :: DA_AV, DA_AV_N, DA_Cal, DA_Cal_N                      ! MCT defined type
  type(Rearranger) :: DA_Re, DA_Re_back
  type(Router),allocatable,dimension(:) :: R_DA2ATM             ! MCT defined type
  !type(Router)          :: R_DA2ATM                      ! MCT defined type

  integer, allocatable,dimension(:) :: start,length
  integer, allocatable,dimension(:) :: start_t,length_t
  integer, allocatable,dimension(:) :: start_n,length_n
  axb_type, dimension(:), pointer :: avdata, area_count, area_count_no_overlap
  integer :: avsize, av_nsize, size_nooverlap
  !type(c_ptr) :: da_rc
  integer :: myid

  contains

subroutine init_da_coupler(MPI_DA_GROUP)

  implicit none
  integer, intent(in) :: MPI_DA_GROUP
  integer :: i, j, k, mem_num, mem_num_local, ij
  integer :: ierr
  integer :: halo2
  integer :: gnum
  integer :: ylength_max
  integer :: groupid

  call log_notice("da_couper init")

  call MPI_Comm_rank(MPI_DA_GROUP, myid, ierr)

  halo2 = da_halo * 2

  call MCTWorld_init(ncomps, dida_comm, MPI_DA_GROUP, DA_Id)

  ylength_max = yed-yst+1+halo2
  allocate(start_t(num_lev*ylength_max*3))
  allocate(length_t(num_lev*ylength_max*3))

  mem_num = atm_ensemble_total
  mem_num_local = atm_ensemble_total / atm_ensemble_group

  gnum = num_lev*num_lat*num_lon*mem_num_local

  ij = 0
  if (xst - da_halo >= 1 .and. xed + da_halo <= num_lon) then
    do j = 1, num_lev
      do i = yst-da_halo, yed+da_halo
        if (i >= 1 .and. i <= num_lat) then
          ij = ij + 1
          length_t(ij) = (xed - xst + 1 + halo2) * mem_num_local
          start_t(ij) = INDEX(j, i, xst-da_halo, num_lat, num_lon) * mem_num_local + 1
        end if
      end do
    end do
  end if
  if (xst - da_halo < 1 .and. xed + da_halo <= num_lon) then
    do j = 1, num_lev
      do i = yst-da_halo, yed+da_halo
        if (i >= 1 .and. i <= num_lat) then
          ij = ij + 1
          length_t(ij) = da_halo * mem_num_local
          start_t(ij) = INDEX(j, i, num_lon-da_halo+1, num_lat, num_lon) * mem_num_local + 1
          ij = ij + 1
          length_t(ij) = (xed - xst + 1 + da_halo) * mem_num_local
          start_t(ij) = INDEX(j, i, xst, num_lat, num_lon) * mem_num_local + 1
        end if
      end do
    end do
  end if
  if (xed + da_halo > num_lon .and. xst - da_halo >= 1) then
    do j = 1, num_lev
      do i = yst-da_halo, yed+da_halo
        if (i >= 1 .and. i <= num_lat) then
          ij = ij + 1
          length_t(ij) = (xed - xst + 1 + da_halo) * mem_num_local
          start_t(ij) = INDEX(j, i, xst-da_halo, num_lat, num_lon) * mem_num_local + 1
          ij = ij + 1
          length_t(ij) = da_halo * mem_num_local
          start_t(ij) = INDEX(j, i, 1, num_lat, num_lon) * mem_num_local + 1
        end if
      end do
    end do
  end if
  if (xst - da_halo < 1 .and. xed + da_halo > num_lon) then
    do j = 1, num_lev
      do i = yst-da_halo, yed+da_halo
        if (i >= 1 .and. i <= num_lat) then
          ij = ij + 1
          length_t(ij) = da_halo * mem_num_local
          start_t(ij) = INDEX(j, i, num_lon-da_halo+1, num_lat, num_lon) * mem_num_local + 1
          ij = ij + 1
          length_t(ij) = (xed - xst + 1) * mem_num_local
          start_t(ij) = INDEX(j, i, xst, num_lat, num_lon) * mem_num_local + 1
          ij = ij + 1
          length_t(ij) = da_halo * mem_num_local
          start_t(ij) = INDEX(j, i, 1, num_lat, num_lon) * mem_num_local + 1
        end if
      end do
    end do
  end if

  allocate(start(ij))
  allocate(length(ij))
  start(1:ij) = start_t(1:ij)
  length(1:ij) = length_t(1:ij)

  call GlobalSegMap_init(DA_GSMap, start, length, 0, MPI_DA_GROUP, DA_Id, gsize=gnum)
  deallocate(start)
  deallocate(length)
  call log_notice("segmap1 init")

  ij = 0
  do j = 1, num_lev
    do i = yst, yed
      ij = ij + 1
      length_t(ij) = (xed - xst + 1) * mem_num_local
      start_t(ij) = INDEX(j, i, xst, num_lat, num_lon) * mem_num_local + 1
    end do
  end do
  
  allocate(start_n(ij))
  allocate(length_n(ij))
  start_n(1:ij) = start_t(1:ij)
  length_n(1:ij) = length_t(1:ij)
  call GlobalSegMap_init(DA_N_GSMap, start_n, length_n, 0, MPI_DA_GROUP, DA_Id, gsize=gnum)
  call log_notice("segmap2 init")

  call Rearranger_init(DA_GSMap, DA_N_GSMAP, MPI_DA_GROUP, DA_Re)
  call Rearranger_init(DA_N_GSMAP, DA_GSMap, MPI_DA_GROUP, DA_Re_back)
  call log_notice("rearranger init")

  deallocate(start_t)
  deallocate(length_t)
  deallocate(start_n)
  deallocate(length_n)

  avsize = GlobalSegMap_lsize(DA_GSMap,MPI_DA_GROUP)
  av_nsize = GlobalSegMap_lsize(DA_N_GSMAP,MPI_DA_GROUP)

  !print *, da_var_name

  call AttrVect_init(DA_AV,rList="ps:u:v:pt",lsize=avsize)
  call AttrVect_init(DA_Cal, rList="num", lsize=avsize)
  call AttrVect_init(DA_AV_N,rList="ps:u:v:pt",lsize=av_nsize)
  call AttrVect_init(DA_Cal_N, rList="num", lsize=av_nsize)

  !print *, da_var_name

  allocate(R_DA2ATM(1:atm_ensemble_group))
  do groupid = 1, atm_ensemble_group
    call Router_init(atm_subid(groupid),DA_GSMap,MPI_DA_GROUP,R_DA2ATM(groupid))
  end do
  call log_notice("router init")

  allocate(area_count(avsize))
  do i = 1, avsize
    area_count(i) = 1.0
  end do
  call AttrVect_importRAttr(DA_Cal, "num", area_count)
  call Rearranger_rearrange(DA_Cal, DA_Cal_N, DA_Re, sum=.true.) 
  call AttrVect_exportRAttr(DA_Cal_N, "num", area_count_no_overlap, size_nooverlap)

  deallocate(area_count)
  call AttrVect_clean(DA_Cal)
  call AttrVect_clean(DA_Cal_N)

end subroutine 


subroutine dafatm_coupling(MPI_DA_GROUP)

  use datetime
  use coupler_config, only:da_start_time_array=>da_start_time, end_time_array=>end_time, &
                           da_in_seconds

  implicit none
  integer, intent(in) :: MPI_DA_GROUP
  integer :: myrank
  integer :: ierr
  integer :: i, ls, j, ij, ts, da_time
  integer :: halo2
  integer :: mem_num, mem_num_local

  ! x : longitude, y : latitude
  integer :: nx, ny
  integer :: lxst, lxed, lyst, lyed
  integer :: ix, iy, iz, debug_int, niz
  integer :: mem_i, mem_st, mem_ed, mem_n

  axb_type, dimension(:), pointer :: recvdata, recv_nooverlap
  ! N = num of vars
  integer :: n_var
  ! n = num of n dimension (2d + 3d * lev)
  integer :: n_max
  integer :: n_start, ens_start, obs_num_2d, obs_num_3d
  integer, allocatable, dimension(:) :: var_lev
  character(len = 20), dimension(10) :: var_name
  character(len = 20), allocatable, dimension(:) :: var_name_ordered
  character(len = 20) :: tname
  integer :: groupid
  character(12) :: curr_time_str
  
  axb_type, allocatable, dimension(:, :, :, :) :: da_data_t
  axb_type, allocatable, dimension(:, :, :, :) :: da_data
  axb_type, allocatable, dimension(:, :, :, :) :: result_data
  axb_type, allocatable, dimension(:, :, :) :: da_data_mean
  axb_type, allocatable, dimension(:, :, :) :: result_data_mean

  real(8), allocatable, dimension(:, :, :) :: spread_b
  real(8), allocatable, dimension(:, :, :) :: spread_a

  real(8), allocatable, dimension(:, :, :, :) :: data_middle
  real(8), allocatable, dimension(:, :, :, :) :: data_to_svr
  real(8), allocatable, dimension(:, :, :, :) :: data_from_svr
  integer :: data_info(10)
  integer :: send_data_sz, recv_data_sz, global_rank, da_size
  integer :: istat(mpi_status_size)

  type(datetime_type) da_start_time
  type(datetime_type) end_time
  type(datetime_type) curr_time
  type(timedelta_type) delta_time
  

  nx = xed - xst + 1
  ny = yed - yst + 1

  lxst = 1 - da_halo; lxed = nx + da_halo
  lyst = 1 - da_halo; lyed = ny + da_halo
  if (yst == 1) lyst = 1
  if (yed == num_lat) lyed = ny

  mem_num = atm_ensemble_total
  mem_num_local = atm_ensemble_total / atm_ensemble_group

  ! !!!!!!!!!!!!!
  ! vars parsing
  call stringsplit(da_var_name,',',var_name,n_var)
  allocate(var_name_ordered(n_var))
  allocate(var_lev(n_var))
  call var_info(var_name, n_var, var_name_ordered, obs_num_2d, obs_num_3d)
  n_max = 0
  do i = 1, obs_num_2d
    var_lev(i) = 1
    n_max = n_max + var_lev(i)
  end do
  do i = obs_num_2d + 1, n_var
    var_lev(i) = num_lev
    n_max = n_max + var_lev(i)
  end do
  ! !!!!!!!!!!!!!
  ! !!!!!!!!!!!!!

  ! var_name(1) = "u"; var_name(2) = "v"; var_name(3) = "gz"
  ! var_lev(1) = num_lev; var_lev(2) = num_lev; var_lev(3) = num_lev
  ! n_max = 3 * num_lev

  allocate(da_data_t(mem_num_local, 1-da_halo:nx+da_halo, &
  1-da_halo:ny+da_halo, num_lev))
  allocate(da_data(svr_ensemble, n_max, &
  1-da_halo:ny+da_halo, 1-da_halo:nx+da_halo))
  allocate(result_data(svr_ensemble, n_max, &
  1-da_halo:ny+da_halo, 1-da_halo:nx+da_halo))
  allocate(da_data_mean(n_max, 1:ny, 1:nx))
  allocate(result_data_mean(n_max, 1:ny, 1:nx))
  allocate(spread_b(n_max, 1:ny, 1:nx))
  allocate(spread_a(n_max, 1:ny, 1:nx))

  allocate(data_middle(mem_num,n_max,1-da_halo:ny+da_halo,1-da_halo:nx+da_halo))
  allocate(data_to_svr(mem_num,n_max,1:ny,1:nx))
  allocate(data_from_svr(svr_ensemble,n_max,1:ny,1:nx))

  allocate(recvdata(1:avsize))
  allocate(recv_nooverlap(1:av_nsize))

  halo2 = 2 * da_halo

  da_start_time = create_datetime(year=da_start_time_array(1),  &
                                 month=da_start_time_array(2), &
                                 day=da_start_time_array(3),   &
                                 hour=da_start_time_array(4),  &
                                 minute=da_start_time_array(5))

  end_time      = create_datetime(year=end_time_array(1),  &
                                 month=end_time_array(2), &
                                 day=end_time_array(3),   &
                                 hour=end_time_array(4),  &
                                 minute=end_time_array(5))

  delta_time = end_time - da_start_time
  !call MCT_Recv(DA_AV, R_DA2ATM)
  !print *, "Before receiving"
  
  da_time = ceiling((delta_time%days*86400+delta_time%hours*3600+delta_time%minutes*60+delta_time%seconds)/da_in_seconds)
  !da_time = da_time + 1
  !da_time = 3 !Warning !
  !print *, "before recving"
  call MPI_COMM_RANK(MPI_DA_GROUP, myrank, ierr)
  if (myrank == 0) then
    call log_notice("da_time " // to_str(da_time))
  endif

  do ts = 1, da_time

    ens_start = 1

    do groupid = 1, atm_ensemble_group
      !print *, n_var, groupid
      call MCT_Recv(DA_AV, R_DA2ATM(groupid))
      !print *, mem_num_local

      n_start = 1
      do i = 1, n_var

        !print *, trim_len(var_name_ordered(i))
        tname = trim(var_name_ordered(i))

        call log_notice(tname)

        !print *, myrank, "len", len_trim(tname)
        call AttrVect_exportRAttr(DA_AV, &
        static_string(tname, len_trim(tname)), &
        recvdata, ls)
        !call AttrVect_exportRAttr(DA_AV, "u", recvdata, ls)

        !da_data_t(lxst:lxed, lyst:lyed, :) = recvdata(1:ls)
        !print *, ls, mem_num_local*(lxed-lxst+1)*(lyed-lyst+1)*var_lev(i)

        !if (i == 1) then
        !  open(12,file='test_out_'//to_str(myrank)//'.txt')
        !end if

        debug_int = 0
        do iz = 1, var_lev(i)
          do iy = lyst, lyed
            do ix = lxst, lxed
              da_data_t(1:mem_num_local, ix, iy, iz) = recvdata(debug_int + 1:debug_int + mem_num_local)
              if (isnan(da_data_t(1, ix, iy, iz))) print *, rank, ix, iy, iz, "recv data nan"
              debug_int = debug_int + mem_num_local
              !if (i == 1) then
              !  write (12,*), iz, iy, ix, da_data_t(1,ix,iy,iz)
              !end if
            end do
          end do
        end do
        !if (i == 1) then
        !  close(12)
        !end if

        ! da_data_t(1:mem_num, lxst:lxed, lyst:lyed, 1:var_lev(i)) = &
        ! reshape(recvdata(1:ls), (/mem_num,lxed-lxst+1,lyed-lyst+1,var_lev(i)/))

        ! (ens_num, lon_num, lat_num, n_2d + lev_num * n_3d) 
        ! transform to 
        !(ens_num, n_2d + lev_num * n_3d, lat_num, lon_num)
        call log_notice(to_str(mem_num_local)//' '//to_str(mem_num))
        call struct_atm2endnn(mem_num_local, var_lev(i), ny+halo2, nx+halo2, mem_num, n_max, ens_start, n_start, &
        da_data_t(:, :, :, :), data_middle(:, :, :, :))
        !result_data(:, :, :, :) = data_middle(:, :, :, :)

        n_start = n_start + var_lev(i)

      end do

      ens_start = ens_start + mem_num_local
      
    end do


    ! gmcore[ens, lon_with_halo, lat_with_halo, lev]
    ! to 
    ! LETKF [ens, lev, lat_with_halo, lon_with_halo] to svr
    ! ens = atm_ensemble_total, need to change, only 1 gmcore is running
    ! mem_num = atm_ensemble_total(also ensemble_num generated from svr)
    ! only one ensemble send to svr
    !open(12,file='out_'//to_str(myrank)//'.txt')
    ! [ps,u,v,pt] to [ps,pt,u,v]
    do iz = 1, n_max
      do iy = 1, ny
        do ix = 1, nx
          if (iz == 1) then !ps -> ps
            niz = iz
          else if (iz >= 2 .and. iz <= 33) then ! pt(svr) 
            niz = iz + 64
          else !u,v(svr)
            niz = iz - 32
          end if
          data_to_svr(1:mem_num,iz,iy,ix) = data_middle(1:mem_num,niz,iy,ix)
          !if (iz == 1) then
          !  write(12,*), iz, iy, ix, data_to_svr(iz,iy,ix)
          !end if
        end do 
      end do
    end do
    !close(12)
    send_data_sz = n_max * ny * nx * mem_num
    call MPI_Comm_rank(MPI_COMM_WORLD, global_rank, ierr)
    call MPI_Comm_size(MPI_DA_GROUP, da_size, ierr)

    data_info(9) = mem_num !atm_ensemble_total
    data_info(1) = svr_ensemble
    data_info(2) = n_max
    data_info(3) = ny
    data_info(4) = nx
    data_info(5) = yst
    data_info(6) = yed
    data_info(7) = xst
    data_info(8) = xed
    data_info(10) = da_time
    
    call log_notice("send data_info to svr")
    !Send data_info
    call MPI_Send(data_info, 10, MPI_INT, global_rank+da_size, 0, MPI_COMM_WORLD, ierr)
    ! Send to svr
    call log_notice("send to svr")
    call MPI_Send(data_to_svr, send_data_sz, MPI_DOUBLE, global_rank+da_size, 0, MPI_COMM_WORLD, ierr)
    ! Recv from svr
    recv_data_sz = svr_ensemble * n_max * ny * nx
    call log_notice("wait for svr")
    call MPI_Recv(data_from_svr, recv_data_sz, MPI_DOUBLE, global_rank+da_size, 0, MPI_COMM_WORLD, istat, ierr)
    
    ! svr to LETKF [ens, lev, lat_with_halo, lon_with_halo] 
    call log_notice("recv from svr")
    ! [ps,pt,u,v] to [ps,u,v,pt]
    do iz = 1, n_max
      do iy = 1, ny
        do ix = 1, nx
          if (iz == 1) then
            niz = iz
          else if (iz >= 2 .and. iz <= 33) then
            niz = iz + 64
          else 
            niz = iz - 32
          end if
          da_data(1:svr_ensemble, niz,iy,ix) = data_from_svr(1:svr_ensemble,iz,iy,ix)
        end do 
      end do
    end do

    call AttrVect_zero(DA_AV_N, .true., .true.)

    !open(12,file="out"//to_str(myrank)//'.txt')

    mem_n = svr_ensemble / mem_num_local
    do mem_i = 1, mem_n 
      recv_nooverlap(:) = 0
      mem_st = (mem_i - 1) * mem_num_local + 1
      mem_ed = mem_st + mem_num_local - 1
    ! gmcore[ens, lon_with_halo, lat_with_halo, lev]
    ! letkf[ens, lev, lat, lon]
      n_start = 0
      do i = 1, n_var
        tname = trim(var_name_ordered(i))
        debug_int = 0
        do iz = 1, var_lev(i)
          do iy = 1, ny
            do ix = 1, nx
              recv_nooverlap(debug_int+1:debug_int+mem_num_local) = da_data(mem_st:mem_ed,n_start+iz,iy,ix)
              debug_int = debug_int + mem_num_local
            end do
          end do
        end do
        n_start = n_start + var_lev(i)
        call AttrVect_importRAttr(DA_AV_N, &
        static_string(tname, len_trim(tname)), &
        recv_nooverlap)
      end do
      call Rearranger_rearrange(DA_AV_N, DA_AV, DA_Re_back)
      n_start = 0
      do i = 1, n_var
        tname = trim(var_name_ordered(i))
        call AttrVect_exportRAttr(DA_AV, &
        static_string(tname, len_trim(tname)), &
        recvdata, ls)
        debug_int = 0
        do iz = 1, var_lev(i)
          do iy = lyst, lyed
            do ix = lxst, lxed
              da_data(mem_st:mem_ed, n_start+iz, iy, ix) = recvdata(debug_int + 1:debug_int + mem_num_local)
              debug_int = debug_int + mem_num_local
     !         write (12,*) mem_st, n_start+iz, iy, ix, da_data(mem_st, n_start+iz, iy, ix)
              if (isnan(da_data(mem_st, n_start+iz, iy, ix))) then
                call log_error("receive nan data from svr!")
              end if
            end do
          end do
        end do
        n_start = n_start + var_lev(i)
      end do

    end do
    !close(12)
    ! backto gmcore
    ! 

    !result_data(:,:,:,:)=da_data(:,:,:,:)
    call run_da_parallel(n_max, ny, nx, ts-1, da_data, result_data, spread_b, spread_a)
    !print *, da_data(1, 1, 1, 1), result_data(1, 1, 1, 1)
    data_middle(1:mem_num,:,:,:) = result_data(1:mem_num,:,:,:)

    !call log_notice("run da parallel done!")

    ens_start = 1

    do groupid = 1, atm_ensemble_group
      n_start = 1

      do i = 1, n_var

        da_data_t(:, :, :, :) = 0

        call struct_endnn2atm(mem_num_local, var_lev(i), ny+halo2, nx+halo2, mem_num, n_max, ens_start, n_start, &
        data_middle(:, :, :, :), da_data_t(:, :, :, :))

        !da_data_t(:, lxst:lxed, lyst:lyed, :) = 0

        debug_int = 0
        do iz = 1, var_lev(i)
          do iy = lyst, lyed
            do ix = lxst, lxed
              recvdata(debug_int + 1:debug_int + mem_num_local) = da_data_t(1:mem_num_local, ix, iy, iz)
              debug_int = debug_int + mem_num_local
            end do
          end do
        end do

        !recvdata(1:ls) = &
        !reshape(da_data_t(1:mem_num, lxst:lxed, lyst:lyed, 1:var_lev(i)), (/ls/))
        tname = trim(var_name_ordered(i))
        call AttrVect_importRAttr(DA_AV, &
        static_string(tname, len_trim(tname)), &
        recvdata)

        n_start = n_start + var_lev(i)

      end do

      ! call log_notice("before_rearrange_1")
      ! call Rearranger_rearrange(DA_AV, DA_AV_N, DA_Re, sum=.true.) 
      ! call log_notice("rearrange_1_ccomplete")
      ! do i = 1, n_var

      !   !recvdata(1:ls) = &
      !   !reshape(da_data_t(1:mem_num, lxst:lxed, lyst:lyed, 1:var_lev(i)), (/ls/))
      !   tname = trim(var_name_ordered(i))
      !   call AttrVect_exportRAttr(DA_AV_N, &
      !   static_string(tname, len_trim(tname)), &
      !   recv_nooverlap, ls)

      !   ! do j = 1, ls
      !   !   recv_nooverlap(j) = recv_nooverlap(j) / area_count_no_overlap(j)
      !   ! end do

      !   call AttrVect_importRAttr(DA_AV_N, &
      !   static_string(tname, len_trim(tname)), &
      !   recv_nooverlap)

      ! end do
      ! call log_notice("before_rearrange_2")
      ! call Rearranger_rearrange(DA_AV_N, DA_AV, DA_Re_back) 
      ! call log_notice("rearrange_2_ccomplete")


      call MCT_Send(DA_AV, R_DA2ATM(groupid))

      ens_start = ens_start + mem_num_local

    end do

    ! (ens_num, n_2d + lev_num * n_3d, lat_num, lon_num)
    do ix = 1, nx
      do iy = 1, ny
        do iz = 1, n_max
          da_data_mean(iz, iy, ix) = sum(da_data(1:mem_num, iz, iy, ix)) / mem_num
        end do
      end do
    end do
    do ix = 1, nx
      do iy = 1, ny
        do iz = 1, n_max
          result_data_mean(iz, iy, ix) = sum(result_data(1:mem_num, iz, iy, ix)) / mem_num
        end do
      end do
    end do

   ! print *, "haha", result_data_mean(26,61,85)
   !print *, "haha", result_data_mean(2,61,85)
   ! print *, "haha", result_data_mean(72,61,85)

    curr_time = da_start_time+create_timedelta(seconds=(ts-1)*da_in_seconds)
    curr_time_str = curr_time%format('%Y%m%d%H%M')
    !print *, 'curr time:', curr_time_str
    !print *,'var_name_ordered',var_name_ordered
#ifdef DAOUT_NC
    call da_output_nc(myrank, MPI_DA_GROUP,da_rc, yst, yed, xst, xed, obs_num_2d, obs_num_3d, var_name_ordered, &
     curr_time_str, da_data_mean, result_data_mean, spread_b, spread_a)
#else
    call da_output_redis(myrank, MPI_DA_GROUP,da_rc, yst, yed, xst, xed, obs_num_2d, obs_num_3d, var_name_ordered, &
     curr_time_str, da_data_mean, result_data_mean, spread_b, spread_a)
#endif
     !print *, myid, "redis done!"

    
    ! open(51, FILE='test_da_'//to_str(myrank)//'.txt')
    ! write (51, *) myrank, ls, xst, xed, yst, yed
    ! close(51)

  end do

  call log_notice("Out da run")

  deallocate(da_data)
  deallocate(da_data_t)
  deallocate(result_data)
  deallocate(da_data_mean)
  deallocate(result_data_mean)
  deallocate(var_name_ordered)
  deallocate(var_lev)
  deallocate(spread_b)
  deallocate(spread_a)
  deallocate(data_to_svr)
  deallocate(data_from_svr)

end subroutine 

subroutine clean_da_coupler()

  implicit none
  integer :: groupid

  deallocate(area_count_no_overlap)
  do groupid = 1, atm_ensemble_group
     call Router_clean(R_DA2ATM(groupid))
  end do
  deallocate(R_DA2ATM)
  !call Router_clean(R_DA2ATM)
  call Rearranger_clean(DA_Re)
  call Rearranger_clean(DA_Re_back)
  call AttrVect_clean(DA_AV)
  call AttrVect_clean(DA_AV_N)
  call GlobalSegMap_clean(DA_GSMap)
  call GlobalSegMap_clean(DA_N_GSMap)
  call MCTWorld_clean()

end subroutine

#include "data_trans.F90"

! Fortran string
function static_string(f_str, f_len)
  character(20), intent(in) :: f_str
  integer, intent(in) :: f_len
  character(len=f_len) static_string
  integer :: i
  character :: ts

  static_string = ""
  do i = 1, f_len
    ts = f_str(i:i)
    static_string(i:i) = ts
  end do
end function static_string

end module da_coupler
