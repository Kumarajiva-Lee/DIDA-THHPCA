module history_mod

  use fiona
  use flogger
  use string
  use const_mod
  use namelist_mod
  use parallel_mod
  use allocator_mod
  use time_mod, dt => dt_in_seconds
  use mesh_mod
  use state_mod
  use static_mod
  use tend_mod
  use block_mod
  use diag_state_mod
  use process_mod
  use pa_mod
  use member_mod

  implicit none

  private

  public history_init
  public history_final
  public history_write_state
  public history_write_debug

contains

  subroutine history_init( group_id_in)

    !integer, intent(in) :: member_num
    integer , intent(in) , optional :: group_id_in

    character(10) time_value, time_units
    character(4) cell_dims(4), cell_dims_2d(3)
    character(4) lon_dims(4), lon_dims_2d(3)
    character(4) lat_dims(4), lat_dims_2d(3)
    character(4) vtx_dims(4), vtx_dims_2d(3)
    character(4) lon_lev_dims(4), lat_lev_dims(4), lev_dims(4)
    character(4) diag_dims(4)
    real(8) seconds, months

    integer im
    character*20 fid0 , fid1
    integer group_id
    character*10 s_id 
    character*50 io_file_name;

    if (present(group_id_in)) then
      group_id = group_id_in
    else
      group_id = 0
    endif


#ifdef Detail_io_Time 
    write(s_id,"(i4.4)") myid
    io_file_name(1:8)  = 'io_time_'
    io_file_name(9:12) = s_id
    io_file_name(13:16) = '.txt'
    open(unit=(23333+myid),POSITION='APPEND',file=io_file_name)
#endif

    if (history_interval(1) == 'N/A') call log_error('Parameter history_interval is not set!')
    if (case_name == 'N/A') call log_error('Parameter case_name is not set!')

    time_value = split_string(history_interval(1), ' ', 1)
    time_units = split_string(history_interval(1), ' ', 2)
    read(time_value, *) seconds
    select case (time_units)
    case ('days')
      seconds = seconds * 86400
    case ('hours')
      seconds = seconds * 3600
    case ('minutes')
      seconds = seconds * 60
    case ('seconds')
      seconds = seconds
    case default
      if (is_root_proc()) call log_error('Invalid history interval ' // trim(history_interval(1)) // '!')
    end select

    call time_add_alert('history_write', seconds=seconds , first_alert=.true.)

        cell_dims(1) =  'lon';     cell_dims(2) =  'lat';     cell_dims(3) =  'lev';     cell_dims(4) = 'time'
         lon_dims(1) = 'ilon';      lon_dims(2) =  'lat';      lon_dims(3) =  'lev';      lon_dims(4) = 'time'
     lon_lev_dims(1) = 'ilon';  lon_lev_dims(2) =  'lat';  lon_lev_dims(3) = 'ilev';  lon_lev_dims(4) = 'time'
         lat_dims(1) =  'lon';      lat_dims(2) = 'ilat';      lat_dims(3) =  'lev';      lat_dims(4) = 'time'
     lat_lev_dims(1) =  'lon';  lat_lev_dims(2) = 'ilat';  lat_lev_dims(3) = 'ilev';  lat_lev_dims(4) = 'time'
         vtx_dims(1) = 'ilon';      vtx_dims(2) = 'ilat';      vtx_dims(3) =  'lev';      vtx_dims(4) = 'time'
         lev_dims(1) =  'lon';      lev_dims(2) =  'lat';      lev_dims(3) = 'ilev';      lev_dims(4) = 'time'
     cell_dims_2d(1) =  'lon';  cell_dims_2d(2) =  'lat';  cell_dims_2d(3) = 'time'
      lon_dims_2d(1) = 'ilon';   lon_dims_2d(2) =  'lat';   lon_dims_2d(3) = 'time'
      lat_dims_2d(1) =  'lon';   lat_dims_2d(2) = 'ilat';   lat_dims_2d(3) = 'time'
      vtx_dims_2d(1) = 'ilon';   vtx_dims_2d(2) = 'ilat';   vtx_dims_2d(3) = 'time'

    call fiona_init(time_units, start_time_str)

    do im = 1 , member_total

      write(fid0, '("h0_", i6.6)') im + group_id * member_total
      write(fid1, '("h1_", i6.6)') im + group_id * member_total

      call fiona_create_dataset(fid0, desc=case_desc, file_prefix=trim(case_name), mpi_comm=proc%comm)
      call fiona_add_att(fid0, 'time_step_size', dt)
      call fiona_add_dim(fid0, 'time' , add_var=.true.)
      call fiona_add_dim(fid0, 'lon'  , size=global_mesh%num_full_lon, add_var=.true., decomp=.true.)
      call fiona_add_dim(fid0, 'lat'  , size=global_mesh%num_full_lat, add_var=.true., decomp=.true.)
      call fiona_add_dim(fid0, 'ilon' , size=global_mesh%num_half_lon, add_var=.true., decomp=.true.)
      call fiona_add_dim(fid0, 'ilat' , size=global_mesh%num_half_lat, add_var=.true., decomp=.true.)
      ! if (.not. baroclinic .or. .not. diag_state(1)%is_init()) then
      !   call fiona_add_dim(fid0, 'ilon' , size=global_mesh%num_half_lon, add_var=.true., decomp=.true.)
      !   call fiona_add_dim(fid0, 'ilat' , size=global_mesh%num_half_lat, add_var=.true., decomp=.true.)
      ! end if
      ! Common or barotropic variables
      if (baroclinic) then
        call fiona_add_dim(fid0, 'lev'  , size=global_mesh%num_full_lev, add_var=.true., decomp=.false.)
        call fiona_add_dim(fid0, 'ilev' , size=global_mesh%num_half_lev, add_var=.true., decomp=.false.)
        call fiona_add_var(fid0, 't'    , long_name='temperature'                 , units='K'      , dim_names=cell_dims, data_type='real(8)')
        call fiona_add_var(fid0, 'pt'   , long_name='potential temperature'       , units='K'      , dim_names=cell_dims, data_type='real(8)')
        call fiona_add_var(fid0, 'phs'  , long_name='surface hydrostatic pressure', units='Pa'     , dim_names=cell_dims_2d, data_type='real(8)')
        call fiona_add_var(fid0, 'ph'   , long_name='hydrostatic pressure'        , units='Pa'     , dim_names=lev_dims, data_type='real(8)')
        call fiona_add_var(fid0, 'u'    , long_name='u wind component'            , units='m s-1'  , dim_names=lon_dims, data_type='real(8)')
        call fiona_add_var(fid0, 'v'    , long_name='v wind component'            , units='m s-1'  , dim_names=lat_dims, data_type='real(8)')
        if (output_Agrid) then
          call fiona_add_var(fid0, 'u_Agrid', long_name='u_Agrid'                   , units='K'      , dim_names=cell_dims, data_type='real(8)')
          call fiona_add_var(fid0, 'v_Agrid', long_name='v_Agrid'                   , units='K'      , dim_names=cell_dims, data_type='real(8)')
        end if
        !call fiona_add_var(fid0, 'z'    , long_name='height'                      , units='m'      , dim_names=cell_lev_dims, data_type='real(8)')
        call fiona_add_var(fid0, 'zs'   , long_name='surface height'              , units='m'      , dim_names=cell_dims_2d, data_type='real(8)')
        call fiona_add_var(fid0, 'pv'   , long_name='potential vorticity'         , units='m-1 s-1', dim_names=vtx_dims, data_type='real(8)')
        call fiona_add_var(fid0, 'vor'  , long_name='relative vorticity'          , units='s-1'    , dim_names=vtx_dims, data_type='real(8)')
        call fiona_add_var(fid0, 'div'  , long_name='divergence'                  , units='s-1'    , dim_names=cell_dims, data_type='real(8)')
        call fiona_add_var(fid0, 'tm'   , long_name='total mass'                  , units='m'      , dim_names=['time'], data_type='real(8)')
        call fiona_add_var(fid0, 'te'   , long_name='total energy'                , units='m4 s-4' , dim_names=['time'], data_type='real(8)')
        call fiona_add_var(fid0, 'tpe'  , long_name='total potential enstrophy'   , units='m2 s-5' , dim_names=['time'], data_type='real(8)')
        if (nonhydrostatic) then
          call fiona_add_var(fid0, 'w', long_name='vertical speed', units='m s-1', dim_names=lev_dims, data_type='real(8)')
          call fiona_add_var(fid0, 'p', long_name='full pressure' , units='Pa'   , dim_names=lev_dims, data_type='real(8)')
          call fiona_add_var(fid0, 'rhod', long_name='dry air density', units='', dim_names=cell_dims, data_type='real(8)')
        end if
      else
        call fiona_add_var(fid0, 'u'    , long_name='u wind component'            , units='m s-1'  , dim_names=lon_dims_2d, data_type='real(8)')
        call fiona_add_var(fid0, 'v'    , long_name='v wind component'            , units='m s-1'  , dim_names=lat_dims_2d, data_type='real(8)')
        !call fiona_add_var(fid0, 'z'    , long_name='height'                      , units='m'      , dim_names=cell_dims_2d, data_type='real(8)'))
        call fiona_add_var(fid0, 'zs'   , long_name='surface height'              , units='m'      , dim_names=cell_dims_2d, data_type='real(8)')
        call fiona_add_var(fid0, 'pv'   , long_name='potential vorticity'         , units='m-1 s-1', dim_names=vtx_dims_2d, data_type='real(8)')
        call fiona_add_var(fid0, 'vor'  , long_name='relative vorticity'          , units='s-1'    , dim_names=vtx_dims_2d, data_type='real(8)')
        call fiona_add_var(fid0, 'div'  , long_name='divergence'                  , units='s-1'    , dim_names=cell_dims_2d, data_type='real(8)')
        call fiona_add_var(fid0, 'tm'   , long_name='total mass'                  , units='m'      , dim_names=['time'], data_type='real(8)')
        call fiona_add_var(fid0, 'te'   , long_name='total energy'                , units='m4 s-4' , dim_names=['time'], data_type='real(8)')
        call fiona_add_var(fid0, 'tpe'  , long_name='total potential enstrophy'   , units='m2 s-5' , dim_names=['time'], data_type='real(8)')
      end if

      call fiona_create_dataset(fid1, desc=case_desc, file_prefix=trim(case_name), mpi_comm=proc%comm)
      call fiona_add_att(fid1, 'time_step_size', dt)
      call fiona_add_dim(fid1, 'time' , add_var=.true.)
      call fiona_add_dim(fid1, 'lon'  , size=global_mesh%num_full_lon, add_var=.true., decomp=.true.)
      call fiona_add_dim(fid1, 'lat'  , size=global_mesh%num_full_lat, add_var=.true., decomp=.true.)
      call fiona_add_dim(fid1, 'ilon' , size=global_mesh%num_half_lon, add_var=.true., decomp=.true.)
      call fiona_add_dim(fid1, 'ilat' , size=global_mesh%num_half_lat, add_var=.true., decomp=.true.)
      if (baroclinic) then
        call fiona_add_dim(fid1, 'lev'  , size=global_mesh%num_full_lev, add_var=.true., decomp=.false.)
        call fiona_add_dim(fid1, 'ilev' , size=global_mesh%num_half_lev, add_var=.true., decomp=.false.)
        call fiona_add_var(fid1, 'dudt'         , long_name='u wind component tendency'                     , units='', dim_names=lon_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'dvdt'         , long_name='v wind component tendency'                     , units='', dim_names=lat_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'dphsdt'       , long_name='surface hydrostatic pressure tendency'         , units='', dim_names=cell_dims_2d, data_type='real(8)')
        call fiona_add_var(fid1, 'dptdt'        , long_name='potential temperature tendency'                , units='', dim_names=cell_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'dptfdlon'     , long_name='zonal potential temperature flux gradient'     , units='', dim_names=cell_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'dptfdlat'     , long_name='meridional potential temperature flux gradient', units='', dim_names=cell_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'dptfdlev'     , long_name='vertical potential temperature flux gradient'  , units='', dim_names=cell_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'wedphdlev_lev', long_name='vertical coordinate velocity'                  , units='', dim_names=lev_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'pgf_lon'      , long_name='zonal pressure gradient force'                 , units='', dim_names=lon_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'wedudlev'     , long_name='vertical advection of u'                       , units='', dim_names=lon_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'pgf_lat'      , long_name='meridional pressure gradient force'            , units='', dim_names=lat_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'wedvdlev'     , long_name='vertical advection of v'                       , units='', dim_names=lat_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'qhv'          , long_name='nonliear zonal Coriolis force'                 , units='', dim_names=lon_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'qhu'          , long_name='nonliear meridional Coriolis force'            , units='', dim_names=lat_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'dkedlon'      , long_name='zonal kinetic energy gradient force'           , units='', dim_names=lon_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'dkedlat'      , long_name='meridional kinetic energy gradient force'      , units='', dim_names=lat_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'dmfdlon'      , long_name='zonal mass flux divergence'                    , units='', dim_names=cell_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'dmfdlat'      , long_name='meridional mass flux divergence'               , units='', dim_names=cell_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'mf_lon_n'     , long_name='normal mass flux on U grid'                    , units='', dim_names=lon_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'mf_lon_t'     , long_name='tangent mass flux on U grid'                   , units='', dim_names=lon_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'mf_lat_n'     , long_name='normal mass flux on V grid'                    , units='', dim_names=lat_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'mf_lat_t'     , long_name='tangent mass flux on V grid'                   , units='', dim_names=lat_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'm'            , long_name='dph on full levels'                            , units='', dim_names=cell_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'm_lon'        , long_name='dph on lon edges on full levels'               , units='', dim_names=lon_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'pt_lon'       , long_name='pt on lon edges on full levels'                , units='', dim_names=lon_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'pt_lat'       , long_name='pt on lat edges on full levels'                , units='', dim_names=lat_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'm_lat'        , long_name='dph on lat edges on full levels'               , units='', dim_names=lat_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'm_vtx'        , long_name='dph on vertices on full levels'                , units='', dim_names=vtx_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'pv_lon'       , long_name='pv on U grid'                                  , units='', dim_names=lon_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'pv_lat'       , long_name='pv on V grid'                                  , units='', dim_names=lat_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'ke'           , long_name='kinetic energy on cell grid'                   , units='', dim_names=cell_dims, data_type='real(8)')
        if (nonhydrostatic) then
          call fiona_add_var(fid1, 'adv_gz_lon', long_name='advection of gz', units='', dim_names=lev_dims, data_type='real(8)')
          call fiona_add_var(fid1, 'adv_gz_lat', long_name='advection of gz', units='', dim_names=lev_dims, data_type='real(8)')
          call fiona_add_var(fid1, 'adv_gz_lev', long_name='advection of gz', units='', dim_names=lev_dims, data_type='real(8)')
          call fiona_add_var(fid1, 'adv_w_lon' , long_name='advection of w' , units='', dim_names=lev_dims, data_type='real(8)')
          call fiona_add_var(fid1, 'adv_w_lat' , long_name='advection of w' , units='', dim_names=lev_dims, data_type='real(8)')
          call fiona_add_var(fid1, 'adv_w_lev' , long_name='advection of w' , units='', dim_names=lev_dims, data_type='real(8)')
          call fiona_add_var(fid1, 'dzsdlon', long_name='zonal zs gradient', units='', dim_names=cell_dims_2d, data_type='real(8)')
          call fiona_add_var(fid1, 'dzsdlat', long_name='meridional zs gradient', units='', dim_names=cell_dims_2d, data_type='real(8)')
        end if   
      else
        call fiona_add_var(fid1, 'dudt'         , long_name='u wind component tendency'                     , units='', dim_names=lon_dims_2d, data_type='real(8)')
        call fiona_add_var(fid1, 'dvdt'         , long_name='v wind component tendency'                     , units='', dim_names=lat_dims_2d, data_type='real(8)')
        call fiona_add_var(fid1, 'dgzdt'        , long_name='geopotential tendency'                         , units='', dim_names=cell_dims_2d, data_type='real(8)')
        call fiona_add_var(fid1, 'qhv'          , long_name='nonliear zonal Coriolis force'                 , units='', dim_names=lon_dims_2d, data_type='real(8)')
        call fiona_add_var(fid1, 'qhu'          , long_name='nonliear meridional Coriolis force'            , units='', dim_names=lat_dims_2d, data_type='real(8)')
        call fiona_add_var(fid1, 'pgf_lon'      , long_name='zonal geopotential energy gradient force'      , units='', dim_names=lon_dims_2d, data_type='real(8)')
        call fiona_add_var(fid1, 'dkedlon'      , long_name='zonal kinetic energy gradient force'           , units='', dim_names=lon_dims_2d, data_type='real(8)')
        call fiona_add_var(fid1, 'pgf_lat'      , long_name='meridional geopotential energy gradient force' , units='', dim_names=lat_dims_2d, data_type='real(8)')
        call fiona_add_var(fid1, 'dkedlat'      , long_name='meridional kinetic energy gradient force'      , units='', dim_names=lat_dims_2d, data_type='real(8)')
        call fiona_add_var(fid1, 'dmfdlon'      , long_name='zonal mass flux divergence'                    , units='', dim_names=cell_dims_2d, data_type='real(8)')
        call fiona_add_var(fid1, 'dmfdlat'      , long_name='meridional mass flux divergence'               , units='', dim_names=cell_dims_2d, data_type='real(8)')
        call fiona_add_var(fid1, 'mf_lon_n'     , long_name='normal mass flux on U grid'                    , units='', dim_names=lon_dims_2d, data_type='real(8)')
        call fiona_add_var(fid1, 'mf_lon_t'     , long_name='tangent mass flux on U grid'                   , units='', dim_names=lon_dims_2d, data_type='real(8)')
        call fiona_add_var(fid1, 'mf_lat_n'     , long_name='normal mass flux on V grid'                    , units='', dim_names=lat_dims_2d, data_type='real(8)')
        call fiona_add_var(fid1, 'mf_lat_t'     , long_name='tangent mass flux on V grid'                   , units='', dim_names=lat_dims_2d, data_type='real(8)')
        call fiona_add_var(fid1, 'm_lon'        , long_name='dph on lon edges on full levels'               , units='', dim_names=lon_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'm_lat'        , long_name='dph on lat edges on full levels'               , units='', dim_names=lat_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'pt_lat'       , long_name='pt on lat edges on full levels'                , units='', dim_names=lat_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'm_vtx'        , long_name='dph on vertices on full levels'                , units='', dim_names=vtx_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'pv_lon'       , long_name='pv on U grid'                                  , units='', dim_names=lon_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'pv_lat'       , long_name='pv on V grid'                                  , units='', dim_names=lat_dims, data_type='real(8)')
        call fiona_add_var(fid1, 'ke'           , long_name='kinetic energy on cell grid'                   , units='', dim_names=cell_dims_2d, data_type='real(8)')
      end if
    end do

  end subroutine history_init

  subroutine history_final()

#ifdef Detail_io_Time 
    close((23333+myid))
#endif

  end subroutine history_final

  subroutine history_write_state(blocks, itime, group_id_in , new_file)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime
    !integer, intent(in) :: member_num
    integer , intent(in) , optional :: group_id_in
    logical , intent(in) , optional :: new_file

    type(mesh_type), pointer :: mesh
    type(state_type), pointer :: state
    type(static_type), pointer :: static
    integer iblk, is, ie, js, je, ks, ke
    integer i,j,k;
    integer start(3), count(3)

    integer im
    character*20 fid0 , fid1
    integer group_id 

    if (present(group_id_in)) then
      group_id = group_id_in
    else
      group_id = 0
    endif

    do im = 1 , member_num

      if ( (ivector - 1) * member_num + im > member_total ) exit

      write(fid0, '("h0_", i6.6)') (ivector - 1) * member_num + im + group_id * member_total

      call fiona_start_output(fid0, elapsed_seconds, new_file = new_file)
      call fiona_output(fid0, 'lon' , global_mesh%full_lon_deg(1:global_mesh%num_full_lon))
      call fiona_output(fid0, 'lat' , global_mesh%full_lat_deg(1:global_mesh%num_full_lat))
      call fiona_output(fid0, 'ilon', global_mesh%half_lon_deg(1:global_mesh%num_half_lon))
      call fiona_output(fid0, 'ilat', global_mesh%half_lat_deg(1:global_mesh%num_half_lat))
      if (baroclinic) then
        call fiona_output(fid0, 'lev' , global_mesh%full_lev)
        call fiona_output(fid0, 'ilev', global_mesh%half_lev)
      end if

      do iblk = 1, size(blocks)
        mesh => blocks(iblk)%mesh
        state => blocks(iblk)%state(itime,ivector)
        static => blocks(iblk)%static(ivector)

        if (output_Agrid) then
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do j = mesh%full_lat_ibeg, mesh%full_lat_iend
              do i = mesh%full_lon_ibeg, mesh%full_lon_iend
                state%u_Agrid(im,i,j,k) = (state%u(im,i,j,k) + state%u(im,i-1,j,k)) / 2;
                state%v_Agrid(im,i,j,k) = 2 * (state%v(im,i,j-1,k) * mesh%area_subcell(1,j) + state%v(im,i,j,k) * mesh%area_subcell(2,j)) / mesh%area_cell(j) 
              end do
            end do
          end do
        end if  


        is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
        js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
        ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
        start = [is,js,ks]
        count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]

        !call fiona_output(fid0, 'zs' , static%gzs(im,is:ie,js:je      ) / g, start=start, count=count)
        call fiona_output(fid0, 'div', state%div (im,is:ie,js:je,ks:ke)    , start=start, count=count)

        if (baroclinic) then
          call fiona_output(fid0, 't'     , state%t     (im,is:ie,js:je,ks:ke), start=start, count=count)
          call fiona_output(fid0, 'pt'    , state%pt_f    (im,is:ie,js:je,ks:ke), start=start, count=count)
          call fiona_output(fid0, 'phs'   , state%phs_f   (im,is:ie,js:je      ), start=start, count=count)
         
          if (output_Agrid) then
            call fiona_output(fid0, 'u_Agrid'    , state%u_Agrid    (im,is:ie,js:je,ks:ke), start=start, count=count)
            call fiona_output(fid0, 'v_Agrid'    , state%v_Agrid    (im,is:ie,js:je,ks:ke), start=start, count=count)
          end if
        end if
        
        

        is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
        js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
        ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
        start = [is,js,ks]
        count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_full_lev]

        call fiona_output(fid0, 'u'   , state%u_f   (im,is:ie,js:je,ks:ke), start=start, count=count)

        is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
        js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
        ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
        start = [is,js,ks]
        count = [mesh%num_full_lon,mesh%num_half_lat,mesh%num_full_lev]

        call fiona_output(fid0, 'v'   , state%v_f   (im,is:ie,js:je,ks:ke), start=start, count=count)

        is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
        js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
        ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
        start = [is,js,ks]
        count = [mesh%num_half_lon,mesh%num_half_lat,mesh%num_full_lev]

        call fiona_output(fid0, 'pv' , state %pv (im,is:ie,js:je,ks:ke), start=start, count=count)
        call fiona_output(fid0, 'vor', state %vor(im,is:ie,js:je,ks:ke), start=start, count=count)

        is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
        js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
        ks = mesh%half_lev_ibeg; ke = mesh%half_lev_iend
        start = [is,js,ks]
        count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_half_lev]

        call fiona_output(fid0, 'ph'    , state%ph_lev    (im,is:ie,js:je,ks:ke), start=start, count=count)

        !call fiona_output(fid0, 'z'  , state%gz_lev(im,is:ie,js:je,ks:ke) / g, start=start, count=count)
        if (nonhydrostatic) then
          call fiona_output(fid0, 'w', state%w_lev(im,is:ie,js:je,ks:ke), start=start, count=count)
          call fiona_output(fid0, 'p', state%p_lev(im,is:ie,js:je,ks:ke), start=start, count=count)
          call fiona_output(fid0, 'rhod', state%rhod(im,is:ie,js:je,ks:ke), start=start, count=count)
        end if

        call fiona_output(fid0, 'tm' , state %tm(im))
        call fiona_output(fid0, 'te' , state %te(im))
        call fiona_output(fid0, 'tpe', state %tpe(im))
      end do
      call fiona_end_output(fid0)
    end do

  end subroutine history_write_state

  subroutine history_write_debug(blocks, itime , group_id_in, new_file)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime
    !integer, intent(in) :: member_num
    integer , intent(in) , optional :: group_id_in
    logical , intent(in) , optional :: new_file


    type(mesh_type), pointer :: mesh
    type(state_type), pointer :: state
    type(tend_type), pointer :: tend

    integer is, ie, js, je, ks, ke
    integer start(3), count(3)

    integer im
    character*20 fid0 , fid1
    integer group_id

    if (present(group_id_in)) then
      group_id = group_id_in
    else
      group_id = 0
    endif

    mesh => blocks(1)%mesh
    state => blocks(1)%state(itime,ivector)
    tend => blocks(1)%tend(itime,ivector)

    do im = 1 , member_num

      if ( (ivector - 1) * member_num + im > member_total ) exit

      write(fid1, '("h1_", i6.6)') (ivector - 1) * member_num + im + group_id * member_total

      call fiona_start_output(fid1, elapsed_seconds, new_file = new_file)
      
      call fiona_output(fid1, 'lon'   , global_mesh%full_lon_deg(1:global_mesh%num_full_lon))
      call fiona_output(fid1, 'lat'   , global_mesh%full_lat_deg(1:global_mesh%num_full_lat))
      call fiona_output(fid1, 'ilon'  , global_mesh%half_lon_deg(1:global_mesh%num_half_lon))
      call fiona_output(fid1, 'ilat'  , global_mesh%half_lat_deg(1:global_mesh%num_half_lat))
      if (baroclinic) then
        call fiona_output(fid1, 'lev' , global_mesh%full_lev)
        call fiona_output(fid1, 'ilev', global_mesh%half_lev)
      end if

      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]

      call fiona_output(fid1, 'dmfdlon' , tend%dmfdlon  (im,is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output(fid1, 'dmfdlat' , tend%dmfdlat  (im,is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output(fid1, 'ke'      , state%ke      (im,is:ie,js:je,ks:ke), start=start, count=count)
      if (baroclinic) then
        call fiona_output(fid1, 'm'       , state%m      (im,is:ie,js:je,ks:ke), start=start, count=count)
        call fiona_output(fid1, 'dphsdt'  , tend%dphs    (im,is:ie,js:je      ), start=start, count=count)
        call fiona_output(fid1, 'dptdt'   , tend%dpt     (im,is:ie,js:je,ks:ke), start=start, count=count)
        call fiona_output(fid1, 'dptfdlon', tend%dptfdlon(im,is:ie,js:je,ks:ke), start=start, count=count)
        call fiona_output(fid1, 'dptfdlat', tend%dptfdlat(im,is:ie,js:je,ks:ke), start=start, count=count)
        call fiona_output(fid1, 'dptfdlev', tend%dptfdlev(im,is:ie,js:je,ks:ke), start=start, count=count)
      else
        call fiona_output(fid1, 'dgzdt'   , tend%dgz     (im,is:ie,js:je,ks:ke), start=start, count=count)
      end if

      is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_full_lev]

      call fiona_output(fid1, 'qhv'     , tend%qhv      (im,is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output(fid1, 'pgf_lon' , tend%pgf_lon  (im,is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output(fid1, 'dkedlon' , tend%dkedlon  (im,is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output(fid1, 'dudt   ' , tend%du       (im,is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output(fid1, 'm_lon'   , state%m_lon   (im,is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output(fid1, 'pt_lon'  , state%pt_lon  (im,is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output(fid1, 'mf_lon_n', state%mf_lon_n(im,is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output(fid1, 'mf_lon_t', state%mf_lon_t(im,is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output(fid1, 'pv_lon'  , state%pv_lon  (im,is:ie,js:je,ks:ke), start=start, count=count)

      if (baroclinic) then
        call fiona_output(fid1, 'wedudlev', tend%wedudlev(im,is:ie,js:je,ks:ke), start=start, count=count)
      end if

      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_half_lat,mesh%num_full_lev]

      call fiona_output(fid1, 'qhu'     , tend%qhu      (im,is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output(fid1, 'pgf_lat' , tend%pgf_lat  (im,is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output(fid1, 'dkedlat' , tend%dkedlat  (im,is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output(fid1, 'dvdt'    , tend%dv       (im,is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output(fid1, 'm_lat'   , state%m_lat   (im,is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output(fid1, 'pt_lat'  , state%pt_lat  (im,is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output(fid1, 'mf_lat_n', state%mf_lat_n(im,is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output(fid1, 'mf_lat_t', state%mf_lat_t(im,is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output(fid1, 'pv_lat'  , state%pv_lat  (im,is:ie,js:je,ks:ke), start=start, count=count)

      if (baroclinic) then
        call fiona_output(fid1, 'wedvdlev', tend%wedvdlev(im,is:ie,js:je,ks:ke), start=start, count=count)
      end if

      is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
      js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_half_lon,mesh%num_half_lat,mesh%num_full_lev]

      call fiona_output(fid1, 'm_vtx', state%m_vtx(im,is:ie,js:je,ks:ke), start=start, count=count)

      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%half_lev_ibeg; ke = mesh%half_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_half_lev]
  
      if (baroclinic) then
        call fiona_output(fid1, 'wedphdlev_lev', state%wedphdlev_lev(im,is:ie,js:je,ks:ke), start=start, count=count)
      end if

      if (nonhydrostatic) then
        call fiona_output(fid1, 'adv_gz_lon', tend%adv_gz_lon(im,is:ie,js:je,ks:ke), start=start, count=count)
        call fiona_output(fid1, 'adv_gz_lat', tend%adv_gz_lat(im,is:ie,js:je,ks:ke), start=start, count=count)
        call fiona_output(fid1, 'adv_gz_lev', tend%adv_gz_lev(im,is:ie,js:je,ks:ke), start=start, count=count)
        call fiona_output(fid1, 'adv_w_lon' , tend%adv_w_lon (im,is:ie,js:je,ks:ke), start=start, count=count)
        call fiona_output(fid1, 'adv_w_lat' , tend%adv_w_lat (im,is:ie,js:je,ks:ke), start=start, count=count)
        call fiona_output(fid1, 'adv_w_lev' , tend%adv_w_lev (im,is:ie,js:je,ks:ke), start=start, count=count)
        call fiona_output(fid1, 'dzsdlon', blocks(1)%static(ivector)%dzsdlon(im,is:ie,js:je), start=start, count=count)
        call fiona_output(fid1, 'dzsdlat', blocks(1)%static(ivector)%dzsdlat(im,is:ie,js:je), start=start, count=count)
      end if

      call fiona_end_output(fid1)
    end do

  end subroutine history_write_debug

end module history_mod
