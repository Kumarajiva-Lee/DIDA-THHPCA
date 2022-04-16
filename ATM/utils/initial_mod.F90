module initial_mod

    use fiona
    use string
    use flogger
    use datetime
    use const_mod
    use namelist_mod
    use time_mod
    use block_mod
    use parallel_mod
    use parallel_types_mod
    use process_mod
    use member_mod

    implicit none
  
    private

    public initial_write
    public initial_read
    public date_advance
  
  contains

    subroutine date_advance(year , month , day , hour , minute , advance_time)

      integer, intent(inout) :: year
      integer, intent(inout) :: month
      integer, intent(inout) :: day
      integer, intent(inout) :: hour
      integer, intent(inout) :: minute
      integer, intent(in)    :: advance_time

      integer month_day(12) 
      
      data month_day /31,28,31,30,31,30,31,31,30,31,30,31/

      if (mod(year , 400) == 0) then 
        month_day(2) = 29
      else if (mod(year , 100) /= 0 .and. mod(year , 4) ==0 ) then
        month_day(2) = 29
      end if

      minute = minute + advance_time 
      if (minute > 59) then
        hour = hour + ((minute - mod(minute , 60)) / 60)
        minute = mod(minute , 60)
        if (hour > 23) then
          day = day + ((hour - mod(hour , 24)) / 24)
          hour = mod(hour , 24)
          if (day > month_day(month)) then
            day = 1
            month = month + 1
            if (month > 12) then
              month = 1
              year = year + 1
            end if
          end if
        end if
      end if
    end subroutine date_advance
  
    subroutine initial_write(initial_file_)
      
      character(*), intent(in) :: initial_file_


      character(4) cell_dims(4), cell_dims_2d(3)
      character(4) lon_dims(4)
      character(4) lat_dims(4)
      integer iblk, is, ie, js, je, ks, ke
      integer start(3), count(3)

      integer im , length
      character*10 s_id
      character*256 filename

      filename = initial_file_

        cell_dims(1) =  'lon';     cell_dims(2) =  'lat';     cell_dims(3) =  'lev';     cell_dims(4) = 'time'
          lon_dims(1) = 'ilon';      lon_dims(2) =  'lat';      lon_dims(3) =  'lev';      lon_dims(4) = 'time'
          lat_dims(1) =  'lon';      lat_dims(2) = 'ilat';      lat_dims(3) =  'lev';      lat_dims(4) = 'time'
      cell_dims_2d(1) =  'lon';  cell_dims_2d(2) =  'lat';  cell_dims_2d(3) = 'time'
      
      im = 1

      call fiona_create_dataset('i0', file_path=filename, start_time='1970-01-01', time_units='hours', mpi_comm=proc%comm)
      call fiona_add_dim('i0', 'time', add_var=.true.)
      call fiona_add_dim('i0', 'lon'  , size=global_mesh%num_full_lon, add_var=.true., decomp=.true.)
      call fiona_add_dim('i0', 'lat'  , size=global_mesh%num_full_lat, add_var=.true., decomp=.true.)
      call fiona_add_dim('i0', 'ilon' , size=global_mesh%num_half_lon, add_var=.true., decomp=.true.)
      call fiona_add_dim('i0', 'ilat' , size=global_mesh%num_half_lat, add_var=.true., decomp=.true.)
      if (baroclinic) then
        call fiona_add_dim('i0', 'lev'  , size=global_mesh%num_full_lev, add_var=.true., decomp=.false.)
        call fiona_add_dim('i0', 'ilev' , size=global_mesh%num_half_lev, add_var=.true., decomp=.false.)
        call fiona_add_var('i0', 'pt'   , long_name='potential temperature'       , units='K'      , dim_names=cell_dims)
        call fiona_add_var('i0', 'phs'  , long_name='surface hydrostatic pressure', units='Pa'     , dim_names=cell_dims_2d)
        call fiona_add_var('i0', 'u'    , long_name='u wind component'            , units='m s-1'  , dim_names=lon_dims)
        call fiona_add_var('i0', 'v'    , long_name='v wind component'            , units='m s-1'  , dim_names=lat_dims)
      end if
      call fiona_add_var('i0', 'zs', long_name='surface height', units='m' , dim_names=cell_dims_2d)

      call fiona_start_output('i0', 0.0d0)
      call fiona_output('i0', 'lon' , global_mesh%full_lon_deg(1:global_mesh%num_full_lon))
      call fiona_output('i0', 'lat' , global_mesh%full_lat_deg(1:global_mesh%num_full_lat))
      call fiona_output('i0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%num_half_lon))
      call fiona_output('i0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%num_half_lat))
      if (baroclinic) then
        call fiona_output('i0', 'lev' , global_mesh%full_lev)
        call fiona_output('i0', 'ilev', global_mesh%half_lev)
      end if

      do iblk = 1, size(proc%blocks)
        associate (mesh   => proc%blocks(iblk)%mesh,     &
                  state  => proc%blocks(iblk)%state(1,1), &
                  static => proc%blocks(iblk)%static(1))
          is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
          js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
          ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
          start = [is,js,ks]
          count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]

          if (baroclinic) then
            call fiona_output('i0', 'pt' , state%pt (im,is:ie,js:je,ks:ke), start=start, count=count)
            call fiona_output('i0', 'phs', state%phs(im,is:ie,js:je      ), start=start, count=count)
          end if
          call fiona_output('i0', 'zs' , static%gzs(im,is:ie,js:je) / g, start=start, count=count)

          is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
          js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
          ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
          start = [is,js,ks]
          count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_full_lev]

          call fiona_output('i0', 'u', state%u(im,is:ie,js:je,ks:ke), start=start, count=count)

          is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
          js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
          ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
          start = [is,js,ks]
          count = [mesh%num_full_lon,mesh%num_half_lat,mesh%num_full_lev]

          call fiona_output('i0', 'v', state%v(im,is:ie,js:je,ks:ke), start=start, count=count)
        end associate
      end do

      call fiona_end_output('i0')



    end subroutine initial_write


    subroutine initial_read(initial_file_ , group_id_in)
  
      !integer , intent(in) :: member
      character(*), intent(in)  :: initial_file_
      integer , intent(in) , optional :: group_id_in

  
      integer iblk, is, ie, js, je, ks, ke
      integer start(3), count(3)
      integer im , length
      character*10 s_id
      character*256 filename
      integer ini_year , ini_month , ini_day ,ini_hour , ini_minute
      integer group_id

      if (present(group_id_in)) then
        group_id = group_id_in
      else
        group_id = 0
      endif
  
      filename = initial_file_
      length = len_trim(initial_file_)

      if (initial_file_type == 'N/A') then
        if (is_root_proc()) call log_error('Missing initial_file_type!')
      else if (initial_file_type == 'time' .and. initial_interval == -1) then
        if (is_root_proc()) call log_error('Missing initial_interval!')
      end if

      if (initial_file_type == 'time') then
        read(initial_file_(length - 21 : length - 18) , *) ini_year
        read(initial_file_(length - 16 : length - 15) , *) ini_month
        read(initial_file_(length - 13 : length - 12) , *) ini_day
        read(initial_file_(length - 10 : length - 9)  , *) ini_hour
        read(initial_file_(length - 7  : length - 6)  , *) ini_minute
        
        call date_advance(ini_year , ini_month , ini_day , ini_hour , ini_minute , initial_interval * group_id * member_num)

        write(s_id,"(i4.4)") ini_year
        filename(length - 21 : length - 18) = s_id
        s_id = 'N/A'
    
        write(s_id,"(i2.2)") ini_month
        filename(length - 16 : length - 15) = s_id
        s_id = 'N/A'
        
        write(s_id,"(i2.2)") ini_day
        filename(length - 13 : length - 12) = s_id
        s_id = 'N/A'
    
        write(s_id,"(i2.2)") ini_hour
        filename(length - 10 : length - 9)  = s_id
        s_id = 'N/A'

        write(s_id,"(i2.2)") ini_minute
        filename(length - 7  : length - 6)  = s_id
        s_id = 'N/A'

      end if

      do im = 1 , member_num
        
        if ( (ivector - 1) * member_num + im > member_total ) exit

        if (initial_file_type == 'number') then
          write(s_id,"(i6.6)") (ivector - 1) * member_num + im + group_id * member_total
          filename(length-4 : length - 3) = s_id(1:2)
        end if

        if (initial_file_ /= 'N/A') then
          call fiona_open_dataset('i0', file_path=filename, mpi_comm=proc%comm)
          if (is_root_proc()) call log_notice('Read initial data from ' // trim(filename) // '.')
        else
          call fiona_open_dataset('i0', file_path=filename, mpi_comm=proc%comm)
          if (is_root_proc()) call log_notice('Read initial data from ' // trim(filename) // '.')
        end if
        call fiona_start_input('i0')
    
        do iblk = 1, size(proc%blocks)
          associate (block  => proc%blocks(iblk)                    , &
                    mesh   => proc%blocks(iblk)%mesh               , &
                    state  => proc%blocks(iblk)%state(old_time_idx,ivector), &
                    static => proc%blocks(iblk)%static(ivector))
            is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
            js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
            ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
            start = [is,js,ks]
            count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]
    
            call fiona_input('i0', 'zs', static%gzs(im , is:ie,js:je), start=start, count=count)
            static%gzs(im,:,:) = static%gzs(im,:,:) * g
            call  fill_halo_member(block, static%gzs, full_lon=.true., full_lat=.true., async=state%async(async_gzs))
            if (baroclinic) then
              call fiona_input('i0', 'phs', state%phs(im , is:ie,js:je      ), start=start, count=count)
              call  fill_halo_member(block, state%phs, full_lon=.true., full_lat=.true., async=state%async(async_phs))
              call fiona_input('i0', 'pt' , state%pt (im , is:ie,js:je,ks:ke), start=start, count=count)
              call  fill_halo_member(block, state%pt, full_lon=.true., full_lat=.true., full_lev=.true., async=state%async(async_pt))
            else
              call fiona_input('i0', 'z' , state%gz (im , is:ie,js:je,ks:ke), start=start, count=count)
              state%gz(im,:,:,:) = state%gz(im,:,:,:) * g
              call  fill_halo_member(block, state%gz, full_lon=.true., full_lat=.true., full_lev=.true., async=state%async(async_gz))
            end if
    
            is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
            js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
            ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
            start = [is,js,ks]
            count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_full_lev]
    
            call fiona_input('i0', 'u'  , state%u(im , is:ie,js:je,ks:ke), start=start, count=count)
            call  fill_halo_member(block, state%u, full_lon=.false., full_lat=.true., full_lev=.true., async=state%async(async_u))
    
            is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
            js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
            ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
            start = [is,js,ks]
            count = [mesh%num_full_lon,mesh%num_half_lat,mesh%num_full_lev]
    
            call fiona_input('i0', 'v'  , state%v(im , is:ie,js:je,ks:ke), start=start, count=count)
            call fill_halo_member(block, state%v, full_lon=.true., full_lat=.false., full_lev=.true., async=state%async(async_v))

            call state%async(async_gzs)%wait()
            call state%async(async_phs)%wait()
            call state%async(async_pt)%wait()
            call state%async(async_gz)%wait()
            call state%async(async_u)%wait()
            call state%async(async_v)%wait()

          end associate
        end do
        call fiona_end_input('i0')

        if (initial_file_type == 'time') then
          call date_advance(ini_year , ini_month , ini_day , ini_hour , ini_minute , initial_interval)

          write(s_id,"(i4.4)") ini_year
          filename(length - 21 : length - 18) = s_id
          s_id = 'N/A'
      
          write(s_id,"(i2.2)") ini_month
          filename(length - 16 : length - 15) = s_id
          s_id = 'N/A'
          
          write(s_id,"(i2.2)") ini_day
          filename(length - 13 : length - 12) = s_id
          s_id = 'N/A'
      
          write(s_id,"(i2.2)") ini_hour
          filename(length - 10 : length - 9)  = s_id
          s_id = 'N/A'

          write(s_id,"(i2.2)") ini_minute
          filename(length - 7  : length - 6)  = s_id
          s_id = 'N/A'

        end if


      end do !end do im
      if (is_root_proc()) call log_print('Finish initial_read')  
    end subroutine initial_read
  
  end module initial_mod