module filter_mod

    use const_mod
    use block_mod
    use parallel_mod
    use filter_types_mod
    use namelist_mod
    use pa_mod
  
    implicit none
  
    private
  
    public filter_type
    public filter_on_cell
    public filter_on_vtx
    public filter_on_lon_edge
    public filter_on_lat_edge
    public filter_on_lev_edge
    public filter_vector_on_cell
    public filter_vector_on_edges
  
    interface filter_on_cell
      module procedure filter_on_cell_2d
      module procedure filter_on_cell_3d
    end interface filter_on_cell
  
  contains
  
    subroutine filter_on_cell_2d(block, x, y)
  
      type(block_type), intent(in) :: block
      real(r8), intent(inout) :: x(member_num ,&
                                   block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub)
      real(r8), intent(out), optional :: y(member_num ,&
                                           block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                           block%mesh%full_lat_lb:block%mesh%full_lat_ub)
  
      real(r8) tmp(member_num ,block%mesh%full_lon_lb:block%mesh%full_lon_ub)
      integer i, j, n, hn, im
  
#ifdef Detail_Time
      call Get_Start_Time(cal_time_start)
#endif
      associate (mesh => block%mesh, filter => block%filter)
      do im = 1 , member_num
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            if (filter%ngrid_lon(j) > 1) then
            n  = filter%ngrid_lon(j)
            hn = (n - 1) / 2
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
                tmp(im,i) = sum(filter%wgt_lon(:n,j) * x(im,i-hn:i+hn,j))
            end do
            if (present(y)) then
                y(im,:,j) = tmp(im,:)
            else
                x(im,:,j) = tmp(im,:)
            end if
            else if (present(y)) then
            y(im,:,j) = x(im,:,j)
            end if
        end do
      end do
      end associate
#ifdef Detail_Time
      call Get_End_Time(cal_time_end)
      cal_time = cal_time + cal_time_end - cal_time_start
#endif
  
    end subroutine filter_on_cell_2d
  
    subroutine filter_on_cell_3d(block, x, y)
  
      type(block_type), intent(in) :: block
      real(r8), intent(inout) :: x(member_num , &
                                   block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out), optional :: y(member_num , &
                                           block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                           block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                           block%mesh%full_lev_lb:block%mesh%full_lev_ub)
  
      real(r8) tmp(member_num,block%mesh%full_lon_lb:block%mesh%full_lon_ub)
      integer i, j, k, n, hn, im
  
#ifdef Detail_Time
      call Get_Start_Time(cal_time_start)
#endif
      associate (mesh => block%mesh, filter => block%filter)
      do im = 1 , member_num
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            if (filter%ngrid_lon(j) > 1) then
                n  = filter%ngrid_lon(j)
                hn = (n - 1) / 2
                do i = mesh%full_lon_ibeg, mesh%full_lon_iend
                tmp(im,i) = sum(filter%wgt_lon(:n,j) * x(im,i-hn:i+hn,j,k))
                end do
                if (present(y)) then
                y(im,:,j,k) = tmp(im,:)
                else
                x(im,:,j,k) = tmp(im,:)
                end if
            else if (present(y)) then
                y(im,:,j,k) = x(im,:,j,k)
            end if
            end do
        end do
      end do
      end associate
#ifdef Detail_Time
      call Get_End_Time(cal_time_end)
      cal_time = cal_time + cal_time_end - cal_time_start
#endif
  
    end subroutine filter_on_cell_3d
  
    subroutine filter_on_vtx(block, x, y)
  
      type(block_type), intent(in) :: block
      real(r8), intent(inout) :: x(member_num ,&
                                   block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out), optional :: y(member_num ,&
                                           block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                           block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                           block%mesh%full_lev_lb:block%mesh%full_lev_ub)
  
      real(r8) tmp(member_num,block%mesh%half_lon_lb:block%mesh%half_lon_ub)
      integer i, j, k, n, hn, im
  
#ifdef Detail_Time
      call Get_Start_Time(cal_time_start)
#endif
      associate (mesh => block%mesh, filter => block%filter)
      do im = 1 , member_num
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do j = mesh%half_lat_ibeg, mesh%half_lat_iend
            if (filter%ngrid_lat(j) > 1) then
                n  = filter%ngrid_lat(j)
                hn = (n - 1) / 2
                do i = mesh%half_lon_ibeg, mesh%half_lon_iend
                tmp(im,i) = sum(filter%wgt_lat(:n,j) * x(im,i-hn:i+hn,j,k))
                end do
                if (present(y)) then
                y(im,:,j,k) = tmp(im,:)
                else
                x(im,:,j,k) = tmp(im,:)
                end if
            else if (present(y)) then
                y(im,:,j,k) = x(im,:,j,k)
            end if
            end do
        end do
      end do
      end associate
#ifdef Detail_Time
      call Get_End_Time(cal_time_end)
      cal_time = cal_time + cal_time_end - cal_time_start
#endif
  
    end subroutine filter_on_vtx
  
    subroutine filter_on_lon_edge(block, x, y)
  
      type(block_type), intent(in) :: block
      real(r8), intent(inout) :: x(member_num ,&
                                   block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out), optional :: y(member_num ,&
                                           block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                           block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                           block%mesh%full_lev_lb:block%mesh%full_lev_ub)
  
      real(r8) tmp(member_num ,block%mesh%half_lon_lb:block%mesh%half_lon_ub)
      integer i, j, k, n, hn, im
  
      associate (mesh => block%mesh, filter => block%filter)
#ifdef Detail_Time
      call Get_Start_Time(cal_time_start)
#endif
      do im = 1 , member_num
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            if (filter%ngrid_lon(j) > 1) then
                n  = filter%ngrid_lon(j)
                hn = (n - 1) / 2
                do i = mesh%half_lon_ibeg, mesh%half_lon_iend
                tmp(im,i) = sum(filter%wgt_lon(:n,j) * x(im,i-hn:i+hn,j,k))
                end do
                if (present(y)) then
                y(im,:,j,k) = tmp(im,:)
                else
                x(im,:,j,k) = tmp(im,:)
                end if
            else if (present(y)) then
                y(im,:,j,k) = x(im,:,j,k)
            end if
          end do
        end do
      end do
#ifdef Detail_Time
      call Get_End_Time(cal_time_end)
      cal_time = cal_time + cal_time_end - cal_time_start
#endif
      end associate
  
    end subroutine filter_on_lon_edge
  
    subroutine filter_on_lat_edge(block, x, y)
  
      type(block_type), intent(in) :: block
      real(r8), intent(inout) :: x(member_num ,&
                                   block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out), optional :: y(member_num ,&
                                           block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                           block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                           block%mesh%full_lev_lb:block%mesh%full_lev_ub)
  
      real(r8) tmp(member_num,block%mesh%full_lon_lb:block%mesh%full_lon_ub)
      integer i, j, k, n, hn, im
  
      associate (mesh => block%mesh, filter => block%filter)
#ifdef Detail_Time
      call Get_Start_Time(cal_time_start)
#endif
      do im = 1 , member_num
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%half_lat_ibeg, mesh%half_lat_iend
            if (filter%ngrid_lat(j) > 1) then
                n  = filter%ngrid_lat(j)
                hn = (n - 1) / 2
                do i = mesh%full_lon_ibeg, mesh%full_lon_iend
                tmp(im,i) = sum(filter%wgt_lat(:n,j) * x(im,i-hn:i+hn,j,k))
                end do
                if (present(y)) then
                y(im,:,j,k) = tmp(im,:)
                else
                x(im,:,j,k) = tmp(im,:)
                end if
            else if (present(y)) then
                y(im,:,j,k) = x(im,:,j,k)
            end if
          end do
        end do
      end do
#ifdef Detail_Time
      call Get_End_Time(cal_time_end)
      cal_time = cal_time + cal_time_end - cal_time_start
#endif
      end associate
  
    end subroutine filter_on_lat_edge
  
    subroutine filter_on_lev_edge(block, x, y)
  
      type(block_type), intent(in) :: block
      real(r8), intent(inout) :: x(member_num ,&
                                   block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%half_lev_lb:block%mesh%half_lev_ub)
      real(r8), intent(out), optional :: y(member_num ,&
                                           block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                           block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                           block%mesh%half_lev_lb:block%mesh%half_lev_ub)
  
      real(r8) tmp(member_num ,block%mesh%full_lon_lb:block%mesh%full_lon_ub)
      integer i, j, k, n, hn, im
  
#ifdef Detail_Time
      call Get_Start_Time(cal_time_start)
#endif
      associate (mesh => block%mesh, filter => block%filter)
      do im = 1 , member_num
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            if (filter%ngrid_lon(j) > 1) then
                n  = filter%ngrid_lon(j)
                hn = (n - 1) / 2
                do i = mesh%full_lon_ibeg, mesh%full_lon_iend
                tmp(im,i) = sum(filter%wgt_lon(:n,j) * x(im,i-hn:i+hn,j,k))
                end do
                if (present(y)) then
                y(im,:,j,k) = tmp(im,:)
                else
                x(im,:,j,k) = tmp(im,:)
                end if
            else if (present(y)) then
                y(im,:,j,k) = x(im,:,j,k)
            end if
          end do
        end do
      end do
      end associate
#ifdef Detail_Time
      call Get_End_Time(cal_time_end)
      cal_time = cal_time + cal_time_end - cal_time_start
#endif
  
    end subroutine filter_on_lev_edge
  
    subroutine filter_vector_on_cell(block, u, v)
  
      type(block_type), intent(in) :: block
      real(r8), intent(inout) :: u(member_num ,&
                                   block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(inout) :: v(member_num ,&
                                   block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
  
      real(r8) us(member_num ,block%mesh%full_lon_lb:block%mesh%full_lon_ub)
      real(r8) vs(member_num ,block%mesh%full_lon_lb:block%mesh%full_lon_ub)
      real(r8) s
      real(r8) tmp(member_num ,block%mesh%full_lon_lb:block%mesh%full_lon_ub)
      integer i, j, k, n, hn, im
  
      associate (mesh => block%mesh, filter => block%filter)
      
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          if (filter%ngrid_lon(j) > 1) then
            s = sign(1.0_r8, mesh%full_lat(j))
            ! Transform onto polar plane.
#ifdef Detail_Time
            call Get_Start_Time(cal_time_start)
#endif
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              us(:,i) = -s * u(:,i,j,k) * mesh%full_sin_lon(i) - v(:,i,j,k) * mesh%full_cos_lon(i)
              vs(:,i) =  s * u(:,i,j,k) * mesh%full_cos_lon(i) - v(:,i,j,k) * mesh%full_sin_lon(i)
            end do
#ifdef Detail_Time
            call Get_End_Time(cal_time_end)
            cal_time = cal_time + cal_time_end - cal_time_start
#endif
#ifdef Detail_Time
            call Get_Start_Time(tran_time_start)
#endif
            call fill_zonal_halo_member(block, mesh%lon_halo_width, us)
            call fill_zonal_halo_member(block, mesh%lon_halo_width, vs)
#ifdef Detail_Time
            call Get_End_Time(tran_time_end)
            tran_time = tran_time + tran_time_end - tran_time_start
#endif
#ifdef Detail_Time
            call Get_Start_Time(cal_time_start)
#endif
            n  = filter%ngrid_lon(j)
            hn = (n - 1) / 2
            do im = 1 , member_num
              do i = mesh%full_lon_ibeg, mesh%full_lon_iend
                tmp(im,i) = sum(filter%wgt_lon(:n,j) * us(im,i-hn:i+hn))
              end do
              us(im,:) = tmp(im,:)
              do i = mesh%half_lon_ibeg, mesh%half_lon_iend
                tmp(im,i) = sum(filter%wgt_lon(:n,j) * vs(im,i-hn:i+hn))
              end do
              vs(im,:) = tmp(im,:)
            end do
            ! Transform back.
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              u(:,i,j,k) = s * (-us(:,i) * mesh%full_sin_lon(i) + vs(:,i) * mesh%full_cos_lon(i))
              v(:,i,j,k) =      -us(:,i) * mesh%full_cos_lon(i) - vs(:,i) * mesh%full_sin_lon(i)
            end do
#ifdef Detail_Time
            call Get_End_Time(cal_time_end)
            cal_time = cal_time + cal_time_end - cal_time_start
#endif
          end if
        end do
      end do
      end associate
  
    end subroutine filter_vector_on_cell
  
    subroutine filter_vector_on_edges(block, u, v, u_f, v_f)
  
      type(block_type), intent(in) :: block
      real(r8), intent(inout) :: u(member_num ,&
                                   block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(inout) :: v(member_num ,&
                                   block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(inout), optional :: u_f(member_num,&
                                               block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                               block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                               block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(inout), optional :: v_f(member_num,&
                                               block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                               block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                               block%mesh%full_lev_lb:block%mesh%full_lev_ub)
  
      real(r8) v_lon(member_num,block%mesh%half_lon_lb:block%mesh%half_lon_ub)
      real(r8) u_lat(member_num,block%mesh%full_lon_lb:block%mesh%full_lon_ub)
      real(r8) us(member_num,block%mesh%full_lon_lb:block%mesh%full_lon_ub)
      real(r8) vs(member_num,block%mesh%full_lon_lb:block%mesh%full_lon_ub)
      real(r8) s
      real(r8) tmp(member_num,block%mesh%full_lon_lb:block%mesh%full_lon_ub)
      integer i, j, k, n, hn, im
  
      associate (mesh => block%mesh, filter => block%filter)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          if (filter%ngrid_lon(j) > 1) then
#ifdef Detail_Time
            call Get_Start_Time(cal_time_start)
#endif
            s = sign(1.0_r8, mesh%full_lat(j))
            ! Transform onto polar plane.
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              v_lon(:,i) = mesh%full_tangent_wgt(1,j) * (v(:,i,j-1,k) + v(:,i+1,j-1,k)) + &
                           mesh%full_tangent_wgt(2,j) * (v(:,i,j  ,k) + v(:,i+1,j  ,k))
              us(:,i) = -s * u(:,i,j,k) * mesh%half_sin_lon(i) - v_lon(:,i) * mesh%half_cos_lon(i)
              vs(:,i) =  s * u(:,i,j,k) * mesh%half_cos_lon(i) - v_lon(:,i) * mesh%half_sin_lon(i)
            end do
#ifdef Detail_Time
            call Get_End_Time(cal_time_end)
            cal_time = cal_time + cal_time_end - cal_time_start
#endif
#ifdef Detail_Time
            call Get_Start_Time(tran_time_start)
#endif
            call fill_zonal_halo_member(block, mesh%lon_halo_width, us)
            call fill_zonal_halo_member(block, mesh%lon_halo_width, vs)
#ifdef Detail_Time
            call Get_End_Time(tran_time_end)
            tran_time = tran_time + tran_time_end - tran_time_start
#endif
#ifdef Detail_Time
            call Get_Start_Time(cal_time_start)
#endif
            n  = filter%ngrid_lon(j)
            hn = (n - 1) / 2
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              do im = 1 , member_num  
                tmp(im,i) = sum(filter%wgt_lon(:n,j) * us(im,i-hn:i+hn))
              end do
            end do
            us = tmp
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              do im = 1 , member_num
                tmp(im,i) = sum(filter%wgt_lon(:n,j) * vs(im,i-hn:i+hn))
              end do
            end do
            vs = tmp
            ! Transform back.
            if (present(u_f)) then
              do i = mesh%half_lon_ibeg, mesh%half_lon_iend
                u_f(:,i,j,k) = s * (-us(:,i) * mesh%half_sin_lon(i) + vs(:,i) * mesh%half_cos_lon(i))
              end do
            else
              do i = mesh%half_lon_ibeg, mesh%half_lon_iend
                u(:,i,j,k) = s * (-us(:,i) * mesh%half_sin_lon(i) + vs(:,i) * mesh%half_cos_lon(i))
              end do
            end if
#ifdef Detail_Time
            call Get_End_Time(cal_time_end)
            cal_time = cal_time + cal_time_end - cal_time_start
#endif
          else if (present(u_f)) then
#ifdef Detail_Time
            call Get_Start_Time(cal_time_start)
#endif
            u_f(:,:,j,k) = u(:,:,j,k)
#ifdef Detail_Time
            call Get_End_Time(cal_time_end)
            cal_time = cal_time + cal_time_end - cal_time_start
#endif
          end if
        end do
      end do
  
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          if (filter%ngrid_lat(j) > 1) then
#ifdef Detail_Time
            call Get_Start_Time(cal_time_start)
#endif
            s = sign(1.0_r8, mesh%half_lat(j))
            ! Transform onto polar plane.
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              u_lat(:,i) = mesh%half_tangent_wgt(1,j) * (u(:,i-1,j  ,k) + u(:,i,j  ,k)) + &
                           mesh%half_tangent_wgt(2,j) * (u(:,i-1,j+1,k) + u(:,i,j+1,k))
              us(:,i) = -s * u_lat(:,i) * mesh%full_sin_lon(i) - v(:,i,j,k) * mesh%full_cos_lon(i)
              vs(:,i) =  s * u_lat(:,i) * mesh%full_cos_lon(i) - v(:,i,j,k) * mesh%full_sin_lon(i)
            end do
#ifdef Detail_Time
            call Get_End_Time(cal_time_end)
            cal_time = cal_time + cal_time_end - cal_time_start
#endif
#ifdef Detail_Time
            call Get_Start_Time(tran_time_start)
#endif
            call fill_zonal_halo_member(block, mesh%lon_halo_width, us)
            call fill_zonal_halo_member(block, mesh%lon_halo_width, vs)
#ifdef Detail_Time
            call Get_End_Time(tran_time_end)
            tran_time = tran_time + tran_time_end - tran_time_start
#endif
#ifdef Detail_Time
            call Get_Start_Time(cal_time_start)
#endif
            n  = filter%ngrid_lat(j)
            hn = (n - 1) / 2
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              do im = 1 , member_num
                tmp(im,i) = sum(filter%wgt_lat(:n,j) * us(im,i-hn:i+hn))
              end do
            end do
            us = tmp
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              do im = 1 , member_num
                tmp(im,i) = sum(filter%wgt_lat(:n,j) * vs(im,i-hn:i+hn))
              end do
            end do
            vs = tmp
            ! Transform back.
            if (present(v_f)) then
              do i = mesh%full_lon_ibeg, mesh%full_lon_iend
                v_f(:,i,j,k) = -us(:,i) * mesh%full_cos_lon(i) - vs(:,i) * mesh%full_sin_lon(i)
              end do
            else
              do i = mesh%full_lon_ibeg, mesh%full_lon_iend
                v(:,i,j,k) = -us(:,i) * mesh%full_cos_lon(i) - vs(:,i) * mesh%full_sin_lon(i)
              end do
            end if
#ifdef Detail_Time
            call Get_End_Time(cal_time_end)
            cal_time = cal_time + cal_time_end - cal_time_start
#endif
          else if (present(v_f)) then
#ifdef Detail_Time
            call Get_Start_Time(cal_time_start)
#endif
            v_f(:,:,j,k) = v(:,:,j,k)
#ifdef Detail_Time
            call Get_End_Time(cal_time_end)
            cal_time = cal_time + cal_time_end - cal_time_start
#endif
          end if
        end do
      end do
      end associate
  
    end subroutine filter_vector_on_edges
  
  end module filter_mod
  