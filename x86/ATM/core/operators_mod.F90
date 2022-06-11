module operators_mod

  use const_mod
  use vert_coord_mod
  use block_mod
  use parallel_mod
  use parallel_types_mod
  use formula_mod
  use namelist_mod
  use time_mod
  use log_mod
  use pgf_mod
  use nh_mod
  use pgf_mod
  use interp_mod 
  use upwind_mod
  use pa_mod
  use member_mod

  implicit none

  private

  public operators_prepare
  public diag_ph
  public diag_m
  public diag_t
  public diag_gz_lev
  public calc_wedphdlev_lev
  public calc_div
  public calc_vor
  public calc_qhu_qhv
  public calc_dkedlon_dkedlat
  public calc_dmfdlon_dmfdlat
  public calc_dptfdlon_dptfdlat
  public calc_dptfdlev
  public calc_dphs
  public calc_wedudlev_wedvdlev
  public nh_prepare
  public nh_solve
  public interp_gz

  interface operators_prepare
    module procedure operators_prepare_1
    module procedure operators_prepare_2
  end interface operators_prepare

contains

  subroutine operators_prepare_1(blocks, itime, dt)

    type(block_type), intent(inout) :: blocks(:)
    integer, intent(in) :: itime
    real(8), intent(in) :: dt

    integer iblk

    do iblk = 1, size(blocks)
      if (baroclinic) then
        call diag_ph                  (blocks(iblk), blocks(iblk)%state(itime,ivector))
        call interp_pt                (blocks(iblk), blocks(iblk)%state(itime,ivector))
        call diag_t                   (blocks(iblk), blocks(iblk)%state(itime,ivector))
      end if
      call diag_m                     (blocks(iblk), blocks(iblk)%state(itime,ivector))
      if (nonhydrostatic) then
        !call diag_m_lev               (blocks(iblk), blocks(iblk)%state(itime,ivector))
      end if
      call interp_m_vtx               (blocks(iblk), blocks(iblk)%state(itime,ivector))
      call calc_mf                    (blocks(iblk), blocks(iblk)%state(itime,ivector))
      call calc_ke                    (blocks(iblk), blocks(iblk)%state(itime,ivector))
      call diag_pv                    (blocks(iblk), blocks(iblk)%state(itime,ivector))
      call interp_pv                  (blocks(iblk), blocks(iblk)%state(itime,ivector))
      call calc_div                   (blocks(iblk), blocks(iblk)%state(itime,ivector))
      if (hydrostatic) then
        call diag_gz_lev              (blocks(iblk), blocks(iblk)%state(itime,ivector))
      end if
      call pgf_prepare                (blocks(iblk), blocks(iblk)%state(itime,ivector))
    end do

  end subroutine operators_prepare_1

  subroutine operators_prepare_2(block, state, dt, pass)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state
    real(8), intent(in) :: dt
    integer, intent(in) :: pass

    select case (pass)
      ! --------------------------------------------------------------------------
      case (all_pass)
        if (baroclinic) then
          call diag_ph                  (block, state)
          call interp_pt                (block, state)
          call diag_t                   (block, state)
        end if
        call diag_m                     (block, state)
        if (nonhydrostatic) then
          !call diag_m_lev               (block, state)
        end if
        call calc_mf                    (block, state)
        call calc_ke                    (block, state)
        call calc_div                   (block, state)
        call interp_m_vtx               (block, state)
        call diag_pv                    (block, state)
        call interp_pv                  (block, state)
        if (hydrostatic) then
          call diag_gz_lev              (block, state)
        end if
        call pgf_prepare                (block, state)
      ! --------------------------------------------------------------------------
      case (forward_pass)
        if (baroclinic) then
          call interp_pt                (block, state)
        end if
        call calc_mf                    (block, state)
        call calc_ke                    (block, state)
        call calc_div                   (block, state)
        call interp_m_vtx               (block, state)
        call diag_pv                    (block, state)
        call interp_pv                  (block, state)
      ! --------------------------------------------------------------------------
      case (backward_pass)
        if (hydrostatic) then
          call diag_ph                  (block, state)
          call diag_t                   (block, state)
          call diag_gz_lev              (block, state)
        end if
        call diag_m                     (block, state)
        call pgf_prepare                (block, state)
      ! --------------------------------------------------------------------------
      case (no_wind_pass)
        if (baroclinic) then
          call diag_ph                  (block, state)
          call interp_pt                (block, state)
          call diag_t                   (block, state)
          call diag_gz_lev              (block, state)
        end if
        call diag_m                     (block, state)
        if (nonhydrostatic) then
          !call diag_m_lev               (block, state)
        end if
        call pgf_prepare                (block, state)
      end select
    
  end subroutine operators_prepare_2

  subroutine diag_ph(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k , im


    associate (mesh   => block%mesh  , &
               phs    => state%phs_f , & ! in
               ph_lev => state%ph_lev, & ! out
               ph     => state%ph)       ! out
               

    call state%async(async_phs)%wait()

    do k = mesh%half_lev_ibeg, mesh%half_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          do im = 1 , member_num
            ph_lev(im,i,j,k) = vert_coord_calc_ph_lev(k, phs(im,i,j))
          end do  
        end do
      end do
    end do



    call fill_halo_member(block, state%ph_lev, full_lon=.true., full_lat=.true., full_lev=.false.,async=state%async(async_ph_lev))

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          ph(:,i,j,k) = 0.5_r8 * (ph_lev(:,i,j,k) + ph_lev(:,i,j,k+1))
        end do
      end do
    end do

    call fill_halo_member(block, state%ph, full_lon=.true., full_lat=.true., full_lev=.true., async=state%async(async_ph))

    end associate
  end subroutine diag_ph

  subroutine diag_t(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k ,im

    associate (mesh => block%mesh, &
               pt   => state%pt_f, & ! in
               ph   => state%ph  , & ! in
               t    => state%t)      ! out

    call state%async(async_ph)%wait()
    call state%async(async_pt_f)%wait()

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          do im = 1 , member_num
            t(im,i,j,k) = temperature(pt(im,i,j,k), ph(im,i,j,k))
          end do
        end do
      end do
    end do

    call fill_halo_member(block, state%t, full_lon=.true., full_lat=.true., full_lev=.true., async=state%async(async_t))

    end associate

  end subroutine diag_t

  subroutine calc_wedphdlev_lev(block, state, tend, dt)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(in) :: tend
    real(8), intent(in) :: dt

    integer i, j, k, l, im
    real(r8) mf


    associate(mesh              => block%mesh             , &
              dmfdlon           => tend%dmfdlon           , & ! in
              dmfdlat           => tend%dmfdlat           , & ! in
              dphs              => tend%dphs              , & ! in
              wedphdlev_lev     => state%wedphdlev_lev    , & ! out
              wedphdlev_lev_lon => state%wedphdlev_lev_lon, & ! out
              wedphdlev_lev_lat => state%wedphdlev_lev_lat)   ! out
    
    !fix l loop !!!!!!!!!            
    do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          do im = 1 , member_num
            mf = 0.0_r8
            do l = 1, k - 1
              mf = mf + dmfdlon(im,i,j,l) + dmfdlat(im,i,j,l)
            end do
            wedphdlev_lev(im,i,j,k) = - vert_coord_calc_dphdt_lev(k, dphs(im,i,j)) - mf
          end do
        end do
      end do
    end do



    ! Set vertical boundary conditions.
    wedphdlev_lev(:,:,:,mesh%half_lev_ibeg) = 0.0_r8
    wedphdlev_lev(:,:,:,mesh%half_lev_iend) = 0.0_r8
      
    call fill_halo_member(block, state%wedphdlev_lev, full_lon=.true., full_lat=.true., full_lev=.false., south_halo=.false., async=state%async(async_wedphdlev_lev))

    call state%async(async_wedphdlev_lev)%wait()
      
    call interp_lev_edge_to_lev_lon_edge(mesh, state%wedphdlev_lev, state%wedphdlev_lev_lon)
    call interp_lev_edge_to_lev_lat_edge(mesh, state%wedphdlev_lev, state%wedphdlev_lev_lat)

    end associate
  
  end subroutine calc_wedphdlev_lev

  subroutine calc_ke(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    integer i, j, k
    real(r8) ke_vtx(member_num , 4)
    real(r8) work(member_num , state%mesh%full_lon_ibeg:state%mesh%full_lon_iend,state%mesh%num_full_lev)
    real(r8) pole(member_num , state%mesh%num_full_lev)
    real(r8) order_reduce(member_num , state%mesh%full_lon_ibeg : state%mesh%full_lon_iend)


    call state%async(async_u_f)%wait()
    call state%async(async_v_f)%wait()

    associate (mesh => block%mesh, &
                u    => state%u_f , & ! in
                v    => state%v_f , & ! in
                ke   => state%ke)     ! out

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend + 1
          ke(:,i,j,k) = (mesh%area_lon_west (j  ) * u(:,i-1,j  ,k)**2  + &
                          mesh%area_lon_east (j  ) * u(:,i  ,j  ,k)**2 + &
                          mesh%area_lat_north(j-1) * v(:,i  ,j-1,k)**2 + &
                          mesh%area_lat_south(j  ) * v(:,i  ,j  ,k)**2   &
                        ) / mesh%area_cell(j)
        end do
      end do
    end do

    if (ke_scheme == 2) then

    !
    !      ________u_________________u________
    !     |     i-1,j+1     |       i,j+1     |
    !     |                 |                 |
    !     |                 |                 |
    !     |        1        |        4        |
    !     v        o--------v--------o        v
    !  i-1,j    i-1,j      i,j      i,j    i+1,j
    !     |        |        |        |        |
    !     |        |        |        |        |
    !     |________u________|________u________|
    !     |     i-1,j      i,j      i,j       |
    !     |        |        |        |        |
    !     |        |        |        |        |
    !     |        |        |        |        |
    !     v        o--------v--------o        v
    !  i-1,j-1  i-1,j-1    i,j-1    i,j-1  i+1,j-1
    !     |        2        |        3        |
    !     |                 |                 |
    !     |________u________|________u________|
    !           i-1,j-1             i,j-1
    !

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole + merge(0, 1, mesh%has_north_pole())
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend +1
            ke_vtx(:,1) = (                                        &
              mesh%area_lat_east (j  ) * v(:,i-1,j  ,k)**2 + &
              mesh%area_lat_west (j  ) * v(:,i  ,j  ,k)**2 + &
              mesh%area_lon_north(j  ) * u(:,i-1,j  ,k)**2 + &
              mesh%area_lon_south(j+1) * u(:,i-1,j+1,k)**2   &
            ) / mesh%area_vtx(j)
            ke_vtx(:,2) = (                                        &
              mesh%area_lat_east (j-1) * v(:,i-1,j-1,k)**2 + &
              mesh%area_lat_west (j-1) * v(:,i  ,j-1,k)**2 + &
              mesh%area_lon_north(j-1) * u(:,i-1,j-1,k)**2 + &
              mesh%area_lon_south(j  ) * u(:,i-1,j  ,k)**2   &
            ) / mesh%area_vtx(j-1)
            ke_vtx(:,3) = (                                        &
              mesh%area_lat_east (j-1) * v(:,i  ,j-1,k)**2 + &
              mesh%area_lat_west (j-1) * v(:,i+1,j-1,k)**2 + &
              mesh%area_lon_north(j-1) * u(:,i  ,j-1,k)**2 + &
              mesh%area_lon_south(j  ) * u(:,i  ,j  ,k)**2   &
            ) / mesh%area_vtx(j-1)
            ke_vtx(:,4) = (                                        &
              mesh%area_lat_east (j  ) * v(:,i  ,j  ,k)**2 + &
              mesh%area_lat_west (j  ) * v(:,i+1,j  ,k)**2 + &
              mesh%area_lon_north(j  ) * u(:,i  ,j  ,k)**2 + &
              mesh%area_lon_south(j+1) * u(:,i  ,j+1,k)**2   &
            ) / mesh%area_vtx(j)
            ke(:,i,j,k) = (1.0_r8 - ke_cell_wgt) * (               &
              (ke_vtx(:,1) + ke_vtx(:,4)) * mesh%area_subcell(2,j) + &
              (ke_vtx(:,2) + ke_vtx(:,3)) * mesh%area_subcell(1,j)   &
            ) / mesh%area_cell(j) + ke_cell_wgt * ke(:,i,j,k)
          end do
        end do
      end do
    end if

    ! Note: area_lat_south and area_lat_north at the Poles is the same as area_cell.
    if (mesh%has_south_pole()) then
      j = mesh%full_lat_ibeg
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          work(:,i,k) = v(:,i,j,k)**2
        end do
      end do
      call zonal_sum_member(proc%zonal_circle, work, pole)
      pole = pole / global_mesh%num_full_lon
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          ke(:,i,j,k) = pole(:,k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_lat_iend
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          work(:,i,k) = v(:,i,j-1,k)**2
        end do
      end do
      call zonal_sum_member(proc%zonal_circle, work, pole)
      pole = pole / global_mesh%num_full_lon
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          ke(:,i,j,k) = pole(:,k)
        end do
      end do
    end if

    !call fill_halo_member(block, ke, full_lon=.true., full_lat=.true., full_lev=.true., async=state%async(async_ke))

    end associate

  end subroutine calc_ke

  subroutine calc_div(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    real(r8) pole(member_num , global_mesh%num_full_lev)
    integer i, j, k
    real(r8) order_reduce(member_num ,  state%mesh%full_lon_ibeg : state%mesh%full_lon_iend)
    real(r8) work(member_num , state%mesh%full_lon_ibeg:state%mesh%full_lon_iend,state%mesh%num_full_lev)


    associate (mesh => block%mesh, &
               u    => state%u   , & ! in
               v    => state%v   , & ! in
               div  => state%div , & ! out
               div2 => state%div2)   ! out

    call state%async(async_u)%wait()
    call state%async(async_v)%wait()

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          div(:,i,j,k) = (                                                          &
            (u(:,i,j,k) * mesh%le_lon(j) - u(:,i-1,  j,k) * mesh%le_lon(j  )) + &
            (v(:,i,j,k) * mesh%le_lat(j) - v(:,i  ,j-1,k) * mesh%le_lat(j-1))   &
          ) / mesh%area_cell(j)
        end do
      end do
    end do



    if (mesh%has_south_pole()) then
      j = mesh%full_lat_ibeg
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          work(:,i,k) = state%v(:,i,j,k)
        end do
      end do
      call zonal_sum_member(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%num_full_lon / mesh%area_cell(j)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%div(:,i,j,k) = pole(:,k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_lat_iend
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          work(:,i,k) = - state%v(:,i,j-1,k)
        end do
      end do
      call zonal_sum_member(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%num_full_lon / mesh%area_cell(j)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%div(:,i,j,k) = pole(:,k)
        end do
      end do
    end if

    if (div_damp_order == 4) then
      call fill_halo_member(block, div, full_lon=.true., full_lat=.true., full_lev=.true., async=state%async(async_div))
    else
      call fill_halo_member(block, div, full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false., async=state%async(async_div))
    end if

    call state%async(async_div)%wait()

    if (div_damp_order == 4) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            div2(:,i,j,k) = (                                                                 &
              div(:,i+1,j,k) - 2 * div(:,i,j,k) + div(:,i-1,j,k)                              &
            ) / mesh%de_lon(j)**2 + (                                                         &
              (div(:,i,j+1,k) - div(:,i,j  ,k)) * mesh%half_cos_lat(j  ) / mesh%de_lat(j  ) - &
              (div(:,i,j  ,k) - div(:,i,j-1,k)) * mesh%half_cos_lat(j-1) / mesh%de_lat(j-1)   &
            ) / mesh%le_lon(j) / mesh%full_cos_lat(j)
          end do
        end do
      end do
      call fill_halo_member(block, state%div2, full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false., async=state%async(async_div2))
    end if
    end associate

  end subroutine calc_div

  subroutine diag_gz_lev(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    integer i, j, k, l
    !real(r8) dgz(member_num, block%mesh%full_lon_ibeg:block%mesh%full_lon_iend, block%mesh%full_lat_ibeg:block%mesh%full_lat_iend)
    real(r8) dgz(member_num)

    call state%async(async_t)%wait()
    call state%async(async_ph_lev)%wait()
    associate (mesh   => block%mesh      , &
               t      => state%t         , &          ! in
               ph_lev => state%ph_lev    , &          ! in
               gzs    => block%static(ivector)%gzs, & ! in
               gz_lev => state%gz_lev)                ! out

    ! !fix dong
    do k = mesh%half_lev_ibeg, mesh%half_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          dgz = 0.0_r8
          do l = k, mesh%num_full_lev
            dgz(:) = dgz(:) + Rd * t(:,i,j,l) * log(ph_lev(:,i,j,l+1) / ph_lev(:,i,j,l))
          end do
          gz_lev(:,i,j,k) = gzs(:,i,j) + dgz(:)
        end do
      end do
    end do


    
    !!!!! Maybe problem
    ! dgz(:,:,:) = 0.0_r8
    ! gz_lev(:,:,:,mesh%half_lev_iend) = gzs(:,:,:)
    ! do k = mesh%half_lev_iend - 1, mesh%half_lev_ibeg, -1
    !   do j = mesh%full_lat_ibeg, mesh%full_lat_iend
    !     do i = mesh%full_lon_ibeg, mesh%full_lon_iend
    !       dgz(:,i,j) = dgz(:,i,j) + Rd * t(:,i,j,k) * log(ph_lev(:,i,j,k+1) / ph_lev(:,i,j,k))
    !       gz_lev(:,i,j,k) = gzs(:,i,j) + dgz(:,i,j)
    !     end do
    !   end do
    ! end do
          
    call fill_halo_member(block, state%gz_lev, full_lon=.true., full_lat=.true., full_lev=.false., async=state%async(async_gz_lev))

  end associate
   
  end subroutine diag_gz_lev

  subroutine diag_m(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k

    associate (mesh   => block%mesh      , &
               ph_lev => state%ph_lev    , & ! in
               gz     => state%gz_f      , & ! in
               gzs    => block%static(ivector)%gzs, & ! in
               m      => state%m         , & ! out
               m_lon  => state%m_lon     , & ! out
               m_lat  => state%m_lat)        ! out


    if (baroclinic) then

      call state%async(async_ph_lev)%wait()

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend + merge(0, 1, mesh%has_north_pole())
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend + 1
            m(:,i,j,k) = ph_lev(:,i,j,k+1) - ph_lev(:,i,j,k)
          end do
        end do
      end do
    else
      call state%async(async_gz)%wait()
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend + 1
          m(:,i,j,1) = (gz(:,i,j,1) - gzs(:,i,j)) / g
        end do
      end do
    end if
    call average_cell_to_lon_edge(mesh, m, m_lon)
    call fill_halo_member(block, state%m_lon, full_lon=.false., full_lat=.true., full_lev=.true., async=state%async(async_m_lon))
    call average_cell_to_lat_edge(mesh, m, m_lat)
    call fill_halo_member(block, state%m_lat, full_lon=.true., full_lat=.false., full_lev=.true., async=state%async(async_m_lat))
    
    end associate

  end subroutine diag_m

  subroutine interp_pt(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh

#ifdef Detail_Time
  call Get_Time_Pa(wait_time_start)
#endif

    call state%async(async_pt_f)%wait()
    call state%async(async_u_f)%wait()
    call state%async(async_v_f)%wait()

    associate (u      => state%u_f          , & ! in
               v      => state%v_f          , & ! in
               w      => state%wedphdlev_lev, & ! in
               pt     => state%pt_f         , & ! in
               pt_lon => state%pt_lon       , & ! out
               pt_lat => state%pt_lat       , & ! out
               pt_lev => state%pt_lev)          ! out
    call interp_cell_to_lon_edge(state%mesh, pt, pt_lon, reversed_area=.true., u=u, upwind_wgt_=upwind_wgt_pt)
    call interp_cell_to_lat_edge(state%mesh, pt, pt_lat, reversed_area=.true., v=v, upwind_wgt_=upwind_wgt_pt)   
    call fill_halo_member(block, state%pt_lon, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false., async=state%async(async_pt_lon))
    call fill_halo_member(block, state%pt_lat, full_lon=.true., full_lat=.false., full_lev=.true., west_halo=.false.,  east_halo=.false., north_halo=.false., async=state%async(async_pt_lat))  
    call interp_cell_to_lev_edge(state%mesh, pt, pt_lev, w=w, upwind_wgt_=upwind_wgt_pt)
    end associate
  end subroutine interp_pt

  subroutine interp_m_vtx(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    call interp_cell_to_vtx(state%mesh, state%m, state%m_vtx)

  end subroutine interp_m_vtx

  subroutine calc_mf(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k

    associate (mesh     => block%mesh    , &
               m_lon    => state%m_lon   , & ! in
               m_lat    => state%m_lat   , & ! in
               u        => state%u_f     , & ! in
               v        => state%v_f     , & ! in
               mf_lon_n => state%mf_lon_n, & ! out
               mf_lat_n => state%mf_lat_n, & ! out
               mf_lon_t => state%mf_lon_t, & ! out
               mf_lat_t => state%mf_lat_t)   ! out

    call state%async(async_m_lon)%wait()
    call state%async(async_u_f)%wait()

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole + merge(0, 1, mesh%has_north_pole())
        do i = mesh%half_lon_ibeg - 1, mesh%half_lon_iend
          mf_lon_n(:,i,j,k) = m_lon(:,i,j,k) * u(:,i,j,k)
        end do
      end do
    end do



    call state%async(async_m_lat)%wait()
    call state%async(async_v_f)%wait()

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg - merge(0, 1, mesh%has_south_pole()), mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend + 1
          mf_lat_n(:,i,j,k) = m_lat(:,i,j,k) * v(:,i,j,k)
        end do
      end do
    end do

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          mf_lat_t(:,i,j,k) = mesh%half_tangent_wgt(1,j) * (mf_lon_n(:,i-1,j  ,k) + mf_lon_n(:,i,j  ,k)) + &
                              mesh%half_tangent_wgt(2,j) * (mf_lon_n(:,i-1,j+1,k) + mf_lon_n(:,i,j+1,k))
        end do
      end do
    end do

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          mf_lon_t(:,i,j,k) = mesh%full_tangent_wgt(1,j) * (mf_lat_n(:,i,j-1,k) + mf_lat_n(:,i+1,j-1,k)) + &
                              mesh%full_tangent_wgt(2,j) * (mf_lat_n(:,i,j  ,k) + mf_lat_n(:,i+1,j  ,k))
        end do
      end do
    end do

    end associate
  end subroutine calc_mf

  subroutine calc_vor(block, state, u, v)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state
    real(r8), intent(in) :: u(member_num ,&
                              block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                              block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                              block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in) :: v(member_num ,&
                              block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                              block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                              block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    real(r8) pole(member_num , state%mesh%num_full_lev)
    integer i, j, k
    real(r8) order_reduce(member_num , state%mesh%half_lon_ibeg : state%mesh%half_lon_iend)
    real(r8) work(member_num , state%mesh%half_lon_ibeg:state%mesh%half_lon_iend,state%mesh%num_full_lev)
     

    associate (mesh     => block%mesh    , &
               mf_lat_t => state%mf_lat_t, & ! in
               m_lat    => state%m_lat   , & ! in
               vor      => state%vor)        ! out

    call state%async(async_u)%wait()
    call state%async(async_v)%wait()

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg - merge(0, 1, mesh%has_south_pole()), mesh%half_lat_iend
        do i = mesh%half_lon_ibeg - 1, mesh%half_lon_iend
          vor(:,i,j,k) = (                                                                 &
            u(:,i  ,j  ,k) * mesh%de_lon(j  ) - u(:,i  ,j+1,k) * mesh%de_lon(j+1) + &
            v(:,i+1,j  ,k) * mesh%de_lat(j  ) - v(:,i  ,j  ,k) * mesh%de_lat(j  )   &
          ) / mesh%area_vtx(j)
        end do
      end do
    end do

    if (pv_pole_stokes) then
      ! Special treatment of vorticity around Poles
      if (mesh%has_south_pole()) then
        j = mesh%half_lat_ibeg
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            work(:,i,k) = -mf_lat_t(:,i,j,k) / m_lat(:,i,j,k) * mesh%le_lat(j)
          end do
        end do
        call zonal_sum_member(proc%zonal_circle, work, pole)
        pole = pole / global_mesh%num_full_lon / mesh%area_cell(j)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            vor(:,i,j,k) = vor(:,i,j+1,k) / 3.0_r8 + pole(:,k) * 2.0_r8 / 3.0_r8
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%half_lat_iend
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            work(:,i,k) = mf_lat_t(:,i,j,k) / m_lat(:,i,j,k) * mesh%le_lat(j)
          end do
        end do
        call zonal_sum_member(proc%zonal_circle, work, pole)
        pole = pole / global_mesh%num_full_lon / mesh%area_cell(j+1)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            vor(:,i,j,k) = vor(:,i,j-1,k) / 3.0_r8 + pole(:,k) * 2.0_r8 / 3.0_r8
          end do
        end do
      end if
    end if
    end associate
  end subroutine calc_vor

  subroutine diag_pv(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    integer i, j, k

    associate (mesh  => block%mesh , &
               m_vtx => state%m_vtx, & ! in
               u     => state%u_f  , & ! in
               v     => state%v_f  , & ! in
               vor   => state%vor  , & ! in
               pv    => state%pv)      ! out

    call calc_vor(block, state, u, v)

    call state%async(async_vor)%wait()

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          pv(:,i,j,k) = (vor(:,i,j,k) + mesh%half_f(j)) / m_vtx(:,i,j,k)
        end do
      end do
    end do

    call fill_halo_member(block, state%pv, full_lon=.false., full_lat=.false., full_lev=.true., async=state%async(async_pv))
    
    end associate
  end subroutine diag_pv

  subroutine interp_pv_midpoint(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    integer i, j, k

    associate (mesh   => block%mesh  , &
               pv     => state%pv    , & ! in
               pv_lon => state%pv_lon, & ! out
               pv_lat => state%pv_lat)   ! out

    call state%async(async_pv)%wait()

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
         pv_lat(:,i,j,k) = 0.5_r8 * (pv(:,i-1,j,k) + pv(:,i,j,k))
        end do
      end do
    end do
    
    call fill_halo_member(block, state%pv_lat, full_lon=.true., full_lat=.false., full_lev=.true., north_halo=.false., async=state%async(async_pv_lat))

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          pv_lon(:,i,j,k) = 0.5_r8 * (pv(:,i,j,k) + pv(:,i,j-1,k))
        end do
      end do
    end do

    call fill_halo_member(block, state%pv_lon, full_lon=.false., full_lat=.true., full_lev=.true., async=state%async(async_pv_lon))
    end associate
  end subroutine interp_pv_midpoint

  subroutine interp_pv_upwind(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    real(r8) ut, vt, b
    integer i, j, k, im

    associate (mesh     => block%mesh    , &
               un       => state%u_f       , & ! in
               vn       => state%v_f       , & ! in
               m_lon    => state%m_lon   , & ! in
               m_lat    => state%m_lat   , & ! in
               mf_lon_t => state%mf_lon_t, & ! in
               mf_lat_t => state%mf_lat_t, & ! in
               pv       => state%pv      , & ! in
               pv_lon   => state%pv_lon  , & ! out
               pv_lat   => state%pv_lat)     ! out


      call state%async(async_u_f)%wait()
      call state%async(async_v_f)%wait()
      call state%async(async_m_lon)%wait()
      call state%async(async_m_lat)%wait()
      call state%async(async_pv)%wait()
         
      select case (upwind_order_pv)
      case (1)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              do im = 1 , member_num
                vt = mf_lon_t(im,i,j,k) / m_lon(im,i,j,k)
                b  = abs(vt) / (sqrt(un(im,i,j,k)**2 + vt**2) + eps)
                pv_lon(im,i,j,k) = b * upwind1(sign(1.0_r8, vt), upwind_wgt_pv, pv(im,i,j-1:j,k)) + &
                                (1 - b) * 0.5_r8 * (pv(im,i,j-1,k) + pv(im,i,j,k))
              end do 
            end do
          end do
        end do
        call fill_halo_member(block, pv_lon, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false.,async=state%async(async_pv_lon))
  
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              do im = 1 , member_num
                ut = mf_lat_t(im,i,j,k) / m_lat(im,i,j,k)
                b  = abs(ut) / (sqrt(ut**2 + vn(im,i,j,k)**2) + eps)
                pv_lat(im,i,j,k) = b * upwind1(sign(1.0_r8, ut), upwind_wgt_pv, pv(im,i-1:i,j,k)) + &
                                (1 - b) * 0.5_r8 * (pv(im,i-1,j,k) + pv(im,i,j,k))
              end do
            end do
          end do
        end do
        call fill_halo_member(block, pv_lat, full_lon=.true., full_lat=.false., full_lev=.true., north_halo=.false.,async=state%async(async_pv_lat))
      case (3)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            if (mesh%is_full_lat_next_to_pole(j)) then
              do i = mesh%half_lon_ibeg, mesh%half_lon_iend
                do im = 1 , member_num
                  vt = mf_lon_t(im,i,j,k) / m_lon(im,i,j,k)
                  b  = abs(vt) / (sqrt(un(im,i,j,k)**2 + vt**2) + eps)
                  pv_lon(im,i,j,k) = b * upwind1(sign(1.0_r8, vt), upwind_wgt_pv, pv(im,i,j-1:j,k)) + &
                                    (1 - b) * 0.5_r8 * (pv(im,i,j-1,k) + pv(im,i,j,k))
                end do
              end do
            else
              do i = mesh%half_lon_ibeg, mesh%half_lon_iend
                do im = 1 , member_num
                  vt = mf_lon_t(im,i,j,k) / m_lon(im,i,j,k)
                  b  = abs(vt) / (sqrt(un(im,i,j,k)**2 + vt**2) + eps)
                  pv_lon(im,i,j,k) = b * upwind3(sign(1.0_r8, vt), upwind_wgt_pv, pv(im,i,j-2:j+1,k)) + &
                                    (1 - b) * 0.5_r8 * (pv(im,i,j-1,k) + pv(im,i,j,k))
                end do
              end do
            end if
          end do
        end do
        
        call fill_halo_member(block, pv_lon, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false., async=state%async(async_pv_lon))

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              do im = 1 , member_num
                ut = mf_lat_t(im,i,j,k) / m_lat(im,i,j,k)
                b  = abs(ut) / (sqrt(ut**2 + vn(im,i,j,k)**2) + eps)
                pv_lat(im,i,j,k) = b * upwind3(sign(1.0_r8, ut), upwind_wgt_pv, pv(im,i-2:i+1,j,k)) + &
                                  (1 - b) * 0.5_r8 * (pv(im,i-1,j,k) + pv(im,i,j,k))
              end do
            end do
          end do
        end do

        call fill_halo_member(block, pv_lat, full_lon=.true., full_lat=.false., full_lev=.true., north_halo=.false., async=state%async(async_pv_lat))

      end select
    end associate

  end subroutine interp_pv_upwind

  subroutine interp_pv(block, state)
    
    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state 

#ifdef Detail_Time
  call Get_Time_Pa(wait_time_start)
#endif

    call state%async(async_pv)%wait()

#ifdef Detail_Time
  call Get_Time_Pa(wait_time_end)
  wait_time = wait_time + wait_time_end - wait_time_start
#endif

    select case (pv_scheme)
    case (1)
      call interp_pv_midpoint(block, state)
    case (2)
      call interp_pv_upwind(block, state)
    case default
      if (is_root_proc()) call log_error('Unknown PV scheme!')
    end select


  end subroutine interp_pv

  subroutine calc_qhu_qhv(block, state, tend, dt)

    type(block_type), intent(inout), target :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real(8), intent(in) :: dt

    integer i, j, k, move, rf


#ifdef Detail_Time
  call Get_Time_Pa(wait_time_start)
#endif

    call state%async(async_mf_lon_n)%wait()
    call state%async(async_mf_lat_n)%wait()
    call state%async(async_pv_lon)%wait()
    call state%async(async_pv_lat)%wait()

#ifdef Detail_Time
  call Get_Time_Pa(wait_time_end)
  wait_time = wait_time + wait_time_end - wait_time_start
#endif

    associate (mesh          => block%mesh         , &
               mf_lon_n      => state%mf_lon_n     , & ! in
               mf_lat_n      => state%mf_lat_n     , & ! in
               mf_lon_t      => state%mf_lon_t     , & ! in
               mf_lat_t      => state%mf_lat_t     , & ! in
               pv_lon        => state%pv_lon       , & ! in
               pv_lat        => state%pv_lat       , & ! in
               qhu           => tend%qhu           , & ! out
               qhv           => tend%qhv)              ! out
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          if (coriolis_scheme == 1) then
            qhu(:,i,j,k) = (                                                  &
              mesh%half_tangent_wgt(1,j) * (                                &
                mf_lon_n(:,i-1,j  ,k) * (pv_lat(:,i,j,k) + pv_lon(:,i-1,j  ,k)) + &
                mf_lon_n(:,i  ,j  ,k) * (pv_lat(:,i,j,k) + pv_lon(:,i  ,j  ,k))   &
              ) +                                                           &
              mesh%half_tangent_wgt(2,j) * (                                &
                mf_lon_n(:,i-1,j+1,k) * (pv_lat(:,i,j,k) + pv_lon(:,i-1,j+1,k)) + &
                mf_lon_n(:,i  ,j+1,k) * (pv_lat(:,i,j,k) + pv_lon(:,i  ,j+1,k))   &
              )                                                             &
            ) * 0.5_r8
          else if (coriolis_scheme == 2) then
            qhu(:,i,j,k) = mf_lat_t(:,i,j,k) * pv_lat(:,i,j,k)
          end if
        end do
      end do
    end do
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          if (coriolis_scheme == 1) then
            qhv(:,i,j,k) = (                                                  &
              mesh%full_tangent_wgt(1,j) * (                                &
                mf_lat_n(:,i  ,j-1,k) * (pv_lon(:,i,j,k) + pv_lat(:,i  ,j-1,k)) + &
                mf_lat_n(:,i+1,j-1,k) * (pv_lon(:,i,j,k) + pv_lat(:,i+1,j-1,k))   &
              ) +                                                           &
              mesh%full_tangent_wgt(2,j) * (                                &
                mf_lat_n(:,i  ,j  ,k) * (pv_lon(:,i,j,k) + pv_lat(:,i  ,j  ,k)) + &
                mf_lat_n(:,i+1,j  ,k) * (pv_lon(:,i,j,k) + pv_lat(:,i+1,j  ,k))   &
              )                                                             &
            ) * 0.5_r8
          else if (coriolis_scheme == 2) then
            qhv(:,i,j,k) = mf_lon_t(:,i,j,k) * pv_lon(:,i,j,k)
          end if
        end do
      end do
    end do
    end associate

  end subroutine calc_qhu_qhv

  subroutine calc_dkedlon_dkedlat(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real(8), intent(in) :: dt

    integer i, j, k

    associate (mesh    => block%mesh  , &
               ke      => state%ke    , & ! in
               dkedlon => tend%dkedlon, & ! out
               dkedlat => tend%dkedlat)   ! out
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          dkedlon(:,i,j,k) = (ke(:,i+1,j,k) - ke(:,i,j,k)) / mesh%de_lon(j)
        end do
      end do
    end do
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          dkedlat(:,i,j,k) = (ke(:,i,j+1,k) - ke(:,i,j,k)) / mesh%de_lat(j)
        end do
      end do
    end do
    end associate

  end subroutine calc_dkedlon_dkedlat

  subroutine calc_dmfdlon_dmfdlat(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real(8), intent(in) :: dt
    integer i, j, k, move, rf
    real(r8) work(member_num , state%mesh%full_lon_ibeg:state%mesh%full_lon_iend,state%mesh%num_full_lev)
    real(r8) pole(member_num , state%mesh%num_full_lev)

    associate (mesh     => block%mesh    , &
               mf_lon_n => state%mf_lon_n, & ! in
               mf_lat_n => state%mf_lat_n, & ! in
               dmfdlon  => tend%dmfdlon  , & ! out
               dmfdlat  => tend%dmfdlat)     ! out
    
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          dmfdlon(:,i,j,k) = (                    &
            mf_lon_n(:,i,j,k) - mf_lon_n(:,i-1,j,k) &
          ) * mesh%le_lon(j) / mesh%area_cell(j)
        end do
      end do
    end do
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          dmfdlat(:,i,j,k) = (                       &
            mf_lat_n(:,i,j  ,k) * mesh%le_lat(j  ) - &
            mf_lat_n(:,i,j-1,k) * mesh%le_lat(j-1)   &
          ) / mesh%area_cell(j)
        end do
      end do
    end do

    if (mesh%has_south_pole()) then
      j = mesh%full_lat_ibeg
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          work(:,i,k) = state%mf_lat_n(:,i,j,k)
        end do
      end do
      call zonal_sum_member(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%num_full_lon / mesh%area_cell(j)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%dmfdlat(:,i,j,k) = pole(:,k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_lat_iend
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          work(:,i,k) = - state%mf_lat_n(:,i,j-1,k)
        end do
      end do
      call zonal_sum_member(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%num_full_lon / mesh%area_cell(j)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%dmfdlat(:,i,j,k) = pole(:,k)
        end do
      end do
    end if
    end associate
  end subroutine calc_dmfdlon_dmfdlat

  subroutine calc_dptfdlon_dptfdlat(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real(8), intent(in) :: dt

    integer i, j, k, move, rf
    real(r8) work(member_num , state%mesh%full_lon_ibeg:state%mesh%full_lon_iend,state%mesh%num_full_lev)
    real(r8) pole(member_num , state%mesh%num_full_lev)

    call state%async(async_pt_lon)%wait()
    call state%async(async_pt_lat)%wait()

    associate (mesh     => block%mesh    , &
               mf_lon_n => state%mf_lon_n, & ! in
               mf_lat_n => state%mf_lat_n, & ! in
               pt_lon   => state%pt_lon  , & ! in
               pt_lat   => state%pt_lat  , & ! in
               dptfdlon => tend%dptfdlon , & ! out
               dptfdlat => tend%dptfdlat)    ! out

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          dptfdlon(:,i,j,k) = (                     &
            mf_lon_n(:,i  ,j,k) * pt_lon(:,i  ,j,k) - &
            mf_lon_n(:,i-1,j,k) * pt_lon(:,i-1,j,k)   &
          ) * mesh%le_lon(j) / mesh%area_cell(j)
        end do
      end do
    end do
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          dptfdlat(:,i,j,k) = (                                        &
            mf_lat_n(:,i,j  ,k) * pt_lat(:,i,j  ,k) * mesh%le_lat(j  ) - &
            mf_lat_n(:,i,j-1,k) * pt_lat(:,i,j-1,k) * mesh%le_lat(j-1)   &
          ) / mesh%area_cell(j)
        end do
      end do
    end do

    if (mesh%has_south_pole()) then
      j = mesh%full_lat_ibeg
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          work(:,i,k) = state%mf_lat_n(:,i,j,k) * state%pt_lat(:,i,j,k)
        end do
      end do
      call zonal_sum_member(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%num_full_lon / mesh%area_cell(j)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%dptfdlat(:,i,j,k) = pole(:,k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_lat_iend
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          work(:,i,k) = - state%mf_lat_n(:,i,j-1,k) * state%pt_lat(:,i,j-1,k)
        end do
      end do
      call zonal_sum_member(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%num_full_lon / mesh%area_cell(j)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%dptfdlat(:,i,j,k) = pole(:,k)
        end do
      end do
    end if

    end associate
  end subroutine calc_dptfdlon_dptfdlat

  subroutine calc_dptfdlev(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real(8), intent(in) :: dt

    integer i, j, k

    associate (mesh          => block%mesh         , &
                wedphdlev_lev => state%wedphdlev_lev, &
                pt_lev        => state%pt_lev       , &
                dptfdlev      => tend%dptfdlev)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            dptfdlev(:,i,j,k) = wedphdlev_lev(:,i,j,k+1) * pt_lev(:,i,j,k+1) - &
                                wedphdlev_lev(:,i,j,k  ) * pt_lev(:,i,j,k  )
          end do
        end do
      end do
    end associate

  end subroutine calc_dptfdlev

  subroutine calc_dphs(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt


    integer i, j, k

    associate (mesh    => block%mesh  , &
      dmfdlon => tend%dmfdlon, & ! in
      dmfdlat => tend%dmfdlat, & ! in
      dphs    => tend%dphs)      ! out

    tend%dphs = 0.0_r8
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          dphs(:,i,j) = dphs(:,i,j) - dmfdlon(:,i,j,k) - dmfdlat(:,i,j,k)
        end do
      end do
    end do
    end associate
  end subroutine calc_dphs

  subroutine calc_wedudlev_wedvdlev(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real(8), intent(in) :: dt

    integer i, j, k

#ifdef Detail_Time
  call Get_Time_Pa(wait_time_start)
#endif

    call state%async(async_u)%wait()
    call state%async(async_v)%wait()
    call state%async(async_m_lon)%wait()
    call state%async(async_m_lat)%wait()

#ifdef Detail_Time
  call Get_Time_Pa(wait_time_end)
  wait_time = wait_time + wait_time_end - wait_time_start
#endif

    associate (mesh              => block%mesh             , &
              u                 => state%u_f              , & ! in
              v                 => state%v_f              , & ! in
              m_lon             => state%m_lon            , & ! in
              m_lat             => state%m_lat            , & ! in
              wedphdlev_lev_lon => state%wedphdlev_lev_lon, & ! in
              wedphdlev_lev_lat => state%wedphdlev_lev_lat, & ! in
              wedudlev          => tend%wedudlev          , & ! out
              wedvdlev          => tend%wedvdlev)             ! out
  
    do k = mesh%full_lev_ibeg + 1, mesh%full_lev_iend - 1
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          wedudlev(:,i,j,k) = (                                          &
              wedphdlev_lev_lon(:,i,j,k+1) * (u(:,i,j,k+1) - u(:,i,j,k  )) + &
              wedphdlev_lev_lon(:,i,j,k  ) * (u(:,i,j,k  ) - u(:,i,j,k-1))   &
            ) / m_lon(:,i,j,k) / 2.0_r8
        end do
      end do
    end do
    k = mesh%full_lev_ibeg
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        wedudlev(:,i,j,k) = (wedphdlev_lev_lon(:,i,j,k+1) * &
          (u(:,i,j,k+1) - u(:,i,j,k  ))) / m_lon(:,i,j,k) / 2.0_r8
      end do
    end do
    k = mesh%full_lev_iend
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        wedudlev(:,i,j,k) = (wedphdlev_lev_lon(:,i,j,k  ) * &
          (u(:,i,j,k  ) - u(:,i,j,k-1))) / m_lon(:,i,j,k) / 2.0_r8
      end do
    end do

    do k = mesh%full_lev_ibeg + 1, mesh%full_lev_iend - 1
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          wedvdlev(:,i,j,k) = (                                        &
            wedphdlev_lev_lat(:,i,j,k+1) * (v(:,i,j,k+1) - v(:,i,j,k  )) + &
            wedphdlev_lev_lat(:,i,j,k  ) * (v(:,i,j,k  ) - v(:,i,j,k-1))   &
          ) / m_lat(:,i,j,k) / 2.0_r8
        end do
      end do
    end do
    k = mesh%full_lev_ibeg
    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        wedvdlev(:,i,j,k) = (wedphdlev_lev_lat(:,i,j,k+1) * &
          (v(:,i,j,k+1) - v(:,i,j,k  ))) / m_lat(:,i,j,k) / 2.0_r8
      end do
    end do
    k = mesh%full_lev_iend
    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        wedvdlev(:,i,j,k) = (wedphdlev_lev_lat(:,i,j,k  ) * &
          (v(:,i,j,k  ) - v(:,i,j,k-1))) / m_lat(:,i,j,k) / 2.0_r8
      end do
    end do
    end associate

  end subroutine calc_wedudlev_wedvdlev

end module operators_mod
