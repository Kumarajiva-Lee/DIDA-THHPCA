module parallel_types_mod

  use mpi

  implicit none

  private

  public async_type

  ! operators_mod中用到的
  integer, public, parameter :: async_u                       = 1
  integer, public, parameter :: async_u_f                     = 2
  integer, public, parameter :: async_v                       = 3
  integer, public, parameter :: async_v_f                     = 4
  integer, public, parameter :: async_wedphdlev_lev           = 5
  integer, public, parameter :: async_wedphdlev_lev_lon       = 6
  integer, public, parameter :: async_wedphdlev_lev_lat       = 7
  integer, public, parameter :: async_gz                      = 8
  integer, public, parameter :: async_gz_f                    = 9
  integer, public, parameter :: async_gz_lev                  = 10
  integer, public, parameter :: async_m                       = 11
  integer, public, parameter :: async_m_vtx                   = 12
  integer, public, parameter :: async_m_lon                   = 13
  integer, public, parameter :: async_m_lat                   = 14
  integer, public, parameter :: async_mf_lon_n                = 15
  integer, public, parameter :: async_mf_lat_n                = 16
  integer, public, parameter :: async_mf_lat_t                = 17
  integer, public, parameter :: async_mf_lon_t                = 18
  integer, public, parameter :: async_pv                      = 19
  integer, public, parameter :: async_pv_lon                  = 20
  integer, public, parameter :: async_pv_lat                  = 21
  integer, public, parameter :: async_ke                      = 22    
  integer, public, parameter :: async_pt                      = 23
  integer, public, parameter :: async_pt_f                    = 91
  integer, public, parameter :: async_pt_lon                  = 24  
  integer, public, parameter :: async_pt_lat                  = 25
  integer, public, parameter :: async_pt_lev                  = 26
  integer, public, parameter :: async_t                       = 27
  integer, public, parameter :: async_ph                      = 28
  integer, public, parameter :: async_ph_lev                  = 29  
  integer, public, parameter :: async_phs                     = 30
  integer, public, parameter :: async_phs_f                   = 31
  integer, public, parameter :: async_div                     = 32  
  integer, public, parameter :: async_div2                    = 33  
  integer, public, parameter :: async_vor                     = 34
  integer, public, parameter :: async_m_lev                   = 35  
  integer, public, parameter :: async_wedphdlev               = 36
  integer, public, parameter :: async_w                       = 37
  integer, public, parameter :: async_w_lev                   = 38
  integer, public, parameter :: async_w_lev_lon               = 39
  integer, public, parameter :: async_w_lev_lat               = 40
  integer, public, parameter :: async_gz_lev_lon              = 41
  integer, public, parameter :: async_gz_lev_lat              = 42
  integer, public, parameter :: async_rhod                    = 43
  integer, public, parameter :: async_rhod_lon                = 44 
  integer, public, parameter :: async_rhod_lat                = 45
  integer, public, parameter :: async_p                       = 46
  integer, public, parameter :: async_p_lev                   = 47  
  integer, public, parameter :: async_p_lev_lon               = 48
  integer, public, parameter :: async_p_lev_lat               = 49
  integer, public, parameter :: async_u_lev_lon               = 50
  integer, public, parameter :: async_v_lev_lat               = 51
  integer, public, parameter :: async_mf_lev_lon_n            = 52  
  integer, public, parameter :: async_mf_lev_lat_n            = 53
  integer, public, parameter :: async_tension_h               = 54
  integer, public, parameter :: async_shear_h                 = 55
  integer, public, parameter :: async_kmh                     = 56
  integer, public, parameter :: async_kmh_vtx                 = 57
  integer, public, parameter :: async_kmh_lon                 = 58
  integer, public, parameter :: async_kmh_lat                 = 59
  integer, public, parameter :: async_landmask                = 60
  integer, public, parameter :: async_gzs                     = 61
  integer, public, parameter :: async_zs_std                  = 62
  integer, public, parameter :: async_dzsdlon                 = 63
  integer, public, parameter :: async_dzsdlat                 = 64
  integer, public, parameter :: async_du                      = 65
  integer, public, parameter :: async_dv                      = 66  
  integer, public, parameter :: async_dgz                     = 67
  integer, public, parameter :: async_dpt                     = 68 
  integer, public, parameter :: async_dphs                    = 69
  integer, public, parameter :: async_qhv                     = 70
  integer, public, parameter :: async_qhu                     = 71  
  integer, public, parameter :: async_dkedlon                 = 72
  integer, public, parameter :: async_dkedlat                 = 73
  integer, public, parameter :: async_dmfdlon                 = 74
  integer, public, parameter :: async_dmfdlat                 = 75
  integer, public, parameter :: async_dptfdlon                = 76
  integer, public, parameter :: async_dptfdlat                = 77  
  integer, public, parameter :: async_dptfdlev                = 78
  integer, public, parameter :: async_pgf_lon                 = 79  
  integer, public, parameter :: async_pgf_lat                 = 80
  integer, public, parameter :: async_wedudlev                = 81
  integer, public, parameter :: async_wedvdlev                = 82
  integer, public, parameter :: async_smag_dpt                = 83
  integer, public, parameter :: async_adv_gz_lon              = 84    
  integer, public, parameter :: async_adv_gz_lat              = 85
  integer, public, parameter :: async_adv_gz_lev              = 86
  integer, public, parameter :: async_adv_w_lon               = 87  
  integer, public, parameter :: async_adv_w_lat               = 88
  integer, public, parameter :: async_adv_w_lev               = 89
  integer, public, parameter :: async_tmp                     = 90

  integer, public, parameter :: async_total_num               = 91

  type async_type
    integer, allocatable :: send_req(:)
    integer, allocatable :: recv_req(:)
  contains
    procedure :: init => async_init
    procedure :: wait => async_wait
    final :: async_final
  end type async_type

contains

  subroutine async_init(this, num_send_ngb , num_recv_ngb)

    class(async_type), intent(inout) :: this
    integer, intent(in) :: num_send_ngb
    integer, intent(in) :: num_recv_ngb

    if (allocated(this%send_req)) deallocate(this%send_req)
    if (allocated(this%recv_req)) deallocate(this%recv_req)

    allocate(this%send_req(num_send_ngb)); this%send_req = MPI_REQUEST_NULL
    allocate(this%recv_req(num_recv_ngb)); this%recv_req = MPI_REQUEST_NULL

  end subroutine async_init

  subroutine async_wait(this)

    class(async_type), intent(inout) :: this

    integer i, ierr

    ! do i = 1, size(this%send_req)
    !   if (this%send_req(i) /= MPI_REQUEST_NULL) then
    !     call MPI_WAIT(this%send_req(i), MPI_STATUS_IGNORE, ierr)
    !   end if
    !   if (this%recv_req(i) /= MPI_REQUEST_NULL) then
    !     call MPI_WAIT(this%recv_req(i), MPI_STATUS_IGNORE, ierr)
    !   end if
    ! end do

    do i = 1, size(this%send_req)
      call MPI_WAIT(this%send_req(i), MPI_STATUS_IGNORE, ierr)
    end do
    do i = 1, size(this%recv_req)
      call MPI_WAIT(this%recv_req(i), MPI_STATUS_IGNORE, ierr)
    end do

  end subroutine async_wait

  subroutine async_final(this)

    type(async_type), intent(inout) :: this

    if (allocated(this%send_req)) deallocate(this%send_req)
    if (allocated(this%recv_req)) deallocate(this%recv_req)

  end subroutine async_final

end module parallel_types_mod
