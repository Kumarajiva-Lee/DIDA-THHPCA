! DiDA 2.0
! letkf interface

#include "../../utils/da_config.inc"

#define DEBUG 1
#define DEBUG_PRINT 1
#define DEBUG_BARRIER 0
#define SINGLE_POINT 0

#define pi 3.1415926536

module letkf_mod

use obsmaker_mod
use redis_mod
use pro_info
use flogger
#if (SINGLE_POINT == 1)
use da_namelist_mod, only : onepoint_z
#endif
use da_namelist_mod, only : ps_pos
use namelist_mod, only : num_lev
use da_vert_coord_mod

implicit none

! ens        : ensemble number
! nx, ny     : block size in a process, nx, ny horizonal
! halo       : size of halo region
! n          : number of prognostic variables with level
integer :: ens
integer :: nx, ny, nz, halo
integer :: n

! obs_1l : first observation position (index start from 1) in a process
! obs_1g : first observation position (index start from 1) in global observation mesh 
! obs_dl : delta of grid points on latitude and longitude of observations
integer :: obs_1lat, obs_1lon
! integer :: obs_1gox, obs_1goy
integer :: obs_dlat, obs_dlon

! polar : If lat is in polar area, we will give the background back
real(8) :: polar

! inflation_factor : inflation coefficient in computation
real(8) :: inflation_factor

! pid      : process id
! global_b : start index of global grid block, x(latitude) and y(longitude)
! global_n : number of global grid blocks, x(latitude) and y(longitude)
integer :: pid
integer :: global_bx, global_by
integer :: global_nx, global_ny

! cube_n : letkf cube size
integer :: cube_n

! overlap : overlap of letkf cubes
integer :: overlap

! gc_enable : whether we will use gc filter
! neighbour_halo : whether we will use halo observations in north east / ... processes
logical :: gc_enable
logical :: neighbour_halo

! da_local_distance : localization distance
real(8) :: da_local_distance

! OBS_DENSE : whether we will use a dense strategy to choose observations
integer :: OBS_DENSE

! ndist      : number of localization distance
! local_dist(ndist) : localization distance (grid points)
integer :: ndist
integer, allocatable, dimension(:) :: local_dist

! rinv(n) : observation error
real(4), allocatable, dimension(:) :: rinv

! nlat/nlon : number of observations of current process on latitude/lontitude 
! total_obs : total observation number of current process
integer :: nlat, nlon
integer :: total_obs

! n_obs : actual size of n dimension in observation
integer :: obs_nz

! obs_l   : number of observations on a vertical strike
! obs_n   : number of observations paticipating in yb matrix
! delta_l : delta on vertical grid points
! integer :: obs_l
integer :: obs_n
integer, allocatable, dimension(:) :: obs_var_start
integer, allocatable, dimension(:) :: obs_var_lev
integer :: delta_l

! ps_pos : ps variable in vertical dimension
!integer :: ps_pos

! remote_comm : whether there should be remote communications among processes
integer :: remote_comm

! remote_obsn : remote observation number of a remote process
integer :: remote_obsn

! comm_letkf : communication group
integer :: comm_letkf

! obs_info : information of an obs item
type :: obs_info
    obs_type :: lat, lon
end type

! obs_struct : a structure for observation data
! different obs_struct for different processes
type :: obs_struct

    ! obs_id : horizonal
    ! obs_l  : vertical

    ! 2d information of a observation point (obs_id)
    type(obs_info), allocatable, dimension(:) :: info

    ! 3d information of a observation point (obs_l, obs_id)
    ayb_type, allocatable, dimension(:, :) :: ph

    ! projected observation Yb (ens, obs_nz, obs_id)
    yb_type, allocatable, dimension(:, :, :) :: pro_yb

    ! yo - ayb (obs_nz, obs_id)
    ayb_type, allocatable, dimension(:, :) :: dis_yo

    ! observation number on latitude/longitude
    ! obs_list(5) stores its total number of observations in lat_n, and lon_n = 1
    integer :: lat_n, lon_n
end type obs_struct

! 0 current process, 1 ~ 4 east/west/south/north, 5 for discrete observations
! 6 ~ 9 north east / south west / north west / south east
type(obs_struct) :: obs_list(0: 9)

! yb/ayb_for_cpy: obs data array which combines obs_list(0: 4, 6: 9)
! for better performance on sw
yb_type,  allocatable, dimension(:, :, :, :) :: yb_for_cpy
ayb_type, allocatable, dimension(:, :, :)    :: ayb_for_cpy

! neighbour id, 1 ~ 4 east/west/south/north
integer :: neighbour(0:4)

! r_neighbour : remote neighbour id (ndist, 4), 1 ~ 4 east/west/south/north
! r_sn        : start observation number of remote neighbour(1, 2 lon_n, 3, 4 lat_n)
! r_en        : end   observation number of remote neighbour(1, 2 lon_n, 3, 4 lat_n)
! rns         : remote neighbours
integer, allocatable, dimension(:, :) :: r_neighbour
integer, allocatable, dimension(:, :) :: r_sn
integer, allocatable, dimension(:, :) :: r_en
integer :: rns(4)

! !!!!!!!!!!!!!!!!!!!!!!!!!
! vars for obs_maker
type(c_ptr) :: letkf_rc
integer :: obs_vertex_id(4)
real(4), allocatable, dimension(:) :: obs_lat_t
real(4), allocatable, dimension(:) :: obs_lon_t
! !!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!

private ens, nx, ny, nz, halo, n
private obs_1lat, obs_1lon, obs_dlat, obs_dlon!, obs_1gox, obs_1goy
private polar
private pid, global_bx, global_by, global_nx, global_ny
private cube_n, overlap, gc_enable, neighbour_halo
private da_local_distance, ndist, local_dist, rinv
private OBS_DENSE
private obs_l, delta_l
private nlat, nlon, total_obs, remote_comm
private obs_list, neighbour
private r_neighbour, r_sn, r_en, rns
private letkf_rc, obs_vertex_id, obs_lat_t, obs_lon_t
!private obs_vertex_rid, obs_lat_r, obs_lon_r

contains

! #define MAX_OBS (2 * ens)
#define MAX_OBS 100000

#include "letkf_init.F90"

subroutine letkf_run(ntimes, xb, axb, xa, axa)
    
#if (DEBUG_BARRIER == 1)
    use mpi
#endif

    ! ntimes : number of da time
    integer, intent(in) :: ntimes

    ! xb  : background variation
    ! axb : average background

    xb_type,  intent(in) :: xb(ens, n, 1-halo:nx+halo, 1-halo:ny+halo)
    axb_type, intent(in) :: axb(n, 1-halo:nx+halo, 1-halo:ny+halo)

    ! xb  : analyses variation
    ! axb : average analyses

    xb_type,  intent(out) :: xa(ens, n, 1-halo:nx+halo, 1-halo:ny+halo)
    axb_type, intent(out) :: axa(n, 1-halo:nx+halo, 1-halo:ny+halo)

    integer :: i, j, k, indx
    yb_type :: debug_sum
#if (DEBUG_BARRIER == 1)
    integer :: ierr

    call mpi_barrier(comm_letkf, ierr)
#endif

#if (DEBUG_PRINT == 1)
    call log_notice("Letkf run start!")
#endif

#if (DEBUG == 1)
    call GPTLstart("Letkf run")
    call GPTLstart("Obslist generation")
#endif
    if (ps_pos == 0) then
        call gen_obslist(ntimes, xb, axb)
    else
        call gen_obslist_3d(ntimes, xb, axb)
    end if
#if (DEBUG == 1)
    call GPTLstop("Obslist generation")
#endif

#if (DEBUG_BARRIER == 1)
    call mpi_barrier(comm_letkf, ierr)
#endif

#if (DEBUG == 1)
    call GPTLstart("Obslist communication")
#endif
    call comm_obslist()
    call gen_obs_for_cpy()
#if (DEBUG == 1)
    call GPTLstop("Obslist communication")
#endif

#if (DEBUG_PRINT == 1)
    call log_notice("Obslist generation & communication done!")
#endif

#if (DEBUG_BARRIER == 1)
    call mpi_barrier(comm_letkf, ierr)
#endif

#if (DEBUG == 1)
    call GPTLstart("Letkf compute")
#endif
    if (overlap == 0) then
        call letkf_compute_nooverlap(xb, axb, xa, axa)
    else
        call letkf_compute(xb, axb, xa, axa)
    end if
#if (DEBUG == 1)
    call GPTLstop("Letkf compute")
#endif

#if (DEBUG_PRINT == 1)
    call log_notice("Letkf run done!")
#endif

#if (DEBUG == 1)
    call GPTLstop("Letkf run")
    call GPTLpr(pid)
    call GPTLpr_summary(comm_letkf)
#endif

#if (DEBUG_BARRIER == 1)
    call mpi_barrier(comm_letkf, ierr)
#endif

end subroutine

#include "letkf_compute.F90"

#include "obs_yb.F90"

#include "obs_list.F90"

#include "obs_make.F90"

#undef MAX_OBS

end module

#undef SINGLE_POINT
#undef DEBUG
#undef pi
