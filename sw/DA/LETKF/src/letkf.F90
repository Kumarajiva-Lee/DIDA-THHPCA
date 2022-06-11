! LETKF 3.0

#include "../../utils/da_config.inc"

#define DEBUG 1
#define DEBUG_PRINT 1
#define DEBUG_BARRIER 1
#define SINGLE_POINT 0

module letkf_mod

use read_obs_mod
use pro_info
use flogger
use mpal

implicit none

private

! ens        : ensemble number
! nx, ny     : block size in a process, nx(lat), ny(lon) horizonal
! halo       : size of halo region
! n          : number of prognostic variables with level
integer :: ens
integer :: nx, ny, nz, halo
integer :: n

! inflation_factor : inflation coefficient in computation
real(8) :: inflation_factor

! local_dist : localization distances
real(8) :: local_dist

! vert_dist : vertical localization distances
real(8) :: vert_dist

! pid      : process id
! global_b : start index of global grid block, x(latitude) and y(longitude)
! global_n : number of global grid blocks, x(latitude) and y(longitude)
! global_g : number of global grids, x(latitude) and y(longitude)
integer :: pid
integer :: global_bx, global_by
integer :: global_nx, global_ny
integer :: global_gx, global_gy

! south_n    : south direction processes
! north_n    : north direction processes
integer :: south_n, north_n

! obs_type_num : number of types of observations
! now in read_obs_mod.F90
! integer :: obs_type_num

type :: obs_latcircle

    integer :: east_n, west_n

    ! obs_data : array of obs_struct (obs_type_num, -east_n: west_n)
    type(obs_list), allocatable, dimension(:, :) :: obs_data

    ! neighbour_id : communication id (-east_n: west_n)
    integer, allocatable, dimension(:) :: nid

end type

! obs_lc : observation data in different lat circle (-north_n, south_n)
type(obs_latcircle), allocatable, dimension(:) :: obs_lc

! gc_filter : whether we will use gc filter
logical :: gc_filter

! ps_indx : variable position in n
integer :: ps_indx

! obs_polar : discard obs in polar
real(8) :: obs_polar

! n_var :: number of da variables
integer :: n_var
! da_var_list (n_var): da variable name
character(len = 20),  dimension(10) :: da_var_list
! da_var_start(n_var): da variable start point in n dimension
! da_var_level(n_var): da variable level 
integer, allocatable, dimension(:) :: da_var_start
integer, allocatable, dimension(:) :: da_var_lev

! cube_n : letkf horizonal cube size
integer :: cube_n
! cube_nz : letkf vertical cube size
integer :: cube_nz
! overlap : overlap of letkf cubes
integer :: overlap

! precision : computing precision in letkf core
integer :: precision

! ! obs_info (obs_type_num)
! type(obs_information), allocatable, dimension(:) :: obs_info

! ! observation information
! type :: obs_information

!     ! zonal level of observation nodes
!     integer :: obs_l
!     ! number of observation variables
!     integer :: obs_v
!     ! obs_nz = obs_l * obs_v
!     ! actual size of n dimension in observation
!     integer :: obs_nz

!     ! observation error rinv(obs_v)
!     real(8), allocatable, dimension(:) :: rinv

!     ! obs_vid(obs_v)
!     ! observation variable id in background variables
!     integer, allocatable, dimension(:) :: obs_vid

! end type

! ! obs_struct : structure for observation data
! ! one obs_struct corresponds to one observation data type
! type :: obs_struct

!     ! use obs_id to get obs_info
!     integer :: obs_id

!     ! number of observation nodes
!     integer :: obs_num

!     ! horizonal position lat/lon(obs_num)
!     real(8), allocatable, dimension(:) :: lat, lon

!     ! ph value of a observation point (obs_nz, obs_num)
!     real(8), allocatable, dimension(:, :) :: ph

!     ! projected observation Yb (ens, obs_nz, obs_num)
!     yb_type, allocatable, dimension(:, :, :) :: pro_yb

!     ! yo - ayb (obs_nz, obs_num)
!     ayb_type, allocatable, dimension(:, :) :: dis_yo

! end type

! comm_letkf : communication group
integer :: comm_letkf

! redis rc
type(c_ptr) :: letkf_rc

public letkf_init, letkf_run

contains

#define MAX_OBS 8192

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
        call gen_obs_struct(ntimes, xb, axb)
    
#if (DEBUG_BARRIER == 1)
        call mpi_barrier(comm_letkf, ierr)
#endif

#if (DEBUG == 1)
        call GPTLstop("Obslist generation")
#endif
    
#if (DEBUG == 1)
        call GPTLstart("Communication and computation")
        call GPTLstart("Obslist communication")
#endif
        call comm_obs_data()
        ! call comm_obslist()
        ! call gen_obs_for_cpy()
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
            if (vert_dist <= 0d0) then
                call letkf_compute_nooverlap(xb, axb, xa, axa)
            else
                call letkf_compute_nooverlap_vertical(xb, axb, xa, axa)
            end if
        else
            call letkf_compute(xb, axb, xa, axa)
        end if
#if (DEBUG == 1)
        call GPTLstop("Letkf compute")
        call GPTLstop("Communication and computation")
#endif

        call letkf_final()
    
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
#include "letkf_init.F90"

#include "letkf_obs.F90"

#include "obs_make.F90"

#include "letkf_yb.F90"

#include "letkf_compute.F90"

#include "letkf_final.F90"

#undef MAX_OBS

end module
