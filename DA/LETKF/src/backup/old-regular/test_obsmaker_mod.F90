#include "letkf_config.inc"

module test_obsmaker_mod
use obsmaker_mod
use redis_mod
use coupler_config, only:redis_address
implicit none

contains
  subroutine test_obsmaker(lonids, lonide, latids, latide)
    type(c_ptr) :: rc
    integer, intent(in) :: lonids, lonide, latids, latide
    integer :: obs_lat_num, obs_lon_num
    integer, dimension(4) :: obs_vertex_id
    real(4), allocatable :: obs_lon(:), obs_lat(:)
    real(4), allocatable :: obs_r(:)
    real(4), allocatable :: obs_p(:,:)
    real(4), allocatable :: yo(:,:)
    integer :: ii, jj
    character(12) :: date_char

    rc = RedisConnect(redis_address)
    call obs_init()

    call obs_get_num(latids, latide, lonids, lonide, obs_lat_num, obs_lon_num, obs_vertex_id)

    allocate(obs_lon(obs_lon_num))
    allocate(obs_lat(obs_lat_num))

    call obs_get_loc(rc, obs_lon_num, obs_lat_num, obs_vertex_id, obs_lon, obs_lat)

    allocate(obs_r(obs_num_2d+obs_num_3d*obs_l))

    allocate(obs_p(obs_l, obs_lat_num*obs_lon_num))
    allocate(yo(obs_num_2d+obs_num_3d*obs_l, obs_lat_num*obs_lon_num))

    call obs_get_r(rc, obs_r)
    date_char='200001010000'
    call obs_get_yo(rc, obs_vertex_id, obs_lat_num, obs_lon_num, date_char, obs_r, obs_p, yo)

    open(1,file='obs_p.dat')
    open(2,file='yo.dat')

130     format(20f16.8)
!
    write(1,130) ((obs_p(ii,jj),jj=1,obs_lat_num*obs_lon_num),ii=1,obs_l)
    write(2,130) ((yo(ii,jj),jj=1,obs_lat_num*obs_lon_num),ii=1,obs_l*obs_num_3d+obs_num_2d)

    close(1)
    close(2)

    deallocate(obs_lon)
    deallocate(obs_lat)
    deallocate(obs_r)
    deallocate(obs_p)
    deallocate(yo)
  end subroutine test_obsmaker
end module test_obsmaker_mod
