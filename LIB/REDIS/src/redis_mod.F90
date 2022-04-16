module redis_mod
  use iso_c_binding

  implicit none

  character(len=1), parameter :: nullchar=char(0)

  interface RedisHget
      module procedure RedisHgeti
      module procedure RedisHgetf
      module procedure RedisHgetd
      module procedure RedisHgets
  end interface RedisHget 

  interface RedisHset
      module procedure RedisHseti
      module procedure RedisHsetf
      module procedure RedisHsetd
      module procedure RedisHsets
  end interface RedisHset 

  !interface RedisHmsetf
  !    module procedure RedisHmsetf3d
  !    module procedure RedisHmsetf2d
  !end interface RedisHmsetf
  !interface RedisCmd
  !    module procedure RedisCmdreturn
  !    module procedure RedisCmdreturnf
  !end interface RedisCmd

  interface

     subroutine err_print() bind(c, name='err_to_stdout') 
        use iso_c_binding
     end subroutine err_print
     
     function usleep(n) bind(c)
       use iso_c_binding
       integer(c_int), value :: n
       integer(c_int) usleep
     end function usleep

     !--------------------------------------------------------
     !  init part
     !--------------------------------------------------------
     subroutine redis_test_set(c) bind(c, name='redisTestSet')
       use iso_c_binding
       type(c_ptr), value :: c
     end subroutine redis_test_set

     type(c_ptr) function RedisConnect(redis_address) &
             bind(c, name='setupConnection')
       use iso_c_binding
       character(c_char)     :: redis_address(*)
     end function RedisConnect

     subroutine RedisDisconnect(c) bind(c, name='redisClusterFree')
       use iso_c_binding
       type(c_ptr), value :: c
     end subroutine RedisDisconnect

     integer(c_int) function set_redis_socket_timeout_f(cc, sec) &
             bind(c, name='set_redis_socket_timeout')
       use iso_c_binding
       type(c_ptr), value    :: cc
       integer(c_int)        :: sec
     end function set_redis_socket_timeout_f

     integer(c_int) function redis_cmd_f(c, cmd) bind(c, name='redis_cmd')
       use iso_c_binding
       type(c_ptr), value    :: c
       character(c_char)     :: cmd(*)
     end function redis_cmd_f

     integer(c_int) function redis_cmdreturnf_f(c, cmd, buf) bind(c, name='redis_cmdreturnf')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: cmd(*)
     end function redis_cmdreturnf_f

     integer(c_int) function redis_cmdreturnd_f(c, cmd, buf) bind(c, name='redis_cmdreturnd')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: cmd(*)
     end function redis_cmdreturnd_f

     integer(c_int) function redis_cmdreturni_f(c, cmd, buf) bind(c, name='redis_cmdreturni')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: cmd(*)
     end function redis_cmdreturni_f

     !--------------------------------------------------------------------------
     !  string part
     !-------------------------------------------------------------------------

     integer(c_int) function redis_set_f(c, key, val) bind(c, name='redis_set')
       use iso_c_binding
       type(c_ptr), value    :: c
       character(c_char)     :: val(*), key(*)
     end function redis_set_f

     integer(c_int) function redis_setb_f(c, key, buf, len) bind(c, name='redis_setb')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: key(*)
       integer(c_int), intent(in) ::len
     end function redis_setb_f

     !----------------------------------------------------------------------------
     !  hash part
     !----------------------------------------------------------------------------
     integer(c_int) function RedisHgeti_f(c, hash, key, buf) bind(c, name='RedisHgeti')
       use iso_c_binding
       type(c_ptr), value    :: c
       character(c_char)     :: hash(*), key(*)
       integer(c_int)        :: buf
     end function RedisHgeti_f

     integer(c_int) function RedisHgetf_f(c, hash, key, buf) bind(c, name='RedisHgetf')
       use iso_c_binding
       type(c_ptr), value    :: c
       character(c_char)     :: hash(*), key(*)
       real(c_float)                :: buf
     end function RedisHgetf_f

     integer(c_int) function RedisHgetd_f(c, hash, key, buf) bind(c, name='RedisHgetd')
       use iso_c_binding
       type(c_ptr), value    :: c
       character(c_char)     :: hash(*), key(*)
       real(c_double)                :: buf
     end function RedisHgetd_f

     integer(c_int) function RedisHgets_f(c, hash, key, buf, len1, len2, len3) bind(c, name='RedisHgets')
       use iso_c_binding
       type(c_ptr), value    :: c
       integer(c_int)        :: len1, len2, len3
       character(c_char)     :: hash(*), key(*), buf(*)
     end function RedisHgets_f

     integer(c_int) function RedisHseti_f(c, hash, key, buf) bind(c, name='RedisHseti')
       use iso_c_binding
       type(c_ptr), value    :: c
       character(c_char)     :: hash(*), key(*)
       integer(c_int)        :: buf
     end function RedisHseti_f

     integer(c_int) function RedisHsetf_f(c, hash, key, buf) bind(c, name='RedisHsetf')
       use iso_c_binding
       type(c_ptr), value    :: c
       character(c_char)     :: hash(*), key(*)
       real(c_float)                :: buf
     end function RedisHsetf_f

     integer(c_int) function RedisHsetd_f(c, hash, key, buf) bind(c, name='RedisHsetd')
       use iso_c_binding
       type(c_ptr), value    :: c
       character(c_char)     :: hash(*), key(*)
       real(c_double)                :: buf
     end function RedisHsetd_f

     integer(c_int) function RedisHsets_f(c, hash, key, buf) bind(c, name='RedisHsets')
       use iso_c_binding
       type(c_ptr), value    :: c
       character(c_char)     :: hash(*), key(*), buf(*)
     end function RedisHsets_f

     !------------------------------------------------------------------
     !  zset part
     !------------------------------------------------------------------

      !--------------------------------------------------------------------------
      !  interface for non_regular observation part
      !-------------------------------------------------------------------------

      !zadd & hadd(4)
      ! zset <attr>:hkey index score
      ! hset <attr>:hkey 
      ! n是变量个数
      integer(c_int) function redis_setobserve_f(c, hkey, varlist, lon_buf, lat_buf, &
              lev_buf, val_buf, obs_num, n) bind(c, name='redis_setobservef')
        use iso_c_binding
        type(c_ptr), value :: c, lon_buf, lat_buf, lev_buf, val_buf
        character(c_char)  :: hkey(*), varlist(*)
        integer(c_int), intent(in) :: obs_num, n
      end function redis_setobserve_f

      integer(c_int) function redis_setobserve_d(c, hkey, varlist, lon_buf, lat_buf, &
              lev_buf, val_buf, obs_num, n) bind(c, name='redis_setobserved')
        use iso_c_binding
        type(c_ptr), value :: c, lon_buf, lat_buf, lev_buf, val_buf
        character(c_char)  :: hkey(*), varlist(*)
        integer(c_int), intent(in) :: obs_num, n
      end function redis_setobserve_d

      ! zcount <attr>:hkey
      integer(c_int) function redis_obscount(c, hkey, its, ite, jts, jte, buf) bind(c, name='redis_obscount')
        use iso_c_binding
        type(c_ptr), value :: c
        character(c_char) :: hkey(*)
        real(c_float), intent(in) :: its, ite, jts, jte
        integer(c_int) :: buf
      end function redis_obscount

      ! zrange <attr>:hkey
      ! hget
      integer(c_int) function redis_getobserve_f(c, hkey, varlist, its, ite, jts, jte, &
              lon_buf, lat_buf, lev_buf, val_buf, n) bind(c, name='redis_getobservef')
        use iso_c_binding
        type(c_ptr), value :: c, lon_buf, lat_buf, lev_buf, val_buf
        character(c_char) :: hkey(*), varlist(*)
        real(c_float), intent(in) :: its, ite, jts, jte
        integer(c_int), intent(in) :: n
      end function redis_getobserve_f

      integer(c_int) function redis_getobserve_d(c, hkey, varlist, its, ite, jts, jte, &
              lon_buf, lat_buf, lev_buf, val_buf, n) bind(c, name='redis_getobserved')
        use iso_c_binding
        type(c_ptr), value :: c, lon_buf, lat_buf, lev_buf, val_buf
        character(c_char) :: hkey(*), varlist(*)
        real(c_double), intent(in) :: its, ite, jts, jte
        integer(c_int), intent(in) :: n
      end function redis_getobserve_d

     !------------------------------------------------------------------
     !  custum part
     !-----------------------------------------------------------------
     integer(c_int) function redis_hmsetf_f(c, hkey, varlist, its, ite, istride, jts, jte, &
             jstride, kms, kme, kstride, n, buf) bind(c, name='redis_hmsetf')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: its, ite, jts, jte, kms, kme, n
       integer(c_int), intent(in) :: istride, jstride, kstride
     end function redis_hmsetf_f

    integer(c_int) function redis_hmsetf2d_f(c, hkey, varlist, its, ite, istride, jts, jte, &
             jstride, n, buf) bind(c, name='redis_hmsetf2d')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: its, ite, jts, jte, n
       integer(c_int), intent(in) :: istride, jstride
     end function redis_hmsetf2d_f

    integer(c_int) function redis_hmsetf_f_nonblock(c, hkey, varlist, its, ite, istride, jts, jte, &
             jstride, kms, kme, kstride, n, buf) bind(c, name='redis_hmsetf_nonblock')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: its, ite, jts, jte, kms, kme, n
       integer(c_int), intent(in) :: istride, jstride, kstride
     end function redis_hmsetf_f_nonblock

    integer(c_int) function redis_hmsetf2d_f_nonblock(c, hkey, varlist, its, ite, istride, jts, jte, &
             jstride, n, buf) bind(c, name='redis_hmsetf2d_nonblock')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: its, ite, jts, jte, n
       integer(c_int), intent(in) :: istride, jstride
     end function redis_hmsetf2d_f_nonblock

     integer(c_int) function redis_hmsetf1d_f(c, hkey, varlist, its, ite, istride, n, buf) &
             bind(c, name='redis_hmsetf1d')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: its, ite, n
       integer(c_int), intent(in) :: istride
     end function redis_hmsetf1d_f

     integer(c_int) function redis_hmgetf1d_f(c, hkey, varlist,  its, ite, istride, n, buf) &
             bind(c, name='redis_hmgetf1d')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: its, ite, n
       integer(c_int), intent(in) :: istride
     end function redis_hmgetf1d_f

     integer(c_int) function redis_hmgetf2d_f(c, hkey, varlist,  its, ite, istride, jts, jte, &
             jstride, n, buf) bind(c, name='redis_hmgetf2d')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: its, ite, jts, jte, n
       integer(c_int), intent(in) :: istride, jstride
     end function redis_hmgetf2d_f

     integer(c_int) function redis_hmgetf_f(c, hkey, varlist,  its, ite, istride, jts, jte, &
             jstride, kms, kme, kstride, n, buf) bind(c, name='redis_hmgetf')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: its, ite, jts, jte, kms, kme, n
       integer(c_int), intent(in) :: istride, jstride, kstride
     end function redis_hmgetf_f

     integer(c_int) function redis_hmgetf2d_f_nonblock(c, hkey, varlist,  its, ite, istride, jts, jte, &
             jstride, n, buf) bind(c, name='redis_hmgetf2d_nonblock')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: its, ite, jts, jte, n
       integer(c_int), intent(in) :: istride, jstride
     end function redis_hmgetf2d_f_nonblock

     integer(c_int) function redis_hmgetf_f_nonblock(c, hkey, varlist,  its, ite, istride, jts, jte, &
             jstride, kms, kme, kstride, n, buf) bind(c, name='redis_hmgetf_nonblock')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: its, ite, jts, jte, kms, kme, n
       integer(c_int), intent(in) :: istride, jstride, kstride
     end function redis_hmgetf_f_nonblock

     integer(c_int) function redis_hmsetd1d_f(c, hkey, varlist,  its, ite, istride, n, buf) &
             bind(c, name='redis_hmsetd1d')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: its, ite, n
       integer(c_int), intent(in) :: istride
     end function redis_hmsetd1d_f

     integer(c_int) function redis_hmsetd2d_f(c, hkey, varlist,  its, ite, istride, jts, jte, &
             jstride, n, buf) bind(c, name='redis_hmsetd2d')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: its, ite, jts, jte, n
       integer(c_int), intent(in) :: istride, jstride
     end function redis_hmsetd2d_f

     integer(c_int) function redis_hmsetd_f(c, hkey, varlist,  its, ite, istride, jts, jte,&
             jstride, kms, kme, kstride, n, buf) bind(c, name='redis_hmsetd')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: its, ite, jts, jte, kms, kme, n
       integer(c_int), intent(in) :: istride, jstride, kstride
     end function redis_hmsetd_f

     integer(c_int) function redis_hmsetd2d_f_nonblock(c, hkey, varlist,  its, ite, istride, jts, jte, &
             jstride, n, buf) bind(c, name='redis_hmsetd2d_nonblock')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: its, ite, jts, jte, n
       integer(c_int), intent(in) :: istride, jstride
     end function redis_hmsetd2d_f_nonblock

     integer(c_int) function redis_hmsetd_f_nonblock(c, hkey, varlist,  its, ite, istride, jts, jte,&
             jstride, kms, kme, kstride, n, buf) bind(c, name='redis_hmsetd_nonblock')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: its, ite, jts, jte, kms, kme, n
       integer(c_int), intent(in) :: istride, jstride, kstride
     end function redis_hmsetd_f_nonblock

     integer(c_int) function redis_hmgetd1d_f(c, hkey, varlist,  its, ite, istride, n, buf) &
             bind(c, name='redis_hmgetd1d')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: its, ite, n
       integer(c_int), intent(in) :: istride
     end function redis_hmgetd1d_f

     integer(c_int) function redis_hmgetd2d_f(c, hkey, varlist,  its, ite, istride, jts, jte, &
             jstride, n, buf) bind(c, name='redis_hmgetd2d')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: its, ite, jts, jte, n
       integer(c_int), intent(in) :: istride, jstride
     end function redis_hmgetd2d_f

     integer(c_int) function redis_hmgetd_f(c, hkey, varlist,  its, ite, istride, jts, jte, &
             jstride, kms, kme, kstride, n, buf) bind(c, name='redis_hmgetd')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: its, ite, jts, jte, kms, kme, n
       integer(c_int), intent(in) :: istride, jstride, kstride
     end function redis_hmgetd_f

     integer(c_int) function redis_hmgetd2d_f_nonblock(c, hkey, varlist,  its, ite, istride, jts, jte, &
             jstride, n, buf) bind(c, name='redis_hmgetd2d_nonblock')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: its, ite, jts, jte, n
       integer(c_int), intent(in) :: istride, jstride
     end function redis_hmgetd2d_f_nonblock

     integer(c_int) function redis_hmgetd_f_nonblock(c, hkey, varlist,  its, ite, istride, jts, jte, &
             jstride, kms, kme, kstride, n, buf) bind(c, name='redis_hmgetd_nonblock')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: its, ite, jts, jte, kms, kme, n
       integer(c_int), intent(in) :: istride, jstride, kstride
     end function redis_hmgetd_f_nonblock

     !------------------------------------------------------------------------------
     ! da_letkf_customize
     !-----------------------------------------------------------------------------

     integer(c_int) function redis_da_outputd_f(c, hkey, varlist,  lonids, lonide, lonstep, latids, latide, latstep, &
             mpas_num_lev, num_lev, zstep, num_2d, num_3d, buf) bind(c, name='redis_da_outputd')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: lonids, lonide, lonstep, latids, latide, &
           latstep, mpas_num_lev, num_lev, zstep, num_2d, num_3d
     end function

     integer(c_int) function redis_da_outputf_f(c, hkey, varlist,  lonids, lonide, lonstep, latids, latide, latstep, &
             mpas_num_lev, num_lev, zstep, num_2d, num_3d, buf) bind(c, name='redis_da_outputf')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: lonids, lonide, lonstep, latids, latide, &
           latstep, mpas_num_lev, num_lev, zstep, num_2d, num_3d
     end function

     integer(c_int) function redis_da_outputd_f_nonblock(c, hkey, varlist,  lonids, lonide, lonstep, latids, latide, latstep, &
             mpas_num_lev, num_lev, zstep, num_2d, num_3d, buf) bind(c, name='redis_da_outputd_nonblock')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: lonids, lonide, lonstep, latids, latide, &
           latstep, mpas_num_lev, num_lev, zstep, num_2d, num_3d
     end function

     integer(c_int) function redis_da_outputf_f_nonblock(c, hkey, varlist,  lonids, lonide, lonstep, latids, latide, latstep, &
             mpas_num_lev, num_lev, zstep, num_2d, num_3d, buf) bind(c, name='redis_da_outputf_nonblock')
       use iso_c_binding
       type(c_ptr), value    :: c, buf
       character(c_char)     :: hkey(*) , varlist(*)
       integer(c_int), intent(in) :: lonids, lonide, lonstep, latids, latide, &
           latstep, mpas_num_lev, num_lev, zstep, num_2d, num_3d
     end function

  end interface
contains

  !------------------------------------------------------------------------------------
  !  init part
  !------------------------------------------------------------------------------------
  subroutine RedisCmd(c, cmd)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    character(len=*), intent(in)             :: cmd

    errno = redis_cmd_f(c, c_str(cmd))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine RedisCmdreturn(c, cmd, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    character(len=*), intent(in)             :: cmd
    type(*), target :: buf(*)

    errno = redis_cmdreturnf_f(c, c_str(cmd), c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine RedisCmdreturnd(c, cmd, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    character(len=*), intent(in)             :: cmd
    type(*), target :: buf(*)

    errno = redis_cmdreturnd_f(c, c_str(cmd), c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine RedisCmdreturni(c, cmd, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    character(len=*), intent(in)             :: cmd
    type(*), target :: buf(*)

    errno = redis_cmdreturni_f(c, c_str(cmd), c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  !------------------------------------------------------------------------------------
  !  string part
  !------------------------------------------------------------------------------------
  subroutine redis_set(c, key, val)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    character(len=*), intent(in)             :: key, val

    errno = redis_set_f(c, c_str(key), c_str(val))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine redis_set

  subroutine redis_setb(c, key, buf, len)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    character(len=*), intent(in)             :: key
    type(*), target :: buf(*)
    integer, intent(in) :: len

    errno = redis_setb_f(c, c_str(key), c_loc(buf), len)
    if (errno /= 0) then
        call err_print
    end if
  end subroutine redis_setb

  !------------------------------------------------------------------------------------
  !  hash part
  !------------------------------------------------------------------------------------
  subroutine RedisHgets(c, hash, key, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    type(c_ptr)                              :: c_string
    character(len=*), intent(in)             :: hash, key
    character(len=*), intent(out)            :: buf
    integer                                  :: len1, len2, len3

    errno = RedisHgets_f(c, c_str(hash), c_str(key), buf, len1, len2, len3)
    if (errno /= 0) then
        call err_print
    end if

  end subroutine RedisHgets

  subroutine RedisHgeti(c, hash, key, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    character(len=*), intent(in)             :: hash, key
    integer(c_int)  :: buf

    errno = RedisHgeti_f(c, c_str(hash), c_str(key), buf)
    if (errno /= 0) then
        call err_print
    end if
  end subroutine RedisHgeti

  subroutine RedisHgetf(c, hash, key, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    character(len=*), intent(in)             :: hash, key
    real(c_float)  :: buf

    errno = RedisHgetf_f(c, c_str(hash), c_str(key), buf)
    if (errno /= 0) then
        call err_print
    end if
  end subroutine RedisHgetf

  subroutine RedisHgetd(c, hash, key, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    character(len=*), intent(in)             :: hash, key
    real(c_double)  :: buf

    errno = RedisHgetd_f(c, c_str(hash), c_str(key), buf)
    if (errno /= 0) then
        call err_print
    end if
  end subroutine RedisHgetd

  subroutine RedisHsets(c, hash, key, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    type(c_ptr)                              :: c_string
    character(len=*), intent(in)             :: hash, key
    character(len=*), intent(in)            :: buf

    errno = RedisHsets_f(c, c_str(hash), c_str(key), c_str(buf))
    if (errno /= 0) then
        call err_print
    end if

  end subroutine RedisHsets

  subroutine RedisHseti(c, hash, key, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    character(len=*), intent(in)             :: hash, key
    integer(c_int)  :: buf

    errno = RedisHseti_f(c, c_str(hash), c_str(key), buf)
    if (errno /= 0) then
        call err_print
    end if
  end subroutine RedisHseti

  subroutine RedisHsetf(c, hash, key, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    character(len=*), intent(in)             :: hash, key
    real(c_float)  :: buf

    errno = RedisHsetf_f(c, c_str(hash), c_str(key), buf)
    if (errno /= 0) then
        call err_print
    end if
  end subroutine RedisHsetf

  subroutine RedisHsetd(c, hash, key, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    character(len=*), intent(in)             :: hash, key
    real(c_double)  :: buf

    errno = RedisHsetd_f(c, c_str(hash), c_str(key), buf)
    if (errno /= 0) then
        call err_print
    end if
  end subroutine RedisHsetd

  !------------------------------------------------------------------------------------
  !  zset part
  !------------------------------------------------------------------------------------

  ! non_regular observation interface

  subroutine RedisSetobserveF(c, hkey, varlist, lon_buf, lat_buf, lev_buf, &
          val_buf, obs_num, n)
    type(c_ptr), intent(in)                  :: c
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: lon_buf(*), lat_buf(*), lev_buf(*), val_buf(*)
    integer(c_int), intent(in)               :: obs_num, n
    integer(c_int)                           :: errno

    errno = redis_setobserve_f(c, c_str(hkey), c_str(varlist), c_loc(lon_buf), &
        c_loc(lat_buf), c_loc(lev_buf), c_loc(val_buf), obs_num, n)
    if (errno /= 0) then
      call err_print
    end if
  end subroutine

  subroutine RedisSetobserveD(c, hkey, varlist, lon_buf, lat_buf, lev_buf, &
    val_buf, obs_num, n)
    type(c_ptr), intent(in)                  :: c
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: lon_buf(*), lat_buf(*), lev_buf(*), val_buf(*)
    integer(c_int), intent(in)               :: obs_num, n
    integer(c_int)                           :: errno

    errno = redis_setobserve_d(c, c_str(hkey), c_str(varlist), c_loc(lon_buf), &
      c_loc(lat_buf), c_loc(lev_buf), c_loc(val_buf), obs_num, n)
    if (errno /= 0) then
      call err_print
    end if
  end subroutine

  subroutine RedisObscount(c, hkey, its, ite, jts, jte, buf)
    type(c_ptr), intent(in)                  :: c
    character(len=*), intent(in)             :: hkey
    real(c_float), intent(in)               :: its, ite, jts, jte
    integer(c_int)                           :: buf
    integer(c_int)                           :: errno
    
    errno = redis_obscount(c, c_str(hkey), its, ite, jts, jte, buf)
    if (errno /= 0) then
      call err_print
    end if
  end subroutine

  subroutine RedisGetobserveF(c, hkey, varlist, its, ite, jts, jte, &
    lon_buf, lat_buf, lev_buf, val_buf, n)
    type(c_ptr), intent(in)                  :: c
    character(len=*), intent(in)             :: hkey, varlist
    real(c_float), intent(in)               :: its, ite, jts, jte
    integer(c_int), intent(in)               :: n
    type(*), target :: lon_buf(*), lat_buf(*), lev_buf(*), val_buf(*)
    integer(c_int)                           :: errno

    errno = redis_getobserve_f(c, c_str(hkey), c_str(varlist), its, ite, jts, jte, &
      c_loc(lon_buf), c_loc(lat_buf), c_loc(lev_buf), c_loc(val_buf), n)
    if (errno /= 0) then
      call err_print
    end if
  end subroutine

  subroutine RedisGetobserveD(c, hkey, varlist, its, ite, jts, jte, &
    lon_buf, lat_buf, lev_buf, val_buf, n)
    type(c_ptr), intent(in)                  :: c
    character(len=*), intent(in)             :: hkey, varlist
    real(c_double), intent(in)               :: its, ite, jts, jte
    integer(c_int), intent(in)               :: n
    type(*), target :: lon_buf(*), lat_buf(*), lev_buf(*), val_buf(*)
    integer(c_int)                           :: errno

    errno = redis_getobserve_d(c, c_str(hkey), c_str(varlist), its, ite, jts, jte, &
      c_loc(lon_buf), c_loc(lat_buf), c_loc(lev_buf), c_loc(val_buf), n)
    if (errno /= 0) then
      call err_print
    end if
  end subroutine

  !------------------------------------------------------------------------------------
  !  custom part
  !------------------------------------------------------------------------------------

  subroutine RedisHmsetf3d(c, hkey, varlist, its, ite, istride, jts, jte, &
          jstride, kms, kme, kstride, n, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    integer(c_int), intent(in) :: its, ite, jts, jte, kms, kme, n
    integer(c_int), intent(in) :: istride, jstride, kstride
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: buf(*)

    errno = redis_hmsetf_f(c, c_str(hkey), c_str(varlist), its, ite, istride, jts, jte,&
        jstride, kms, kme, kstride, n, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine RedisHmsetf2d(c, hkey, varlist, its, ite, istride, jts, jte, jstride, n, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    integer(c_int), intent(in) :: its, ite, jts, jte, n
    integer(c_int), intent(in) :: istride, jstride
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: buf(*)

    errno = redis_hmsetf2d_f(c, c_str(hkey), c_str(varlist), its, ite, istride, &
        jts, jte, jstride, n, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine RedisHmsetf1d(c, hkey, varlist, its, ite, istride, n, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    integer(c_int), intent(in) :: its, ite, n
    integer(c_int), intent(in) :: istride
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: buf(*)

    errno = redis_hmsetf1d_f(c, c_str(hkey), c_str(varlist), its, ite, istride, n, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine RedisHmgetf3d(c, hkey, varlist, its, ite, istride, jts, jte, jstride, &
          kms, kme, kstride, n, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    integer(c_int), intent(in) :: its, ite, jts, jte, kms, kme, n
    integer(c_int), intent(in) :: istride, jstride, kstride
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: buf(*)

    errno = redis_hmgetf_f(c, c_str(hkey), c_str(varlist), its, ite, istride, jts, jte, &
        jstride, kms, kme, kstride, n, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine RedisHmgetf2d(c, hkey, varlist, its, ite, istride, jts, jte, jstride, n, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    integer(c_int), intent(in) :: its, ite, jts, jte, n
    integer(c_int), intent(in) :: istride, jstride
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: buf(*)

    errno = redis_hmgetf2d_f(c, c_str(hkey), c_str(varlist), its, ite, istride, jts, jte, &
        jstride, n, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine RedisHmgetf1d(c, hkey, varlist, its, ite, istride, n, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    integer(c_int), intent(in) :: its, ite, n
    integer(c_int), intent(in) :: istride
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: buf(*)

    errno = redis_hmgetf1d_f(c, c_str(hkey), c_str(varlist), its, ite, istride, n, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine RedisHmsetd1d(c, hkey, varlist, its, ite, istride, n, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    integer(c_int), intent(in) :: its, ite, n
    integer(c_int), intent(in) :: istride
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: buf(*)

    errno = redis_hmsetd1d_f(c, c_str(hkey), c_str(varlist), its, ite, istride, n, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine RedisHmsetd2d(c, hkey, varlist, its, ite, istride, jts, jte, jstride, n, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    integer(c_int), intent(in) :: its, ite, jts, jte, n
    integer(c_int), intent(in) :: istride, jstride
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: buf(*)

    errno = redis_hmsetd2d_f(c, c_str(hkey), c_str(varlist), its, ite, istride, jts, jte,&
        jstride, n, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine RedisHmsetd3d(c, hkey, varlist, its, ite, istride, jts, jte, jstride,&
          kms, kme, kstride, n, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    integer(c_int), intent(in) :: its, ite, jts, jte, kms, kme, n
    integer(c_int), intent(in) :: istride, jstride, kstride
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: buf(*)

    errno = redis_hmsetd_f(c, c_str(hkey), c_str(varlist), its, ite, istride, jts, jte, &
        jstride, kms, kme, kstride, n, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine RedisHmgetd1d(c, hkey, varlist, its, ite, istride, n, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    integer(c_int), intent(in) :: its, ite, n
    integer(c_int), intent(in) :: istride
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: buf(*)

    errno = redis_hmgetd1d_f(c, c_str(hkey), c_str(varlist), its, ite, istride, n, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine RedisHmgetd2d(c, hkey, varlist, its, ite, istride, jts, jte, jstride, n, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    integer(c_int), intent(in) :: its, ite, jts, jte, n
    integer(c_int), intent(in) :: istride, jstride
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: buf(*)

    errno = redis_hmgetd2d_f(c, c_str(hkey), c_str(varlist), its, ite, istride, jts, jte,&
        jstride, n, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine RedisHmgetd3d(c, hkey, varlist, its, ite, istride, jts, jte, jstride, &
          kms, kme, kstride, n, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    integer(c_int), intent(in) :: its, ite, jts, jte, kms, kme, n
    integer(c_int), intent(in) :: istride, jstride, kstride
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: buf(*)

    errno = redis_hmgetd_f(c, c_str(hkey), c_str(varlist), its, ite, istride, &
        jts, jte, jstride, kms, kme, kstride, n, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  !-----------------------------------------------------------------------------
  ! nonblock interface
  !----------------------------------------------------------------------------- 
  subroutine RedisHmsetf3dNonblock(c, hkey, varlist, its, ite, istride, jts, jte, &
          jstride, kms, kme, kstride, n, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    integer(c_int), intent(in) :: its, ite, jts, jte, kms, kme, n
    integer(c_int), intent(in) :: istride, jstride, kstride
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: buf(*)

    errno = redis_hmsetf_f_nonblock(c, c_str(hkey), c_str(varlist), its, ite, istride, jts, jte,&
        jstride, kms, kme, kstride, n, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine RedisHmsetf2dNonblock(c, hkey, varlist, its, ite, istride, jts, jte, jstride, n, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    integer(c_int), intent(in) :: its, ite, jts, jte, n
    integer(c_int), intent(in) :: istride, jstride
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: buf(*)

    errno = redis_hmsetf2d_f_nonblock(c, c_str(hkey), c_str(varlist), its, ite, istride, &
        jts, jte, jstride, n, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine RedisHmgetf3dNonblock(c, hkey, varlist, its, ite, istride, jts, jte, jstride, &
          kms, kme, kstride, n, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    integer(c_int), intent(in) :: its, ite, jts, jte, kms, kme, n
    integer(c_int), intent(in) :: istride, jstride, kstride
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: buf(*)

    errno = redis_hmgetf_f_nonblock(c, c_str(hkey), c_str(varlist), its, ite, istride, jts, jte, &
        jstride, kms, kme, kstride, n, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine RedisHmgetf2d_nonblock(c, hkey, varlist, its, ite, istride, jts, jte, jstride, n, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    integer(c_int), intent(in) :: its, ite, jts, jte, n
    integer(c_int), intent(in) :: istride, jstride
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: buf(*)

    errno = redis_hmgetf2d_f_nonblock(c, c_str(hkey), c_str(varlist), its, ite, istride, jts, jte, &
        jstride, n, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine


  subroutine RedisHmsetd2dNonblock(c, hkey, varlist, its, ite, istride, jts, jte, jstride, n, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    integer(c_int), intent(in) :: its, ite, jts, jte, n
    integer(c_int), intent(in) :: istride, jstride
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: buf(*)

    errno = redis_hmsetd2d_f_nonblock(c, c_str(hkey), c_str(varlist), its, ite, istride, jts, jte,&
        jstride, n, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine RedisHmsetd3dNonblock(c, hkey, varlist, its, ite, istride, jts, jte, jstride,&
          kms, kme, kstride, n, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    integer(c_int), intent(in) :: its, ite, jts, jte, kms, kme, n
    integer(c_int), intent(in) :: istride, jstride, kstride
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: buf(*)

    errno = redis_hmsetd_f_nonblock(c, c_str(hkey), c_str(varlist), its, ite, istride, jts, jte, &
        jstride, kms, kme, kstride, n, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine
  subroutine RedisHmgetd2dNonblock(c, hkey, varlist, its, ite, istride, jts, jte, jstride, n, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    integer(c_int), intent(in) :: its, ite, jts, jte, n
    integer(c_int), intent(in) :: istride, jstride
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: buf(*)

    errno = redis_hmgetd2d_f_nonblock(c, c_str(hkey), c_str(varlist), its, ite, istride, jts, jte,&
        jstride, n, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine RedisHmgetd3dNonblock(c, hkey, varlist, its, ite, istride, jts, jte, jstride, &
          kms, kme, kstride, n, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    integer(c_int), intent(in) :: its, ite, jts, jte, kms, kme, n
    integer(c_int), intent(in) :: istride, jstride, kstride
    character(len=*), intent(in)             :: hkey, varlist
    type(*), target :: buf(*)

    errno = redis_hmgetd_f_nonblock(c, c_str(hkey), c_str(varlist), its, ite, istride, &
        jts, jte, jstride, kms, kme, kstride, n, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  !------------------------------------------------------------------------------
  ! da_letkf_customize
  !-----------------------------------------------------------------------------

  subroutine redis_da_outputd(c, hkey, varlist,  lonids, lonide, lonstep, &
          latids, latide, latstep, mpas_num_lev, num_lev, zstep, num_2d, num_3d, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    type(*), target :: buf(*)
    character(len=*), intent(in)             :: hkey, varlist
    integer(c_int), intent(in) :: lonids, lonide, lonstep, latids, latide, &
        latstep, mpas_num_lev, num_lev, zstep, num_2d, num_3d
    errno = redis_da_outputd_f(c, c_str(hkey), c_str(varlist), lonids, lonide, lonstep, &
          latids, latide, latstep, mpas_num_lev, num_lev, zstep, num_2d, num_3d, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine redis_da_outputf(c, hkey, varlist,  lonids, lonide, lonstep, &
          latids, latide, latstep, mpas_num_lev, num_lev, zstep, num_2d, num_3d, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    type(*), target :: buf(*)
    character(len=*), intent(in)             :: hkey, varlist
    integer(c_int), intent(in) :: lonids, lonide, lonstep, latids, latide, &
        latstep, mpas_num_lev, num_lev, zstep, num_2d, num_3d
    errno = redis_da_outputf_f(c, c_str(hkey), c_str(varlist), lonids, lonide, lonstep, &
          latids, latide, latstep, mpas_num_lev, num_lev, zstep, num_2d, num_3d, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine
  
  subroutine redis_da_outputd_nonblock(c, hkey, varlist,  lonids, lonide, lonstep, &
          latids, latide, latstep, mpas_num_lev, num_lev, zstep, num_2d, num_3d, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    type(*), target :: buf(*)
    character(len=*), intent(in)             :: hkey, varlist
    integer(c_int), intent(in) :: lonids, lonide, lonstep, latids, latide, &
        latstep, mpas_num_lev, num_lev, zstep, num_2d, num_3d
    errno = redis_da_outputd_f_nonblock(c, c_str(hkey), c_str(varlist), lonids, lonide, lonstep, &
          latids, latide, latstep, mpas_num_lev, num_lev, zstep, num_2d, num_3d, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine redis_da_outputf_nonblock(c, hkey, varlist,  lonids, lonide, lonstep, &
          latids, latide, latstep, mpas_num_lev, num_lev, zstep, num_2d, num_3d, buf)
    type(c_ptr), intent(in)                  :: c
    integer(c_int)                           :: errno
    type(*), target :: buf(*)
    character(len=*), intent(in)             :: hkey, varlist
    integer(c_int), intent(in) :: lonids, lonide, lonstep, latids, latide, &
        latstep, mpas_num_lev, num_lev, zstep, num_2d, num_3d
    errno = redis_da_outputf_f_nonblock(c, c_str(hkey), c_str(varlist), lonids, lonide, lonstep, &
          latids, latide, latstep, mpas_num_lev, num_lev, zstep, num_2d, num_3d, c_loc(buf))
    if (errno /= 0) then
        call err_print
    end if
  end subroutine

  subroutine RedisSetSocketTimeout(cc, sec) 
    type(c_ptr), intent(in)                  :: cc
    integer(c_int)                           :: errno
    integer(c_int), intent(in)               :: sec

    errno = set_redis_socket_timeout_f(cc, sec)
    if (errno /= 0) then
        call err_print
    end if
  end subroutine
    
! 此处注释是因为该函数依赖
! string，而这里没有相关依赖，在模块中使用时，可以取消注释 
  !function RedisConnect_Balanced(group_size, group_rank, redis_address)
  !  use string
  !  implicit none
  !  type(c_ptr) :: RedisConnect_Balanced
  !  integer, allocatable, dimension(:) :: addnum
  !  type(string_type), allocatable :: redis_add(:)
  !  integer,intent(in) :: group_rank, group_size
  !  character(len=*), intent(in)  :: redis_address

  !  redis_add = split_string(redis_address, ',')
  !  allocate(addnum(0:group_size - 1))
  !  addnum(group_rank) = mod(group_rank, size(redis_add)) + 1

  !  RedisConnect_Balanced = RedisConnect(c_str(trim(redis_add(addnum(group_rank))%value)))

  !end function RedisConnect_Balanced
  !------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------
  !  function part
  !------------------------------------------------------------------------------------
  function c_str(f_str)

    character(*), intent(in) :: f_str
    character(len=len_trim(f_str)+1,kind=c_char) c_str

    c_str = trim(f_str) // c_null_char

  end function c_str

  function static_s(f_str, f_len)
    character, allocatable, intent(in) :: f_str(:)
    integer, intent(in) :: f_len
    character(len=f_len) static_s
    integer :: i
    character :: ts

    static_s = ""
    do i = 1, f_len
      ts = f_str(i)
      static_s(i:i) = ts
    end do
  end function static_s

end module redis_mod
