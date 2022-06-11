module atm_config
  implicit none

  !namelist total
  character(30)  :: atm_split_scheme
  integer        :: atm_dt_in_seconds
  integer        :: atm_reduce_factors
  integer        :: atm_damp_orders
  integer        :: atm_fast_cycles
  integer        :: atm_pv_scheme
  integer        :: atm_coarse_polar_lats_exp
  integer        :: atm_coarse_polar_lats_mul

contains
  subroutine init_config_atm()
    use coupler_config
    use redis_mod

    implicit none

    type(c_ptr) :: rds_cc

    rds_cc = RedisConnect(redis_address)

    call RedisHget(rds_cc, case_name, 'atm_split_scheme', atm_split_scheme)
    call RedisHget(rds_cc, case_name, 'atm_dt_in_seconds', atm_dt_in_seconds)
    call RedisHget(rds_cc, case_name, 'atm_reduce_factors', atm_reduce_factors )
    call RedisHget(rds_cc, case_name, 'atm_damp_orders',  atm_damp_orders)
    call RedisHget(rds_cc, case_name, 'atm_fast_cycles', atm_fast_cycles)
    call RedisHget(rds_cc, case_name, 'atm_pv_scheme', atm_pv_scheme)
    call RedisHget(rds_cc, case_name, 'atm_coarse_polar_lats_exp', atm_coarse_polar_lats_exp)
    call RedisHget(rds_cc, case_name, 'atm_coarse_polar_lats_mul', atm_coarse_polar_lats_mul)

  end subroutine
end module atm_config 

