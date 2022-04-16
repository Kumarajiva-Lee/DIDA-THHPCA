! Letkf computation core

! nanals     : number of ensember members
! nobsl      : number of observations
! hxens      : transpose of Yb matrix
! rdiaginv   : observation error
! dep        : y_o - ayb
! rloc       : localization function
! trans      : W^a
! overline_w : overline w^a

#include "../../../utils/da_config.inc"
#define DEBUG 1

subroutine letkf_core_half(nanals,nobsl,hxens,rdiaginv,dep,rloc,trans,overline_w,inflation_factor,pflag)

    implicit none
    integer, intent(in ) :: nanals
    integer, intent(in ) :: nobsl
    yb_type, dimension(nanals,nobsl ),intent(in) :: hxens
    real(4), dimension(nobsl        ),intent(in ) :: rdiaginv
    ayb_type,dimension(nobsl        ),intent(in ) :: dep
    yb_type, dimension(nobsl        ),intent(in ) :: rloc
    yb_type,dimension(nanals,nanals),intent(out) :: trans
    yb_type, dimension(nanals)       ,intent(out) :: overline_w
    real(8), intent(in) :: inflation_factor
    logical, intent(in) :: pflag

    call letkf_core_half_c(nanals,nobsl,hxens,rdiaginv,dep,rloc,trans,overline_w,inflation_factor,pflag)
   
    return

end subroutine
