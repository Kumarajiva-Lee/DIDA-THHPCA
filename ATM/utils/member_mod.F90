module member_mod
    use namelist_mod

    implicit none

    private

    public calc_member_id

    integer, public :: ivector

    contains 

  subroutine calc_member_id(member, member_id, vector_group)
    integer, intent(in) :: member
    integer, intent(out):: member_id
    integer, intent(out):: vector_group

    vector_group = (member + member_num - 1) / member_num
    member_id = member - (vector_group-1)*member_num
 
  end subroutine calc_member_id

end module member_mod