module var_info_mod
implicit none

contains
  subroutine var_info(var_name_list, var_num, var_name_list_sort, num_2d, num_3d)
    use coupler_config,only:num_lev

    implicit none

    integer, intent(in):: var_num
    character(20), dimension(var_num), intent(in) :: var_name_list
    character(20), dimension(var_num), intent(out) :: var_name_list_sort
    integer, dimension(var_num) :: var_lev_list
    integer, intent(out) :: num_2d, num_3d
  
    integer :: i, j, k
    integer, dimension(var_num) :: var_lev_list_tmp

    num_2d = 0
    num_3d = 0
    do i=1, var_num
        if (trim(adjustl(var_name_list(i)))=='ps') then
            num_2d = num_2d + 1
            var_lev_list(i) = 1
        else
            num_3d = num_3d + 1
            var_lev_list(i) = num_lev
        endif
    end do

    if (num_2d .ne. 0) then
        do i=1, num_2d
            k=0
            do j=1, num_2d+num_3d
                if (var_lev_list(j) .eq. 1) then
                    var_name_list_sort(i) = var_name_list(j)
                    var_lev_list_tmp(i) = var_lev_list(j)
                else
                    if (i .eq. 1) then
                        k=k+1
                        var_name_list_sort(num_2d+k) = var_name_list(j)
                        var_lev_list_tmp(num_2d+k) = var_lev_list(j)
                    end if
                endif
            end do
        end do
    else
        var_name_list_sort = var_name_list
    end if
  end subroutine var_info
end module var_info_mod
