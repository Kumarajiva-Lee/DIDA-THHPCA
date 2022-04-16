! deallocate obs struct
subroutine letkf_final()

    integer :: i, j, k

    do i = -north_n, south_n
        do j = -obs_lc(i)%west_n, obs_lc(i)%east_n
            do k = 1, obs_type_num
                call obs_final(obs_lc(i)%obs_data(k, j))
            end do
        end do
    end do

end subroutine