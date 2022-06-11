#include "../src/letkf_config.inc"

#define pi 3.1415926536

module toy_model

    implicit none

    integer :: ens
    integer :: n
    integer :: nx, ny, nz
    integer :: nx_s, ny_s
    integer :: wx, wy
    integer :: halo

    real(8) :: ens_phase = 0.1
    real(8) :: t_phase = 0.2

    private ens, n, nx, ny, nz, halo, nx_s, ny_s, wx, wy
    private ens_phase, t_phase

contains
    subroutine toy_model_init(ens_in, n_in, nx_in, ny_in, nz_in, halo_in, &
        nx_s_in, ny_s_in, wx_in, wy_in, xb, axb)
        integer, intent(in) :: ens_in
        integer, intent(in) :: n_in
        integer, intent(in) :: nx_in, ny_in, nz_in
        integer, intent(in) :: halo_in
        integer, intent(in) :: nx_s_in, ny_s_in
        integer, intent(in) :: wx_in, wy_in
        xb_type, intent(out) :: xb(ens_in, n_in*nz_in, 1-halo_in:nx_in+halo_in, 1-halo_in:ny_in+halo_in)
        axb_type,intent(out) :: axb(n_in*nz_in,1-halo_in:nx_in+halo_in, 1-halo_in:ny_in+halo_in)

        integer :: i, j, l, ens_num
        axb_type :: t_sum

        ens = ens_in
        n = n_in
        nx = nx_in
        ny = ny_in
        nz = nz_in
        halo = halo_in
        nx_s = nx_s_in
        ny_s = ny_s_in
        wx = wx_in
        wy = wy_in

        do i = 1, ny
            do j = 1, nx
                !do k = 1, n
                    do l = 1, nz*n
                        t_sum = 0.0
                        do ens_num = 1, ens
                            xb(ens_num, l, j, i) = 100 + &
                            100*sin(2*pi*((nx_s+j-1)*1.0/wx+(ny_s+i-1)/2.0/wy+ &
                            ens_phase*ens_num/ens))
                            t_sum = t_sum + xb(ens_num, l, j, i)
                        end do
                        axb(l, j, i) = t_sum / ens
                        do ens_num = 1, ens
                            xb(ens_num, l, j, i) = xb(ens_num, l, j, i) - axb(l, j, i)
                        end do
                    end do
                !end do
            end do
        end do

    end subroutine

    function toy_model_get(lat, lon) result(ans)
        real(8), intent(in) :: lat, lon
        real(8) :: ans
        ans = 100 + 100*sin(2*pi*(lat*1.0/360+lon/2.0/180))
    end function

    subroutine toy_model_advance(xb, axb)
        xb_type, intent(inout) :: xb
        axb_type,intent(inout) :: axb

    end subroutine
end module

#undef pi
