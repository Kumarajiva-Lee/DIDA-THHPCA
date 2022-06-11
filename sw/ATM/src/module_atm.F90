module atm_mod

  use atm_parallel

  implicit none

  real*4, allocatable, dimension(:,:) :: u, v, t, ps, w

  contains

  subroutine init_atm()

    integer :: i, j, k

    allocate(u(yst:yed, xst:xed))
    allocate(v(yst:yed, xst:xed))
    allocate(t(yst:yed, xst:xed))
    allocate(ps(yst:yed, xst:xed))
    allocate(w(yst:yed, xst:xed))

    do i = xst, xed
      do j = yst, yed
          u(j, i) = 10.002
          v(j, i) = 20.002
          t(j, i) = 30.002
          w(j, i) = 40.002
          ps(j, i) = 50.002
      end do
    end do

  end subroutine init_atm

  subroutine atm_run()

    implicit none

  end subroutine atm_run

  subroutine atm_final()

    implicit none

    deallocate(u)
    deallocate(v)
    deallocate(t)
    deallocate(ps)
    deallocate(w)

  end subroutine atm_final
end module atm_mod

