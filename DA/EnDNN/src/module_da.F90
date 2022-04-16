module da_mod

  use da_parallel
 ! use coupler_config

  implicit none

  real*4, allocatable, dimension(:,:) :: u, v, t, ps, w
  
  contains

  subroutine init_da()

    integer :: i, j, k

    allocate(u(yst:yed, xst:xed))
    allocate(v(yst:yed, xst:xed))
    allocate(t(yst:yed, xst:xed))
    allocate(w(yst:yed, xst:xed))
    allocate(ps(yst:yed, xst:xed))

    do i = xst, xed
      do j = yst, yed
          u(j, i) = 10.003
          v(j, i) = 20.003
          t(j, i) = 30.003
          w(j, i) = 40.003
          ps(j, i) = 50.003
      end do
    end do
    
  end subroutine init_da

  subroutine da_run()

    implicit none

  end subroutine da_run

  subroutine da_final()

    implicit none

    deallocate(u)
    deallocate(v)
    deallocate(t)
    deallocate(ps)
    deallocate(w)

  end subroutine da_final
end module da_mod

