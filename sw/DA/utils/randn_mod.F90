module randn_mod
  Implicit None

  double precision, parameter :: pi = 3.141592653589793239
! ran return a uniform random number between 0-1  
! norma return a normal distribution  
contains 
  function ran()   !returns random number between 0 - 1  
    implicit none 
    integer , save :: flag = 0
    double precision :: ran 
    if(flag==0) then 
      call random_seed()
      flag = 1 
    endif 
    call random_number(ran)     ! built in fortran 90 random number
  end function ran
  
  function randn(n) result(res)
    implicit none 
    integer, intent(in) :: n
    integer :: flag, i
    double precision :: u1, u2
    double precision :: res(n)
    save flag 
    data flag /0/ 
    do i =1, n
        u1 = ran(); u2 = ran()
        if (flag.eq.0) then 
          res(i) = sqrt(-2.0d0*log(u1))*cos(2.0d0*pi*u2) 
          flag = 1 
        else 
          res(i) = sqrt(-2.0d0*log(u1))*sin(2.0d0*pi*u2) 
          flag = 0 
        endif  
    enddo
  end function randn
  
  function randn_fix(n) result(res)
      implicit none
      integer, intent(in) :: n
      real :: res(n)
      real, allocatable :: r(:,:)
      real, allocatable :: zeta(:,:)

      if (mod(n,2) .ne. 0) then 
          allocate(r((n+1)/2,2))
          allocate(zeta((n+1)/2,2))
      else
          allocate(r(n/2,2))
          allocate(zeta(n/2,2))
      end if

      call random_number(zeta)
      r(:,1) = sqrt(-2.0*log(zeta(:,1)))*cos(2.0*pi*zeta(:,2))
      r(:,2) = sqrt(-2.0*log(zeta(:,1)))*sin(2.0*pi*zeta(:,2))

      if (mod(n,2) .ne. 0) then
          res(1:(n+1)/2) = r(:,1)
          res((n+1)/2+1:n) = r(1:(n+1)/2-1,2)
      else
          res(1:n/2) = r(:,1)
          res(n/2:n) = r(:,2)
      end if

      deallocate(r)
      deallocate(zeta)
  end function randn_fix

end module randn_mod

