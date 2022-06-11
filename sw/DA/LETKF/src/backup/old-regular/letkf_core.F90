! Letkf computation core

! nanals     : number of ensember members
! nobsl      : number of observations
! hxens      : transpose of Yb matrix
! rdiaginv   : observation error
! dep        : y_o - ayb
! rloc       : localization function
! trans      : W^a
! overline_w : overline w^a

#include "../../utils/da_config.inc"

subroutine letkf_core(nanals,nobsl,hxens,rdiaginv,dep,rloc,trans,&
   overline_w,inflation_factor,pflag,trans_rate)

   !use flogger

   implicit none
   integer, intent(in ) :: nanals
   integer, intent(in ) :: nobsl
   yb_type, dimension(nanals,nobsl ),intent(inout) :: hxens
   real(4), dimension(nobsl        ),intent(in ) :: rdiaginv
   ayb_type,dimension(nobsl        ),intent(in ) :: dep
   yb_type, dimension(nobsl        ),intent(in ) :: rloc
   yb_type,dimension(nanals,nanals),intent(out) :: trans
   yb_type, dimension(nanals)       ,intent(out) :: overline_w
   real(8), intent(in) :: inflation_factor
   integer, intent(in) :: pflag
   real(8), intent(in) :: trans_rate
   yb_type, allocatable, dimension(:,:) :: work1, eivec, pa, work2
   yb_type, allocatable, dimension(:) :: rrloc,eival
   yb_type :: rho
   integer :: i,j,nob,nanal,ierr,lwork,eignums
   integer :: debugi, debugj

   !for dsyevr
   ! integer iwork(10*nanals),isuppz(2*nanals)
   ! yb_type vl,vu,work(70*nanals)
   !for dsyevd
   integer iwork(5*nanals+3)
   yb_type work(2*nanals*nanals+6*nanals+1)

   allocate(work2(nanals,nobsl))
   allocate(eivec(nanals,nanals),pa(nanals,nanals))
   allocate(work1(nanals,nanals),eival(nanals),rrloc(nobsl))

   !yb_size = kind(hxens)

   ! hxens sqrt(Rinv)
   ! rrloc(1:nobsl) = rdiaginv(1:nobsl) * rloc(1:nobsl)
   ! rho = tiny(rrloc)
   ! where (rrloc < rho) rrloc = rho
   ! rrloc = sqrt(rrloc)

   ! work1 = hdxb^T Rinv hdxb + (m-1) I
   do nob=1, nobsl
      ! rho = 0
      ! do nanal=1,nanals
      !    rho = rho + hxens(nanal, nob) * hxens(nanal, nob)
      ! end do
      ! if (rho == 0) rho = 1
      rrloc(nob) = rdiaginv(nob) * rloc(nob)
      !rrloc(nob) = rdiaginv(nob) * rloc(nob) / rho
      rrloc(nob) = sqrt(rrloc(nob))
      hxens(1:nanals,nob) = hxens(1:nanals,nob) * rrloc(nob)
      !if (pflag) write(*, *) hxens(:, nob), rdiaginv(nob), rloc(nob)
   end do
   ! do nanal=1,nanals
   !    hxens(nanal,1:nobsl) = hxens(nanal,1:nobsl) * rrloc(1:nobsl)
   ! end do
   if(yb_size == kind(1.d0)) then
      call dgemm('n','t',nanals,nanals,nobsl,1.d0,hxens,nanals, &
         hxens,nanals,0.d0,work1,nanals)
   else
      call sgemm('n','t',nanals,nanals,nobsl,1.e0,hxens,nanals, &
         hxens,nanals,0.e0,work1,nanals)
   end if
   do nanal=1,nanals
      work1(nanal,nanal) = work1(nanal,nanal)
   end do

   ! eig decomposition
#if (DEBUG == 1)
   call GPTLstart("eig decomposition")
#endif
   if(yb_size == kind(1.d0)) then
      call dsyevd('V','L',nanals,work1,nanals,eival,work,size(work),iwork,size(iwork),ierr)
      ! call dsyevr('V','A','L',nanals,work1,nanals,vl,vu,1,nanals,-1.d0,eignums,eival,eivec, &
      !             nanals,isuppz,work,size(work),iwork,size(iwork),ierr)
   else
      call ssyevd('V','L',nanals,work1,nanals,eival,work,size(work),iwork,size(iwork),ierr)
      ! call ssyevr('V','A','L',nanals,work1,nanals,vl,vu,1,nanals,-1.d0,eignums,eival,eivec, &
      !             nanals,isuppz,work,size(work),iwork,size(iwork),ierr)
   end if

   rho = 1
   !if (eival(nanals) > nanals-1) rho = (nanals - 1) / eival(nanals) 

   do nanal = 1, nanals
      eival(nanal) = eival(nanal) * rho + real(nanals-1, 8) / inflation_factor
      !if (eival(nanal) > 2*real(nanals-1, 8)/inflation_factor) eival(nanal) = 2*real(nanals-1, 8)/inflation_factor
   end do

#if (DEBUG == 1)
   call GPTLstop("eig decomposition")
#endif
   !if (eival(nanals) > 50d0) call log_error("eival!!!")
   if (ierr .ne. 0) then
      ! call log_error("dsyev* failed, ierr=", to_str(ierr))
      print *, "dsyev* failed, ierr=", ierr
   end if

   if (pflag) print *, "eival: ", eival
   if (pflag) print *, "nobsl", nobsl
   !if (pflag) print *, "dep", dep

   ! pa = [ hdxb^T Rinv hdxb + (m-1) I ]inv
   do j=1,nanals
      do i=1,nanals
         ! work1(i,j) = eivec(i,j) / eival(j)
         eivec(i, j) = work1(i, j)
         work1(i,j) = work1(i,j) / eival(j)
      end do
   end do

   if(yb_size == kind(1.d0)) then
      call dgemm('n','t',nanals,nanals,nanals,1.d0,work1,nanals,eivec,&
         nanals,0.d0,pa,nanals)
   else
      call sgemm('n','t',nanals,nanals,nanals,1.e0,work1,nanals,eivec,&
         nanals,0.e0,pa,nanals)
   end if

   ! convert hxens * Rinv^T from hxens * sqrt(Rinv)^T
   do nob=1,nobsl
      hxens(1:nanals,nob) = hxens(1:nanals,nob) * rrloc(nob)
   end do
   ! do nanal=1,nanals
   !    hxens(nanal,1:nobsl) = hxens(nanal,1:nobsl) * rrloc(1:nobsl)
   ! end do

   ! if (pflag == 1) then
   !    write (*,*) 'print hxens: '
   !    do nob = 1, nob
   !       do nanal = 1, nanals
   !           write(*, 668) hxens(nanal, nob)
   !       end do
   !   end do
   ! end if

   ! Pa hdxb_rinv^T
   if(yb_size == kind(1.d0)) then
      call dgemm('n','n',nanals,nobsl,nanals,1.d0,pa,nanals,hxens,&
         nanals,0.d0,work2,nanals)
   else
      call sgemm('n','n',nanals,nobsl,nanals,1.e0,pa,nanals,hxens,&
         nanals,0.e0,work2,nanals)
   end if

   ! if (pflag == 1) then
   !    write (*,*) 'print work2: '
   !    do nob = 1, nobsl
   !       do nanal = 1, nanals
   !          write(*, 668) work2(nanal, nob)
   !       end do
   !    end do
   ! end if
   ! over-write hxens with Pa hdxb_rinv
   ! (pre-multiply with ensemble perts to compute Kalman gain - 
   !  eqns 20-23 in Hunt et al 2007 paper)
   ! hxens = work2
   ! overline_w = Pa hdxb_rinv^T dep
   ! do nanal=1,nanals
   !    overline_w(nanal) = work2(nanal,1) * dep(1)
   !    do nob=2,nobsl
   !       overline_w(nanal) = overline_w(nanal) + work2(nanal,nob) * dep(nob)
   !    end do
   ! end do
   overline_w(:) = 0d0
   do nob=1,nobsl
      do nanal=1,nanals
         overline_w(nanal) = overline_w(nanal) + work2(nanal,nob) * dep(nob)
      end do
   end do

   ! trans = sqrt[(m-1)Pa]
   do j=1,nanals
      rho = sqrt( real(nanals-1, 8) / eival(j) )
      do i=1,nanals
         work1(i,j) = eivec(i,j) * rho
      end do
   end do
   if(yb_size == kind(1.d0)) then
      call dgemm('n','t',nanals,nanals,nanals,1.d0,work1,nanals,eivec,&
         & nanals,0.d0,trans,nanals)
   else
      call sgemm('n','t',nanals,nanals,nanals,1.e0,work1,nanals,eivec,&
         & nanals,0.e0,trans,nanals)
   end if
   
   !do j=1,nanals
   !   do i=1,nanals
   !      trans(i,j) = work1(i,1) * eivec(j,1)
   !      do k=2,nanals
   !         trans(i,j) = trans(i,j) + work1(i,k) * eivec(j,k)
   !      end do
   !   end do
   !end do
   ! T + Pa hdxb_rinv^T dep
   do j=1,nanals
      do i=1,nanals
         trans(i, j) = trans_rate * (trans(i, j) + overline_w(i))
         if (i == j) trans(i, j) = trans(i, j) + (1 - trans_rate)
      end do
   end do
   deallocate(work2,eivec,pa,work1,rrloc,eival)

   if (pflag) print *, "overline_w", overline_w

   668 FORMAT(F10.5\)
   
   return
end subroutine letkf_core