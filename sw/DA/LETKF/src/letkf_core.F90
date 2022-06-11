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
#define DEBUG 1

subroutine letkf_core(nanals,nobsl,hxens_in,rdiaginv,dep,rloc,trans,overline_w,inflation_factor,pflag)

   implicit none
   integer, intent(in ) :: nanals
   integer, intent(in ) :: nobsl
   yb_type, dimension(nanals,nobsl ),intent(in) :: hxens_in
   real(4), dimension(nobsl        ),intent(in ) :: rdiaginv
   ayb_type,dimension(nobsl        ),intent(in ) :: dep
   yb_type, dimension(nobsl        ),intent(in ) :: rloc
   yb_type,dimension(nanals,nanals),intent(out) :: trans
   yb_type, dimension(nanals)       ,intent(out) :: overline_w
   real(8), intent(in) :: inflation_factor
   logical, intent(in) :: pflag
   !real(8), intent(in) :: trans_rate
   yb_type, allocatable, dimension(:,:) :: hxens, work1, eivec, pa, work2
   yb_type, allocatable, dimension(:) :: rrloc,eival
   yb_type :: rho
   integer :: i,j,nob,nanal,ierr,lwork,eignums

   ! for acc tiling
   integer :: acc_tile
   integer :: get_tile

   !for dsyevr
   ! integer iwork(10*nanals),isuppz(2*nanals)
   ! yb_type vl,vu,work(70*nanals)

   !for dsyevd
   integer iwork(5*nanals+3)
   yb_type work(2*nanals*nanals+6*nanals+1)

   allocate(hxens(nanals,nobsl),work2(nanals,nobsl))
   allocate(eivec(nanals,nanals),pa(nanals,nanals))
   allocate(work1(nanals,nanals),eival(nanals),rrloc(nobsl))

   !if (pflag) open(1, file="data.out", form="unformatted")

   !yb_size = kind(hxens)

   ! hxens sqrt(Rinv)
   ! rrloc(1:nobsl) = rdiaginv(1:nobsl) * rloc(1:nobsl)
   ! rho = tiny(rrloc)
   ! where (rrloc < rho) rrloc = rho
   ! rrloc = sqrt(rrloc)

   ! work1 = hdxb^T Rinv hdxb + (m-1) I
   acc_tile = get_tile(nobsl, nanals*yb_size)
   !$ACC PARALLEL LOOP TILE(acc_tile) COPYIN(hxens_in(*,nob)) COPYOUT(hxens(*,nob)) COPYIN(rdiaginv,rloc,nobsl,nanals,acc_tile) COPYOUT(rrloc)
   do nob=1,nobsl
      ! rho = 0
      ! do nanal = 1, nanals
      !   rho = rho + hxens(nanal, nob) * hxens(nanal, nob)
      ! end do
      ! if (rho > 0) then
      !   rrloc(nob) = 1 / rho * rloc(nob)
      ! else
      !   rrloc(nob) = 0
      ! end if
      rrloc(nob) = rdiaginv(nob) * rloc(nob)
      rrloc(nob) = sqrt(rrloc(nob))! / 15
      !if (pflag) print *, "rrloc hxens", rrloc(nob), hxens(1:nanals,nob)
      hxens(1:nanals,nob) = hxens_in(1:nanals,nob) * rrloc(nob)
   end do
   !$ACC END PARALLEL LOOP
   !$ACC SCHEDYIELD
   ! do nanal=1,nanals
   !    hxens(nanal,1:nobsl) = hxens(nanal,1:nobsl) * rrloc(1:nobsl)
   ! end do

#if (DEBUG == 1)
   call GPTLstart("big gemm")
#endif
   if(yb_size == kind(1.d0)) then
      call dgemm('n','t',nanals,nanals,nobsl,1.d0,hxens,nanals, &
         hxens,nanals,0.d0,work1,nanals)
   else
      call sgemm('n','t',nanals,nanals,nobsl,1.e0,hxens,nanals, &
         hxens,nanals,0.e0,work1,nanals)
   end if
#if (DEBUG == 1)
   call GPTLstop("big gemm")
#endif


   if (pflag) then
        write(*, *) "work1", work1
    end if

   !if (pflag) then
   ! write(*, *) "writing"
   ! do nob = 1, nanals
   !  do nanal = 1, nanals
   !     write(1) work1(nanal, nob)
   !  end do
   ! end do
   !end if

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
#if (DEBUG == 1)
   call GPTLstop("eig decomposition")
#endif
   if (ierr .ne. 0) then
      ! call log_error("dsyev* failed, ierr=", to_str(ierr))
      print *, "dsyev* failed, ierr=", ierr
   end if

   do nanal=1,nanals
      !work1(nanal,nanal) = work1(nanal,nanal) + real(nanals-1, 8) / inflation_factor
      eival(nanal) = eival(nanal) + real(nanals-1, 8) / inflation_factor
   end do

   if (pflag) print *, "nobsl", nobsl
   if (pflag) print *, "eival", eival

   ! pa = [ hdxb^T Rinv hdxb + (m-1) I ]inv
   acc_tile = get_tile(nanals, 2*nanals*yb_size)
   !$ACC PARALLEL LOOP TILE(acc_tile) COPYIN(eival(j),nanals,acc_tile) COPY(work1(*,j)) COPYOUT(eivec(*,j))
   do j=1,nanals
      do i=1,nanals
         ! work1(i,j) = eivec(i,j) / eival(j)
         eivec(i, j) = work1(i, j)
         work1(i,j) = work1(i,j) / eival(j)
      end do
   end do
   !$ACC END PARALLEL LOOP
   !$ACC SCHEDYIELD

#if (DEBUG == 1)
   call GPTLstart("gemm")
#endif
   if(yb_size == kind(1.d0)) then
      call dgemm('n','t',nanals,nanals,nanals,1.d0,work1,nanals,eivec,&
         nanals,0.d0,pa,nanals)
   else
      call sgemm('n','t',nanals,nanals,nanals,1.e0,work1,nanals,eivec,&
         nanals,0.e0,pa,nanals)
   end if
#if (DEBUG == 1)
   call GPTLstop("gemm")
#endif

   ! convert hxens * Rinv^T from hxens * sqrt(Rinv)^T
   acc_tile = get_tile(nobsl, nanals*yb_size)
   !$ACC PARALLEL LOOP TILE(acc_tile) COPY(hxens(*, nob)) COPYIN(rrloc(nob),nobsl,nanals,acc_tile)
   do nob=1,nobsl
      hxens(1:nanals,nob) = hxens(1:nanals,nob) * rrloc(nob)
   end do
   !$ACC END PARALLEL LOOP
   !$ACC SCHEDYIELD
   ! do nanal=1,nanals
   !    hxens(nanal,1:nobsl) = hxens(nanal,1:nobsl) * rrloc(1:nobsl)
   ! end do

   ! Pa hdxb_rinv^T
#if (DEBUG == 1)
   call GPTLstart("medium gemm")
#endif
   if(yb_size == kind(1.d0)) then
      call dgemm('n','n',nanals,nobsl,nanals,1.d0,pa,nanals,hxens,&
         nanals,0.d0,work2,nanals)
   else
      call sgemm('n','n',nanals,nobsl,nanals,1.e0,pa,nanals,hxens,&
         nanals,0.e0,work2,nanals)
   end if
#if (DEBUG == 1)
   call GPTLstop("medium gemm")
#endif

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
   !overline_w(:) = 0d0
   acc_tile = get_tile(nanals,(2*nobsl+nanals)*yb_size)
   !$ACC PARALLEL LOOP COPYIN(dep,nanals,nobsl) COPYOUT(overline_w(nanal)) SWAPIN(work2(dimension order:2,1)) ANNOTATE(ENTIRE(dep))
   do nanal=1,nanals
      overline_w(nanal) = 0
      do nob=1,nobsl
          overline_w(nanal) = overline_w(nanal) + work2(nanal,nob) * dep(nob)
      end do
   end do
   !$ACC END PARALLEL LOOP

   ! trans = sqrt[(m-1)Pa]
   acc_tile = get_tile(nanals,2*nanals*yb_size)
   !$ACC PARALLEL LOOP TILE(acc_tile) COPYIN(eival(j),eivec(*,j),nanals,acc_tile) COPYOUT(work1(*,j)) LOCAL(rho)
   do j=1,nanals
      rho = sqrt( real(nanals-1, 8) / eival(j) )
      do i=1,nanals
         work1(i,j) = eivec(i,j) * rho
      end do
   end do
   !$ACC END PARALLEL LOOP
   !$ACC SCHEDYIELD

#if (DEBUG == 1)
   call GPTLstart("gemm")
#endif
   if(yb_size == kind(1.d0)) then
      call dgemm('n','t',nanals,nanals,nanals,1.d0,work1,nanals,eivec,&
         & nanals,0.d0,trans,nanals)
   else
      call sgemm('n','t',nanals,nanals,nanals,1.e0,work1,nanals,eivec,&
         & nanals,0.e0,trans,nanals)
   end if
#if (DEBUG == 1)
   call GPTLstop("gemm")
#endif
   
   !do j=1,nanals
   !   do i=1,nanals
   !      trans(i,j) = work1(i,1) * eivec(j,1)
   !      do k=2,nanals
   !         trans(i,j) = trans(i,j) + work1(i,k) * eivec(j,k)
   !      end do
   !   end do
   !end do
   ! T + Pa hdxb_rinv^T dep
   acc_tile = get_tile(nanals,2*nanals*yb_size)
   !$ACC PARALLEL LOOP TILE(acc_tile) COPY(trans(*,j)) COPYIN(overline_w,nanals,acc_tile) ANNOTATE(ENTIRE(overline_w))
   do j=1,nanals
      do i=1,nanals
         trans(i,j) = trans(i,j) + overline_w(i)
         !if (i == j) trans(i,j) = trans(i,j) + (1 - trans_rate)
      end do
   end do
   !$ACC END PARALLEL LOOP
   !$ACC SCHEDYIELD
   deallocate(hxens,work2,eivec,pa,work1,rrloc,eival)
   
   return
end subroutine letkf_core

function get_tile(tasks, cpy_length)

  implicit none
  
  integer, intent(in) :: tasks, cpy_length

  integer :: get_tile

  integer :: cpes = 64

  get_tile = min((tasks + cpes - 1) / cpes, 131072/cpy_length)

  return

end function get_tile
