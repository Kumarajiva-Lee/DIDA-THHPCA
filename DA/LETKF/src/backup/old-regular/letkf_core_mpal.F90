#include "../../utils/da_config.inc"

subroutine letkf_core2(nanals,nobsl,hxens,rdiaginv,dep,rloc,trans,overline_w,inflation_factor, pflag,trans_rate)
   
    use mpal
    use flogger
 
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
    yb_type :: min_eigval
    integer :: i,j,nob,nanal,ierr,lwork,eignums
    integer :: debugi, debugj
 
    integer iwork(5*nanals+3)
    yb_type work(2*nanals*nanals+6*nanals+1)

    type(mpal_st) :: rrloc_mp(nobsl), hxens_mp(nanals, nobsl)
    real(8) :: hxens_test(nanals, nobsl)
    type(mpal_st) :: work1_mp(nanals, nanals)
    real(8) :: work1_test(nanals, nanals)
    type(mpal_st) :: eivec_mp(nanals, nanals)
    type(mpal_st) :: pa_mp(nanals, nanals)
    type(mpal_st) :: work2_mp(nanals, nobsl)
    type(mpal_st) :: overline_w_mp(nanals)
    type(mpal_st) :: trans_mp(nanals, nanals)
    type(mpal_st) :: work2_dep
    type(mpal_st) :: hxens_mp_acc

    allocate(work2(nanals,nobsl))
    allocate(eivec(nanals,nanals),pa(nanals,nanals))
    allocate(work1(nanals,nanals),eival(nanals),rrloc(nobsl))

    do nob=1,nobsl
        rrloc(nob) = sqrt(rdiaginv(nob) * rloc(nob))
        rrloc_mp(nob) = rrloc(nob)
        hxens_mp_acc = 0d0
        do nanal = 1, nanals
            hxens_mp(nanal,nob) = rrloc(nob) * hxens(nanal,nob)
            hxens_test(nanal,nob) = hxens_mp(nanal,nob)%rpe_st%val
            !hxens_mp_acc = hxens_mp_acc + hxens_mp(nanal,nob)
        end do
        !hxens_mp(nanals,nob) = hxens_mp_acc * (-1d0)
        !hxens_mp_acc = hxens_mp_acc + hxens_mp(nanals,nob)
        !if (pflag .and. hxens_mp_acc%rpe_st%val .ne. 0) write(*, *) "warning for", nob, hxens_mp_acc
    end do

    ! do nob = 1, nanals
    !     do nanal = 1, nanals
    !         work1_mp(nanal, nob)%rpe_st%sbits = 23
    !     end do
    ! end do
    
    call mpal_gemm('n','t',nanals,nanals,nobsl,hxens_mp,nanals, &
        hxens_mp,nanals,work1_mp)
    
    do nanal=1,nanals
        work2_dep = real(nanals-1, 8) / inflation_factor
        work1_mp(nanal,nanal) = work1_mp(nanal,nanal) + work2_dep
    end do

    if (pflag) write(*, *) "work1_mp2", work1_mp(1, 1)%rpe_st%sbits

    if (pflag) write(*, *) "work1_mp"
    do i = 1, nanals
        do j = 1, nanals
            work1(j, i) = work1_mp(j, i)%rpe_st%val
        end do
    end do

    if (pflag == 1) then
        write (*,*) 'print work1: '
        do debugi = 1, nanals
           do debugj = 1, nanals
              write(*,668) work1(debugj, debugi)
           end do
           write(*, *) " "
        end do
    end if

    if (pflag == 1) then
        write (*,*) 'print work1_test: '
        do debugi = 1, nanals
           do debugj = 1, nanals
              write(*,668) work1_test(debugj, debugi)
           end do
           write(*, *) " "
        end do
    end if

    if(yb_size == kind(1.d0)) then
        call dsyevd('V','L',nanals,work1,nanals,eival,work,size(work),iwork,size(iwork),ierr)
    else
        call ssyevd('V','L',nanals,work1,nanals,eival,work,size(work),iwork,size(iwork),ierr)
    end if

    if (eival(1) < nanals-1) then
        min_eigval = eival(1) - nanals + 1
        do i = 1, nanals
            eival(i) = eival(i) - min_eigval
        end do
    end if

    if (pflag == 1) then
        write (*,*) 'print syevd: '
        do debugi = 1, nanals
        !    do debugj = 1, nanals
        !       write(*,668) work1(debugj, debugi)
        !    end do
           write(*, 668) eival(debugi)
           write(*, *) " "
        end do
    end if

    if (ierr .ne. 0) then
        print *, "dsyev* failed, ierr=", ierr
    end if

    do j=1,nanals
        do i=1,nanals
            eivec_mp(i, j) = work1(i, j)
            work1(i, j) = work1(i, j) / eival(j)
            work1_mp(i, j)%rpe_st%sbits = RPE_DEFAULT_SBITS
            work1_mp(i, j) = work1(i, j)
        end do
    end do

    ! if (pflag == 1) then
    !     write (*,*) 'print work1 hahha: '
    !     do nob = 1, nanals
    !        do nanal = 1, nanals
    !           write(*, 668) work1(nanal, nob)
    !        end do
    !     end do
    !  end if
    !  if (pflag == 1) then
    !     write (*,*) 'print eivec: '
    !     do nob = 1, nanals
    !        do nanal = 1, nanals
    !           write(*, 668) eivec_mp(nanal, nob)%rpe_st%val
    !        end do
    !     end do
    !  end if
    
    call mpal_gemm('n','t',nanals,nanals,nanals,work1_mp,nanals,eivec_mp,&
        nanals,pa_mp)
    
    do nob=1,nobsl
        do nanal = 1, nanals
            hxens_mp(nanal,nob) = hxens_mp(nanal,nob) * rrloc_mp(nob)
        end do
    end do

    ! if (pflag == 1) then
    !     write (*,*) 'print pa: '
    !     do nob = 1, nanals
    !        do nanal = 1, nanals
    !           write(*, 668) pa_mp(nanal, nob)%rpe_st%val
    !        end do
    !     end do
    !  end if

    call mpal_gemm('n','n',nanals,nobsl,nanals,pa_mp,nanals,hxens_mp,&
        nanals,work2_mp)

    do nanal = 1, nanals
        overline_w_mp(nanal) = 0d0
    end do
    do nob=1,nobsl
        do nanal=1,nanals
            ! if (work2_mp(nanal,nob)%rpe_st%val .ne. work2_mp(nanal,nob)%svalue) then
            !     print *, work2_mp(nanal,nob)%rpe_st%val, work2_mp(nanal,nob)%svalue
            !     call log_error("before overline_w!")
            ! end if
            work2_dep = dep(nob)
            overline_w_mp(nanal) = overline_w_mp(nanal) + work2_mp(nanal,nob) * work2_dep
            ! if (overline_w_mp(nanal)%rpe_st%val .ne. overline_w_mp(nanal)%svalue) then
            !     print *, overline_w_mp(nanal)%rpe_st%val, overline_w_mp(nanal)%svalue
            !     print *, work2_dep%rpe_st%sbits
            !     print *, overline_w_mp(nanal)%rpe_st%sbits
            !     call log_error("after overline_w!")
            ! end if
        end do
    end do

    ! if (pflag == 1) then
    !     write (*,*) 'print overline_w: '
    !     do nanal = 1, nanals
    !         write(*, 668) overline_w_mp(nanal)%rpe_st%val
    !     end do
    !     write (*,*) " "
    ! end if

    do j=1,nanals
        rho = sqrt( real(nanals-1, 8) / eival(j) )
        do i=1,nanals
            work1_mp(i,j) = eivec_mp(i,j) * rho
        end do
    end do

    call mpal_gemm('n','t',nanals,nanals,nanals,work1_mp,nanals,eivec_mp,&
        nanals,trans_mp)

    ! if (pflag == 1) then
    !     write (*,*) 'print trans: '
    !     do nob = 1, nanals
    !         do nanal = 1, nanals
    !            write(*, 668) trans_mp(nanal, nob)%rpe_st%val
    !         end do
    !         write (*,*) " "
    !     end do
    ! end if

    do j=1,nanals
        do i=1,nanals
            trans_mp(i,j) = trans_mp(i,j) + overline_w_mp(i)
            trans(i, j) = trans_mp(i, j)%rpe_st%val
            ! if (trans(i, j) .ne. trans_mp(i, j)%svalue) then
            !     print *, trans(i, j), trans_mp(i, j)%svalue
            !     print *, trans_mp(i, j)%svalue
            !     print *, trans_mp(i, j)%rpe_st%sbits
            !     call log_error("warning")
            ! end if
        end do
        overline_w(j) = overline_w_mp(j)%rpe_st%val
    end do

    do j=1,nanals
        do i=1,nanals
           trans(i, j) = trans_rate * trans(i, j)
           if (i == j) trans(i, j) = trans(i, j) + (1 - trans_rate)
        end do
    end do

    deallocate(work2,eivec,pa,work1,rrloc,eival)

    668 FORMAT(F15.6\)
 
 end subroutine