! Test letkf

#include "../src/letkf_config.inc"
#include "model_config.inc"

program main

    use mpi
    use toy_model
    use letkf_mod
    use pro_info
    implicit none 

    integer :: pn, sqrtpn
    integer :: bni, bnj
    integer :: rowi, colj
    integer :: nx, ny
    integer :: nx_s, ny_s
    integer :: pi, pj
    integer :: iter
    real(8) :: maxerror
    ! real(8) :: dlat, dlon
    real(8) :: tlat, tlon, wlat, wlon
    character*10 :: file1= "actual.out"
    character*8  :: file2= "back.out"
    character*7  :: file3= "ana.out"
    character*9  :: file4= "param.out"
    character*1  :: fileid

    character * (MPI_MAX_PROCESSOR_NAME) processor_name
    integer myid, numprocs, namelen, ierr
    
    xb_type, allocatable, dimension(:,:,:,:,:) :: xb
    axb_type,allocatable, dimension(:,:,:,:)   :: axb

    xb_type, allocatable, dimension(:,:,:,:,:) :: xa
    axb_type,allocatable, dimension(:,:,:,:)   :: axa

    axb_type,allocatable, dimension(:,:,:,:)   :: actual

    integer :: i, j, k, nn, ens_num
    integer :: debugi, debugj, debugk
    integer :: marki, markj
    integer :: halo
    real(8) :: debug_sum

    call mpi_init(ierr)
    call mpi_comm_rank( MPI_COMM_WORLD, myid, ierr)
    call mpi_comm_size( MPI_COMM_WORLD, numprocs, ierr)
    call mpi_get_processor_name(processor_name, namelen, ierr)

    pn = numprocs / 8
    sqrtpn = int(sqrt(REAL(pn)))
    if (sqrtpn * sqrtpn == pn) then
        bni = sqrtpn
        bnj = sqrtpn
    else
        bni = sqrtpn * 2
        bnj = sqrtpn 
    end if
    rowi = bni * 4
    colj = bnj * 2
    nx = world_x / rowi
    ny = world_y / colj
    call id2pos(myid, pi, pj)
    nx_s = pi * nx + 1
    ny_s = pj * ny + 1
    halo = 2
    call pro_info_init(nx, ny, 360d0 / world_x, 180d0 / world_y, pi + 1, pj + 1)
    ! dlat = 360d0 / world_x
    ! dlon = 180d0 / world_y

    allocate(xb(ens, n, lz, 1-halo:nx+halo, 1-halo:ny+halo))
    allocate(axb(n, lz, 1-halo:nx+halo, 1-halo:ny+halo))
    allocate(xa(ens, n, lz, 1-halo:nx+halo, 1-halo:ny+halo))
    allocate(axa(n, lz, 1-halo:nx+halo, 1-halo:ny+halo))
    allocate(actual(n, lz, 1-halo:nx+halo, 1-halo:ny+halo))

    write(fileid,"(I1)") myid

    open(100,file='./result/actual/'//trim(file1)//fileid,form='formatted',status='replace')
    open(101,file='./result/back/'//trim(file2)//fileid,form='formatted',status='replace')
    open(102,file='./result/ana/'//trim(file3)//fileid,form='formatted',status='replace')
    open(103,file='./result/param/'//trim(file4)//fileid,form='formatted',status='replace')


    write(103,'(i6,i6,i6,i6,i6,i6,f8.2,f8.2)') ny,nx,lz,n,nx_s,ny_s,dlat,dlon 

    !write(*, *) "hello", myid
    !write(*, *) "pij", pi, pj, "nxy_s", nx_s, ny_s

    call letkf_init(ens, nx, ny, lz, halo, n, &
    10, 10, 1, myid, pi + 1, pj + 1, rowi, colj, &
    4, 0, pos2id)
    call toy_model_init(ens, n, nx, ny, lz, halo, nx_s, ny_s, world_x, world_y, xb, axb)

    do i = 1, ny
        do j = 1, nx
            tlat = (nx_s + j - 2) * dlat
            tlon = (ny_s + i - 2) * dlon
            do k = 1, lz
                do nn = 1, n
                    actual(nn, k, j, i) = toy_model_get(tlat, tlon)
                    wlat = (nx_s+j-1)*dlat
                    wlon = (ny_s+i-1)*dlon
                    write(100,'(f10.1,f10.1,f10.4)') wlon,wlat,actual(nn, k, j, i) 
                end do
            end do
        end do
    end do

    do iter = 1, 20

        call letkf_run(xb, axb, xa, axa)

        maxerror = 0
        do i = 1, ny
            do j = 1, nx
                do k = 1, lz
                    do nn = 1, n
                        if (abs(axb(nn, k, j, i) - actual(nn, k, j, i)) > maxerror) then
                            maxerror = abs(axb(nn, k, j, i) - actual(nn, k, j, i))
                            marki = i
                            markj = j
                        end if
                    end do
                end do
                ! if (myid == 1 .and. j == 90) then
                !     write(*, *) actual(1, 1, j, i)
                ! end if
            end do
        end do

        !if (myid == 0) write(*, *) "ij: ", marki, markj
        !if (myid == 0) write(*, *)"maxerr: ", maxerror, axb(1, 1, markj, marki), actual(1, 1, markj, marki)
        write(*, 300) myid, iter, maxerror
        300 format("myid ", I3, " iter ", I3, " maxerr ", F10.2)

        ! if (myid == 0) then
        !     write (*,*) myid, 'xb'
        !     debug_sum = 0
        !     do debugi = 1, ens
        !         write(*,101) xb(debugi, 1, 1, 1, 1)
        !         debug_sum = debug_sum + xb(debugi, 1, 1, 1, 1)
        !     end do
        !     write(*, *) axb(1, 1, 1, 1), " ", debug_sum/ens
        ! end if
        ! if (myid == 0) then
        !     write (*,*) myid, 'xa'
        !     debug_sum = 0
        !     do debugi = 1, ens
        !         write(*,101) xa(debugi, 1, 1, 1, 1)
        !         debug_sum = debug_sum + xa(debugi, 1, 1, 1, 1)
        !     end do
        !     write(*, *) axa(1, 1, 1, 1), " ", debug_sum/ens
        !     101 FORMAT(F10.2\)
        ! end if

        do i = 1, ny
            do j = 1, nx
                do k = 1, lz
                    do nn = 1, n
                       !write(101,'(i6,f10.4)') myid,axb(nn, k, j, i) 
                       wlat = (nx_s+j-1)*dlat
                       wlon = (ny_s+i-1)*dlon
                       write(101,'(i3,f10.1,f10.1,f10.4)') iter,wlon,wlat,axb(nn, k, j, i) 
                    end do
                end do
            end do
        end do

        xb(:, :, :, :, :) = xa(:, :, :, :, :)
        axb(:, :, :, :) = axa(:, :, :, :)

        do i = 1, ny
            do j = 1, nx
                do k = 1, lz
                    do nn = 1, n
                       !write(102,'(i6,f10.4)') myid,axb(nn, k, j, i) 
                       wlat = (nx_s+j-1)*dlat
                       wlon = (ny_s+i-1)*dlon
                       write(102,'(i3,f10.1,f10.1,f10.4)') iter,wlon,wlat,axb(nn, k, j, i) 
                    end do
                end do
            end do
        end do

        ! do debugi = 1, nx
        !     do debugj = 1, ny
        !         do debugk = 1, ens
        !             if (isnan(xa(debugk, 1, 1, debugj, debugi))) then
        !                 print *, "iter ", iter, " myid ", myid, debugk, debugi, debugj
        !                 stop
        !             end if
        !         end do
        !     end do
        ! end do

    end do 

    deallocate(xb)
    deallocate(axb)
    deallocate(xa)
    deallocate(axa)
    deallocate(actual)

contains
    subroutine id2pos(id, pi, pj)
        integer, intent(in)  :: id
        integer, intent(out) :: pi, pj
        integer bid, boff

        bid = id / pn
        boff = mod(id, pn)
        pi = mod(bid, 4) * bni + mod(boff, bni)
        pj = bid / 4 * bnj + boff / bni

    end subroutine

    subroutine pos2id(pi, pj, id)
        integer, intent(in) :: pi, pj
        integer, intent(out):: id
        integer bid, boff
        integer bis, bjs
        integer bi, bj
        
        bid = pj / bnj * 4 + pi / bni
        bis = mod(bid, 4) * bni
        bjs = bid / 4 * bnj
        boff = (pj - bjs) * bni + pi - bis

        id = bid * pn + boff
    
    end subroutine

end program
