#include "../../mp_config/mct_config.inc"

module atm_coupler

  use m_MCTWorld,only: MCTWorld_init => init
  use m_MCTWorld,only: MCTWorld_clean => clean
  use m_GlobalSegMap,only: GlobalSegMap
  use m_GlobalSegMap,only: GlobalSegMap_init => init
  use m_GlobalSegMap,only: GlobalSegMap_lsize => lsize
  use m_GlobalSegMap,only: GlobalSegMap_clean => clean
  use m_GlobalSegMap,only: GlobalSegMap_Ordpnts => OrderedPoints
  use m_AttrVect,only    : AttrVect
  use m_AttrVect,only    : AttrVect_init => init
  use m_AttrVect,only    : AttrVect_clean => clean
  use m_AttrVect,only    : AttrVect_indxR => indexRA
  use m_AttrVect,only    : AttrVect_importRAttr => importRAttr
  use m_Router,only: Router
  use m_Router,only: Router_init => init
  use m_Router,only: Router_clean => clean
  use m_Transfer,only : MCT_Send => send
  use m_Transfer,only : MCT_Recv => recv

  use mpi
  use atm_parallel
  use coupler_config
  use atm_mod

  implicit none

  type(GlobalSegMap)   :: ATM_GSMap             ! MCT defined type
  type(AttrVect) :: ATM_AV                      ! MCT defined type
  type(Router) :: R_ATM2DA                      ! MCT defined type

  integer, allocatable,dimension(:) :: start,length
  axb_type, dimension(:), pointer :: avdata
  integer :: avsize

  contains

subroutine init_atm_coupler(MPI_ATM_GROUP)

  implicit none
  integer, intent(in) :: MPI_ATM_GROUP
  integer :: ierr
  integer :: myrank
  integer i

  call MPI_COMM_RANK(MPI_ATM_GROUP, myrank, ierr)
  call MCTWorld_init(ncomps, MPI_COMM_WORLD, MPI_ATM_GROUP, ATM_Id)

  allocate(start(xst:xed))
  allocate(length(xst:xed))

  do i = xst, xed
    length(i) = yed - yst + 1
    start(i) = (i - 1) * atm_num_lat + yst
  end do

  call GlobalSegMap_init(ATM_GSMap, start, length, 0, MPI_ATM_GROUP, ATM_Id)

  deallocate(start)
  deallocate(length)

  avsize = GlobalSegMap_lsize(ATM_GSMap,MPI_ATM_GROUP)

  print *, "atmid:", xst, xed, yst, yed, ATM_Id, avsize

  call AttrVect_init(ATM_AV,rList="u:v:t:ps:w",lsize=avsize)

  call Router_init(DA_Id,ATM_GSMap,MPI_ATM_GROUP,R_ATM2DA)

end subroutine 

subroutine atm2da_coupling(MPI_ATM_GROUP)

  implicit none
  integer, intent(in) :: MPI_ATM_GROUP
  integer i, j, ij

  allocate(avdata(1:avsize))

  ij=0
  do i = xst, xed - 1
    do j = yst, yed - 1
      ij = ij + 1
      avdata(ij) = u(j, i)
    end do
  end do
  call AttrVect_importRAttr(ATM_AV,"u",avdata)

  ij=0
  do i = xst, xed - 1
    do j = yst, yed - 1
      ij = ij + 1
      avdata(ij) = v(j, i)
    end do
  end do
  call AttrVect_importRAttr(ATM_AV,"v",avdata)

  ij=0
  do i = xst, xed - 1
    do j = yst, yed - 1
      ij = ij + 1
      avdata(ij) = t(j, i)
    end do
  end do
  call AttrVect_importRAttr(ATM_AV,"t",avdata)

  ij=0
  do i = xst, xed - 1
    do j = yst, yed - 1
      ij = ij + 1
      avdata(ij) = w(j, i)
    end do
  end do
  call AttrVect_importRAttr(ATM_AV,"w",avdata)

  ij=0
  do i = xst, xed - 1
    do j = yst, yed - 1
      ij = ij + 1
      avdata(ij) = ps(j, i)
    end do
  end do
  call AttrVect_importRAttr(ATM_AV,"ps",avdata)

  deallocate(avdata)
  call MCT_Send(ATM_AV, R_ATM2DA)

end subroutine 

subroutine clean_atm_coupler()

  implicit none

  call Router_clean(R_ATM2DA)
  call AttrVect_clean(ATM_AV)
  call GlobalSegMap_clean(ATM_GSMap)
  call MCTWorld_clean()

end subroutine 

end module atm_coupler
