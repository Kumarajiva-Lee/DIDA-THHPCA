!#include "letkf_config.inc"

module splitchar

  implicit none

contains

subroutine StringSplit(InStr,delimiter,StrArray,nsize)
!----------------------------------------------
!---将字符串InStr进行分割,结果放入StrArray中
!---delimiter::分隔符号,例如';,,' 使用;和,分割字符串
!---nsize:分割数目
!----------------------------------------------
character(len = *) , Intent( IN ) :: InStr
character(len = *)  , Intent( IN ) :: delimiter
character(len = LEN(InStr)),dimension(LEN(InStr)),Intent( OUT ) :: StrArray
integer, Intent( OUT ) :: nsize ! Effective Size of StrArray
integer:: i,j ! loop variable
integer:: istart ! split index for Start Position
!write(*,*)LEN(InStr)
nsize=0
istart=1
do i=1,LEN(InStr)
  do j=1,LEN(delimiter)
    if (InStr(i:i) == delimiter(j:j)) then
      if (istart == i) then
        istart=i+1 ! ---可防止分隔符相连的情况
      end if
      if (istart<i) then
        nsize=nsize+1
        StrArray(nsize)=InStr(istart:i-1)
        istart=i+1
      end if
    end if
  end do
end do
! ---匹配最后一个子字符串
if (nsize>0) then
  if (istart<LEN(InStr)) then
    nsize=nsize+1
    StrArray(nsize)=InStr(istart:LEN(InStr))
  end if
end if
! ---如果无可分割的子字符串,则包含整个字符串为数组的第一元素
if ( (nsize<1) .AND. (LEN(TRIM(InStr)) > 0 )) then
  nsize=1
  StrArray(1)=InStr
end if
end subroutine StringSplit
end module splitchar
