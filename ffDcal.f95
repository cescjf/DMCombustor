!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the subroutine is coded to model friction
!      term in dual-mode acramjet simularion
!      4*Cf/Dh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Coded by       : Black Ghost
! Created Date   : 2016-01-13
! Completed Date : 2016-01-13
! Revised Date   :
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note :
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ffDcal(i)
   use VariableDef
   implicit none
   integer::i
   
!   4*cf/Dh
   Dh=sqrt(4*PR(i-1,13)/pi)
   ffD=4*cf/Dh

end subroutine ffDcal
