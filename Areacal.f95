!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the subroutine is coded for modeling area for
!    properties calculation of dual mode scramjet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Coded by       : Black Ghost
! Created Date   : 2016-01-07
! Completed Date : 2016-01-13
! Revised Date   :
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note :
!   1. all angles should be in rad
!   2. assume only contant divergent section occurs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Areacal(i)
   use VariableDef
   implicit none
   integer::i,k
   real*8::xA
   
   xA=PR(i,1)                              ! current axial location
   Hd=sqrt(4*PR(i-1,13)/pi)                ! effective(inviscid) height of duct,
                                           !    assume it is circular section
   if (i == 2) then
      dA=tan(AREA(1,2))
   end if
   if (NA >= 2) then
      do k=2,NA
         if (abs(xA-AREA(k,1)) < 0.5*h) then
            dA=tan(AREA(k,2))
            exit
         end if
      end do
   end if
   
   PR(i,13)=PR(i-1,13)+dA*h
   dAA=dA/PR(i,13)
   Hd=sqrt(4*PR(i-1,13)/pi)
   
end subroutine Areacal
