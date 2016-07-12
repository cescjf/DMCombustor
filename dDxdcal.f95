!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the subroutine is coded for modeling drag for
!    properties calculation of dual mode scramjet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Coded by       : Black Ghost
! Created Date   : 2016-01-07
! Completed Date : 2016-01-13
! Revised Date   :
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the drag is the sum of : 
!    1. the drag of stationary bodies immersed in the stream
!       within the control-surface boundaries;
!    2. the drag of liquid droplets anf filaments traveling
!       more sloely than the main gas stream, and,
!    3. the component of body or gravity forces acting on
!       the material within the control surface in the 
!       direction opposite to that of the velocity vector.
! ref. to Ascher H. Shapiro, the dynamics and aerodynamics 
!    of compressible fluid flow, vol. I, pp. 224-225
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note :
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dDxdcal(i)
   use VariableDef
   implicit none
   integer::i
   real*8::p,A,Ma
   real*8::dDx
!   real*8::
   dDx=0.0
   p=PR(i-1,2)
   A=PR(i-1,13)
   Ma=PR(i-1,3)
!   dDx/(gamma*p*A*Ma*Ma/2)
   dDxd=dDx/(gamma*p*A*Ma*Ma/2)
   
end subroutine dDxdcal
