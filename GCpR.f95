!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the subroutine is coded for modeling specific
!    heat ratio for properties calculation of dual
!    mode scramjet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Coded by       : Black Ghost
! Created Date   : 2016-01-07
! Completed Date : 2016-01-13
! Revised Date   :
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ref. 徐旭, 陈兵, 徐大军. 冲压发动机原理及技术. 北航, 2014.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note :
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GCpR(i)
   use VariableDef
   implicit none
   integer::i
   real*8::b0,b1,b2,b3,b4,b5
   real*8::T,G1
   
   b0=1.817160e3
   b1=9.890504e2
   b2=-9.595592e-3
   b3=1.041469e-4
   b4=-4.433065e-8
   b5=5.879263e-12
   
   G1=gamma
   T=PR(i-1,6)
   
   ! Cp, gamma and R calculation
   
   ! specific heat ratio，ref. 2
   Cp=b1+2*b2*T+3*b3*T**2+4*b4*T**3+5*b5*T**4
   ! specific heat (pressure)，ref. 2
   gamma=(b1+2*b2*T+3*b3*T**2+4*b4*T**3+5*b5*T**4)/(b1-Rg+2*b2*T+3*b3*T**2+4*b4*T**3+5*b5*T**4)
   ! gas constant
   Rg=Cp/(gamma/(gamma-1))
   
!   dGG and gamma storing
   dGG=(gamma-G1)/(PR(i,1)-PR(i-1,1))
   dGG=dGG/gamma
   PR(i,14)=gamma
!   added if there a constant gammma required
   dGG=0
   PR(i,14)=G1
   
end subroutine GCpR
