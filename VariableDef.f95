!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the module defines the variables used in main and 
!    subroutines of DualMode program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Coded by       : Black Ghost
! Created Date   : 2016-01-07
! Completed Date : 2016-01-19
! Revised Date   : 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Nomenclature :
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module VariableDef
   implicit none
   integer::Nmax,NA,Nm,Nend,ST
   real*8::pi
   real*8::L,cf,Dh,gamma,Cp,Rg,dGG,h,Dinj,fst,dmdotm,dWW,y,dA,dAA,ffD,dDxd,Hd
   real*8::Tw,Cpf,dmfdotgm,dmfdotlm,Xg,Xl,Xa,Vl,entampre,dQWxH,Tf,pf,time_begin,time_end,mfdot
   real*8,allocatable::PR(:,:),AREA(:,:),MASS(:,:)

end module VariableDef
