!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the subroutine is coded for calculate entropy after every
!    step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Coded by       : Black Ghost
! Created Date   : 2016-01-19
! Completed Date : 2016-01-19
! Revised Date   :
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ref. Maurice J. Zucrow and Joe D. Hoffman. Gas Dynamics V1.
!    New York: John Wiley & Sons, Inc., 1976. pp. 51-52
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note :
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine entro(i)
   use VariableDef
   implicit none
   integer::i
   real*8::T0,Ti,p0,pii,ds
   
   T0=PR(i-1,6)
   Ti=PR(i,6)
   p0=PR(i-1,2)
   pii=PR(i,2)
   ds=Cp*Xa*log(Ti/T0)+Cpf*Xg*log(Ti/Tf)+Cpf*Xl*log(Ti/Tf)
   ds=ds-Rg*(Xa*log(pii/p0)+Xg*log(pii/pf)+Xl*log(pii/pf))
!   ln(0) is singular point of ln(x)
   if ((Xg == 0) .and. (Xl == 0)) then
      ds=ds-Rg*Xa*log(Xa)
   else if (Xl == 0) then
      ds=ds-Rg*(Xa*log(Xa)+Xg*log(Xg))
   else
      ds=ds-Rg*(Xa*log(Xa)+Xg*log(Xg)+Xl*log(Xl))
   end if
   PR(i,10)=PR(i-1,10)+ds

end subroutine entro
