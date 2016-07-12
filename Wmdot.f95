!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the subroutine is coded for modeling molar mass
!    and injected mass through injector when
!    properties calculation of dual mode scramjet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Coded by       : Black Ghost
! Created Date   : 2016-01-07
! Completed Date : 2016-01-15
! Revised Date   :
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ref. William H. Heiser, David T. Pratt. Hypersonic
!   Airbreathing Propulsion AIAA 1994.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note :
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Wmdot(i)
   use VariableDef
   implicit none
   integer::i,k
   real*8::xm,mdot,W,T0,p0
   real*8::dmdot,mfdotg,mfdotl,Wnew,dW
   
   xm=PR(i,1)
   mdot=PR(i-1,15)
   W=PR(i-1,16)
   T0=PR(i-1,6)
   p0=PR(i-1,2)

   do k=1,Nm
      if ((xm > (MASS(k,1)+0.1*h)) .and. (xm <= (MASS(k,1)+Dinj+0.1*h))) then
         mfdotg=MASS(k,2)*mdot*fst                  ! gas fuel mass flow rate
         mfdotl=0                                   ! liquid fuel mass flow rate
         Wnew=(mdot+mfdotg)/(mdot/W+mfdot/2)        ! total added molecular weight
         dW=(Wnew-W)/Dinj
         W=W+dW*h                                   ! current molecular weight
         dWW=dW/W
         mfdot=mfdotg+mfdotl                        ! dmdot of fuel
         dmdot=mfdot/Dinj
         mdot=mdot+dmdot*h
         dmfdotgm=mfdotg/(Dinj*mdot)
         dmfdotlm=mfdotl/(Dinj*mdot)
         dmdotm=dmdot/mdot
         y=MASS(k,3)
         Xa=mdot*1000/W
         Xg=(mfdotg*1000/2)/Xa
         Xl=(mfdotl*1000/2)/Xa
         Xa=(PR(i-1,15)*1000/PR(i-1,16))/Xa
         exit                                       ! exit for keep from excute contents in else statement
      else
         dWW=0.0
         dmfdotgm=0.0
         dmfdotlm=0.0
         dmdotm=0.0
         y=0.0
         Xg=0.0
         Xl=0.0
         Xa=1.0
      end if
   end do
   Vl=y*PR(i-1,5)                                   !!!!!!!!!!!!! evaporated liquid fuel velocity
   PR(i,15)=mdot
   PR(i,16)=W

end subroutine Wmdot
