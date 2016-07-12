!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the subroutine is coded for calculate aerodynamic 
!    parameters across the pre-combustion shock train
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Coded by       : Black Ghost
! Created Date   : 2016-01-25
! Completed Date : 2016-01-27
! Revised Date   :
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ref. William H. Heiser, David T. Pratt. Hypersonic
!    Airbreathing Propulsion. Washington D. C.: AIAA, 1994.
!    pp. 342-346.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note :
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine PCST(i)
   use VariableDef
   implicit none
   integer::i
   real*8::x2,p2,Ma2,A2,gamma2,gamma3a,Ma3a,p3a,Ac3a
   real*8::b,c,l,omega2,omega3,omegax,x
   
   x2=PR(i-1,1)
   p2=PR(i-1,2)
   Ma2=PR(i-1,3)
   A2=PR(i-1,13)
   gamma2=PR(i-1,14)
   
   open(unit=102,file='point3.inp',status='old')
   read(102,*) gamma3a                           ! specific heat ratio at point 3a
   read(102,*) Ma3a                              ! Mach number at point 3
   read(102,*) p3a                               ! static pressure at point 3
   close(102)
   
!   assume :
!   1. no mass addition or extraction within PCST section, constant Rg
!   2. adiabatic wall, constant Tt
!   3. there is no seperation at point 2
   if ((Ma3a == 0.0) .and. (p3a /= 0.0)) then
      Ma3a=gamma2*gamma2*Ma2*Ma2*(1+(gamma2-1)*Ma2*Ma2/2)/((1+gamma2*Ma2*Ma2-p3a/p2)**2)
      Ma3a=(gamma3a*Ma3a/gamma2-(gamma3a-1)/2)**(-0.5)
   else if ((p3a == 0.0) .and. (Ma3a /= 0.0)) then
      p3a=sqrt(gamma3a*(1+(gamma2-1)*Ma2*Ma2/2)/(gamma2*(1+(gamma3a-1)*Ma3a*Ma3a/2)))
      p3a=p2*(1+gamma2*Ma2*Ma2-gamma2*Ma2*Ma3a*p3a)
   else
      write(*,*) 'ERROR : input error at point 3 !'
      read(*,*)
      stop
   end if
   Ac3a=A2*((1+gamma2*Ma2*Ma2)/(p3a/p2)-1)/(gamma3a*Ma3a*Ma3a)
   if ((Ac3a/A2) > 1.0) then
      write(*,*) 'The Mach number at point 3 lower than that through a normal shock from Ma2'
      stop
   end if
!   assume only one wall divergent angle within PCST
   
   c=0.114
   omega2=sqrt(((gamma2-1)*Ma2*Ma2/2)/(1+(gamma2-1)*Ma2*Ma2/2))
   omega3=sqrt(((gamma3a-1)*Ma3a*Ma3a/2)/(1+(gamma3a-1)*Ma3a*Ma3a/2))
   b=c*acosh(omega2/omega3)/(log(omega2/omega3))
   l=(log(omega2/omega3))/c
   
!   jump to point 3a
   PR(i,1)=PR(i-1,1)+l
   PR(i,2)=p3a
   PR(i,3)=Ma3a
   PR(i,6)=PR(i-1,6)*(1+(gamma2-1)*Ma2*Ma2/2)/(1+(gamma3a-1)*Ma3a*Ma3a/2)
   PR(i,4)=p3a/(Rg*PR(i,6))
   PR(i,11)=sqrt(gamma3a*Rg*PR(i,6))
   PR(i,5)=Ma3a*PR(i,11)
   PR(i,13)=A2+dA*l
   PR(i,7)=p3a*PR(i,13)+p3a*Ac3a*gamma3a*Ma3a*Ma3a
   PR(i,8)=PR(i-1,8)
   PR(i,9)=p3a*(1+(gamma3a-1)*Ma3a*Ma3a/2)**(gamma3a/(gamma3a-1))
   PR(i,12)=Ac3a
   PR(i,14)=gamma3a
   PR(i,15)=PR(i-1,15)
   PR(i,16)=PR(i-1,16)
   PR(i,17)=2*(Ac3a-A2)/(l*PR(i,12))-gamma3a*PR(i,3)*PR(i,3)*(ffD+dDxd)
   PR(i,17)=PR(i,17)*(1+(gamma3a-1)*Ma3a*Ma3a/2)-(Ma3a*Ma3a-1)*(gamma3a-gamma2)/(l*gamma3a)
   PR(i,17)=PR(i,17)*Ma3a*Ma3a
   call entro(i)
   
!   x=PR(i-1,1)+l/100
!   
!   do while (x <= (x2+l))
!      PR(i,1)=x
!      omegax=omega3*cosh(b*(l+x2-x))
!      PR(i,2)=(omega2-omegax)*(p3a-p2)/(omega2-omega3)+p2
!      PR(i,3)=PR(i-1,3)+(Ma3a-Ma2)/l
!      PR(i,6)=PR(i-1,6)*(1+(gamma-1)*PR(i-1,3)*PR(i-1,3)/2)/(1+(gamma-1)*PR(i,3)*PR(i,3)/2)
!      PR(i,4)=PR(i,2)/(Rg*PR(i,6))
!      PR(i,11)=sqrt(gamma*Rg*PR(i,6))
!      PR(i,5)=PR(i,3)*PR(i,11)
!!      call Areacal(i+1)                          ! dAA
!      PR(i,13)=PR(i-1,13)+dA*(PR(i,1)-PR(i-1,1))
!      PR(i,12)=
!      PR(i,7)=PR(i,2)*PR(i,13)+PR(i,2)*PR(i,12)*gamma*PR(i,3)*PR(i,3)
!      PR(i,8)=PR(i-1,8)
!      PR(i,9)=PR(i,2)*(1+(gamma-1)*PR(i,3)*PR(i,3)/2)**(gamma/(gamma-1))
!      PR(i,10)=
!      PR(i,15)=PR(i-1,15)
!      PR(i,16)=PR(i-1,16)
!      PR(i,17)=2*(Ac3a-A2)/(l*PR(i,12))-gamma3a*PR(i,3)*PR(i,3)*(ffD+dDxd)
!      PR(i,17)=PR(i,17)*(1+(gamma3a-1)*Ma3a*Ma3a/2)-(Ma3a*Ma3a-1)*(gamma3a-gamma2)/(l*gamma3a)
!      PR(i,17)=PR(i,17)*Ma3a*Ma3a
!      call entro(i)
!      
!      call ffDcal(i+1)                           ! calculate 4Cf/Dh
!      call GCpR(i+1)                             ! Cp, gamma, dGG
!      call Wmdot(i+1)                            ! dWW, dmdotm, y, dmfdotgm, dmfdotlm, Vl, sg, sl
!      call dDxdcal(i+1)                          ! dDxd
!      call dQWxHcal(i+1)                         ! dQWxH
!      x=PR(i,1)+l/100
!      i=i+1
!   end do
   
end subroutine PCST
