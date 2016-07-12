!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The subroutine is the integrator for step-advancing of
!    program DualMode, Runge-Kutta 4-4 is chosen as the 
!    integrate method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Coded by       : Black Ghost
! Created Date   : 2016-01-07
! Completed Date : 2016-01-28
! Revised Date   :
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note :
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Integrate()
   use VariableDef
   implicit none
   integer::i,j
   real*8::Ma4,sh,hmax,hmin
   
   hmax=0.01                                   ! maximum step size
   hmin=0.0001                                 ! minimum step size
   
!  justify if this is case with PSCT#######################################################
   j=2
   PR(j,1)=PR(j-1,1)+h
   call ffDcal(j)                           ! calculate 4Cf/Dh
   call GCpR(j)                             ! Cp, gamma, dGG
   call Wmdot(j)                            ! dWW, dmdotm, y, dmfdotgm, dmfdotlm, Vl, sg, sl
   call Areacal(j)                          ! dAA
   call dDxdcal(j)                          ! dDxd
   call dQWxHcal(j)                         ! dQWxH
   PR(1,17)=-gamma*PR(1,3)*PR(1,3)*(ffD+dDxd)
   PR(1,17)=PR(1,17)*(1+(gamma-1)*PR(1,3)*PR(1,3)/2)
   PR(1,17)=PR(1,17)*PR(1,3)*PR(1,3)
   if ((ST == 1)) then
      call PCST(j)
   else
      call RK44(j,PR(j-1,2),PR(j-1,3),PR(j-1,4),PR(j-1,5),PR(j-1,6),PR(j-1,7),&
         &PR(j-1,8),PR(j-1,9),PR(j-1,11),PR(j-1,12)/PR(j-1,13))     
!!      call for variable step Runge-Kutta method
!      call RK54(j,PR(j-1,2),PR(j-1,3),PR(j-1,4),PR(j-1,5),PR(j-1,6),PR(j-1,7),&
!         &PR(j-1,8),PR(j-1,9),PR(j-1,11),PR(j-1,12)/PR(j-1,13),Ma4)
!!      Ma determines step size
!      if (abs(PR(j,3)-Ma4) <= 1e-3) then
!         sh=0.0
!      else
!         sh=0.84*(1e-3*h/abs(PR(j,3)-Ma4))**0.25
!      end if
!      if ((sh < 0.75) .and. (h > 2*hmin)) then
!         h=h/2
!      else if ((sh > 1.25) .and. (h*2 < hmax)) then
!         h=h*2
!      end if
!      determine the subsonic flow
   end if
   call entro(j)                            ! calculate entropy
   
!   general ODE####################################################################
   do i=j+1,Nmax
!      account for the accumulation error
      if ((PR(i-1,1)+h) > (L+PR(1,1)+0.1*h)) then
         Nend=i-1                              ! steps excuted
         exit                                  ! length exceeded, exit
      end if
      PR(i,1)=PR(i-1,1)+h                      ! current axial location
      
      call ffDcal(i)                           ! calculate 4Cf/Dh
      call GCpR(i)                             ! Cp, gamma, dGG
      call Wmdot(i)                            ! dWW, dmdotm, y, dmfdotgm, dmfdotlm, Vl, sg, sl
      call Areacal(i)                          ! dAA
      call dDxdcal(i)                          ! dDxd
      call dQWxHcal(i)                         ! dQWxH

!      for debug
      if (i == 109) then
         read(*,*)
      end if
      
!      only for cases from subsonic to supersonic
      if (((PR(i-1,17)-PR(i-2,17)) > 0) .and. (abs(PR(i-1,17)) <= 2e-5) .and. (PR(i-1,3) > 0.99)) then
         call jump(i)
         cycle
      end if
      
!      call constant step Runge-Kutta integrator
      call RK44(i,PR(i-1,2),PR(i-1,3),PR(i-1,4),PR(i-1,5),PR(i-1,6),PR(i-1,7),&
         &PR(i-1,8),PR(i-1,9),PR(i-1,11),PR(i-1,12)/PR(i-1,13))
      
!!      call for variable step Runge-Kutta method
!      call RK54(i,PR(i-1,2),PR(i-1,3),PR(i-1,4),PR(i-1,5),PR(i-1,6),PR(i-1,7),&
!         &PR(i-1,8),PR(i-1,9),PR(i-1,11),PR(i-1,12)/PR(i-1,13),Ma4)
!!      Ma determines step size
!      if (abs(PR(i,3)-Ma4) <= 1e-3) then
!         sh=0.0
!      else
!         sh=0.84*(1e-3*h/abs(PR(i,3)-Ma4))**0.25
!      end if
!      if ((sh < 0.75) .and. (h > 2*hmin)) then
!         h=h/2
!      else if ((sh > 1.25) .and. (h*2 < hmax)) then
!         h=h*2
!      end if
!      determine the subsonic flow
      
      call entro(i)                            ! calculate entropy
      
!      if ((PR(i,3) < (1+h)) .and. (PR(i,3) > (1-h))) then
!         write(*,*) 'transonic encountered !'
!         Nend=i
!         read(*,*)
!         exit
!      end if
      
   end do

end subroutine Integrate

! transonic jump
subroutine jump(i)
   use VariableDef
   implicit none
   integer::i
   PR(i,1)=PR(i-1,1)+10*h
   PR(i,2)=PR(i-1,2)
   PR(i,3)=2-PR(i-1,3)
   PR(i,4)=PR(i-1,4)
   PR(i,5)=PR(i,3)*PR(i-1,11)
   PR(i,6)=PR(i-1,6)
   PR(i,7)=PR(i-1,7)
   PR(i,8)=PR(i-1,8)
   PR(i,9)=PR(i-1,9)
   PR(i,10)=PR(i-1,10)
   PR(i,11)=PR(i-1,11)
   PR(i,12)=PR(i-1,12)
   PR(i,13)=PR(i-1,13)
   PR(i,14)=PR(i-1,14)
   PR(i,15)=PR(i-1,15)
   PR(i,16)=PR(i-1,16)
   PR(i,17)=-gamma*PR(1,3)*PR(1,3)*(ffD+dDxd)
   PR(i,17)=PR(1,17)*(1+(gamma-1)*PR(1,3)*PR(1,3)/2)
   PR(i,17)=PR(1,17)*PR(1,3)*PR(1,3)
end subroutine jump
