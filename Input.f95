!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the subroutine is used ti read-in initial values
!    and initiation for the program of DualMode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Coded by       : Black Ghost
! Created Date   : 2016-01-07
! Completed Date : 2016-01-13
! Revised Date   : 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note :
!   1. input angle in degree, but should be convert
!      to in rad in all calaulations.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Input()
   use VariableDef
   implicit none
   integer::k
!   real*8::
   
   pi=3.141592653
   
   open(unit=101,file='initC.inp',status='old')
   read(101,*) Nmax                        ! maximum step number
   read(101,*) NA                          ! duct section number
   read(101,*) Nm                          ! injector number
   
!   allocate output array
!   second dimension represents in sequence: 
!      PR : x, p, Ma, rho, V, T, F, Tt, Pt, s, a, Ac, A, gamma, mdot, W, G
!      AREA : Ax, theta(slope of area change)
!      MASS : mx, psi, y, theta(injection angle, nornal - 90, parallel - 0)
   allocate(PR(Nmax,17))
   allocate(AREA(NA,2))
   allocate(MASS(Nm,4))
   
   read(101,*) PR(1,1)                     ! axial location of entry plane
   read(101,*) PR(1,2)                     ! main flow pressure enter the duct
   read(101,*) PR(1,3)                     ! main flow Mach number enter the duct
   read(101,*) PR(1,4)                     ! main flow desity enter the duct
   read(101,*) PR(1,5)                     ! main flow velocity enter the duct
   read(101,*) PR(1,6)                     ! main flow temperature enter the duct
   read(101,*) PR(1,7)                     ! main flow axial force act on the entry of the duct
   read(101,*) PR(1,8)                     ! main flow stagnation temperature enter the duct
   read(101,*) PR(1,9)                     ! main flow stagnation pressure enter the duct
   read(101,*) PR(1,10)                    ! main flow entropy enter the duct
   read(101,*) PR(1,11)                    ! main flow sound speed enter the duct
   read(101,*) PR(1,12)                    ! main flow core flow area enter the duct
   read(101,*) PR(1,13)                    ! geometrical entry area of duct
   read(101,*) PR(1,14)                    ! ratio of specific heat of entry flow
   
   read(101,*) cf                          ! coefficient of friction
   read(101,*) Dinj	                      ! injector diameter
   read(101,*) L                           ! duct length
   
!   input axial changes of duct area, area section - 1
!   x(A), theta
   do k=1,NA
      read(101,*) AREA(k,1),AREA(k,2)
      AREA(k,2)=AREA(k,2)*pi/180           ! all angles convert to in rad
   end do
   
!   input inject flow, flow section - 1
!   x(m), psi, y, angle
   do k=1,Nm
      read(101,*) MASS(k,1),MASS(k,2),MASS(k,3),MASS(k,4)
   end do
   
   close(101)
!   initiation
   gamma=PR(1,14)
   Rg=(PR(1,5)/PR(1,3))**2/(PR(1,14)*PR(1,6))
   Cp=Rg*PR(1,14)/(PR(1,14)-1)
   h=0.001
   fst=6.72/(8*28.92)                             ! stochiometric fuel / air ratio
   PR(1,15)=PR(1,4)*PR(1,13)*PR(1,5)              ! mdot=rho*A*V
   dmdotm=0.0
   PR(1,16)=28*0.78+32*0.21+36*0.01               ! molecule weight of air
   entampre=0.0
!   dWW=0.0
!   dA=0.0                                         ! assume constant area section at duct entry
   Cpf=14300                                      ! Cp of fuel, H2
   Tw=800                                         ! wall temperature
   Tf=500                                         ! fuel temperature
   pf=101325                                      ! fuel pressure
   write(*,*) "Is this with a PCST ? yes for '1', otherwise for '0'."
   read(*,*) ST                                   ! if it is not PCST free
   call cpu_time(time_begin)
   write (*,*) time_begin
   
end subroutine Input
