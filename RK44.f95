!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The subroutine is the constant step Runge-Kutta method 
!    kernel used in the integrate method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Coded by       : Black Ghost
! Created Date   : 2016-01-12
! Completed Date : 2016-01-27
! Revised Date   :
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ref. 胡健伟, 汤怀民. 微分方程数值方法, 第二版. 北京：科学出版社, 2007
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note :
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine RK44(i,p1,Ma1,rho1,V1,T1,F1,Tt1,Pt1,a1,AcA1)
   use VariableDef
   implicit none
   integer::i,k
   real*8::p1,Ma1,rho1,V1,T1,F1,Tt1,Pt1,a1,AcA1,Ac1
   real*8::kp1,kMa1,krho1,kV1,kT1,kF1,kTt1,kPt1,ka1,kAcA1,kAc1
   real*8::kp2,kMa2,krho2,kV2,kT2,kF2,kTt2,kPt2,ka2,kAcA2,kAc2
   real*8::kp3,kMa3,krho3,kV3,kT3,kF3,kTt3,kPt3,ka3,kAcA3,kAc3
   real*8::kp4,kMa4,krho4,kV4,kT4,kF4,kTt4,kPt4,ka4,kAcA4,kAc4
   real*8::kp,kMa,krho,kV,kT,kTt,kAc
   real*8::funp44,funMa44,funrho44,funV44,funT44,funF44,funTt44,funPt44,funa44,funAcA44,funGx44
   
   Ac1=AcA1*PR(i-1,13)
   if ((i == 2) .or. (ST == 1)) then
      kp=0.0
      kMa=0.0
      krho=0.0
      kV=0.0
      kT=0.0
      kTt=0.0
      kAc=0.0 
      ST=0
   else
      kp=(PR(i-1,2)-PR(i-2,2))/h
      kMa=(PR(i-1,3)-PR(i-2,3))/h
      krho=(PR(i-1,4)-PR(i-2,4))/h
      kV=(PR(i-1,5)-PR(i-2,5))/h
      kT=(PR(i-1,6)-PR(i-2,6))/h
      kTt=(PR(i-1,8)-PR(i-2,8))/h
      kAc=(PR(i-1,12)-PR(i-2,12))/h
   end if
!   do this for stable after injection
   do k=1,Nm
      if (((PR(i,1)-MASS(k,1)-Dinj) > 0.9*h) .and. ((PR(i,1)-MASS(k,1)-Dinj) < 4.1*h)) then
         kp=(PR(10,2)-PR(9,2))/h
         kMa=(PR(10,3)-PR(9,3))/h
         krho=(PR(10,4)-PR(9,4))/h
         kV=(PR(10,5)-PR(9,5))/h
         kT=(PR(10,6)-PR(9,6))/h
         kTt=(PR(10,8)-PR(9,8))/h
         kAc=(PR(10,12)-PR(9,12))/h 
      end if
   end do
   
   kp1=funp44(p1,rho1,T1,krho,kT)
   kMa1=funMa44(Ma1,V1,T1,kV,kT)
   krho1=funrho44(rho1,Ac1,V1,kAc,kV)
   kV1=funV44(V1,Ma1,p1,kp)
   kT1=funT44(T1,Ma1,V1,kV)
   kF1=funF44(F1,Ma1,p1,kMa,kp)
   kTt1=funTt44(Tt1,Ma1,T1,kMa,kT)
   kPt1=funPt44(Pt1,Ma1,p1,kMa,kp)
   ka1=funa44(a1,T1,kT)
   kAcA1=funAcA44(AcA1,Ma1,p1,Tt1,Ac1,kp,kTt,kAc1)
   
   kp2=funp44(p1+h*kp1/2,rho1,T1,krho,kT)
   kMa2=funMa44(Ma1+h*kMa1/2,V1,T1,kV,kT)
   krho2=funrho44(rho1+h*krho1/2,Ac1,V1,kAc,kV)
   kV2=funV44(V1+h*kV1/2,Ma1,p1,kp)
   kT2=funT44(T1+h*kT1/2,Ma1,V1,kV)
   kF2=funF44(F1+h*kF1/2,Ma1,p1,kMa,kp)
   kTt2=funTt44(Tt1+h*kTt1/2,Ma1,T1,kMa,kT)
   kPt2=funPt44(Pt1+h*kPt1/2,Ma1,p1,kMa,kp)
   ka2=funa44(a1+h*ka1/2,T1,kT)
   kAcA2=funAcA44(AcA1+h*kAcA1/2,Ma1,p1,Tt1,Ac1,kp,kTt,kAc2)
   
   kp3=funp44(p1+h*kp2/2,rho1,T1,krho,kT)
   kMa3=funMa44(Ma1+h*kMa2/2,V1,T1,kV,kT)
   krho3=funrho44(rho1+h*krho2/2,Ac1,V1,kAc,kV)
   kV3=funV44(V1+h*kV2/2,Ma1,p1,kp)
   kT3=funT44(T1+h*kT2/2,Ma1,V1,kV)
   kF3=funF44(F1+h*kF2/2,Ma1,p1,kMa,kp)
   kTt3=funTt44(Tt1+h*kTt2/2,Ma1,T1,kMa,kT)
   kPt3=funPt44(Pt1+h*kPt2/2,Ma1,p1,kMa,kp)
   ka3=funa44(a1+h*ka2/2,T1,kT)
   kAcA3=funAcA44(AcA1+h*kAcA2/2,Ma1,p1,Tt1,Ac1,kp,kTt,kAc3)
   
   kp4=funp44(p1+h*kp3,rho1,T1,krho,kT)
   kMa4=funMa44(Ma1+h*kMa3,V1,T1,kV,kT)
   krho4=funrho44(rho1+h*krho3,Ac1,V1,kAc,kV)
   kV4=funV44(V1+h*kV3,Ma1,p1,kp)
   kT4=funT44(T1+h*kT3,Ma1,V1,kV)
   kF4=funF44(F1+h*kF3,Ma1,p1,kMa,kp)
   kTt4=funTt44(Tt1+h*kTt3,Ma1,T1,kMa,kT)
   kPt4=funPt44(Pt1+h*kPt3,Ma1,p1,kMa,kp)
   ka4=funa44(a1+h*ka3,T1,kT)
   kAcA4=funAcA44(AcA1+h*kAcA3,Ma1,p1,Tt1,Ac1,kp,kTt,kAc4)
   
   PR(i,2)=p1+h*(kp1+2*kp2+2*kp3+kp4)/6
   PR(i,3)=Ma1+h*(kMa1+2*kMa2+2*kMa3+kMa4)/6
   PR(i,4)=rho1+h*(krho1+2*krho2+2*krho3+krho4)/6
   PR(i,5)=V1+h*(kV1+2*kV2+2*kV3+kV4)/6
   PR(i,6)=T1+h*(kT1+2*kT2+2*kT3+kT4)/6
   PR(i,7)=F1+h*(kF1+2*kF2+2*kF3+kF4)/6
   PR(i,8)=Tt1+h*(kTt1+2*kTt2+2*kTt3+kTt4)/6
   PR(i,9)=Pt1+h*(kPt1+2*kPt2+2*kPt3+kPt4)/6
   PR(i,11)=a1+h*(ka1+2*ka2+2*ka3+ka4)/6
   PR(i,12)=Ac1+h*(kAc1+2*kAc2+2*kAc3+kAc4)/6
   
   ! if the dlow condition can not support the BL seperation
   if (PR(i,12) > PR(i,13)) then
      PR(i,12)=PR(i,13)
   end if
   
!   calculate G(x)
   PR(i,17)=funGx44(Ma1,Tt1,Ac1,kTt,kAc)
   
end subroutine RK44

! function funp44
function funp44(p,rho,T,krho,kT)
   use VariableDef
   implicit none
   real*8::p,rho,T,krho,kT,funp44
   funp44=krho/rho+kT/T-dWW
   funp44=p*funp44
end function funp44

! function funMa44
function funMa44(Ma,V,T,kV,kT)
   use VariableDef
   implicit none
   real*8::Ma,V,T,kV,kT,funMa44
   funMa44=2*kV/V+dWW-dGG-kT/T
   funMa44=Ma*funMa44/2
end function funMa44

! function funrho44
function funrho44(rho,Ac,V,kAc,kV)
   use VariableDef
   implicit none
   real*8::rho,Ac,V,kAc,kV,funrho44
   funrho44=dmdotm-kAc/Ac-kV/V
   funrho44=rho*funrho44
end function funrho44

! function funV44
function funV44(V,Ma,p,kp)
   use VariableDef
   implicit none
   real*8::V,Ma,p,kp,funV44
   real*8::temp
   temp=gamma*Ma*Ma
   funV44=-kp/p-temp*(ffD+dDxd)/2-temp*(1-y)*dmdotm
   funV44=V*funV44/temp
end function funV44

! function funT44
function funT44(T,Ma,V,kV)
   use VariableDef
   implicit none
   real*8::T,Ma,V,kV,funT44
   real*8::temp
   temp=(gamma-1)*Ma*Ma*kV/V
   funT44=T*(dQWxH/(Cp*T)-temp)
end function funT44

! function funF44
function funF44(F,Ma,p,kMa,kp)
   use VariableDef
   implicit none
   real*8::F,Ma,p,kMa,kp,funF44
   real*8::temp
   temp=gamma*Ma*Ma
   temp=temp/(1+temp)
   funF44=dAA+kp/p+temp*2*kMa/Ma+temp*dGG
   funF44=funF44*F
end function funF44

! function funTt44
function funTt44(Tt,Ma,T,kMa,kT)
   use VariableDef
   implicit none
   real*8::Tt,Ma,T,kMa,kT,funTt44
   real*8::temp
   temp=(gamma-1)*Ma*Ma/2
   temp=temp/(1+temp)
   funTt44=kT/T+temp*2*kMa/Ma
   funTt44=Tt*funTt44
end function funTt44

! function funPt44
function funPt44(Pt,Ma,p,kMa,kp)
   use VariableDef
   implicit none
   real*8::Pt,Ma,p,kMa,kp,funPt44
   real*8::temp
   temp=(gamma-1)*Ma*Ma/2
   temp=(temp+Ma*Ma/2)/(1+temp)
   funPt44=kp/p+temp*2*kMa/Ma
   funPt44=Pt*funPt44
end function funPt44

! function funa44
function funa44(a,T,kT)
   use VariableDef
   implicit none
   real*8::a,T,kT,funa44
   funa44=(dGG+kT/T-dWW)/2
   funa44=a*funa44
end function funa44

! function funAc
function funAcA44(AcA,Ma,p,Tt,Ac,kp,kTt,kAc)
   use VariableDef
   implicit none
   real*8::AcA,Ma,p,Tt,Ac,kp,kTt,kAc,funAcA44
   real*8::temp0,temp1,temp2,temp3
   real*8::tempp,tempD,tempTt,tempm
   temp0=Ma*Ma
   temp1=gamma*temp0
   temp2=temp1-temp0
   temp3=temp2/2
   tempp=(1-temp0*(1-gamma*(1-AcA)))/(temp1*AcA)
   tempD=(1+temp2)/(2*AcA)
   tempTt=1+temp3
   tempm=1+(1+temp2)*(1-y)
   funAcA44=tempp*kp/p+tempD*(ffD+dDxd)+tempTt*kTt/Tt+tempm*dmdotm+temp3*dGG-tempTt*dWW
   kAc=funAcA44*Ac
   funAcA44=AcA*funAcA44
end function funAcA44

! function funGx
function funGx44(Ma,Tt,Ac,kTt,kAc)
   use VariableDef
   implicit none
   real*8::Ma,Tt,Ac,kTt,kAc,funGx44
   real*8::temp0,temp1,temp2
   temp0=Ma*Ma
   temp1=gamma*temp0
   temp2=(temp1-temp0)/2
   funGx44=2*kAc/Ac-(1+temp1)*kTt/Tt-temp1*(ffD+dDxd)-2*(1+temp1*(1-y))*dmdotm
!   funGx44=2*dAA-(1+temp1)*kTt/Tt-temp1*(ffD+dDxd)-2*(1+temp1*(1-y))*dmdotm
   funGx44=funGx44*(1+temp2)+(1+temp1)*dWW-(temp0-1)*dGG
   funGx44=temp0*funGx44
end function funGx44
