!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The subroutine is the variable step Runge-Kutta-Ferlberg
!    method kernel used in the integrate method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Coded by       : Black Ghost
! Created Date   : 2016-01-08
! Completed Date : 2016-01-27
! Revised Date   :
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ref. John H. Mathews and Kurtis D Fink. Numerical Methods
!   Using MATLAB, 4th Edition. ÖÜè´£¬ ³ÂÓå£¬ Ç®·½ µÈÒë. 2012.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note :
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine RK54(i,p1,Ma1,rho1,V1,T1,F1,Tt1,Pt1,a1,AcA1,Ma4)
   use VariableDef
   implicit none
   integer::i,k
   real*8::p1,Ma1,rho1,V1,T1,F1,Tt1,Pt1,a1,AcA1,Ac1
   real*8::kp1,kMa1,krho1,kV1,kT1,kF1,kTt1,kPt1,ka1,kAcA1,kAc1
   real*8::kp2,kMa2,krho2,kV2,kT2,kF2,kTt2,kPt2,ka2,kAcA2,kAc2
   real*8::kp3,kMa3,krho3,kV3,kT3,kF3,kTt3,kPt3,ka3,kAcA3,kAc3
   real*8::kp4,kMa4,krho4,kV4,kT4,kF4,kTt4,kPt4,ka4,kAcA4,kAc4
   real*8::kp5,kMa5,krho5,kV5,kT5,kF5,kTt5,kPt5,ka5,kAcA5,kAc5
   real*8::kp6,kMa6,krho6,kV6,kT6,kF6,kTt6,kPt6,ka6,kAcA6,kAc6
   real*8::p4,Ma4,rho4,V4,T4,F4,Tt4,Pt4,a4,AcA4
   real*8::kp,kMa,krho,kV,kT,kTt,kAc
   real*8::funp54,funMa54,funrho54,funV54,funT54,funF54,funTt54,funPt54,funa54,funAcA54,funGx54
   
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
   
   kp1=h*funp54(p1,rho1,T1,krho,kT)
   kMa1=h*funMa54(Ma1,V1,T1,kV,kT)
   krho1=h*funrho54(rho1,Ac1,V1,kAc,kV)
   kV1=h*funV54(V1,Ma1,p1,kp)
   kT1=h*funT54(T1,Ma1,V1,kV)
   kF1=h*funF54(F1,Ma1,Ac1,p1,kMa,kAc,kp)
   kTt1=h*funTt54(Tt1,Ma1,T1,kMa,kT)
   kPt1=h*funPt54(Pt1,Ma1,p1,kMa,kp)
   ka1=h*funa54(a1,T1,kT)
   kAcA1=h*funAcA54(AcA1,Ma1,p1,Tt1,Ac1,kp,kTt,kAc1)
   kAc1=kAc1*h
   
   kp2=h*funp54(p1+kp1/4,rho1,T1,krho,kT)
   kMa2=h*funMa54(Ma1+kMa1/4,V1,T1,kV,kT)
   krho2=h*funrho54(rho1+krho1/4,Ac1,V1,kAc,kV)
   kV2=h*funV54(V1+kV1/4,Ma1,p1,kp)
   kT2=h*funT54(T1+kT1/4,Ma1,V1,kV)
   kF2=h*funF54(F1+kF1/4,Ma1,Ac1,p1,kMa,kAc,kp)
   kTt2=h*funTt54(Tt1+kTt1/4,Ma1,T1,kMa,kT)
   kPt2=h*funPt54(Pt1+kPt1/4,Ma1,p1,kMa,kp)
   ka2=h*funa54(a1+ka1/4,T1,kT)
   kAcA2=h*funAcA54(AcA1+kAcA1/4,Ma1,p1,Tt1,Ac1,kp,kTt,kAc2)
   kAc2=kAc2*h
   
   kp3=h*funp54(p1+3*kp1/32+9*kp2/32,rho1,T1,krho,kT)
   kMa3=h*funMa54(Ma1+3*kMa1/32+9*kMa2/32,V1,T1,kV,kT)
   krho3=h*funrho54(rho1+3*krho1/32+9*krho2/32,Ac1,V1,kAc,kV)
   kV3=h*funV54(V1+3*kV1/32+9*kV2/32,Ma1,p1,kp)
   kT3=h*funT54(T1+3*kT1/32+9*kT2/32,Ma1,V1,kV)
   kF3=h*funF54(F1+3*kF1/32+9*kF2/32,Ma1,Ac1,p1,kMa,kAc,kp)
   kTt3=h*funTt54(Tt1+3*kTt1/32+9*kTt2/32,Ma1,T1,kMa,kT)
   kPt3=h*funPt54(Pt1+3*kPt1/32+9*kPt2/32,Ma1,p1,kMa,kp)
   ka3=h*funa54(a1+3*ka1/32+9*ka2/32,T1,kT)
   kAcA3=h*funAcA54(AcA1+3*kAcA1/32+9*kAcA2/32,Ma1,p1,Tt1,Ac1,kp,kTt,kAc3)
   kAc3=kAc3*h
   
   kp4=h*funp54(p1+1932*kp1/2197-7200*kp2/2197+7296*kp3/2197,rho1,T1,krho,kT)
   kMa4=h*funMa54(Ma1+1932*kMa1/2197-7200*kMa2/2197+7296*kMa3/2197,V1,T1,kV,kT)
   krho4=h*funrho54(rho1+1932*krho1/2197-7200*krho2/2197+7296*krho3/2197,Ac1,V1,kAc,kV)
   kV4=h*funV54(V1+1932*kV1/2197-7200*kV2/2197+7296*kV3/2197,Ma1,p1,kp)
   kT4=h*funT54(T1+1932*kT1/2197-7200*kT2/2197+7296*kT3/2197,Ma1,V1,kV)
   kF4=h*funF54(F1+1932*kF1/2197-7200*kF2/2197+7296*kF3/2197,Ma1,Ac1,p1,kMa,kAc,kp)
   kTt4=h*funTt54(Tt1+1932*kTt1/2197-7200*kTt2/2197+7296*kTt3/2197,Ma1,T1,kMa,kT)
   kPt4=h*funPt54(Pt1+1932*kPt1/2197-7200*kPt2/2197+7296*kPt3/2197,Ma1,p1,kMa,kp)
   ka4=h*funa54(a1+1932*ka1/2197-7200*ka2/2197+7296*ka3/2197,T1,kT)
   kAcA4=h*funAcA54(AcA1+1932*kAcA1/2197-7200*kAcA2/2197+7296*kAcA3/2197,Ma1,p1,Tt1,Ac1,kp,kTt,kAc4)
   kAc4=kAc4*h
   
   kp5=h*funp54(p1+439*kp1/216-8*kp2+3680*kp3/513-845*kp4/4104,rho1,T1,krho,kT)
   kMa5=h*funMa54(Ma1+439*kMa1/216-8*kMa2+3680*kMa3/513-845*kMa4/4104,V1,T1,kV,kT)
   krho5=h*funrho54(rho1+439*krho1/216-8*krho2+3680*krho3/513-845*krho4/4104,Ac1,V1,kAc,kV)
   kV5=h*funV54(V1+439*kV1/216-8*kV2+3680*kV3/513-845*kV4/4104,Ma1,p1,kp)
   kT5=h*funT54(T1+439*kT1/216-8*kT2+3680*kT3/513-845*kT4/4104,Ma1,V1,kV)
   kF5=h*funF54(F1+439*kF1/216-8*kF2+3680*kF3/513-845*kF4/4104,Ma1,Ac1,p1,kMa,kAc,kp)
   kTt5=h*funTt54(Tt1+439*kTt1/216-8*kTt2+3680*kTt3/513-845*kTt4/4104,Ma1,T1,kMa,kT)
   kPt5=h*funPt54(Pt1+439*kPt1/216-8*kPt2+3680*kPt3/513-845*kPt4/4104,Ma1,p1,kMa,kp)
   ka5=h*funa54(a1+439*ka1/216-8*ka2+3680*ka3/513-845*ka4/4104,T1,kT)
   kAcA5=h*funAcA54(AcA1+439*kAcA1/216-8*kAcA2+3680*kAcA3/513-845*kAcA4/4104,Ma1,p1,Tt1,Ac1,kp,kTt,kAc5)
   kAc5=kAc5*h
   
   kp6=h*funp54(p1-8*kp1/27+2*kp2-3544*kp3/2565+1859*kp4/4104-11*kp5/40,rho1,T1,krho,kT)
   kMa6=h*funMa54(Ma1-8*kMa1/27+2*kMa2-3544*kMa3/2565+1859*kMa4/4104-11*kMa5/40,V1,T1,kV,kT)
   krho6=h*funrho54(rho1-8*krho1/27+2*krho2-3544*krho3/2565+1859*krho4/4104-11*krho5/40,Ac1,V1,kAc,kV)
   kV6=h*funV54(V1-8*kV1/27+2*kV2-3544*kV3/2565+1859*kV4/4104-11*kV5/40,Ma1,p1,kp)
   kT6=h*funT54(T1-8*kT1/27+2*kT2-3544*kT3/2565+1859*kT4/4104-11*kT5/40,Ma1,V1,kV)
   kF6=h*funF54(F1-8*kF1/27+2*kF2-3544*kF3/2565+1859*kF4/4104-11*kF5/40,Ma1,Ac1,p1,kMa,kAc,kp)
   kTt6=h*funTt54(Tt1-8*kTt1/27+2*kTt2-3544*kTt3/2565+1859*kTt4/4104-11*kTt5/40,Ma1,T1,kMa,kT)
   kPt6=h*funPt54(Pt1-8*kPt1/27+2*kPt2-3544*kPt3/2565+1859*kPt4/4104-11*kPt5/40,Ma1,p1,kMa,kp)
   ka6=h*funa54(a1-8*ka1/27+2*ka2-3544*ka3/2565+1859*ka4/4104-11*ka5/40,T1,kT)
   kAcA6=h*funAcA54(AcA1-8*kAcA1/27+2*kAcA2-3544*kAcA3/2565+1859*kAcA4/4104-11*kAcA5/40,Ma1,p1,Tt1,Ac1,kp,kTt,kAc6)
   kAc6=kAc6*h
   
   p4=p1+25*kp1/216+1408*kp3/2565+2197*kp4/4104-kp5/5
   Ma4=Ma1+25*kMa1/216+1408*kMa3/2565+2197*kMa4/4104-kMa5/5
   rho4=rho1+25*krho1/216+1408*krho3/2565+2197*krho4/4104-krho5/5
   V4=V1+25*kV1/216+1408*kV3/2565+2197*kV4/4104-kV5/5
   T4=T1+25*kT1/216+1408*kT3/2565+2197*kT4/4104-kT5/5
   F4=F1+25*kF1/216+1408*kF3/2565+2197*kF4/4104-kF5/5
   Tt4=Tt1+25*kTt1/216+1408*kTt3/2565+2197*kTt4/4104-kTt5/5
   Pt4=Pt1+25*kPt1/216+1408*kPt3/2565+2197*kPt4/4104-kPt5/5
   a4=a1+25*ka1/216+1408*ka3/2565+2197*ka4/4104-ka5/5
   AcA4=AcA1+25*kAcA1/216+1408*kAcA3/2565+2197*kAcA4/4104-kAcA5/5
!   AcA4=AcA4*PR(i,13)
   
   PR(i,2)=p1+16*kp1/135+6656*kp3/12825+28561*kp4/56430-9*kp5/50+2*kp6/55
   PR(i,3)=Ma1+16*kMa1/135+6656*kMa3/12825+28561*kMa4/56430-9*kMa5/50+2*kMa6/55
   PR(i,4)=rho1+16*krho1/135+6656*krho3/12825+28561*krho4/56430-9*krho5/50+2*krho6/55
   PR(i,5)=V1+16*kV1/135+6656*kV3/12825+28561*kV4/56430-9*kV5/50+2*kV6/55
   PR(i,6)=T1+16*kT1/135+6656*kT3/12825+28561*kT4/56430-9*kT5/50+2*kT6/55
   PR(i,7)=F1+16*kF1/135+6656*kF3/12825+28561*kF4/56430-9*kF5/50+2*kF6/55
   PR(i,8)=Tt1+16*kTt1/135+6656*kTt3/12825+28561*kTt4/56430-9*kTt5/50+2*kTt6/55
   PR(i,9)=Pt1+16*kPt1/135+6656*kPt3/12825+28561*kPt4/56430-9*kPt5/50+2*kPt6/55
   PR(i,11)=a1+16*ka1/135+6656*ka3/12825+28561*ka4/56430-9*ka5/50+2*ka6/55
   PR(i,12)=Ac1+16*kAc1/135+6656*kAc3/12825+28561*kAc4/56430-9*kAc5/50+2*kAc6/55
   
   ! if the dlow condition can not support the BL seperation
   if (PR(i,12) > PR(i,13)) then
      PR(i,12)=PR(i,13)
   end if
   
!   calculate G(x)
   PR(i,17)=funGx54(Ma1,Tt1,Ac1,kTt,kAc)
   
end subroutine RK54

! function funp54
function funp54(p,rho,T,krho,kT)
   use VariableDef
   implicit none
   real*8::p,rho,T,krho,kT,funp54
   funp54=krho/rho+kT/T-dWW
   funp54=p*funp54
end function funp54

! function funMa54
function funMa54(Ma,V,T,kV,kT)
   use VariableDef
   implicit none
   real*8::Ma,V,T,kV,kT,funMa54
   funMa54=2*kV/V+dWW-dGG-kT/T
   funMa54=Ma*funMa54/2
end function funMa54

! function funrho54
function funrho54(rho,Ac,V,kAc,kV)
   use VariableDef
   implicit none
   real*8::rho,Ac,V,kAc,kV,funrho54
   funrho54=dmdotm-kAc/Ac-kV/V
   funrho54=rho*funrho54
end function funrho54

! function funV54
function funV54(V,Ma,p,kp)
   use VariableDef
   implicit none
   real*8::V,Ma,p,kp,funV54
   real*8::temp
   temp=gamma*Ma*Ma
   funV54=-kp/p-temp*(ffD+dDxd)/2-temp*(1-y)*dmdotm
   funV54=V*funV54/temp
end function funV54

! function funT54
function funT54(T,Ma,V,kV)
   use VariableDef
   implicit none
   real*8::T,Ma,V,kV,funT54
   real*8::temp
   temp=(gamma-1)*Ma*Ma*kV/V
   funT54=T*(dQWxH/(Cp*T)-temp)
end function funT54

! function funF54
function funF54(F,Ma,Ac,p,kMa,kAc,kp)
   use VariableDef
   implicit none
   real*8::F,Ma,Ac,p,kMa,kAc,kp,funF54
   real*8::temp
   temp=gamma*Ma*Ma
   temp=temp/(1+temp)
   funF54=dAA+kp/p+temp*2*kMa/Ma+temp*dGG
   funF54=funF54*F
end function funF54

! function funTt54
function funTt54(Tt,Ma,T,kMa,kT)
   use VariableDef
   implicit none
   real*8::Tt,Ma,T,kMa,kT,funTt54
   real*8::temp
   temp=(gamma-1)*Ma*Ma/2
   temp=temp/(1+temp)
   funTt54=kT/T+temp*2*kMa/Ma
   funTt54=Tt*funTt54
end function funTt54

! function funPt54
function funPt54(Pt,Ma,p,kMa,kp)
   use VariableDef
   implicit none
   real*8::Pt,Ma,p,kMa,kp,funPt54
   real*8::temp
   temp=(gamma-1)*Ma*Ma/2
   temp=(temp+Ma*Ma/2)/(1+temp)
   funPt54=kp/p+temp*2*kMa/Ma
   funPt54=Pt*funPt54
end function funPt54

! function funa54
function funa54(a,T,kT)
   use VariableDef
   implicit none
   real*8::a,T,kT,funa54
   funa54=(dGG+kT/T-dWW)/2
   funa54=a*funa54
end function funa54

! function funAc
function funAcA54(AcA,Ma,p,Tt,Ac,kp,kTt,kAc)
   use VariableDef
   implicit none
   real*8::AcA,Ma,p,Tt,Ac,kp,kTt,kAc,funAcA54
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
   funAcA54=tempp*kp/p+tempD*(ffD+dDxd)+tempTt*kTt/Tt+tempm*dmdotm+temp3*dGG-tempTt*dWW
   kAc=funAcA54*Ac
   funAcA54=AcA*funAcA54
end function funAcA54

! function funGx
function funGx54(Ma,Tt,Ac,kTt,kAc)
   use VariableDef
   implicit none
   real*8::Ma,Tt,Ac,kTt,kAc,funGx54
   real*8::temp0,temp1,temp2
   temp0=Ma*Ma
   temp1=gamma*temp0
   temp2=(temp1-temp0)/2
   funGx54=2*kAc/Ac-(1+temp1)*kTt/Tt-temp1*(ffD+dDxd)-2*(1+temp1*(1-y))*dmdotm
!   funGx=2*dAA-(1+temp1)*kTt/Tt-temp1*(ffD+dDxd)-2*(1+temp1*(1-y))*dmdotm
   funGx54=funGx54*(1+temp2)+(1+temp1)*dWW-(temp0-1)*dGG
   funGx54=temp0*funGx54
end function funGx54
