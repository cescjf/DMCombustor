!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the subroutine is coded for modeling net heat
!    transfered to main flow, the work extracted
!    from main flow and dH which defined by 
!    Eq. 8.15 of Shapiro's book when properties
!    calculation of dual mode scramjet
!    dQ-dWx+dH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Coded by       : Black Ghost
! Created Date   : 2016-01-07
! Completed Date : 2016-01-13
! Revised Date   : 
!    2016-02-02      add empirical ignition delay time from
!       Heiser & Pratt's book, pp. 326
!    2016-02-25      add empirical mixing efficiency from the
!       report of G. B. Northam & G. Y. Anderson
!    2016-06-08      add emperical mixing length from the 
!       paper of Cristian Birzer and Con J. Doolan
!    2016-06-08      add  emperical mixing length from G. B. 
!       Northam & G. Y. Anderson and M .V. Pulsonetti, J. 
!       Erdos & K. Early
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note :
!   1. assume ignition delayed length is 10*Dinj for first
!      injection and 1*Dinj for others
!   2. assume the length of mixing is equals to the length 
!      of reaction
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! a function of exothermic chemical reaction is required, 
!    which stretching from ignition to infinite and arrives
!    its maximum at some distance downstream ignition
 
subroutine dQWxHcal(i)
   use VariableDef
   implicit none
   integer::i,k
   real*8::entam,entaC,lambda,Dw,Lind,Lrea,xphi,entam0,entam90,entamcur
   real*8::xReact,Q,dQ,Wx,T0,dWx,V0,hl,hv,dH
   real*8::gammaf,af,uf,Mc,fc,rhof,xbar
   
   entam=1.0                                 ! mixing efficiency
   entaC=0.9                                 ! combustion efficiency
   lambda=24                                 ! thermo-conductivity
   Dw=0.01                                   ! wall depth
   Lind=4.5e-9*(101325/PR(i-1,3))*exp(10000/PR(i-1,6))
   Lind=PR(i-1,5)*Lind                       ! ignition delay length
   Lrea=170*Dinj	                           ! length of reaction region
   
   T0=PR(i-1,6)
   xReact=PR(i,1)
   V0=PR(i-1,5)
   
!   assume ideal gas ^H=Cp*^T, h=^H/mass
!   if mfdotl not equals to zero, then hl-hv = Cpf*(T0-Tw)/mfdotl
   hl=0.0                                    ! liquid fuel enthapy
   hv=0.0                                    ! evaporated liquid fuel enthapy
   dQ=0.0                                    ! initiate dQ
   
!   dQ, combustion heat release
!   assume gaseous hydrogen is fuel
!   assume combustion occurs between Lind and Lind+Lrea for first injection
!      and between Dinj and Dinj+Lrea for others (51 cm length for combustion)
!#  require the second injector locate after x=Lind from first one #
!   2H2(g)+O2(g)-->2H2O(g), ^H=-483.6 kJ/mol(O2)
   do k=1,Nm
!      mixing length, Cristian Birzer and Con J. Doolan
!      gammaf=28.82/20.18
!      af=sqrt(gammaf*8314.4598*Tf/2)   !sound speed of fuel H2
!      uf=PR(i-1,5)*MASS(k,3)/cosd(MASS(k,4))
!      Mc=(uf-PR(i-1,5))/(af-PR(i-1,11))
!      fc=0.25+0.75*exp(-3.0*Mc*Mc)
!      rhof=4*mfdot/(af*Dinj*Dinj)
!      Lrea=Dinj*390*sqrt(rhof*uf/(PR(i-1,4)*PR(i-1,5)))/fc   !390 for strut injector
!      xbar=(PR(i,1)-MASS(k,1))/Lrea
!      entam=1.06492*(1-exp(-(3.69639*xbar)**0.80586))

      ! mixing length of G. B. Northam & G. Y. Anderson and M .V. Pulsonetti, J. Erdos & K. Early
      Lrea=50*Hd
      
!      entam specifying, G. B. Northam & G. Y. Anderson
      if (MASS(k,2) <= 1.0) then
         xphi=Lrea*0.179*exp(1.72*MASS(k,2))
      else
         xphi=Lrea*3.33*exp(-1.204*MASS(k,2))
      end if
      entam0=(PR(i,1)-MASS(k,1))/Lrea
      entam90=1.01+0.176*log((PR(i,1)-MASS(k,1))/xphi)
      entamcur=entam0+(entam90-entam0)*MASS(k,4)/90
      entam=entamcur-entampre
      entampre=entamcur

      if ((xReact > (MASS(1,1)+Lind+0.1*h)) .and. (xReact <= (MASS(1,1)+Lind+Lrea+0.1*h)) .and. (k == 1)) then
         Q=entam*entaC*MASS(1,2)*PR(1,15)*fst*483600/4
         dQ=Q/(PR(i,1)-PR(i-1,1))
         cycle
      else if ((xReact > (MASS(k,1)+Dinj+0.1*h)) .and. (xReact <= (MASS(k,1)+Dinj+Lrea+0.1*h)) .and. (Nm >= 2) .and. &
         &(k >= 2)) then
         Q=entam*entaC*MASS(k,2)*PR(1,15)*fst*483600/4
         dQ=dQ+Q/(PR(i,1)-PR(i-1,1))
         cycle
      else if (xReact <= (MASS(1,1)+Lind+0.1*h)) then
         dQ=dQ+0.0
         exit
      else if ((Nm == 1) .and. (xReact > (MASS(1,1)+Lind+Lrea+0.1*h))) then
         dQ=dQ+0.0
         exit
      else if ((xReact > (MASS(1,1)+Lind+Lrea+0.1*h)) .and. (xReact <= (MASS(2,1)+Dinj+0.1*h)) .and. (Nm == 2)) then
         dQ=dQ+0.0
         cycle
      else if ((xReact > (MASS(k,1)+Dinj+Lrea+0.1*h)) .and. (xReact <= (MASS(k+1,1)+Dinj+0.1*h)) .and. (Nm > 2) .and. &
         &(k >= 2) .and. ((k+1) <= Nm)) then
         dQ=dQ+0.0
         cycle
      else if ((k == Nm) .and. (xReact > (MASS(k,1)+Dinj+Lrea+0.1*h))) then
         dQ=dQ+0.0
         exit
      end if
      
   end do
   
!   dWx, work and heat transfer
   Wx=0.0+lambda*Dh*(PR(i,1)-PR(i-1,1))*(T0-Tw)/Dw
   dWx=Wx/(PR(i,1)-PR(i-1,1))
   
!   dH, heat changer caused by fuel injection, which is defined by the second 
!      and third terms of Eq. 8.15 (Shapiro)
!   Cp(H2)=14.30 kJ/(kg*K) according to "百度百科"
!   assume the temperature of fuel equals to the temperature of wall Tw
   dH=-(Cpf*(T0-Tf)+V0*V0/2)*dmfdotgm-(hl-hv+(V0*V0-Vl*Vl)/2)*dmfdotlm
   
   dQWxH=dQ-dWx+dH
   
end subroutine dQWxHcal
