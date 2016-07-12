!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the subroutine is coded for output paramenters
!    generated in the program DualMode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Coded by       : Black Ghost
! Created Date   : 2016-01-07
! Completed Date : 2016-01-27
! Revised Date   :
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note :
!   1. all angles in rad in calculations, but 
!      should be convert to in degree for output
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Output()
   use VariableDef
   implicit none
   integer::i,k

  ! output
   open(unit=201,file='Results.txt')
   write(201,200) Nend,NA,Nm,cf,Dinj,L
   write(201,201)
   do k=1,NA
      write(201,202) AREA(k,1),AREA(k,2)*180/pi
   end do
   write(201,203)
   do k=1,Nm
      write(201,204) MASS(k,1),MASS(k,2),MASS(k,3),MASS(k,4)
   end do
   write(201,205)
   do i=1,Nend
      write(201,206) PR(i,1),PR(i,2),PR(i,3),PR(i,4),PR(i,5),PR(i,6),PR(i,7),PR(i,8),&
         &PR(i,9),PR(i,10),PR(i,11),PR(i,12),PR(i,13),PR(i,14),PR(i,15),PR(i,16),PR(i,17)
   end do
   
   200 format(/,T10,'Dual Mode Scramjet',//,2x,'N',T11,'NA',T17,'Nm',T23,'cf',T34,'Dinj',&
      &T45,'L',/,2x,I5,T11,I1,T17,I1,T23,F7.5,T34,F7.5,T45,F5.3)
   201 format(T10,'Area changes along axis',/,2x,'xA',T11,'theta')
   202 format(2x,F5.3,T11,F5.2)
   203 format(T10,'Injections along axis',/,2x,'xm',T11,'psi',T19,'y',T27,'angle')
   204 format(2x,F5.3,T11,F4.2,T19,F4.2,T27,F5.2)
   205 format(//,2x,'x',T12,'p',T24,'Ma',T33,'rho',T45,'V',T56,'T',T66,'F',T76,'Tt',T86,&
      &'Pt',T98,'s',T109,'a',T119,'Ac',T132,'A',T145,'gamma',T157,'mdot',T168,'W',T178,'G',/)
   206 format(2x,F7.4,T12,F9.1,T24,F6.4,T33,F9.6,T45,F8.3,T56,F7.2,T66,F7.2,T76,F7.2,T86,&
      &F9.1,T98,F8.2,T109,F7.3,T119,F10.8,T132,F10.8,T145,F9.7,T157,F8.5,T168,F7.4,T178,F10.5)
   close(201)
   call cpu_time(time_end)
   write(*,*) time_end
   write(*,*) time_end-time_begin

end subroutine Output
