!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program is coded for one-dimensionally simulates
!     dual-mode scramjet combustor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Coded by       : Black Ghost
! Created Date   : 2016-01-06
! Completed Date : 2016-01-11
! Revised Date   :
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ref. :
!   1. Ascher H. Shapiro. The Dynamics and Thermodynamics of Compressible Fluid Flow, VOL. I, Chapter 8.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Version :
!    0.0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Contains :
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note :
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program DualMode
   use VariableDef
   implicit none
   
!## input initial values ##!
   call Input()

!## call for simulation ##!
   call Integrate()

!## output the calculated parameters ##!
   call Output()

end program DualMode
