!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

   module current_meters

!BOP
! !MODULE: current_meters
!
! !DESCRIPTION:
! This module collects data to compare with current meter data.
! NOTE: this module currently does not work. old CM-5 routines
! are appended but must be re-done.
!
! !REVISION HISTORY:
! CVS:$Id: current_meters.F90,v 1.4 2002/04/16 15:22:53 pwjones Exp $
! CVS:$Name: POP_2_0_1 $

! !USES:

   use kinds_mod

   implicit none
   private
   save

!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_current_meters
! !INTERFACE:

 subroutine init_current_meters

! !DESCRIPTION:
! Initializes all necessary variables for current meter diagnostics.
!
! !REVISION HISTORY:
! same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!EOC

 end subroutine init_current_meters

!***********************************************************************

 end module current_meters

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
