!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

   module drifters

!BOP
! !MODULE: drifters
!
! !DESCRIPTION:
! This module contains routines necessary for moving drifters.
! NOTE: this module currently does not really exist - this version
! has old CM-5 code and will be replaced with some new code when
! Mat is ready to add it in...
!
! !REVISION HISTORY:
! CVS:$Id: drifters.F90,v 1.3 2002/02/26 22:49:04 pwjones Exp $
! CVS:$Name: POP_2_0_1 $
!
! !USES

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
! !IROUTINE: init_drifters
! !INTERFACE:

 subroutine init_drifters

! !DESCRIPTION:
! Initializes all variables for drifter diagnostics.
!
! !REVISION HISTORY:
! same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!EOC

 end subroutine init_drifters

!***********************************************************************

 end module drifters

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
