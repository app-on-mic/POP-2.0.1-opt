!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module hydro_sections

!BOP
! !MODULE: hydro_sections
!
! !DESCRIPTION:
! This module computes data along hydrographic sections to compare
! with cruise data.
! NOTE: currently does not work. routines appended below are from old
! CM-5 version and must be re-done.
!
! !REVISION HISTORY:
! CVS:$Id: hydro_sections.F90,v 1.4 2002/04/16 15:22:56 pwjones Exp $
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
! !IROUTINE: init_hydro_sections
! !INTERFACE:

 subroutine init_hydro_sections

! !DESCRIPTION:
! Initializes all variables to be used for hydrographic sections.
!
! !REVISION HISTORY:
! same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!EOC

 end subroutine init_hydro_sections

!***********************************************************************

 end module hydro_sections

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
