!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module forcing_coupled

!MODULE: forcing_coupled
! !DESCRIPTION:
! This module contains all the routines necessary for coupling POP to
! atmosphere and sea ice models using the NCAR CCSM flux coupler. To
! enable the routines in this module, the coupled ifdef option must
! be specified during the make process.
!
! !REVISION HISTORY:
! CVS:$Id: forcing_coupled.F90,v 1.9 2003/02/24 16:18:43 pwjones Exp $
! CVS:$Name: POP_2_0_1 $
!
! !USES:

   use kinds_mod
   use domain_size
   use domain
   use communicate
   use constants
   use broadcast
   use io
   use time_management
   use grid
   use prognostic
   use exit_mod
   use ice, only: tfreez, get_ice_flux
   use forcing_shf
   use forcing_sfwf
   !NCAR use qflux_mod
   !NCAR use ms_balance
   use timers
   use gather_scatter, only: scatter_global, gather_global
   use global_reductions, only: global_sum_prod
   !NCAR use shr_sys_mod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_coupled, &
             set_coupled_forcing

! !PUBLIC DATA MEMBERS:

   logical (log_kind), public :: &
      lcoupled, &! flag for coupled forcing
      ldiag_cpl ! flag for computing coupled diagnostics

   integer (int_kind), public :: &
      coupled_freq_iopt, &! coupler frequency option
      coupled_freq ! frequency of coupling

!EOP
!BOC
!-----------------------------------------------------------------------
!
! module variables
!
!-----------------------------------------------------------------------
!EOC
!***********************************************************************
      contains
!***********************************************************************
!BOP
! !IROUTINE: init_coupled
! !INTERFACE:
 subroutine init_coupled(SMF, SMFT, STF, SHF_QSW, lsmft_avail)
! !DESCRIPTION:
! This routine sets up everything necessary for coupling with
! the NCAR flux coupler.
!
! !REVISION HISTORY:
! same as coupled
! !INPUT/OUTPUT PARAMETERS:
   real (r8), dimension(nx_block,ny_block,2,max_blocks_clinic), &
      intent(inout) :: &
      SMF, &! surface momentum fluxes (wind stress)
      SMFT ! surface momentum fluxes at T points
   real (r8), dimension(nx_block,ny_block,nt,max_blocks_clinic), &
      intent(inout) :: &
      STF ! surface tracer fluxes
   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      intent(inout) :: &
      SHF_QSW ! penetrative solar heat flux
   logical (log_kind), intent(inout) :: &
      lsmft_avail ! true if SMFT is an available field
!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------
   character (char_len) :: &
      coupled_freq_opt
   namelist /coupled_nml/ coupled_freq_opt, coupled_freq
   integer (int_kind) :: &
      i,j,k, &! dummy loop index
      ncouple_per_day, &! num of coupler comms per day
      nml_error ! namelist i/o error flag
!-----------------------------------------------------------------------
!
! read coupled_nml namelist to start coupling and determine
! coupling frequency
!
!-----------------------------------------------------------------------
   lcoupled = .false.
   coupled_freq_opt = 'never'
   coupled_freq_iopt = freq_opt_never
   coupled_freq = 100000
   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error = 1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=coupled_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif
   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR: reading coupled_nml')
   endif
   if (my_task == master_task) then
      write(stdout, delim_fmt)
      write(stdout,blank_fmt)
      write(stdout,'(a22)') 'Model coupling options'
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
      select case (coupled_freq_opt)
      case ('nyear')
         coupled_freq_iopt = -1000
      case ('nmonth')
         coupled_freq_iopt = -1000
      case ('nday')
         if (coupled_freq == 1) then
            lcoupled = .true.
            coupled_freq_iopt = freq_opt_nday
            ncouple_per_day = 1
         else
            coupled_freq_iopt = -1000
         endif
      case ('nhour')
         if (coupled_freq <= 24) then
            lcoupled = .true.
            coupled_freq_iopt = freq_opt_nhour
            ncouple_per_day = 24/coupled_freq
         else
            coupled_freq_iopt = -1000
         endif
      case ('nsecond')
         if (coupled_freq <= seconds_in_day) then
            lcoupled = .true.
            coupled_freq_iopt = freq_opt_nsecond
            ncouple_per_day = seconds_in_day/coupled_freq
         else
            coupled_freq_iopt = -1000
         endif
      case ('nstep')
         if (coupled_freq <= nsteps_per_day) then
            lcoupled = .true.
            coupled_freq_iopt = freq_opt_nstep
            ncouple_per_day = nsteps_per_day/coupled_freq
         else
            coupled_freq_iopt = -1000
         endif
      case ('never')
         lcoupled = .false.
      case default
         coupled_freq_iopt = -2000
      end select
   endif
   call broadcast_scalar(lcoupled, master_task)
   call broadcast_scalar(coupled_freq_iopt, master_task)
   call broadcast_scalar(coupled_freq , master_task)
   if (coupled_freq_iopt == -1000) then
      call exit_POP(sigAbort, &
             'ERROR: Coupling frequency must be at least once per day')
   else if (coupled_freq_iopt == -2000) then
      call exit_POP(sigAbort, &
                    'ERROR: Unknown option for coupling frequency')
   endif
!-----------------------------------------------------------------------
!EOC
 end subroutine init_coupled
!***********************************************************************
!BOP
! !IROUTINE: set_coupled_forcing
! !INTERFACE:
 subroutine set_coupled_forcing(SMF,SMFT,STF,SHF_QSW,FW,TFW,IFRAC)
! !DESCRIPTION:
! This routine calls communicates with NCAR flux coupler to set
! surface forcing data
!
! !REVISION HISTORY:
! same as module
! !INPUT/OUTPUT PARAMETERS:
   real (r8), dimension(nx_block,ny_block,2,max_blocks_clinic), &
      intent(inout) :: &
      SMF, &! surface momentum fluxes (wind stress)
      SMFT ! surface momentum fluxes at T points
   real (r8), dimension(nx_block,ny_block,nt,max_blocks_clinic), &
      intent(inout) :: &
      STF, &! surface tracer fluxes
      TFW ! tracer concentration in water flux
   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      intent(inout) :: &
      SHF_QSW, &! penetrative solar heat flux
      FW, &! fresh water flux
      IFRAC ! fractional ice coverage
!EOP
!BOC
!-----------------------------------------------------------------------
!
! exit if coupling not enabled
!
!-----------------------------------------------------------------------
   if (.not. lcoupled) return
!-----------------------------------------------------------------------
!EOC
 end subroutine set_coupled_forcing
!***********************************************************************
 end module forcing_coupled
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!!***********************************************************************
!
! subroutine set_combined_forcing (STF)
!
!!-----------------------------------------------------------------------
!!
!! this routine combines terms when the "partially-coupled"
!! has been selected
!!
!!-----------------------------------------------------------------------
!
! real (kind=dbl_kind), dimension(nx,ny,nt) ::
! & STF ! surface tracer fluxes
!
! real (kind=dbl_kind), dimension(nx,ny) :: WORK1, WORK2
!
!
!#if coupled
!
! if ( shf_formulation == 'partially-coupled' ) then
! STF(:,:,1) = SHF_COMP(:,:,shf_comp_wrest)
! & + SHF_COMP(:,:,shf_comp_srest)
! & + SHF_COMP(:,:,shf_comp_cpl)
! endif
!
! if ( sfwf_formulation == 'partially-coupled' ) then
! if (sfc_layer_type == sfc_layer_varthick .and.
! & .not. lfw_as_salt_flx) then
! STF(:,:,2) = SFWF_COMP(:,:,sfwf_comp_wrest)
! & + SFWF_COMP(:,:,sfwf_comp_srest)
! FW = SFWF_COMP(:,:,sfwf_comp_cpl)
! & + SFWF_COMP(:,:,sfwf_comp_flxio)
! TFW = TFW_COMP(:,:,:,tfw_comp_cpl)
! & + TFW_COMP(:,:,:,tfw_comp_flxio)
! else
! if ( lms_balance ) then
!
! WORK1 = SFWF_COMP(:,:,sfwf_comp_flxio) / salinity_factor
! WORK2 = SFWF_COMP(:,:,sfwf_comp_cpl)
!
! call ms_balancing ( WORK2, EVAP_F, PREC_F, MELT_F,
! & ROFF_F, SALT_F, 'salt',
! & ICEOCN_F=WORK1 )
!
! STF(:,:,2) = SFWF_COMP(:,:,sfwf_comp_wrest)
! & + SFWF_COMP(:,:,sfwf_comp_srest)
! & + WORK2
! & + SFWF_COMP(:,:,sfwf_comp_flxio) * MASK_SR
! else
!
! STF(:,:,2) = SFWF_COMP(:,:,sfwf_comp_wrest)
! & + SFWF_COMP(:,:,sfwf_comp_srest)
! & + SFWF_COMP(:,:,sfwf_comp_cpl)
! & + SFWF_COMP(:,:,sfwf_comp_flxio)
!
! endif
! endif
! endif
!
!#endif
!!-----------------------------------------------------------------------
!
! end subroutine set_combined_forcing
!
