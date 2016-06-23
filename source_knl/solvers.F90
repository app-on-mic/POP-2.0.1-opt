!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module solvers

!BOP
! !MODULE: solvers
!
! !DESCRIPTION:
!  This module contains routines and operators for solving the elliptic
!  system for surface pressure in the barotropic mode.
!
! !REVISION HISTORY:
!  CVS:$Id: solvers.F90,v 1.14 2002/12/02 13:45:10 pwjones Exp $
!  CVS:$Name: POP_2_0_1 $

! !USES:

   use kinds_mod
   use blocks
   use distribution
   use domain
   use domain_size
   use constants
   use boundary
   use global_reductions
   use gather_scatter
   use broadcast
   use grid
   use io
   use time_management
   use exit_mod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_solvers, &
             elliptic_solver

! !PUBLIC DATA MEMBERS:

   real (r8), dimension (nx_block,ny_block,max_blocks_clinic), &
      public :: & 
      AC,                &! time-independent part of center 9pt weight
      A0_CLINIC           ! time-dependent center weight of 9pt operator
                          !   in baroclinic block distribution

   integer (int_kind), public :: &
      solv_sum_iters      ! accumulated no of iterations (diagnostic)

   real (r8), public ::  &
      rms_residual        ! residual (also a diagnostic)

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  other operator and preconditioner weights for barotropic operator
!
!-----------------------------------------------------------------------

   real (r8), dimension (nx_block,ny_block,max_blocks_tropic) :: & 
      A0,AN,AE,ANE,         &! barotropic (9pt) operator coefficients
      RCALCT_B               ! land mask in barotropic distribution 

   ! with thickened halo
   real (r8), dimension (block_size_x+2*cacg_matlevel,block_size_y+2*cacg_matlevel,max_blocks_tropic) :: & 
      A0L,ANL,AEL,ANEL,      &! barotropic (9pt) operator coefficients
      RCALCT_BL               ! land mask in barotropic distribution 

   real (r8), dimension (:,:,:), allocatable :: & 
      PCC,PCN,PCS,PCE,PCW,  &! preconditioner coefficients
      PCNE,PCSE,PCNW,PCSW

!-----------------------------------------------------------------------
!
!  scalar convergence-related variables
!
!-----------------------------------------------------------------------

   logical (log_kind) :: &
      lprecond            ! true if computed preconditioner to be used

   real (r8) ::          &
      solv_convrg,       &! convergence error criterion
      sor,               &! for jacobi solver
      resid_norm          ! residual normalization

   integer (int_kind), parameter :: &
      solv_pcg = 1,      &! predefined solver types
      solv_cgr = 2,      &
      solv_jac = 3,      &
      solv_capcg = 4      ! add ca-pcg, by Junmin, 2015.11.8

   integer (int_kind) :: &
      solv_itype,        &! integer solver method (1=pcg, 2=cgr, 3=jac)
      solv_max_iters,    &! max number of iterations for solver
      solv_ncheck         ! check convergence every ncheck iterations

   integer (int_kind), parameter :: &
      !cacg_sstep = 8,      &! s-step parameter for ca-pcg Alg, moved to "blocks.F90"
      cacg_smat = 2*cacg_sstep+1,      &! s-step Gram matrix size
      !cacg_matlevel = 4, &!cacg_sstep,     &! number of successvie matrix powers per boundary update, moved to "blocks.F90"
      cacg_block_x = 90,   &! blocking for matrix powers operator
      cacg_block_y = 50     ! blocking for matrix powers operator

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: elliptic_solver
! !INTERFACE:

 subroutine elliptic_solver(PRESS, RHS)

! !DESCRIPTION:
!  Solves the elliptic equation for surface pressure by calling
!  the requested solver routine.  Also redistributes necessary
!  array to the barotropic distribution of blocks for better performance
!  of the solver.
!  The elliptic equation is
!  \begin{equation}
!     AF = B
!  \end{equation}
!  where $F$ is a field (eg surface pressure), $B$ is the right hand side
!  and $A$ is the operator defined as
!  \begin{equation}
!     AF = a \nabla\cdot(H \nabla F)
!  \end{equation}
!  where $a$ is the cell area.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      intent(in) :: &
      RHS                  ! right-hand-side of linear system
                           !  for blocks in baroclinic distribution

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
     intent(inout) :: &
     PRESS              ! on input,  initial guess
                        ! on output, final solution for sfc pressure

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
      P_TROPIC,         &! surface pressure in barotropic distribution
      RHS_TROPIC         ! right hand side  in barotropic distribution

!-----------------------------------------------------------------------
!
!  switch to the barotropic distribution for iterative solvers
!
!-----------------------------------------------------------------------
   ! Added by Junmin Sep 20,2015
   if (solv_itype == solv_capcg) then
      ! For CA-PCG algorithm
      call update_ghost_cells(A0_CLINIC, bndy_tropic, field_loc_center, &
                                           field_type_scalar)
   endif

   call redistribute_blocks(A0,         distrb_tropic, &
                            A0_CLINIC,  distrb_clinic)
!      if (my_task == nproc_cpu_pn*nnode) then
!         write(6,*) 'elliptic_solver:A0,A0_CLINIC:',my_task,sum(A0(nghost+1:nx_block-nghost,nghost+1:ny_block-nghost,1:nx_bkfactor*ny_bkfactor)),sum(A0_CLINIC(nghost+1:nx_block-nghost,nghost+1:ny_block-nghost,1:nx_bkfactor*ny_bkfactor))
!         write(6,*) 'elliptic_solver:A0,A0_CLINIC:',my_task,sum(A0(1:nghost+90,1:nghost,1))+sum(A0(1:nghost,nghost+1:nghost+25,1)),sum(A0_CLINIC(1:nghost+90,1:nghost,1))+sum(A0_CLINIC(1:nghost,nghost+1:nghost+25,1))
!      endif
   call redistribute_blocks(P_TROPIC,   distrb_tropic, &
                            PRESS,      distrb_clinic)
   call redistribute_blocks(RHS_TROPIC, distrb_tropic, &
                            RHS,        distrb_clinic)

!-----------------------------------------------------------------------
!
!  call proper routine based on user choice of solver
!
!-----------------------------------------------------------------------

   if (my_task < distrb_tropic%nprocs) then
      select case(solv_itype)
      case (solv_pcg)
         call pcg(P_TROPIC,RHS_TROPIC)      ! precond conjg grad solver
      case (solv_capcg)
         call capcg(P_TROPIC,RHS_TROPIC)  ! communication-avoiding precond conjg grad solver
      case (solv_cgr)
         call cgr(P_TROPIC,RHS_TROPIC)      ! conjugate residual solver
      case (solv_jac)
         call jacobi(P_TROPIC,RHS_TROPIC)   ! simple jacobi solver
      end select
   endif

!-----------------------------------------------------------------------
!
!  switch solution back to the baroclinic distribution
!
!-----------------------------------------------------------------------

   call redistribute_blocks(PRESS,    distrb_clinic, &
                            P_TROPIC, distrb_tropic)

!-----------------------------------------------------------------------
!EOC

 end subroutine elliptic_solver

!***********************************************************************
!BOP
! !IROUTINE: init_solvers
! !INTERFACE:

 subroutine init_solvers

! !DESCRIPTION:
!  This routine initializes choice of solver, calculates the 
!  coefficients of the 9-point stencils for the barotropic operator and
!  reads in a preconditioner if requested.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables:
!
!       {X,Y}{NE,SE,NW,SW} = contribution to {ne,se,nw,sw} coefficients 
!         from {x,y} components of divergence
!       HU = depth at U points
!
!-----------------------------------------------------------------------

   real (r8) ::         &
      xne,xse,xnw,xsw,  &! contribution to coefficients from x,y
      yne,yse,ynw,ysw,  &!   components of divergence
      ase,anw,asw

   character (char_len) :: &
      solv_type,           &! user choice of solver method
      precond_file          ! file containing preconditioner

   namelist /solver_nml/ solv_convrg, solv_max_iters, solv_ncheck, &
                         lprecond, solv_type, precond_file

   integer (int_kind) :: &
      i,j,n,             &! dummy counter
      iblock,            &! block counter
      ncheck,            &! scalar for checking PC/mask compatibility
      nu,                &! I/O unit number and status flag
      nml_error           ! namelist i/o error flag

   integer (int_kind),dimension(:), allocatable :: &
      icheck              ! check for PC/mask compatibility in block

   real (r8), dimension(:,:,:), allocatable :: &
      WORK0,WORKC,WORKN,  &! temp space for computing operator and
      WORKS,WORKE,WORKW,  &! preconditioner coefficients
      WORKNE,WORKNW,WORKSE,WORKSW, &
      RCALC_TMP

   logical (log_kind) :: &
      mlandne, mlandnw, mlandse, mlandsw ! land mask at nbr points

   type (block) ::      &
      this_block         ! block information for current block

!-----------------------------------------------------------------------
!
!  read solver choice and solver constants from namelist input
!  (namelist input file opened in initial.F)
!
!-----------------------------------------------------------------------

   solv_convrg    = eps
   solv_max_iters = 1000
   solv_ncheck    = 10
   lprecond       = .false.
   solv_type      = 'pcg'
   precond_file   = 'empty'

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=solver_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR reading solver_nml')
   endif

   if (my_task == master_task) then
      write(stdout,delim_fmt)
      write(stdout,blank_fmt)
      write(stdout,'(a35)') ' Solver options (barotropic solver)'
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)

      write(stdout,'(a13)') ' Solver type:'
      ! Add capcg by Junmin, 2015.11.08
      select case(solv_type)
      case('pcg')
         write(stdout,'(a35)') '  Preconditioned Conjugate Gradient'
         solv_itype = solv_pcg
      case('capcg')
         write(stdout,'(a28)') '  Communication-avoiding PCG'
         solv_itype = solv_capcg
      case('cgr')
         write(stdout,'(a29)') '  Conjugate Gradient Residual'
         solv_itype = solv_cgr
      case('jac')
         write(stdout,'(a8)')  '  Jacobi'
         solv_itype = solv_jac
      case default
         solv_itype = -1000
      end select

      write(stdout,'(a28,1pe12.5)') ' Solver converged for err < ', &
                                      solv_convrg
      write(stdout,'(a29,i6)') ' Solver maximum iterations = ', &
                                 solv_max_iters
      write(stdout,'(a35,i6,a11)') ' Solver convergence checked every ',&
                                     solv_ncheck, ' iterations'

   endif

   call broadcast_scalar(solv_convrg,    master_task)
   call broadcast_scalar(solv_max_iters, master_task)
   call broadcast_scalar(solv_ncheck,    master_task)
   call broadcast_scalar(lprecond,       master_task)
   call broadcast_scalar(solv_itype,     master_task)
   call broadcast_scalar(precond_file,   master_task)

   if (solv_itype == -1000) then
      call exit_POP(sigAbort, &
                    'unknown solver type: must be pcg, cgr or jacobi')
   endif

!-----------------------------------------------------------------------
!
!  set sor for jacobi solver
!
!-----------------------------------------------------------------------

   sor = p25     ! should be < 1/2

!-----------------------------------------------------------------------
!
!  compute nine point operator coefficients: compute on baroclinic
!  decomposition first where grid info defined and redistribute
!  to barotropic distribution
!  leave A0,AC in baroclinic distribution to facilitate easy
!  time-dependent changes in barotropic routine
!
!-----------------------------------------------------------------------

   allocate(WORK0 (nx_block,ny_block,nblocks_clinic), &
            WORKC (nx_block,ny_block,nblocks_clinic), &
            WORKN (nx_block,ny_block,nblocks_clinic), &
            WORKE (nx_block,ny_block,nblocks_clinic), &
            WORKNE(nx_block,ny_block,nblocks_clinic), &
            RCALC_TMP(nx_block,ny_block,nblocks_clinic))

   do iblock = 1,nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)  

      WORK0    (:,:,iblock) = c0
      WORKC    (:,:,iblock) = c0
      WORKN    (:,:,iblock) = c0
      WORKE    (:,:,iblock) = c0
      WORKNE   (:,:,iblock) = c0
      RCALC_TMP(:,:,iblock) = c0

      do j=2,ny_block
      do i=2,nx_block

         xne = p25*HU(i  ,j  ,iblock)*DXUR(i  ,j  ,iblock)* &
                                      DYU (i  ,j  ,iblock)
         xse = p25*HU(i  ,j-1,iblock)*DXUR(i  ,j-1,iblock)* &
                                      DYU (i  ,j-1,iblock)
         xnw = p25*HU(i-1,j  ,iblock)*DXUR(i-1,j  ,iblock)* &
                                      DYU (i-1,j  ,iblock)
         xsw = p25*HU(i-1,j-1,iblock)*DXUR(i-1,j-1,iblock)* &
                                      DYU (i-1,j-1,iblock)

         yne = p25*HU(i  ,j  ,iblock)*DYUR(i  ,j  ,iblock)* &
                                      DXU (i  ,j  ,iblock)
         yse = p25*HU(i  ,j-1,iblock)*DYUR(i  ,j-1,iblock)* &
                                      DXU (i  ,j-1,iblock)
         ynw = p25*HU(i-1,j  ,iblock)*DYUR(i-1,j  ,iblock)* &
                                      DXU (i-1,j  ,iblock)
         ysw = p25*HU(i-1,j-1,iblock)*DYUR(i-1,j-1,iblock)* &
                                      DXU (i-1,j-1,iblock)

         WORKNE(i,j,iblock) = xne + yne
         ase                = xse + yse
         anw                = xnw + ynw
         asw                = xsw + ysw
 
         WORKE(i,j,iblock)  = xne + xse - yne - yse
         WORKN(i,j,iblock)  = yne + ynw - xne - xnw

         AC(i,j,iblock)  = -(WORKNE(i,j,iblock) + ase + anw + asw)

         WORK0(i,j,iblock)  = TAREA(i,j,iblock)**2
         RCALC_TMP(i,j,iblock) = RCALCT(i,j,iblock)

      end do
      end do
   end do

   A0_CLINIC  = AC

   call redistribute_blocks(AN , distrb_tropic, WORKN,  distrb_clinic)
   call redistribute_blocks(AE , distrb_tropic, WORKE,  distrb_clinic)
   call redistribute_blocks(ANE, distrb_tropic, WORKNE, distrb_clinic)

!-----------------------------------------------------------------------
!
!  calculate normalization constant (darea,darea) for rms_residual
!  in cgr routine.
!
!-----------------------------------------------------------------------

   resid_norm = c1/global_sum(WORK0, distrb_clinic, &
                              field_loc_center, RCALC_TMP)
   solv_convrg = solv_convrg**2/resid_norm

   call redistribute_blocks(RCALCT_B,   distrb_tropic, &
                            RCALC_TMP,  distrb_clinic)

   deallocate(RCALC_TMP)

   ! Added by Junmin Sep 20,2015
   ! For CA-PCG algorithm
   if (solv_itype == solv_capcg) then
   call update_ghost_cells(AN, bndy_tropic, field_loc_center, &
                                           field_type_scalar)
   call update_ghost_cells(AE, bndy_tropic, field_loc_center, &
                                           field_type_scalar)
   call update_ghost_cells(ANE, bndy_tropic, field_loc_center, &
                                           field_type_scalar)
   call update_ghost_cells(RCALCT_B, bndy_tropic, field_loc_center, &
                                           field_type_scalar)

   ANL(cacg_matlevel+1:cacg_matlevel+block_size_x,cacg_matlevel+1:cacg_matlevel+block_size_y,:) = AN(nghost+1:nghost+block_size_x,nghost+1:nghost+block_size_y,:)
   call update_matpow_halo(ANL, bndy_capcg, field_loc_center, &
                                           field_type_scalar)
   AEL(cacg_matlevel+1:cacg_matlevel+block_size_x,cacg_matlevel+1:cacg_matlevel+block_size_y,:) = AE(nghost+1:nghost+block_size_x,nghost+1:nghost+block_size_y,:)
   call update_matpow_halo(AEL, bndy_capcg, field_loc_center, &
                                           field_type_scalar)
   ANEL(cacg_matlevel+1:cacg_matlevel+block_size_x,cacg_matlevel+1:cacg_matlevel+block_size_y,:) = ANE(nghost+1:nghost+block_size_x,nghost+1:nghost+block_size_y,:)
   call update_matpow_halo(ANEL, bndy_capcg, field_loc_center, &
                                           field_type_scalar)
   RCALCT_BL(cacg_matlevel+1:cacg_matlevel+block_size_x,cacg_matlevel+1:cacg_matlevel+block_size_y,:) = RCALCT_B(nghost+1:nghost+block_size_x,nghost+1:nghost+block_size_y,:)
   call update_matpow_halo(RCALCT_BL, bndy_capcg, field_loc_center, &
                                           field_type_scalar)
   endif

!-----------------------------------------------------------------------
!
!  setup preconditioner if required
!
!-----------------------------------------------------------------------

   if (lprecond) then

      call exit_POP(sigAbort,'This option not currently supported')
!     if (my_task == master_task) then
!       write(stdout,*) ' Preconditioner read from file: ', &
!                         trim(precond_file)
!     endif
!
!     allocate(WORKS (nx_block,ny_block,nblocks_clinic), &
!              WORKW (nx_block,ny_block,nblocks_clinic), &
!              WORKNW(nx_block,ny_block,nblocks_clinic), &
!              WORKSE(nx_block,ny_block,nblocks_clinic), &
!              WORKSW(nx_block,ny_block,nblocks_clinic))
!
!     allocate(PCC     (nx_block,ny_block,nblocks_tropic), &
!              PCN     (nx_block,ny_block,nblocks_tropic), &
!              PCS     (nx_block,ny_block,nblocks_tropic), &
!              PCE     (nx_block,ny_block,nblocks_tropic), &
!              PCW     (nx_block,ny_block,nblocks_tropic), &
!              PCNE    (nx_block,ny_block,nblocks_tropic), &
!              PCSE    (nx_block,ny_block,nblocks_tropic), &
!              PCNW    (nx_block,ny_block,nblocks_tropic), &
!              PCSW    (nx_block,ny_block,nblocks_tropic))
!
!     allocate(icheck(nblocks_clinic))
!
!-----------------------------------------------------------------------
!
!    read preconditioner and check that it is consistent with
!    KMU field
!
!-----------------------------------------------------------------------
!
!     call open_parallel_file(nu,precond_file,recl_dbl)
!     call read_array(nu,WORKC)
!     call read_array(nu,WORKN)
!     call read_array(nu,WORKS)
!     call read_array(nu,WORKE)
!     call read_array(nu,WORKW)
!     call read_array(nu,WORKNE)
!     call read_array(nu,WORKNW)
!     call read_array(nu,WORKSE)
!     call read_array(nu,WORKSW)
!     call close_parallel_file(nu)
!
!     if (my_task == master_task) then
!       write(stdout,blank_fmt)
!       write(stdout,*) ' file read: ', trim(precond_file)
!     endif
!
!-----------------------------------------------------------------------
!
!    check that PC is consistent with KMU field
!
!-----------------------------------------------------------------------
!
!     do iblock = 1,nblocks_clinic
!
!       this_block = get_block(blocks_clinic(iblock),iblock)  
!
!       icheck(iblock) = 0
!
!       do j=this_block%jb,this_block%je
!       do i=this_block%ib,this_block%ie
!
!         mlandne = .false.
!         mlandnw = .false.
!         mlandse = .false.
!         mlandsw = .false.
!         if (KMU(i  ,j  ,iblock) == 0) mlandne = .true.
!         if (KMU(i-1,j  ,iblock) == 0) mlandnw = .true.
!         if (KMU(i  ,j-1,iblock) == 0) mlandse = .true.
!         if (KMU(i-1,j-1,iblock) == 0) mlandsw = .true.
!
!         if (mlandne .and. WORKNE(i,j,iblock) /= c0)  &
!                           icheck(iblock) = icheck(iblock) + 1
!
!         if (mlandnw .and. WORKNW(i,j,iblock) /= c0)  &
!                           icheck(iblock) = icheck(iblock) + 1
!
!         if (mlandse .and. WORKSE(i,j,iblock) /= c0)  &
!                           icheck(iblock) = icheck(iblock) + 1
!
!         if (mlandsw .and. WORKSW(i,j,iblock) /= c0)  &
!                           icheck(iblock) = icheck(iblock) + 1
!      
!         if (mlandne .and. mlandnw .and. (WORKN(i,j,iblock) /= c0)) &
!                           icheck(iblock) = icheck(iblock) + 1
!         if (mlandne .and. mlandse .and. (WORKE(i,j,iblock) /= c0)) &
!                           icheck(iblock) = icheck(iblock) + 1
!         if (mlandnw .and. mlandsw .and. (WORKW(i,j,iblock) /= c0)) &
!                           icheck(iblock) = icheck(iblock) + 1
!         if (mlandse .and. mlandsw .and. (WORKS(i,j,iblock) /= c0)) &
!                           icheck(iblock) = icheck(iblock) + 1
!         if (mlandne .and. mlandse .and.                            &
!             mlandnw .and. mlandsw .and. (WORKC(i,j,iblock) /= c0)) &
!                           icheck(iblock) = icheck(iblock) + 1
!       end do
!       end do
!     end do
!
!     ncheck = sum(icheck)
!     if (global_sum(ncheck, distrb_clinic) /= 0) then
!       call exit_POP(sigAbort,'PC and KMU are incompatible')
!     endif
!
!     deallocate(icheck)
!
!     call redistribute_blocks(PCC ,distrb_tropic,WORKC ,distrb_clinic)
!     call redistribute_blocks(PCN ,distrb_tropic,WORKN ,distrb_clinic)
!     call redistribute_blocks(PCE ,distrb_tropic,WORKE ,distrb_clinic)
!     call redistribute_blocks(PCS ,distrb_tropic,WORKS ,distrb_clinic)
!     call redistribute_blocks(PCW ,distrb_tropic,WORKW ,distrb_clinic)
!     call redistribute_blocks(PCNE,distrb_tropic,WORKNE,distrb_clinic)
!     call redistribute_blocks(PCNW,distrb_tropic,WORKNW,distrb_clinic)
!     call redistribute_blocks(PCSE,distrb_tropic,WORKSE,distrb_clinic)
!     call redistribute_blocks(PCSW,distrb_tropic,WORKSW,distrb_clinic)
!
!     deallocate(WORKS, WORKW, WORKNW, WORKSE, WORKSW)
!
   else ! no preconditioner

      if (my_task == master_task) then
         write(stdout,'(a18)') ' No preconditioner'
      endif

   endif

   deallocate(WORK0, WORKC, WORKN, WORKE, WORKNE)

!-----------------------------------------------------------------------
!EOC

 end subroutine init_solvers

!***********************************************************************
!BOP
! !IROUTINE: pcg
! !INTERFACE:

 subroutine pcg(X,B)

! !DESCRIPTION:
!  This routine uses a preconditioned conjugate-gradient solver to
!  solve the equation $Ax=b$.  Both the operator $A$ and preconditioner
!  are nine-point stencils. If no preconditioner has been supplied,
!  a diagonal preconditioner is applied.  Convergence is checked
!  every {\em ncheck} steps.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: &
      B                         ! right hand side of linear system

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(inout) :: &
      X                  ! on input,  an initial guess for the solution
                         ! on output, solution of the linear system

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      m,                 &! local iteration counter
      ierr,              &! local MPI error flag
      iblock              ! local block     counter
   integer (int_kind) :: &
      i,j                 ! local iteration counter

   real (r8) ::          &
      eta0,eta1,rr        ! scalar inner product results

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
      A0R,               &! reciprocal of A0
      R,                 &! residual (b-Ax)
      S,                 &! conjugate direction vector
      Q,WORK0,WORK1       ! various cg intermediate results

   character (char_len) :: & 
      noconvrg           ! error message for no convergence

   type (block) ::      &
      this_block         ! block information for current block

!-----------------------------------------------------------------------
!
!  compute initial residual and initialize S
!
!-----------------------------------------------------------------------
   where (A0 /= c0)
      A0R = c1/A0
   elsewhere
      A0R = c0
   endwhere

   !$OMP PARALLEL DO PRIVATE(iblock,this_block,i,j) 

   do iblock=1,nblocks_tropic
      this_block = get_block(blocks_tropic(iblock),iblock)  

!      call btrop_operator(R,X,this_block,iblock)
      R(:,1:nghost,iblock) = B(:,1:nghost,iblock)
   do j=this_block%jb,this_block%je
      R(1:nghost,j,iblock) = B(1:nghost,j,iblock)
   do i=this_block%ib,this_block%ie
      R(i,j,iblock) = A0 (i  ,j  ,iblock)*X(i  ,j  ,iblock) + &
                    AN (i  ,j  ,iblock)*X(i  ,j+1,iblock) + &
                    AN (i  ,j-1,iblock)*X(i  ,j-1,iblock) + &
                    AE (i  ,j  ,iblock)*X(i+1,j  ,iblock) + &
                    AE (i-1,j  ,iblock)*X(i-1,j  ,iblock) + &
                    ANE(i  ,j  ,iblock)*X(i+1,j+1,iblock) + &
                    ANE(i  ,j-1,iblock)*X(i+1,j-1,iblock) + &
                    ANE(i-1,j  ,iblock)*X(i-1,j+1,iblock) + &
                    ANE(i-1,j-1,iblock)*X(i-1,j-1,iblock)

      R(i,j,iblock) = B(i,j,iblock) - R(i,j,iblock)
   end do
      R(nx_block-nghost+1:nx_block,j,iblock) = B(nx_block-nghost+1:nx_block,j,iblock)
   end do
      R(:,ny_block-nghost+1:ny_block,iblock) = B(:,ny_block-nghost+1:ny_block,iblock)
      S(:,:,iblock) = c0
   end do ! block loop

   !$OMP END PARALLEL DO

!      if (my_task == nproc_cpu_pn*nnode) then
!         write(6,*) 'pcg-init:R:',my_task,sum(R(nghost+1:nx_block-nghost,nghost+1:ny_block-nghost,1:nx_bkfactor*ny_bkfactor))
!         write(6,*) 'pcg-init:R_halo1:',my_task,sum(R(1:nghost+90,1:nghost,1))+sum(R(1:nghost,nghost+1:nghost+25,1))
!         write(6,*) 'pcg-init:R_h2:',my_task,sum(R(1:nghost+92,nghost+26:nghost+27,1))
!      endif
!-----------------------------------------------------------------------
!
!  initialize fields and scalars
!
!-----------------------------------------------------------------------

   call update_ghost_cells(R, bndy_tropic, field_loc_center, &
                                           field_type_scalar)
   eta0 =c1 
   solv_sum_iters = solv_max_iters
 
!   !ljm tuning
!   do j=1,distrb_tropic%local_block_num
!      i = distrb_tropic%local_block_ids(j)
!      if (i>nblocks_tot) then ! tiny-block
!        this_block = get_block(i,j) 
!        write(6,*) 'pcg-init:R:',this_block%block_id,i-nblocks_tot,sum(R(:,:,j))
!      else
!        write(6,*) 'pcg-init:R:',i,0,sum(R(:,:,j))
!      endif
!   enddo
!   call MPI_BARRIER(MPI_COMM_OCN, ierr)
!-----------------------------------------------------------------------
!
!  iterate
!
!-----------------------------------------------------------------------

   iter_loop: do m = 1, solv_max_iters

   !ljm tuning
   !if (m==3) &
   !   call exit_POP(sigAbort,'DEBUG here.')
!-----------------------------------------------------------------------
!
!     calculate (PC)r 
!     diagonal preconditioner if preconditioner not specified
!
!-----------------------------------------------------------------------

      if (lprecond) then
      !$OMP PARALLEL DO PRIVATE(iblock,this_block) 

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)  

            call preconditioner(WORK1,R,this_block,iblock)
         WORK0(:,:,iblock) = R(:,:,iblock)*WORK1(:,:,iblock)
      end do ! block loop

      !$OMP END PARALLEL DO
      else
      !$OMP PARALLEL DO PRIVATE(iblock,j) 

      do iblock=1,nblocks_tropic
!        !DIR$ UNROLL_AND_JAM  ! cause step 1 not converge
        do j=1,ny_block
!            where (A0(:,j,iblock) /= c0)
!               WORK1(:,j,iblock) = R(:,j,iblock)*A0R(:,j,iblock)
!            elsewhere
!               WORK1(:,j,iblock) = c0
!            endwhere
         WORK1(:,j,iblock) = R(:,j,iblock)*A0R(:,j,iblock)

         WORK0(:,j,iblock) = R(:,j,iblock)*WORK1(:,j,iblock)
        end do ! j
      end do ! block loop

      !$OMP END PARALLEL DO
      endif

!-----------------------------------------------------------------------
!
!     update conjugate direction vector s
!
!-----------------------------------------------------------------------

!      if (my_task == nproc_cpu_pn*nnode .and. &
!          m <= 3 ) then
!         write(6,*) 'pcg_iter:R,WORK1,WORK0:',my_task,m,sum(R(nghost+1:nx_block-nghost,nghost+1:ny_block-nghost,1:nx_bkfactor*ny_bkfactor)),sum(WORK1(nghost+1:nx_block-nghost,nghost+1:ny_block-nghost,1:nx_bkfactor*ny_bkfactor)),sum(WORK0(nghost+1:nx_block-nghost,nghost+1:ny_block-nghost,1:nx_bkfactor*ny_bkfactor))
!         write(6,*) 'pcg_iter:R_h2,WORK1_h2,WORK0_h2:',my_task,m,&
!         sum(R(1:nghost+92,nghost+26:nghost+27,1)),&
!         sum(WORK1(1:nghost+92,nghost+26:nghost+27,1)),&
!         sum(WORK0(1:nghost+92,nghost+26:nghost+27,1))
!      endif
      if (lprecond) &
         call update_ghost_cells(WORK1,bndy_tropic, field_loc_center,&
                                                    field_type_scalar)
      !*** (r,(PC)r)
      eta1 = global_sum(WORK0, distrb_tropic, field_loc_center, RCALCT_B)

!      if (my_task == nproc_cpu_pn*nnode .and. &
!          m <= 3 ) then
!         write(6,*) 'pcg_iter:S,A0,AN,AE,ANE:',my_task,m,sum(S(nghost+1:nx_block-nghost,nghost+1:ny_block-nghost,1:nx_bkfactor*ny_bkfactor)),&
!        sum(A0(nghost+0:nx_block-nghost+1,nghost+0:ny_block-nghost+1,1:nx_bkfactor*ny_bkfactor)),&
!        sum(AN(nghost+0:nx_block-nghost+1,nghost+0:ny_block-nghost+1,1:nx_bkfactor*ny_bkfactor)),&
!        sum(AE(nghost+0:nx_block-nghost+1,nghost+0:ny_block-nghost+1,1:nx_bkfactor*ny_bkfactor)),&
!        sum(ANE(nghost+0:nx_block-nghost+1,nghost+0:ny_block-nghost+1,1:nx_bkfactor*ny_bkfactor))
!         write(6,*) 'pcg_iter:S_h1,A0_h1,AN_h1,AE_h1,ANE_h1:',my_task,m,sum(S(1:nghost+90,1:nghost,1))+sum(S(1:nghost,nghost+1:nghost+25,1)),&
!        sum(A0(1:nghost+90,1:nghost,1))+sum(A0(1:nghost,nghost+1:nghost+25,1)),&
!        sum(AN(1:nghost+90,1:nghost,1))+sum(AN(1:nghost,nghost+1:nghost+25,1)),&
!        sum(AE(1:nghost+90,1:nghost,1))+sum(AE(1:nghost,nghost+1:nghost+25,1)),&
!        sum(ANE(1:nghost+90,1:nghost,1))+sum(ANE(1:nghost,nghost+1:nghost+25,1))
!      endif
      !$OMP PARALLEL DO PRIVATE(iblock,this_block,i,j) 

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)  

         S(:,:,iblock) = WORK1(:,:,iblock) + S(:,:,iblock)*(eta1/eta0) 

!-----------------------------------------------------------------------
!
!        compute As
!
!-----------------------------------------------------------------------

!         call btrop_operator(Q,S,this_block,iblock)
      Q(:,1:nghost,iblock) = c0
      WORK0(:,1:nghost,iblock) = c0
!   do j=this_block%jb,this_block%je
!      Q(1:nghost,j,iblock) = c0
!      WORK0(1:nghost,j,iblock) = c0
!   do i=this_block%ib,this_block%ie
!      Q(i,j,iblock) = A0 (i  ,j  ,iblock)*S(i  ,j  ,iblock) + &
!                    AN (i  ,j  ,iblock)*S(i  ,j+1,iblock) + &
!                    AN (i  ,j-1,iblock)*S(i  ,j-1,iblock) + &
!                    AE (i  ,j  ,iblock)*S(i+1,j  ,iblock) + &
!                    AE (i-1,j  ,iblock)*S(i-1,j  ,iblock) + &
!                    ANE(i  ,j  ,iblock)*S(i+1,j+1,iblock) + &
!                    ANE(i  ,j-1,iblock)*S(i+1,j-1,iblock) + &
!                    ANE(i-1,j  ,iblock)*S(i-1,j+1,iblock) + &
!                    ANE(i-1,j-1,iblock)*S(i-1,j-1,iblock)
! 
!         WORK0(i,j,iblock) = Q(i,j,iblock)*S(i,j,iblock)
!   end do
!      Q(nx_block-nghost+1:nx_block,j,iblock) = c0
!      WORK0(nx_block-nghost+1:nx_block,j,iblock) = c0
!   end do

!-----------------------------------------------------------------------
!
!     2-dimensional to 1-dimensional transformation wth stencil code
!
!-----------------------------------------------------------------------
      call pcg_iter_1d(Q(:,:,iblock),WORK0(:,:,iblock), &
                       this_block,                      &
                       S(:,:,iblock),A0(:,:,iblock),    &
                       AN(:,:,iblock),AE(:,:,iblock),   &
                       ANE(:,:,iblock))

      Q(:,ny_block-nghost+1:ny_block,iblock) = c0
      WORK0(:,ny_block-nghost+1:ny_block,iblock) = c0
      end do ! block loop

      !$OMP END PARALLEL DO

!      if (my_task == nproc_cpu_pn*nnode .and. &
!          m == 2 ) then
!        if (lmic_proc) then
!         write(6,*) 'pcg_iter:Q_1,Q_2,Q_3,Q_4,Q_5,Q_6,Q_7,Q_8,S:',&
!         my_task,m,&
!         sum(Q(nghost+2:nghost+89,nghost+1:nghost+1,1)),&
!         sum(Q(nghost+2:nghost+89,nghost+25:nghost+25,1)),&
!         sum(A0(nghost+2:nghost+89,nghost+25:nghost+25,1)),&
!         sum(S(nghost+1:nghost+90,nghost+24:nghost+26,1)),&
!         sum(AN(nghost+2:nghost+89,nghost+24:nghost+25,1)),&
!         sum(AE(nghost+1:nghost+89,nghost+25:nghost+25,1)),&
!         sum(ANE(nghost+1:nghost+89,nghost+24:nghost+25,1)),&
!         sum(Q(nghost+1:nghost+1,nghost+2:nghost+24,1)),&
!         sum(Q(nghost+90:nghost+90,nghost+2:nghost+24,1)),&
!         sum(Q(nghost+2:nghost+89,nghost+2:nghost+24,2)),&
!         sum(Q(nghost+2:nghost+89,nghost+2:nghost+24,3)),&
!         sum(Q(nghost+2:nghost+89,nghost+2:nghost+24,4)),&
!         sum(Q(nghost+2:nghost+89,nghost+2:nghost+24,5)),&
!         sum(Q(nghost+2:nghost+89,nghost+2:nghost+24,6)),&
!         sum(Q(nghost+2:nghost+89,nghost+2:nghost+24,7)),&
!         sum(Q(nghost+2:nghost+89,nghost+2:nghost+24,8)),&
!         sum(S(nghost+1:nx_block-nghost,nghost+1:ny_block-nghost,1:nx_bkfactor*ny_bkfactor))
!        else
!         write(6,*) 'pcg_iter:Q_1,Q_2,Q_3,Q_4,Q_5,Q_6,Q_7,Q_8,S:',&
!         my_task,m,&
!         sum(Q(nghost+2:nghost+89,nghost+1:nghost+1,1)),&
!         sum(Q(nghost+2:nghost+89,nghost+25:nghost+25,1)),&
!         sum(A0(nghost+2:nghost+89,nghost+25:nghost+25,1)),&
!         sum(S(nghost+1:nghost+90,nghost+24:nghost+26,1)),&
!         sum(AN(nghost+2:nghost+89,nghost+24:nghost+25,1)),&
!         sum(AE(nghost+1:nghost+89,nghost+25:nghost+25,1)),&
!         sum(ANE(nghost+1:nghost+89,nghost+24:nghost+25,1)),&
!         sum(Q(nghost+1:nghost+1,nghost+2:nghost+24,1)),&
!         sum(Q(nghost+90:nghost+90,nghost+2:nghost+24,1)),&
!         sum(Q(nghost+92:nghost+179,nghost+2:nghost+24,1)),&
!         sum(Q(nghost+2:nghost+89,nghost+27:nghost+49,1)),&
!         sum(Q(nghost+92:nghost+179,nghost+27:nghost+49,1)),&
!         sum(Q(nghost+2:nghost+89,nghost+52:nghost+74,1)),&
!         sum(Q(nghost+92:nghost+179,nghost+52:nghost+74,1)),&
!         sum(Q(nghost+2:nghost+89,nghost+77:nghost+99,1)),&
!         sum(Q(nghost+92:nghost+179,nghost+77:nghost+99,1)),&
!         sum(S(nghost+1:nx_block-nghost,nghost+1:ny_block-nghost,1:nx_bkfactor*ny_bkfactor))
!        endif
!         write(6,*) 'pcg_iter:Q_h1,S_h1:',my_task,m,sum(Q(nx_block-90-nghost+1:nx_block-nghost,ny_block-nghost+1:ny_block,nx_bkfactor*ny_bkfactor))+sum(Q(nx_block-nghost+1:nx_block,ny_block-25-nghost+1:ny_block-nghost,nx_bkfactor*ny_bkfactor)),sum(S(1:nghost+90,1:nghost,1))+sum(S(1:nghost,nghost+1:nghost+25,1))
!      endif
!-----------------------------------------------------------------------
!
!     compute next solution and residual
!
!-----------------------------------------------------------------------

      call update_ghost_cells(Q, bndy_tropic, field_loc_center, &
                                              field_type_scalar)
!   if (my_task == nproc_cpu_pn*nnode) &
!      write(6,*) 'update_ghost_cells:pcg:Q_h2',my_task,&
!         sum(Q(1:nghost+92,nghost+26:nghost+27,1))

!      !ljm tuning
!      do j=1,distrb_tropic%local_block_num
!         i = distrb_tropic%local_block_ids(j)
!         if (i>nblocks_tot) then ! tiny-block
!           this_block = get_block(i,j) 
!           write(6,*) 'pcg-iter:Q:',m,this_block%block_id,i-nblocks_tot,sum(Q(:,:,j))
!         else
!           write(6,*) 'pcg-iter:Q:',m,i,0,sum(Q(:,:,j))
!         endif
!      enddo
!      call MPI_BARRIER(MPI_COMM_OCN, ierr)

      eta0 = eta1
      eta1 = eta0/global_sum(WORK0, distrb_tropic, &
                             field_loc_center, RCALCT_B)
!   if (my_task == nproc_cpu_pn*nnode) &
!      write(6,*) 'after_global_sum:pcg:eta1',my_task,eta1

      if (mod(m,solv_ncheck) == 0) then
      !$OMP PARALLEL DO PRIVATE(iblock,this_block,i,j) 

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)  

         X(:,:,iblock) = X(:,:,iblock) + eta1*S(:,:,iblock)
!         R(:,:,iblock) = R(:,:,iblock) - eta1*Q(:,:,iblock)

!            call btrop_operator(R,X,this_block,iblock)
      R(:,1:nghost,iblock) = B(:,1:nghost,iblock)
      WORK0(:,1:nghost,iblock) = R(:,1:nghost,iblock)**2
   do j=this_block%jb,this_block%je
      R(1:nghost,j,iblock) = B(1:nghost,j,iblock)
      WORK0(1:nghost,j,iblock) = R(1:nghost,j,iblock)**2
   do i=this_block%ib,this_block%ie
      R(i,j,iblock) = A0 (i  ,j  ,iblock)*X(i  ,j  ,iblock) + &
                    AN (i  ,j  ,iblock)*X(i  ,j+1,iblock) + &
                    AN (i  ,j-1,iblock)*X(i  ,j-1,iblock) + &
                    AE (i  ,j  ,iblock)*X(i+1,j  ,iblock) + &
                    AE (i-1,j  ,iblock)*X(i-1,j  ,iblock) + &
                    ANE(i  ,j  ,iblock)*X(i+1,j+1,iblock) + &
                    ANE(i  ,j-1,iblock)*X(i+1,j-1,iblock) + &
                    ANE(i-1,j  ,iblock)*X(i-1,j+1,iblock) + &
                    ANE(i-1,j-1,iblock)*X(i-1,j-1,iblock)

            R(i,j,iblock) = B(i,j,iblock) - R(i,j,iblock)
            WORK0(i,j,iblock) = R(i,j,iblock)**2
   end do
      R(nx_block-nghost+1:nx_block,j,iblock) = B(nx_block-nghost+1:nx_block,j,iblock)
      WORK0(nx_block-nghost+1:nx_block,j,iblock) = R(nx_block-nghost+1:nx_block,j,iblock)**2
   end do
      R(:,ny_block-nghost+1:ny_block,iblock) = B(:,ny_block-nghost+1:ny_block,iblock)
      WORK0(:,ny_block-nghost+1:ny_block,iblock) = R(:,ny_block-nghost+1:ny_block,iblock)**2

      end do ! block loop

      !$OMP END PARALLEL DO
      else
      !$OMP PARALLEL DO PRIVATE(iblock,j) 

      do iblock=1,nblocks_tropic
        do j=1,ny_block
         X(:,j,iblock) = X(:,j,iblock) + eta1*S(:,j,iblock)
         R(:,j,iblock) = R(:,j,iblock) - eta1*Q(:,j,iblock)
        enddo
      end do ! block loop

      !$OMP END PARALLEL DO
      endif

!-----------------------------------------------------------------------
!
!     test for convergence
!
!-----------------------------------------------------------------------

      if (mod(m,solv_ncheck) == 0) then

         call update_ghost_cells(R, bndy_tropic, field_loc_center,&
                                                 field_type_scalar)
!         !ljm tuning
!         do j=1,distrb_tropic%local_block_num
!            i = distrb_tropic%local_block_ids(j)
!            if (i>nblocks_tot) then ! tiny-block
!              this_block = get_block(i,j) 
!              write(6,*) 'pcg-check:R:',m,this_block%block_id,i-nblocks_tot,sum(R(:,:,j))
!            else
!              write(6,*) 'pcg-check:R:',m,i,0,sum(R(:,:,j))
!            endif
!         enddo
!         call MPI_BARRIER(MPI_COMM_OCN, ierr)

         rr = global_sum(WORK0, distrb_tropic, &
                         field_loc_center, RCALCT_B)   ! (r,r)

         if (rr < solv_convrg) then
            ! ljm tuning
            if (my_task == master_task) &
            write(6,*) 'pcg_iter_loop:iter#=',m
            solv_sum_iters = m
            exit iter_loop
         endif

      endif

   enddo iter_loop

   rms_residual = sqrt(rr*resid_norm)

   if (solv_sum_iters == solv_max_iters) then
      if (solv_convrg /= c0) then
         write(noconvrg,'(a45,i11)') & 
           'Barotropic solver not converged at time step ', nsteps_total
         call exit_POP(sigAbort,noconvrg)
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine pcg

!DIR$ ATTRIBUTES FORCEINLINE :: pcg_iter_1d
 subroutine pcg_iter_1d(Q,WORK0,this_block, &
                        S,A0,AN,AE,ANE)

! !DESCRIPTION:
!
! !INPUT PARAMETERS:

   type (block), intent(in) :: &
      this_block             ! block info for this sub block

   real (r8), dimension(1:nx_block*ny_block), intent(in) :: &
     S,A0,AN,AE,ANE

! !OUTPUT PARAMETERS:

   real (r8), dimension(1:nx_block*ny_block), intent(out) :: &
     Q,WORK0

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j                 ! dummy indices

   do j=this_block%jb,this_block%je
      Q(1+(j-1)*nx_block:nghost+(j-1)*nx_block) = c0
      WORK0(1+(j-1)*nx_block:nghost+(j-1)*nx_block) = c0
      do i=this_block%ib,this_block%ie
         Q((i)+(j-1)*nx_block) = A0 ((i  )+(j  -1)*nx_block)*S((i  )+(j  -1)*nx_block) + &
                    AN ((i  )+(j  -1)*nx_block)*S((i  )+(j+1-1)*nx_block) + &
                    AN ((i  )+(j-1-1)*nx_block)*S((i  )+(j-1-1)*nx_block) + &
                    AE ((i  )+(j  -1)*nx_block)*S((i+1)+(j  -1)*nx_block) + &
                    AE ((i-1)+(j  -1)*nx_block)*S((i-1)+(j  -1)*nx_block) + &
                    ANE((i  )+(j  -1)*nx_block)*S((i+1)+(j+1-1)*nx_block) + &
                    ANE((i  )+(j-1-1)*nx_block)*S((i+1)+(j-1-1)*nx_block) + &
                    ANE((i-1)+(j  -1)*nx_block)*S((i-1)+(j+1-1)*nx_block) + &
                    ANE((i-1)+(j-1-1)*nx_block)*S((i-1)+(j-1-1)*nx_block)
 
         WORK0((i)+(j-1)*nx_block) = Q((i)+(j-1)*nx_block)*S((i)+(j-1)*nx_block)
      end do
      Q(nx_block-nghost+1+(j-1)*nx_block:nx_block+(j-1)*nx_block) = c0
      WORK0(nx_block-nghost+1+(j-1)*nx_block:nx_block+(j-1)*nx_block) = c0
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine pcg_iter_1d

!***********************************************************************
!BOP
! !IROUTINE: cgr
! !INTERFACE:

 subroutine cgr(X,B)

! !DESCRIPTION:
!  This routine uses a conjugate residual solver to
!  solve the equation $Ax=b$.  Both the operator $A$ and preconditioner
!  are nine-point stencils. If no preconditioner has been supplied,
!  a diagonal preconditioner is applied.  Convergence is checked
!  every {\em ncheck} steps.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: &
      B                         ! right hand side of linear system

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(inout) :: &
      X                 ! on input,  an initial guess for the solution
                        ! on output, solution of the linear system

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::   & 
      m,                   &! local iteration counter
      iblock                ! local block     counter

   real (r8) ::            &
      eps0,eps1,eta0,eta1, &! scalar inner product results
      rr,rrold              ! residual norms

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
      R,                   &! residual (b-Ax)
      S,                   &! conjugate direction vector
      Q,                   &! As (operator acting on s)
      AR,                  &! Ar (operator acting on residual)
      WORK                  ! various cg intermediate results

   character (char_len) ::  & 
      noconvrg              ! error message for no convergence

   type (block) ::          &
      this_block            ! block information for current block

!-----------------------------------------------------------------------
!
!  compute initial residual and its norm
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,this_block) 

   do iblock=1,nblocks_tropic
      this_block = get_block(blocks_tropic(iblock),iblock)  

      call btrop_operator(S,X,this_block,iblock)
      R(:,:,iblock) = B(:,:,iblock) - S(:,:,iblock)
      WORK(:,:,iblock) = R(:,:,iblock)*R(:,:,iblock)
      S(:,:,iblock) = c0
      Q(:,:,iblock) = c0
   end do ! block loop

   !$OMP END PARALLEL DO

   call update_ghost_cells(R, bndy_tropic, field_loc_center, &
                                           field_type_scalar)
   ! initial (r,r)
   rr = global_sum(WORK, distrb_tropic, field_loc_center, RCALCT_B)
   rrold = rr
 
!-----------------------------------------------------------------------
!
!  initialize scalars
!
!-----------------------------------------------------------------------

   eps1 = c1 
   solv_sum_iters = solv_max_iters
 
!-----------------------------------------------------------------------
!
!  iterate
!
!-----------------------------------------------------------------------

   iter_loop: do m = 1, solv_max_iters

!-----------------------------------------------------------------------
!
!     calculate Ar
!
!-----------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblock,this_block) 

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)  

         call btrop_operator(AR,R,this_block,iblock)
         WORK(:,:,iblock) = R(:,:,iblock)*AR(:,:,iblock)
      end do ! block loop

      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!     update conjugate direction vector s, and q = As
!
!-----------------------------------------------------------------------

      call update_ghost_cells(AR, bndy_tropic, field_loc_center, &
                                               field_type_scalar)
      ! (r,Ar)
      eps0 = global_sum(WORK, distrb_tropic, field_loc_center, RCALCT_B)
      eta1 = eps0/eps1

      !$OMP PARALLEL DO PRIVATE(iblock,this_block) 

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)  

         S(:,:,iblock) =  R(:,:,iblock) +  eta1*S(:,:,iblock) 
         Q(:,:,iblock) = AR(:,:,iblock) +  eta1*Q(:,:,iblock) 
         WORK(:,:,iblock) = Q(:,:,iblock)*Q(:,:,iblock)
      end do ! block loop

      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!     compute next solution and residual, update (r,r).
!     every ncheck steps recalculate (r,r) = (b - Ax,b - Ax) and
!     exit if it is not decreasing (due to roundoff error).
!
!-----------------------------------------------------------------------
      ! (As,As)
      eps1 = global_sum(WORK, distrb_tropic, field_loc_center, RCALCT_B)
      eta0 = eps0/eps1

      !$OMP PARALLEL DO PRIVATE(iblock,this_block) 

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)  

         X(:,:,iblock) = X(:,:,iblock) + eta0*S(:,:,iblock)

         if (mod(m,solv_ncheck) == 0) then
            call btrop_operator(R,X,this_block,iblock)
            R(:,:,iblock) = B(:,:,iblock) - R(:,:,iblock)
            WORK(:,:,iblock) = R(:,:,iblock)*R(:,:,iblock)
         else
            R(:,:,iblock) = R(:,:,iblock) - eta0*Q(:,:,iblock)
         endif
      end do ! block loop

      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!     test for convergence
!     every ncheck steps the residual norm is recalculated as
!     r = b - Ax to avoid roundoff error in the accumulated
!     residual (r,r) = old (r,r) - eta0*q.  if the recalculated
!     (r,r) is not less than its previously calculated value, 
!     then the solution is not converging due to machine roundoff
!     error, and the routine is exited.
!
!-----------------------------------------------------------------------

      if (mod(m,solv_ncheck) == 0) then
         call update_ghost_cells(R, bndy_tropic, field_loc_center, &
                                                 field_type_scalar)

         rr = global_sum(WORK, distrb_tropic, field_loc_center, RCALCT_B)
         if (rr > rrold) then
            solv_sum_iters = m
            exit iter_loop
         endif
         rrold = rr
      else
         rr = rr - eta0**2*eps1
      endif

      eps1 = eps0  ! update for next pass

      if (rr < solv_convrg) then
         solv_sum_iters = m
         exit iter_loop
      endif

   enddo iter_loop

   rms_residual = sqrt(abs(rr)*resid_norm)

   if (solv_sum_iters == solv_max_iters) then
      if (solv_convrg /= c0) then
        write(noconvrg,'(a45,i11)') & 
           'Barotropic solver not converged at time step ', nsteps_total
        call exit_POP(sigAbort,noconvrg)
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine cgr

!***********************************************************************
!BOP
! !IROUTINE: jacobi
! !INTERFACE:

 subroutine jacobi(X,B)

! !DESCRIPTION:
!  This routine uses a simple Richardson-Jacobi solver to
!  solve the equation $Ax=b$.  Both the operator $A$ and preconditioner
!  are nine-point stencils. If no preconditioner has been supplied,
!  a diagonal preconditioner is applied.  Convergence is checked
!  every {\em ncheck} steps.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: &
      B                         ! right hand side of linear system

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(inout) :: &
      X                 ! on input,  an initial guess for the solution
                        ! on output, solution of the linear system

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  & 
      m,                  &! local iteration counter
      iblock               ! local block     counter

   real (r8) ::           &
      rr                   ! scalar inner product results

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
      R,                  &! residual (b-Ax)
      WORK0,WORK1          ! various intermediate results

   character (char_len) :: & 
      noconvrg             ! error message for no convergence

   type (block) ::         &
      this_block           ! block information for current block

!-----------------------------------------------------------------------
!
!  iterate
!
!-----------------------------------------------------------------------

   solv_sum_iters = solv_max_iters

   iter_loop: do m = 1, solv_max_iters

!-----------------------------------------------------------------------
!
!     calculate residual r = b - Ax
!
!-----------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblock,this_block) 

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)  

         call btrop_operator(WORK0,X,this_block,iblock)
         R(:,:,iblock) = B(:,:,iblock) - WORK0(:,:,iblock)
      end do ! block loop

      !$OMP END PARALLEL DO

      call update_ghost_cells(R, bndy_tropic, field_loc_center, &
                                              field_type_scalar)

!-----------------------------------------------------------------------
!
!     calculate (PC)r 
!
!-----------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblock,this_block) 

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)  

         if (lprecond) then
            call preconditioner(WORK1,R,this_block,iblock)
         else
            where (A0(:,:,iblock) /= c0)
               WORK1(:,:,iblock) = R(:,:,iblock)/A0(:,:,iblock)
            elsewhere
               WORK1(:,:,iblock) = c0
            endwhere
         endif

!-----------------------------------------------------------------------
!
!        compute next solution
!
!-----------------------------------------------------------------------

         X(:,:,iblock) = X(:,:,iblock) + sor*WORK1(:,:,iblock)
      end do ! block loop

      !$OMP END PARALLEL DO

      if (lprecond) &
         call update_ghost_cells(X, bndy_tropic, field_loc_center, &
                                                 field_type_scalar)

!-----------------------------------------------------------------------
!
!     test for convergence
!
!-----------------------------------------------------------------------

      if (mod(m,solv_ncheck) == 0) then

         WORK0 = R*R
         ! (r,r)
         rr = global_sum(WORK0, distrb_tropic, field_loc_center, RCALCT_B)
          
         if (rr < solv_convrg) then
            solv_sum_iters = m
            exit iter_loop
         endif

      endif

   enddo iter_loop

   rms_residual = sqrt(rr*resid_norm)

   if (solv_sum_iters == solv_max_iters) then
      if (solv_convrg /= c0) then
         write(noconvrg,'(a45,i11)') & 
           'Barotropic solver not converged at time step ', nsteps_total
         call exit_POP(sigAbort,noconvrg)
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine jacobi

!***********************************************************************
!BOP
! !IROUTINE: preconditioner
! !INTERFACE:

 subroutine preconditioner(PX,X,this_block,bid)

! !DESCRIPTION:
!  This function applies a precomputed preconditioner as a nine-point
!  stencil operator.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: & 
      X                     ! array to be operated on 

   type (block), intent(in) :: &
      this_block             ! block info for this block

   integer (int_kind), intent(in) :: &
      bid                    ! local block address for this block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(out) :: &
      PX                  ! nine point operator result

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j                ! dummy counters

!-----------------------------------------------------------------------

   PX(:,:,bid) = c0

   do j=this_block%jb,this_block%je
   do i=this_block%ib,this_block%ie
      PX(i,j,bid) = PCNE(i,j,bid)*X(i+1,j+1,bid) + &
                    PCNW(i,j,bid)*X(i-1,j+1,bid) + &
                    PCSE(i,j,bid)*X(i+1,j-1,bid) + &
                    PCSW(i,j,bid)*X(i-1,j-1,bid) + &
                    PCN (i,j,bid)*X(i  ,j+1,bid) + &
                    PCS (i,j,bid)*X(i  ,j-1,bid) + &
                    PCE (i,j,bid)*X(i+1,j  ,bid) + &
                    PCW (i,j,bid)*X(i-1,j  ,bid) + &
                    PCC (i,j,bid)*X(i  ,j  ,bid)
   end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine preconditioner

!***********************************************************************
!BOP
! !IROUTINE: btrop_operator
! !INTERFACE:

! !DIR$ ATTRIBUTES FORCEINLINE :: btrop_operator
 subroutine btrop_operator(AX,X,this_block,bid)

! !DESCRIPTION:
!  This routine applies the nine-point stencil operator for the
!  barotropic solver.  It takes advantage of some 9pt weights being 
!  shifted versions of others.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: & 
      X                  ! array to be operated on 

   type (block), intent(in) :: &
      this_block             ! block info for this block

   integer (int_kind), intent(in) :: &
      bid                    ! local block address for this block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(out) :: &
      AX                     ! nine point operator result (Ax)

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j                ! dummy counters

!-----------------------------------------------------------------------

   AX(:,:,bid) = c0

   do j=this_block%jb,this_block%je
   do i=this_block%ib,this_block%ie
      AX(i,j,bid) = A0 (i  ,j  ,bid)*X(i  ,j  ,bid) + &
                    AN (i  ,j  ,bid)*X(i  ,j+1,bid) + &
                    AN (i  ,j-1,bid)*X(i  ,j-1,bid) + &
                    AE (i  ,j  ,bid)*X(i+1,j  ,bid) + &
                    AE (i-1,j  ,bid)*X(i-1,j  ,bid) + &
                    ANE(i  ,j  ,bid)*X(i+1,j+1,bid) + &
                    ANE(i  ,j-1,bid)*X(i+1,j-1,bid) + &
                    ANE(i-1,j  ,bid)*X(i-1,j+1,bid) + &
                    ANE(i-1,j-1,bid)*X(i-1,j-1,bid)
   end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine btrop_operator

!***********************************************************************
!BOP
! !IROUTINE: capcg
! !INTERFACE:

 subroutine capcg(X,B)

! !DESCRIPTION:
!  This routine uses a communication avoiding conjugate-gradient solver to
!  solve the equation $Ax=b$.  Both the operator $A$ and preconditioner
!  are nine-point stencils. If no preconditioner has been supplied,
!  a diagonal preconditioner is applied.  Convergence is checked
!  every {\em ncheck} steps.
!
!  References:
!  [1] M. Hoemmen, Communication-avoiding Krylov subspace methods, PhD thesis,
!  University of California, Berkeley, 2010.  
!  [2] Erin Carson, Communication-Avoiding Krylov Subspace Methods in Theory and
!  Practice, PhD thesis, University of California, Berkeley, 2015.
!  [3] Erin Carson and James Demmel, A Residual Replacement Strategy for Improving
!  the Maximum Attainable Accuracy of s-step Krylov Subspace Methods. 
!  SIAM J. Matrix Anal. Appl., 35(1), pp. 22-43.
!  [4] J. Demmel, M. Hoemmen, M. Mohiyuddin, and K. Yelick, Avoiding communication
!  in computing Krylov subspaces, Tech. Report UCB/EECS-2007-123, EECS Dept., U.C.
!  Berkeley, Oct 2007.
!
!  CA-PCG algorithm: (derived by Junmin on Sep 6,2015)
! 1: r_1=b−Ax_1, z_1=M^(−1) r_1, p_1=z_1
! 2: for k=0,1,..., until convergence do
! 3:  compute V_k=[P_k,Z_k ] from p_(k,1), z_(k,1)
! 4:  compute G_k=V_k^T MV_k
! 5:  assemble B_k such that (3) holds
! 6:  p′_(k,j)=[1,0_2s ]^T, z′_(k,j)=[0_(s+1),1,0_(s−1) ]^T,
!     x′_(k,j)=[0_(2s+1) ]^T
! 7:  for j=1:s do
! 8:    alpha_(k,j)=(z′_(k,j)^T G_k z'_(k,j))/(p′_(k,j)^T G_k B_k p'_(k,j))
! 9:    x′_(k,j+1)=x′_(k,j)+alpha_(k,j) p′_(k,j)
! 10:   z′_(k,j+1)=z′_(k,j)−alpha_(k,j) B_k p′_(k,j)
! 11:   beta_(k,j)=(z′_(k,j+1)^T G_k z'_(k,j+1))/(z′_(k,j)^T G_k z'_(k,j))
! 12:   p′_(k,j+1)=z′_(k,j+1)+beta_(k,j) p′_(k,j)
! 13: end for
! 14: p_(k,j)=V_k  p′_(k,j), z_(k,j)=V_k  z′_(k,j), x_(k,j)=V_k x′_(k,j)+x_(k,1)
! 15:end for
! 
! Assemble B_k from polynomial rho_j(z), rho_j(z) satisfies the three-term recurrence:
! rho_0(z)=1, rho_1(z)=(z−theta_0)rho_0(z)/gamma_0,
! rho_j(z)=((z−theta_(j−1))rho_(j−1)(z)−sigma_(j−2)rho_(j−2)(z)/gamma_(j−1).  (3)

! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: &
      B                         ! right hand side of linear system

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(inout) :: &
      X                  ! on input,  an initial guess for the solution
                         ! on output, solution of the linear system

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      k,s,cb,ce,i,j,     &! local iteration counter
      iblock              ! local block     counter

   real (r8) ::          &
      rr,                &! scalar inner product results
      local_sum,         &! local_sum for inner product
      local_sum2,        &! local_sum for inner product
      work_zgz,          &! intermediate: Z_p^T.G_k.Z_p
      work_zgzn,         &! intermediate: Z_p^T.G_k.Z_p for next iteration(j+1)
      work_pgbp           ! intermediate: P_p^T.G_k.B_k.P_p

   real (r8), dimension(block_size_x+2*cacg_matlevel,block_size_y+2*cacg_matlevel,max_blocks_tropic) :: &
      MASK_A0L           ,&! RCALCT_BL*A0L
      A0RL                 ! reciprocal of A0

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
      R,                 &! residual (b-Ax)
      Z,                 &! M^(-1) R
      P,                 &! conjugate direction vector
      WORK0               ! various cg intermediate results

   real (r8), dimension(block_size_x+2*cacg_matlevel,block_size_y+2*cacg_matlevel,max_blocks_tropic,cacg_smat) :: &
      V_k               ! basis vectors for Krylov subspace of degree s, [P,AP,...,A^sP,R,...,A^(s-1)R]

   real (r8), dimension(cacg_smat,cacg_smat) :: &
      G_k,             &! Gram matrix for s-step: G_k=V_k^T MV_k
      B_k               ! change basis matrix: tridiagonal, covers monomial,Newton,Chebyshev,etc

   real (r8), dimension(6*cacg_sstep+1) :: &
      grwork,          &! intermediate: local sum vector for G_k and/or rr
      grsumN            ! intermediate: global sum vector for G_k and/or rr

   real (r8), dimension(6*cacg_sstep+1,max_blocks_tropic) :: &
      grmpbuf           ! intermediate: local sum per block for G_k and/or rr

   integer (int_kind) :: &
      grcount           ! number of inner product to sum

   real (r8) ::          &
      alpha, beta        ! coefficients for CG vectors update

   real (r8), dimension(cacg_smat) :: &
      P_p,Z_p,X_p,      &! short CG vector p',z',x'
      work_bp,          &! intermediate: B_k.P_p
      work_v0            ! intermediate: temporary vector

   character (char_len) :: & 
      noconvrg           ! error message for no convergence

   type (block) ::      &
      this_block         ! block information for current block

!-----------------------------------------------------------------------

   !if (my_task == master_task) &
   !if (my_task == nproc_cpu_pn*nnode) &
      !write(6,"(4(A,I3))") 'capcg_param: sstep=',cacg_sstep,', matlevel=',cacg_matlevel,', block_x=', cacg_block_x,', block_y=', cacg_block_y

   A0L(cacg_matlevel+1:cacg_matlevel+block_size_x,cacg_matlevel+1:cacg_matlevel+block_size_y,:) = A0(nghost+1:nghost+block_size_x,nghost+1:nghost+block_size_y,:)
   call update_matpow_halo(A0L, bndy_capcg, field_loc_center, &
                                           field_type_scalar)
   !$OMP PARALLEL DO PRIVATE(iblock,i,j)
   do iblock=1,nblocks_tropic
      do j=1,block_size_y+2*cacg_matlevel,1
      do i=1,block_size_x+2*cacg_matlevel,1
         if (A0L(i,j,iblock) /= c0) then
            A0RL(i,j,iblock) = c1/A0L(i,j,iblock)
            MASK_A0L(i,j,iblock) = RCALCT_BL(i,j,iblock)*A0L(i,j,iblock)
         else
            A0RL(i,j,iblock) = c0
            MASK_A0L(i,j,iblock) = c0
         endif
      end do
      end do
   end do ! block loop
   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!  step 1: r_1=b−Ax_1, z_1=M^(−1) r_1, p_1=r_1
!  compute initial residual and conjugate vector
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,this_block)
   do iblock=1,nblocks_tropic
      this_block = get_block(blocks_tropic(iblock),iblock)  

      call btrop_operator(WORK0,X,this_block,iblock)
      R(:,:,iblock) = B(:,:,iblock) - WORK0(:,:,iblock)
      ! use diagonal preconditioner
      Z(:,:,iblock) = R(:,:,iblock)*A0RL(cacg_matlevel-nghost+1:cacg_matlevel+block_size_x+nghost,cacg_matlevel-nghost+1:cacg_matlevel+block_size_y+nghost,iblock)
      P(:,:,iblock) = Z(:,:,iblock)
   end do ! block loop
   !$OMP END PARALLEL DO
!   call update_ghost_cells(R, bndy_tropic, field_loc_center, &
!                                           field_type_scalar)

!   call update_ghost_cells(Z, bndy_tropic, field_loc_center, &
!                                           field_type_scalar)
 
!-----------------------------------------------------------------------
!
!  initialize fields and scalars
!
!-----------------------------------------------------------------------

   solv_sum_iters = solv_max_iters
 
!-----------------------------------------------------------------------
!
!  step 2: for k=0, 1 , …, until convergence do
!  iterate
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!  step 3:  compute V_k=[P_k,Z_k ] from p_(k,1), z_(k,1)
!     P_k = [p_(k,1),M^(-1)Ap_(k,1),...,(M^(-1)A)^s p_(k,1)]
!     Z_k = [z_(k,1),M^(-1)Az_(k,1),...,(M^(-1)A)^s z_(k,1)]
!     use diagonal preconditioner if preconditioner not specified
!
!-----------------------------------------------------------------------
      grmpbuf = c0
      call btrop_matpower_operator(V_k,P,Z,A0RL,A0L,MASK_A0L,grmpbuf)

!   if (my_task == master_task) then
!      do j = 1,cacg_smat
!         write(6,"(A,I1,A,E)") 'VS(',j,'):', sum(V_k(cacg_matlevel+1:cacg_matlevel+block_size_x,cacg_matlevel+1:cacg_matlevel+block_size_y,:,j))
!      end do
!   endif

   iter_loop: do k = 1, solv_max_iters/cacg_sstep

!-----------------------------------------------------------------------
!
!  step 4:  compute G_k=V_k^T MV_k
!    [ P^T MP, P^T MR ]
!    [ R^T MP, R^T MR ]
!    For monomial bases and symmetric A and M, 6*sstep dot products are
!    calculated. 
!-----------------------------------------------------------------------
      grwork = c0
!    [ P^T MP ]
!     r_(k,j) = M z_(k,j)
!     rr = r_(k,j)^T r_(k,j)
      grcount = 0
      ce = 1
      cb = 1
            ! G(cb,ce) = V_k(:,cb)*A0R(:)*V_k(:,ce)
            local_sum = c0
            local_sum2 = c0
            !$OMP PARALLEL DO PRIVATE(iblock,this_block,i,j) &
            !$OMP             REDUCTION(+:local_sum,local_sum2)
            do iblock=1,nblocks_tropic
               this_block = get_block(blocks_tropic(iblock),iblock)  
               do j=this_block%jb+cacg_matlevel-nghost,this_block%je+cacg_matlevel-nghost
               !DIR$ SIMD REDUCTION(+:local_sum,local_sum2)
               do i=this_block%ib+cacg_matlevel-nghost,this_block%ie+cacg_matlevel-nghost
! Empty loop
!                  local_sum = local_sum + V_k(i,j,iblock,cb)*MASK_A0L(i,j,iblock)*V_k(i,j,iblock,ce)
!     rr = r_(k,j)^T r_(k,j)
!                  local_sum2 = local_sum2 + Z(i-(cacg_matlevel-nghost),j-(cacg_matlevel-nghost),iblock)*Z(i-(cacg_matlevel-nghost),j-(cacg_matlevel-nghost),iblock)*MASK_A0L(i,j,iblock)*A0L(i,j,iblock)
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
            grcount = grcount + 1
!            grwork(grcount) = local_sum

      !do ce=1,cacg_sstep+1,1
      !   do cb=ce,min(ce+1,cacg_sstep+1),1
            ! G(cb,ce) = V_k(:,cb)*A0R(:)*V_k(:,ce)
         !$OMP PARALLEL DO PRIVATE(iblock,this_block,i,j,cb)
         do iblock=1,nblocks_tropic
            this_block = get_block(blocks_tropic(iblock),iblock)  
            do j=this_block%jb+cacg_matlevel-nghost,this_block%je+cacg_matlevel-nghost
               do cb=2,cacg_sstep+1,1
               !DIR$ SIMD 
               do i=this_block%ib+cacg_matlevel-nghost,this_block%ie+cacg_matlevel-nghost
! Empty loop
!    [ P^T MP ] except the first element P^T(M^(-1)A)P
!                  grmpbuf(cb,iblock) = grmpbuf(cb,iblock) + V_k(i,j,iblock,cb)*MASK_A0L(i,j,iblock)*V_k(i,j,iblock,1)
!                  grmpbuf(cb+cacg_sstep,iblock) = grmpbuf(cb+cacg_sstep,iblock) + V_k(i,j,iblock,cb)*MASK_A0L(i,j,iblock)*V_k(i,j,iblock,cacg_sstep+1)
!    [ R^T MP ]
!                  grmpbuf(cb+2*cacg_sstep,iblock) = grmpbuf(cb+2*cacg_sstep,iblock) + V_k(i,j,iblock,cb+cacg_sstep)*MASK_A0L(i,j,iblock)*V_k(i,j,iblock,1)
!                  grmpbuf(cb+3*cacg_sstep,iblock) = grmpbuf(cb+3*cacg_sstep,iblock) + V_k(i,j,iblock,cb+cacg_sstep)*MASK_A0L(i,j,iblock)*V_k(i,j,iblock,cacg_sstep+1)
!    [ R^T MR ]
!                  grmpbuf(cb+4*cacg_sstep,iblock) = grmpbuf(cb+4*cacg_sstep,iblock) + V_k(i,j,iblock,cb+cacg_sstep)*MASK_A0L(i,j,iblock)*V_k(i,j,iblock,cacg_sstep+2)
!                  grmpbuf(cb+5*cacg_sstep,iblock) = grmpbuf(cb+5*cacg_sstep,iblock) + V_k(i,j,iblock,cb+cacg_sstep)*MASK_A0L(i,j,iblock)*V_k(i,j,iblock,2*cacg_sstep+1)
               enddo
               enddo
            enddo
         enddo
         !$OMP END PARALLEL DO

         !$OMP PARALLEL DO PRIVATE(iblock,cb)
         !do cb=2,6*cacg_sstep+1,1
         do cb=1,6*cacg_sstep+1,1
         do iblock=1,nblocks_tropic
            grwork(cb) = grwork(cb) + grmpbuf(cb,iblock)
         enddo
         enddo
         !$OMP END PARALLEL DO
     
     !if (my_task == nproc_cpu_pn*nnode) then
        !write(6,"(A,I3,17E)") 'V_k(1):',my_task,sum(V_k(cacg_matlevel+1:cacg_matlevel+block_size_x,cacg_matlevel+1:cacg_matlevel+block_size_y,:,1))
        !write(6,"(A,I3,17E)") 'V_k(2):',my_task,sum(V_k(cacg_matlevel+1:cacg_matlevel+block_size_x,cacg_matlevel+1:cacg_matlevel+block_size_y,:,2))
        !write(6,"(A,I3,17E)") 'RCALCT_BL:',my_task,sum(RCALCT_BL(cacg_matlevel+1:cacg_matlevel+block_size_x,cacg_matlevel+1:cacg_matlevel+block_size_y,:))
        !write(6,"(A,I3,17E)") 'A0L:',my_task,sum(A0L(cacg_matlevel+1:cacg_matlevel+block_size_x,cacg_matlevel+1:cacg_matlevel+block_size_y,:))
     !endif
!    [ R^T MP ]

!    [ R^T MR ]

      grcount = 6*cacg_sstep+1
!     Adjust last column of [ R^T MR ]
!      do ce=5*cacg_sstep+1,6*cacg_sstep,1
!         grwork(ce) = grwork(ce+1)
!      enddo
      grcount = 6*cacg_sstep

!      G_k = global_sum(grwork, distrb_tropic)

!     merged with [P^T MP]_(0,0)
!     test for convergence
!     r_(k,j) = M z_(k,j)
!     rr = r_(k,j)^T r_(k,j)
!            local_sum = c0
!            do iblock=1,nblocks_tropic
!               this_block = get_block(blocks_tropic(iblock),iblock)  
!               do j=this_block%jb,this_block%je
!               !DIR$ SIMD REDUCTION(+:local_sum)
!               do i=this_block%ib,this_block%ie
!                  local_sum = local_sum + Z(i,j,iblock)*Z(i,j,iblock)*RCALCT_B(i,j,iblock)*A0(i,j,iblock)*A0(i,j,iblock)
!               enddo
!               enddo
!            enddo
            grcount = grcount + 1
!            grwork(grcount) = local_sum2  ! computed ahead 

!     Global reduction for G_k and rr
      grsumN = global_sum(grwork, distrb_tropic)

      rr = grsumN(grcount)
      !if (my_task == master_task) &
      !write(6,"(A,I5,A,E)") 'capcg_iter#:',k,', rr=',rr
      if (rr < solv_convrg) then
         if (my_task == master_task) &
         write(6,*) 'cacg_iter_loop:iter#=',k * cacg_sstep
         solv_sum_iters = k * cacg_sstep
         exit iter_loop
      endif

!     G_k
!      G_k = reshape(grsumN(1:cacg_smat*cacg_smat), (/cacg_smat,cacg_smat/))
!    [ P^T MP ]
      do ce=1,cacg_sstep+1,1
         G_k(1:cacg_sstep+1,ce) = grsumN(ce:ce+cacg_sstep)
      enddo
!    [ R^T MP ]
      do cb=cacg_sstep+2,cacg_smat,1
         G_k(cb,1:cacg_sstep+1) = grsumN(cacg_sstep+cb:2*cacg_sstep+cb)
      enddo
!    [ P^T MR ]
      do ce=cacg_sstep+2,cacg_smat,1
         G_k(1:cacg_sstep+1,ce) = grsumN(cacg_sstep+ce:2*cacg_sstep+ce)
      enddo
!    [ R^T MR ]
      do ce=cacg_sstep+2,cacg_smat,1
         G_k(cacg_sstep+2:cacg_smat,ce) = grsumN(3*cacg_sstep+ce:4*cacg_sstep-1+ce)
      enddo
      ! Debug
!      if (my_task == master_task) then
!      !if (my_task == nproc_cpu_pn*nnode) then
!      do j = 1,cacg_smat
!         print "(A,I5,A,I5,A,17E)",'capcg_iter#:',k,', G_k(',j,')', (G_k(i,j), i = 1,cacg_smat)
!      end do
!      endif
!-----------------------------------------------------------------------
!
!  step 5:  assemble B_k such that (3) holds
!-----------------------------------------------------------------------
      ! skipped, use T operator instead
      ! choose monomial basis: T= [D_s,0,D_(s-1),0]^T

!-----------------------------------------------------------------------
!
!  step 6:  p′_(k,j)=[1,0_2s ]^T, z′_(k,j)=[0_(s+1),1,0_(s−1) ]^T,
!       Init coefficient P_p,Z_p,X_p
!-----------------------------------------------------------------------
      P_p = c0
      P_p(1) = c1
      Z_p = c0
      Z_p(cacg_sstep+2) = c1
      X_p = c0

!-----------------------------------------------------------------------
!
!  step 7:  for j=1:s do
!     Inner loop
!-----------------------------------------------------------------------
         ! work_zgz = Z_p^T.G_k.Z_p
         work_v0 = c0
         do i=1,cacg_smat
           do j=1,cacg_smat
            work_v0(j) = work_v0(j) + Z_p(i)*G_k(i,j)
           end do
         end do
         work_zgzn = sum(work_v0*Z_p)
 
      do s=1,cacg_sstep
!-----------------------------------------------------------------------
!
!  step 8:    alpha_(k,j)=(z′_(k,j)^T G_k z'_(k,j))/(p′_(k,j)^T G_k B_k p'_(k,j))
!-----------------------------------------------------------------------
         ! reuse work_zgz from previous(prologue) iteration
         work_zgz = work_zgzn
         ! compute B_k p'_(k,j)
         call btrop_t_operator(work_bp,P_p)
         ! work_pgbp = P_p^T.G_k.B_k.P_p
         do i=1,cacg_smat
            work_v0(i) = sum(P_p*G_k(:,i))
         end do
         work_pgbp = sum(work_v0*work_bp)
         alpha = work_zgz / work_pgbp
!-----------------------------------------------------------------------
!
!  step 9:    x′_(k,j+1)=x′_(k,j)+alpha_(k,j) p′_(k,j)
!  step 10:   z′_(k,j+1)=z′_(k,j)−alpha_(k,j) B_k p′_(k,j)
!-----------------------------------------------------------------------
         X_p = X_p + alpha*P_p
         Z_p = Z_p - alpha*work_bp
!-----------------------------------------------------------------------
!
!  step 11:   beta_(k,j)=(z′_(k,j+1)^T G_k z'_(k,j+1))/(z′_(k,j)^T G_k z'_(k,j))
!-----------------------------------------------------------------------
         ! work_zgzn = Z_p^T.G_k.Z_p for next iteration(j+1)
         do i=1,cacg_smat
            work_v0(i) = sum(Z_p*G_k(:,i))
         end do
         work_zgzn = sum(work_v0*Z_p)
         beta = work_zgzn / work_zgz
!-----------------------------------------------------------------------
!
!  step 12:   p′_(k,j+1)=z′_(k,j+1)+beta_(k,j) p′_(k,j)
!-----------------------------------------------------------------------
         P_p = Z_p + beta*P_p
      end do ! inner loop

!-----------------------------------------------------------------------
!
! 14: p_(k,j)=V_k  p′_(k,j), z_(k,j)=V_k  z′_(k,j), x_(k,j)=V_k x′_(k,j)+x_(k,1)
!     Recover X,R,P vectors from the coordinates results of inner loop
!-----------------------------------------------------------------------
      !$OMP PARALLEL DO PRIVATE(iblock,i,j,s)
      do iblock=1,nblocks_tropic
         Z(:,:,iblock) = c0
         P(:,:,iblock) = c0
      do j=1,ny_block,1
         do s=1,cacg_smat
         do i=1,nx_block,1
            X(i,j,iblock) = X(i,j,iblock) + V_k(i+cacg_matlevel-nghost,j+cacg_matlevel-nghost,iblock,s)*X_p(s)
            Z(i,j,iblock) = Z(i,j,iblock) + V_k(i+cacg_matlevel-nghost,j+cacg_matlevel-nghost,iblock,s)*Z_p(s)
            P(i,j,iblock) = P(i,j,iblock) + V_k(i+cacg_matlevel-nghost,j+cacg_matlevel-nghost,iblock,s)*P_p(s)
         end do
      end do ! 2s+1
      end do
      end do ! block loop
      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!  step 3:  compute V_k=[P_k,Z_k ] from p_(k,1), z_(k,1)
!     P_k = [p_(k,1),M^(-1)Ap_(k,1),...,(M^(-1)A)^s p_(k,1)]
!     Z_k = [z_(k,1),M^(-1)Az_(k,1),...,(M^(-1)A)^s z_(k,1)]
!     use diagonal preconditioner if preconditioner not specified
!
!-----------------------------------------------------------------------
      grmpbuf = c0
      call btrop_matpower_operator(V_k,P,Z,A0RL,A0L,MASK_A0L,grmpbuf)

   enddo iter_loop

   rms_residual = sqrt(rr*resid_norm)

   if (solv_sum_iters == solv_max_iters) then
      if (solv_convrg /= c0) then
         write(noconvrg,'(a45,i11)') & 
           'Barotropic solver not converged at time step ', nsteps_total
         call exit_POP(sigAbort,noconvrg)
      endif
   endif

   call update_ghost_cells(X, bndy_tropic, field_loc_center, &
                                           field_type_scalar)
!-----------------------------------------------------------------------
!EOC

 end subroutine capcg

!***********************************************************************
!BOP
! !IROUTINE: btrop_matpower_operator
! !INTERFACE:

! !DIR$ ATTRIBUTES FORCEINLINE :: btrop_matpower_operator
 subroutine btrop_matpower_operator(VS,P,Z,A0RL,A0L,MASK_A0L,grmpbuf)

! !DESCRIPTION:
!  This routine continuously applies the nine-point stencil operator for the
!  barotropic solver. 
!  It calculate V = [P_k, Z_k]
!     P_k = [p_(k,1),M^(-1)Ap_(k,1),...,(M^(-1)A)^s p_(k,1)]
!     Z_k = [z_(k,1),M^(-1)Az_(k,1),...,(M^(-1)A)^s z_(k,1)]
!  a matrix power by a vector
!
!  Naive version: iterate on interleaved Axpy and boundary exchange;
!  PA2 version as in [1],etc: partial redundant computation and one-time communication
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: & 
      P,                 &! array to be operated on 
      Z                   ! array to be operated on 
 
   real (r8), dimension(block_size_x+2*cacg_matlevel,block_size_y+2*cacg_matlevel,max_blocks_tropic),  &
      intent(in) :: & 
      A0RL              ,&! reciprocal of A0
      A0L               ,&! reciprocal of A0
      MASK_A0L           ! reciprocal of A0
 
! !OUTPUT PARAMETERS:

   real (r8), dimension(block_size_x+2*cacg_matlevel,block_size_y+2*cacg_matlevel,max_blocks_tropic,cacg_smat),  &
      intent(out) :: &
      VS                     ! nine point operator result (V)

   real (r8), dimension(6*cacg_sstep+1,max_blocks_tropic), &
      intent(out) :: &
      grmpbuf           ! intermediate: local sum per block for G_k and/or rr

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   type (block) :: &
      this_block             ! block info for this block

   integer (int_kind) :: &
      bid,                  &! local block address for this block
      s,i,j                  ! dummy counters

   integer (int_kind) :: &
      ss,ii,jj               ! dummy counters

   real (r8) :: &
      ap1,ap2,az1,az2       ,&! intermedium variable
      p1_s,p2_s,p3_s,pz1_s,pz2_s,z_s,r_s ! intermedium variable for local sum

   !integer (int_kind), save :: call_count = 0

   !call_count = call_count + 1
!-----------------------------------------------------------------------

   VS = c0

   VS(cacg_matlevel+1:cacg_matlevel+block_size_x,cacg_matlevel+1:cacg_matlevel+block_size_y,:,1) = P(nghost+1:nghost+block_size_x,nghost+1:nghost+block_size_y,:)
   call update_matpow_halo(VS(:,:,:,1), bndy_capcg, field_loc_center, &
                                           field_type_scalar)
   VS(cacg_matlevel+1:cacg_matlevel+block_size_x,cacg_matlevel+1:cacg_matlevel+block_size_y,:,cacg_sstep+2) = Z(nghost+1:nghost+block_size_x,nghost+1:nghost+block_size_y,:)
   call update_matpow_halo(VS(:,:,:,cacg_sstep+2), bndy_capcg, field_loc_center, &
                                           field_type_scalar)

   do s = 1,cacg_sstep-1,cacg_matlevel

   !$OMP PARALLEL DO PRIVATE(bid,this_block,i,j,ss,ii,jj,ap1,ap2,az1,az2,&
   !$OMP             p1_s,p2_s,p3_s,z_s,pz1_s,pz2_s)
   do bid=1,nblocks_tropic
      this_block = get_block(blocks_tropic(bid),bid)  
   do j=cacg_matlevel+1,cacg_matlevel+block_size_y,cacg_block_y
   do i=cacg_matlevel+1,cacg_matlevel+block_size_x,cacg_block_x
     do ss = s,min(s+cacg_matlevel-1,cacg_sstep-1),1
      p1_s = c0
      p2_s = c0
      p3_s = c0
      pz1_s = c0
      pz2_s = c0
      z_s  = c0
      do jj=j-(cacg_matlevel-1-(ss-s)),min(j+cacg_block_y-1,cacg_matlevel+block_size_y)+(cacg_matlevel-1-(ss-s)),1
      !DIR$ SIMD ASSERT PRIVATE(ap1,ap2,az1,az2) REDUCTION(+:p1_s,p2_s,p3_s,z_s,pz1_s,pz2_s)
      do ii=i-(cacg_matlevel-1-(ss-s)),min(i+cacg_block_x-1,cacg_matlevel+block_size_x)+(cacg_matlevel-1-(ss-s)),1
      !VS(ii,jj,bid,ss+1) = VS(ii  ,jj  ,bid,ss) + &
      ap1 = VS(ii  ,jj  ,bid,ss)
      ap2 = ap1 + &
                   (ANL (ii  ,jj  ,bid)*VS(ii  ,jj+1,bid,ss) + &
                    ANL (ii  ,jj-1,bid)*VS(ii  ,jj-1,bid,ss) + &
                    AEL (ii  ,jj  ,bid)*VS(ii+1,jj  ,bid,ss) + &
                    AEL (ii-1,jj  ,bid)*VS(ii-1,jj  ,bid,ss) + &
                    ANEL(ii  ,jj  ,bid)*VS(ii+1,jj+1,bid,ss) + &
                    ANEL(ii  ,jj-1,bid)*VS(ii+1,jj-1,bid,ss) + &
                    ANEL(ii-1,jj  ,bid)*VS(ii-1,jj+1,bid,ss) + &
                    ANEL(ii-1,jj-1,bid)*VS(ii-1,jj-1,bid,ss) &
                   ) * A0RL(ii,jj,bid)
      !  use diagonal preconditioner
      VS(ii,jj,bid,ss+1) = ap2
      !VS(ii,jj,bid,ss+cacg_sstep+2) = VS(ii  ,jj  ,bid,ss+cacg_sstep+1) + &
      az1 = VS(ii  ,jj  ,bid,ss+cacg_sstep+1)
      az2 = az1 + &
                   (ANL (ii  ,jj  ,bid)*VS(ii  ,jj+1,bid,ss+cacg_sstep+1) + &
                    ANL (ii  ,jj-1,bid)*VS(ii  ,jj-1,bid,ss+cacg_sstep+1) + &
                    AEL (ii  ,jj  ,bid)*VS(ii+1,jj  ,bid,ss+cacg_sstep+1) + &
                    AEL (ii-1,jj  ,bid)*VS(ii-1,jj  ,bid,ss+cacg_sstep+1) + &
                    ANEL(ii  ,jj  ,bid)*VS(ii+1,jj+1,bid,ss+cacg_sstep+1) + &
                    ANEL(ii  ,jj-1,bid)*VS(ii+1,jj-1,bid,ss+cacg_sstep+1) + &
                    ANEL(ii-1,jj  ,bid)*VS(ii-1,jj+1,bid,ss+cacg_sstep+1) + &
                    ANEL(ii-1,jj-1,bid)*VS(ii-1,jj-1,bid,ss+cacg_sstep+1) &
                   ) * A0RL(ii,jj,bid)
      !  use diagonal preconditioner
      VS(ii,jj,bid,ss+cacg_sstep+2) = az2
      ! calc the inner product of V_k(ss) and V_k(ss+1)
      if (jj >= j .and. jj <= min(j+cacg_block_y-1,cacg_matlevel+block_size_y) .and.&
          ii >= i .and. ii <= min(i+cacg_block_x-1,cacg_matlevel+block_size_x)) then
!    [ P^T MP ] except the first element P^T(M^(-1)A)P
      p1_s = p1_s + ap1*ap1*MASK_A0L(ii,jj,bid)
      p2_s = p2_s + ap2*ap1*MASK_A0L(ii,jj,bid)
!    [ R^T MP ]
      pz1_s = pz1_s + ap1*az1*MASK_A0L(ii,jj,bid)
      pz2_s = pz2_s + ap2*az1*MASK_A0L(ii,jj,bid)
!    [ R^T MR ]
      p3_s = p3_s + az1*az1*MASK_A0L(ii,jj,bid)
      z_s = z_s + az2*az1*MASK_A0L(ii,jj,bid)
      endif
      end do ! ii
      end do ! jj
!    [ P^T MP ] except the first element P^T(M^(-1)A)P
      grmpbuf(ss+ss-1,bid) = grmpbuf(ss+ss-1,bid) + p1_s
      grmpbuf(ss+ss,bid) = grmpbuf(ss+ss,bid) + p2_s
!    [ R^T MP ]
      grmpbuf(ss+ss+2*cacg_sstep,bid) = grmpbuf(ss+ss+2*cacg_sstep,bid) + pz1_s
      grmpbuf(ss+ss+1+2*cacg_sstep,bid) = grmpbuf(ss+ss+1+2*cacg_sstep,bid) + pz2_s
!    [ R^T MR ]
      grmpbuf(ss+ss+4*cacg_sstep,bid) = grmpbuf(ss+ss+4*cacg_sstep,bid) + p3_s
      grmpbuf(ss+ss+1+4*cacg_sstep,bid) = grmpbuf(ss+ss+1+4*cacg_sstep,bid) + z_s
     end do ! do ss
   end do ! i
   end do ! j
   end do ! block loop
   !$OMP END PARALLEL DO
   !if (my_task == nproc_cpu_pn*nnode .and. call_count == 1 .and. s==1) &
   !write(6,"(A,I3,17E)") 'VS(2):',my_task,sum(VS(cacg_matlevel+1:cacg_matlevel+block_size_x,cacg_matlevel+1:cacg_matlevel+block_size_y,:,2))
   !do ss = s,min(s+cacg_matlevel-1,cacg_sstep-1),1
   ss = min(s+cacg_matlevel-1,cacg_sstep-1)
   call update_matpow_halo(VS(:,:,:,ss+1), bndy_capcg, field_loc_center, &
                                           field_type_scalar)
   call update_matpow_halo(VS(:,:,:,ss+cacg_sstep+2), bndy_capcg, field_loc_center, &
                                           field_type_scalar)

   end do ! do s matrix powers

   !$OMP PARALLEL DO PRIVATE(ap1,ap2,az1,az2,p1_s,p2_s,p3_s,z_s,r_s,pz1_s,pz2_s,i,j)
   do bid=1,nblocks_tropic
      p1_s = c0
      p2_s = c0
      p3_s = c0
      pz1_s = c0
      pz2_s = c0
      z_s = c0
      r_s = c0
   do j=cacg_matlevel+1,cacg_matlevel+block_size_y,1
   !DIR$ SIMD PRIVATE(ap1,ap2) REDUCTION(+:p1_s,p2_s,p3_s,z_s,r_s)
   do i=cacg_matlevel+1,cacg_matlevel+block_size_x,1
      !VS(i,j,bid,cacg_sstep+1) = VS(i  ,j  ,bid,cacg_sstep) + &
      ap1 = VS(i  ,j  ,bid,cacg_sstep)
      ap2 = ap1 + &
                   (ANL (i  ,j  ,bid)*VS(i  ,j+1,bid,cacg_sstep) + &
                    ANL (i  ,j-1,bid)*VS(i  ,j-1,bid,cacg_sstep) + &
                    AEL (i  ,j  ,bid)*VS(i+1,j  ,bid,cacg_sstep) + &
                    AEL (i-1,j  ,bid)*VS(i-1,j  ,bid,cacg_sstep) + &
                    ANEL(i  ,j  ,bid)*VS(i+1,j+1,bid,cacg_sstep) + &
                    ANEL(i  ,j-1,bid)*VS(i+1,j-1,bid,cacg_sstep) + &
                    ANEL(i-1,j  ,bid)*VS(i-1,j+1,bid,cacg_sstep) + &
                    ANEL(i-1,j-1,bid)*VS(i-1,j-1,bid,cacg_sstep) &
                   ) * A0RL(i,j,bid)
      !  use diagonal preconditioner 
      VS(i,j,bid,cacg_sstep+1) = ap2
!    [ P^T MP ] except the first element P^T(M^(-1)A)P
      ! calc the inner product of V_k(ss) and V_k(ss+1)
      p1_s = p1_s + ap1*ap1*MASK_A0L(i,j,bid)
      p2_s = p2_s + ap2*ap1*MASK_A0L(i,j,bid)
      p3_s = p3_s + ap2*ap2*MASK_A0L(i,j,bid)
!    [ R^T MP ]
      pz1_s = pz1_s + ap1*VS(i,j,bid,2*cacg_sstep+1)*MASK_A0L(i,j,bid)
      pz2_s = pz2_s + ap2*VS(i,j,bid,2*cacg_sstep+1)*MASK_A0L(i,j,bid)
!    [ R^T MR ]
      z_s = z_s + VS(i,j,bid,2*cacg_sstep+1)*VS(i,j,bid,2*cacg_sstep+1)*MASK_A0L(i,j,bid)
!     rr = r_(k,j)^T r_(k,j)
      r_s = r_s + VS(i,j,bid,cacg_sstep+2)*VS(i,j,bid,cacg_sstep+2)*MASK_A0L(i,j,bid)*A0L(i,j,bid)
   end do
   end do
!    [ P^T MP ] 
      grmpbuf(2*cacg_sstep-1,bid) = p1_s
      grmpbuf(2*cacg_sstep,bid) = p2_s
      grmpbuf(2*cacg_sstep+1,bid) = p3_s
!    [ R^T MP ]
      grmpbuf(4*cacg_sstep,bid) = grmpbuf(4*cacg_sstep,bid) + pz1_s
      grmpbuf(4*cacg_sstep+1,bid) = grmpbuf(4*cacg_sstep+1,bid) + pz2_s
!    [ R^T MR ]
      grmpbuf(6*cacg_sstep,bid) = z_s
!     rr = r_(k,j)^T r_(k,j)
      grmpbuf(6*cacg_sstep+1,bid) = r_s

   end do ! block loop
   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

 end subroutine btrop_matpower_operator

!***********************************************************************
!BOP
! !IROUTINE: btrop_t_operator
! !INTERFACE:

! !DIR$ ATTRIBUTES FORCEINLINE :: btrop_t_operator
 subroutine btrop_t_operator(bp,P_p)

! !DESCRIPTION:
!  This routine calculate bp = B_k p'_(k,j)
!  B_k is the change basis matrix, with all 0 in column s+1 and 2s+1.
!  Assemble B_k from rho_j(z), rho_j(z) satisfies the three-term recurrence
!  rho_0(z)=1, rho_1(z)=(z−theta_0)rho_0(z)/gamma_0,
!  rho_j(z)=((z−theta_(j−1))rho_(j−1)(z)−sigma_(j−2)rho_(j−2)(z)/gamma_(j−1).  (3)
!  For monomial bases, gamma_j=1, theta_j=sigma_j=0.
!  For Newton bases, ...
!  For Chebyshev bases, ...
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(cacg_smat), intent(in) :: & 
      P_p                   ! short vector

! !OUTPUT PARAMETERS:

   real (r8), dimension(cacg_smat), intent(out) :: &
      bp                    ! B_k p'_(k,j)

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i                   ! dummy counters

!-----------------------------------------------------------------------
   bp(1) = c0
   do i=2,2*cacg_sstep+1,1
      bp(i) = P_p(i-1)
   end do
   bp(cacg_sstep+2) = c0

!-----------------------------------------------------------------------
!EOC

 end subroutine btrop_t_operator

!***********************************************************************

 end module solvers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
