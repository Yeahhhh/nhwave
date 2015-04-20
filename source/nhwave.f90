     program NH_WAVE
!--------------------------------------------------------
!
!    Nonhydrostatic WAVE dynamics
!    
!    Code developer: Gangfeng Ma, University of Delaware
!    Last update: 14/04/2011
!    Last update: 14/08/2011, parallel implementation
!
!-------------------------------------------------------
     use global
     implicit none
     integer :: j,Istage
     real(SP) :: tbegin,tend


     call MPI_INIT(ier)
     call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ier)
     call MPI_COMM_SIZE(MPI_COMM_WORLD,NumP,ier)


     ! record wall time
     call wall_time_secs(tbegin)

     ! read parameter input
     call read_input

     ! work index
     call index
 
     ! allocate variables
     call allocate_variables

     ! generate grids
     call generate_grid

     ! read bathymetry
     call read_bathymetry

     ! initialize model run
     call initial

     ! time integration
     do while (TIME<TOTAL_TIME.and.RUN_STEP<SIM_STEPS)

       ! time step     
       call estimate_dt








       ! update boundary conditions       
       call update_wave_bc

       ! update mask
       call update_mask






       ! update wind
       call update_wind

       ! update vars
       call update_vars

       ! SSP Runge-Kutta time stepping
       do Istage = 1,It_Order

         ! well-balanced source terms
         call source_terms

         ! fluxes at cell faces
         call fluxes

         ! update all variables
         call eval_duvw(Istage)


         ! sponge layer
         if(SPONGE_ON) then
           call sponge_damping
         endif

         ! turbulence model
         if(VISCOUS_FLOW) call eval_turb(Istage)





       enddo


       ! wave average quantities
       if(WAVE_AVERAGE_ON) then
         call wave_average
       endif

       ! screen output
       Screen_Count = Screen_Count+dt
       if(Screen_Count>=Screen_Intv) then
         Screen_Count = Screen_Count-Screen_Intv
         call statistics
       endif

       ! probe output to files
       if(NSTAT>0) then
         Plot_Count_Stat = Plot_Count_Stat+dt
         if(Plot_Count_Stat>=Plot_Intv_Stat) then
           Plot_Count_Stat=Plot_Count_Stat-Plot_Intv_Stat
           call probes
         endif
       endif

       ! field output to files
       if(TIME>=Plot_Start) then
         Plot_Count = Plot_Count+dt
         if(Plot_Count>=Plot_Intv) then
           Plot_Count=Plot_Count-Plot_Intv
           call preview
         endif
       endif
      
     end do

     ! write out wave height and setup
     if(WAVE_AVERAGE_ON) then
       call print_wh_setup
     endif

     if(myid.eq.0) write(*,*) 'Normal Termination!'
     if(myid.eq.0) write(3,*) 'Normal Termination!'

     ! wall time at the end
     call wall_time_secs(tend)

     if(myid.eq.0) write(*,*) 'Simulation takes',tend-tbegin,'seconds'
     if(myid.eq.0) write(3,*) 'Simulation takes',tend-tbegin,'seconds'

     call MPI_FINALIZE(ier)

     end


     subroutine projection_corrector
!-------------------------------------------
!    Correct the velocity field
!    Called by
!       eval_duvw
!    Last update: 25/03/2011, Gangfeng Ma
!-------------------------------------------
     use global
     implicit none
     integer :: i,j,k
     real(SP), dimension(:,:,:),allocatable :: DelxP,DelyP,DelzP
 
     allocate(DelxP(Mloc,Nloc,Kloc))
     allocate(DelyP(Mloc,Nloc,Kloc))
     allocate(DelzP(Mloc,Nloc,Kloc))

     DelxP = Zero
     DelyP = Zero
     DelzP = Zero
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       DelxP(i,j,k) = 0.5*((P(i+1,j,k)-P(i-1,j,k))/(2.*dx)+  &
              (P(i+1,j,k+1)-P(i-1,j,k+1))/(2.*dx))
       DelyP(i,j,k) = 0.5*((P(i,j+1,k)-P(i,j-1,k))/(2.*dy)+  &
              (P(i,j+1,k+1)-P(i,j-1,k+1))/(2.*dy))
       DelzP(i,j,k) = (P(i,j,k+1)-P(i,j,k))/dsig(k)
     enddo
     enddo
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==0) cycle

       DU(i,j,k) = DU(i,j,k)-D(i,j)*dt/Rho0*(DelxP(i,j,k)+DelzP(i,j,k)*DelxSc(i,j,k))
       DV(i,j,k) = DV(i,j,k)-D(i,j)*dt/Rho0*(DelyP(i,j,k)+DelzP(i,j,k)*DelySc(i,j,k))
       DW(i,j,k) = DW(i,j,k)-dt/Rho0*DelzP(i,j,k)
     enddo
     enddo
     enddo

     deallocate(DelxP)
     deallocate(DelyP)
     deallocate(DelzP)

     return
     end subroutine projection_corrector
     

     subroutine poisson_solver
!--------------------------------------------
!    Solve poisson equation for dynamic pressure
!    Called by
!       eval_duvw
!    Last update: 24/03/2011, Gangfeng Ma
!----------------------------------------------
     use global
     implicit none
     integer :: i,j,k,imask

     ! generate coefficient matrix and rhs
     call generate_coef_rhs

     ! use HYPRE package for parallel computation
     call hypre_pres_solver
  
     ! fyshi gave boundary condition for dry cells
     ! set zero for dry set is inaccurate
     ! dry cells
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==0) then
         P(i,j,k) = Zero
         
         ! south boundary 
         if(Mask(i,j+1)==1)then
           do imask=1,Nghost
             P(i,j-imask+1,k)=P(i,j+imask,k)
           enddo
         ! north boundary
         elseif(Mask(i,j-1)==1)then
           do imask=1,Nghost
             P(i,j+imask-1,k)=P(i,j-imask,k)
           enddo
         ! west boundary
         elseif(Mask(i+1,j)==1)then
           do imask=1,Nghost
             P(i-imask+1,j,k)=P(i+imask,j,k)
           enddo
         ! east boundary
         elseif(Mask(i-1,j)==1)then
           do imask=1,Nghost
             P(i+imask-1,j,k)=P(i-imask,j,k)
           enddo
         endif
       endif 
     enddo
     enddo
     enddo

     ! collect into ghost cells
     if(n_west.eq.MPI_PROC_NULL) then
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       do i = 1,Nghost
         P(Ibeg-i,j,k) = P(Ibeg+i-1,j,k)
       enddo
     enddo
     enddo
     endif

     if(n_east.eq.MPI_PROC_NULL) then
     do k = Kbeg,Kend
     do j = Jbeg,Jend 
       do i = 1,Nghost     
         P(Iend+i,j,k) = P(Iend-i+1,j,k)
       enddo
     enddo
     enddo
     endif

     if(n_suth.eq.MPI_PROC_NULL) then
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       do j = 1,Nghost
         P(i,Jbeg-j,k) = P(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
     endif

     if(n_nrth.eq.MPI_PROC_NULL) then
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       do j = 1,Nghost
         P(i,Jend+j,k) = P(i,Jend-j+1,k)
       enddo
     enddo
     enddo
     endif

     call phi_3D_exch(P)

     end subroutine poisson_solver


     subroutine hypre_pres_solver
!---------------------------------------------
!    solve pressure using hypre package
!    called by
!       poisson_solver
!    Last update: 22/08/2011, Gangfeng Ma
!---------------------------------------------
     use global
     implicit none
     integer, parameter :: ndim=3
     integer, parameter :: nentries=15
     integer :: i,j,k,n,ivalues,nvalues,neq,ientry,num_iterations,  &
                precond_id,n_pre,n_post,ierr
     integer*8 :: grid,stencil,matrix,vec_b,vec_x,solver,precond
     integer :: i_glob(Mloc),j_glob(Nloc),k_glob(Kloc)
     integer :: ilower(ndim),iupper(ndim),offsets(nentries,ndim),stencil_indices(nentries), &
                periodic_shift(ndim)
     real(SP) :: final_res_norm
     real(SP), dimension(:), allocatable :: values,Phi
     integer, dimension(:,:,:), allocatable :: indx 
     data ((offsets(i,j),j=1,ndim),i=1,nentries)/0,0,0,1,0,0,0,1,0,0,-1,1,-1,0,1,  &
             0,0,1,1,0,1,0,1,1,-1,0,0,0,-1,0,  &
             0,1,-1,1,0,-1,0,0,-1,-1,0,-1,0,-1,-1/

     ! set up a three dimensional grid
     call HYPRE_StructGridCreate(MPI_COMM_WORLD,ndim,grid,ierr)

     ! global indices
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       i_glob(i) = npx*(Iend-Ibeg+1)+i-Nghost
       j_glob(j) = npy*(Jend-Jbeg+1)+j-Nghost
       k_glob(k) = k-Nghost
     enddo
     enddo
     enddo

     ilower(1) = i_glob(Ibeg)
     ilower(2) = j_glob(Jbeg)
     ilower(3) = k_glob(Kbeg)
     iupper(1) = i_glob(Iend)
     iupper(2) = j_glob(Jend)
     iupper(3) = k_glob(Kend)

     call HYPRE_StructGridSetExtents(grid,ilower,iupper,ierr)

     if(PERIODIC_X.or.PERIODIC_Y) then
       if(PERIODIC_X) then
         periodic_shift(1) = Mglob
       else
         periodic_shift(1) = 0
       endif
       if(PERIODIC_Y) then
         periodic_shift(2) = Nglob
       else
         periodic_shift(2) = 0
       endif
       periodic_shift(3) = 0
       call HYPRE_StructGridSetPeriodic(grid,periodic_shift,ierr)
     endif

     call HYPRE_StructGridAssemble(grid,ierr)

     ! define the discretization stencil
     call HYPRE_StructStencilCreate(ndim,nentries,stencil,ierr)

     do ientry = 1,nentries
       call HYPRE_StructStencilSetElement(stencil,(ientry-1),offsets(ientry,:),ierr)
     enddo

     ! create matrix object
     call HYPRE_StructMatrixCreate(MPI_COMM_WORLD,grid,stencil,matrix,ierr)

     call HYPRE_StructMatrixInitialize(matrix,ierr)

     ! set the matrix coefficient
     do i = 1,nentries
       stencil_indices(i) = i-1
     enddo

     allocate(indx(Mloc,Nloc,Kloc))
 
     neq = 0
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       neq = neq+1
       indx(i,j,k) = neq
     enddo
     enddo
     enddo
    
     nvalues = (Iend-Ibeg+1)*(Jend-Jbeg+1)*(Kend-Kbeg+1)*nentries
     allocate(values(nvalues))

     ivalues = 0
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       do n = 1,nentries
         ivalues = ivalues+1
         values(ivalues) = Coef(indx(i,j,k),n)
       enddo
     enddo
     enddo
     enddo

     call HYPRE_StructMatrixSetBoxValues(matrix,ilower,iupper,nentries,  &
                                  stencil_indices,values,ierr) 
     call HYPRE_StructMatrixAssemble(matrix,ierr)
     !call HYPRE_StructMatrixPrint(matrix,zero,ierr)

     ! set up struct vectors for b and x
     call HYPRE_StructVectorCreate(MPI_COMM_WORLD,grid,vec_b,ierr)
     call HYPRE_StructVectorCreate(MPI_COMM_WORLD,grid,vec_x,ierr)

     call HYPRE_StructVectorInitialize(vec_b,ierr)
     call HYPRE_StructVectorInitialize(vec_x,ierr)

     ! set the vector coefficients
     call HYPRE_StructVectorSetBoxValues(vec_b,ilower,iupper,Rhs,ierr)   
     call HYPRE_StructVectorAssemble(vec_b,ierr)     
     !call HYPRE_StructVectorPrint(vec_b,zero,ierr)

     ! initial guess
     allocate(Phi(neqns))
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Phi(indx(i,j,k)) = P(i,j,k)
     enddo
     enddo
     enddo
     
     call HYPRE_StructVectorSetBoxValues(vec_x,ilower,iupper,Phi,ierr)
     call HYPRE_StructVectorAssemble(vec_x,ierr)
     !call HYPRE_StructVectorPrint(vec_x,zero,ierr)

     ! set up and use a solver
     call HYPRE_StructGMRESCreate(MPI_COMM_WORLD,solver,ierr)
     call HYPRE_StructGMRESSetMaxIter(solver,itmax,ierr)
     call HYPRE_StructGMRESSetTol(solver,tol,ierr)
     call HYPRE_StructGMRESSetPrintLevel(solver,0,ierr)
     call HYPRE_StructGMRESSetLogging(solver,0,ierr)

     ! use symmetric SMG as preconditioner
     n_pre = 1; n_post = 1
     call HYPRE_StructSMGCreate(MPI_COMM_WORLD,precond,ierr)
     call HYPRE_StructSMGSetMemoryUse(precond,0,ierr)
     call HYPRE_StructSMGSetMaxIter(precond,1,ierr)
     call HYPRE_StructSMGSetTol(precond,0.0,ierr)
     call HYPRE_StructSMGSetNumPreRelax(precond,n_pre,ierr)
     call HYPRE_StructSMGSetNumPostRelax(precond,n_post,ierr)
     call HYPRE_StructSMGSetLogging(precond,0,ierr)

     ! set up preconditioner
     precond_id = 0
     call HYPRE_StructGMRESSetPrecond(solver,precond_id,precond,ierr)
     
     ! do the setup
     call HYPRE_StructGMRESSetup(solver,matrix,vec_b,vec_x,ierr)
 
     ! do the solve
     call HYPRE_StructGMRESSolve(solver,matrix,vec_b,vec_x,ierr)

     ! get results
     call HYPRE_StructVectorGetBoxValues(vec_x,ilower,iupper,Phi,ierr)

     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       P(i,j,k) = Phi(indx(i,j,k))
     enddo
     enddo
     enddo

     ! get some info
     !call HYPRE_StructGMRESGetFinalRelati(solver,final_res_norm,ierr)
     !call HYPRE_StructGMRESGetNumIteratio(solver,num_iterations,ierr);
     !
     !if(myid.eq.0) then
     !  write(*,*)'Iterations = ',num_iterations
     !  write(*,*)'Final Relative Residual Norm = ',final_res_norm
     !endif

     ! free memory
     call HYPRE_StructGridDestroy(grid,ierr)
     call HYPRE_StructStencilDestroy(stencil,ierr)
     call HYPRE_StructMatrixDestroy(matrix,ierr)
     call HYPRE_StructVectorDestroy(vec_b,ierr)
     call HYPRE_StructVectorDestroy(vec_x,ierr)
     call HYPRE_StructGMRESDestroy(solver,ierr)
     call HYPRE_StructSMGDestroy(precond,ierr)

     deallocate(indx)
     deallocate(values)
     deallocate(Phi)

     return
     end subroutine hypre_pres_solver


     subroutine generate_coef_rhs
!---------------------------------------------
!    Generate coefficient matrix and rhs
!    Called by 
!       poisson_solver
!    Last update: 24/03/2011, Gangfeng Ma
!--------------------------------------------
     use global
     implicit none
     integer :: i,j,k,neq,n,ic
     real(SP), dimension(:,:,:), allocatable :: DelxS,DelyS,DelzS,A1
     integer,  dimension(:,:,:), allocatable :: indx

     allocate(DelxS(Mloc,Nloc,Kloc1))
     allocate(DelyS(Mloc,Nloc,Kloc1))
     allocate(DelzS(Mloc,Nloc,Kloc1))
     allocate(A1(Mloc,Nloc,Kloc1))
     allocate(indx(Mloc,Nloc,Kloc))

     DelxS = Zero
     DelyS = Zero
     DelzS = Zero
     A1 = Zero
     do k = Kbeg,Kend1
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       DelxS(i,j,k) = (1.-sig(k))/D(i,j)*DelxH(i,j)-sig(k)/D(i,j)*DelxEta(i,j)
       DelyS(i,j,k) = (1.-sig(k))/D(i,j)*DelyH(i,j)-sig(k)/D(i,j)*DelyEta(i,j)
       DelzS(i,j,k) = 1./D(i,j)

       A1(i,j,k) = DelxS(i,j,k)*DelxS(i,j,k)+DelyS(i,j,k)*DelyS(i,j,k)+  &
            DelzS(i,j,k)*DelzS(i,j,k)
     enddo
     enddo
     enddo
   
     ! generate coefficient matrix
     neq = 0
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       neq = neq+1
       indx(i,j,k) = neq
     enddo
     enddo 
     enddo

     ! generate source term 
     Rhs = Zero
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Rhs(indx(i,j,k)) = -((Uf(i+1,j,k)-Uf(i-1,j,k))/(2.0*dx)+(U(i,j,k)-U(i,j,k-1))/(0.5*(dsig(k)+dsig(k-1)))*  &
              DelxS(i,j,k)+(Vf(i,j+1,k)-Vf(i,j-1,k))/(2.0*dy)+(V(i,j,k)-V(i,j,k-1))/(0.5*(dsig(k)+dsig(k-1)))*  &
              DelyS(i,j,k)+(W(i,j,k)-W(i,j,k-1))/(0.5*(dsig(k)+dsig(k-1)))*DelzS(i,j,k)-SourceC(i,j))*Rho0/dt
     enddo
     enddo
     enddo

     Coef = Zero
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Coef(indx(i,j,k),1) = (2./(dx*dx)+2./(dy*dy)+A1(i,j,k)/(0.5*(dsig(k)+dsig(k-1))*dsig(k))+  &
                A1(i,j,k)/(0.5*(dsig(k)+dsig(k-1))*dsig(k-1)))
       Coef(indx(i,j,k),2) = -1./(dx*dx)
       Coef(indx(i,j,k),3) = -1./(dy*dy)
       Coef(indx(i,j,k),4) = (DelyS(i,j-1,k)/(2.*dy*(dsig(k)+dsig(k-1)))+DelyS(i,j,k)/(2.*dy*(dsig(k)+dsig(k-1))))   
       Coef(indx(i,j,k),5) = (DelxS(i-1,j,k)/(2.*dx*(dsig(k)+dsig(k-1)))+DelxS(i,j,k)/(2.*dx*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),6) = -A1(i,j,k)/(0.5*(dsig(k)+dsig(k-1))*dsig(k))
       Coef(indx(i,j,k),7) = -(DelxS(i+1,j,k)/(2.*dx*(dsig(k)+dsig(k-1)))+DelxS(i,j,k)/(2.*dx*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),8) = -(DelyS(i,j+1,k)/(2.*dy*(dsig(k)+dsig(k-1)))+DelyS(i,j,k)/(2.*dy*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),9) = -1./(dx*dx)
       Coef(indx(i,j,k),10) = -1./(dy*dy)
       Coef(indx(i,j,k),11) = (DelyS(i,j+1,k)/(2.*dy*(dsig(k)+dsig(k-1)))+DelyS(i,j,k)/(2.*dy*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),12) = (DelxS(i+1,j,k)/(2.*dx*(dsig(k)+dsig(k-1)))+DelxS(i,j,k)/(2.*dx*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),13) = -A1(i,j,k)/(0.5*(dsig(k)+dsig(k-1))*dsig(k-1))
       Coef(indx(i,j,k),14) = -(DelxS(i-1,j,k)/(2.*dx*(dsig(k)+dsig(k-1)))+DelxS(i,j,k)/(2.*dx*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),15) = -(DelyS(i,j-1,k)/(2.*dy*(dsig(k)+dsig(k-1)))+DelyS(i,j,k)/(2.*dy*(dsig(k)+dsig(k-1))))
     enddo
     enddo
     enddo

     ! fyshi added boundary conditions at masks face 02/15/2013
     do i = Ibeg+1,Iend-1
     do j = Jbeg+1,Jend-1
     do k = Kbeg,Kend
       if(mask(i,j)==0) then
         ! left 
         if(mask(i+1,j)==1) then
           ic = indx(I+1,j,k)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,9)
           Coef(ic,6) = Coef(ic,6)+Coef(ic,5)
           Coef(ic,13) = Coef(ic,13)+Coef(ic,14)
           Coef(ic,9) = Zero
           Coef(ic,5) = Zero
           Coef(ic,14) = Zero
         ! right 
         elseif(mask(i-1,j)==1) then
           ic = indx(I-1,j,k)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,2)
           Coef(ic,6) = Coef(ic,6)+Coef(ic,7)
           Coef(ic,13) = Coef(ic,13)+Coef(ic,12)
           Coef(ic,2) = Zero
           Coef(ic,7) = Zero
           Coef(ic,12) = Zero
         ! south
         elseif(mask(i,j+1)==1) then
           ic = indx(i,J+1,k)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,10)
           Coef(ic,6) = Coef(ic,6)+Coef(ic,4)
           Coef(ic,13) = Coef(ic,13)+Coef(ic,15)
           Coef(ic,10) = Zero
           Coef(ic,4) = Zero
           Coef(ic,15) = Zero
         ! north
         elseif(mask(i,j-1)==1) then
           ic = indx(i,J-1,k)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,3)
           Coef(ic,6) = Coef(ic,6)+Coef(ic,8)
           Coef(ic,13) = Coef(ic,13)+Coef(ic,11)
           Coef(ic,3) = Zero
           Coef(ic,8) = Zero
           Coef(ic,11) = Zero
         endif ! end mask+1=1 
       endif ! end mask=0
     enddo
     enddo
     enddo


     ! boundary conditions
     ! left side
     if(n_west.eq.MPI_PROC_NULL) then
     i = Ibeg
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       ic = indx(i,j,k)
       Coef(ic,1) = Coef(ic,1)+Coef(ic,9)
       Coef(ic,6) = Coef(ic,6)+Coef(ic,5)
       Coef(ic,13) = Coef(ic,13)+Coef(ic,14)
       Coef(ic,9) = Zero
       Coef(ic,5) = Zero
       Coef(ic,14) = Zero
     enddo
     enddo
     endif

     ! right side
     if(n_east.eq.MPI_PROC_NULL) then
     i = Iend
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       ic = indx(i,j,k)
       Coef(ic,1) = Coef(ic,1)+Coef(ic,2)
       Coef(ic,6) = Coef(ic,6)+Coef(ic,7)
       Coef(ic,13) = Coef(ic,13)+Coef(ic,12)
       Coef(ic,2) = Zero
       Coef(ic,7) = Zero
       Coef(ic,12) = Zero
     enddo
     enddo
     endif

     ! front side
     if(n_suth.eq.MPI_PROC_NULL) then
     j = Jbeg
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       ic = indx(i,j,k)         
       Coef(ic,1) = Coef(ic,1)+Coef(ic,10)
       Coef(ic,6) = Coef(ic,6)+Coef(ic,4)
       Coef(ic,13) = Coef(ic,13)+Coef(ic,15)
       Coef(ic,10) = Zero
       Coef(ic,4) = Zero
       Coef(ic,15) = Zero
     enddo
     enddo
     endif

     ! back side
     if(n_nrth.eq.MPI_PROC_NULL) then
     j = Jend
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       ic = indx(i,j,k)
       Coef(ic,1) = Coef(ic,1)+Coef(ic,3)
       Coef(ic,6) = Coef(ic,6)+Coef(ic,8)
       Coef(ic,13) = Coef(ic,13)+Coef(ic,11)
       Coef(ic,3) = Zero
       Coef(ic,8) = Zero
       Coef(ic,11) = Zero
     enddo
     enddo
     endif

     ! bottom side
     k = Kbeg
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       ic = indx(i,j,k)

!# if defined (TWOLAYERSLIDE)
!       Rhs(ic) = Rhs(ic)+Rho0*(dsig(Kbeg)+dsig(Kbeg-1))*(Coef(ic,13)*D(i,j)*Delt2H(i,j)+ &             
!          Coef(ic,12)*D(i+1,j)*Delt2H(i+1,j)+Coef(ic,11)*D(i,j+1)*Delt2H(i,j+1)+ &
!          Coef(ic,14)*D(i-1,j)*Delt2H(i-1,j)+Coef(ic,15)*D(i,j-1)*Delt2H(i,j-1))
!# endif

       Coef(ic,6) = Coef(ic,6)+Coef(ic,13)
       Coef(ic,7) = Coef(ic,7)+Coef(ic,12)
       Coef(ic,8) = Coef(ic,8)+Coef(ic,11)
       Coef(ic,5) = Coef(ic,5)+Coef(ic,14)
       Coef(ic,4) = Coef(ic,4)+Coef(ic,15)
       Coef(ic,13) = Zero
       Coef(ic,12) = Zero
       Coef(ic,11) = Zero
       Coef(ic,14) = Zero
       Coef(ic,15) = Zero
     enddo
     enddo

     ! top side (Dirichlet boundary)
     k = Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       ic = indx(i,j,k)
       Coef(ic,4) = Zero
       Coef(ic,5) = Zero
       Coef(ic,6) = Zero
       Coef(ic,7) = Zero
       Coef(ic,8) = Zero
     enddo
     enddo

     ! take (i=2,j=2,k=2) to obtain the diagonal information
     JCoef(1) = indx(Ibeg+1,Jbeg+1,Kbeg+1)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i,j,k)
     JCoef(2) = indx(Ibeg+2,Jbeg+1,Kbeg+1)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i+1,j,k) 
     JCoef(3) = indx(Ibeg+1,Jbeg+2,Kbeg+1)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i,j+1,k)
     JCoef(4) = indx(Ibeg+1,Jbeg,Kbeg+2)-indx(Ibeg+1,Jbeg+1,Kbeg+1)    ! (i,j-1,k+1)
     JCoef(5) = indx(Ibeg,Jbeg+1,Kbeg+2)-indx(Ibeg+1,Jbeg+1,Kbeg+1)    ! (i-1,j,k+1)
     JCoef(6) = indx(Ibeg+1,Jbeg+1,Kbeg+2)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i,j,k+1)
     JCoef(7) = indx(Ibeg+2,Jbeg+1,Kbeg+2)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i+1,j,k+1)
     JCoef(8) = indx(Ibeg+1,Jbeg+2,Kbeg+2)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i,j+1,k+1)
     JCoef(9) = indx(Ibeg,Jbeg+1,Kbeg+1)-indx(Ibeg+1,Jbeg+1,Kbeg+1)    ! (i-1,j,k)
     JCoef(10) = indx(Ibeg+1,Jbeg,Kbeg+1)-indx(Ibeg+1,Jbeg+1,Kbeg+1)   ! (i,j-1,k)
     JCoef(11) = indx(Ibeg+1,Jbeg+2,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+1)   ! (i,j+1,k-1)
     JCoef(12) = indx(Ibeg+2,Jbeg+1,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+1)   ! (i+1,j,k-1)
     JCoef(13) = indx(Ibeg+1,Jbeg+1,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+1)   ! (i,j,k-1)
     JCoef(14) = indx(Ibeg,Jbeg+1,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+1)     ! (i-1,j,k-1)
     JCoef(15) = indx(Ibeg+1,Jbeg,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+1)     ! (i,j-1,k-1)

     deallocate(DelxS)
     deallocate(DelyS)
     deallocate(DelzS)
     deallocate(A1) 
     deallocate(indx)

     return
     end subroutine generate_coef_rhs


     subroutine eval_duvw(ISTEP)
!-----------------------------------------------
!    Update all variables D,U,V,W,Omega
!    Called by
!       main
!    Last update: 25/12/2010, Gangfeng Ma
!----------------------------------------------
     use global
     implicit none
     integer,intent(in) :: ISTEP
     real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0
     real(SP),dimension(:,:),allocatable :: R1,qz
     real(SP),dimension(:,:,:),allocatable :: R2,R3,R4,RhsU,RhsV,RhsW
     real(SP) :: dedt,Umag,Dz1,Cdrag,Ustar2,Wtop,Wbot
     integer :: i,j,k,n,Ista,Nlen

     Nlen = Kend-Kbeg+1

     allocate(qz(Mloc,Nloc))
     allocate(R1(Mloc,Nloc))
     allocate(R2(Mloc,Nloc,Kloc))
     allocate(R3(Mloc,Nloc,Kloc))
     allocate(R4(Mloc,Nloc,Kloc))
     allocate(Acoef(Nlen))
     allocate(Bcoef(Nlen))
     allocate(Ccoef(Nlen))
     allocate(Xsol(Nlen))
     allocate(Rhs0(Nlen))

     ! calculate baroclinic pressure gradient
     if(.not.BAROTROPIC) call baropg_z

     ! estimate horizontal diffusion terms
     if(VISCOUS_FLOW) call diffusion

     ! external forcing
     if(EXTERNAL_FORCING) call driver(ExtForceX,ExtForceY)


     ! solve total water depth D
     R1 = Zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==0) cycle
       do k = Kbeg,Kend
         R1(i,j) = R1(i,j)-1.0/dx*(Ex(i+1,j,k)-Ex(i,j,k))*dsig(k)  &
                        -1.0/dy*(Ey(i,j+1,k)-Ey(i,j,k))*dsig(k)
       enddo
       ! internal wavemaker
       R1(i,j) = R1(i,j)+D(i,j)*SourceC(i,j)
       D(i,j) = ALPHA(ISTEP)*D0(i,j)+BETA(ISTEP)*(D(i,j)+dt*R1(i,j)) 
     enddo
     enddo

     ! update D and Eta          
     D = max(D,MinDep)
     call wl_bc
     Eta = D-Hc

     call delxFun_2D(Eta,DelxEta)
     call delyFun_2D(Eta,DelyEta)

     ! prepare right-hand side terms
     R2 = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Mask(i,j)==0) cycle  
       do k = Kbeg,Kend
         R2(i,j,k) = -1.0/dx*(Fx(i+1,j,k)-Fx(i,j,k))-1.0/dy*(Fy(i,j+1,k)-Fy(i,j,k)) &
                -1.0/dsig(k)*(Fz(i,j,k+1)-Fz(i,j,k))+fcor*DV(i,j,k)+DRhoX(i,j,k)+  &
                SourceX(i,j)+Diffxx(i,j,k)+Diffxy(i,j,k)+ExtForceX(i,j,k)
       enddo
     enddo
     enddo

     R3 = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Mask(i,j)==0) cycle
       do k = Kbeg,Kend
         R3(i,j,k) = -1.0/dx*(Gx(i+1,j,k)-Gx(i,j,k))-1.0/dy*(Gy(i,j+1,k)-Gy(i,j,k)) &   
                        -1.0/dsig(k)*(Gz(i,j,k+1)-Gz(i,j,k))-fcor*DU(i,j,k)+DRhoY(i,j,k)  &                                                
                        +SourceY(i,j)+Diffyx(i,j,k)+Diffyy(i,j,k)+ExtForceY(i,j,k)                                                           
       enddo
     enddo
     enddo

     R4 = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Mask(i,j)==0) cycle
       do k = Kbeg,Kend
         R4(i,j,k) = -1.0/dx*(Hx(i+1,j,k)-Hx(i,j,k))-1.0/dy*(Hy(i,j+1,k)-Hy(i,j,k)) &  
                        -1.0/dsig(k)*(Hz(i,j,k+1)-Hz(i,j,k))+Diffzx(i,j,k)+Diffzy(i,j,k)
       enddo
     enddo
     enddo


     ! solve DU
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Mask(i,j)==0) cycle

       if(VISCOUS_FLOW) then
         Nlen = 0
         do k = Kbeg,Kend
           Nlen = Nlen+1
           if(k==Kbeg) then
             Acoef(Nlen) = 0.0
           else
             Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+  &
                     Cmu(i,j,k)+CmuR(i,j,k-1)+CmuR(i,j,k))+  &
                     0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schmidt)/  &
                     (0.5*dsig(k)*(dsig(k)+dsig(k-1)))
           endif

           if(k==Kend) then
             Ccoef(Nlen) = 0.0
           else
             Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1)+  &
                     CmuR(i,j,k)+CmuR(i,j,k+1))+  &
                     0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schmidt)/  &
                     (0.5*dsig(k)*(dsig(k)+dsig(k+1)))
           endif
        
           if(k==Kbeg.and.Bc_Z0==2) then  ! no-slip
             Bcoef(Nlen) = 1.0-Ccoef(Nlen)+  &
                dt/D(i,j)**2*(Cmu(i,j,k)+CmuR(i,j,k)+  &
                CmuVt(i,j,k)/Schmidt)/(0.5*dsig(k)*dsig(k))
           else
             Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)
           endif


           if(k==Kbeg.and.Bc_Z0==5) then  ! friction law
             Dz1 = 0.5*D(i,j)*dsig(Kbeg)
             if(ibot==1) then
               Cdrag = Cd0
             else
               Cdrag = 1./(1./Kappa*log(30.0*Dz1/Zob))**2
             endif
             Ustar2 = Cdrag*sqrt(U(i,j,k)**2+V(i,j,k)**2)*U(i,j,k)
             Rhs0(Nlen) = DU(i,j,k)+dt*R2(i,j,k)-dt*Ustar2/dsig(k)
           elseif(k==Kend) then
             Rhs0(Nlen) = DU(i,j,k)+dt*R2(i,j,k)+dt*Wsx(i,j)/dsig(k) 
           else
             Rhs0(Nlen) = DU(i,j,k)+dt*R2(i,j,k)
           endif

         enddo
      
         call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

         Nlen = 0
         do k = Kbeg,Kend
           Nlen = Nlen+1
           DU(i,j,k) = Xsol(Nlen)
         enddo
       else
         do k = Kbeg,Kend
           DU(i,j,k) = DU(i,j,k)+dt*R2(i,j,k)
         enddo
       endif
     enddo
     enddo

     ! solve DV
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Mask(i,j)==0) cycle

       if(VISCOUS_FLOW) then
         Nlen = 0
         do k = Kbeg,Kend
           Nlen = Nlen+1
           if(k==Kbeg) then
             Acoef(Nlen) = 0.0
           else
             Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k)+  &
                  CmuR(i,j,k-1)+CmuR(i,j,k))+  &
                  0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schmidt)/  &
                  (0.5*dsig(k)*(dsig(k)+dsig(k-1)))  
           endif

           if(k==Kend) then
             Ccoef(Nlen) = 0.0
           else
             Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+  &
                  Cmu(i,j,k+1)+CmuR(i,j,k)+CmuR(i,j,k+1))+  &
                  0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schmidt)/  &
                  (0.5*dsig(k)*(dsig(k)+dsig(k+1)))   
           endif

           if(k==Kbeg.and.Bc_Z0==2) then  ! no-slip                                             
             Bcoef(Nlen) = 1.0-Ccoef(Nlen)+  &
                dt/D(i,j)**2*(Cmu(i,j,k)+CmuR(i,j,k)+  &
                CmuVt(i,j,k)/Schmidt)/(0.5*dsig(k)*dsig(k))
           else
             Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)
           endif


           if(k==Kbeg.and.Bc_Z0==5) then
             Dz1 = 0.5*D(i,j)*dsig(Kbeg)
             if(ibot==1) then
               Cdrag = Cd0
             else
               Cdrag = 1./(1./Kappa*log(30.0*Dz1/Zob))**2
             endif
             Ustar2 = Cdrag*sqrt(U(i,j,k)**2+V(i,j,k)**2)*V(i,j,k)

             Rhs0(Nlen) = DV(i,j,k)+dt*R3(i,j,k)-dt*Ustar2/dsig(k)
           elseif(k==Kend) then
             Rhs0(Nlen) = DV(i,j,k)+dt*R3(i,j,k)+dt*Wsy(i,j)/dsig(k)
           else
             Rhs0(Nlen) = DV(i,j,k)+dt*R3(i,j,k)
           endif


         enddo

         call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

         Nlen = 0
         do k = Kbeg,Kend
           Nlen = Nlen+1
           DV(i,j,k) = Xsol(Nlen)
         enddo
       else
         do k = Kbeg,Kend
           DV(i,j,k) = DV(i,j,k)+dt*R3(i,j,k)
         enddo
       endif
     enddo
     enddo

     ! solve DW
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(Mask(i,j)==0) cycle

       if(VISCOUS_FLOW) then
         Nlen = 0
         do k = Kbeg,Kend
           Nlen = Nlen+1
           if(k==Kbeg) then
             Acoef(Nlen) = 0.0
           else
             Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+  &
                  Cmu(i,j,k)+CmuR(i,j,k-1)+CmuR(i,j,k))+  &
                  0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schmidt)/  &
                  (0.5*dsig(k)*(dsig(k)+dsig(k-1))) 
           endif
 
           if(k==Kend) then
             Ccoef(Nlen) = 0.0
           else
             Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+  &
                 Cmu(i,j,k+1)+CmuR(i,j,k)+CmuR(i,j,k+1))+  &
                 0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schmidt)/  &
                 (0.5*dsig(k)*(dsig(k)+dsig(k+1)))    
           endif

           if(k==Kbeg) then
             Bcoef(Nlen) = 1.0-Ccoef(Nlen)+  &
                 dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k)+CmuR(i,j,k-1)+CmuR(i,j,k))+  &
                 0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schmidt)/(0.25*dsig(k)*(dsig(k)+dsig(k-1)))
           elseif(k==Kend) then
             Bcoef(Nlen) = 1.0-Acoef(Nlen)+  &
                 dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1)+CmuR(i,j,k)+CmuR(i,j,k+1))+  &
                 0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schmidt)/(0.25*dsig(k)*(dsig(k)+dsig(k+1)))
           else
             Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)
           endif


           if(k==Kbeg) then
             Wbot = -DeltH(i,j)-U(i,j,Kbeg)*DelxH(i,j)-V(i,j,Kbeg)*DelyH(i,j) 
             Rhs0(Nlen) = DW(i,j,k)+dt*R4(i,j,k)+  &
                dt/D(i,j)*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k)+CmuR(i,j,k-1)+  &
                CmuR(i,j,k))+0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schmidt)/  &
                (0.25*dsig(k)*(dsig(k)+dsig(k-1)))*Wbot
           elseif(k==Kend) then
             Wtop = (Eta(i,j)-Eta0(i,j))/dt+U(i,j,Kend)*DelxEta(i,j)+V(i,j,Kend)*DelyEta(i,j) 
             Rhs0(Nlen) = DW(i,j,k)+dt*R4(i,j,k)+  &
                dt/D(i,j)*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1)+CmuR(i,j,k)+  &
                CmuR(i,j,k+1))+0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schmidt)/  &
                (0.25*dsig(k)*(dsig(k)+dsig(k+1)))*Wtop
           else
             Rhs0(Nlen) = DW(i,j,k)+dt*R4(i,j,k)
           endif


         enddo

         call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

         Nlen = 0
         do k = Kbeg,Kend
           Nlen = Nlen+1
           DW(i,j,k) = Xsol(Nlen)
         enddo
       else
         do k = Kbeg,Kend
           DW(i,j,k) = DW(i,j,k)+dt*R4(i,j,k)
         enddo 
       endif
     enddo
     enddo

     ! sigma transformation coefficient 
     call sigma_transform

     ! run non-hydrostatic simulation  
     if(NON_HYDRO) then
       ! obtain hydrostatic velocity
       call get_UVW

       ! interpolate velocity into vertical faces
       call interpolate_velocity_to_faces

       ! solve dynamic pressure 
       call poisson_solver

       ! correct velocity field  
       call projection_corrector
     endif

     ! SSP Runge-Kutta time stepping
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       DU(i,j,k) = ALPHA(ISTEP)*DU0(i,j,k)+BETA(ISTEP)*DU(i,j,k)
       DV(i,j,k) = ALPHA(ISTEP)*DV0(i,j,k)+BETA(ISTEP)*DV(i,j,k)
       DW(i,j,k) = ALPHA(ISTEP)*DW0(i,j,k)+BETA(ISTEP)*DW(i,j,k)

       if(Mask(i,j)==0) then
         DU(i,j,k) = Zero
         DV(i,j,k) = Zero
         DW(i,j,k) = Zero
       endif
     enddo
     enddo
     enddo

     ! boundary conditions and final velocity
     call get_UVW  

     ! update Omega
     call get_Omega(R1)

     ! if running hydrostatic mode, replace vertical velocity
     ! in fact, W is useless. Only for output
     if(.not.NON_HYDRO) then
       do i = Ibeg,Iend
       do j = Jbeg,Jend
       do k = Kbeg,Kend
         W(i,j,k) = 0.5*(Omega(i,j,k)+Omega(i,j,k+1))-  &
              DeltH(i,j)+sigc(k)*R1(i,j)-U(i,j,k)*  &
              ((1.0-sigc(k))*DelxH(i,j)+sigc(k)*DelxEta(i,j))-  &
              V(i,j,k)*((1.0-sigc(k))*DelyH(i,j)+sigc(k)*DelyEta(i,j))
         if(Mask(i,j)==0) W(i,j,k) = Zero
         DW(i,j,k) = D(i,j)*W(i,j,k)         
       enddo
       enddo
       enddo

       ! update velocity field
       call get_UVW
     endif

     deallocate(qz)
     deallocate(R1)
     deallocate(R2)
     deallocate(R3)
     deallocate(R4)
     deallocate(Acoef)
     deallocate(Bcoef)
     deallocate(Ccoef)
     deallocate(Xsol)
     deallocate(Rhs0)

     end subroutine eval_duvw


     subroutine driver(ExtForceX,ExtForceY)
!--------------------------------------------------------------------------
!    Specify external forcing
!    Called by 
!       eval_duvw
!    Last Update: 18/07/2012, Gangfeng Ma
!--------------------------------------------------------------------------
     use global, only: SP,Zero,Mloc,Nloc,Kloc,Ibeg,Iend,Jbeg,Jend,Kbeg,Kend,D
     implicit none
     real(SP), dimension(Mloc,Nloc,Kloc), intent(inout) :: ExtForceX,ExtForceY
     integer :: i,j,k
    
     ExtForceX = Zero
     ExtForceY = Zero

     ! specify energy slope for open channel flow
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       ExtForceX(i,j,k) = 0.02*D(i,j)
     enddo
     enddo
     enddo

     return
     end subroutine driver


     subroutine interpolate_velocity_to_faces
!------------------------------------------------                                        
!    Interpolate U,V,W to vertical faces                                                 
!    Called by                                                                           
!       main                                                                             
!    Last Update: 19/03/2011, Gangfeng Ma                                                
!------------------------------------------------                                        
     use global, only: SP,U,V,W,Uf,Vf,Wf,dsig,sig,sigc,  &
                       Mloc,Nloc,Kloc,Kloc1,Nghost
     implicit none
     integer  :: i,j,k
     real(SP) :: Int_factor,Int_factor1,Int_factor2,Int_factor3
     logical  :: Linear_Interp=.True.

     if(Linear_Interp) then
       ! first-order linear interpolation
       do k = 2,Kloc
       do j = 1,Nloc
       do i = 1,Mloc
         Int_factor = dsig(k)/(dsig(k)+dsig(k-1))
         Uf(i,j,k) = (1.0-Int_factor)*U(i,j,k)+Int_factor*U(i,j,k-1)
         Vf(i,j,k) = (1.0-Int_factor)*V(i,j,k)+Int_factor*V(i,j,k-1)
         Wf(i,j,k) = (1.0-Int_factor)*W(i,j,k)+Int_factor*W(i,j,k-1)
       enddo
       enddo
       enddo

       do j = 1,Nloc
       do i = 1,Mloc
         Uf(i,j,1) = U(i,j,1)
         Vf(i,j,1) = V(i,j,1)
         Wf(i,j,1) = W(i,j,1)
         Uf(i,j,Kloc1) = U(i,j,Kloc)
         Vf(i,j,Kloc1) = V(i,j,Kloc)
         Wf(i,j,Kloc1) = W(i,j,Kloc)
       enddo
       enddo
     else
       ! second-order lagrange interpolation
       do k = 3,Kloc
       do j = 1,Nloc
       do i = 1,Mloc
         Int_factor1 = (sig(k)-sigc(k-1))*(sig(k)-sigc(k))/  &
             ((sigc(k-2)-sigc(k-1))*(sigc(k-2)-sigc(k)))
         Int_factor2 = (sig(k)-sigc(k-2))*(sig(k)-sigc(k))/  &
             ((sigc(k-1)-sigc(k-2))*(sigc(k-1)-sigc(k)))
         Int_factor3 = (sig(k)-sigc(k-2))*(sig(k)-sigc(k-1))/  &
             ((sigc(k)-sigc(k-2))*(sigc(k)-sigc(k-1)))
         Uf(i,j,k) = Int_factor1*U(i,j,k-2)+Int_factor2*U(i,j,k-1)+Int_factor3*U(i,j,k)
         Vf(i,j,k) = Int_factor1*V(i,j,k-2)+Int_factor2*V(i,j,k-1)+Int_factor3*V(i,j,k)
         Wf(i,j,k) = Int_factor1*W(i,j,k-2)+Int_factor2*W(i,j,k-1)+Int_factor3*W(i,j,k)
       enddo
       enddo
       enddo

       do j = 1,Nloc
       do i = 1,Mloc
         Int_factor1 = (sig(2)-sigc(2))*(sig(2)-sigc(3))/  &
             ((sigc(1)-sigc(2))*(sigc(1)-sigc(2)))
         Int_factor2 = (sig(2)-sigc(1))*(sig(2)-sigc(3))/  &
             ((sigc(2)-sigc(1))*(sigc(2)-sigc(3)))
         Int_factor3 = (sig(2)-sigc(1))*(sig(2)-sigc(2))/  &
             ((sigc(3)-sigc(1))*(sigc(3)-sigc(2)))
         Uf(i,j,2) = Int_factor1*U(i,j,1)+Int_factor2*U(i,j,2)+Int_factor3*U(i,j,3)
         Vf(i,j,2) = Int_factor1*V(i,j,1)+Int_factor2*V(i,j,2)+Int_factor3*V(i,j,3)
         Wf(i,j,2) = Int_factor1*W(i,j,1)+Int_factor2*W(i,j,2)+Int_factor3*W(i,j,3)
         Uf(i,j,1) = U(i,j,1)
         Vf(i,j,1) = V(i,j,1)
         Wf(i,j,1) = W(i,j,1)
         Uf(i,j,Kloc1) = U(i,j,Kloc)
         Vf(i,j,Kloc1) = V(i,j,Kloc)
         Wf(i,j,Kloc1) = W(i,j,Kloc)
       enddo
       enddo
     endif

     end subroutine interpolate_velocity_to_faces


     subroutine get_Omega(R1)
!-----------------------------------------------
!    Obtain vertical velocity in sigma-coord.
!    Called by 
!       eval_duvw
!    Last update: 30/08/2013, Gangfeng Ma
!-----------------------------------------------
     use global
     implicit none
     integer :: i,j,k
     real(SP), dimension(Mloc,Nloc), intent(in) :: R1
     real(SP), dimension(:,:,:), allocatable :: D3xL,D3xR,D3yL,D3yR

     ! reconstruct flux using new velocities
     call delxFun_2D(Eta,DelxEta)
     call delyFun_2D(Eta,DelyEta)
     call delxFun_3D(U,DelxU)
     call delxFun_3D(V,DelxV)
     call delyFun_3D(U,DelyU)
     call delyFun_3D(V,DelyV)
     call delxFun_3D(DU,DelxDU)
     call delxFun_3D(DV,DelxDV)
     call delyFun_3D(DU,DelyDU)
     call delyFun_3D(DV,DelyDV)

     call construct_2D_x(Eta,DelxEta,EtaxL,EtaxR)
     call construct_2D_y(Eta,DelyEta,EtayL,EtayR)
     call construct_3D_x(U,DelxU,UxL,UxR)
     call construct_3D_x(V,DelxV,VxL,VxR)
     call construct_3D_y(U,DelyU,UyL,UyR)
     call construct_3D_y(V,DelyV,VyL,VyR)
     call construct_3D_x(DU,DelxDU,DUxL,DUxR)
     call construct_3D_x(DV,DelxDV,DVxL,DVxR)
     call construct_3D_y(DU,DelyDU,DUyL,DUyR)
     call construct_3D_y(DV,DelyDV,DVyL,DVyR)

     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend1
       DxL(i,j) = EtaxL(i,j)+hfx(i,j)
       DxR(i,j) = EtaxR(i,j)+hfx(i,j)
       ExL(i,j,k) = DUxL(i,j,k)
       ExR(i,j,k) = DUxR(i,j,k)
     enddo
     enddo
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend1
     do i = Ibeg,Iend
       DyL(i,j) = EtayL(i,j)+hfy(i,j)
       DyR(i,j) = EtayR(i,j)+hfy(i,j) 
       EyL(i,j,k) = DVyL(i,j,k)
       EyR(i,j,k) = DVyR(i,j,k)
     enddo
     enddo
     enddo     

     allocate(D3xL(Mloc1,Nloc,Kloc))
     allocate(D3xR(Mloc1,Nloc,Kloc))
     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc1
       D3xL(i,j,k) = EtaxL(i,j)
       D3xR(i,j,k) = EtaxR(i,j)
     enddo
     enddo
     enddo

     allocate(D3yL(Mloc,Nloc1,Kloc))
     allocate(D3yR(Mloc,Nloc1,Kloc))
     do k = 1,Kloc
     do j = 1,Nloc1
     do i = 1,Mloc
       D3yL(i,j,k) = EtayL(i,j)
       D3yR(i,j,k) = EtayR(i,j)
     enddo
     enddo
     enddo

     call wave_speed

     call HLL(Mloc1,Nloc,Kloc,SxL,SxR,ExL,ExR,D3xL,D3xR,Ex)
     call HLL(Mloc,Nloc1,Kloc,SyL,SyR,EyL,EyR,D3yL,D3yR,Ey)

     ! left and right side
     if(n_west.eq.MPI_PROC_NULL) then
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_X0==1.or.Bc_X0==2) then
         Ex(Ibeg,j,k) = Zero
       elseif(Bc_X0==3) then
         Ex(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)
       endif
     enddo
     enddo
     endif

     if(n_east.eq.MPI_PROC_NULL) then
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_Xn==1.or.Bc_Xn==2) then
         Ex(Iend1,j,k) = Zero
       elseif(Bc_Xn==3) then
         Ex(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)
       endif
     enddo
     enddo
     endif

     ! front and back side  
     if(n_suth.eq.MPI_PROC_NULL) then
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Y0==1.or.Bc_Y0==2) then
         Ey(i,Jbeg,k) = Zero
       endif
     enddo
     enddo
     endif

     if(n_nrth.eq.MPI_PROC_NULL) then
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Yn==1.or.Bc_Yn==2) then
         Ey(i,Jend1,k) = Zero
       endif
     enddo
     enddo
     endif

     ! update Omega
     Omega = zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
     do k = Kbeg+1,Kend1
       Omega(i,j,k) = Omega(i,j,k-1)-dsig(k-1)*  &
          (R1(i,j)+(Ex(i+1,j,k-1)-Ex(i,j,k-1))/dx+(Ey(i,j+1,k-1)-Ey(i,j,k-1))/dy)                      
       if(Mask(i,j)==0) Omega(i,j,k) = zero
     enddo
     enddo
     enddo

     ! adjust omega to make omega=zero at free surface
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       if(abs(Omega(i,j,Kend1))>1.e-8) then
         do k = Kbeg+1,Kend1
           Omega(i,j,k) = Omega(i,j,k)-  &
              float(k-Kbeg)/float(Kend-Kbeg+1)*Omega(i,j,Kend1)                                        
         enddo
       endif
     enddo
     enddo

     deallocate(D3xL)
     deallocate(D3xR)
     deallocate(D3yL)
     deallocate(D3yR)

     return
     end subroutine get_Omega


     subroutine get_UVW
!------------------------------------------------
!    Obtain U,V,W
!    Called by
!       eval_duvw
!    Last update: 25/12/2010, Gangfeng Ma
!-----------------------------------------------
     use global
     implicit none
     integer :: i,j,k

     do j = Jbeg,Jend
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       U(i,j,k) = DU(i,j,k)/D(i,j)
       V(i,j,k) = DV(i,j,k)/D(i,j)
       W(i,j,k) = DW(i,j,k)/D(i,j)
     enddo
     enddo
     enddo

     ! collect data into ghost cells
     call vel_bc
     call phi_3D_exch(U)
     call phi_3D_exch(V)
     call phi_3D_exch(W)
     call phi_3D_exch(DU)
     call phi_3D_exch(DV)
     call phi_3D_exch(DW)

     end subroutine get_UVW


     subroutine fluxes
!------------------------------------------------
!    This subroutine is used to calculate fluxes 
!    at cell faces
!    Called by
!       main
!    Last update: 23/12/2010, Gangfeng Ma
!------------------------------------------------
     use global
     implicit none

     ! second order construction
     call delxyzFun
     call construction  

     ! calculate wave speed
     call wave_speed

     ! calculate fluxes at faces
     if(ADV_HLLC) then
       call fluxes_at_faces_HLLC
     else
       call fluxes_at_faces_HLL
     endif

     ! impose boundary conditions
     call flux_bc

     end subroutine fluxes


     subroutine construction
!------------------------------------------
!    Second-order construction
!    Called by 
!       fluxes
!    Last update: 04/01/2011, Gangfeng Ma
!-----------------------------------------
     use global
     implicit none
     integer :: i,j,k

     call construct_2D_x(Eta,DelxEta,EtaxL,EtaxR)
     call construct_3D_x(U,DelxU,UxL,UxR)
     call construct_3D_x(V,DelxV,VxL,VxR)
     call construct_3D_x(W,DelxW,WxL,WxR)
     call construct_3D_x(DU,DelxDU,DUxL,DUxR)
     call construct_3D_x(DV,DelxDV,DVxL,DVxR)
     call construct_3D_x(DW,DelxDW,DWxL,DWxR)

     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend1
       DxL(i,j) = EtaxL(i,j)+hfx(i,j)
       DxR(i,j) = EtaxR(i,j)+hfx(i,j)
       ExL(i,j,k) = DUxL(i,j,k)
       ExR(i,j,k) = DUxR(i,j,k)
       FxL(i,j,k) = DUxL(i,j,k)*UxL(i,j,k)+0.5*Grav*(EtaxL(i,j)*EtaxL(i,j)+2.0*EtaxL(i,j)*hfx(i,j))
       FxR(i,j,k) = DUxR(i,j,k)*UxR(i,j,k)+0.5*Grav*(EtaxR(i,j)*EtaxR(i,j)+2.0*EtaxR(i,j)*hfx(i,j))
       GxL(i,j,k) = DxL(i,j)*UxL(i,j,k)*VxL(i,j,k)
       GxR(i,j,k) = DxR(i,j)*UxR(i,j,k)*VxR(i,j,k)
       HxL(i,j,k) = DxL(i,j)*UxL(i,j,k)*WxL(i,j,k)
       HxR(i,j,k) = DxR(i,j)*UxR(i,j,k)*WxR(i,j,k)
     enddo
     enddo
     enddo

     call construct_2D_y(Eta,DelyEta,EtayL,EtayR)
     call construct_3D_y(U,DelyU,UyL,UyR)
     call construct_3D_y(V,DelyV,VyL,VyR)
     call construct_3D_y(W,DelyW,WyL,WyR)
     call construct_3D_y(DU,DelyDU,DUyL,DUyR)
     call construct_3D_y(DV,DelyDV,DVyL,DVyR)
     call construct_3D_y(DW,DelyDW,DWyL,DWyR)

     do k = Kbeg,Kend
     do j = Jbeg,Jend1
     do i = Ibeg,Iend
       DyL(i,j) = EtayL(i,j)+hfy(i,j)
       DyR(i,j) = EtayR(i,j)+hfy(i,j)
       EyL(i,j,k) = DVyL(i,j,k)
       EyR(i,j,k) = DVyR(i,j,k)
       FyL(i,j,k) = DyL(i,j)*UyL(i,j,k)*VyL(i,j,k)
       FyR(i,j,k) = DyR(i,j)*UyR(i,j,k)*VyR(i,j,k)
       GyL(i,j,k) = DVyL(i,j,k)*VyL(i,j,k)+0.5*Grav*(EtayL(i,j)*EtayL(i,j)+2.0*EtayL(i,j)*hfy(i,j))
       GyR(i,j,k) = DVyR(i,j,k)*VyR(i,j,k)+0.5*Grav*(EtayR(i,j)*EtayR(i,j)+2.0*EtayR(i,j)*hfy(i,j))
       HyL(i,j,k) = DyL(i,j)*VyL(i,j,k)*WyL(i,j,k)
       HyR(i,j,k) = DyR(i,j)*VyR(i,j,k)*WyR(i,j,k) 
     enddo
     enddo
     enddo

     call construct_3D_z(U,DelzU,UzL,UzR)
     call construct_3D_z(V,DelzV,VzL,VzR)
     call construct_3D_z(W,DelzW,WzL,WzR)

     end subroutine construction


     subroutine construct_2D_x(Vin,Din,OutL,OutR)
!-------------------------------------------------
!    Construct 2D variables in x-direction
!    Called by
!       construction
!    Last update: 04/01/2011, Gangfeng Ma
!------------------------------------------------
     use global, only: SP,Zero,dx,Mloc,Nloc,Mloc1, &
                       Ibeg,Iend,Jbeg,Jend,Iend1
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)   :: Vin,Din
     real(SP),intent(out),dimension(Mloc1,Nloc) :: OutL,OutR      
     integer :: i,j

     OutL = Zero
     OutR = Zero
     do i = Ibeg,Iend1
     do j = Jbeg,Jend
       OutL(i,j) = Vin(i-1,j)+0.5*dx*Din(i-1,j)
       OutR(i,j) = Vin(i,j)-0.5*dx*Din(i,j)
     enddo
     enddo

     end subroutine construct_2D_x


     subroutine construct_3D_x(Vin,Din,OutL,OutR)
!-------------------------------------------------
!    Construct 3D variables in x-direction
!    Called by 
!       construction 
!    Last update: 04/01/2011, Gangfeng Ma 
!------------------------------------------------
     use global, only: SP,Zero,dx,Mloc,Nloc,Kloc,Mloc1, &
                       Ibeg,Iend,Jbeg,Jend,Iend1,Kbeg,Kend
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)   :: Vin,Din
     real(SP),intent(out),dimension(Mloc1,Nloc,Kloc) :: OutL,OutR      
     integer :: i,j,k

     OutL = Zero
     OutR = Zero
     do i = Ibeg,Iend1
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       OutL(i,j,k) = Vin(i-1,j,k)+0.5*dx*Din(i-1,j,k)
       OutR(i,j,k) = Vin(i,j,k)-0.5*dx*Din(i,j,k)
     enddo
     enddo
     enddo

     end subroutine construct_3D_x


     subroutine construct_2D_y(Vin,Din,OutL,OutR)
!-------------------------------------------------
!    Construct 2D variables in y-direction
!    Called by
!       construction 
!    Last update: 04/01/2011, Gangfeng Ma 
!------------------------------------------------ 
     use global, only: SP,Zero,dy,Mloc,Nloc,Nloc1, &
                       Ibeg,Iend,Jbeg,Jend,Jend1
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)   :: Vin,Din
     real(SP),intent(out),dimension(Mloc,Nloc1) :: OutL,OutR      
     integer :: i,j

     OutL = Zero
     OutR = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend1
       OutL(i,j) = Vin(i,j-1)+0.5*dy*Din(i,j-1)
       OutR(i,j) = Vin(i,j)-0.5*dy*Din(i,j)
     enddo
     enddo

     end subroutine construct_2D_y


     subroutine construct_3D_y(Vin,Din,OutL,OutR)
!-------------------------------------------------
!    Construct 3D variables in y-direction 
!    Called by
!       construction 
!    Last update: 04/01/2011, Gangfeng Ma 
!------------------------------------------------
     use global, only: SP,Zero,dy,Mloc,Nloc,Kloc,Nloc1, &
                       Ibeg,Iend,Jbeg,Jend,Jend1,Kbeg,Kend
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)   :: Vin,Din
     real(SP),intent(out),dimension(Mloc,Nloc1,Kloc) :: OutL,OutR      
     integer :: i,j,k

     OutL = Zero
     OutR = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend1
     do k = Kbeg,Kend
       OutL(i,j,k) = Vin(i,j-1,k)+0.5*dy*Din(i,j-1,k)
       OutR(i,j,k) = Vin(i,j,k)-0.5*dy*Din(i,j,k)
     enddo
     enddo
     enddo

     end subroutine construct_3D_y


     subroutine construct_3D_z(Vin,Din,OutL,OutR)
!-------------------------------------------------
!    Construct 3D variables in z-direction
!    Called by
!       construction 
!    Last update: 04/01/2011, Gangfeng Ma
!------------------------------------------------
     use global, only: SP,Zero,dsig,Mloc,Nloc,Kloc,Kloc1, &
                       Ibeg,Iend,Jbeg,Jend,Kbeg,Kend,Kend1
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)   :: Vin,Din
     real(SP),intent(out),dimension(Mloc,Nloc,Kloc1) :: OutL,OutR      
     integer :: i,j,k

     OutL = Zero
     OutR = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
     do k = Kbeg,Kend1
       OutL(i,j,k) = Vin(i,j,k-1)+0.5*dsig(k-1)*Din(i,j,k-1)
       OutR(i,j,k) = Vin(i,j,k)-0.5*dsig(k)*Din(i,j,k)
     enddo
     enddo
     enddo

     end subroutine construct_3D_z


     subroutine delxyzFun 
!-------------------------------------------
!    Calculate variable derivatives 
!    Called by 
!       fluxes 
!    Last update: 04/01/2011, Gangfeng Ma
!------------------------------------------
     use global
     implicit none
     integer :: i
     
     call delxFun_2D(Eta,DelxEta)
     call delxFun_3D(U,DelxU)
     call delxFun_3D(V,DelxV)
     call delxFun_3D(W,DelxW)
     call delxFun_3D(DU,DelxDU)
     call delxFun_3D(DV,DelxDV)
     call delxFun_3D(DW,DelxDW)

     call delyFun_2D(Eta,DelyEta)
     call delyFun_3D(U,DelyU)
     call delyFun_3D(V,DelyV)
     call delyFun_3D(W,DelyW)
     call delyFun_3D(DU,DelyDU)
     call delyFun_3D(DV,DelyDV)
     call delyFun_3D(DW,DelyDW)

     call delzFun_3D(U,DelzU)
     call delzFun_3D(V,DelzV)
     call delzFun_3D(W,DelzW)

     end subroutine delxyzFun

     
     subroutine delxFun_2D(Din,Dout)
!-------------------------------------------
!    Second-order derivative in x
!    Called by
!       delxyzFun
!    Last update: 04/01/2011, Gangfeng Ma
!------------------------------------------
     use global, only: SP,Small,Zero,dx,Mloc,Nloc,Mask,Brks
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc) :: Dout
     real(SP) :: TMP1,TMP2,LIMITER
     integer :: i,j
    
     do i = 2,Mloc-1
     do j = 1,Nloc
       if(Mask(i,j)==0) then
         Dout(i,j) = Zero
       else
         TMP1 = (Din(i+1,j)-Din(i,j))/dx
         TMP2 = (Din(i,j)-Din(i-1,j))/dx

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j) = Zero
         else
           Dout(i,j) = LIMITER(TMP1,TMP2)
         endif
       endif
     enddo
     enddo

     do j = 1,Nloc
       Dout(1,j) = (Din(2,j)-Din(1,j))/dx
       Dout(Mloc,j) = (Din(Mloc,j)-Din(Mloc-1,j))/dx
     enddo  

     return
     end subroutine delxFun_2D


     subroutine delxFun_3D(Din,Dout)
!------------------------------------------
!    Second-order derivative in x
!    Called by 
!       delxyzFun
!    Last update: 04/01/2011, Gangfeng Ma
!------------------------------------------
     use global, only: SP,Small,Zero,dx,Mloc,Nloc,Kloc,Mask
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc,Kloc) :: Dout
     real(SP) :: TMP1,TMP2,LIMITER
     integer :: i,j,k

     do i = 2,Mloc-1
     do j = 1,Nloc
     do k = 1,Kloc
       if(Mask(i,j)==0) then
         Dout(i,j,k) = Zero
       else
         TMP1 = (Din(i+1,j,k)-Din(i,j,k))/dx
         TMP2 = (Din(i,j,k)-Din(i-1,j,k))/dx

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j,k) = Zero
         else
           Dout(i,j,k) = LIMITER(TMP1,TMP2)
         endif
       endif
     enddo
     enddo
     enddo

     do j = 1,Nloc
     do k = 1,Kloc
       Dout(1,j,k) = (Din(2,j,k)-Din(1,j,k))/dx
       Dout(Mloc,j,k) = (Din(Mloc,j,k)-Din(Mloc-1,j,k))/dx
     enddo
     enddo

     return
     end subroutine delxFun_3D


     subroutine delyFun_2D(Din,Dout)
!-----------------------------------------
!    Second-order derivative in y
!    Called by 
!       delxyzFun  
!    Last update: 04/01/2011, Gangfeng Ma  
!------------------------------------------ 
     use global, only: SP,Small,Zero,dy,Mloc,Nloc,Mask
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc) :: Dout
     real(SP) :: TMP1,TMP2,LIMITER
     integer :: i,j

     do i = 1,Mloc
     do j = 2,Nloc-1
       if(Mask(i,j)==0) then 
         Dout(i,j) = Zero
       else
         TMP1 = (Din(i,j+1)-Din(i,j))/dy
         TMP2 = (Din(i,j)-Din(i,j-1))/dy

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j) = Zero
         else
           Dout(i,j) = LIMITER(TMP1,TMP2)
         endif
       endif
     enddo
     enddo

     do i = 1,Mloc
       Dout(i,1) = (Din(i,2)-Din(i,1))/dy
       Dout(i,Nloc) = (Din(i,Nloc)-Din(i,Nloc-1))/dy
     enddo

     return
     end subroutine delyFun_2D


     subroutine delyFun_3D(Din,Dout)
!-------------------------------------------
!    Second-order derivative in y
!    Called by
!       delxyzFun 
!    Last update: 04/01/2011, Gangfeng Ma 
!-------------------------------------------
     use global, only: SP,Small,Zero,dy,Mloc,Nloc,Kloc,Mask
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc,Kloc) :: Dout
     real(SP) :: TMP1,TMP2,LIMITER
     integer :: i,j,k

     do i = 1,Mloc
     do j = 2,Nloc-1
     do k = 1,Kloc
       if(Mask(i,j)==0) then
         Dout(i,j,k) = Zero
       else
         TMP1 = (Din(i,j+1,k)-Din(i,j,k))/dy
         TMP2 = (Din(i,j,k)-Din(i,j-1,k))/dy

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j,k) = Zero
         else
           Dout(i,j,k) = LIMITER(TMP1,TMP2)
         endif
       endif
     enddo
     enddo
     enddo

     do i = 1,Mloc
     do k = 1,Kloc
       Dout(i,1,k) = (Din(i,2,k)-Din(i,1,k))/dy
       Dout(i,Nloc,k) = (Din(i,Nloc,k)-Din(i,Nloc-1,k))/dy
     enddo
     enddo

     return 
     end subroutine delyFun_3D

     
     subroutine delzFun_3D(Din,Dout)
!-------------------------------------------
!    Second-order derivative in z
!    Called by
!       delxyzFun
!    Last update: 04/01/2011, Gangfeng Ma
!-------------------------------------------
     use global, only: SP,Small,Zero,dsig,sigc,Mloc,Nloc,Kloc
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc,Kloc) :: Dout
     real(SP) :: TMP1,TMP2,LIMITER
     integer :: i,j,k

     do i = 1,Mloc
     do j = 1,Nloc
     do k = 2,Kloc-1
       TMP1 = (Din(i,j,k+1)-Din(i,j,k))/(sigc(k+1)-sigc(k))
       TMP2 = (Din(i,j,k)-Din(i,j,k-1))/(sigc(k)-sigc(k-1))

       if((abs(TMP1)+abs(TMP2))<Small) then
         Dout(i,j,k) = Zero
       else
         Dout(i,j,k) = LIMITER(TMP1,TMP2)
       endif
     enddo
     enddo
     enddo

     do i = 1,Mloc
     do j = 1,Nloc
       Dout(i,j,1) = (Din(i,j,2)-Din(i,j,1))/(0.5*(dsig(1)+dsig(2)))
       Dout(i,j,Kloc) = (Din(i,j,Kloc)-Din(i,j,Kloc-1))/(0.5*(dsig(Kloc-1)+dsig(Kloc)))
     enddo
     enddo

     return
     end subroutine delzFun_3D


     subroutine flux_bc
!--------------------------------------------
!    This is subroutine to provide boundary conditions
!    Called by
!       fluxes
!    Last update: 25/12/2010, Gangfeng Ma
!--------------------------------------------
     use global
     implicit none
     integer :: i,j,k

     ! left and right side
     if(n_west.eq.MPI_PROC_NULL) then
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_X0==1.or.Bc_X0==2) then
         Ex(Ibeg,j,k) = Zero
         Fx(Ibeg,j,k) = 0.5*Grav*(EtaxR(Ibeg,j)*EtaxR(Ibeg,j)+2.0*EtaxR(Ibeg,j)*hfx(Ibeg,j))
         Gx(Ibeg,j,k) = Zero
         Hx(Ibeg,j,k) = Zero
       elseif(Bc_X0==3) then
         Ex(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)
         Fx(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)*Uin_X0(j,k)+  &
                  0.5*Grav*(Ein_X0(j)*Ein_X0(j)+2.0*Ein_X0(j)*hfx(Ibeg,j))
         Gx(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)*Vin_X0(j,k)
         Hx(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)*Win_X0(j,k)
       endif
     enddo
     enddo
     endif

     if(n_east.eq.MPI_PROC_NULL) then
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_Xn==1.or.Bc_Xn==2) then
         Ex(Iend1,j,k) = Zero
         Fx(Iend1,j,k) = 0.5*Grav*(EtaxL(Iend1,j)*EtaxL(Iend1,j)+  &
                  2.0*EtaxL(Iend1,j)*hfx(Iend1,j))
         Gx(Iend1,j,k) = Zero
         Hx(Iend1,j,k) = Zero
       elseif(Bc_Xn==3) then
         Ex(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)
         Fx(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)*Uin_Xn(j,k)+  &
                 0.5*Grav*(Ein_Xn(j)*Ein_Xn(j)+2.0*Ein_Xn(j)*hfx(Iend1,j)) 
         Gx(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)*Vin_Xn(j,k)
         Hx(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)*Win_Xn(j,k)
       endif
     enddo
     enddo
     endif

     ! front and back side
     if(n_suth.eq.MPI_PROC_NULL) then
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Y0==1.or.Bc_Y0==2) then
         Ey(i,Jbeg,k) = Zero
         Fy(i,Jbeg,k) = Zero
         Gy(i,Jbeg,k) = 0.5*Grav*(EtayR(i,Jbeg)*EtayR(i,Jbeg)+  &
                2.0*EtayR(i,Jbeg)*hfy(i,Jbeg))
         Hy(i,Jbeg,k) = Zero
       endif
     enddo
     enddo
     endif

     if(n_nrth.eq.MPI_PROC_NULL) then
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Yn==1.or.Bc_Yn==2) then
         Ey(i,Jend1,k) = Zero
         Fy(i,Jend1,k) = Zero
         Gy(i,Jend1,k) = 0.5*Grav*(EtayL(i,Jend1)*EtayL(i,Jend1)+  &
                 2.0*EtayL(i,Jend1)*hfy(i,Jend1))
         Hy(i,Jend1,k) = Zero
       endif
     enddo
     enddo
     endif

     ! upper and bottom
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Fz(i,j,Kbeg) = Zero
       Gz(i,j,Kbeg) = Zero
       Hz(i,j,Kbeg) = Zero
       Fz(i,j,Kend1) = Zero
       Gz(i,j,Kend1) = Zero
       Hz(i,j,Kend1) = Zero
     enddo
     enddo

     do k = Kbeg,Kend
     do j = Jbeg-1,Jend+1
     do i = Ibeg-1,Iend+1
       if(Mask(i,j)==0) then
         Ex(i,j,k) = Zero
         if(i==Ibeg.and.n_west.eq.MPI_PROC_NULL) then
           Fx(i,j,k) = Zero
         else
           Fx(i,j,k) = 0.5*Grav*(EtaxL(i,j)*EtaxL(i,j)+  &
                2.0*EtaxL(i,j)*hfx(i,j))*Mask(i-1,j)
         endif
         Gx(i,j,k) = Zero
         Hx(i,j,k) = Zero

         Ex(i+1,j,k) = Zero
         if(i==Iend.and.n_east.eq.MPI_PROC_NULL) then
           Fx(i+1,j,k) = Zero
         else
           Fx(i+1,j,k) = 0.5*Grav*(EtaxR(i+1,j)*EtaxR(i+1,j)+  &
                  2.0*EtaxR(i+1,j)*hfx(i+1,j))*Mask(i+1,j)
         endif
         Gx(i+1,j,k) = Zero
         Hx(i+1,j,k) = Zero

         Ey(i,j,k) = Zero
         Fy(i,j,k) = Zero
         if(j==Jbeg.and.n_suth.eq.MPI_PROC_NULL) then
           Gy(i,j,k) = Zero
         else
           Gy(i,j,k) = 0.5*Grav*(EtayL(i,j)*EtayL(i,j)+  &
                2.0*EtayL(i,j)*hfy(i,j))*Mask(i,j-1)
         endif
         Hy(i,j,k) = Zero

         Ey(i,j+1,k) = Zero
         Fy(i,j+1,k) = Zero
         if(j==Jend.and.n_nrth.eq.MPI_PROC_NULL) then
           Gy(i,j+1,k) = Zero
         else
           Gy(i,j+1,k) = 0.5*Grav*(EtayR(i,j+1)*EtayR(i,j+1)+  &
                2.0*EtayR(i,j+1)*hfy(i,j+1))*Mask(i,j+1)
         endif
         Hy(i,j+1,k) = Zero
       endif
     enddo
     enddo
     enddo

     end subroutine flux_bc


     subroutine fluxes_at_faces_HLLC
!---------------------------------------------
!    Fluxes at cell faces estimated by HLLC approximation
!    Called by 
!       fluxes
!    Last update: 24/12/2010, Gangfeng Ma
!---------------------------------------------
     use global
     implicit none
     integer  :: i,j,k
     real(SP), dimension(:,:,:), allocatable :: D3xL,D3xR,D3yL,D3yR,D3xLS,D3xRS,  &
                 D3yLS,D3yRS,DUxLS,DUxRS,DVxLS,DVxRS,DWxLS,DWxRS,DUyLS,DUyRS,DVyLS, &
                 DVyRS,DWyLS,DWyRS

     ! temporary arrays
     allocate(D3xL(Mloc1,Nloc,Kloc))
     allocate(D3xR(Mloc1,Nloc,Kloc))
     allocate(D3xLS(Mloc1,Nloc,Kloc))
     allocate(D3xRS(Mloc1,Nloc,Kloc))
     allocate(DUxLS(Mloc1,Nloc,Kloc))
     allocate(DUxRS(Mloc1,Nloc,Kloc))
     allocate(DVxLS(Mloc1,Nloc,Kloc))
     allocate(DVxRS(Mloc1,Nloc,Kloc))
     allocate(DWxLS(Mloc1,Nloc,Kloc))
     allocate(DWxRS(Mloc1,Nloc,Kloc))
     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc1
       D3xL(i,j,k) = DxL(i,j)
       D3xR(i,j,k) = DxR(i,j)
       D3xLS(i,j,k) = DxL(i,j)*(SxL(i,j,k)-UxL(i,j,k)+Small)/(SxL(i,j,k)-SxS(i,j,k)+Small)
       D3xRS(i,j,k) = DxR(i,j)*(SxR(i,j,k)-UxR(i,j,k)+Small)/(SxR(i,j,k)-SxS(i,j,k)+Small)       
       DUxLS(i,j,k) = DxL(i,j)*(SxL(i,j,k)-UxL(i,j,k)+Small)/(SxL(i,j,k)-SxS(i,j,k)+Small)*SxS(i,j,k)
       DUxRS(i,j,k) = DxR(i,j)*(SxR(i,j,k)-UxR(i,j,k)+Small)/(SxR(i,j,k)-SxS(i,j,k)+Small)*SxS(i,j,k)
       DVxLS(i,j,k) = DxL(i,j)*(SxL(i,j,k)-UxL(i,j,k)+Small)/(SxL(i,j,k)-SxS(i,j,k)+Small)*VxL(i,j,k)
       DVxRS(i,j,k) = DxR(i,j)*(SxR(i,j,k)-UxR(i,j,k)+Small)/(SxR(i,j,k)-SxS(i,j,k)+Small)*VxR(i,j,k)
       DWxLS(i,j,k) = DxL(i,j)*(SxL(i,j,k)-UxL(i,j,k)+Small)/(SxL(i,j,k)-SxS(i,j,k)+Small)*WxL(i,j,k)
       DWxRS(i,j,k) = DxR(i,j)*(SxR(i,j,k)-UxR(i,j,k)+Small)/(SxR(i,j,k)-SxS(i,j,k)+Small)*WxR(i,j,k)
     enddo
     enddo
     enddo

     allocate(D3yL(Mloc,Nloc1,Kloc))
     allocate(D3yR(Mloc,Nloc1,Kloc))
     allocate(D3yLS(Mloc,Nloc1,Kloc))
     allocate(D3yRS(Mloc,Nloc1,Kloc))
     allocate(DUyLS(Mloc,Nloc1,Kloc))
     allocate(DUyRS(Mloc,Nloc1,Kloc))
     allocate(DVyLS(Mloc,Nloc1,Kloc))
     allocate(DVyRS(Mloc,Nloc1,Kloc))
     allocate(DWyLS(Mloc,Nloc1,Kloc))
     allocate(DWyRS(Mloc,Nloc1,Kloc))
     do k = 1,Kloc
     do j = 1,Nloc1
     do i = 1,Mloc
       D3yL(i,j,k) = DyL(i,j)
       D3yR(i,j,k) = DyR(i,j)
       D3yLS(i,j,k) = DyL(i,j)*(SyL(i,j,k)-VyL(i,j,k)+Small)/(SyL(i,j,k)-SyS(i,j,k)+Small)
       D3yRS(i,j,k) = DyR(i,j)*(SyR(i,j,k)-VyR(i,j,k)+Small)/(SyR(i,j,k)-SyS(i,j,k)+Small)
       DUyLS(i,j,k) = DyL(i,j)*(SyL(i,j,k)-VyL(i,j,k)+Small)/(SyL(i,j,k)-SyS(i,j,k)+Small)*UyL(i,j,k)
       DUyRS(i,j,k) = DyR(i,j)*(SyR(i,j,k)-VyR(i,j,k)+Small)/(SyR(i,j,k)-SyS(i,j,k)+Small)*UyR(i,j,k)
       DVyLS(i,j,k) = DyL(i,j)*(SyL(i,j,k)-VyL(i,j,k)+Small)/(SyL(i,j,k)-SyS(i,j,k)+Small)*SyS(i,j,k)
       DVyRS(i,j,k) = DyR(i,j)*(SyR(i,j,k)-VyR(i,j,k)+Small)/(SyR(i,j,k)-SyS(i,j,k)+Small)*SyS(i,j,k)
       DWyLS(i,j,k) = DyL(i,j)*(SyL(i,j,k)-VyL(i,j,k)+Small)/(SyL(i,j,k)-SyS(i,j,k)+Small)*WyL(i,j,k)
       DWyRS(i,j,k) = DyR(i,j)*(SyR(i,j,k)-VyR(i,j,k)+Small)/(SyR(i,j,k)-SyS(i,j,k)+Small)*WyR(i,j,k)
     enddo
     enddo
     enddo

     ! horizontal fluxes
     call HLLC(Mloc1,Nloc,Kloc,SxL,SxR,SxS,ExL,ExR,D3xL,D3xLS,D3xR,D3xRS,Ex)
     call HLLC(Mloc,Nloc1,Kloc,SyL,SyR,SyS,EyL,EyR,D3yL,D3yLS,D3yR,D3yRS,Ey)
     call HLLC(Mloc1,Nloc,Kloc,SxL,SxR,SxS,FxL,FxR,DUxL,DUxLS,DUxR,DUxRS,Fx)
     call HLLC(Mloc,Nloc1,Kloc,SyL,SyR,SyS,FyL,FyR,DUyL,DUyLS,DUyR,DUyRS,Fy)
     call HLLC(Mloc1,Nloc,Kloc,SxL,SxR,SxS,GxL,GxR,DVxL,DVxLS,DVxR,DVxRS,Gx)
     call HLLC(Mloc,Nloc1,Kloc,SyL,SyR,SyS,GyL,GyR,DVyL,DVyLS,DVyR,DVyRS,Gy)
     call HLLC(Mloc1,Nloc,Kloc,SxL,SxR,SxS,HxL,HxR,DWxL,DWxLS,DWxR,DWxRS,Hx)
     call HLLC(Mloc,Nloc1,Kloc,SyL,SyR,SyS,HyL,HyR,DWyL,DWyLS,DWyR,DWyRS,Hy)     

     ! vertical fluxes
     do k = Kbeg+1,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Fz(i,j,k) = 0.5*(Omega(i,j,k)*(UzL(i,j,k)+UzR(i,j,k))-abs(Omega(i,j,k))*(UzR(i,j,k)-UzL(i,j,k)))
       Gz(i,j,k) = 0.5*(Omega(i,j,k)*(VzL(i,j,k)+VzR(i,j,k))-abs(Omega(i,j,k))*(VzR(i,j,k)-VzL(i,j,k)))
       Hz(i,j,k) = 0.5*(Omega(i,j,k)*(WzL(i,j,k)+WzR(i,j,k))-abs(Omega(i,j,k))*(WzR(i,j,k)-WzL(i,j,k)))
!       Fz(i,j,k) = Omega(i,j,k)*(dsig(k)*U(i,j,k-1)+dsig(k-1)*U(i,j,k))/(dsig(k)+dsig(k-1))
!       Gz(i,j,k) = Omega(i,j,k)*(dsig(k)*V(i,j,k-1)+dsig(k-1)*V(i,j,k))/(dsig(k)+dsig(k-1))                        
!       Hz(i,j,k) = Omega(i,j,k)*(dsig(k)*W(i,j,k-1)+dsig(k-1)*W(i,j,k))/(dsig(k)+dsig(k-1))
     enddo
     enddo
     enddo

     deallocate(D3xL)
     deallocate(D3xR)
     deallocate(D3yL)
     deallocate(D3yR)
     deallocate(D3xLS)
     deallocate(D3xRS)
     deallocate(D3yLS)
     deallocate(D3yRS)
     deallocate(DUxLS)
     deallocate(DUxRS)
     deallocate(DUyLS)
     deallocate(DUyRS)
     deallocate(DVxLS)
     deallocate(DVxRS)
     deallocate(DVyLS)
     deallocate(DVyRS)
     deallocate(DWxLS)
     deallocate(DWxRS)
     deallocate(DWyLS)
     deallocate(DWyRS)

     return
     end subroutine fluxes_at_faces_HLLC


     subroutine fluxes_at_faces_HLL
!---------------------------------------------
!    Fluxes at cell faces estimated by HLL approximation
!    Called by 
!       fluxes
!    Last update: 24/12/2010, Gangfeng Ma
!---------------------------------------------
     use global
     implicit none
     integer  :: i,j,k
     real(SP), dimension(:,:,:), allocatable :: D3xL,D3xR,D3yL,D3yR

     ! temporary arrays
     allocate(D3xL(Mloc1,Nloc,Kloc))
     allocate(D3xR(Mloc1,Nloc,Kloc))
     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc1
       D3xL(i,j,k) = EtaxL(i,j)
       D3xR(i,j,k) = EtaxR(i,j)
     enddo
     enddo
     enddo

     allocate(D3yL(Mloc,Nloc1,Kloc))
     allocate(D3yR(Mloc,Nloc1,Kloc))
     do k = 1,Kloc
     do j = 1,Nloc1
     do i = 1,Mloc
       D3yL(i,j,k) = EtayL(i,j)
       D3yR(i,j,k) = EtayR(i,j)
     enddo
     enddo
     enddo

     ! horizontal fluxes
     call HLL(Mloc1,Nloc,Kloc,SxL,SxR,ExL,ExR,D3xL,D3xR,Ex)
     call HLL(Mloc,Nloc1,Kloc,SyL,SyR,EyL,EyR,D3yL,D3yR,Ey)
     call HLL(Mloc1,Nloc,Kloc,SxL,SxR,FxL,FxR,DUxL,DUxR,Fx)
     call HLL(Mloc,Nloc1,Kloc,SyL,SyR,FyL,FyR,DUyL,DUyR,Fy)
     call HLL(Mloc1,Nloc,Kloc,SxL,SxR,GxL,GxR,DVxL,DVxR,Gx)
     call HLL(Mloc,Nloc1,Kloc,SyL,SyR,GyL,GyR,DVyL,DVyR,Gy)
     call HLL(Mloc1,Nloc,Kloc,SxL,SxR,HxL,HxR,DWxL,DWxR,Hx)
     call HLL(Mloc,Nloc1,Kloc,SyL,SyR,HyL,HyR,DWyL,DWyR,Hy)     

     ! vertical fluxes
     do k = Kbeg+1,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Fz(i,j,k) = 0.5*(Omega(i,j,k)*(UzL(i,j,k)+UzR(i,j,k))-abs(Omega(i,j,k))*(UzR(i,j,k)-UzL(i,j,k)))
       Gz(i,j,k) = 0.5*(Omega(i,j,k)*(VzL(i,j,k)+VzR(i,j,k))-abs(Omega(i,j,k))*(VzR(i,j,k)-VzL(i,j,k)))
       Hz(i,j,k) = 0.5*(Omega(i,j,k)*(WzL(i,j,k)+WzR(i,j,k))-abs(Omega(i,j,k))*(WzR(i,j,k)-WzL(i,j,k)))
!       Fz(i,j,k) = Omega(i,j,k)*(dsig(k)*U(i,j,k-1)+dsig(k-1)*U(i,j,k))/(dsig(k)+dsig(k-1))
!       Gz(i,j,k) = Omega(i,j,k)*(dsig(k)*V(i,j,k-1)+dsig(k-1)*V(i,j,k))/(dsig(k)+dsig(k-1))
!       Hz(i,j,k) = Omega(i,j,k)*(dsig(k)*W(i,j,k-1)+dsig(k-1)*W(i,j,k))/(dsig(k)+dsig(k-1))
     enddo
     enddo
     enddo

     deallocate(D3xL)
     deallocate(D3xR)
     deallocate(D3yL)
     deallocate(D3yR)

     return
     end subroutine fluxes_at_faces_HLL


     subroutine HLL(M,N,L,SL,SR,FL,FR,UL,UR,FOUT)
!----------------------------------------------
!    HLLC reconstruction 
!    Called by
!       fluxes_at_faces_HLL
!    Last update: 24/12/2010, Gangfeng Ma
!---------------------------------------------
     use global, only: SP,ZERO,SMALL
     implicit none
     INTEGER,INTENT(IN)::M,N,L
     REAL(SP),INTENT(IN),DIMENSION(M,N,L)::SL,SR,FL,FR,UL,UR
     REAL(SP),INTENT(OUT),DIMENSION(M,N,L)::FOUT
     INTEGER :: I,J,K

     DO K = 1,L
     DO J = 1,N
     DO I = 1,M
       IF(SL(I,J,K)>=ZERO) THEN
         FOUT(I,J,K) = FL(I,J,K)
       ELSEIF(SR(I,J,K)<=ZERO) THEN
         FOUT(I,J,K) = FR(I,J,K)
       ELSE
         FOUT(I,J,K) = SR(I,J,K)*FL(I,J,K)-SL(I,J,K)*FR(I,J,K)+  &
               SL(I,J,K)*SR(I,J,K)*(UR(I,J,K)-UL(I,J,K))
         IF((ABS(SR(I,J,K)-SL(I,J,K)))<SMALL)THEN
           FOUT(I,J,K) = FOUT(I,J,K)/SMALL
         ELSE
           FOUT(I,J,K) = FOUT(I,J,K)/(SR(I,J,K)-SL(I,J,K))
         ENDIF
       ENDIF
     ENDDO
     ENDDO
     ENDDO

     return
     end subroutine HLL

   
     subroutine HLLC(M,N,L,SL,SR,SS,FL,FR,UL,ULS,UR,URS,FOUT)
!----------------------------------------------
!    HLLC reconstruction 
!    Called by
!       fluxes_at_faces_HLLC
!    Last update: 24/12/2010, Gangfeng Ma
!---------------------------------------------
     use global, only: SP,ZERO,SMALL
     implicit none
     INTEGER,INTENT(IN)::M,N,L
     REAL(SP),INTENT(IN),DIMENSION(M,N,L)::SL,SR,SS,FL,FR,UL,ULS,UR,URS
     REAL(SP),INTENT(OUT),DIMENSION(M,N,L)::FOUT
     INTEGER :: I,J,K

     DO K = 1,L
     DO J = 1,N
     DO I = 1,M
       IF(SL(I,J,K)>=ZERO) THEN
         FOUT(I,J,K) = FL(I,J,K)
       ELSEIF(SR(I,J,K)<=ZERO) THEN
         FOUT(I,J,K) = FR(I,J,K)
       ELSEIF(SS(I,J,K)>=ZERO) THEN
         FOUT(I,J,K) = FL(I,J,K)+SL(I,J,K)*(ULS(I,J,K)-UL(I,J,K))
       ELSE
         FOUT(I,J,K) = FR(I,J,K)+SR(I,J,K)*(URS(I,J,K)-UR(I,J,K))
       ENDIF
     ENDDO
     ENDDO
     ENDDO

     return
     end subroutine HLLC


     subroutine wave_speed
!----------------------------------------------
!    This subroutine is used to calculate wave speeds
!    Called by
!       fluxes
!    Last update: 24/12/2010, Gangfeng Ma
!    Last update: 12/04/2011, Gangfeng Ma, wetting-drying
!-----------------------------------------------
     use global, only: SP,Ibeg,Iend,Iend1,Jbeg,Jend,Jend1,Kbeg,Kend, &
                       DxL,DxR,DyL,DyR,UxL,UxR,VyL,VyR, &
                       SxL,SxR,SxS,SyL,SyR,SyS,Grav,Mask
     implicit none
     integer  :: i,j,k
     real(SP) :: SQR_PHI_L,SQR_PHI_R,SQR_PHI_S,U_S
     
     ! x-faces
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend1
       if(Mask(i-1,j)==1.and.Mask(i,j)==1) then
         SQR_PHI_L = sqrt(Grav*abs(DxL(i,j)))
         SQR_PHI_R = sqrt(Grav*abs(DxR(i,j)))
         SQR_PHI_S = 0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(UxL(i,j,k)-UxR(i,j,k))
         U_S = 0.5*(UxL(i,j,k)+UxR(i,j,k))+SQR_PHI_L-SQR_PHI_R
         SxL(i,j,k) = min(UxL(i,j,k)-SQR_PHI_L,U_S-SQR_PHI_S)
         SxR(i,j,k) = max(UxR(i,j,k)+SQR_PHI_R,U_S+SQR_PHI_S)
         SxS(i,j,k) = U_S
       elseif(Mask(i-1,j)==0.and.Mask(i,j)==1) then
         ! left-side dry case
         SQR_PHI_R = sqrt(Grav*abs(DxR(i,j)))
         SxL(i,j,k) = UxR(i,j,k)-2.0*SQR_PHI_R
         SxR(i,j,k) = UxR(i,j,k)+SQR_PHI_R
         SxS(i,j,k) = SxL(i,j,k)
       elseif(Mask(i-1,j)==1.and.Mask(i,j)==0) then
         ! right-side dry case
         SQR_PHI_L = sqrt(Grav*abs(DxL(i,j)))
         SxL(i,j,k) = UxL(i,j,k)-SQR_PHI_L
         SxR(i,j,k) = UxL(i,j,k)+2.0*SQR_PHI_L
         SxS(i,j,k) = SxR(i,j,k)
       endif
     enddo
     enddo
     enddo

     ! y-faces
     do k = Kbeg,Kend
     do j = Jbeg,Jend1
     do i = Ibeg,Iend
       if(Mask(i,j-1)==1.and.Mask(i,j)==1) then
         SQR_PHI_L = sqrt(Grav*abs(DyL(i,j)))
         SQR_PHI_R = sqrt(Grav*abs(DyR(i,j)))
         SQR_PHI_S = 0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(VyL(i,j,k)-VyR(i,j,k))
         U_S = 0.5*(VyL(i,j,k)+VyR(i,j,k))+SQR_PHI_L-SQR_PHI_R
         SyL(i,j,k) = min(VyL(i,j,k)-SQR_PHI_L,U_S-SQR_PHI_S)
         SyR(i,j,k) = max(VyR(i,j,k)+SQR_PHI_R,U_S+SQR_PHI_S)
         SyS(i,j,k) = U_S
       elseif(Mask(i,j-1)==0.and.Mask(i,j)==1) then
         ! left-side dry case
         SQR_PHI_R = sqrt(Grav*abs(DyR(i,j)))
         SyL(i,j,k) = VyR(i,j,k)-2.0*SQR_PHI_R
         SyR(i,j,k) = VyR(i,j,k)+SQR_PHI_R
         SyS(i,j,k) = SyL(i,j,k)
       elseif(Mask(i,j-1)==1.and.Mask(i,j)==0) then
         ! right-side dry case
         SQR_PHI_L = sqrt(Grav*abs(DyL(i,j)))
         SyL(i,j,k) = VyL(i,j,k)-SQR_PHI_L
         SyR(i,j,k) = VyL(i,j,k)+2.0*SQR_PHI_L
         SyS(i,j,k) = SyR(i,j,k)
       endif
     enddo
     enddo
     enddo

     end subroutine wave_speed


     FUNCTION LIMITER(A,B)
     use global, only: SP,Zero,One,Small
     IMPLICIT NONE
     REAL(SP),INTENT(IN) :: A,B
     REAL(SP) :: LIMITER

!     ! minmod limiter
!     LIMITER=max(Zero,min(A,B))

!     ! van Leer limiter
     LIMITER=(A*ABS(B)+ABS(A)*B)/(ABS(A)+ABS(B))

!     ! superbee limiter
!     LIMITER=SIGN(One,B)*MAX(Zero,MIN(2.0*ABS(B),SIGN(One,B)*A),  &
!          MIN(ABS(B),2.0*SIGN(One,B)*A))

     RETURN
     END FUNCTION LIMITER


     subroutine source_terms
!------------------------------------------------
!    This subroutine is used to evaluate source
!    Called by
!       main
!    Last update: 23/12/2010, Gangfeng Ma
!------------------------------------------------
     use global
     implicit none
     integer :: i,j,k,Iter,nn,ndir,nfreq,nk
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Source_Area,myvar, &
                 WnumX,WnumY,Phs_lag,dfreq,ddir,Angle,tmp1,tmp2,tmp3,tmp4,Umag,Ubar,Vbar
     real(SP) :: Ytrough,Mod1,Zup,Zlow,Zmid,Xstart,Zero1,cnoidal_cn,cnoidal_ck,Atmp

     ! internal wavemaker for linear wave
     if(WaveMaker(1:7)=='INT_LIN') then
       ! Find wave number for linear wave (Newton-Ralphson Method)
       Segma = 2.0*pi/Per_Wave
       Celerity = sqrt(Grav*Dep_Wave)
       Wave_Length = Celerity*Per_Wave
       Wave_Number = 2.0*pi/Wave_Length
     
       Iter = 0
 55    Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
       if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 65
       Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &
           Grav*tanh(Wave_Number*Dep_Wave) 
       Wave_Number = Wave_Number-Fk/Fkdif
       Iter = Iter+1
       goto 55
 65    continue
       Wave_Length = 2.0*pi/Wave_Number
       Celerity = Wave_Length/Per_Wave
       WnumX = Wave_Number*cos(Theta_Wave*pi/180.)
       WnumY = Wave_Number*sin(Theta_Wave*pi/180.)       

       Source_Area = 0.0    
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East) then
           Source_Area = Source_Area+dx*D(i,Jbeg+1)
         endif
       enddo

       call MPI_ALLREDUCE(Source_Area,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)
       Source_Area = myvar/float(PY)

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
           Phs_lag = (j-Jbeg)*dy*WnumY
           SourceC(i,j) = Celerity*Amp_Wave/Source_Area*cos(pi/2-Segma*time+Phs_lag)
         endif
       enddo
       enddo
     endif

     ! internal wavemaker for random waves 
     if(WaveMaker(1:7)=='INT_SPC') then
       Source_Area = 0.0
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East) then
           Source_Area = Source_Area+dx*D(i,Jbeg+1)
         endif
       enddo

       call MPI_ALLREDUCE(Source_Area,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)  
       Source_Area = myvar/float(PY)

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         SourceC(i,j) = Zero
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
              yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
           do nfreq = 1,NumFreq
           do ndir = 1,NumDir
             Per_Wave = 1.0/Freq(nfreq)
             Segma = 2.0*pi/Per_Wave
             Celerity = sqrt(Grav*Dep_Wave)
             Wave_Length = Celerity*Per_Wave
             Wave_Number = 2.0*pi/Wave_Length
       
             Iter = 0
   75        Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
             if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 85
             Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  & 
               Grav*tanh(Wave_Number*Dep_Wave)
             Wave_Number = Wave_Number-Fk/Fkdif
             Iter = Iter+1
             goto 75
   85        continue
             Wave_Length = 2.0*pi/Wave_Number
             Celerity = Wave_Length/Per_Wave

             ! adjust wave direction for periodic bc
             Angle = Dire(ndir)*pi/180.
             if(Angle>zero) then
               tmp3 = zero
               tmp1 = Wave_Number
               nk = 0
               do while (tmp3<Angle)
                 nk = nk+1
                 tmp2 = nk*2.0*pi/(Nglob*dy)
                 if(tmp2>=tmp1) then
                   tmp3 = 0.5*pi-small
                 else
                   tmp3 = asin(tmp2/tmp1)
                 endif
               enddo

               ! judge between nk-1 and nk which is closer                                          
               tmp4 = asin((nk-1)*2.0*pi/(Nglob*dy)/tmp1)
               if(abs(tmp4-Angle)<abs(Angle-tmp3)) then
                 Angle = tmp4
               else
                 Angle = tmp3
               endif
             else
               tmp3 = zero
               tmp1 = Wave_Number
               nk = 0
               do while (tmp3>Angle)
                 nk = nk+1
                tmp2 = nk*2.0*pi/(Nglob*dy)
                 if(tmp2>=tmp1) then
                   tmp3 = -0.5*pi+small
                 else
                   tmp3 = -asin(tmp2/tmp1)
                 endif
               enddo

               ! judge between nk-1 and nk which is closer                                          
               tmp4= asin((nk-1)*2.0*pi/(Nglob*dy)/tmp1)
               if(abs(tmp4-Angle)<abs(Angle-tmp3)) then
                 Angle = tmp4
               else
                 Angle = tmp3
               endif
             endif

             WnumX = Wave_Number*cos(Angle)
             WnumY = Wave_Number*sin(Angle)

             ! calculate root-mean-squre wave height for each component
             if(nfreq==1) then
               dfreq = Freq(2)-Freq(1)
             elseif(nfreq==NumFreq) then
               dfreq = Freq(NumFreq)-Freq(NumFreq-1)
             else
               dfreq = 0.5*(Freq(nfreq+1)-Freq(nfreq-1))
             endif
             dfreq = abs(dfreq)

             if(ndir==1) then
               ddir = Dire(2)-Dire(1)
             elseif(ndir==NumDir) then
               ddir = Dire(NumDir)-Dire(NumDir-1)
             else
               ddir = 0.5*(Dire(ndir+1)-Dire(ndir-1))
             endif
             ddir = abs(ddir)
         
             Amp_Wave = 2.0*sqrt(2.0*Wave_Spc2d(ndir,nfreq)*ddir*dfreq)

             Phs_lag = (dy/2.0+(j-Jbeg)*dy)*WnumY
             SourceC(i,j) = SourceC(i,j)+Celerity*Amp_Wave/Source_Area*  &
                  cos(pi/2-Segma*time+Phs_lag+Random_Phs(ndir,nfreq)) 
           enddo
           enddo
         endif
       enddo
       enddo
     endif

     ! internal wavemaker for jonswap spectrum
     if((WaveMaker(1:7)=='INT_JON').or.(WaveMaker(1:7)=='INT_TMA')) then
       Source_Area = 0.0
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East) then
           Source_Area = Source_Area+dx*D(i,Jbeg+1)
         endif
       enddo

       call MPI_ALLREDUCE(Source_Area,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)  
       Source_Area = myvar/float(PY)

       dfreq = (Freq_Max-Freq_Min)/float(NumFreq)

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         SourceC(i,j) = Zero
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
              yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
           do nfreq = 1,NumFreq
             Per_Wave = 1.0/Freq(nfreq)
             Segma = 2.0*pi/Per_Wave
             Celerity = sqrt(Grav*Dep_Wave)
             Wave_Length = Celerity*Per_Wave
             Wave_Number = 2.0*pi/Wave_Length
       
             Iter = 0
 76          Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
             if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 86
             Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  & 
               Grav*tanh(Wave_Number*Dep_Wave)
             Wave_Number = Wave_Number-Fk/Fkdif
             Iter = Iter+1
             goto 76
 86          continue
             Wave_Length = 2.0*pi/Wave_Number
             Celerity = Wave_Length/Per_Wave

             ! root-mean-square wave height
             Amp_Wave = 2.0*sqrt(2.0*Jon_Spc(nfreq)*DFreq)

             SourceC(i,j) = SourceC(i,j)+Celerity*Amp_Wave/Source_Area*  &
                cos(pi/2-Segma*time+RanPhs(nfreq))
           enddo
         endif
       enddo
       enddo
     endif


     ! internal wavemaker for cnoidal wave
     if(WaveMaker(1:7)=='INT_CON') then
       call cnoidal(Amp_Wave,Dep_Wave,Per_Wave,Wave_Length,Celerity,Ytrough,Mod1)

       ! wave number
       Wave_Number = 2.0*pi/Wave_Length

!# if defined(1)
!       if(myid.eq.0) write(*,*) 'Mod=',Mod1,'Ytrough=',Ytrough, &
!            'Wave_Number=',wave_number
!# endif      

       ! find zero start
       Zup = 1.0
       Zlow = 0.0
       Zmid= (Zup+Zlow)/2.0
       nn = 0
 200   nn = nn+1
       Zero1 = Ytrough+Amp_Wave*cnoidal_cn(Zmid*0.5*cnoidal_ck(Mod1),Mod1)**2                            

       if(abs(Zero1)<=1.0e-6) goto 210
       if(nn>1000) then
         write(*,*)'too many iterations; stop'
         stop
       endif
       if(Zero1<0.0) then
         Zup = Zmid
         Zmid = (Zup+Zlow)/2.0
         goto 200
       else
         Zlow = Zmid
         Zmid = (Zup+Zlow)/2.0
         goto 200
       endif
 210   continue
       Xstart = Zmid

       Source_Area = 0.0
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East) then
           Source_Area = Source_Area+dx*D(i,Jbeg+1)
         endif
       enddo

       call MPI_ALLREDUCE(Source_Area,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier) 
       Source_Area = myvar/float(PY)

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
           SourceC(i,j) = 2.0*Celerity/Source_Area*(Ytrough+Amp_Wave*cnoidal_cn(  &
               Xstart*0.5*cnoidal_ck(Mod1)+2.0*cnoidal_ck(Mod1)*(-TIME/Per_Wave),Mod1)**2)
         endif
       enddo
       enddo
     endif

     ! internal wavemaker for solitary wave
     if(WaveMaker(1:7)=='INT_SOL') then
       Celerity = sqrt(Grav*Dep_Wave*(1.0+Amp_Wave/Dep_Wave))
       Atmp = sqrt(0.75*Amp_Wave/Dep_Wave**3)
       Xstart = 4.0*Dep_Wave/sqrt(Amp_Wave/Dep_Wave)
       
       Source_Area = 0.0
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East) then
           Source_Area = Source_Area+dx*D(i,Jbeg+1)
         endif
       enddo

       call MPI_ALLREDUCE(Source_Area,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier) 
       Source_Area = myvar/float(PY)

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
           SourceC(i,j) = 2.0*Celerity/Source_Area*  &
               Amp_Wave/cosh(Atmp*(Xstart-Celerity*TIME))**2
         endif
       enddo
       enddo
     endif

     ! source terms for momentum eqs.
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       SourceX(i,j) = Grav*Eta(i,j)*DelxH(i,j)*Mask(i,j)
       SourceY(i,j) = Grav*Eta(i,j)*DelyH(i,j)*Mask(i,j)
     enddo
     enddo

     end subroutine source_terms
 

     subroutine wave_average
!---------------------------------------------------
!    Estimate wave averaged quantities
!    Called by                        
!       main 
!    Last update: 13/01/2012, Gangfeng Ma
!--------------------------------------------------
     use global
     implicit none
     real(SP), dimension(:,:), allocatable :: U_Dep_Ave,V_Dep_Ave,S_Dep_Ave
     integer :: i,j,k,n
     real(SP) :: Tmp,Tmp_0,Dz,Zk,U_zk,V_zk,W_zk,S_zk,Zbeg,Zend,Zn,Zn1,Sinterp

     allocate(U_Dep_Ave(Mloc,Nloc))
     allocate(V_Dep_Ave(Mloc,Nloc))

     if(TIME>Wave_Ave_Start.and.TIME<=Wave_Ave_End) then
       ! Lagrangian mean velocity
       do k = Kbeg,Kend
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         Lag_Umean(i,j,k) = Lag_Umean(i,j,k)+U(i,j,k)*dt/(Wave_Ave_End-Wave_Ave_Start)
         Lag_Vmean(i,j,k) = Lag_Vmean(i,j,k)+V(i,j,k)*dt/(Wave_Ave_End-Wave_Ave_Start)
         Lag_Wmean(i,j,k) = Lag_Wmean(i,j,k)+W(i,j,k)*dt/(Wave_Ave_End-Wave_Ave_Start)
       enddo
       enddo
       enddo

       ! Eulerian mean velocity
       do k = Kbeg,Kend
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         Dz = Hc(i,j)/float(Kglob)
         Zk = (k-Kbeg)*Dz+Dz/2.0

         U_zk = Zero; V_zk = Zero; W_zk = Zero
         Zbeg = sigc(Kbeg)*(Eta(i,j)+Hc(i,j))
         Zend = sigc(Kend)*(Eta(i,j)+Hc(i,j))
         if(Zk<=Zbeg) then
           U_zk = U(i,j,Kbeg)
           V_zk = V(i,j,Kbeg)
           W_zk = W(i,j,Kbeg)
         elseif(Zk>=Zend) then
           U_zk = U(i,j,Kend)
           V_zk = V(i,j,Kend)
           W_zk = W(i,j,Kend)
         else
           do n = Kbeg,Kend-1
             Zn = sigc(n)*(Eta(i,j)+Hc(i,j)) 
             Zn1 = sigc(n+1)*(Eta(i,j)+Hc(i,j))
             if(Zk>=Zn.and.Zk<Zn1) then
               Sinterp = (Zk-Zn)/(Zn1-Zn)
               U_zk = U(i,j,n)*(1.0-Sinterp)+U(i,j,n+1)*Sinterp
               V_zk = V(i,j,n)*(1.0-Sinterp)+V(i,j,n+1)*Sinterp
               W_zk = W(i,j,n)*(1.0-Sinterp)+W(i,j,n+1)*Sinterp
             endif 
           enddo
         endif
         Euler_Umean(i,j,k) = Euler_Umean(i,j,k)+U_zk*dt/(Wave_Ave_End-Wave_Ave_Start)
         Euler_Vmean(i,j,k) = Euler_Vmean(i,j,k)+V_Zk*dt/(Wave_Ave_End-Wave_Ave_Start)
         Euler_Wmean(i,j,k) = Euler_Wmean(i,j,k)+W_Zk*dt/(Wave_Ave_End-Wave_Ave_Start)
       enddo
       enddo
       enddo         

       ! depth-averaged velocity
       U_Dep_Ave = Zero
       V_Dep_Ave = Zero
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         do k = Kbeg,Kend
           U_Dep_Ave(i,j) = U_Dep_Ave(i,j)+U(i,j,k)/float(Kend-Kbeg+1)
           V_Dep_Ave(i,j) = V_Dep_Ave(i,j)+V(i,j,k)/float(Kend-Kbeg+1)
         enddo
       enddo
       enddo

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         Setup(i,j) = Setup(i,j)+Eta(i,j)*dt/(Wave_Ave_End-Wave_Ave_Start)
         Umean(i,j) = Umean(i,j)+U_Dep_Ave(i,j)*dt/(Wave_Ave_End-Wave_Ave_Start)
         Vmean(i,j) = Vmean(i,j)+V_Dep_Ave(i,j)*dt/(Wave_Ave_End-Wave_Ave_Start)

         if(Eta(i,j)>Emax(i,j)) Emax(i,j) = Eta(i,j)
         if(Eta(i,j)<Emin(i,j)) Emin(i,j) = Eta(i,j)
         
         Tmp = Eta(i,j)
         Tmp_0 = Eta0(i,j)
         if(Tmp>Tmp_0.and.Tmp*Tmp_0<=Zero) then
           Num_Zero_Up(i,j) = Num_Zero_Up(i,j)+1
           if(Num_Zero_Up(i,j)>=2) then
             if(WaveheightID==1) then  ! Average wave height
               WaveHeight(i,j) = WaveHeight(i,j)+Emax(i,j)-Emin(i,j)
             elseif(WaveheightID==2) then  ! RMS wave height
               WaveHeight(i,j) = WaveHeight(i,j)+(Emax(i,j)-Emin(i,j))**2
             endif
           endif

           ! reset Emax and Emin to find next wave
           Emax(i,j) = -1000.
           Emin(i,j) = 1000.
         endif  
       enddo
       enddo
     endif

     deallocate(U_Dep_Ave)
     deallocate(V_Dep_Ave)

     end subroutine wave_average

     
     subroutine print_wh_setup
!---------------------------------------------------
!    Estimate wave averaged quantities
!    Called by
!       main 
!    Last update: 13/01/2012, Gangfeng Ma
!--------------------------------------------------
     use global
     implicit none
     integer :: i,j
     character(len=80) :: FDIR,file

     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Num_Zero_Up(i,j)>=2) then
         if(WaveheightID==1) then
           WaveHeight(i,j) = WaveHeight(i,j)/float(Num_Zero_Up(i,j)-1)
         elseif(WaveheightID==2) then
           WaveHeight(i,j) = sqrt(WaveHeight(i,j))/float(Num_Zero_Up(i,j)-1)
         endif
       else
         WaveHeight(i,j) = Zero
       endif
     enddo
     enddo

     ! results directory
     FDIR = TRIM(RESULT_FOLDER)

     file = TRIM(FDIR)//'waveheight'
     call putfile2D(file,WaveHeight)

     file = TRIM(FDIR)//'setup'
     call putfile2D(file,Setup)
     
     file = TRIM(FDIR)//'umean'
     call putfile2D(file,Umean)

     file = TRIM(FDIR)//'vmean'
     call putfile2D(file,Vmean)


     file = TRIM(FDIR)//'lag_umean'
     call putfile3D(file,Lag_Umean)

     file = TRIM(FDIR)//'lag_vmean'
     call putfile3D(file,Lag_Vmean)

     file = TRIM(FDIR)//'lag_wmean'
     call putfile3D(file,Lag_Wmean)

     file = TRIM(FDIR)//'euler_umean'
     call putfile3D(file,Euler_Umean)

     file = TRIM(FDIR)//'euler_vmean'
     call putfile3D(file,Euler_Vmean)

     file = TRIM(FDIR)//'euler_wmean'
     call putfile3D(file,Euler_Wmean)


     end subroutine print_wh_setup


     subroutine statistics
!---------------------------------------------------
!    This subroutine is used to show statistics
!    Called by
!       main
!    Last update: 23/12/2010, Gangfeng Ma
!--------------------------------------------------
     use global
     implicit none
     real(SP) :: MassVolume,CellMass,Energy,MaxEta,MinEta,MaxU, &
                 MaxV,MaxW,MaxS,MinS
     integer :: i,j,k
     real(SP) :: myvar

     ! Vol = sum(D*dx*dy)
     ! Energy = sum(m*g*h+0.5*m*u^2), reference is at z = 0
     MassVolume = Zero
     Energy = Zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       MassVolume = MassVolume+D(i,j)*dx*dy
       do k = Kbeg,Kend
         CellMass = Rho0*dsig(k)*D(i,j)*dx*dy
         Energy = Energy+CellMass*Grav*(D(i,j)*sigc(k)-Hc(i,j))+  &
                    0.5*CellMass*(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
       enddo
     enddo
     enddo

     MaxEta = MAXVAL(Eta(Ibeg:Iend,Jbeg:Jend))
     MinEta = MINVAL(Eta(Ibeg:Iend,Jbeg:Jend))
     MaxU = MAXVAL(abs(U(Ibeg:Iend,Jbeg:Jend,Kbeg:Kend)))
     MaxV = MAXVAL(abs(V(Ibeg:Iend,Jbeg:Jend,Kbeg:Kend)))
     MaxW = MAXVAL(abs(W(Ibeg:Iend,Jbeg:Jend,Kbeg:Kend))) 


     call MPI_ALLREDUCE(MassVolume,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)        
     MassVolume = myvar
     call MPI_ALLREDUCE(Energy,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)            
     Energy = myvar
     call MPI_ALLREDUCE(MaxEta,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)            
     MaxEta = myvar
     call MPI_ALLREDUCE(MinEta,myvar,1,MPI_SP,MPI_MIN,MPI_COMM_WORLD,ier)            
     MinEta = myvar
     call MPI_ALLREDUCE(MaxU,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     MaxU = myvar
     call MPI_ALLREDUCE(MaxV,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     MaxV = myvar
     call MPI_ALLREDUCE(MaxW,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     MaxW = myvar

     if(myid.eq.0) then
     ! print screen
     WRITE(*,*),'----------------- STATISTICS ----------------'
     WRITE(*,*),' TIME        DT         DT_CONSTRAINT'
     WRITE(*,102) TIME,dt,TRIM(dt_constraint)
     WRITE(*,103) ' MassVolume  Energy      MaxEta      MinEta      MaxU       MaxV       MaxW       MaxS       MinS'
     WRITE(*,101) MassVolume,Energy,MaxEta,MinEta,MaxU,MaxV,MaxW,MaxS,MinS

     ! print log file 
     WRITE(3,*),'----------------- STATISTICS ----------------'
     WRITE(3,*),' TIME        DT         DT_CONSTRAINT'
     WRITE(3,102) TIME,dt,TRIM(dt_constraint)
     WRITE(3,103) ' MassVolume  Energy      MaxEta      MinEta      MaxU       MaxV       MaxW       MaxS       MinS'
     WRITE(3,101), MassVolume,Energy,MaxEta,MinEta,MaxU,MaxV,MaxW,MaxS,MinS
     endif

101  FORMAT(10E12.4)
102  FORMAT(2E12.4,A8)
103  FORMAT(A97)
 
     end subroutine statistics

 
     subroutine probes
!--------------------------------------------------
!    This subroutine is used to output probes
!    Called by
!       main
!    Last update: 16/11/2011, Gangfeng Ma
!--------------------------------------------------
     use global
     implicit none
     integer :: n,iu,i,j,k
     character(len=80) :: STAT_FILE,FDIR,FILE_NUM

     FDIR = TRIM(RESULT_FOLDER)
     
     do n = 1,NSTAT
       iu = 100+n
       write(FILE_NUM(1:4),'(I4.4)') n
       STAT_FILE = TRIM(FDIR)//'probe_'//TRIM(FILE_NUM)
       open(iu,file=TRIM(STAT_FILE),access='APPEND')
       
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         if(xstat(n)>=x(i).and.xstat(n)<=x(i+1).and.  &
            ystat(n)>=y(j).and.ystat(n)<=y(j+1)) then
           write(iu,'(100E12.4)') time,eta(i,j),(u(i,j,k),v(i,j,k),w(i,j,k),k=Kbeg,Kend)
         endif
       enddo
       enddo
       close(iu)
     enddo

     end subroutine probes


     subroutine preview
!--------------------------------------------------- 
!    This subroutine is used to preview
!    Called by                         
!       main 
!    Last update: 23/12/2010, Gangfeng Ma 
!--------------------------------------------------
     use global
     implicit none
     integer :: i,j,k,I1,I2,I3,I4,I5
     character(len=80) :: FDIR=''
     character(len=80) :: FILE_NAME=''
     character(len=80) :: file=''

     ! file number
     Icount = Icount+1
   
     ! results directory
     FDIR = TRIM(RESULT_FOLDER)

     if(myid.eq.0) write(*,102) 'Printing file No.',Icount,' TIME/TOTAL: ',TIME,'/',TOTAL_TIME
     if(myid.eq.0) write(3,102) 'Printing file No.',Icount,' TIME/TOTAL: ',TIME,'/',TOTAL_TIME     
102  FORMAT(A20,I5,A14,F8.3,A2,F8.3)
100  FORMAT(5000E16.6)

     I1 = mod(Icount/10000,10)
     I2 = mod(Icount/1000,10)
     I3 = mod(Icount/100,10)
     I4 = mod(Icount/10,10)
     I5 = mod(Icount,10)

     write(FILE_NAME(1:1),'(I1)') I1
     write(FILE_NAME(2:2),'(I1)') I2
     write(FILE_NAME(3:3),'(I1)') I3
     write(FILE_NAME(4:4),'(I1)') I4
     write(FILE_NAME(5:5),'(I1)') I5

     if(myid.eq.0) then
     open(5,file=TRIM(FDIR)//'time',position="append")
     write(5,*) TIME
     close(5)
     endif

     if(Icount==1) then
       if(OUT_H) then
         file=TRIM(FDIR)//'depth'
         call putfile2D(file,Hc0)
       endif
     endif

     if(OUT_E) then
       file=TRIM(FDIR)//'eta_'//TRIM(FILE_NAME)
       call putfile2D(file,Eta)
     endif

     if(OUT_U) then
       file=TRIM(FDIR)//'u_'//TRIM(FILE_NAME)
       call putfile3D(file,U)
     endif

     if(OUT_V) then
       file=TRIM(FDIR)//'v_'//TRIM(FILE_NAME)
       call putfile3D(file,V)
     endif

     if(OUT_W) then
       file=TRIM(FDIR)//'w_'//TRIM(FILE_NAME)
       call putfile3D(file,W)
     endif

     if(OUT_P) then
       file=TRIM(FDIR)//'p_'//TRIM(FILE_NAME)
       call putfile3D(file,P)
     endif

     if(OUT_K) then
       file=TRIM(FDIR)//'k_'//TRIM(FILE_NAME)
       call putfile3D(file,Tke)
     endif

     if(OUT_D) then
       file=TRIM(FDIR)//'d_'//TRIM(FILE_NAME)
       call putfile3D(file,Eps)
     endif

     if(OUT_S) then
       file=TRIM(FDIR)//'s_'//TRIM(FILE_NAME)
       call putfile3D(file,Prod_s)
     endif

     if(OUT_C) then
       file=TRIM(FDIR)//'c_'//TRIM(FILE_NAME)
       call putfile3D(file,CmuVt)
     endif

     if(OUT_A) then
       file=TRIM(FDIR)//'upwp_'//TRIM(FILE_NAME)
       call putfile3D(file,UpWp)
     endif







     end subroutine preview


    subroutine putfile2D(file,phi)
    use global
    implicit none
    real(SP),dimension(Mloc,Nloc),intent(in) :: phi
    character(len=80) :: file
    integer,dimension(NumP) :: npxs,npys
    integer,dimension(1) :: req
    real(SP),dimension(Mloc,Nloc) :: xx
    real(SP),dimension(Mglob,Nglob) :: phiglob
    integer,dimension(MPI_STATUS_SIZE,1) :: status
    integer :: i,j,iglob,jglob,len,n

    call MPI_GATHER(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,  &
           0,MPI_COMM_WORLD,ier)
    call MPI_GATHER(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,  &
           0,MPI_COMM_WORLD,ier)

    ! put the data in master processor into the global var
    if(myid==0) then
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        iglob = i-Nghost
        jglob = j-Nghost
        phiglob(iglob,jglob) = Phi(i,j)
      enddo
      enddo
    endif

    ! collect data from other processors into the master processor
    len = Mloc*Nloc

    do n = 1,NumP-1
      if(myid==0) then
        call MPI_IRECV(xx,len,MPI_SP,n,0,MPI_COMM_WORLD,req(1),ier)
        call MPI_WAITALL(1,req,status,ier)
        do j = Jbeg,Jend
        do i = Ibeg,Iend
          iglob = npxs(n+1)*(Iend-Ibeg+1)+i-Nghost
          jglob = npys(n+1)*(Jend-Jbeg+1)+j-Nghost
          phiglob(iglob,jglob) = xx(i,j)
        enddo
        enddo
      endif

      if(myid==n) then
        call MPI_SEND(phi,len,MPI_SP,0,0,MPI_COMM_WORLD,ier)
      endif
    enddo       

    if(myid==0) then
      open(5,file=TRIM(file))
      do j = 1,Nglob
        write(5,100) (phiglob(i,j),i=1,Mglob)
      enddo
      close(5)
    endif
100 FORMAT(5000f15.6)

    end subroutine putfile2D


    subroutine putfile3D(file,phi)
    use global
    implicit none
    real(SP),dimension(Mloc,Nloc,Kloc),intent(in) :: phi
    character(len=80) :: file
    integer,dimension(NumP) :: npxs,npys
    integer,dimension(1) :: req
    real(SP),dimension(:,:),allocatable :: xx,philoc
    real(SP),dimension(Mglob,Nglob,Kglob) :: phiglob
    integer,dimension(MPI_STATUS_SIZE,1) :: status
    integer :: i,j,k,jk,iglob,jglob,kk,n,len,nreq,NKloc

    call MPI_GATHER(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,  &
          0,MPI_COMM_WORLD,ier)
    call MPI_GATHER(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,  &
          0,MPI_COMM_WORLD,ier)

    NKloc = Nloc*Kloc

    ! put the data in master processor into the global var
    if(myid==0) then
      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        iglob = i-Nghost
        jglob = j-Nghost
        kk = k-Nghost
        phiglob(iglob,jglob,kk) = Phi(i,j,k)
      enddo
      enddo
      enddo
    endif

    allocate(philoc(Mloc,NKloc))
    allocate(xx(Mloc,NKloc))

    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      jk = (k-1)*Nloc+j
      philoc(i,jk) = phi(i,j,k)
    enddo
    enddo
    enddo

    ! collect data from other processors into the master processor
    len = Mloc*NKloc

    do n = 1,NumP-1
      if(myid==0) then
        call MPI_IRECV(xx,len,MPI_SP,n,0,MPI_COMM_WORLD,req(1),ier)
        call MPI_WAITALL(1,req,status,ier)
        do k = Kbeg,Kend
        do j = Jbeg,Jend
        do i = Ibeg,Iend
          iglob = npxs(n+1)*(Iend-Ibeg+1)+i-Nghost
          jglob = npys(n+1)*(Jend-Jbeg+1)+j-Nghost
          kk = k-Nghost
          jk = (k-1)*Nloc+j
          phiglob(iglob,jglob,kk) = xx(i,jk)
        enddo
        enddo
        enddo
      endif

      if(myid==n) then
        call MPI_SEND(philoc,len,MPI_SP,0,0,MPI_COMM_WORLD,ier)
      endif
    enddo

    if(myid.eq.0) then
      open(5,file=TRIM(file))
      do k = 1,Kglob
      do j = 1,Nglob
        write(5,100) (phiglob(i,j,k),i=1,Mglob)
      enddo
      enddo
      close(5)
    endif
100 FORMAT(5000f15.6)

    deallocate(philoc)
    deallocate(xx)

    end subroutine putfile3D

 
     subroutine estimate_dt
!----------------------------------------------------
!    This subroutine is used to estimate dt
!    Called by
!       main
!    Last update: 22/12/2010, Gangfeng Ma
!---------------------------------------------------
     use global
     implicit none
     integer :: i,j,k
     real(SP) :: tmp1,tmp2,dxonu,dyonv,dzonw,dt_growth,dt_courant,dt_viscous
     real(SP) :: myvar

     ! save previous time step
     dt_old = dt
     dt_growth = 1.05*dt_old     

     tmp2 = Large
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       tmp1 = abs(U(i,j,k))+sqrt(Grav*D(i,j))
       tmp1 = max(tmp1,Small)
       dxonu = dx/tmp1
       if(dxonu<tmp2) tmp2=dxonu

       tmp1 = abs(V(i,j,k))+sqrt(Grav*D(i,j))
       tmp1 = max(tmp1,Small)
       dyonv = dy/tmp1
       if(dyonv<tmp2) tmp2=dyonv

       tmp1 = max(abs(W(i,j,k)),Small)
       dzonw = dsig(k)*D(i,j)/tmp1
       if(dzonw<tmp2) tmp2=dzonw
     enddo
     enddo
     enddo


     call MPI_ALLREDUCE(tmp2,myvar,1,MPI_SP,MPI_MIN,MPI_COMM_WORLD,ier)
     tmp2 = myvar
     dt_courant = CFL*tmp2

     ! time step limit due to explicit viscous stress terms
     dt_viscous = Large
     if(VISCOUS_FLOW) then
       tmp2 = Large
       do k = Kbeg,Kend
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         tmp1 = dx**2/(abs(CmuHt(i,j,k))+1.e-16)
         if(tmp1<tmp2) tmp2 = tmp1

         tmp1 = dy**2/(abs(CmuHt(i,j,k))+1.e-16)
         if(tmp1<tmp2) tmp2 = tmp1
       enddo
       enddo
       enddo
       call MPI_ALLREDUCE(tmp2,myvar,1,MPI_SP,MPI_MIN,MPI_COMM_WORLD,ier)
       tmp2 = myvar
       dt_viscous = VISCOUS_NUMBER*tmp2
     endif 

     ! get dt    
     dt = min(dt_growth,dt_courant,dt_viscous,dt_max)
     if(dt<dt_min) then
       if(myid.eq.0) then
         write(3,*) 'time step too small !!',dt,dt_courant,dt_viscous
         stop
       endif
     endif
     TIME = TIME+dt
     RUN_STEP = RUN_STEP+1 
     if(myid.eq.0) write(3,*) RUN_STEP,dt,TIME

     if(dt==dt_growth) then
       dt_constraint = 'GROWTH'
     elseif(dt==dt_courant) then
       dt_constraint = 'COURANT'
     elseif(dt==dt_viscous) then
       dt_constraint = 'VISCOUS'
     elseif(dt==dt_max) then
       dt_constraint = 'MAXIMUM'
     endif  

     end subroutine estimate_dt


     subroutine vel_bc
!----------------------------------------------------
!    Boundary conditions for velocity
!    Called by 
!       main and get_UVW
!    Last update: 01/02/2011, Gangfeng Ma
!---------------------------------------------------
     use global
     implicit none
     integer :: i,j,k,imask
     real(SP) :: Wtop,Wbot,Cdrag,Phi,Dz1,Cg

     ! left and right boundary
     if(n_west.eq.MPI_PROC_NULL) then
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       if(Bc_X0==1) then  ! free-slip wall
         do i = 1,Nghost
           U(Ibeg-i,j,k) = -U(Ibeg+i-1,j,k)
           V(Ibeg-i,j,k) = V(Ibeg+i-1,j,k)
           W(Ibeg-i,j,k) = W(Ibeg+i-1,j,k)
           DU(Ibeg-i,j,k) = -DU(Ibeg+i-1,j,k)
           DV(Ibeg-i,j,k) = DV(Ibeg+i-1,j,k)
           DW(Ibeg-i,j,k) = DW(Ibeg+i-1,j,k)
         enddo
       elseif(Bc_X0==2) then ! no-slip wall
         do i =1,Nghost
           U(Ibeg-i,j,k) = -U(Ibeg+i-1,j,k)
           V(Ibeg-i,j,k) = -V(Ibeg+i-1,j,k)
           W(Ibeg-i,j,k) = -W(Ibeg+i-1,j,k)
           DU(Ibeg-i,j,k) = -DU(Ibeg+i-1,j,k)
           DV(Ibeg-i,j,k) = -DV(Ibeg+i-1,j,k)
           DW(Ibeg-i,j,k) = -DW(Ibeg+i-1,j,k)
         enddo
       elseif(Bc_X0==3) then ! inflow and outflow
         if(WaveMaker(1:7)=='LEF_TID') then ! for long-wave
           do i =1,Nghost
             U(Ibeg-i,j,k) = U(Ibeg+i-1,j,k)
             V(Ibeg-i,j,k) = V(Ibeg+i-1,j,k)
             W(Ibeg-i,j,k) = W(Ibeg+i-1,j,k)
             DU(Ibeg-i,j,k) = DU(Ibeg+i-1,j,k)
             DV(Ibeg-i,j,k) = DV(Ibeg+i-1,j,k)
             DW(Ibeg-i,j,k) = DW(Ibeg+i-1,j,k)
           enddo
         else
           do i = 1,Nghost
             U(Ibeg-i,j,k) = 2.0*Uin_X0(j,k)-U(Ibeg+i-1,j,k)
             V(Ibeg-i,j,k) = 2.0*Vin_X0(j,k)-V(Ibeg+i-1,j,k)
             W(Ibeg-i,j,k) = 2.0*Win_X0(j,k)-W(Ibeg+i-1,j,k)
             DU(Ibeg-i,j,k) = 2.0*Din_X0(j)*Uin_X0(j,k)-DU(Ibeg+i-1,j,k)
             DV(Ibeg-i,j,k) = 2.0*Din_X0(j)*Vin_X0(j,k)-DV(Ibeg+i-1,j,k)
             DW(Ibeg-i,j,k) = 2.0*Din_X0(j)*Win_X0(j,k)-DW(Ibeg+i-1,j,k)
           enddo
         endif
       elseif(Bc_X0==4) then
         do i =1,Nghost
           U(Ibeg-i,j,k) = U(Ibeg+i-1,j,k)
           V(Ibeg-i,j,k) = V(Ibeg+i-1,j,k)
           W(Ibeg-i,j,k) = W(Ibeg+i-1,j,k)
           DU(Ibeg-i,j,k) = DU(Ibeg+i-1,j,k)
           DV(Ibeg-i,j,k) = DV(Ibeg+i-1,j,k)
           DW(Ibeg-i,j,k) = DW(Ibeg+i-1,j,k)
         enddo
       elseif(Bc_X0==8) then  ! specify u,v,w at ghost cells
         do i = 1,Nghost
           U(Ibeg-i,j,k) = Uin_X0(j,k)
           V(Ibeg-i,j,k) = Vin_X0(j,k)
           W(Ibeg-i,j,k) = Win_X0(j,k)
           DU(Ibeg-i,j,k) = Din_X0(j)*Uin_X0(j,k)
           DV(Ibeg-i,j,k) = Din_X0(j)*Vin_X0(j,k)
           DW(Ibeg-i,j,k) = Din_X0(j)*Win_X0(j,k)
         enddo
       endif
     enddo
     enddo
     endif

     if(n_east.eq.MPI_PROC_NULL) then
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       if(Bc_Xn==1) then  ! free-slip wall 
         do i = 1,Nghost
           U(Iend+i,j,k) = -U(Iend-i+1,j,k)
           V(Iend+i,j,k) = V(Iend-i+1,j,k)
           W(Iend+i,j,k) = W(Iend-i+1,j,k)
           DU(Iend+i,j,k) = -DU(Iend-i+1,j,k)
           DV(Iend+i,j,k) = DV(Iend-i+1,j,k)
           DW(Iend+i,j,k) = DW(Iend-i+1,j,k)
         enddo
       elseif(Bc_Xn==2) then ! no-slip wall
         do i = 1,Nghost
           U(Iend+i,j,k) = -U(Iend-i+1,j,k)
           V(Iend+i,j,k) = -V(Iend-i+1,j,k)
           W(Iend+i,j,k) = -W(Iend-i+1,j,k)
           DU(Iend+i,j,k) = -DU(Iend-i+1,j,k)
           DV(Iend+i,j,k) = -DV(Iend-i+1,j,k)
           DW(Iend+i,j,k) = -DW(Iend-i+1,j,k)
         enddo
       elseif(Bc_Xn==3) then
         do i = 1,Nghost
           U(Iend+i,j,k) = 2.0*Uin_Xn(j,k)-U(Iend-i+1,j,k)
           V(Iend+i,j,k) = 2.0*Vin_Xn(j,k)-V(Iend-i+1,j,k)
           W(Iend+i,j,k) = 2.0*Win_Xn(j,k)-W(Iend-i+1,j,k)
           DU(Iend+i,j,k) = 2.0*Din_Xn(j)*Uin_Xn(j,k)-DU(Iend-i+1,j,k)
           DV(Iend+i,j,k) = 2.0*Din_Xn(j)*Vin_Xn(j,k)-DV(Iend-i+1,j,k)
           DW(Iend+i,j,k) = 2.0*Din_Xn(j)*Win_Xn(j,k)-DW(Iend-i+1,j,k)
         enddo
       elseif(Bc_Xn==4) then 
         do i = 1,Nghost
           U(Iend+i,j,k) = U(Iend-i+1,j,k)
           V(Iend+i,j,k) = V(Iend-i+1,j,k)
           W(Iend+i,j,k) = W(Iend-i+1,j,k)
           DU(Iend+i,j,k) = DU(Iend-i+1,j,k)
           DV(Iend+i,j,k) = DV(Iend-i+1,j,k)
           DW(Iend+i,j,k) = DW(Iend-i+1,j,k)
         enddo
       elseif(Bc_Xn==6) then
         do i = 1,Nghost
!           Cg = -(U0(Iend+i-1,j,k)-U00(Iend+i-1,j,k)+1.e-16)/  &
!                 (U00(Iend+i-1,j,k)-U00(Iend+i-2,j,k)+1.e-16)
           Cg = sqrt(Grav*D(Iend,j))*dt/dx
           Cg = max(min(Cg,1.0),0.0)
           U(Iend+i,j,k) = Cg*U0(Iend+i-1,j,k)+(1.0-Cg)*U0(Iend+i,j,k)

!           Cg =-(V0(Iend+i-1,j,k)-V00(Iend+i-1,j,k)+1.e-16)/  &
!                 (V00(Iend+i-1,j,k)-V00(Iend+i-2,j,k)+1.e-16)
!           Cg= max(min(Cg,1.0),0.0)
           V(Iend+i,j,k) = Cg*V0(Iend+i-1,j,k)+(1.0-Cg)*V0(Iend+i,j,k)
           
!           Cg =-(W0(Iend+i-1,j,k)-W00(Iend+i-1,j,k)+1.e-16)/  &
!                 (W00(Iend+i-1,j,k)-W00(Iend+i-2,j,k)+1.e-16)
!           Cg= max(min(Cg,1.0),0.0)
           W(Iend+i,j,k) = Cg*W0(Iend+i-1,j,k)+(1.0-Cg)*W0(Iend+i,j,k)

           DU(Iend+i,j,k) = D(Iend+i,j)*U(Iend+i,j,k)
           DV(Iend+i,j,k) = D(Iend+i,j)*V(Iend+i,j,k)
           DW(Iend+i,j,k) = D(Iend+i,j)*W(Iend+i,j,k)
         enddo 
       elseif(Bc_Xn==8) then
         do i = 1,Nghost
           U(Iend+i,j,k) = Uin_Xn(j,k)
           V(Iend+i,j,k) = Vin_Xn(j,k)
           W(Iend+i,j,k) = Win_Xn(j,k)
           DU(Iend+i,j,k) = Din_Xn(j)*Uin_Xn(j,k)
           DV(Iend+i,j,k) = Din_Xn(j)*Vin_Xn(j,k)
           DW(Iend+i,j,k) = Din_Xn(j)*Win_Xn(j,k)
         enddo
       endif
     enddo
     enddo
     endif

     if(n_suth.eq.MPI_PROC_NULL) then
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       if(Bc_Y0==1) then  ! free-slip wall 
         do j = 1,Nghost
           U(i,Jbeg-j,k) = U(i,Jbeg+j-1,k)
           V(i,Jbeg-j,k) = -V(i,Jbeg+j-1,k)
           W(i,Jbeg-j,k) = W(i,Jbeg+j-1,k)
           DU(i,Jbeg-j,k) = DU(i,Jbeg+j-1,k)
           DV(i,Jbeg-j,k) = -DV(i,Jbeg+j-1,k)
           DW(i,Jbeg-j,k) = DW(i,Jbeg+j-1,k)
         enddo
       elseif(Bc_Y0==2) then ! no-slip wall 
         do j = 1,Nghost
           U(i,Jbeg-j,k) = -U(i,Jbeg+j-1,k)
           V(i,Jbeg-j,k) = -V(i,Jbeg+j-1,k)
           W(i,Jbeg-j,k) = -W(i,Jbeg+j-1,k)
           DU(i,Jbeg-j,k) = -DU(i,Jbeg+j-1,k)
           DV(i,Jbeg-j,k) = -DV(i,Jbeg+j-1,k)
           DW(i,Jbeg-j,k) = -DW(i,Jbeg+j-1,k)
         enddo
       elseif(Bc_Y0==4) then
         do j = 1,Nghost
           U(i,Jbeg-j,k) = U(i,Jbeg+j-1,k)
           V(i,Jbeg-j,k) = V(i,Jbeg+j-1,k)
           W(i,Jbeg-j,k) = W(i,Jbeg+j-1,k)
           DU(i,Jbeg-j,k) = DU(i,Jbeg+j-1,k)
           DV(i,Jbeg-j,k) = DV(i,Jbeg+j-1,k)
           DW(i,Jbeg-j,k) = DW(i,Jbeg+j-1,k)
         enddo
       endif
     enddo
     enddo
     endif

     if(n_nrth.eq.MPI_PROC_NULL) then
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       if(Bc_Yn==1) then  ! free-slip wall 
         do j = 1,Nghost
           U(i,Jend+j,k) = U(i,Jend-j+1,k)
           V(i,Jend+j,k) = -V(i,Jend-j+1,k)
           W(i,Jend+j,k) = W(i,Jend-j+1,k)
           DU(i,Jend+j,k) = DU(i,Jend-j+1,k)
           DV(i,Jend+j,k) = -DV(i,Jend-j+1,k)
           DW(i,Jend+j,k) = DW(i,Jend-j+1,k)
         enddo
       elseif(Bc_Yn==2) then ! no-slip wall 
         do j = 1,Nghost
           U(i,Jend+j,k) = -U(i,Jend-j+1,k)
           V(i,Jend+j,k) = -V(i,Jend-j+1,k)
           W(i,Jend+j,k) = -W(i,Jend-j+1,k)
           DU(i,Jend+j,k) = -DU(i,Jend-j+1,k)
           DV(i,Jend+j,k) = -DV(i,Jend-j+1,k)
           DW(i,Jend+j,k) = -DW(i,Jend-j+1,k)
         enddo
       elseif(Bc_Yn==4) then
         do j = 1,Nghost
           U(i,Jend+j,k) = U(i,Jend-j+1,k)
           V(i,Jend+j,k) = V(i,Jend-j+1,k)
           W(i,Jend+j,k) = W(i,Jend-j+1,k)
           DU(i,Jend+j,k) = DU(i,Jend-j+1,k)
           DV(i,Jend+j,k) = DV(i,Jend-j+1,k)
           DW(i,Jend+j,k) = DW(i,Jend-j+1,k)
         enddo
       endif
     enddo
     enddo
     endif

     ! top and bottom
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Dz1 = 0.5*D(i,j)*dsig(Kbeg)
       if(ibot==1) then
         Cdrag = Cd0
       else
         Cdrag = 1./(1./Kappa*log(30.0*Dz1/Zob))**2
       endif
       Phi = 2.0*Dz1*Cdrag*sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)/(Cmu(i,j,Kbeg)+CmuVt(i,j,Kbeg))
       Phi = dmin1(Phi,2.0)

       if(Bc_Z0==1) then  ! free-slip
         Wbot = -DeltH(i,j)-U(i,j,Kbeg)*DelxH(i,j)-V(i,j,Kbeg)*DelyH(i,j)
         do k = 1,Nghost
           U(i,j,Kbeg-k) = U(i,j,Kbeg+k-1)
           V(i,j,Kbeg-k) = V(i,j,Kbeg+k-1)
           W(i,j,Kbeg-k) = 2.0*Wbot-W(i,j,Kbeg+k-1)
           DU(i,j,Kbeg-k) = D(i,j)*U(i,j,Kbeg-k)
           DV(i,j,Kbeg-k) = D(i,j)*V(i,j,Kbeg-k)
           DW(i,j,Kbeg-k) = D(i,j)*W(i,j,Kbeg-k)
         enddo
       elseif(Bc_Z0==2) then  ! no-slip
         Wbot = -DeltH(i,j)
         do k = 1,Nghost
           U(i,j,Kbeg-k) = -U(i,j,Kbeg+k-1)
           V(i,j,Kbeg-k) = -V(i,j,Kbeg+k-1)
           W(i,j,Kbeg-k) = 2.0*Wbot-W(i,j,Kbeg+k-1)
           DU(i,j,Kbeg-k) = D(i,j)*U(i,j,Kbeg-k)
           DV(i,j,Kbeg-k) = D(i,j)*V(i,j,Kbeg-k)
           DW(i,j,Kbeg-k) = D(i,j)*W(i,j,Kbeg-k)
         enddo
       elseif(Bc_Z0==5) then
         do k = 1,Nghost
           U(i,j,Kbeg-k) = (1.0-Phi)*U(i,j,Kbeg+k-1)
           V(i,j,Kbeg-k) = (1.0-Phi)*V(i,j,Kbeg+k-1)
           Wbot = -DeltH(i,j)-0.5*(U(i,j,Kbeg)+U(i,j,Kbeg-1))*DelxH(i,j)-  &
                    0.5*(V(i,j,Kbeg)+V(i,j,Kbeg-1))*DelyH(i,j)
           W(i,j,Kbeg-k) = 2.0*Wbot-W(i,j,Kbeg+k-1)
           DU(i,j,Kbeg-k) = D(i,j)*U(i,j,Kbeg-k)
           DV(i,j,Kbeg-k) = D(i,j)*V(i,j,Kbeg-k)
           DW(i,j,Kbeg-k) = D(i,j)*W(i,j,Kbeg-k)
         enddo
       endif

       ! at the surface (no stress)
       Wtop = (Eta(i,j)-Eta0(i,j))/dt+U(i,j,Kend)*DelxEta(i,j)+V(i,j,Kend)*DelyEta(i,j)
       do k = 1,Nghost
         U(i,j,Kend+k) = U(i,j,Kend-k+1)
         V(i,j,Kend+k) = V(i,j,Kend-k+1)
         W(i,j,Kend+k) = 2.0*Wtop-W(i,j,Kend-k+1)
         DU(i,j,Kend+k) = D(i,j)*U(i,j,Kend+k)
         DV(i,j,Kend+k) = D(i,j)*V(i,j,Kend+k)
         DW(i,j,Kend+k) = D(i,j)*W(i,j,Kend+k)
       enddo
     enddo
     enddo

     ! fyshi added boundary conditions at masks 02/15/2013
     DO K=Kbeg,Kend
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
       IF(Mask(i,j)==0) THEN
         ! south boundary 
         IF(Mask(i,j+1)==1)then
           if(Bc_X0==1) then  ! free-slip wall 
             do imask = 1,Nghost
               U(i,j-imask+1,k) = U(i,j+imask,k)
               V(i,j-imask+1,k) = -V(i,j+imask,k)
               W(i,j-imask+1,k) = W(i,j+imask,k)
               DU(i,j-imask+1,k) = DU(i,j+imask,k)
               DV(i,j-imask+1,k) = -DV(i,j+imask,k)
               DW(i,j-imask+1,k) = DW(i,j+imask,k)
             enddo
           elseif(Bc_X0==2) then ! no-slip wall 
             do imask =1,Nghost
               U(i,j-imask+1,k) = -U(i,j+imask,k)
               V(i,j-imask+1,k) = -V(i,j+imask,k)
               W(i,j-imask+1,k) = -W(i,j+imask,k)
               DU(i,j-imask+1,k) = -DU(i,j+imask,k)
               DV(i,j-imask+1,k) = -DV(i,j+imask,k)
               DW(i,j-imask+1,k) = -DW(i,j+imask,k)
             enddo
           endif
         ! north  
         ELSEIF(Mask(i,j-1)==1)then
           if(Bc_X0==1) then  ! free-slip wall 
             do imask = 1,Nghost
               U(i,j+imask-1,k) = U(i,j-imask,k)
               V(i,j+imask-1,k) = -V(i,j-imask,k)
               W(i,j+imask-1,k) = W(i,j-imask,k)
               DU(i,j+imask-1,k) = DU(i,j-imask,k)
               DV(i,j+imask-1,k) = -DV(i,j-imask,k)
               DW(i,j+imask-1,k) = DW(i,j-imask,k)
             enddo
           elseif(Bc_X0==2) then ! no-slip wall 
             do imask =1,Nghost
               U(i,j+imask-1,k) = -U(i,j-imask,k)
               V(i,j+imask-1,k) = -V(i,j-imask,k)
               W(i,j+imask-1,k) = -W(i,j-imask,k)
               DU(i,j+imask-1,k) = -DU(i,j-imask,k)
               DV(i,j+imask-1,k) = -DV(i,j-imask,k)
               DW(i,j+imask-1,k) = -DW(i,j-imask,k)
             enddo
           endif
         ! west
         ELSEIF(Mask(i+1,j)==1)THEN
           if(Bc_X0==1) then  ! free-slip wall 
             do imask = 1,Nghost
               U(I-imask+1,j,k) = -U(I+imask,j,k)
               V(I-imask+1,j,k) = V(I+imask,j,k)
               W(I-imask+1,j,k) = W(I+imask,j,k)
               DU(I-imask+1,j,k) = -DU(I+imask,j,k)
               DV(I-imask+1,j,k) = DV(I+imask,j,k)
               DW(I-imask+1,j,k) = DW(I+imask,j,k)
             enddo
           elseif(Bc_X0==2) then ! no-slip wall
             do imask =1,Nghost
               U(I-imask+1,j,k) = -U(I+imask,j,k)
               V(I-imask+1,j,k) = -V(I+imask,j,k)
               W(I-imask+1,j,k) = -W(I+imask,j,k)
               DU(I-imask+1,j,k) = -DU(I+imask,j,k)
               DV(I-imask+1,j,k) = -DV(I+imask,j,k)
               DW(I-imask+1,j,k) = -DW(I+imask,j,k)
             enddo
           endif
         ! east 
         ELSEIF(Mask(i-1,j)==1)THEN
           if(Bc_X0==1) then  ! free-slip wall  
             do imask = 1,Nghost
               U(i+imask-1,j,k) = -U(i-imask,j,k)
               V(i+imask-1,j,k) = V(i-imask,j,k)
               W(i+imask-1,j,k) = W(i-imask,j,k)
               DU(i+imask-1,j,k) = -DU(i-imask,j,k)
               DV(i+imask-1,j,k) = DV(i-imask,j,k)
               DW(i+imask-1,j,k) = DW(i-imask,j,k)
             enddo
           elseif(Bc_X0==2) then ! no-slip wall 
             do imask =1,Nghost
               U(i+imask-1,j,k) = -U(i-imask,j,k)
               V(i+imask-1,j,k) = -V(i-imask,j,k)
               W(i+imask-1,j,k) = -W(i-imask,j,k)
               DU(i+imask-1,j,k) = -DU(i-imask,j,k)
               DV(i+imask-1,j,k) = -DV(i-imask,j,k)
               DW(i+imask-1,j,k) = -DW(i-imask,j,k)
             enddo
           endif
         ENDIF ! end mask+1=1 
       ENDIF ! end mask=0 
     ENDDO
     ENDDO
     ENDDO

     end subroutine vel_bc

     subroutine wl_bc
!-----------------------------------------------------------
!    Boundary condition for surface elevation or water depth
!    Called by
!       eval_duvw
!    Last update: 14/06/2012, Gangfeng Ma
!-----------------------------------------------------------
     use global
     implicit none
     real(SP) :: Cg
     integer :: i,j

     ! left and right boundary
     if(n_west.eq.MPI_PROC_NULL) then
     do j = Jbeg,Jend
       if(Bc_X0==1.or.Bc_X0==2) then  ! free/no-slip wall
         do i = 1,Nghost
           D(Ibeg-i,j) = D(Ibeg+i-1,j)
         enddo
       elseif(Bc_X0==3) then ! inflow
         if(WaveMaker(1:7)=='LEF_TID') then
           do i = 0,Nghost
             D(Ibeg-i,j) = Din_X0(j)
           enddo
         else
           do i = 1,Nghost
             D(Ibeg-i,j) = 2.0*Din_X0(j)-D(Ibeg+i-1,j)
           enddo
         endif
       elseif(Bc_X0==4) then ! outflow
         do i = 1,Nghost
           D(Ibeg-i,j) = Din_X0(j)
         enddo
       elseif(Bc_X0==8) then
         do i = 1,Nghost
           D(Ibeg-i,j) = Din_X0(j)
         enddo
       endif
     enddo
     endif

     if(n_east.eq.MPI_PROC_NULL) then
     do j = Jbeg,Jend
       if(Bc_Xn==1.or.Bc_Xn==2) then 
         do i = 1,Nghost
           D(Iend+i,j) = D(Iend-i+1,j)
         enddo
       elseif(Bc_Xn==3) then 
         do i = 1,Nghost
           D(Iend+i,j) = 2.0*Din_Xn(j)-D(Iend-i+1,j)
         enddo
       elseif(Bc_Xn==4) then
         do i = 1,Nghost
           D(Iend+i,j) = Din_Xn(j)
         enddo
       elseif(Bc_Xn==6) then 
         do i = 1,Nghost
!           Cg =-(Eta0(Iend+i-1,j)-Eta00(Iend+i-1,j)+1.e-16)/  &
!                 (Eta00(Iend+i-1,j)-Eta00(Iend+i-2,j)+1.e-16)
!           Cg= max(min(Cg,1.0),0.0)
           Cg = sqrt(Grav*D(Iend,j))*dt/dx
           Eta(Iend+i,j) = Cg*Eta0(Iend+i-1,j)+(1.0-Cg)*Eta0(Iend+i,j)
           D(Iend+i,j) = Hc(Iend+i,j)+Eta(Iend+i,j)
         enddo
       elseif(Bc_Xn==8) then
         do i = 1,Nghost
           D(Iend+i,j) = Din_Xn(j)
         enddo
       endif
     enddo
     endif

! y-direction and corners                                                                                                     
     if(n_suth.eq.MPI_PROC_NULL) then
       do i = 1,Mloc
       do j = 1,Nghost
         D(i,j) = D(i,Jbeg+Nghost-j)
       enddo
       enddo
     endif

     if(n_nrth.eq.MPI_PROC_NULL) then
       do i = 1,Mloc
       do j = 1,Nghost
         D(i,Jend+j) = D(i,Jend-j+1)
       enddo
       enddo
     endif

     call phi_2D_exch(D)
     
     return
     end subroutine wl_bc

     subroutine phi_2D_coll(phi)
!-----------------------------------------------------
!    This subroutine is used to collect data into ghost cells
!    Called by
!       eval_duvw
!    Last update: 22/12/2010, Gangfeng Ma
!-----------------------------------------------------
     use global
     implicit none
     real(SP), intent(inout) :: phi(Mloc,Nloc)
     integer :: i,j

     ! x-direction
     if(n_west.eq.MPI_PROC_NULL) then
       do j = Jbeg,Jend
       do i = 1,Nghost
         phi(i,j) = phi(Ibeg+Nghost-i,j)
       enddo
       enddo
     endif

     if(n_east.eq.MPI_PROC_NULL) then
       do j = Jbeg,Jend
       do i = 1,Nghost
         phi(Iend+i,j) = phi(Iend-i+1,j)
       enddo
       enddo
     endif
 
     ! y-direction and corners
     if(n_suth.eq.MPI_PROC_NULL) then
       do i = 1,Mloc
       do j = 1,Nghost
         phi(i,j) = phi(i,Jbeg+Nghost-j)
       enddo
       enddo
     endif

     if(n_nrth.eq.MPI_PROC_NULL) then
       do i = 1,Mloc
       do j = 1,Nghost
         phi(i,Jend+j) = phi(i,Jend-j+1)
       enddo
       enddo
     endif

     call phi_2D_exch(phi)

     end subroutine phi_2D_coll


     subroutine update_wind
!--------------------------------------------------------
!    Update wind speed at current time step
!    Called by
!       main
!    Last update: 09/07/2013, Gangfeng Ma
!--------------------------------------------------------
     use global
     implicit none
     integer :: i,j,iglob,jglob
     real(SP),dimension(Mglob,Nglob) :: WUG,WVG
     real(SP) :: Wds,Cds,Cdsmin,Cdsmax

     ! update wind speed
     if(Iws==1) then
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         WdU(i,j) = WindU
         WdV(i,j) = WindV
       enddo
       enddo
     elseif(Iws==2) then
       open(7,file='wind.txt')

       do j = 1,Nglob
       do i = 1,Mglob
         read(7,*) WUG(i,j),WVG(i,j)
       enddo
       enddo

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         iglob = npx*(Mloc-2*Nghost)+i-Nghost
         jglob = npy*(Nloc-2*Nghost)+j-Nghost
         WdU(i,j) = WUG(iglob,jglob)
         WdV(i,j) = WVG(iglob,jglob)
       enddo
       enddo

     endif

     ! wind drag coefficient
     Cdsmin = 1.e-3*(0.61+0.063*6.0)
     Cdsmax = 1.e-3*(0.61+0.063*50.0)
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Wds = sqrt(WdU(i,j)**2+WdV(i,j)**2)
       Cds = 1.e-3*(0.61+0.063*Wds)
       Cds = dmin1(dmax1(Cds,Cdsmin),Cdsmax)
       
       Wsx(i,j) = 0.001293*Cds*WdU(i,j)*Wds
       Wsy(i,j) = 0.001293*Cds*WdV(i,j)*Wds
     enddo
     enddo

     return
     end subroutine


     subroutine update_mask
!------------------------------------------------------  
!    This subroutine is used to update mask for wetting-drying
!    Called by                                                
!       main
!    Last update: 22/12/2010, Gangfeng Ma 
!-----------------------------------------------------
     use global, only: Ibeg,Iend,Jbeg,Jend,Eta,Hc,D,MinDep,  &
                       Mask,Mask_Struct,Mask9
     implicit none
     integer :: i,j

     ! Mask at ghost cells keeps no change
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask_Struct(i,j)==0) cycle
       
       ! flooding (dry->wet)
       if(Mask(i,j)==0) then
         if(Mask(i-1,j)==1.and.Eta(i-1,j)>Eta(i,j)) Mask(i,j)=1
         if(Mask(i+1,j)==1.and.Eta(i+1,j)>Eta(i,j)) Mask(i,j)=1
         if(Mask(i,j-1)==1.and.Eta(i,j-1)>Eta(i,j)) Mask(i,j)=1
         if(Mask(i,j+1)==1.and.Eta(i,j+1)>Eta(i,j)) Mask(i,j)=1
       else
         ! drying (wet->dry)
         if(abs(D(i,j)-MinDep)<=1.e-6) then
           Mask(i,j) = 0
           Eta(i,j) = MinDep-Hc(i,j)
           D(i,j) = Eta(i,j)+Hc(i,j)           
         endif
       endif
     enddo
     enddo

     Mask = Mask*Mask_Struct

     ! collect mask into ghost cells
     call phi_int_exch(Mask)    

     do j = Jbeg,Jend
     do i = Ibeg,Iend
      Mask9(i,j) = Mask(i,j)*Mask(i-1,j)*Mask(i+1,j)  &
                *Mask(i+1,j+1)*Mask(i,j+1)*Mask(i-1,j+1) &
                *Mask(i+1,j-1)*Mask(i,j-1)*Mask(i-1,j-1)
     enddo
     enddo

     end subroutine update_mask


     subroutine update_vars
!------------------------------------------------------ 
!    This subroutine is used to save variables at 
!    last time step
!    Called by   
!       main 
!    Last update: 22/12/2010, Gangfeng Ma 
!----------------------------------------------------- 
     use global        
     implicit none

     Eta00 = Eta0
     U00 = U0
     V00 = V0
     W00 = W0

     D0 = D
     Eta0 = Eta
     U0 = U
     V0 = V
     W0 = W
     DU0 = DU
     DV0 = DV
     DW0 = DW
     DTke0 = DTke
     DEps0 = DEps





     end subroutine update_vars

  
     subroutine update_wave_bc
!------------------------------------------------------
!    This subroutine is used to update boundary conditions
!    Called by
!       main
!    Last update: 22/12/2010, Gangfeng Ma 
!-----------------------------------------------------
     use global
     implicit none

     if(n_west.eq.MPI_PROC_NULL) then
     if(WaveMaker(1:7)=='LEF_SOL') then
       call solitary_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_LIN') then
       call linear_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_CON') then
       call cnoidal_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_STK') then
       call stokes_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_SPC') then
       call random_wave_left_boundary
     elseif((WaveMaker(1:7)=='LEF_JON').or.(WaveMaker(1:7)=='LEF_TMA')) then
       call jonswap_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_TID') then
       call tidal_wave_left_boundary
     endif
     endif


     if(n_east.eq.MPI_PROC_NULL) then
     if(WaveMaker(1:7)=='RIG_LIN') then
       call linear_wave_right_boundary
     endif
     endif


     if(WaveMaker(1:7)=='FLUX_LR') then
       call flux_left_right_boundary
     endif

     end subroutine update_wave_bc 


     subroutine flux_left_right_boundary
!-----------------------------------------------------------   
!    This subroutine is used to specify left/right boundary                                         
!    Called by 
!       update_wave_bc 
!    Last update: 14/06/2012, Gangfeng Ma
!-----------------------------------------------------------
     use global
     implicit none
     integer :: j,k,n
     real(SP) :: Zlev1,Zlev2,Uavg_Left,Uavg_Right,Ufric,sintep,UU,  &
                 FluxL,FluxR,myvar,Ramp

     Uavg_Left = 0.116
     Uavg_Right = 0.116

     FluxL = 0.0
     FluxR = 0.0

     if(TRamp>0.0) then
       Ramp = tanh(TIME/TRamp)
     else
       Ramp = 1.0
     endif

     if(n_west.eq.MPI_PROC_NULL) then
     do j = 1,Nloc
       Ein_X0(j) = 1.5*Eta(Ibeg,j)-0.5*Eta(Ibeg+1,j)
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)
     enddo

     ! log-profile of velocity
     Ufric = Uavg_Left*Kappa/(log(30.*Din_X0(Jbeg)/Zob)-1) 

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev1 = sigc(k)*Din_X0(j)
       Uin_X0(j,k) = Ufric/Kappa*log(30.*Zlev1/Zob)*Ramp
       FluxL = FluxL+dsig(k)*Din_X0(j)*dy*Uin_X0(j,k)
       Win_X0(j,k) = 0.0
       Vin_X0(j,k) = 0.0


     enddo
     enddo
     endif

     call MPI_ALLREDUCE(FluxL,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)
     FluxL = myvar

     if(n_east.eq.MPI_PROC_NULL) then

!     if(TIME>TRamp) then
!       call linear_wave_right_boundary
!     endif

     do j = 1,Nloc
!       Ein_Xn(j) = 1.5*Eta(Iend,j)-0.5*Eta(Iend-1,j)
       Ein_Xn(j) = 0.0
       Din_Xn(j) = Ein_Xn(j)+Hfx(Iend+1,j)
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev2 = sigc(k)*Din_Xn(j)
       Uin_Xn(j,k) = FluxL/(Din_Xn(Jbeg)*dy*float(Nglob))
       Win_Xn(j,k) = 0.0
       Vin_Xn(j,k) = 0.0


     enddo
     enddo
     endif
     
     return
     end subroutine flux_left_right_boundary


     subroutine stokes_wave_left_boundary
!-----------------------------------------------------
!    This subroutine is used to specify left boundary
!    Called by
!       update_wave_bc
!    Last update: 26/04/2011, Gangfeng Ma
!-----------------------------------------------------
     use global, only: SP,pi,Zero,Ibeg,Grav,TIME,Nloc,Kloc,Amp_Wave,Per_Wave,Dep_Wave, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend,sigc
     implicit none
     integer  :: j,k,Iter
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Zlev

     ! Find wave number for linear wave (Newton-Ralphson Method)
     Segma = 2.0*pi/Per_Wave
     Celerity = sqrt(Grav*Dep_Wave)
     Wave_Length = Celerity*Per_Wave
     Wave_Number = 2.0*pi/Wave_Length
     
     Iter = 0
55   Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
     if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 65
     Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &
        Grav*tanh(Wave_Number*Dep_Wave) 
     Wave_Number = Wave_Number-Fk/Fkdif
     Iter = Iter+1
     goto 55
65   continue
     
     Wave_Length = 2.0*pi/Wave_Number
     Celerity = Wave_Length/Per_Wave
     
     do j = 1,Nloc
       Ein_X0(j) = 0.5*Amp_Wave*cos(pi/2-Segma*TIME)+  &
                   Amp_Wave**2*Wave_Number/16.*cosh(Wave_Number*Dep_Wave)/  &
                   sinh(Wave_Number*Dep_Wave)**3*(2.0+cosh(2.*Wave_Number*Dep_Wave))* &
                   cos(2.*(pi/2-Segma*TIME))
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)     
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)
       Uin_X0(j,k) = 0.5*Amp_Wave*Segma*cosh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*cos(pi/2-Segma*TIME)+  &
           3./16.*Amp_Wave**2*Segma*Wave_Number*cosh(2.*Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)**4*cos(2.*(pi/2-Segma*TIME))
       Win_X0(j,k) = 0.5*Amp_Wave*Segma*sinh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*sin(pi/2-Segma*TIME)+  &
           3./16.*Amp_Wave**2*Segma*Wave_Number*sinh(2.*Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)**4*sin(2.*(pi/2-Segma*TIME))
       Vin_X0(j,k) = 0.0
     enddo
     enddo

     end subroutine stokes_wave_left_boundary


     subroutine tidal_wave_left_boundary
!-----------------------------------------------------
!    This subroutine is used to specify left boundary
!    Called by
!       update_wave_bc
!    Last update: 06/02/2011, Gangfeng Ma
!-----------------------------------------------------
     use global
     implicit none
     integer :: i,j,k
     real(SP), parameter, dimension(8) :: &
                 !  s2      m2    n2     k2     k1     p1     o1     q1
        period = (/43200.,44712.,45570.,43082.,86164.,86637.,92950.,96726./)

     do j = 1,Nloc
       Ein_X0(j) = 0.5*Amp_Wave*cos(2.0*pi/period(1)*TIME-pi/2.0)
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Uin_X0(j,k) = sqrt(Grav/Hfx(Ibeg,j))*Ein_X0(j)
       Vin_X0(j,k) = 0.0
       Win_X0(j,k) = 0.0
     enddo
     enddo


     end subroutine tidal_wave_left_boundary


     subroutine cnoidal_wave_left_boundary
!-----------------------------------------------------
!    This subroutine is used to specify left boundary 
!    Called by 
!       update_wave_bc 
!    Last update: 06/02/2011, Gangfeng Ma 
!----------------------------------------------------- 
     use global, only: SP,pi,Zero,Ibeg,Grav,TIME,Nloc,Kloc,Amp_Wave,Per_Wave,Dep_Wave, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend,sigc
     implicit none
     integer :: j,k,nn
     real(SP) :: Wave_Length,Celerity,Ytrough,Mod1,Wave_Number,Zup,Zlow,Zmid,Xstart, &
                 Zero1,Zlev,Xcn,Xsn,Xdn,cnoidal_cn,cnoidal_ck,Stokes_Drift,Fact

     call cnoidal(Amp_Wave,Dep_Wave,Per_Wave,Wave_Length,Celerity,Ytrough,Mod1)
     
     ! wave number
     Wave_Number = 2.0*pi/Wave_Length

     ! find zero start
     Zup = 1.0
     Zlow = 0.0
     Zmid= (Zup+Zlow)/2.0
     nn = 0
200  nn = nn+1
     Zero1 = Ytrough+Amp_Wave*cnoidal_cn(Zmid*0.5*cnoidal_ck(Mod1),Mod1)**2                            

     if(abs(Zero1)<=1.0e-6) goto 210
     if(nn>1000) then
       write(*,*)'too many iterations; stop'
       stop
     endif
     if(Zero1<0.0) then
       Zup = Zmid
       Zmid = (Zup+Zlow)/2.0
       goto 200
     else
       Zlow = Zmid
       Zmid = (Zup+Zlow)/2.0
       goto 200
     endif
 210 continue
     Xstart = Zmid
!     write(*,*) "Mod=",Mod1,"Ytrough=",Ytrough,"Wave_Number=",wave_number,  &    
!                   "Xstart=",Xstart

     do j = 1,Nloc
       Ein_X0(j) = Ytrough+Amp_Wave*cnoidal_cn(Xstart*0.5*cnoidal_ck(Mod1)+  &
             2.0*cnoidal_ck(Mod1)*(-TIME/Per_Wave),Mod1)**2
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)
     enddo

     ! mean mass transport
     Stokes_Drift = Grav*Amp_Wave**2/(Wave_Length/Per_Wave)/Dep_Wave/8.0
     Fact = 1.0

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)

       Xcn = cnoidal_cn(Xstart*0.5*cnoidal_ck(Mod1)+2.0*cnoidal_ck(Mod1)*(-TIME/Per_Wave),Mod1)
       Xsn = sqrt(1.0-Xcn**2)
       Xdn = sqrt(1.0-Mod1*(1.0-Xcn**2)) 

       Uin_X0(j,k) = sqrt(Grav*Dep_Wave)*(-5./4.+3.*(Dep_Wave+Ytrough)/2.0/Dep_Wave-  &
            (Dep_Wave+Ytrough)**2/4./Dep_Wave**2+(3.*Amp_Wave/2./Dep_Wave-  &                       
            (Dep_Wave+Ytrough)*Amp_Wave/2./Dep_Wave**2)*Xcn**2-Amp_Wave**2/4./  &                      
            Dep_Wave**2*Xcn**4-8.0*Amp_Wave*cnoidal_ck(Mod1)**2/Wave_Length**2*(Dep_Wave/3-  &            
            Zlev**2/2/Dep_Wave)*(-Mod1**2*Xsn**2*Xcn**2+Xcn**2*Xdn**2-Xsn**2*Xdn**2))-  &
            ! substract mean mass transport
            Stokes_Drift*Fact
       Win_X0(j,k) = sqrt(Grav*Dep_Wave)*Zlev*2*Amp_Wave*cnoidal_ck(Mod1)/  &
            Wave_Length/Dep_Wave*(1+(Dep_Wave+Ytrough)/Dep_Wave+Amp_Wave/  &                        
            Dep_Wave*Xcn**2+32*cnoidal_ck(Mod1)**2/3/Wave_Length**2*(Dep_Wave**2-  &
            Zlev**2/2)*(Mod1**2*Xsn**2-Mod1**2*Xcn**2-Xdn**2))*Xsn*Xcn*Xdn
       Vin_X0(j,k) =0.0
     enddo
     enddo

     end subroutine cnoidal_wave_left_boundary


     subroutine cnoidal(height,depth,period,l,c,x2,mod)
!----------------------------------------------------------------
!    Cnoidal Function 
!    Called by
!       cnoidal_wave_left_boundary 
!    Last update: 06/02/2011, Gangfeng Ma 
!----------------------------------------------------------------
     use global, only: SP,Grav
     implicit none

     real(SP), intent(in) :: height, depth, period
     real(SP), intent(out) :: l, c, x2, mod
     real(SP) :: xa,xb,xtemp,cnoidal_cn,cnoidal_ck,cnoidal_ce
     integer :: n

     mod=0.99999999d0
     n=0
 40  n=n+1

     xa=mod*depth+2.0d0*height-mod*height-3.0d0*height*  &
        cnoidal_ce(mod)/cnoidal_ck(mod)-16.0d0*depth**3*mod**2*cnoidal_ck(mod)**2/  &
        3.0d0/Grav/height/period**2
     if(abs(xa).le.1.0e-8.or.n.gt.1000) goto 50
  
     xb=depth-height+3.0d0*height/2.0d0/mod/(1.0d0-mod)/cnoidal_ck(mod)**2*  &
        ((1.0d0-mod)*cnoidal_ck(mod)**2+cnoidal_ce(mod)**2-2.0d0*(1.0d0-mod)* &
        cnoidal_ck(mod)*cnoidal_ce(mod))-16.0d0*depth**3*mod*cnoidal_ck(mod)/3.0d0/Grav/ &
        (1.0d0-mod)/height/period**2*((1.0d0-mod)*cnoidal_ck(mod)+cnoidal_ce(mod))                          

     mod=mod-xa/xb
     goto 40
 50  continue

     ! sobet el at (1987, J. Waterway)                        
     l=4.0*cnoidal_ck(mod)*depth*sqrt(mod*depth/height/3.0)
     ! mei (1983) or simply c=L/T 
     xtemp=-mod+2.0-3.0*cnoidal_ce(mod)/cnoidal_ck(mod)
     c=sqrt(Grav*depth*(1.0+height/depth/mod*xtemp))
     x2=height/mod*(1.0-mod-cnoidal_ce(mod)/cnoidal_ck(mod))
 
     return
     end subroutine cnoidal


     function cnoidal_cn(u,mod)
!----------------------------------------------------------------
!    Cnoidal Function
!----------------------------------------------------------------
     use global, only: SP
     implicit none

     real(SP), intent(in) :: mod,u
     real(SP) ::  mod1,a0,a1,b0,b1,c0,c1,c(1000),a(1000),y(1000),cnoidal_cn
     integer :: n,i

     mod1 = 1.0-mod
     a0 = 1.0
     b0 = sqrt(mod1)
     c0 = sqrt(mod)
     n = 1
     a(n) = a0
     c(n) = c0

 15  if(abs(c0)<1.0e-15.or.n>1000) then
       goto 30
     else
       n = n+1
       a1 = (a0+b0)/2.0
       b1 = sqrt(a0*b0)
       c1 = (a0-b0)/2.0

       a0 = a1
       b0 = b1
       c0 = c1

       a(n) = a0
       c(n) = c0

       goto 15
     endif

 30  y(n) = 2.0**(n-1)*a(n)*u

     do i = n-1,1,-1
       y(i) = (y(i+1)+asin(c(i+1)/a(i+1)*sin(y(i+1))))/2.0
     enddo

     cnoidal_cn = cos(y(1))

     return 
     end function cnoidal_cn

     function cnoidal_ck(mod)
!----------------------------------------------------------------
!    Cnoidal Function
!----------------------------------------------------------------
     use global, only: SP
     implicit none

     real(SP), intent(in) :: mod
     real(SP) :: mod1,a0,a1,b0,b1,c0,c1,cnoidal_ck
     integer :: n

     mod1 = 1.0-mod
     cnoidal_ck = 0.0
     a0 = 1.0
     b0 = dsqrt(mod1)
     c0 = dsqrt(mod)
     n = 1

 15  if(abs(c0)<1.0e-15.or.n>1000) then
       goto 30
     else
       n = n+1
       a1 = (a0+b0)/2.0
       b1 = sqrt(a0*b0)
       c1 = (a0-b0)/2.0

       a0 = a1
       b0 = b1
       c0 = c1

       goto 15
     endif

 30  cnoidal_ck = 3.1415926535897932384626/2.0/a0

     return
     end function cnoidal_ck

     function cnoidal_ce(mod)
!----------------------------------------------------------------
!    Cnoidal Function
!----------------------------------------------------------------
     use global, only: SP
     implicit none

     real(SP), intent(in) :: mod
     real(SP) :: sum,mod1,a0,a1,b0,b1,c0,c1,ck,c(1000),cnoidal_ce
     integer :: n,k1

     cnoidal_ce = 0.0
     mod1 = 1.0-mod
     a0 = 1.0
     b0 = sqrt(mod1)
     c0 = sqrt(mod)
     n = 1
     c(n) = c0

 15  if(abs(c0)<1.0e-15.or.n>1000) then
       goto 30
     else
       n = n+1
       a1 = (a0+b0)/2.0
       b1 = sqrt(a0*b0)
       c1 = (a0-b0)/2.0

       a0 = a1
       b0 = b1
       c0 = c1
       c(n) = c0

       goto 15
     endif

 30  ck = 3.1415926535897932384626/2.0/a0

     sum = 0.0
     do k1 = 1,n
       sum = sum+2.0**(k1-2)*c(k1)**2
     enddo

     cnoidal_ce = ck*(1.0-sum)

     return
     end function cnoidal_ce


     subroutine random_wave_left_boundary
!-----------------------------------------------------    
!    This subroutine is used to specify left boundary 
!    Called by 
!       update_wave_bc  
!    Last update: 04/11/2011, Gangfeng Ma  
!-----------------------------------------------------  
     use global, only: SP,pi,Zero,small,Ibeg,Grav,TIME,Nloc,Kloc,NumFreq,NumDir,Freq,Dire,Wave_Spc2d, & 
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend,sigc,dy, &
                       MaxNumFreq,MaxNumDir,Nglob,Dep_Wave,Random_Phs
     use global, only: myid
     implicit none
     integer :: nfreq,ndir,Iter,nk,i,j,k
     real(SP) :: Per_Wave,Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Angle,tmp1,tmp2,tmp3,tmp4,Fact, &
                 Wtheta(MaxNumDir,MaxNumFreq),Wh(MaxNumDir,MaxNumFreq),Stokes_Drift(MaxNumDir,MaxNumFreq), &
                 Wnum(MaxNumDir,MaxNumFreq),Phs_Lag,dfreq,ddir,Zlev

     ! find the right angle for periodic bc
     do nfreq = 1,NumFreq
     do ndir = 1,NumDir
       Per_Wave = 1.0/Freq(nfreq)
       Segma = 2.0*pi/Per_Wave
       Celerity = sqrt(Grav*Dep_Wave)
       Wave_Length = Celerity*Per_Wave
       Wave_Number = 2.0*pi/Wave_Length

       Iter = 0
 75    Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
       if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 85
         Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &                
            Grav*tanh(Wave_Number*Dep_Wave)
         Wave_Number = Wave_Number-Fk/Fkdif
         Iter = Iter+1
         goto 75
 85    continue
       Wave_Length = 2.0*pi/Wave_Number

       Angle = Dire(ndir)*pi/180.
       goto 100
       if(Angle>zero) then
         tmp3 = zero
         tmp1 = Wave_Number
         nk = 0
         do while (tmp3<Angle)
           nk = nk+1
           tmp2 = nk*2.0*pi/(Nglob*dy)
           if(tmp2>=tmp1) then
             tmp3 = 0.5*pi-small
           else
             tmp3 = asin(tmp2/tmp1)
           endif
         enddo

         ! judge between nk-1 and nk which is closer
         tmp4 = asin((nk-1)*2.0*pi/(Nglob*dy)/tmp1)
         if(abs(tmp4-Angle)<abs(Angle-tmp3)) then
           Angle = tmp4
         else
           Angle = tmp3
         endif
       else
         tmp3 = zero
         tmp1 = Wave_Number
         nk = 0
         do while (tmp3>Angle)
           nk = nk+1
           tmp2 = nk*2.0*pi/(Nglob*dy)
           if(tmp2>=tmp1) then
             tmp3 = -0.5*pi+small
           else
             tmp3 = -asin(tmp2/tmp1)
           endif
         enddo

         ! judge between nk-1 and nk which is closer 
         tmp4= asin((nk-1)*2.0*pi/(Nglob*dy)/tmp1)
         if(abs(tmp4-Angle)<abs(Angle-tmp3)) then
           Angle = tmp4
         else
           Angle = tmp3
         endif
       endif
 100   continue
       Wtheta(ndir,nfreq) = Angle
 
       if(nfreq==1) then
         dfreq = Freq(2)-Freq(1)
       elseif(nfreq==NumFreq) then
         dfreq = Freq(NumFreq)-Freq(NumFreq-1)
       else
         dfreq = 0.5*(Freq(nfreq+1)-Freq(nfreq-1))
       endif
       dfreq = abs(dfreq)

       if(ndir==1) then
         ddir = Dire(2)-Dire(1)
       elseif(ndir==NumDir) then
         ddir = Dire(NumDir)-Dire(NumDir-1)
       else
         ddir = 0.5*(Dire(ndir+1)-Dire(ndir-1))
       endif
       ddir = abs(ddir)
       
       ! save wave number and wave height
       Wnum(ndir,nfreq) = Wave_Number
       Wh(ndir,nfreq) = 2.0*sqrt(2.0*Wave_Spc2d(ndir,nfreq)*dfreq*ddir)

       ! Stokes Drift
       Stokes_Drift(ndir,nfreq) = Grav*Wh(ndir,nfreq)**2/(Wave_Length/Per_Wave)/Dep_Wave/8.0
!       if(myid.eq.0)write(12,'(8E20.10)') Freq(nfreq),Dire(ndir),dfreq,ddir,Wtheta(ndir,nfreq)*180./pi,  &
!               Wnum(ndir,nfreq),Wh(ndir,nfreq),Stokes_Drift(ndir,nfreq)
     enddo
     enddo

     Fact = 1.0
     do j = 1,Nloc
       Ein_X0(j) = 0.0
       do nfreq = 1,NumFreq
       do ndir = 1,NumDir
         Segma = 2.0*pi*Freq(nfreq)
         Phs_Lag = ((j-Jbeg)*dy)*sin(Wtheta(ndir,nfreq))*Wnum(ndir,nfreq)+Random_Phs(ndir,nfreq)
         Ein_X0(j) = Ein_X0(j)+0.5*Wh(ndir,nfreq)*cos(pi/2-Segma*TIME+Phs_Lag)
       enddo
       enddo
       Ein_X0(j)=Ein_X0(j)*tanh(TIME/20.0)
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)
       !if(myid.eq.0) write(12,*) TIME,Ein_X0(Jbeg+Nglob/2)
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)
       Uin_X0(j,k) = 0.0
       Win_X0(j,k) = 0.0
       Vin_X0(j,k) = 0.0
       do nfreq = 1,NumFreq
       do ndir = 1,NumDir
         Segma = 2.0*pi*Freq(nfreq)
         Phs_Lag = ((j-Jbeg)*dy)*sin(Wtheta(ndir,nfreq))*Wnum(ndir,nfreq)+Random_Phs(ndir,nfreq)
         Uin_X0(j,k) = Uin_X0(j,k)+(0.5*Wh(ndir,nfreq)*Segma*cosh(Wnum(ndir,nfreq)*Zlev)/  &
             sinh(Wnum(ndir,nfreq)*Dep_Wave)*cos(pi/2-Segma*TIME+Phs_Lag)-Stokes_Drift(ndir,nfreq)*Fact)*  &
             cos(Wtheta(ndir,nfreq)) 
         Win_X0(j,k) = Win_X0(j,k)+0.5*Wh(ndir,nfreq)*Segma*sinh(Wnum(ndir,nfreq)*Zlev)/  &
             sinh(Wnum(ndir,nfreq)*Dep_Wave)*sin(pi/2-Segma*TIME+Phs_Lag)
         Vin_X0(j,k) = Vin_X0(j,k)+(0.5*Wh(ndir,nfreq)*Segma*cosh(Wnum(ndir,nfreq)*Zlev)/  & 
             sinh(Wnum(ndir,nfreq)*Dep_Wave)*cos(pi/2-Segma*TIME+Phs_Lag)-Stokes_Drift(ndir,nfreq)*Fact)*  & 
             sin(Wtheta(ndir,nfreq))
       enddo
       enddo
     enddo
     enddo

     end subroutine random_wave_left_boundary


     subroutine jonswap_wave_left_boundary
!-----------------------------------------------------                  
!    This subroutine is used to specify left boundary 
!    Called by 
!       update_wave_bc  
!    Last update: 04/11/2011, Gangfeng Ma 
!----------------------------------------------------- 
     use global, only: SP,pi,Zero,small,Ibeg,Grav,TIME,Nloc,Kloc,NumFreq,Jon_Spc, &
                       Freq,Freq_Max,Freq_Min,Dep_Wave,Per_Wave,Hfx,RanPhs, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend,sigc
     implicit none
     real(SP), dimension(NumFreq) :: Wh,Wnum,Stokes_Drift
     real(SP) :: dfreq,Segma,Celerity,Wave_Length,Wave_NUmber,Fk,Fkdif,Fact,Zlev
     integer :: nfreq,Iter,i,j,k

     dfreq = (Freq_Max-Freq_Min)/float(NumFreq)
     do nfreq = 1,NumFreq
       Per_Wave = 1.0/Freq(nfreq)
       Segma = 2.0*pi/Per_Wave
       Celerity = sqrt(Grav*Dep_Wave)
       Wave_Length = Celerity*Per_Wave
       Wave_Number = 2.0*pi/Wave_Length

       Iter = 0
 77    Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
       if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 87
         Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &                 
              Grav*tanh(Wave_Number*Dep_Wave)
         Wave_Number = Wave_Number-Fk/Fkdif
         Iter = Iter+1
         goto 77
 87    continue
       Wave_Length = 2.0*pi/Wave_Number
       Celerity = Wave_Length/Per_Wave

       ! root-mean-square wave height                                                                   
       Wh(nfreq) = 2.0*sqrt(2.0*Jon_Spc(nfreq)*DFreq)
       Wnum(nfreq) = Wave_Number
       Stokes_Drift(nfreq) = Grav*Wh(nfreq)**2/(Wave_Length/Per_Wave)/Dep_Wave/8.0
     enddo

     Fact = 1.0
     do j = 1,Nloc
       Ein_X0(j) = 0.0
       do nfreq = 1,NumFreq
         Segma = 2.0*pi*Freq(nfreq)
         Ein_X0(j) = Ein_X0(j)+0.5*Wh(nfreq)*cos(pi/2-Segma*TIME+RanPhs(nfreq))                           
       enddo
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)
       Uin_X0(j,k) = 0.0
       Win_X0(j,k) = 0.0
       Vin_X0(j,k) = 0.0
       do nfreq = 1,NumFreq
         Segma = 2.0*pi*Freq(nfreq)
         Uin_X0(j,k) = Uin_X0(j,k)+0.5*Wh(nfreq)*Segma*cosh(Wnum(nfreq)*Zlev)/  &                         
             sinh(Wnum(nfreq)*Dep_Wave)*cos(pi/2-Segma*TIME+RanPhs(nfreq))-Stokes_Drift(nfreq)*Fact          
         Win_X0(j,k) = Win_X0(j,k)+0.5*Wh(nfreq)*Segma*sinh(Wnum(nfreq)*Zlev)/  &                            
             sinh(Wnum(nfreq)*Dep_Wave)*sin(pi/2-Segma*TIME+RanPhs(nfreq))
         Vin_X0(j,k) = 0.0
       enddo
     enddo
     enddo

     end subroutine jonswap_wave_left_boundary


     subroutine linear_wave_right_boundary
!-----------------------------------------------------
!    This subroutine is used to specify left boundary 
!    Called by
!       update_wave_bc 
!    Last update: 05/02/2011, Gangfeng Ma
!-----------------------------------------------------
     use global, only: SP,pi,Zero,Iend1,Grav,TIME,Nloc,Kloc,Amp_Wave,Per_Wave,Dep_Wave,Theta_Wave, &
                       Ein_Xn,Din_Xn,Uin_Xn,Vin_Xn,Win_Xn,Hfx,Jbeg,Jend,Kbeg,Kend,sigc,dy
     implicit none
     integer  :: j,k,Iter
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Zlev,Stokes_Drift,Fact,Phs_Lag(Nloc)

     ! Find wave number for linear wave (Newton-Ralphson Method)
     Segma = 2.0*pi/Per_Wave
     Celerity = sqrt(Grav*Dep_Wave)
     Wave_Length = Celerity*Per_Wave
     Wave_Number = 2.0*pi/Wave_Length
     
     Iter = 0
 55  Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
     if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 65
     Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &
        Grav*tanh(Wave_Number*Dep_Wave) 
     Wave_Number = Wave_Number-Fk/Fkdif
     Iter = Iter+1
     goto 55
 65  continue
     Wave_Length = 2.0*pi/Wave_Number

     Stokes_Drift = Grav*Amp_Wave**2/(Wave_Length/Per_Wave)/Dep_Wave/8.0
     Fact = 1.0     

     do j = 1,Nloc
       Phs_Lag(j) = ((j-Jbeg)*dy)*sin(Theta_Wave*pi/180.)*Wave_Number
       Ein_Xn(j) = 0.5*Amp_Wave*cos(pi/2-Segma*TIME+Phs_Lag(j))
       Din_Xn(j) = Ein_Xn(j)+Hfx(Iend1,j)     
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_Xn(j)
       Uin_Xn(j,k) = (0.5*Amp_Wave*Segma*cosh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*cos(pi/2-Segma*TIME+Phs_Lag(j))-Stokes_Drift*Fact)*cos(Theta_Wave*pi/180.)
       Win_Xn(j,k) = 0.5*Amp_Wave*Segma*sinh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*sin(pi/2-Segma*TIME+Phs_Lag(j))
       Vin_Xn(j,k) = (0.5*Amp_Wave*Segma*cosh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*cos(pi/2-Segma*TIME+Phs_Lag(j))-Stokes_Drift*Fact)*sin(Theta_Wave*pi/180.)
     enddo
     enddo         

     end subroutine linear_wave_right_boundary


     subroutine linear_wave_left_boundary
!-----------------------------------------------------
!    This subroutine is used to specify left boundary 
!    Called by
!       update_wave_bc 
!    Last update: 05/02/2011, Gangfeng Ma
!-----------------------------------------------------
     use global
     implicit none
     integer  :: j,k,Iter
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Zlev,Stokes_Drift,Fact,Phs_Lag(Nloc)

     ! Find wave number for linear wave (Newton-Ralphson Method)
     Segma = 2.0*pi/Per_Wave
     Celerity = sqrt(Grav*Dep_Wave)
     Wave_Length = Celerity*Per_Wave
     Wave_Number = 2.0*pi/Wave_Length
     
     Iter = 0
55   Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
     if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 65
     Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &
        Grav*tanh(Wave_Number*Dep_Wave) 
     Wave_Number = Wave_Number-Fk/Fkdif
     Iter = Iter+1
     goto 55
65   continue
     Wave_Length = 2.0*pi/Wave_Number

     Stokes_Drift = Grav*Amp_Wave**2/(Wave_Length/Per_Wave)/Dep_Wave/8.0
     Fact = 1.0     

     do j = 1,Nloc
       Phs_Lag(j) = ((j-Jbeg)*dy)*sin(Theta_Wave*pi/180.)*Wave_Number
       Ein_X0(j) = 0.5*Amp_Wave*cos(pi/2-Segma*TIME+Phs_Lag(j))
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)     
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)
       Uin_X0(j,k) = (0.5*Amp_Wave*Segma*cosh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*cos(pi/2-Segma*TIME+Phs_Lag(j))-Stokes_Drift*Fact)*cos(Theta_Wave*pi/180.)
       Win_X0(j,k) = 0.5*Amp_Wave*Segma*sinh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*sin(pi/2-Segma*TIME+Phs_Lag(j))
       Vin_X0(j,k) = (0.5*Amp_Wave*Segma*cosh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*cos(pi/2-Segma*TIME+Phs_Lag(j))-Stokes_Drift*Fact)*sin(Theta_Wave*pi/180.)
     enddo
     enddo         


     end subroutine linear_wave_left_boundary


     subroutine solitary_wave_left_boundary
!------------------------------------------------------
!    This subroutine is used to specify left boundary
!    Called by
!       update_bc
!    Last update: 22/12/2010, Gangfeng Ma 
!-----------------------------------------------------
     use global, only: SP,Zero,pi,Ibeg,Grav,TIME,Nloc,Kloc,Amp_Wave,Dep_Wave, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Kbeg,Kend,sigc,BC_X0
     implicit none
     integer  :: j,k
     real(SP) :: Celerity,Atmp,Xstart,C2,D1,D2,D3,Zlev,Xc

     Uin_X0 = Zero
     Vin_X0 = Zero
     Win_X0 = Zero
     
     Celerity = sqrt(Grav*Dep_Wave*(1.0+Amp_Wave/Dep_Wave))
     do j = 1,Nloc
       Atmp = sqrt(0.75*Amp_Wave/Dep_Wave**3)
       Xstart = 4.0*Dep_Wave/sqrt(Amp_Wave/Dep_Wave)
       Xc = Xstart-Celerity*TIME
       Ein_X0(j) = Amp_Wave/cosh(Atmp*Xc)**2
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)
 
       C2 = sqrt(Grav*Dep_Wave)
       D2 = 2.0*Amp_Wave*Atmp**2*(2.0*cosh(Atmp*Xc)**2-3.)/cosh(Atmp*Xc)**4
       D1 = -2.0*Amp_Wave*sinh(Atmp*Xc)*Atmp/cosh(Atmp*Xc)**3
       D3 = -8.0*Amp_Wave*sinh(Atmp*Xc)*Atmp**3*(cosh(Atmp*Xc)**2-3.)/cosh(Atmp*Xc)**5
       do k = Kbeg,Kend
         Zlev = sigc(k)*Din_X0(j)
         Uin_X0(j,k) = C2*Ein_X0(j)/Dep_Wave*(1.0-1.0/4.0*Ein_X0(j)/Dep_Wave+Dep_Wave/3.0*(Dep_Wave/Ein_X0(j))*  &
                           (1.0-3.0/2.0*Zlev**2/Dep_Wave**2)*D2)
         Win_X0(j,k) = -C2*Zlev/Dep_Wave*((1.0-1.0/2.0*Ein_X0(j)/Dep_Wave)*D1+1.0/3.0*Dep_Wave**2*  &
                           (1.0-1.0/2.0*Zlev**2/Dep_Wave**2)*D3)
         Vin_X0(j,k) = 0.0
       enddo
     enddo

     end subroutine solitary_wave_left_boundary

 
     subroutine read_bathymetry
!------------------------------------------------------ 
!    This subroutine is used to read bathymetry                                                
!    Called by  
!       main    
!    Last update: 21/12/2010, Gangfeng Ma      
!-----------------------------------------------------
     use global
     implicit none
     integer :: i,j,m,n,iter,iglob,jglob
     integer :: Maskp(Mglob+1,Nglob+1)
     real(SP), dimension(Mglob+1,Nglob+1) :: HG


     ! read bathymetry at grid points
     if(trim(adjustl(DEPTH_TYPE))=='CELL_GRID') then
       if(ANA_BATHY) then
         do j = 1,Nglob+1
         do i = 1,Mglob+1
           HG(i,j) = 30.0
         enddo
         enddo
       else
         open(5,file='depth.txt',status='old')
         do j = 1,Nglob+1
           read(5,*) (HG(i,j),i=1,Mglob+1)
         enddo
       endif

       ! find pernament dry points
       Maskp = 1
       do j = 1,Nglob+1
       do i = 1,Mglob+1
         if(HG(i,j)<-1000.0) Maskp(i,j) = 0
       enddo
       enddo
 
       ! interpolate depth into cell center
       do j = 1,Nglob
       do i = 1,Mglob
         HCG(i,j) = (HG(i,j)*Maskp(i,j)+HG(i+1,j)*Maskp(i+1,j)+  &
             HG(i,j+1)*Maskp(i,j+1)+HG(i+1,j+1)*Maskp(i+1,j+1))/  &
             (Maskp(i,j)+Maskp(i+1,j)+Maskp(i,j+1)+Maskp(i+1,j+1)+1.e-16)
       enddo
       enddo

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         iglob = npx*(Mloc-2*Nghost)+i-Nghost
         jglob = npy*(Nloc-2*Nghost)+j-Nghost
         Hc(i,j) = HCG(iglob,jglob)
       enddo
       enddo

     elseif(trim(adjustl(DEPTH_TYPE))=='CELL_CENTER') then
       ! read bathymetry at cell center
       if(ANA_BATHY) then
         do j = 1,Nglob
         do i = 1,Mglob
           HCG(i,j) = 0.5
         enddo
         enddo
       else ! not analytical bathymetry, read from depth file
         open(5,file='depth.txt',status='old')
         do j = 1,Nglob
           read(5,*) (HCG(i,j),i=1,Mglob)
         enddo

         do j = Jbeg,Jend
         do i = Ibeg,Iend
           iglob = npx*(Mloc-2*Nghost)+i-Nghost
           jglob = npy*(Nloc-2*Nghost)+j-Nghost
           Hc(i,j) = HCG(iglob,jglob)
         enddo
         enddo
       endif
     endif

     ! collect data into ghost cells 
     call phi_2D_coll(Hc)

     ! save the initial water depth
     Hc0 = Hc



     ! find pernament dry cells
     Mask_Struct = 1
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Hc(i,j)<-10.0) then
         Mask_Struct(i,j) = 0
       endif
     enddo
     enddo

     ! reconstruct depth at x-y faces
     do j = 1,Nloc
     do i = 2,Mloc
       Hfx(i,j) = 0.5*(Hc(i-1,j)+Hc(i,j))
     enddo
     Hfx(1,j) = Hc(1,j)
     Hfx(Mloc1,j) = Hc(Mloc,j)
     enddo

     do i = 1,Mloc
     do j = 2,Nloc
       Hfy(i,j) = 0.5*(Hc(i,j-1)+Hc(i,j))
     enddo
     Hfy(i,1) = Hc(i,1)
     Hfy(i,Nloc1) = Hc(i,Nloc)
     enddo

     ! derivatives of water depth at cell center
     do j = 1,Nloc
     do i = 1,Mloc
       DelxH(i,j) = (Hfx(i+1,j)-Hfx(i,j))/dx
       DelyH(i,j) = (Hfy(i,j+1)-Hfy(i,j))/dy
     enddo
     enddo

     end subroutine read_bathymetry

    
     subroutine generate_grid
!------------------------------------------------------
!    This subroutine is used to generate grids
!    Called by
!       main
!    Last update: 20/12/2010, Gangfeng Ma
!-----------------------------------------------------
     use global
     implicit none
     integer :: i,j,k

     ! horizontal grid
     x(Ibeg) = npx*(Mloc-2*Nghost)*dx
     do i = Ibeg+1,Mloc1
       x(i) = x(i-1)+dx
       xc(i-1) = x(i-1)+0.5*dx
     enddo
     do i = Ibeg-1,Ibeg-Nghost,-1
       x(i) = x(i+1)-dx
       xc(i) = x(i+1)-0.5*dx
     enddo

     y(Jbeg) = npy*(Nloc-2*Nghost)*dy
     do j = Jbeg+1,Nloc1
       y(j) = y(j-1)+dy
       yc(j-1) = y(j-1)+0.5*dy
     enddo
     do j = Jbeg-1,Jbeg-Nghost,-1
       y(j) = y(j+1)-dy
       yc(j) = y(j+1)-0.5*dy
     enddo

     ! vertical grid
     if(Ivgrd==1) then
       do k = 1,Kloc
         dsig(k) = 1.0/float(Kglob)
       enddo
     elseif(Ivgrd==2) then
       dsig(Kbeg) = (Grd_R-1.0)/(Grd_R**float(Kglob)-1.0)
       do k = Kbeg+1,Kend
         dsig(k) = dsig(k-1)*Grd_R
       enddo
       do k = 1,Nghost
         dsig(Kbeg-k) = dsig(Kbeg+k-1)
       enddo
       do k = 1,Nghost
         dsig(Kend+k) = dsig(Kend-k+1)         
       enddo
     endif

     sig(Kbeg) = Zero
     do k = Kbeg+1,Kloc1
       sig(k) = sig(k-1)+dsig(k-1)
       sigc(k-1) = sig(k-1)+0.5*dsig(k-1)
     enddo
     do k = Kbeg-1,1,-1
       sig(k) = sig(k+1)-dsig(k)
       sigc(k) = sig(k+1)-0.5*dsig(k)
     enddo

     end subroutine generate_grid

     
     subroutine initial
!---------------------------------------------------  
!    This subroutine is used to initialize model run 
!    Called by                                       
!       main 
!    Last update: 21/12/2010, Gangfeng Ma  
!---------------------------------------------------
     use global
     implicit none
     integer  :: i,j,k,n,m,nmax,iglob,jglob
     real(SP) :: xsol(80),zsol(80),zmax,xmax,xterp,zterp,tmp,zc(Kloc), &
                 utmp1,wtmp1,utmp2,wtmp2,xk(321,16),zk(321,16),  &
                 uk(321,16),wk(321,16)
     real(SP) :: Ufric,Zlev1,Zlev,mud_dens,Zslide(Mloc,Nloc),conc_slide, &
                 alpha0,L0,T,bl,hs0,ls0,ls1,ls2,lsx,lsx1,hslide,  &                   
                 eslide,xt,yt,zt,kb,kw,Slope,Xslide,SlideX1,SlideX2
     real(SP), dimension(Mglob,Nglob) :: EtaG
     real(SP), dimension(Mglob,Nglob,Kglob) :: UG,VG,WG,SaliG

     ! simulation time
     TIME = Zero
     RUN_STEP = 0
     dt = dt_ini     
     Screen_Count = Zero
     Plot_Count = Zero
     Plot_Count_Stat = Zero
     Icount = 0
     
     ! working arrays
     D = Zero
     U = Zero
     V = Zero
     W = Zero
     P = Zero
     Omega = Zero
     DU = Zero
     DV = Zero
     DW = Zero
     D0 = Zero
     Eta0 = Zero
     DU0 = Zero
     DV0 = Zero
     DW0 = Zero
     Uf = Zero
     Vf = Zero
     Wf = Zero
     Rho = Rho0
     
     ! source terms
     SourceC = Zero
     SourceX = Zero
     SourceY = Zero

     ! fluxes
     DxL = Zero
     DxR = Zero
     DyL = Zero
     DyR = Zero
     UxL = Zero
     UxR = Zero
     UyL = Zero
     UyR = Zero
     UzL = Zero
     UzR = Zero
     VxL = Zero
     VxR = Zero
     VyL = Zero
     VyR = Zero
     VzL = Zero
     VzR = Zero
     WxL = Zero
     WxR = Zero
     WyL = Zero
     WyR = Zero
     WzL = Zero
     WzR = Zero
     DUxL = Zero
     DUxR = Zero
     DUyL = Zero
     DUyR = Zero
     DVxL = Zero
     DVxR = Zero
     DVyL = Zero
     DVyR = Zero
     DWxL = Zero
     DWxR = Zero
     DWyL = Zero
     DWyR = Zero
     OzL = Zero
     OzR = Zero
     SxL = Zero
     SxR = Zero
     SxS = Zero
     SyL = Zero
     SyR = Zero
     SyS = Zero
     Ex = Zero
     Ey = Zero
     Fx = Zero
     Fy = Zero
     Fz = Zero
     Gx = Zero
     Gy = Zero
     Gz = Zero
     Hx = Zero
     Hy = Zero
     Hz = Zero
     EtaxL = Zero
     EtaxR = Zero
     EtayL = Zero
     EtayR = Zero
     DelxEta = Zero
     DelyEta = Zero
     DeltH = Zero
     DeltHo = Zero
     Delt2H = Zero
     DelxD = Zero
     DelyD = Zero
     DelxU = Zero
     DelyU = Zero
     DelzU = Zero
     DelxV = Zero
     DelyV = Zero
     DelzV = Zero
     DelxW = Zero
     DelyW = Zero
     DelzW = Zero
     DelxDU = Zero
     DelyDU = Zero
     DelxDV = Zero
     DelyDV = Zero
     DelzO = Zero
     Sponge = One
     Cmu = Visc
     CmuHt = Zero
     CmuVt = Zero
     CmuR  = Zero
     Richf = Zero

     Diffxx = Zero
     Diffxy = Zero
     Diffxz = Zero
     Diffyx = Zero
     Diffyy = Zero
     Diffyz = Zero
     Diffzx = Zero
     Diffzy = Zero
     Diffzz = Zero

     Uin_X0 = Zero
     Vin_X0 = Zero
     Win_X0 = Zero
     Ein_X0 = Zero
     Din_X0 = Zero
     Uin_Xn = Zero
     Vin_Xn = Zero
     Win_Xn = Zero
     Ein_Xn = Zero
     Din_Xn = Zero   

     Setup = Zero
     WaveHeight = Zero
     Umean = Zero
     Vmean = Zero
     Num_Zero_Up = 0
     Emax = -1000.
     Emin = 1000.

     WdU = Zero
     WdV = Zero

     Lag_Umean = Zero
     Lag_Vmean = Zero
     Lag_Wmean = Zero
     Euler_Umean = Zero
     Euler_Vmean = Zero
     Euler_Wmean = Zero


     ExtForceX = Zero
     ExtForceY = Zero

     ! baroclinic terms
     DRhoX = Zero
     DRhoY = Zero








     Tke = Zero
     Eps = Zero
     DTke = Zero
     DEps = Zero
     DTke0 = Zero
     DEps0 = Zero

     ! wave breaking mask
     Brks = 0

     ! pressure boundary
     Bc_Prs = Zero

     ! initial surface elevation (user-specified)
     Eta = Zero

     if(INITIAL_EUVW) then

        ! initial condition for eta
        open(21,file='eta0.txt',status='old')
        do j = 1,Nglob
          read(21,*) (EtaG(i,j),i=1,Mglob)
        enddo
      
        ! initial condition for U, V, W
        open(22,file='uvw0.txt',status='old')
        do k = 1,Kglob
        do j = 1,Nglob
          read(22,*) (UG(i,j,k),i=1,Mglob)
        enddo
        enddo

        do k = 1,Kglob
        do j = 1,Nglob
          read(22,*) (VG(i,j,k),i=1,Mglob)
        enddo
        enddo

        do k = 1,Kglob
        do j = 1,Nglob
          read(22,*) (WG(i,j,k),i=1,Mglob)
        enddo
        enddo

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         iglob = npx*(Mloc-2*Nghost)+i-Nghost
         jglob = npy*(Nloc-2*Nghost)+j-Nghost
         Eta(i,j) = EtaG(iglob,jglob)
         do k = Kbeg,Kend
           U(i,j,k) = UG(iglob,jglob,k-Nghost)
           V(i,j,k) = VG(iglob,jglob,k-Nghost)
           W(i,j,k) = WG(iglob,jglob,k-Nghost)
         enddo
       enddo
       enddo

!       ! solitary wave from Tanaka solution
!       open(21,file='soliton.dat')
!       do n = 1,80
!         read(21,*) i,xsol(n),zsol(n),tmp,tmp
!       enddo
!       close(21)
!
!       ! find the peak location                                                   
!       zmax = -1.0e+10
!       do n = 1,80
!         if(zsol(n)>zmax) then
!           zmax = zsol(n)
!           xmax = xsol(n)
!         endif
!       enddo
!
!       ! move the peak to x = 3.0m                                                
!       do n = 1,80
!         xsol(n) = xsol(n)+8.0-xmax
!       enddo
!
!       ! interpolate into computational grid
!       do j = Jbeg,Jend
!       do i = Ibeg,Iend
!         if(xc(i)>xsol(80).and.xc(i)<xsol(1)) cycle
!         do n = 2,80
!           if(xc(i)>=xsol(n-1).and.xc(i)<xsol(n)) then
!             xterp = (xc(i)-xsol(n-1))/(xsol(n)-xsol(n-1))
!             Eta(i,j) = (1.0-xterp)*zsol(n-1)+xterp*zsol(n)
!           endif
!         enddo
!       enddo
!       enddo
!
!       open(22,file='plotuv.dat')
!       do n = 1,16
!       do m = 1,321
!         read(22,*) k,xk(m,n),zk(m,n),tmp,uk(m,n),wk(m,n)
!         xk(m,n) = xk(m,n)+8.0
!       enddo
!       enddo
!       close(22)
!       
!       do j = Jbeg,Jend
!       do k = Kbeg,Kend
!       do i = Ibeg,Iend
!         if(xc(i)>xsol(80).and.xc(i)<xsol(1)) cycle
!         zc(k) = (1.0+Eta(i,j))*sigc(k)-1.0
!
!         do n = 2,16
!         do m = 2,321
!           if(xc(i)>=xk(m-1,n).and.xc(i)<xk(m,n).and.  &
!                    zc(k)>=zk(m,n-1).and.zc(k)<zk(m,n)) then
!             xterp = (xc(i)-xk(m-1,n))/(xk(m,n)-xk(m-1,n))
!             zterp = (zc(k)-zk(m,n-1))/(zk(m,n)-zk(m,n-1)) 
!             utmp1 = (1.0-xterp)*uk(m-1,n-1)+xterp*uk(m,n-1)
!             wtmp1 = (1.0-xterp)*wk(m-1,n-1)+xterp*wk(m,n-1)
!             utmp2 = (1.0-xterp)*uk(m-1,n)+xterp*uk(m,n)
!             wtmp2 = (1.0-xterp)*wk(m-1,n)+xterp*wk(m,n)
!             
!             U(i,j,k) = (1.0-zterp)*utmp1+zterp*utmp2
!             W(i,j,k) = (1.0-zterp)*wtmp1+zterp*wtmp2
!           endif
!         enddo
!         enddo
!       enddo
!       enddo
!       enddo
! 100   continue
     endif    

     ! collect data into ghost cells
     call phi_2D_coll(Eta)
     Eta0 = Eta

     ! wetting-drying mask
     ! Mask: 1 - wet; 0 - dry
     ! Mask_Struct: 0 - permanent dry point
     ! Mask9: mask for itself and 8 elements around
     Mask = 1
     do j = 1,Nloc
     do i = 1,Mloc
       if((Eta(i,j)+Hc(i,j))<=MinDep) then
         Mask(i,j) = 0
         Eta(i,j) = MinDep-Hc(i,j)
       else
         Mask(i,j) = 1
       endif
     enddo
     enddo
     Mask = Mask*Mask_Struct

     ! collect mask into ghost cells
     call phi_int_exch(Mask)

     do j = Jbeg,Jend
     do i = Ibeg,Iend
      Mask9(i,j) = Mask(i,j)*Mask(i-1,j)*Mask(i+1,j)  &
                *Mask(i+1,j+1)*Mask(i,j+1)*Mask(i-1,j+1) &
                *Mask(i+1,j-1)*Mask(i,j-1)*Mask(i-1,j-1)
     enddo
     enddo

     ! total water depth and flux
     D = max(Hc+Eta, MinDep)

     call vel_bc
     call phi_3D_exch(U)
     call phi_3D_exch(V)
     call phi_3D_exch(W)

     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc
       DU(i,j,k) = D(i,j)*U(i,j,k)*Mask(i,j)
       DV(i,j,k) = D(i,j)*V(i,j,k)*Mask(i,j)
       DW(i,j,k) = D(i,j)*W(i,j,k)*Mask(i,j)
     enddo
     enddo
     enddo

     if(VISCOUS_FLOW) then
       ! initial seeding values for turbulence
!       Tke_min = 0.5*(1.4e-3)**2
!       Eps_min = 0.09*Tke_min**2/(0.1*Visc)
!       Tke_min = 1.e-12
!       Eps_min = 0.09*Tke_min**2/(1.e-4*Visc)
       Tke_min = 1.e-9
       Eps_min = 1.e-9
       Cmut_min = 9.e-2*Tke_min**2/Eps_min     
       do k = 1,Kloc
       do j = 1,Nloc
       do i = 1,Mloc
         Tke(i,j,k) = Tke_min
         Eps(i,j,k) = Eps_min
         CmuHt(i,j,k) = Cmut_min
         CmuVt(i,j,k) = Cmut_min
         DTke(i,j,k) = D(i,j)*Tke(i,j,k)*Mask(i,j)
         DEps(i,j,k) = D(i,j)*Eps(i,j,k)*Mask(i,j)
       enddo
       enddo
       enddo
     endif



     ! SSP Runge-Kutta method parameters
     if(TIME_ORDER(1:3)=='THI') then
       It_Order = 3
       ALPHA(1) = 0.0
       ALPHA(2) = 3.0/4.0
       ALPHA(3) = 1.0/3.0
       BETA(1) = 1.0
       BETA(2) = 1.0/4.0
       BETA(3) = 2.0/3.0
     elseif(TIME_ORDER(1:3)=='SEC') then
       It_Order = 2
       ALPHA(1) = 0.0
       ALPHA(2) = 1.0/2.0
       BETA(1) = 1.0
       BETA(2) = 1.0/2.0
     else
       It_Order = 1
       ALPHA(1) = 0.0
       BETA(1) = 1.0
     endif

     ! sponge layer
     if(SPONGE_ON) then
       call calculate_sponge
     endif
    
     end subroutine initial

     subroutine allocate_variables 
!--------------------------------------------------- 
!    This subroutine is used to allocate variables
!    Called by                                    
!       main                                                                                  
!    Last update: 23/12/2010, Gangfeng Ma                                                   
!---------------------------------------------------
     use global
     implicit none

     ! one-dimensional vars
     ALLOCATE(x(Mloc1),xc(Mloc),y(Nloc1),yc(Nloc),sig(Kloc1),dsig(Kloc),sigc(Kloc),  &
              Ein_X0(Nloc),Din_X0(Nloc),Ein_Xn(Nloc),Din_Xn(Nloc))

     ! two-dimensional vars
     ALLOCATE(HCG(Mglob,Nglob),Ho(Mloc,Nloc),H(Mloc1,Nloc1),Hc(Mloc,Nloc),Hc0(Mloc,Nloc), &
              Hfx(Mloc1,Nloc),Hfy(Mloc,Nloc1),D(Mloc,Nloc),Eta0(Mloc,Nloc),Eta00(Mloc,Nloc),  &
              D0(Mloc,Nloc),DeltH(Mloc,Nloc),DelxH(Mloc,Nloc),DelyH(Mloc,Nloc),Eta(Mloc,Nloc),Mask(Mloc,Nloc),  &
              Mask_Struct(Mloc,Nloc),Mask9(Mloc,Nloc),SourceC(Mloc,Nloc),SourceX(Mloc,Nloc), &
              SourceY(Mloc,Nloc),DeltHo(Mloc,Nloc),Uin_X0(Nloc,Kloc),Vin_X0(Nloc,Kloc),  &
              Win_X0(Nloc,Kloc),Uin_Xn(Nloc,Kloc),Delt2H(Mloc,Nloc),  &
              Vin_Xn(Nloc,Kloc),Win_Xn(Nloc,Kloc),Bc_Prs(Mloc,Nloc),Brks(Mloc,Nloc))
     ALLOCATE(DxL(Mloc1,Nloc),DxR(Mloc1,Nloc),DyL(Mloc,Nloc1),DyR(Mloc,Nloc1), &
              EtaxL(Mloc1,Nloc),EtaxR(Mloc1,Nloc),EtayL(Mloc,Nloc1),EtayR(Mloc,Nloc1), &
              DelxEta(Mloc,Nloc),DelyEta(Mloc,Nloc),DelxD(Mloc,Nloc),DelyD(Mloc,Nloc),Sponge(Mloc,Nloc), &
              Setup(Mloc,Nloc),WaveHeight(Mloc,Nloc),Umean(Mloc,Nloc),Vmean(Mloc,Nloc),Num_Zero_Up(Mloc,Nloc), &
              Emax(Mloc,Nloc),Emin(Mloc,Nloc),WdU(Mloc,Nloc),WdV(Mloc,Nloc),Wsx(Mloc,Nloc),Wsy(Mloc,Nloc))

     ! three-dimensional vars
     ALLOCATE(U(Mloc,Nloc,Kloc),V(Mloc,Nloc,Kloc),W(Mloc,Nloc,Kloc),Omega(Mloc,Nloc,Kloc1), &
              P(Mloc,Nloc,Kloc1),DU(Mloc,Nloc,Kloc),DV(Mloc,Nloc,Kloc),DW(Mloc,Nloc,Kloc),  &
              U0(Mloc,Nloc,Kloc),V0(Mloc,Nloc,Kloc),W0(Mloc,Nloc,Kloc),  &
              U00(Mloc,Nloc,Kloc),V00(Mloc,Nloc,Kloc),W00(Mloc,Nloc,Kloc),  &
              DU0(Mloc,Nloc,Kloc),DV0(Mloc,Nloc,Kloc),DW0(Mloc,Nloc,Kloc),Uf(Mloc,Nloc,Kloc1), &
              Vf(Mloc,Nloc,Kloc1),Wf(Mloc,Nloc,Kloc1),Cmu(Mloc,Nloc,Kloc),CmuR(Mloc,Nloc,Kloc), &
              Diffxx(Mloc,Nloc,Kloc),Diffxy(Mloc,Nloc,Kloc),Diffxz(Mloc,Nloc,Kloc), &
              Diffyx(Mloc,Nloc,Kloc),Diffyy(Mloc,Nloc,Kloc),Diffyz(Mloc,Nloc,Kloc),Diffzx(Mloc,Nloc,Kloc), &
              Diffzy(Mloc,Nloc,Kloc),Diffzz(Mloc,Nloc,Kloc),DelxSc(Mloc,Nloc,Kloc),DelySc(Mloc,Nloc,Kloc), &
              CmuHt(Mloc,Nloc,Kloc),CmuVt(Mloc,Nloc,Kloc),Rho(Mloc,Nloc,Kloc),Rmean(Mloc,Nloc,Kloc),Tke(Mloc,Nloc,Kloc), &
              Eps(Mloc,Nloc,Kloc),Skl(Mloc,Nloc,Kloc),DTke(Mloc,Nloc,Kloc),DEps(Mloc,Nloc,Kloc),DTke0(Mloc,Nloc,Kloc), &
              DEps0(Mloc,Nloc,Kloc),Prod_s(Mloc,Nloc,Kloc),Prod_b(Mloc,Nloc,Kloc),Lag_Umean(Mloc,Nloc,Kloc), &
              Lag_Vmean(Mloc,Nloc,Kloc),Lag_Wmean(Mloc,Nloc,Kloc),Euler_Umean(Mloc,Nloc,Kloc),Euler_Vmean(Mloc,Nloc,Kloc), &
              Euler_Wmean(Mloc,Nloc,Kloc),DRhoX(Mloc,Nloc,Kloc),DRhoY(Mloc,Nloc,Kloc),ExtForceX(Mloc,Nloc,Kloc), &
              ExtForceY(Mloc,Nloc,Kloc),UpWp(Mloc,Nloc,Kloc),IsMove(Mloc,Nloc,Kloc),Richf(Mloc,Nloc,Kloc), &
              PresForceX(Mloc,Nloc,Kloc),PresForceY(Mloc,Nloc,Kloc),PresForceZ(Mloc,Nloc,Kloc))

     ! fluxes for construction at cell faces    
     ALLOCATE(UxL(Mloc1,Nloc,Kloc),UxR(Mloc1,Nloc,Kloc),VxL(Mloc1,Nloc,Kloc),VxR(Mloc1,Nloc,Kloc), &
              WxL(Mloc1,Nloc,Kloc),WxR(Mloc1,Nloc,Kloc),DUxL(Mloc1,Nloc,Kloc),DUxR(Mloc1,Nloc,Kloc), &
              DVxL(Mloc1,Nloc,Kloc),DVxR(Mloc1,Nloc,Kloc),DWxL(Mloc1,Nloc,Kloc),DWxR(Mloc1,Nloc,Kloc), &
              UyL(Mloc,Nloc1,Kloc),UyR(Mloc,Nloc1,Kloc),VyL(Mloc,Nloc1,Kloc),VyR(Mloc,Nloc1,Kloc), &
              WyL(Mloc,Nloc1,Kloc),WyR(Mloc,Nloc1,Kloc),DUyL(Mloc,Nloc1,Kloc),DUyR(Mloc,Nloc1,Kloc), &
              DVyL(Mloc,Nloc1,Kloc),DVyR(Mloc,Nloc1,Kloc),DWyL(Mloc,Nloc1,Kloc),DWyR(Mloc,Nloc1,Kloc), &
              UzL(Mloc,Nloc,Kloc1),UzR(Mloc,Nloc,Kloc1),VzL(Mloc,Nloc,Kloc1),VzR(Mloc,Nloc,Kloc1), &
              WzL(Mloc,Nloc,Kloc1),WzR(Mloc,Nloc,Kloc1),OzL(Mloc,Nloc,Kloc1),OzR(Mloc,Nloc,Kloc1), &
              SxL(Mloc1,Nloc,Kloc),SxR(Mloc1,Nloc,Kloc),SxS(Mloc1,Nloc,Kloc), &
              SyL(Mloc,Nloc1,Kloc),SyR(Mloc,Nloc1,Kloc),SyS(Mloc,Nloc1,Kloc), &
              ExL(Mloc1,Nloc,Kloc),ExR(Mloc1,Nloc,Kloc),FxL(Mloc1,Nloc,Kloc),FxR(Mloc1,Nloc,Kloc), &
              GxL(Mloc1,Nloc,Kloc),GxR(Mloc1,Nloc,Kloc),HxL(Mloc1,Nloc,Kloc),HxR(Mloc1,Nloc,Kloc), &
              EyL(Mloc,Nloc1,Kloc),EyR(Mloc,Nloc1,Kloc),FyL(Mloc,Nloc1,Kloc),FyR(Mloc,Nloc1,Kloc), &
              GyL(Mloc,Nloc1,Kloc),GyR(Mloc,Nloc1,Kloc),HyL(Mloc,Nloc1,Kloc),HyR(Mloc1,Nloc1,Kloc), &
              Ex(Mloc1,Nloc,Kloc),Ey(Mloc,Nloc1,Kloc),Fx(Mloc1,Nloc,Kloc),Fy(Mloc,Nloc1,Kloc), &
              Gx(Mloc1,Nloc,Kloc),Gy(Mloc,Nloc1,Kloc),Hx(Mloc1,Nloc,Kloc),Hy(Mloc,Nloc1,Kloc), &
              Fz(Mloc,Nloc,Kloc1),Gz(Mloc,Nloc,Kloc1),Hz(Mloc,Nloc,Kloc1),DelxU(Mloc,Nloc,Kloc), &
              DelyU(Mloc,Nloc,Kloc),DelzU(Mloc,Nloc,Kloc),DelxV(Mloc,Nloc,Kloc),DelyV(Mloc,Nloc,Kloc), &
              DelzV(Mloc,Nloc,Kloc),DelxW(Mloc,Nloc,Kloc),DelyW(Mloc,Nloc,Kloc),DelzW(Mloc,Nloc,Kloc), &
              DelxDU(Mloc,Nloc,Kloc),DelyDU(Mloc,Nloc,Kloc),DelxDV(Mloc,Nloc,Kloc),DelyDV(Mloc,Nloc,Kloc), &
              DelxDW(Mloc,Nloc,Kloc),DelyDW(Mloc,Nloc,Kloc),DelzO(Mloc,Nloc,Kloc)) 









     ! poisson solver (for NSPCG use)
     neqns = (Iend-Ibeg+1)*(Jend-Jbeg+1)*(Kend-Kbeg+1)
     ALLOCATE(Coef(5*neqns,5*15),JCoef(5*15),Rhs(neqns))

     end subroutine allocate_variables
     

     subroutine index
!---------------------------------------------------
!    This subroutine is used to creat work index
!    Called by                   
!       main                                                                
!    Last update: 20/12/2010, Gangfeng Ma                                     
!---------------------------------------------------
     use global
     implicit none

     dims(1)=PX
     dims(2)=PY
     periods(1)=.false.
     periods(2)=.false.
     if(PERIODIC_X) periods(1)=.true.
     if(PERIODIC_Y) periods(2)=.true.
     coords(1)=0
     coords(2)=0

     call MPI_CART_CREATE(MPI_COMM_WORLD,ndims,dims, &
         periods,reorder,comm2d,ier)
     call MPI_CART_COORDS(comm2d,myid,2,coords,ier)

     npx=coords(1)
     npy=coords(2)
 
     call MPI_CART_SHIFT(comm2d,0,1,n_west,n_east,ier)
     call MPI_CART_SHIFT(comm2d,1,1,n_suth,n_nrth,ier)

     ! local index
     Mloc = Mglob/PX+2*Nghost
     Nloc = Nglob/PY+2*Nghost
     Kloc = Kglob+2*Nghost
     Mloc1 = Mloc+1
     Nloc1 = Nloc+1
     Kloc1 = Kloc+1

     Ibeg = Nghost+1
     Iend = Mloc-Nghost
     Iend1 = Mloc1-Nghost
     Jbeg = Nghost+1
     Jend = Nloc-Nghost
     Jend1 = Nloc1-Nghost
     Kbeg = Nghost+1
     Kend = Kloc-Nghost
     Kend1 = Kloc1-Nghost

     end subroutine index


     subroutine read_input
!---------------------------------------------------
!    This subroutine is used to read input.txt
!    Called by 
!       main
!    Last update: 20/12/2010, Gangfeng Ma
!---------------------------------------------------
     use global
     use input_util
     implicit none
     character(len=80) :: FILE_NAME
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,  &
                 Fk,Fkdif,Theta_Calc,Wnumy,tmp,tmp1,DFreq, &
                 Freq_Peak,gam,sa,sb,SumInt,A_Jon, rand
     integer :: line,ierr,Iter,i,j
 
     ! log and error file
     open(3,file='log.txt')

     ! read from input.txt
     FILE_NAME='input.txt'

     ! title
     CALL GET_STRING_VAL(TITLE,FILE_NAME,'TITLE',line,ierr)
     IF(ierr==1)THEN
       if(myid.eq.0) write(3,*) 'No TITLE in ', FILE_NAME, 'use default'
       TITLE='---TEST RUN---'
     ENDIF
     if(myid.eq.0) WRITE(3,*)'---- LOG FILE ---'
     if(myid.eq.0) WRITE(3,*)TITLE
     if(myid.eq.0) WRITE(3,*)'--------------input start --------------'

     ! dimension                                             
     CALL GET_INTEGER_VAL(Mglob,FILE_NAME,'Mglob',line)
     CALL GET_INTEGER_VAL(Nglob,FILE_NAME,'Nglob',line)
     CALL GET_INTEGER_VAL(Kglob,FILE_NAME,'Kglob',line)
     if(myid.eq.0) WRITE(3,'(A7,I5)')'Mglob= ',Mglob
     if(myid.eq.0) WRITE(3,'(A7,I5)')'Nglob= ',Nglob
     if(myid.eq.0) WRITE(3,'(A7,I5)')'Kglob= ',Kglob

     ! processor number
     CALL GET_INTEGER_VAL(PX,FILE_NAME,'PX',line)
     CALL GET_INTEGER_VAL(PY,FILE_NAME,'PY',line)
     if(myid.eq.0) WRITE(3,'(A4,I5)')'PX= ',PX
     if(myid.eq.0) WRITE(3,'(A4,I5)')'PY= ',PY
     if(PX*PY.ne.NumP) then
       if(myid.eq.0) WRITE(3,'(A6,I5)') 'NumP= ',NumP
       stop
     endif

     ! grid sizes
     CALL GET_Float_VAL(dx,FILE_NAME,'DX',line)
     CALL GET_Float_VAL(dy,FILE_NAME,'DY',line)
     if(myid.eq.0) WRITE(3,'(A4,F8.4)')'DX= ',dx
     if(myid.eq.0) WRITE(3,'(A4,F8.4)')'DY= ',dy

     ! vertical grid option
     call GET_INTEGER_VAL(Ivgrd,FILE_NAME,'IVGRD',line)
     CALL GET_Float_VAL(Grd_R,FILE_NAME,'GRD_R',line)
     if(myid.eq.0) WRITE(3,'(A7,I3)')'Ivgrd= ',Ivgrd
     if(myid.eq.0) WRITE(3,'(A7,f5.2)')'Grd_R= ',Grd_R

     ! time step
     CALL GET_Float_VAL(dt_ini,FILE_NAME,'DT_INI',line)
     CALL GET_Float_VAL(dt_min,FILE_NAME,'DT_MIN',line)
     CALL GET_Float_VAL(dt_max,FILE_NAME,'DT_MAX',line)
     if(myid.eq.0) WRITE(3,'(A8,F8.4)')'DT_INI= ',dt_ini
     if(myid.eq.0) WRITE(3,'(A8,F8.4)')'DT_MIN= ',dt_min
     if(myid.eq.0) WRITE(3,'(A8,F8.4)')'DT_MAX= ',dt_max

     ! result folder                                     
     CALL GET_STRING_VAL(RESULT_FOLDER,FILE_NAME,'RESULT_FOLDER',line,ierr)
     if(myid.eq.0) WRITE(3,'(A15,A50)')'RESULT_FOLDER= ', RESULT_FOLDER

     ! simulation steps and time
     call GET_INTEGER_VAL(SIM_STEPS,FILE_NAME,'SIM_STEPS',line)
     CALL GET_Float_VAL(TOTAL_TIME,FILE_NAME,'TOTAL_TIME',line)
     CALL GET_Float_VAL(Plot_Start,FILE_NAME,'PLOT_START',line)
     CALL GET_Float_VAL(Plot_Intv,FILE_NAME,'PLOT_INTV',line)
     CALL GET_Float_VAL(Screen_Intv,FILE_NAME,'SCREEN_INTV',line)
     if(myid.eq.0) WRITE(3,'(A11,I12)')'SIM_STEPS= ', SIM_STEPS
     if(myid.eq.0) WRITE(3,'(A12,F8.2)')'TOTAL_TIME= ', TOTAL_TIME
     if(myid.eq.0) WRITE(3,'(A11,F8.2)')'PLOT_START= ', Plot_Start
     if(myid.eq.0) WRITE(3,'(A11,F8.2)')'PLOT_INTV= ', Plot_Intv
     if(myid.eq.0) WRITE(3,'(A13,F8.2)')'SCREEN_INTV= ', Screen_Intv

     ! courant number
     CALL GET_Float_VAL(CFL,FILE_NAME,'CFL',line)
     if(myid.eq.0) WRITE(3,'(A5,F8.3)')'CFL= ',CFL

     ! viscous number
     CALL GET_Float_VAL(VISCOUS_NUMBER,FILE_NAME,'VISCOUS_NUMBER',line)
     if(myid.eq.0) WRITE(3,'(A16,F8.3)')'VISCOUS_NUMBER= ',VISCOUS_NUMBER

     ! minimum depth
     CALL GET_Float_VAL(MinDep,FILE_NAME,'MinDep',line)
     if(myid.eq.0) WRITE(3,'(A8,F8.3)')'MinDep= ',MinDep

     ! laminar viscosity
     CALL GET_LOGICAL_VAL(VISCOUS_FLOW,FILE_NAME,'VISCOUS_FLOW',line)
     CALL GET_INTEGER_VAL(IVturb,FILE_NAME,'IVTURB',line)
     CALL GET_INTEGER_VAL(IHturb,FILE_NAME,'IHTURB',line)
     CALL GET_Float_VAL(Visc,FILE_NAME,'VISCOSITY',line)
     CALL GET_Float_VAL(Schmidt,FILE_NAME,'Schmidt',line)
     CALL GET_Float_VAL(Cvs,FILE_NAME,'Cvs',line)
     CALL GET_Float_VAL(Chs,FILE_NAME,'Chs',line)
     if(myid.eq.0) WRITE(3,'(A14,L4)')'VISCOUS_FLOW= ',VISCOUS_FLOW
     if(myid.eq.0) WRITE(3,'(A8,I2)')'IVTURB= ',IVturb
     if(myid.eq.0) WRITE(3,'(A8,I2)')'IHTURB= ',IHturb
     if(myid.eq.0) WRITE(3,'(A11,F8.3)')'VISCOSITY= ',Visc
     if(myid.eq.0) WRITE(3,'(A9,F8.3)')'Schmidt= ',Schmidt
     if(myid.eq.0) WRITE(3,'(A6,F8.3)')'Cvs= ',Cvs
     if(myid.eq.0) WRITE(3,'(A6,F8.3)')'Chs= ',Chs

     ! bathymetry     
     CALL GET_STRING_VAL(DEPTH_TYPE,FILE_NAME,'DEPTH_TYPE',line,ierr)
     CALL GET_LOGICAL_VAL(ANA_BATHY,FILE_NAME,'ANA_BATHY',line)
     if(myid.eq.0) WRITE(3,'(A12,A50)')'DEPTH_TYPE= ',DEPTH_TYPE
     if(myid.eq.0) WRITE(3,'(A11,L4)')'ANA_BATHY= ',ANA_BATHY

     ! initial conditions
     CALL GET_LOGICAL_VAL(INITIAL_EUVW,FILE_NAME,'INITIAL_EUVW',line)
     if(myid.eq.0) WRITE(3,'(A14,L4)')'INITIAL_EUVW= ',INITIAL_EUVW

     ! bottom roughness
     CALL GET_INTEGER_VAL(Ibot,FILE_NAME,'Ibot',line)
     CALL GET_Float_VAL(Cd0,FILE_NAME,'Cd0',line)
     CALL GET_Float_VAL(Zob,FILE_NAME,'Zob',line)
     if(myid.eq.0) WRITE(3,'(A6,I2)')'Ibot= ',Ibot
     if(myid.eq.0) WRITE(3,'(A5,F8.5)')'Cd0= ',Cd0
     if(myid.eq.0) WRITE(3,'(A5,F8.5)')'Zob= ',Zob


     ! wind speed/stress
     CALL GET_INTEGER_VAL(Iws,FILE_NAME,'Iws',line)
     if(Iws==1) then
       CALL GET_Float_VAL(WindU,FILE_NAME,'WindU',line)
       CALL GET_Float_VAL(WindV,FILE_NAME,'WindV',line)
     endif

     ! Coriolis
     CALL GET_Float_VAL(slat,FILE_NAME,'slat',line)

     slat = slat*pi/180.0
     fcor = 2.0*7.29e-5*sin(slat)

     ! barotropic or baroclinic
     CALL GET_LOGICAL_VAL(BAROTROPIC,FILE_NAME,'BAROTROPIC',line)
     if(myid.eq.0) WRITE(3,'(A12,L4)')'BAROTROPIC= ',BAROTROPIC

     ! numerical scheme
     CALL GET_STRING_VAL(HIGH_ORDER,FILE_NAME,'HIGH_ORDER',line,ierr)
     CALL GET_STRING_VAL(TIME_ORDER,FILE_NAME,'TIME_ORDER',line,ierr)
     CALL GET_STRING_VAL(CONVECTION,FILE_NAME,'CONVECTION',line,ierr)
     CALL GET_LOGICAL_VAL(ADV_HLLC,FILE_NAME,'HLLC',line)
     IF(ierr==1)THEN
       if(myid.eq.0) WRITE(3,'(A12,A50)')'HIGH_ORDER', 'NOT DEFINED, USE DEFAULT'
       HIGH_ORDER='SECOND'
     ENDIF
     if(myid.eq.0) WRITE(3,'(A12,A50)')'HIGH_ORDER= ', HIGH_ORDER
     IF(ierr==1)THEN
       if(myid.eq.0) WRITE(3,'(A12,A50)')'TIME_ORDER', 'NOT DEFINED, USE DEFAULT'
       TIME_ORDER='THIRD'
     ENDIF
     if(myid.eq.0) WRITE(3,'(A12,A50)')'TIME_ORDER= ', TIME_ORDER

     ! ramp up the simulation
     CALL GET_Float_VAL(TRamp,FILE_NAME,'TRAMP',line)
     if(myid.eq.0) WRITE(3,'(A7,E12.3)')'TRAMP= ',TRamp

     ! if non-hydrostatic simulation
     CALL GET_LOGICAL_VAL(NON_HYDRO,FILE_NAME,'NON_HYDRO',line)
     if(myid.eq.0) WRITE(3,'(A11,L4)')'NON_HYDRO= ',NON_HYDRO

     ! poisson solver
     CALL GET_INTEGER_VAL(isolver,FILE_NAME,'ISOLVER',line)
     CALL GET_INTEGER_VAL(itmax,FILE_NAME,'ITMAX',line)
     CALL GET_Float_VAL(tol,FILE_NAME,'TOL',line)
     if(myid.eq.0) WRITE(3,'(A9,I2)')'ISOLVER= ',isolver
     if(myid.eq.0) WRITE(3,'(A7,I5)')'ITMAX= ',itmax
     if(myid.eq.0) WRITE(3,'(A5,E12.3)')'TOL= ',tol

     ! periodic bc
     CALL GET_LOGICAL_VAL(PERIODIC_X,FILE_NAME,'PERIODIC_X',line)
     CALL GET_LOGICAL_VAL(PERIODIC_Y,FILE_NAME,'PERIODIC_Y',line)
     if(myid.eq.0) WRITE(3,'(A12,L4)')'PERIODIC_X= ',PERIODIC_X
     if(myid.eq.0) WRITE(3,'(A12,L4)')'PERIODIC_Y= ',PERIODIC_Y

     ! boundary type
     CALL GET_INTEGER_VAL(Bc_X0,FILE_NAME,'BC_X0',line)
     CALL GET_INTEGER_VAL(Bc_Xn,FILE_NAME,'BC_Xn',line)
     CALL GET_INTEGER_VAL(Bc_Y0,FILE_NAME,'BC_Y0',line)
     CALL GET_INTEGER_VAL(Bc_Yn,FILE_NAME,'BC_Yn',line)
     CALL GET_INTEGER_VAL(Bc_Z0,FILE_NAME,'BC_Z0',line)
     CALL GET_INTEGER_VAL(Bc_Zn,FILE_NAME,'BC_Zn',line)
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_X0= ',Bc_X0
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Xn= ',Bc_Xn
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Y0= ',Bc_Y0
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Yn= ',Bc_Yn
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Z0= ',Bc_Z0
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Zn= ',Bc_Zn

     ! wavemaker 
     CALL GET_STRING_VAL(WaveMaker,FILE_NAME,'WAVEMAKER',line,ierr)
     if(myid.eq.0) WRITE(3,'(A11,A50)')'WAVEMAKER= ', WAVEMAKER
     IF(WaveMaker(1:3)=='LEF'.or.WaveMaker(1:3)=='RIG'  &
         .or.WaveMaker(1:3)=='INT'.or.WaveMaker(1:3)=='FLU')THEN
       CALL GET_Float_VAL(Amp_Wave,FILE_NAME,'AMP',line)
       CALL GET_Float_VAL(Per_Wave,FILE_NAME,'PER',line)
       CALL GET_Float_VAL(Dep_Wave,FILE_NAME,'DEP',line)
       CALL GET_Float_VAL(Theta_Wave,FILE_NAME,'THETA',line)
       if(myid.eq.0) WRITE(3,'(A9,F6.3)')'AMP_WAVE= ', Amp_Wave
       if(myid.eq.0) WRITE(3,'(A9,F6.3)')'PER_WAVE= ', Per_Wave
       if(myid.eq.0) WRITE(3,'(A9,F6.3)')'DEP_WAVE= ', Dep_Wave
       if(myid.eq.0) WRITE(3,'(A12,F6.3)')'THETA_WAVE= ', Theta_Wave

       IF(WaveMaker(1:3)=='INT') then
         CALL GET_Float_VAL(Xsource_West,FILE_NAME,'Xsource_West',line)
         CALL GET_Float_VAL(Xsource_East,FILE_NAME,'Xsource_East',line)
         CALL GET_Float_VAL(Ysource_Suth,FILE_NAME,'Ysource_Suth',line)
         CALL GET_Float_VAL(Ysource_Nrth,FILE_NAME,'Ysource_Nrth',line)
         if(myid.eq.0) WRITE(3,'(A14,F6.3)')'Xsource_West= ',Xsource_West
         if(myid.eq.0) WRITE(3,'(A14,F6.3)')'Xsource_East= ',Xsource_East
         if(myid.eq.0) WRITE(3,'(A14,F6.3)')'Ysource_Suth= ',Ysource_Suth
         if(myid.eq.0) WRITE(3,'(A14,F6.3)')'Ysource_Nrth= ',Ysource_Nrth
       ENDIF

       ! test periodicity
       IF(PERIODIC_Y.and.Theta_Wave.ne.Zero) then
         ! find wave number
         Segma = 2.0*pi/Per_Wave
         Celerity = sqrt(Grav*Dep_Wave)
         Wave_Length = Celerity*Per_Wave
         Wave_Number = 2.0*pi/Wave_Length

         Iter = 0
 75      Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
         if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 85
         Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &                
            Grav*tanh(Wave_Number*Dep_Wave)
         Wave_Number = Wave_Number-Fk/Fkdif
         Iter = Iter+1
         goto 75
 85      continue

         if(Theta_Wave>Zero) then
           ! find right angle for periodic bc
           tmp = Large       
           do Iter = 1,10000
             Wnumy = Iter*2.0*pi/(Nglob*dy)
             if(WnumY<Wave_Number) then
               ! theta based on Ky = K*sin(theta)
               tmp1 = asin(Wnumy/Wave_Number)*180./pi
               if(abs(tmp1-Theta_Wave)<tmp) then
                 tmp = abs(tmp1-Theta_Wave)
                 Theta_Calc = tmp1
               endif
             endif
           enddo
         elseif(Theta_Wave<Zero) then
           ! find right angle for periodic bc 
           tmp = Large
           do Iter = 1,10000
             Wnumy = Iter*2.0*pi/(Nglob*dy)
             if(WnumY<Wave_Number) then
               ! theta based on Ky = K*sin(theta)
               tmp1 = -asin(Wnumy/Wave_Number)*180./pi
               if(abs(tmp1-Theta_Wave)<tmp) then
                 tmp = abs(tmp1-Theta_Wave)
                 Theta_Calc = tmp1
               endif
             endif
           enddo
         endif

         if(myid.eq.0) then
           write(3,'(A20,F6.3)') 'Wave angle you set= ',Theta_Wave         
           write(3,'(A28,F6.3)') 'Wave angle for periodic bc= ',Theta_Calc
         endif
         Theta_Wave = Theta_Calc
       ENDIF
     ENDIF

     ! random wave, read in 2d spectrum
     IF(WaveMaker(5:7)=='SPC') then
       open(14,file='spc2d.txt')
       read(14,*) NumFreq,NumDir
       if(NumFreq>MaxNumFreq) then
         if(myid.eq.0) then
           write(3,'(A)') 'Please set a larger MaxNumFreq in mod_glob.F'
           stop
         endif
       endif
       if(NumDir>MaxNumDir) then
         if(myid.eq.0) then
           write(3,'(A)') 'Please set a larger MaxNumDir in mod_glob.F'
           stop
         endif
       endif
       do i = 1,NumFreq
         read(14,*) Freq(i)
       enddo
       do i = 1,NumDir
         read(14,*) Dire(i)
       enddo
       do j = 1,NumFreq
       do i = 1,NumDir
         read(14,*) Wave_Spc2d(i,j)
       enddo
       enddo
       close(14)

       ! random phase for each component
       do j = 1,NumFreq
       do i = 1,NumDir
         Random_Phs(i,j) = rand(0)*2.0*pi
       enddo
       enddo
     ENDIF

     ! JONSWAP spectrum
     if((WaveMaker(5:7)=='JON').or.(WaveMaker(5:7)=='TMA')) then
       CALL GET_Float_VAL(Hm0,FILE_NAME,'Hm0',line)
       CALL GET_Float_VAL(Tp,FILE_NAME,'Tp',line)
       CALL GET_Float_VAL(Freq_Min,FILE_NAME,'Freq_Min',line)
       CALL GET_Float_VAL(Freq_Max,FILE_NAME,'Freq_Max',line)
       CALL GET_INTEGER_VAL(NumFreq,FILE_NAME,'NumFreq',line) 
       if(myid.eq.0) WRITE(3,'(A5,f6.2)')'Hm0= ', Hm0
       if(myid.eq.0) WRITE(3,'(A4,f6.2)')'Tp= ', Tp
       if(myid.eq.0) WRITE(3,'(A10,f6.2)')'Freq_Min= ', Freq_Min
       if(myid.eq.0) WRITE(3,'(A10,f6.2)')'Freq_Max= ', Freq_Max
       if(myid.eq.0) WRITE(3,'(A9,I5)')'NumFreq= ', NumFreq

       ! jonswap spectrum
       gam = 3.3; sa = 0.07; sb = 0.09
       Freq_Peak = 1.0/Tp

       DFreq = (Freq_Max-Freq_Min)/NumFreq
       do i = 1,NumFreq
         Freq(i) = Freq_Min+0.5*DFreq+(i-1)*DFreq

         Per_Wave = 1.0/Freq(i)
         Segma = 2.0*pi/Per_Wave
         Celerity = sqrt(Grav*Dep_Wave)
         Wave_Length = Celerity*Per_Wave
         Wave_Number = 2.0*pi/Wave_Length
       
         Iter = 0
 76      Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
         if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 86
         Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  & 
                 Grav*tanh(Wave_Number*Dep_Wave)
         Wave_Number = Wave_Number-Fk/Fkdif
         Iter = Iter+1
         goto 76
 86      continue

         if(Freq(i)<Freq_Peak) then
           Jon_Spc(i) = Grav**2/Freq(i)**5*exp(-1.25*(Freq_Peak/Freq(i))**4)*  &
               gam**exp(-0.5*(Freq(i)/Freq_Peak-1.0)**2/sa**2)
         else
           Jon_Spc(i) = Grav**2/Freq(i)**5*exp(-1.25*(Freq_Peak/Freq(i))**4)*  &
               gam**exp(-0.5*(Freq(i)/Freq_Peak-1.0)**2/sb**2)
         endif

         if(WaveMaker(5:7)=='TMA') then
            Jon_Spc(i) = Jon_Spc(i)*tanh(Wave_Number*Dep_Wave)**2/  &
                 (1.0+2*Wave_Number*Dep_Wave/sinh(2.*Wave_Number*Dep_Wave))
         endif
       enddo
         
       ! make sure m0=Hm0**2/16=int S(f)df
       SumInt = Zero
       do i = 1,NumFreq
         SumInt = SumInt+Jon_Spc(i)*DFreq
       enddo
       A_Jon = Hm0**2/16.0/SumInt

       do i = 1,NumFreq
         Jon_Spc(i) = Jon_Spc(i)*A_Jon
         RanPhs(i) = rand()*2.0*pi
       enddo
     endif

     ! sponge layer
     CALL GET_LOGICAL_VAL(SPONGE_ON,FILE_NAME,'SPONGE_ON',line)
     if(myid.eq.0) WRITE(3,'(A11,L4)')'SPONGE_ON= ', SPONGE_ON
     IF(SPONGE_ON)THEN
       CALL GET_Float_VAL(Sponge_West_Width,FILE_NAME,'Sponge_West_Width',line)
       CALL GET_Float_VAL(Sponge_East_Width,FILE_NAME,'Sponge_East_Width',line)
       CALL GET_Float_VAL(Sponge_South_Width,FILE_NAME,'Sponge_South_Width',line)
       CALL GET_Float_VAL(Sponge_North_Width,FILE_NAME,'Sponge_North_Width',line)
       if(myid.eq.0) WRITE(3,'(A19,F6.3)')'Sponge_West_Width= ', Sponge_West_Width
       if(myid.eq.0) WRITE(3,'(A19,F6.3)')'Sponge_East_Width= ', Sponge_East_Width
       if(myid.eq.0) WRITE(3,'(A20,F6.3)')'Sponge_South_Width= ', Sponge_South_Width
       if(myid.eq.0) WRITE(3,'(A20,F6.3)')'Sponge_North_Width= ', Sponge_North_Width
     ENDIF

     ! wave average control
     CALL GET_LOGICAL_VAL(WAVE_AVERAGE_ON,FILE_NAME,'WAVE_AVERAGE_ON',line)
     CALL GET_Float_VAL(Wave_Ave_Start,FILE_NAME,'WAVE_AVERAGE_START',line)
     CALL GET_Float_VAL(Wave_Ave_End,FILE_NAME,'WAVE_AVERAGE_END',line)
     CALL GET_INTEGER_VAL(WaveheightID,FILE_NAME,'WaveheightID',line)





     
     ! whether to consider rheology
     CALL GET_LOGICAL_VAL(RHEOLOGY_ON,FILE_NAME,'RHEOLOGY_ON',line)
     CALL GET_Float_VAL(Yield_Stress,FILE_NAME,'Yield_Stress',line)
     CALL GET_Float_VAL(Plastic_Visc,FILE_NAME,'Plastic_Visc',line)

     ! if there is external forcing
     CALL GET_LOGICAL_VAL(EXTERNAL_FORCING,FILE_NAME,'EXTERNAL_FORCING',line)

     ! probe output
     CALL GET_INTEGER_VAL(NSTAT,FILE_NAME,'NSTAT',line)
     CALL GET_Float_VAL(Plot_Intv_Stat,FILE_NAME,'PLOT_INTV_STAT',line)
       if(myid.eq.0) WRITE(3,'(A7,I3)')'NSTAT= ', NSTAT
       if(myid.eq.0) WRITE(3,'(A19,F6.3)')'Plot_Intv_Stat= ', Plot_Intv_Stat

     if(NSTAT>0) then
       open(15,file='stat.txt',status='old')
       do i = 1,NSTAT
         read(15,*) xstat(i),ystat(i)
       enddo
       close(15)
     endif

     ! output parameters  
     CALL GET_LOGICAL_VAL(OUT_H,FILE_NAME,'OUT_H',line)
     CALL GET_LOGICAL_VAL(OUT_E,FILE_NAME,'OUT_E',line)
     CALL GET_LOGICAL_VAL(OUT_U,FILE_NAME,'OUT_U',line)
     CALL GET_LOGICAL_VAL(OUT_V,FILE_NAME,'OUT_V',line)
     CALL GET_LOGICAL_VAL(OUT_W,FILE_NAME,'OUT_W',line)
     CALL GET_LOGICAL_VAL(OUT_P,FILE_NAME,'OUT_P',line)
     CALL GET_LOGICAL_VAL(OUT_K,FILE_NAME,'OUT_K',line)
     CALL GET_LOGICAL_VAL(OUT_D,FILE_NAME,'OUT_D',line)
     CALL GET_LOGICAL_VAL(OUT_S,FILE_NAME,'OUT_S',line)
     CALL GET_LOGICAL_VAL(OUT_C,FILE_NAME,'OUT_C',line)
     CALL GET_LOGICAL_VAL(OUT_B,FILE_NAME,'OUT_B',line)
     CALL GET_LOGICAL_VAL(OUT_A,FILE_NAME,'OUT_A',line)
     CALL GET_LOGICAL_VAL(OUT_F,FILE_NAME,'OUT_F',line)
     CALL GET_LOGICAL_VAL(OUT_T,FILE_NAME,'OUT_T',line)
     CALL GET_LOGICAL_VAL(OUT_G,FILE_NAME,'OUT_G',line)
     CALL GET_LOGICAL_VAL(OUT_I,FILE_NAME,'OUT_I',line)
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_H= ',OUT_H
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_E= ',OUT_E
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_U= ',OUT_U
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_V= ',OUT_V
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_W= ',OUT_W
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_P= ',OUT_P
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_K= ',OUT_K
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_D= ',OUT_D
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_S= ',OUT_S
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_C= ',OUT_C
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_B= ',OUT_B
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_A= ',OUT_A
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_F= ',OUT_F
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_T= ',OUT_T
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_G= ',OUT_G
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_I= ',OUT_I

     if(myid.eq.0) WRITE(3,*)'--------------input end --------------'

     end subroutine read_input


     subroutine calculate_sponge
!-------------------------------------------------
!    Calculate sponge function
!    Called by
!      initial
!    Last update: 12/02/2011, Gangfeng Ma
!------------------------------------------------
     use global
     implicit none
     integer :: i,j

     if(Sponge_West_Width>Zero) then
       do j = 1,Nloc
       do i = 1,Mloc
         if(xc(i)<=Sponge_West_Width) then
           Sponge(i,j) = sqrt(1.0-  &
              min(((xc(i)-Sponge_West_Width)/Sponge_West_Width)**2,1.0))
         endif
       enddo
       enddo
     endif

     if(Sponge_East_Width>Zero)then
       do j = 1,Nloc
       do i = 1,Mloc
         if(xc(i)>=Mglob*dx-Sponge_East_Width) then
           Sponge(i,j) = sqrt(1.0-  &
             min(((xc(i)-(Mglob*dx-Sponge_East_Width))/Sponge_East_Width)**2,1.0))
         endif
       enddo
       enddo
     endif

     if(Sponge_South_Width>Zero)then
       do j = 1,Nloc
       do i = 1,Mloc
         if(yc(j)<=Sponge_South_Width) then
           Sponge(i,j) = sqrt(1.0-  &
              min(((yc(j)-Sponge_South_Width)/Sponge_South_Width)**2,1.0))
         endif
       enddo
       enddo
     endif

     if(Sponge_North_Width>Zero)then
       do j = 1,Nloc
       do i = 1,Mloc
         if(yc(j)>=Nglob*dy-Sponge_North_Width) then
           Sponge(i,j) = sqrt(1.0-  &
              min(((yc(j)-(Nglob*dy-Sponge_North_Width))/Sponge_North_Width)**2,1.0))
         endif
       enddo
       enddo
     endif

     end subroutine calculate_sponge


     subroutine sponge_damping
!---------------------------------------------------
!    This subroutine is used to damp waves using DHI type
!    sponge layer variables
!    Called by 
!      main
!    Last update: 12/02/2011, Gangfeng Ma
!--------------------------------------------------
     use global, only: Eta,Hc,D,U,V,W,Omega,Sponge,Mask, &
                       Mloc,Nloc,Kloc,DU,DV,DW
     implicit none
     integer :: i,j,k

     do j = 1,Nloc
     do i = 1,Mloc
       if(Mask(i,j)==1) then
         Eta(i,j) = Eta(i,j)*Sponge(i,j)
         D(i,j) = Eta(i,j)+Hc(i,j)

         ! W is estimated from continuity equation
         do k = 1,Kloc
           U(i,j,k) = U(i,j,k)*Sponge(i,j)
           V(i,j,k) = V(i,j,k)*Sponge(i,j)
           DU(i,j,k) = D(i,j)*U(i,j,k)
           DV(i,j,k) = D(i,j)*V(i,j,k)
         enddo
       endif
     enddo
     enddo
     
     end subroutine sponge_damping

 
     subroutine sigma_transform
!--------------------------------------------------- 
!    Calculate sigma transformation coefficient
!    Called by       
!      eval_duvw
!    Last update: 29/03/2011, Gangfeng Ma
!--------------------------------------------------
     use global, only: Zero,DelxSc,DelySc,D,DelxH,DelyH, &
                       DelxEta,DelyEta,sigc,Mloc,Nloc,Kloc
     implicit none
     integer :: i,j,k

     DelxSc = Zero
     DelySc = Zero
     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc
       DelxSc(i,j,k) = (1.0-sigc(k))/D(i,j)*DelxH(i,j)-sigc(k)/D(i,j)*DelxEta(i,j)
       DelySc(i,j,k) = (1.0-sigc(k))/D(i,j)*DelyH(i,j)-sigc(k)/D(i,j)*DelyEta(i,j)
     enddo
     enddo
     enddo

     end subroutine sigma_transform


     subroutine eval_turb(ISTEP)
!---------------------------------------------------
!    This subroutine is used to calculate viscosity
!    Called by                                                                                                             
!      main
!    Last update: 21/06/2011, Gangfeng Ma 
!--------------------------------------------------
     use global
     implicit none
     integer, intent(in) :: ISTEP
     integer :: i,j,k,n,m
     real(SP) :: DelsU,DelsV,Strxx,Stryy,Strxy,StrainMag,Smax
     real(SP), dimension(3,3) :: VelGrad,Stress
     real(SP), dimension(Mloc,Nloc,Kloc) :: StressMag

     ! laminar viscosity
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Cmu(i,j,k) = Visc
     enddo
     enddo
     enddo


     ! vertical turbulent viscosity
     if(IVturb==1) then
       ! constant vertical viscosity
       CmuVt = Cvs
     elseif(IVturb==2) then
       ! subgrid model
       do i = 1,Mloc
       do j = 1,Nloc
       do k = 2,Kloc-1
         DelsU = (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))
         DelsV = (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))
         CmuVt(i,j,k) = Cvs*D(i,j)*dsig(k)**2*sqrt(DelsU**2+DelsV**2)
       enddo
       CmuVt(i,j,1) = CmuVt(i,j,2)
       CmuVt(i,j,Kloc) = CmuVt(i,j,Kloc-1)
       enddo
       enddo
     elseif(IVturb==3) then
       ! k-epsilon turbulence model
       call kepsilon(ISTEP)
     elseif(IVturb==10) then
       ! 3D turbulence model
       call kepsilon_3D(ISTEP)
     elseif(IVturb==20) then
       call les_3D(ISTEP)
     endif
     
     ! horizontal turbulent viscosity
     if(IHturb==1) then
       ! constant viscosity
       CmuHt = Chs
     elseif(IHturb==2) then
       ! subgrid model
       do i = Ibeg,Iend
       do j = Jbeg,Jend
       do k = Kbeg,Kend
         Strxx = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+  &
               (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
         Stryy = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+  &
               (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
         Strxy = 0.5*((U(i,j+1,j)-U(i,j-1,k))/(2.0*dy)+  &
               (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)+  &
               (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+  &
               (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))
         CmuHt(i,j,k) = Chs*dx*dy*sqrt(Strxx**2+2.0*Strxy**2+Stryy**2)
       enddo
       enddo
       enddo

       ! ghost cell
       call phi_3D_exch(CmuHt)

       if(n_west.eq.MPI_PROC_NULL) then
       do j = Jbeg,Jend
       do k = Kbeg,Kend
         do i = 1,Nghost
           CmuHt(Ibeg-i,j,k) = CmuHt(Ibeg+i-1,j,k)
         enddo
       enddo
       enddo
       endif

       if(n_east.eq.MPI_PROC_NULL) then
       do j = Jbeg,Jend
       do k = Kbeg,Kend
         do i = 1,Nghost
           CmuHt(Iend+i,j,k) = CmuHt(Iend-i+1,j,k)
         enddo
       enddo
       enddo
       endif

       if(n_suth.eq.MPI_PROC_NULL) then
       do i = Ibeg,Iend
       do k = Kbeg,Kend
         do j = 1,Nghost
           CmuHt(i,Jbeg-j,k) = CmuHt(i,Jbeg+j-1,k)
         enddo
       enddo
       enddo
       endif

       if(n_nrth.eq.MPI_PROC_NULL) then
       do i = Ibeg,Iend
       do k = Kbeg,Kend
         do j = 1,Nghost
           CmuHt(i,Jend+j,k) = CmuHt(i,Jend-j+1,k)
         enddo
       enddo
       enddo
       endif

       do i = Ibeg,Iend
       do j = Jbeg,Jend
         do k = 1,Nghost
           CmuHt(i,j,Kbeg-k) = CmuHt(i,j,Kbeg+k-1)
         enddo
         do k = 1,Nghost
           CmuHt(i,j,Kend+k) = CmuHt(i,j,Kend-k+1)
         enddo
       enddo
       enddo  
     elseif(IHturb>=10) then
       ! use 3D turbulence model
       ! in this case, the length scales in all directions are
       ! in the same order
       CmuHt = CmuVt
     endif

     end subroutine eval_turb


     subroutine diffusion
!---------------------------------------------------------  
!    This subroutine is used to evaluate diffusion terms
!    in simplified form following FVCOM
!    Called by 
!      eval_duvw 
!    Last update: 20/12/2011, Gangfeng Ma 
!--------------------------------------------------------
     use global
     implicit none
     integer :: i,j,k
   
     Diffxx = zero; Diffxy = zero
     Diffyx = zero; Diffyy = zero
     Diffzx = zero; Diffzy = zero
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask(i,j)==0) cycle

       Diffxx(i,j,k) = ((0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+  &
               0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/Schmidt)* &
               (D(i+1,j)+D(i,j))*(U(i+1,j,k)-U(i,j,k))-  &
               (0.5*(Cmu(i,j,k)+Cmu(i-1,j,k))+0.5*(CmuHt(i,j,k)+  &
               CmuHt(i-1,j,k))/Schmidt)*(D(i,j)+D(i-1,j))* &
               (U(i,j,k)-U(i-1,j,k)))/dx**2
       Diffxy(i,j,k) = 0.5*((0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+  &
               0.5*(CmuHt(i,j+1,k)+CmuHt(i,j,k))/Schmidt)*  &
               (D(i,j+1)+D(i,j))*(U(i,j+1,k)-U(i,j,k))- &
               (0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+0.5*(CmuHt(i,j,k)+  &
               CmuHt(i,j-1,k))/Schmidt)*(D(i,j)+D(i,j-1))*  &
               (U(i,j,k)-U(i,j-1,k)))/dy**2+  &
               ((Cmu(i,j+1,k)+CmuHt(i,j+1,k)/Schmidt)*D(i,j+1)*  &
               (V(i+1,j+1,k)-V(i-1,j+1,k))-(Cmu(i,j-1,k)+  &
               CmuHt(i,j-1,k)/Schmidt)*D(i,j-1)*  &
               (V(i+1,j-1,k)-V(i-1,j-1,k)))/(4.0*dx*dy)

       Diffyx(i,j,k) = 0.5*((0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+  &
               0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/Schmidt)*  &
               (D(i+1,j)+D(i,j))*(V(i+1,j,k)-V(i,j,k))-  & 
               (0.5*(Cmu(i,j,k)+Cmu(i-1,j,k))+0.5*(CmuHt(i,j,k)+  &
               CmuHt(i-1,j,k))/Schmidt)*(D(i,j)+D(i-1,j))*  &
               (V(i,j,k)-V(i-1,j,k)))/dx**2+  &
               ((Cmu(i+1,j,k)+CmuHt(i+1,j,k)/Schmidt)*D(i+1,j)*  &
               (U(i+1,j+1,k)-U(i+1,j-1,k))-(Cmu(i-1,j,k)+  &
               CmuHt(i-1,j,k)/Schmidt)*D(i-1,j)*  &
               (U(i-1,j+1,k)-U(i-1,j-1,k)))/(4.0*dx*dy)
       Diffyy(i,j,k) = ((0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+  &
               0.5*(CmuHt(i,j+1,k)+CmuHt(i,j,k))/Schmidt)*  &
               (D(i,j+1)+D(i,j))*(V(i,j+1,k)-V(i,j,k))-  & 
               (0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+0.5*(CmuHt(i,j,k)+  &
               CmuHt(i,j-1,k))/Schmidt)*(D(i,j)+D(i,j-1))*  &
               (V(i,j,k)-V(i,j-1,k)))/dy**2

       Diffzx(i,j,k) = 0.5*((0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+  &
               0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/Schmidt)*  &
               (D(i+1,j)+D(i,j))*(W(i+1,j,k)-W(i,j,k))-  &
               (0.5*(Cmu(i,j,k)+Cmu(i-1,j,k))+0.5*(CmuHt(i,j,k)+  &
               CmuHt(i-1,j,k))/Schmidt)*(D(i,j)+D(i-1,j))*  &
               (W(i,j,k)-W(i-1,j,k)))/dx**2+  &
               ((Cmu(i+1,j,k)+CmuHt(i+1,j,k)/Schmidt)*  &
               (U(i+1,j,k+1)-U(i+1,j,k-1))/(sigc(k+1)-sigc(k-1))-  &
               (Cmu(i-1,j,k)+CmuHt(i-1,j,k)/Schmidt)*  &
               (U(i-1,j,k+1)-U(i-1,j,k-1))/(sigc(k+1)-sigc(k-1)))/(2.0*dx)
       Diffzy(i,j,k) = 0.5*((0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+  &
               0.5*(CmuHt(i,j+1,k)+CmuHt(i,j,k))/Schmidt)*  &
               (D(i,j+1)+D(i,j))*(W(i,j+1,k)-W(i,j,k))-  & 
               (0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+0.5*(CmuHt(i,j,k)+  &
               CmuHt(i,j-1,k))/Schmidt)*(D(i,j)+D(i,j-1))*  &
               (W(i,j,k)-W(i,j-1,k)))/dy**2+  &
               ((Cmu(i,j+1,k)+CmuHt(i,j+1,k)/Schmidt)*  &
               (V(i,j+1,k+1)-V(i,j+1,k-1))/(sigc(k+1)-sigc(k-1))-  &
               (Cmu(i,j-1,k)+CmuHt(i,j-1,k)/Schmidt)*  & 
               (V(i,j-1,k+1)-V(i,j-1,k-1))/(sigc(k+1)-sigc(k-1)))/(2.0*dy)
     enddo
     enddo
     enddo

     end subroutine diffusion


    subroutine phi_2D_exch(PHI)
    USE GLOBAL
    IMPLICIT NONE
    REAL(SP),INTENT(INOUT) :: PHI(Mloc,Nloc)

    INTEGER,DIMENSION(MPI_STATUS_SIZE,4) :: status
    INTEGER,DIMENSION(4) :: req
    INTEGER :: i,j,nreq,len
    REAL(SP),DIMENSION(Mloc,Nghost) :: rNmsg, sNmsg,rSmsg,sSmsg
    REAL(SP),DIMENSION(Nloc,Nghost) :: rWmsg, sWmsg,rEmsg,sEmsg

! for east-west

    len = Nloc * Nghost

    nreq = 0
    if ( n_west .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rWmsg, len, MPI_SP, &
            n_west, 0, comm2d, req(nreq), ier )
       do j = 1, Nloc
       do i = 1, Nghost
          sWmsg(j,i) = PHI(Ibeg+i-1,j)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sWmsg, len, MPI_SP, &
            n_west, 1, comm2d, req(nreq), ier )
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rEmsg, len, MPI_SP, &
            n_east, 1, comm2d, req(nreq), ier )
       do j = 1, Nloc
       do i = 1, Nghost
          sEmsg(j,i) = PHI(Iend-i+1,j)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sEmsg, len, MPI_SP, &
            n_east, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_west .ne. MPI_PROC_NULL ) then
       do j = 1, Nloc
       do i = 1, Nghost
          PHI(Ibeg-i,j) = rWmsg(j,i)
       enddo
       enddo
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       do j = 1, Nloc
       do i = 1, Nghost
          PHI(Iend+i,j) = rEmsg(j,i)
       enddo
       enddo
    endif

! for nrth-suth

    len = Mloc * Nghost

    nreq = 0
    if ( n_suth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rSmsg, len, MPI_SP, &
            n_suth, 0, comm2d, req(nreq), ier )
       do i = 1, Mloc
       do j = 1, Nghost
          sSmsg(i,j) = PHI(i,Jbeg+j-1)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sSmsg, len, MPI_SP, &
            n_suth, 1, comm2d, req(nreq), ier )
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rNmsg, len, MPI_SP, &
            n_nrth, 1, comm2d, req(nreq), ier )
       do i = 1, Mloc
       do j = 1, Nghost
          sNmsg(i,j) = PHI(i,Jend-j+1)
       enddo
       enddo
       nreq = nreq + 1
       call MPI_ISEND( sNmsg, len, MPI_SP, &
            n_nrth, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_suth .ne. MPI_PROC_NULL ) then
       do i = 1, Mloc
       do j = 1, Nghost
          PHI(i,Jbeg-j) = rSmsg(i,j)
       enddo
       enddo
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       do i = 1, Mloc
       do j = 1, Nghost
          PHI(i,Jend+j) = rNmsg(i,j)
       enddo
       enddo
    endif

    return
    END SUBROUTINE phi_2D_exch


    SUBROUTINE phi_3D_exch(PHI)
    USE GLOBAL
    IMPLICIT NONE
    REAL(SP),INTENT(INOUT) :: PHI(Mloc,Nloc,Kloc)

    INTEGER,DIMENSION(MPI_STATUS_SIZE,4) :: status
    INTEGER,DIMENSION(4) :: req
    INTEGER :: i,j,k,ik,jk,nreq,len
    REAL(SP),DIMENSION(Mloc*Kloc,Nghost) :: rNmsg, sNmsg,rSmsg,sSmsg
    REAL(SP),DIMENSION(Nloc*Kloc,Nghost) :: rWmsg, sWmsg,rEmsg,sEmsg

! for east-west

    len = Nloc * Kloc * Nghost

    nreq = 0
    if ( n_west .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rWmsg, len, MPI_SP, &
            n_west, 0, comm2d, req(nreq), ier )
       do k = 1, Kloc
       do j = 1, Nloc
       do i = 1, Nghost
          jk = (k-1)*Nloc+j
          sWmsg(jk,i) = PHI(Ibeg+i-1,j,k)
       enddo
       enddo
       enddo
       nreq = nreq + 1
       call MPI_ISEND( sWmsg, len, MPI_SP, &
            n_west, 1, comm2d, req(nreq), ier )
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rEmsg, len, MPI_SP, &
            n_east, 1, comm2d, req(nreq), ier )
       do k = 1, Kloc
       do j = 1, Nloc
       do i = 1, Nghost
          jk = (k-1)*Nloc+j
          sEmsg(jk,i) = PHI(Iend-i+1,j,k)
       enddo
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sEmsg, len, MPI_SP, &
            n_east, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_west .ne. MPI_PROC_NULL ) then
       do k = 1, Kloc
       do j = 1, Nloc
       do i = 1, Nghost
          jk = (k-1)*Nloc+j
          PHI(Ibeg-i,j,k) = rWmsg(jk,i)
       enddo
       enddo
       enddo
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       do k = 1, Kloc
       do j = 1, Nloc
       do i = 1, Nghost
          jk = (k-1)*Nloc+j
          PHI(Iend+i,j,k) = rEmsg(jk,i)
       enddo
       enddo
       enddo
    endif

! for nrth-suth

    len = Mloc * Kloc * Nghost

    nreq = 0
    if ( n_suth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rSmsg, len, MPI_SP, &
            n_suth, 0, comm2d, req(nreq), ier )
       do k = 1, Kloc
       do i = 1, Mloc
       do j = 1, Nghost
          ik = (k-1)*Mloc+i
          sSmsg(ik,j) = PHI(i,Jbeg+j-1,k)
       enddo
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sSmsg, len, MPI_SP, &
            n_suth, 1, comm2d, req(nreq), ier )
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rNmsg, len, MPI_SP, &
            n_nrth, 1, comm2d, req(nreq), ier )
       do k = 1, Kloc
       do i = 1, Mloc
       do j = 1, Nghost
          ik = (k-1)*Mloc+i
          sNmsg(ik,j) = PHI(i,Jend-j+1,k)
       enddo
       enddo
       enddo
       nreq = nreq + 1
       call MPI_ISEND( sNmsg, len, MPI_SP, &
            n_nrth, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_suth .ne. MPI_PROC_NULL ) then
       do k = 1, Kloc
       do i = 1, Mloc
       do j = 1, Nghost
          ik = (k-1)*Mloc+i
          PHI(i,Jbeg-j,k) = rSmsg(ik,j)
       enddo
       enddo
       enddo
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       do k = 1, Kloc
       do i = 1, Mloc
       do j = 1, Nghost
          ik = (k-1)*Mloc+i
          PHI(i,Jend+j,k) = rNmsg(ik,j)
       enddo
       enddo
       enddo
    endif

    return
    END SUBROUTINE phi_3D_exch


    ! Jeff added this subroutine to pass mask 02/14/2011
    SUBROUTINE phi_int_exch(PHI)
    USE GLOBAL
    IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: PHI(Mloc,Nloc)

    INTEGER,DIMENSION(MPI_STATUS_SIZE,4) :: status
    INTEGER,DIMENSION(4) :: req
    INTEGER :: i,j,nreq,len
    INTEGER,DIMENSION(Mloc,Nghost) :: rNmsg, sNmsg,rSmsg,sSmsg
    INTEGER,DIMENSION(Nloc,Nghost) :: rWmsg, sWmsg,rEmsg,sEmsg

! for east-west

    len = Nloc * Nghost

    nreq = 0
    if ( n_west .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rWmsg, len, MPI_INTEGER, &
            n_west, 0, comm2d, req(nreq), ier )
       do j = 1, Nloc
       do i = 1, Nghost
          sWmsg(j,i) = PHI(Ibeg+i-1,j)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sWmsg, len, MPI_INTEGER, &
            n_west, 1, comm2d, req(nreq), ier )
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rEmsg, len, MPI_INTEGER, &
            n_east, 1, comm2d, req(nreq), ier )
       do j = 1, Nloc
       do i = 1, Nghost
          sEmsg(j,i) = PHI(Iend-i+1,j)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sEmsg, len, MPI_INTEGER, &
            n_east, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_west .ne. MPI_PROC_NULL ) then
       do j = 1, Nloc
       do i = 1, Nghost
          PHI(Ibeg-i,j) = rWmsg(j,i)
       enddo
       enddo
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       do j = 1, Nloc
       do i = 1, Nghost
          PHI(Iend+i,j) = rEmsg(j,i)
       enddo
       enddo
    endif

! for nrth-suth

    len = Mloc * Nghost

    nreq = 0
    if ( n_suth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rSmsg, len, MPI_INTEGER, &
            n_suth, 0, comm2d, req(nreq), ier )
       do i = 1, Mloc
       do j = 1, Nghost
          sSmsg(i,j) = PHI(i,Jbeg+j-1)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sSmsg, len, MPI_INTEGER, &
            n_suth, 1, comm2d, req(nreq), ier )
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rNmsg, len, MPI_INTEGER, &
            n_nrth, 1, comm2d, req(nreq), ier )
       do i = 1, Mloc
       do j = 1, Nghost
          sNmsg(i,j) = PHI(i,Jend-j+1)
       enddo
       enddo
       nreq = nreq + 1
       call MPI_ISEND( sNmsg, len, MPI_INTEGER, &
            n_nrth, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_suth .ne. MPI_PROC_NULL ) then
       do i = 1, Mloc
       do j = 1, Nghost
          PHI(i,Jbeg-j) = rSmsg(i,j)
       enddo
       enddo
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       do i = 1, Mloc
       do j = 1, Nghost
          PHI(i,Jend+j) = rNmsg(i,j)
       enddo
       enddo
    endif
    END SUBROUTINE phi_int_exch


    subroutine adv_scalar_hlpa(Flx,Fly,Flz,Phi,R5,IVAR)
!--------------------------------------------------------
!   Subroutine for scalar convection and horizontal diffusion  
!   IVAR: indication of different scalars 
!    = 1: turbulent kinetic energy k
!    = 2: dissipation rate epsilon
!    = 3: salinity 
!    = 4: temperature
!    = 5: bubble number density 
!    = 6: sediment concentration 
!   Last update: Gangfeng Ma, 04/04/2012
!-------------------------------------------------------  
    use global
    implicit none
    integer, intent(in) :: IVAR
    real(SP), dimension(Mloc,Nloc,Kloc),  intent(in) :: Phi
    real(SP), dimension(Mloc1,Nloc,Kloc),  intent(in) :: Flx
    real(SP), dimension(Mloc,Nloc1,Kloc),  intent(in) :: Fly
    real(SP), dimension(Mloc,Nloc,Kloc1),  intent(in) :: Flz
    real(SP), dimension(Mloc,Nloc,Kloc), intent(inout) :: R5
    real(SP), dimension(:,:,:), allocatable :: Scalx,Scaly,Scalz,Sdiffx,Sdiffy
    real(SP) :: DUfs,DVfs,Wfs,Fww,Fw,Fp,Fe,hlpa,SchtH
    integer :: i,j,k

    allocate(Scalx(Mloc1,Nloc,Kloc))
    allocate(Scaly(Mloc,Nloc1,Kloc))
    allocate(Scalz(Mloc,Nloc,Kloc1))
    allocate(Sdiffx(MLoc,Nloc,Kloc))
    allocate(Sdiffy(Mloc,Nloc,Kloc))

    ! advection in x direction
    Scalx = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend+1
      DUfs = Flx(i,j,k)
      Fww = Phi(i-2,j,k)
      Fw  = Phi(i-1,j,k)
      Fp  = Phi(i,j,k)
      Fe  = Phi(i+1,j,k)
      Scalx(i,j,k) = DUfs*hlpa(DUfs,Fww,Fw,Fp,Fe)
    enddo
    enddo
    enddo

    ! advection in y direction
    Scaly = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend+1
    do i = Ibeg,Iend      
      DVfs = Fly(i,j,k)
      Fww = Phi(i,j-2,k)
      Fw  = Phi(i,j-1,k)
      Fp  = Phi(i,j,k)
      Fe  = Phi(i,j+1,k)
      Scaly(i,j,k) = DVfs*hlpa(DVfs,Fww,Fw,Fp,Fe)
    enddo
    enddo
    enddo

    ! advection in z direction
    Scalz = Zero
    do k = Kbeg+1,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      Wfs = Flz(i,j,k)
      Fww = Phi(i,j,k-2)
      Fw  = Phi(i,j,k-1)
      Fp  = Phi(i,j,k)
      Fe  = Phi(i,j,k+1)
      Scalz(i,j,k) = Wfs*hlpa(Wfs,Fww,Fw,Fp,Fe)
    enddo
    enddo
    enddo

    ! at boundaries
    call flux_scalar_bc(IVAR,Scalx,Scaly,Scalz)

    ! Schmidt number
    if(IVAR==1) then  ! tke eq.
      SchtH = 1.0
    elseif(IVAR==2) then  ! epsilon eq.
      SchtH = 1.3
    elseif(IVAR==5) then ! bubble
      SchtH = 0.7
    elseif(IVAR==6) then ! sediment
      SchtH = 1.0
    else
      SchtH = 1.0
    endif

    ! estimate horizontal diffusion
    Sdiffx = Zero; Sdiffy = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==0) cycle
      Sdiffx(i,j,k) = 0.5*((0.5*(Cmu(i+1,j,k)+Cmu(i,j,k))+  &
           0.5*(CmuHt(i+1,j,k)+CmuHt(i,j,k))/SchtH)*  &
           (D(i+1,j)+D(i,j))*(Phi(i+1,j,k)-Phi(i,j,k))-  &
           (0.5*(Cmu(i,j,k)+Cmu(i-1,j,k))+  &
           0.5*(CmuHt(i,j,k)+CmuHt(i-1,j,k))/SchtH)*  &
           (D(i,j)+D(i-1,j))*(Phi(i,j,k)-Phi(i-1,j,k)))/dx**2
      Sdiffy(i,j,k) = 0.5*((0.5*(Cmu(i,j+1,k)+Cmu(i,j,k))+  &
           0.5*(CmuHt(i,j+1,k)+CmuHt(i,j,k))/SchtH)*  &
           (D(i,j+1)+D(i,j))*(Phi(i,j+1,k)-Phi(i,j,k))-  &
           (0.5*(Cmu(i,j,k)+Cmu(i,j-1,k))+  &
           0.5*(CmuHt(i,j,k)+CmuHt(i,j-1,k))/SchtH)*  &
           (D(i,j)+D(i,j-1))*(Phi(i,j,k)-Phi(i,j-1,k)))/dy**2
    enddo
    enddo
    enddo

    R5 = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==0) cycle
      R5(i,j,k) = -1.0/dx*(Scalx(i+1,j,k)-Scalx(i,j,k))-  &
                  1.0/dy*(Scaly(i,j+1,k)-Scaly(i,j,k))  &
                  -1.0/dsig(k)*(Scalz(i,j,k+1)-Scalz(i,j,k))+ &
                  Sdiffx(i,j,k)+Sdiffy(i,j,k)
    enddo
    enddo
    enddo

    deallocate(Scalx)
    deallocate(Scaly)
    deallocate(Scalz)
    deallocate(Sdiffx)
    deallocate(Sdiffy)

    end subroutine adv_scalar_hlpa


    function hlpa(Uw,Fww,Fw,Fp,Fe)
!-------------------------------------------------------
!   HLPA scheme
!-------------------------------------------------------
    use global, only: SP,Zero
    implicit none
    real(SP), intent(in)  :: Uw,Fww,Fw,Fp,Fe
    real(SP) :: hlpa,Alpha_pl,Alpha_mn

    if(Uw>=Zero) then
      if(abs(Fp-2.*Fw+Fww)<abs(Fp-Fww)) then
        Alpha_pl = 1.0
      else
        Alpha_pl = 0.0
      endif

      if(abs(Fp-Fww)<=1.e-16) then
        hlpa = Fw
      else
        hlpa = Fw+Alpha_pl*(Fp-Fw)*(Fw-Fww)/(Fp-Fww)
      endif
    endif

    if(Uw<Zero) then
      if(abs(Fw-2.*Fp+Fe)<abs(Fw-Fe)) then
        Alpha_mn = 1.0
      else
        Alpha_mn = 0.0
      endif

      if(abs(Fw-Fe)<=1.e-16) then
        hlpa = Fp
      else
        hlpa = Fp+Alpha_mn*(Fw-Fp)*(Fp-Fe)/(Fw-Fe)
      endif
    endif

    return
    end function hlpa


    subroutine flux_scalar_bc(IVAR,Scalx,Scaly,Scalz)
!--------------------------------------------------------
!   Specify boundary conditions for scalar convection
!   Last update: Gangfeng Ma, 09/02/2011
!-------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: IVAR
    real(SP), dimension(Mloc1,Nloc,Kloc), intent(inout) :: Scalx
    real(SP), dimension(Mloc,Nloc1,Kloc), intent(inout) :: Scaly
    real(SP), dimension(Mloc,Nloc,Kloc1), intent(inout) :: Scalz
    real(SP), dimension(Nloc,Kloc) :: Scal_X0,Scal_Xn
    integer :: i,j,k

    ! temporarily set it here

    ! left and right side
     if(n_west.eq.MPI_PROC_NULL) then
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_X0==1.or.Bc_X0==2) then
        Scalx(Ibeg,j,k) = Zero
       elseif(Bc_X0==3) then
         Scalx(Ibeg,j,k) = Ex(Ibeg,j,k)*Scal_X0(j,k)
       endif
     enddo
     enddo
     endif

     if(n_east.eq.MPI_PROC_NULL) then
     do j = Jbeg,Jend
     do k = Kbeg,Kend
!       if(Bc_Xn==1.or.Bc_Xn==2) then
         Scalx(Iend1,j,k) = Zero
!       elseif(Bc_Xn==3) then
!         Scalx(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)*Scal_Xn(j,k)
!       endif
     enddo
     enddo
     endif

     ! front and back side
     if(n_suth.eq.MPI_PROC_NULL) then
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Y0==1.or.Bc_Y0==2) then
         Scaly(i,Jbeg,k) = Zero
       endif
     enddo
     enddo
     endif


     if(n_nrth.eq.MPI_PROC_NULL) then
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Yn==1.or.Bc_Yn==2) then
         Scaly(i,Jend1,k) = Zero
       endif
     enddo
     enddo
     endif

    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==0) then
        Scalx(i,j,k) = Zero
        Scalx(i+1,j,k) = Zero
        Scaly(i,j,k) = Zero
        Scaly(i,j+1,k) = Zero
      endif
    enddo
    enddo
    enddo

    do j = Jbeg,Jend
    do i = Ibeg,Iend
      Scalz(i,j,Kbeg) = Zero
      Scalz(i,j,Kend1) = Zero
    enddo
    enddo

    return
    end subroutine flux_scalar_bc


    subroutine les_3D(ISTEP)
!------------------------------------------------------
!   large eddy simulation (LES)
!   Last update: Gangfeng Ma, 09/22/2011
!------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: ISTEP
    real(SP), parameter :: Dmin = 0.02
    integer :: i,j,k,Iter
    real(SP) :: S11,S22,S33,S12,S13,S23,SijSij,Filter
    real(SP) :: Umag,Zdis,X0,Xa,Xn,FricU
 
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(Mask9(i,j)==1.and.D(i,j)>Dmin) then
        S11 = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)                                      
        S22 = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)                                      
        S33 = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))
        S12 = 0.5*((U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)+  &                             
              (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)) 
        S13 = 0.5*((U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                                                                     
              (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))
        S23 = 0.5*((V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                                                                     
              (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k))
        SijSij = S11**2+S22**2+S33**2+2.0*(S12**2+S13**2+S23**2)
        Filter = (dx*dy*dsig(k)*D(i,j))**(1./3.)
        CmuVt(i,j,k) = (Cvs*Filter)**2*sqrt(2.0*SijSij)
      endif
    enddo
    enddo
    enddo

    ! ghost cells
    ! at the bottom
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      ! impose wall function
      Umag = sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)
      if(Umag<1.e-6) then
        CmuVt(i,j,Kbeg) = Cmu(i,j,Kbeg)
      else
        Zdis = 0.5*dsig(Kbeg)*D(i,j)
        if(Zob>=0.1) then
          ! rough wall
          FricU = Umag/(1.0/0.41*log(30.0*Zdis/Zob))
        else
          ! smooth wall
          X0 = 0.05
          Iter = 0
       
          Xa = dlog(9.0*Umag*Zdis/Cmu(i,j,Kbeg))
 10       Xn = X0+(0.41-X0*(Xa+dlog(X0)))/(1.0+0.41/X0)
          if(Iter>=20) then
            write(*,*) 'Iteration exceeds 20 steps',i,j,Umag
          endif
          if(dabs((Xn-X0)/X0)>1.e-8.and.Xn>0.0) then
            X0 = Xn
            Iter = Iter+1
            goto 10
          else
            FricU = Xn*Umag
          endif
        endif

        CmuVt(i,j,Kbeg) = 0.41*Zdis*FricU*  &
           (1.0-exp(-Zdis*FricU/Cmu(i,j,Kbeg)/19.0))**2
      endif
 100  continue

      do k = 1,Nghost
        CmuVt(i,j,Kbeg-k) = CmuVt(i,j,Kbeg+k-1)
      enddo
    enddo
    enddo

    ! at the free surface
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        CmuVt(i,j,Kend+k) = CmuVt(i,j,Kend-k+1)
      enddo
    enddo
    enddo

    call phi_3D_exch(CmuVt)

     if(n_west.eq.MPI_PROC_NULL) then
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(WaveMaker(1:3)=='LEF') then
         ! no turbulence at wave generation region
         CmuVt(Ibeg,j,k) = Zero
         do i = 1,Nghost
           CmuVt(Ibeg-i,j,k) = CmuVt(Ibeg+i-1,j,k)
         enddo
       else       
         do i = 1,Nghost
           CmuVt(Ibeg-i,j,k) = CmuVt(Ibeg+i-1,j,k)
         enddo
       endif
     enddo
     enddo
     endif

     if(n_east.eq.MPI_PROC_NULL) then
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       do i = 1,Nghost
         CmuVt(Iend+i,j,k) = CmuVt(Iend-i+1,j,k)
       enddo
     enddo
     enddo
     endif
    
     if(n_suth.eq.MPI_PROC_NULL) then
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         CmuVt(i,Jbeg-j,k) = CmuVt(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
     endif

     if(n_nrth.eq.MPI_PROC_NULL) then
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         CmuVt(i,Jend+j,k) = CmuVt(i,Jend-j+1,k)
       enddo
     enddo
     enddo
     endif

     ! no turbulence in the internal wavemaker region   
    if(WaveMaker(1:3)=='INT') then
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
          CmuVt(i,j,k) = Zero
        endif
      enddo
      enddo
      enddo
    endif

    ! estimate turbulent dissipation rate (Van den Hengel et al., 2005)
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      Filter = (dx*dy*dsig(k)*D(i,j))**(1./3.)
      Eps(i,j,k) = 2.0*CmuVt(i,j,k)**3/(Cvs*Filter)**4
    enddo
    enddo
    enddo

    end subroutine les_3D


    subroutine kepsilon_3D(ISTEP)
!-------------------------------------------------------
!   k-epsilon turbulence model
!   Last update: Gangfeng Ma, 09/07/2011
!-------------------------------------------------------
    use global
    implicit none
    integer,  intent(in) :: ISTEP
    integer,  parameter :: ke_model = 2
    real(SP), parameter :: Dmin = 0.02
    real(SP), dimension(:,:,:), allocatable :: R5,DelzR,Tke_Old,Eps_Old,DUfs,DVfs,Wfs
    real(SP), dimension(:,:), allocatable :: VelGrad,ReynoldStress,Vorticity
    real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0
    real(SP) :: c1e,c2e,c3e,cmiu,cfk,cfe,Umag,Zdis,X0,Xa,Xn,FricU,Sche,Schk
    real(SP) :: smax,dmax,c_d,c_1,c_2,c_3,delta_nm,Tkeb,Epsb
    real(SP) :: S11,S22,S33,S12,S13,S23
    integer :: i,j,k,n,m,l,g,IVAR,Iter,Nlen
    
    allocate(R5(Mloc,Nloc,Kloc))
    allocate(DelzR(Mloc,Nloc,Kloc))
    allocate(Tke_Old(Mloc,Nloc,Kloc))
    allocate(Eps_Old(Mloc,Nloc,Kloc))
    allocate(DUfs(Mloc1,Nloc,Kloc))
    allocate(DVfs(Mloc,Nloc1,Kloc))
    allocate(Wfs(Mloc,Nloc,Kloc1))
    allocate(VelGrad(3,3))
    allocate(ReynoldStress(3,3))
    allocate(Vorticity(3,3))

    ! some parameters
    c1e = 1.44
    c2e = 1.92
!    c3e = -1.4
    c3e = 0.0
    cmiu = 0.09
    cfk = 1.0
    cfe = 1.33
    Sche = 1.3
    Schk = 1.0

    ! save old values
    Tke_Old = Tke
    Eps_Old = Eps

    Prod_s = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(ke_model==1) then
        ! linear model
        if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
          S11 = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          S22 = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          S33 = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1)) 
          S12 = 0.5*((U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)+  &
              (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))                   
          S13 = 0.5*((U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                                            
              (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))                   
          S23 = 0.5*((V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                                            
              (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k))
          Prod_s(i,j,k) = 2.0*CmuVt(i,j,k)*(S11**2+S22**2+S33**2+2.0*(S12**2+S13**2+S23**2))
        endif
      elseif(ke_model==2) then
        ! nonlinear model (Lin and Liu, 1998)
        ! Notice: if using nonlinear model, the initial seeding of tke and epsilon
        !         cannot be zero in order to generate turbulence production. Check 
        !         tke_min and eps_min in subroutine initial. 
        if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
          ! estimate gradient first
          VelGrad = Zero

          VelGrad(1,1) = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+  &
                  (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(1,2) = (U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+  &
                  (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          VelGrad(1,3) = 1./D(i,j)*(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))
          VelGrad(2,1) = (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+  &
                  (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(2,2) = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+  &
                  (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k) 
          VelGrad(2,3) = 1./D(i,j)*(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))
          VelGrad(3,1) = (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+  &
                  (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(3,2) = (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+  &
                  (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          VelGrad(3,3) = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))

          ! estimate Reynolds stress
          ReynoldStress = Zero
          
          smax = zero
          do n = 1,3
            if(abs(VelGrad(n,n))>smax) smax = abs(VelGrad(n,n))
          enddo
          smax = smax*Tke_Old(i,j,k)/Eps_Old(i,j,k)
          c_d = (1./3.)*(1./(3.7+smax))

          dmax = zero
          do n = 1,3
          do m = 1,3
            if(abs(VelGrad(n,m))>dmax) dmax = abs(VelGrad(n,m))
          enddo
          enddo
          dmax = dmax*Tke_Old(i,j,k)/Eps_Old(i,j,k)
          c_1 = 2./3./(123.5+2.0*dmax**2)
          c_2 = -2./3./(39.2+2.0*dmax**2)
          c_3 = 2./3./(246.9+2.0*dmax**2)

          do n = 1,3
          do m = 1,3
            if(n==m) then
              delta_nm = 1.
            else
              delta_nm = 0.
            endif

            ReynoldStress(n,m) = c_d*Tke_Old(i,j,k)**2/Eps_Old(i,j,k)*(VelGrad(n,m)+VelGrad(m,n))-  &
                   (2./3.)*Tke_Old(i,j,k)*delta_nm

            do l = 1,3
              ReynoldStress(n,m) = ReynoldStress(n,m)+  &
                   c_1*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*(VelGrad(n,l)*VelGrad(l,m)+  &
                   VelGrad(m,l)*VelGrad(l,n))
              do g = 1,3
                ReynoldStress(n,m) = ReynoldStress(n,m)-  &
                   (2./3.)*c_1*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*VelGrad(l,g)*VelGrad(g,l)*delta_nm
              enddo
            enddo

            do g = 1,3
              ReynoldStress(n,m) = ReynoldStress(n,m)+  &
                   c_2*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*VelGrad(n,g)*VelGrad(m,g)
              do l = 1,3
                ReynoldStress(n,m) = ReynoldStress(n,m)-  &
                   (1./3.)*c_2*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*VelGrad(l,g)*VelGrad(l,g)*delta_nm
              enddo
            enddo

            do g = 1,3
              ReynoldStress(n,m) = ReynoldStress(n,m)+  &
                   c_3*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*VelGrad(g,n)*VelGrad(g,m)
              do l = 1,3
                ReynoldStress(n,m) = ReynoldStress(n,m)-  &
                   (1./3.)*c_3*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*VelGrad(l,g)*VelGrad(l,g)*delta_nm 
              enddo
            enddo
          enddo
          enddo

          ! estimate shear production
          do n = 1,3
          do m = 1,3
            Prod_s(i,j,k) = Prod_s(i,j,k)+ReynoldStress(n,m)*VelGrad(n,m)
          enddo
          enddo

          !! no negative production at the surface
          !if(k==Kend.and.Prod_s(i,j,k)<0.0) Prod_s(i,j,k) = Zero

          ! Do not allow negative production
          if(Prod_s(i,j,k)<0.0) Prod_s(i,j,k) = Zero

        endif
      elseif(ke_model==3) then
        ! Following Mayer and Madsen (2000), instead of determining the production on the
        ! basis of the strain rate, the production is based on the rotation of the velocity field.
        if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
          ! estimate gradient first                              
          VelGrad = Zero

          VelGrad(1,1) = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+  &
                  (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(1,2) = (U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+  &
                  (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)                  
          VelGrad(1,3) = 1./D(i,j)*(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))                                                        
          VelGrad(2,1) = (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+  &
                  (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)                     
          VelGrad(2,2) = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+  &
                  (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)                     
          VelGrad(2,3) = 1./D(i,j)*(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))                                                        
  
          VelGrad(3,1) = (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+  &
                  (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)                     
          VelGrad(3,2) = (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+  &
                  (W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)                     
          VelGrad(3,3) = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))

          ! vorticity field and production
          do n = 1,3
          do m = 1,3
            Vorticity(n,m) = VelGrad(n,m)-VelGrad(m,n)
            Prod_s(i,j,k) = Prod_s(i,j,k)+CmuVt(i,j,k)*Vorticity(n,m)*Vorticity(n,m)
          enddo
          enddo
        endif
      endif
    enddo
    enddo
    enddo

    ! buoyancy production
    Prod_b = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DelzR(i,j,k) = (Rho(i,j,k+1)-Rho(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)
        Prod_b(i,j,k) = Grav*CmuVt(i,j,k)*DelzR(i,j,k)/Rho0
      endif
    enddo
    enddo
    enddo
 
    ! flux Richardson number
    Richf = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(Prod_s(i,j,k)>Zero) then
        Richf(i,j,k) = -Prod_b(i,j,k)/(Prod_s(i,j,k)+1.0e-16)
        Richf(i,j,k) = dmax1(0.0,dmin1(0.21,Richf(i,j,k)))
      endif
    enddo
    enddo
    enddo


    Nlen = Kend-Kbeg+1
    allocate(Acoef(Nlen))
    allocate(Bcoef(Nlen))
    allocate(Ccoef(Nlen))
    allocate(Xsol(Nlen))
    allocate(Rhs0(Nlen))

    ! transport velocities
    DUfs = Ex
    DVfs = Ey
    Wfs = Omega

    ! solve epsilon equation
    IVAR = 2
    call adv_scalar_hlpa(DUfs,DVfs,Wfs,Eps_Old,R5,IVAR)

    do i = Ibeg,Iend
    do j = Jbeg,Jend
      if(D(i,j)<Dmin.and.Mask(i,j)==0) cycle

      Nlen = 0
      do k = Kbeg,Kend
        R5(i,j,k) = R5(i,j,k)+c1e*D(i,j)*(Prod_s(i,j,k)+  &
                c3e*Prod_b(i,j,k))*Eps_Old(i,j,k)/Tke_Old(i,j,k)                        
        Nlen = Nlen+1
        if(k==Kbeg) then
          Acoef(Nlen) = 0.0
        else
          Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
              0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Sche)/  &
              (0.5*dsig(k)*(dsig(k)+dsig(k-1)))
        endif

        if(k==Kend) then
          Ccoef(Nlen) = 0.0
        else
          Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
              0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Sche)/  &
              (0.5*dsig(k)*(dsig(k)+dsig(k+1)))
        endif
        
        Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)+dt*c2e*Eps_Old(i,j,k)/Tke_Old(i,j,k)

        Rhs0(Nlen) = DEps(i,j,k)+dt*R5(i,j,k)
      enddo
      
      call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

      Nlen = 0
      do k = Kbeg,Kend
        Nlen = Nlen+1
        DEps(i,j,k) = Xsol(Nlen)
      enddo
    enddo
    enddo

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DEps(i,j,k) = ALPHA(ISTEP)*DEps0(i,j,k)+BETA(ISTEP)*DEps(i,j,k)
        DEps(i,j,k) = dmax1(DEps(i,j,k),D(i,j)*Eps_min)
      endif
    enddo
    enddo
    enddo

    ! slove tke equation
    IVAR = 1
    call adv_scalar_hlpa(DUfs,DVfs,Wfs,Tke_Old,R5,IVAR)

    do i = Ibeg,Iend
    do j = Jbeg,Jend
      if(D(i,j)<Dmin.and.Mask(i,j)==0) cycle

      Nlen = 0
      do k = Kbeg,Kend
        R5(i,j,k) = R5(i,j,k)+D(i,j)*(Prod_s(i,j,k)+Prod_b(i,j,k))-DEps(i,j,k)                                                   
        Nlen = Nlen+1
        if(k==Kbeg) then
          Acoef(Nlen) = 0.0
        else
          Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+  &
              0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schk)/  &
              (0.5*dsig(k)*(dsig(k)+dsig(k-1)))
        endif

        if(k==Kend) then
          Ccoef(Nlen) = 0.0
        else
          Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+  &
              0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schk)/  &
              (0.5*dsig(k)*(dsig(k)+dsig(k+1)))
        endif
        
        Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)
        Rhs0(Nlen) = DTke(i,j,k)+dt*R5(i,j,k)
      enddo
      
      call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

      Nlen = 0
      do k = Kbeg,Kend
        Nlen = Nlen+1
        DTke(i,j,k) = Xsol(Nlen)
      enddo
    enddo
    enddo

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DTke(i,j,k) = ALPHA(ISTEP)*DTke0(i,j,k)+BETA(ISTEP)*DTke(i,j,k)
        DTke(i,j,k) = dmax1(DTke(i,j,k),D(i,j)*Tke_min)
      endif
    enddo
    enddo
    enddo

    ! at the bottom
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      if(D(i,j)<Dmin.or.Mask(i,j)==0) cycle

      ! impose wall function 
      Umag = sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)
      if(Umag<1.e-6) then
        Tkeb = Tke_min
        Epsb = Eps_min
  
        DTke(i,j,Kbeg) = D(i,j)*Tkeb
        DEps(i,j,Kbeg) = D(i,j)*Epsb
      else
        Zdis = 0.5*dsig(Kbeg)*D(i,j)

        X0 = 0.05
        Iter = 0

        Xa = dlog(9.0*Umag*Zdis/Visc)
 10     Xn = X0+(0.41-X0*(Xa+dlog(X0)))/(1.0+0.41/X0)
        if(Iter>=20) then
          write(*,*) 'Iteration exceeds 20 steps',i,j,Umag
        endif
        if(dabs((Xn-X0)/X0)>1.e-8.and.Xn>0.0) then
          X0 = Xn
          Iter = Iter+1
          goto 10
        else
          FricU = Xn*Umag
        endif

!        if(Ibot==1) then
!          FricU = sqrt(Cd0)*Umag
!        else
!          FricU = Umag/(1./Kappa*log(30.*Zdis/Zob))
!        endif

        Tkeb = FricU**2/sqrt(cmiu)
        Epsb = FricU**3/(Kappa*Zdis)

        DTke(i,j,Kbeg) = D(i,j)*Tkeb
        DEps(i,j,Kbeg) = D(i,j)*Epsb
      endif

      do k = 1,Nghost
        DTke(i,j,Kbeg-k) = DTke(i,j,Kbeg+k-1)
        DEps(i,j,Kbeg-k) = DEps(i,j,Kbeg+k-1)
      enddo
    enddo
    enddo

    ! at the free surface
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        DTke(i,j,Kend+k) = DTke(i,j,Kend-k+1)
        DEps(i,j,Kend+k) = DEps(i,j,Kend-k+1)
      enddo
    enddo
    enddo

    call phi_3D_exch(DTke)
    call phi_3D_exch(DEps)

     if(n_west.eq.MPI_PROC_NULL) then
     do j = Jbeg,Jend
     do k = 1,Kloc
       do i = 1,Nghost
         DTke(Ibeg-i,j,k) = DTke(Ibeg+i-1,j,k)
         DEps(Ibeg-i,j,k) = DEps(Ibeg+i-1,j,k)
       enddo
     enddo
     enddo
     endif

     if(n_east.eq.MPI_PROC_NULL) then
     do j = Jbeg,Jend
     do k = 1,Kloc
       do i = 1,Nghost
         DTke(Iend+i,j,k) = DTke(Iend-i+1,j,k)
         DEps(Iend+i,j,k) = DEps(Iend-i+1,j,k)
       enddo
     enddo
     enddo
     endif
    
     if(n_suth.eq.MPI_PROC_NULL) then
     do i = 1,Mloc
     do k = 1,Kloc
       do j = 1,Nghost
         DTke(i,Jbeg-j,k) = DTke(i,Jbeg+j-1,k)
         DEps(i,Jbeg-j,k) = DEps(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
     endif

     if(n_nrth.eq.MPI_PROC_NULL) then
     do i = 1,Mloc
     do k = 1,Kloc
       do j = 1,Nghost
         DTke(i,Jend+j,k) = DTke(i,Jend-j+1,k)
         DEps(i,Jend+j,k) = DEps(i,Jend-j+1,k)
       enddo
     enddo
     enddo
     endif

   
    do i = 1,Mloc
    do j = 1,Nloc
    do k = 1,Kloc
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        Tke(i,j,k) = DTke(i,j,k)/D(i,j)
        Eps(i,j,k) = DEps(i,j,k)/D(i,j)
        if(Tke(i,j,k)<1.e-6.or.Eps(i,j,k)<1.e-6) then
          CmuVt(i,j,k) = Cmut_min
        else
          CmuVt(i,j,k) = Cmiu*Tke(i,j,k)**2/Eps(i,j,k)
        endif
      else
        Tke(i,j,k) = Tke_min
        Eps(i,j,k) = Eps_min
        DTke(i,j,k) = D(i,j)*Tke_min
        DEps(i,j,k) = D(i,j)*Eps_min
        CmuVt(i,j,k) = Cmut_min
      endif
    enddo
    enddo
    enddo

    ! no turbulence in the internal wavemaker region
    if(WaveMaker(1:3)=='INT') then
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
          Tke(i,j,k) = Tke_min
          Eps(i,j,k) = Eps_min
          DTke(i,j,k) = D(i,j)*Tke_min
          DEps(i,j,k) = D(i,j)*Eps_min
          CmuVt(i,j,k) = Cmut_min
        endif  
      enddo
      enddo
      enddo
    endif

    ! Reynolds stress (just for output) 
    UpWp = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==1) then
        VelGrad(1,2) = (U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+  &
                  (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
        VelGrad(2,1) = (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+  &
                  (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
        UpWp(i,j,k) = CmuVt(i,j,k)*(VelGrad(1,2)+VelGrad(2,1))
      endif
    enddo
    enddo
    enddo

    deallocate(R5)
    deallocate(DelzR)
    deallocate(Tke_Old)
    deallocate(Eps_Old)
    deallocate(VelGrad)
    deallocate(ReynoldStress)
    deallocate(Acoef)
    deallocate(Bcoef)
    deallocate(Ccoef)
    deallocate(Xsol)
    deallocate(Rhs0)
    deallocate(DUfs)
    deallocate(DVfs)
    deallocate(Wfs)

    end subroutine kepsilon_3D

    subroutine kepsilon(ISTEP)
!-------------------------------------------------------
!   k-epsilon turbulence model
!   Last update: Gangfeng Ma, 09/07/2011
!-------------------------------------------------------
    use global
    implicit none
    integer,  intent(in) :: ISTEP
    integer,  parameter :: ke_model = 2
    real(SP), parameter :: Dmin = 0.02
    real(SP), dimension(:,:,:), allocatable :: R5,DelzR,Tke_Old,Eps_Old,DUfs,DVfs,Wfs
    real(SP), dimension(:), allocatable :: Acoef,Bcoef,Ccoef,Xsol,Rhs0
    real(SP) :: c1e,c2e,c3e,cmiu,cfk,cfe,Umag,Zdis,X0,Xa,Xn,FricU,Sche,Schk
    real(SP) :: smax,dmax,c_d,c_1,c_2,c_3,delta_nm,Tkeb,Epsb,Xlfs,Epsfs
    real(SP) :: S11,S22,S33,S12,S13,S23
    integer :: i,j,k,n,m,l,g,IVAR,Iter,Nlen
    
    allocate(R5(Mloc,Nloc,Kloc))
    allocate(DelzR(Mloc,Nloc,Kloc))
    allocate(Tke_Old(Mloc,Nloc,Kloc))
    allocate(Eps_Old(Mloc,Nloc,Kloc))
    allocate(DUfs(Mloc1,Nloc,Kloc))
    allocate(DVfs(Mloc,Nloc1,Kloc))
    allocate(Wfs(Mloc,Nloc,Kloc1))

    ! some parameters
    c1e = 1.44
    c2e = 1.92
!    c3e = -1.4
    c3e = 0.0
    cmiu = 0.09
    cfk = 1.0
    cfe = cfk*c2e/c1e
!    cfe = 1.25
    Sche = 1.3
    Schk = 1.0

    ! save old values
    Tke_Old = Tke
    Eps_Old = Eps

    Prod_s = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DelzU(i,j,k) = (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)
        DelzV(i,j,k) = (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)
        Prod_s(i,j,k) = CmuVt(i,j,k)*(DelzU(i,j,k)**2+DelzV(i,j,k)**2)   
      endif
    enddo
    enddo
    enddo

    ! buoyancy production
    Prod_b = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DelzR(i,j,k) = (Rho(i,j,k+1)-Rho(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)
        Prod_b(i,j,k) = Grav*CmuVt(i,j,k)*DelzR(i,j,k)/Rho0
     endif
    enddo
    enddo
    enddo


    Nlen = Kend-Kbeg+1
    allocate(Acoef(Nlen))
    allocate(Bcoef(Nlen))
    allocate(Ccoef(Nlen))
    allocate(Xsol(Nlen))
    allocate(Rhs0(Nlen))

    ! transport velocities                                                                                               
    DUfs = Ex
    DVfs = Ey
    Wfs = Omega

    ! solve epsilon equation
    IVAR = 2
    call adv_scalar_hlpa(DUfs,DVfs,Wfs,Eps_Old,R5,IVAR)

    do i = Ibeg,Iend
    do j = Jbeg,Jend
      if(D(i,j)<Dmin.and.Mask(i,j)==0) cycle

      Nlen = 0
      do k = Kbeg,Kend
        R5(i,j,k) = R5(i,j,k)+c1e*D(i,j)*(Prod_s(i,j,k)+c3e*Prod_b(i,j,k))*Eps_Old(i,j,k)/Tke_Old(i,j,k)                        
        Nlen = Nlen+1
        if(k==Kbeg) then
          Acoef(Nlen) = 0.0
        else
          Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Sche)/(0.5*dsig(k)*(dsig(k)+dsig(k-1)))
        endif

        if(k==Kend) then
          Ccoef(Nlen) = 0.0
        else
          Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Sche)/(0.5*dsig(k)*(dsig(k)+dsig(k+1)))
        endif
        
        Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)+dt*c2e*Eps_Old(i,j,k)/Tke_Old(i,j,k)

        Rhs0(Nlen) = DEps(i,j,k)+dt*R5(i,j,k)
      enddo

      call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

      Nlen = 0
      do k = Kbeg,Kend
        Nlen = Nlen+1
        DEps(i,j,k) = Xsol(Nlen)
      enddo
    enddo
    enddo

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DEps(i,j,k) = ALPHA(ISTEP)*DEps0(i,j,k)+BETA(ISTEP)*DEps(i,j,k)
        DEps(i,j,k) = max(DEps(i,j,k),D(i,j)*Eps_min)
      endif
    enddo
    enddo
    enddo

    ! slove tke equation
    IVAR = 1
    call adv_scalar_hlpa(DUfs,DVfs,Wfs,Tke_Old,R5,IVAR)

    do i = Ibeg,Iend
    do j = Jbeg,Jend
      if(D(i,j)<Dmin.and.Mask(i,j)==0) cycle

      Nlen = 0
      do k = Kbeg,Kend
        R5(i,j,k) = R5(i,j,k)+D(i,j)*(Prod_s(i,j,k)+Prod_b(i,j,k))-DEps(i,j,k)                                                   
        Nlen = Nlen+1
        if(k==Kbeg) then
          Acoef(Nlen) = 0.0
        else
          Acoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k-1)+Cmu(i,j,k))+0.5*(CmuVt(i,j,k-1)+CmuVt(i,j,k))/Schk)/(0.5*dsig(k)*(dsig(k)+dsig(k-1)))
        endif

        if(k==Kend) then
          Ccoef(Nlen) = 0.0
        else
          Ccoef(Nlen) = -dt/D(i,j)**2*(0.5*(Cmu(i,j,k)+Cmu(i,j,k+1))+0.5*(CmuVt(i,j,k)+CmuVt(i,j,k+1))/Schk)/(0.5*dsig(k)*(dsig(k)+dsig(k+1)))
        endif
        
        Bcoef(Nlen) = 1.0-Acoef(Nlen)-Ccoef(Nlen)
        Rhs0(Nlen) = DTke(i,j,k)+dt*R5(i,j,k)
      enddo
      
      call trig(Acoef,Bcoef,Ccoef,Rhs0,Xsol,Nlen)

      Nlen = 0
      do k = Kbeg,Kend
        Nlen = Nlen+1
        DTke(i,j,k) = Xsol(Nlen)
      enddo
    enddo
    enddo

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        DTke(i,j,k) = ALPHA(ISTEP)*DTke0(i,j,k)+BETA(ISTEP)*DTke(i,j,k)
        DTke(i,j,k) = max(DTke(i,j,k),D(i,j)*Tke_min)
      endif
    enddo
    enddo
    enddo

    ! at the bottom
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      ! impose wall function 
      Umag = sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)
      if(Umag<1.e-6.or.D(i,j)<Dmin.or.Mask(i,j)==0) cycle

      Zdis = 0.5*dsig(Kbeg)*D(i,j)
!      X0 = 0.05
!      Iter = 0
!
!      Xa = dlog(9.0*Umag*Zdis/Visc)
! 10      Xn = X0+(0.41-X0*(Xa+dlog(X0)))/(1.0+0.41/X0)
!      if(Iter>=20) then
!        write(*,*) 'Iteration exceeds 20 steps',i,j,Umag
!      endif
!      if(dabs((Xn-X0)/X0)>1.e-8.and.Xn>0.0) then
!        X0 = Xn
!        Iter = Iter+1
!        goto 10
!      else
!        FricU = Xn*Umag
!      endif

      if(Ibot==1) then
        FricU = sqrt(Cd0)*Umag
      else
        FricU = Umag/(1./Kappa*log(30.*Zdis/Zob))
      endif

      Tkeb = FricU**2/sqrt(cmiu)
      Epsb = FricU**3/(Kappa*Zdis)

      DTke(i,j,Kbeg) = D(i,j)*Tkeb
      DEps(i,j,Kbeg) = D(i,j)*Epsb

      do k = 1,Nghost
        DTke(i,j,Kbeg-k) = DTke(i,j,Kbeg+k-1)
        DEps(i,j,Kbeg-k) = DEps(i,j,Kbeg+k-1)
      enddo
    enddo
    enddo

    ! at the free surface
    do i = Ibeg,Iend
    do j = Jbeg,Jend
!      Xlfs = 0.5*Kappa*D(i,j)*dsig(Kend)
!      Eps(i,j,Kend) = cmiu**(3./4.)*Tke(i,j,Kend)**(3./2.)/Xlfs
!      DEps(i,j,Kend) = D(i,j)*Eps(i,j,Kend)
      do k = 1,Nghost
        DTke(i,j,Kend+k) = DTke(i,j,Kend-k+1)
        DEps(i,j,Kend+k) = DEps(i,j,Kend-k+1)
      enddo
    enddo
    enddo

    call phi_3D_exch(DTke)
    call phi_3D_exch(DEps)

     if(n_west.eq.MPI_PROC_NULL) then
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(WaveMaker(1:3)=='LEF') then
         do i = Ibeg,Ibeg+5
           Tke(i,j,k) = Tke_min
           Eps(i,j,k) = Eps_min
           DTke(i,j,k) = D(i,j)*Tke_min
           DEps(i,j,k) = D(i,j)*Eps_min
           CmuVt(i,j,k) = Cmut_min
         enddo
       endif

       do i = 1,Nghost
         DTke(Ibeg-i,j,k) = DTke(Ibeg+i-1,j,k)
         DEps(Ibeg-i,j,k) = DEps(Ibeg+i-1,j,k)
       enddo
     enddo
     enddo
     endif

     if(n_east.eq.MPI_PROC_NULL) then
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       do i = 1,Nghost
         DTke(Iend+i,j,k) = DTke(Iend-i+1,j,k)
         DEps(Iend+i,j,k) = DEps(Iend-i+1,j,k)
       enddo
     enddo
     enddo
     endif
    
     if(n_suth.eq.MPI_PROC_NULL) then
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         DTke(i,Jbeg-j,k) = DTke(i,Jbeg+j-1,k)
         DEps(i,Jbeg-j,k) = DEps(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
     endif

     if(n_nrth.eq.MPI_PROC_NULL) then
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         DTke(i,Jend+j,k) = DTke(i,Jend-j+1,k)
         DEps(i,Jend+j,k) = DEps(i,Jend-j+1,k)
       enddo
     enddo
     enddo
     endif

    do i = 1,Mloc
    do j = 1,Nloc
    do k = 1,Kloc
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        Tke(i,j,k) = DTke(i,j,k)/D(i,j)
        Eps(i,j,k) = DEps(i,j,k)/D(i,j)
        CmuVt(i,j,k) = Cmiu*Tke(i,j,k)**2/Eps(i,j,k)
      else
        Tke(i,j,k) = Tke_min
        Eps(i,j,k) = Eps_min
        DTke(i,j,k) = D(i,j)*Tke_min
        DEps(i,j,k) = D(i,j)*Eps_min
        CmuVt(i,j,k) = Cmut_min
      endif
    enddo
    enddo
    enddo

    ! no turbulence in the internal wavemaker region
    if(WaveMaker(1:3)=='INT') then
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
          Tke(i,j,k) = Tke_min
          Eps(i,j,k) = Eps_min
          DTke(i,j,k) = D(i,j)*Tke_min
          DEps(i,j,k) = D(i,j)*Eps_min
          CmuVt(i,j,k) = Cmut_min
        endif  
      enddo
      enddo
      enddo
    endif

    ! Reynolds stress (just for output) 
    UpWp = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==1) then
        UpWp(i,j,k) = CmuVt(i,j,k)*(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j) 
      endif
    enddo
    enddo
    enddo

    deallocate(R5)
    deallocate(DelzR)
    deallocate(Tke_Old)
    deallocate(Eps_Old)
    deallocate(Acoef)
    deallocate(Bcoef)
    deallocate(Ccoef)
    deallocate(Xsol)
    deallocate(Rhs0)
    deallocate(DUfs)
    deallocate(DVfs)
    deallocate(Wfs)

    end subroutine kepsilon


    subroutine wall_time_secs(tcurrent)
!--------------------------------------------------------
!   Calculate current wall time
!   Last update: Gangfeng Ma, 09/12/2011
!--------------------------------------------------------
    use global, only: SP
    implicit none
    integer, dimension(8) :: walltime
    real(SP), intent(out) :: tcurrent
    real(SP) :: msecs,secs,mins,hrs,days,months,mscale,years

    call date_and_time(VALUES=walltime)

    msecs = real(walltime(8))
    secs = real(walltime(7))
    mins = real(walltime(6))
    hrs = real(walltime(5))
    days = real(walltime(3))
    months = real(walltime(2))
    years = real(walltime(1))

    if((months.eq.1).or.(months.eq.3).or.(months.eq.5).or.  &
          (months.eq.7).or.(months.eq.8).or.(months.eq.10).or.  &                                                                                   
          (months.eq.12)) then
      mscale = 31.0
    elseif((months.eq.4).or.(months.eq.6).or.  &
          (months.eq.9).or.(months.eq.11)) then
      mscale = 30.0
    elseif(years.eq.4*int(years/4)) then
      mscale = 29.0
    else
      mscale = 28.0
    endif

    tcurrent = months*mscale*24.0*60.0*60.0+days*24.0*60.0*60.0+  &
         hrs*60.0*60.0+60.0*mins+secs+msecs/1000.0

    return
    end subroutine wall_time_secs







    subroutine eval_dens
!---------------------------------------------------------------------
!
!   equation of state
!
!---------------------------------------------------------------------
    use global
    implicit none
    integer, parameter :: kbb = 101
    integer, parameter :: kbbm1 = kbb-1
    real(SP), dimension(kbbm1) :: phy_z,rhoztmp,rhomean
    real(SP), dimension(Mloc,Nloc,kbbm1) :: rhoz
    real(SP), dimension(Kloc) :: zm,rhos
    integer,dimension(1) :: req
    real(SP),dimension(:,:),allocatable :: xx,rhozloc
    real(SP),dimension(Mglob,Nglob,kbbm1) :: rhozglob
    integer,dimension(MPI_STATUS_SIZE,1) :: status
    integer,dimension(NumP) :: npxs,npys
    real(SP) :: DELTZ,TF,SF,RHOF,HMAX,ETAMAX
    real(SP) :: tmp1,tmp2,myvar
    integer :: I,J,K,isum,jk,iglob,jglob,kk,n,len,nreq,NKloc

    ! calculate density from equation of state
    DO K = 1,KLOC
    DO J = 1,NLOC
    DO I = 1,MLOC
      IF(Mask(I,J)==0) cycle

    ENDDO
    ENDDO
    ENDDO

    ! find maximum water depth
    tmp1 = -large
    tmp2 = -large
    do j = 1,Nloc
    do i = 1,Mloc
      if(Mask(i,j)==0) cycle 
      if(hc(i,j)>tmp1) tmp1 = Hc(i,j)  
      if(eta(i,j)>tmp2) tmp2 = Eta(i,j)
    enddo
    enddo
    call MPI_ALLREDUCE(tmp1,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
    hmax = myvar
    call MPI_ALLREDUCE(tmp2,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
    etamax = myvar

    ! interpolate into physical z levels
    deltz = (hmax+etamax)/float(kbbm1)
    do k = 1,kbbm1
      phy_z(k) = (float(k)-0.5)*deltz-hmax
    enddo

    rhoz = Rho0
    do i = 1,Mloc
    do j = 1,Nloc
      if(Mask(i,j)==0) cycle

      do k = 1,Kloc
        zm(k) = sigc(k)*D(i,j)-Hc(i,j)    
        rhos(k) = Rho(i,j,k)
      enddo

      call sinter(zm,rhos,phy_z,rhoztmp,Kloc,kbbm1)

      do k = 1,kbbm1
        rhoz(i,j,k) = rhoztmp(k)
      enddo
    enddo
    enddo

    call MPI_GATHER(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,  &
           0,MPI_COMM_WORLD,ier)
    call MPI_GATHER(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,  &
           0,MPI_COMM_WORLD,ier)
    
    NKloc = Nloc*kbbm1

    ! put the data in master processor into the global var                                                           
    if(myid==0) then
      do k = 1,kbbm1
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        iglob = i-Nghost
        jglob = j-Nghost
        rhozglob(iglob,jglob,k) = rhoz(i,j,k)
      enddo
      enddo
      enddo
    endif

    allocate(rhozloc(Mloc,NKloc))
    allocate(xx(Mloc,NKloc))

    do k = 1,kbbm1
    do j = 1,Nloc
    do i = 1,Mloc
      jk = (k-1)*Nloc+j
      rhozloc(i,jk) = rhoz(i,j,k)
    enddo
    enddo
    enddo

    ! collect data from other processors into the master processor                                                   
    len = Mloc*NKloc

    do n = 1,NumP-1
      if(myid==0) then
        call MPI_IRECV(xx,len,MPI_SP,n,0,MPI_COMM_WORLD,req(1),ier)
        call MPI_WAITALL(1,req,status,ier)
        do k = 1,kbbm1
        do j = Jbeg,Jend
        do i = Ibeg,Iend
          iglob = npxs(n+1)*(Iend-Ibeg+1)+i-Nghost
          jglob = npys(n+1)*(Jend-Jbeg+1)+j-Nghost
          jk = (k-1)*Nloc+j
          rhozglob(iglob,jglob,k) = xx(i,jk)
        enddo
        enddo
        enddo
      endif

      if(myid==n) then
        call MPI_SEND(rhozloc,len,MPI_SP,0,0,MPI_COMM_WORLD,ier)
      endif
    enddo

    deallocate(rhozloc)
    deallocate(xx)

    if(myid==0) then
      rhomean = zero
      do k = 1,kbbm1
        isum = 0
        do j = 1,Nglob
        do i = 1,Mglob
          if(-HCG(i,j)<=phy_z(k)) then
            isum = isum+1
            rhomean(k) = rhomean(k)+rhozglob(i,j,k)
          endif
        enddo
        enddo
        if(isum>=1) then
          rhomean(k) = rhomean(k)/float(isum)
        else
          rhomean(k) = rhomean(k-1)
        endif
      enddo
    endif

    call MPI_BCAST(rhomean,kbbm1,MPI_SP,0,MPI_COMM_WORLD,ier)


    ! linearly interpolate to obtain density at signa levels
    Rmean = Rho0
    do i = 1,Mloc
    do j = 1,Nloc
      if(Mask(i,j)==0) cycle

      do k = 1,Kloc
        zm(k) = sigc(k)*D(i,j)-Hc(i,j)
      enddo

      call sinter(phy_z,rhomean,zm,rhos,kbbm1,Kloc)          

      Rmean(i,j,1:Kloc) = rhos
    enddo
    enddo
  
    return  
    end subroutine eval_dens

    subroutine sinter(X,A,Y,B,M1,N1)
!------------------------------------------------------------------------------
!                                                                              
!  this subroutine linearly interpolates and extrapolates an                             
!  array b.
!                                                                              
!  x(m1) must be ascending                                                    
!  a(x) given function                                                         
!  b(y) found by linear interpolation and extrapolation                        
!  y(n1) the desired depths                                                    
!  m1   the number of points in x and a                                        
!  n1   the number of points in y and b                                        
!                                                                              
!  a special case of interp ....no extrapolation below data                    
!
!----------------------------------------------------------------------
    use global, only: SP
    implicit none
    INTEGER, INTENT(IN)  :: M1,N1
    REAL(SP),  INTENT(IN)  :: X(M1),A(M1),Y(N1)
    REAL(SP),  INTENT(OUT) :: B(N1)
    INTEGER :: I,J,NM        
!                                                                                                                    
!   EXTRAPOLATION                                                                                                     
!                                                                                                                    
    DO I=1,N1
      IF (Y(I)<X(1 )) B(I) = A(1)
      IF (Y(I)>X(M1)) B(I) = A(M1)
    END DO

!                                                                                                                    
!   INTERPOLATION                                                                                                     
!                                                                                                                    
    NM = M1 - 1
    DO I=1,N1
      DO J=1,NM
        IF (Y(I)>=X(J).AND.Y(I)<=X(J+1)) &
           B(I) = A(J+1) - (A(J+1)- A(J)) * (X(J+1)-Y(I)) / (X(J+1)-X(J))
      END DO
    END DO

    return
    end subroutine sinter


    subroutine sinter_p(X,A,Y,B,M1,N1)
!------------------------------------------------------------------------------
!                                                                              
!  for baroclinic interpolation                                               
!                                                                              
!  this subroutine linearly interpolates and extrapolates an                   
!  array b.                                                                    
!                                                                              
!  x(m1) must be ascending                                                    
!  a(x) given function                                                         
!  b(y) found by linear interpolation and extrapolation                        
!  y(n1) the desired depths                                                    
!  m1   the number of points in x and a                                        
!  n1   the number of points in y and b                                        
!                                                                              
!  a special case of interp ....no extrapolation below data                    
!
!----------------------------------------------------------------------
    use global, only: SP
    implicit none
    INTEGER, INTENT(IN)  :: M1,N1
    REAL(SP),  INTENT(IN)  :: X(M1),A(M1),Y(N1)
    REAL(SP),  INTENT(OUT) :: B(N1)
    INTEGER :: I,J,NM        
!                                                                                                                    
!   EXTRAPOLATION                                                                                                     
!                                                                                                                    
    DO I=1,N1
      IF(Y(I) < X(1 )) B(I) = A(1)-(A(2)-A(1))*(X(1)-Y(I))/(X(2)-X(1))
      IF(Y(I) > X(M1)) B(I) = A(M1)+(A(M1)-A(M1-1))*(Y(I)-X(M1))/(X(M1)-X(M1-1))                            
    END DO

!                                                                                                                    
!   INTERPOLATION                                                                                                     
!                                                                                                                    
    NM = M1 - 1
    DO I=1,N1
      DO J=1,NM
        IF (Y(I)>=X(J).AND.Y(I)<=X(J+1)) &
           B(I) = A(J+1) - (A(J+1)- A(J)) *(X(J+1)-Y(I)) / (X(J+1)-X(J))
      END DO
    END DO

    return
    end subroutine sinter_p


    subroutine baropg_z
!------------------------------------------------------------------------------
!   Calculate baroclinic terms in z levels
!   Called by
!      eval_duvw
!   Last Update: Gangfeng Ma, 07/05/2012
!-----------------------------------------------------------------------------
    use global
    implicit none
    integer, parameter :: kbb = 201
    integer, parameter :: kbbm1 = kbb-1
    real(SP), dimension(Kloc) :: zm,rhos,pbxs,pbys
    real(SP), dimension(kbb) :: phy_z,rhoztmp,pbx,pby
    real(SP), dimension(3,3,kbb) ::pb,rhoz
    real(SP) :: Ramp1,hmax,etamax,tmp1,tmp2,myvar,deltz, &
                Rmean1,Rmean2,dz
    integer :: i,j,k,i2,j2,ic,jc

    if(TRamp==Zero) then
      Ramp1 = 1.0
    else
      Ramp1 = tanh(TIME/TRamp)
    endif

    ! subtract reference density
    Rho = Rho-Rho0

    DRhoX = Zero; DRhoY = Zero
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask9(i,j)==0) cycle

      ! find local maximum water depth 
      hmax = -large
      etamax = -large
      do j2 = j-1,j+1
      do i2 = i-1,i+1
        if(hc(i2,j2)>hmax) hmax = Hc(i2,j2)
        if(eta(i2,j2)>etamax) etamax = Eta(i2,j2)
      enddo
      enddo

      ! interpolate into physical z levels
      if(Ivgrd==1) then  
        deltz = (hmax+etamax)/float(kbbm1)
        do k = 1,kbb
          phy_z(k) = (float(k)-1.0)*deltz-hmax
        enddo
      else
        deltz = (hmax+etamax)*(Grd_R-1.0)/(Grd_R**float(kbbm1)-1.0)
        phy_z(1) = -hmax
        do k = 2,kbb
          phy_z(k) = phy_z(k-1)+deltz
          deltz = deltz*Grd_R
        enddo
      endif

      rhoz = Zero
      do i2 = i-1,i+1
      do j2 = j-1,j+1
        ic = i2-i+2  ! local index
        jc = j2-j+2  
        do k = 1,Kloc
          zm(k) = sigc(k)*D(i2,j2)-Hc(i2,j2)
          rhos(k) = Rho(i2,j2,k)
        enddo

        call sinter(zm,rhos,phy_z,rhoztmp,Kloc,kbb)

        do k = 1,kbb
          rhoz(ic,jc,k) = rhoztmp(k)
        enddo
      enddo
      enddo

      pb = Zero
      do i2 = i-1,i+1
      do j2 = j-1,j+1
        ic = i2-i+2
        jc = j2-j+2
        do k = kbbm1,1,-1
          if(phy_z(k)>=Eta(i2,j2)) then
            pb(ic,jc,k) = 0.0
          else
            if(phy_z(k)>=-Hc(i2,j2)) then
              dz = dmin1(phy_z(k+1)-phy_z(k),Eta(i2,j2)-phy_z(k))
              pb(ic,jc,k) = pb(ic,jc,k+1)+0.5*(Rhoz(ic,jc,k)+Rhoz(ic,jc,k+1))*dz
            elseif(phy_z(k)<-Hc(i2,j2).and.phy_z(k+1)>-Hc(i2,j2)) then
              dz = -Hc(i2,j2)-Phy_z(k)
              pb(ic,jc,k) = pb(ic,jc,k+1)+0.5*(Rhoz(ic,jc,k)+Rhoz(ic,jc,k+1))*dz
            else
              pb(ic,jc,k) = pb(ic,jc,k+1)
            endif
          endif
        enddo
      enddo
      enddo

      pbx = Zero; pby = Zero
      do k = 1,kbb
        if(phy_z(k)<=Eta(i,j).and.phy_z(k)>=-Hc(i,j)) then
          if(phy_z(k)<-Hc(i-1,j).and.phy_z(k)>=-Hc(i+1,j)) then
            pbx(k) = (pb(3,2,k)-pb(2,2,k))/dx
          elseif(phy_z(k)>=-Hc(i-1,j).and.phy_z(k)<-Hc(i+1,j)) then
            pbx(k) = (pb(2,2,k)-pb(1,2,k))/dx
          elseif(phy_z(k)<-Hc(i-1,j).and.phy_z(k)<-Hc(i+1,j)) then
            pbx(k) = Zero
          else
            pbx(k) = (pb(3,2,k)-pb(1,2,k))/(2.0*dx)
          endif

          if(phy_z(k)<-Hc(i,j-1).and.phy_z(k)>=-Hc(i,j+1)) then
            pby(k) = (pb(2,3,k)-pb(2,2,k))/dy
          elseif(phy_z(k)>=-Hc(i,j-1).and.phy_z(k)<-Hc(i,j+1)) then
            pby(k) = (pb(2,2,k)-pb(2,1,k))/dy
          elseif(phy_z(k)<-Hc(i,j-1).and.phy_z(k)<-Hc(i,j+1)) then
            pby(k) = Zero
          else
            pby(k) = (pb(2,3,k)-pb(2,1,k))/(2.0*dy)
          endif
        endif
      enddo

      do k = 1,Kloc
        zm(k) = sigc(k)*D(i,j)-Hc(i,j)
      enddo
 
      call sinter_p(phy_z,pbx,zm,pbxs,kbb,Kloc)
      call sinter_p(phy_z,pby,zm,pbys,kbb,Kloc)

      do k = Kbeg,Kend
        DRhoX(i,j,k) = -pbxs(k)*grav*D(i,j)/Rho0*Ramp1
        DRhoY(i,j,k) = -pbys(k)*grav*D(i,j)/Rho0*Ramp1
      enddo
    enddo
    enddo

    ! Add back reference density
    Rho = Rho+Rho0

    return
    end subroutine baropg_z


    subroutine trig(alpha,beta,gama,b,x,N)
!*************************************************************!
!*                                                           *!
!*         (B1 C1                         )   (x1)     (b1)  *!
!*         (A2 B2 C2                      )   (x2)     (b2)  *!
!*         (   A3 B3 C3                   )   (x3)     (b3)  *!
!*         (      A4 B4 C4                )   (x4)====       *!
!*         (         A5 B5  C5            )   ... ====  ...  *!
!*         (            ... ... ...       )   ...       ...  *!
!*         (                An-1 Bn-1 Cn-1)   (xn-1)   (bn-1)*!
!*         (                     An   Bn  )   (xn)     (bn)  *!
!*                                                           *!
!*                                                           *!
!*************************************************************!
! where A are alpha, B are beta, C are gama
!-----------------------------------------------------------------------------
    use global, only: SP
    implicit none

    integer, intent(in) :: N
    real(SP), dimension(N), intent(in)  :: alpha,beta,gama,b
    real(SP), dimension(N), intent(out) :: x
    real(SP), dimension(N) :: betaPrime,bPrime
    real(SP) :: coeff
    integer :: II
 
    ! Perform forward elimination
    betaPrime(1) = beta(1)
    bPrime(1) = b(1)
 
    do II = 2,N
      coeff = alpha(II)/betaPrime(II-1)
      betaPrime(II) = beta(II)-coeff*gama(II-1)
      bPrime(II) = b(II)-coeff*bPrime(II-1)
    enddo

    ! Perform back substitution
    x(N) = bPrime(N) / betaPrime(N)
    do II = N-1,1,-1
      x(II) = (bPrime(II)-gama(II)*x(II+1))/betaPrime(II)
    enddo

    end subroutine trig




!------------------------The End----------------------------------------------
