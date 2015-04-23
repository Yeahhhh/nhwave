! hackathon profiling self 28.33%
     subroutine hypre_pres_solver
!---------------------------------------------
!    solve pressure using hypre package
!    called by
!       poisson_solver
!    Last update: 22/08/2011, Gangfeng Ma
!---------------------------------------------
# if defined (PARALLEL)
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

# endif
     return
     end subroutine hypre_pres_solver


