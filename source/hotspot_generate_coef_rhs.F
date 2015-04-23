
! hackathon profiling self 7.32%
     subroutine generate_coef_rhs
!---------------------------------------------
!    Generate coefficient matrix and rhs
!    Called by 
!       poisson_solver
!    Last update: 24/03/2011, Gangfeng Ma
!--------------------------------------------
     use global
     use timing

     implicit none
     integer :: i,j,k,neq,n,ic
     real(SP), dimension(:,:,:), allocatable :: DelxS,DelyS,DelzS,A1
     integer,  dimension(:,:,:), allocatable :: indx

     real(kind(0.0d0)) :: t1, t2

     allocate(DelxS(Mloc,Nloc,Kloc1))
     allocate(DelyS(Mloc,Nloc,Kloc1))
     allocate(DelzS(Mloc,Nloc,Kloc1))
     allocate(A1(Mloc,Nloc,Kloc1))
     allocate(indx(Mloc,Nloc,Kloc))

     call TimeNow (t1)

!$acc kernels
     DelxS = Zero
     DelyS = Zero
     DelzS = Zero
     A1 = Zero
!$acc end kernels
!$acc wait


     do k = Kbeg,Kend1
!$acc kernels async(k)
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       DelxS(i,j,k) = (1.-sig(k))/D(i,j)*DelxH(i,j)-sig(k)/D(i,j)*DelxEta(i,j)
       DelyS(i,j,k) = (1.-sig(k))/D(i,j)*DelyH(i,j)-sig(k)/D(i,j)*DelyEta(i,j)
       DelzS(i,j,k) = 1./D(i,j)

       A1(i,j,k) = DelxS(i,j,k)*DelxS(i,j,k)+DelyS(i,j,k)*DelyS(i,j,k)+  &
            DelzS(i,j,k)*DelzS(i,j,k)
     enddo
     enddo
!$acc end kernels
     enddo
!$acc wait
   
     call TimeNow (t2)
     tt(1) = tt(1) + t2 - t1



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
     do k = Kbeg,Kend
     do j = Jbeg+1,Jend-1
     do i = Ibeg+1,Iend-1
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

# if defined (OBSTACLE)
     do k = Kbeg+1,Kend-1
     do j = Jbeg+1,Jend-1
     do i = Ibeg+1,Iend-1
       if(set_flag(i,j,k)==1) then
         ! left 
         if(set_flag(i+1,j,k)==0) then
           ic = indx(I+1,j,k)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,9)
           Coef(ic,9) = Zero
         ! right 
         elseif(set_flag(i-1,j,k)==0) then
           ic = indx(I-1,j,k)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,2)
           Coef(ic,2) = Zero
         ! south
         elseif(set_flag(i,j+1,k)==0) then
           ic = indx(i,J+1,k)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,10)
           Coef(ic,10) = Zero
         ! north
         elseif(set_flag(i,j-1,k)==0) then
           ic = indx(i,J-1,k)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,3)
           Coef(ic,3) = Zero
         ! bottom
         elseif(set_flag(i,j,k+1)==0) then
           ic = indx(i,j,k+1)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,13)
           Coef(ic,13) = Zero
         ! upper
         elseif(set_flag(i,j,k-1)==0) then
           ic = indx(i,j,k-1)
           Coef(ic,1) = Coef(ic,1)+Coef(ic,6)
           Coef(ic,6) = Zero
         endif  
       endif
     enddo
     enddo
     enddo
# endif

     ! boundary conditions
     ! left side
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
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
# if defined (PARALLEL)
     endif
# endif

     ! right side
# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
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
# if defined (PARALLEL)
     endif
# endif

     ! front side
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
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
# if defined (PARALLEL)
     endif
# endif

     ! back side
# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
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
# if defined (PARALLEL)
     endif
# endif

     ! bottom side
     k = Kbeg
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       ic = indx(i,j,k)
# if defined (LANDSLIDE)
       if(trim(SlideType)=='RIGID') then
         Rhs(ic) = Rhs(ic)+Rho0*(dsig(Kbeg)+dsig(Kbeg-1))*(Coef(ic,13)*D(i,j)*Delt2H(i,j)+ &
            Coef(ic,12)*D(i+1,j)*Delt2H(i+1,j)+Coef(ic,11)*D(i,j+1)*Delt2H(i,j+1)+ &
            Coef(ic,14)*D(i-1,j)*Delt2H(i-1,j)+Coef(ic,15)*D(i,j-1)*Delt2H(i,j-1))
       endif
# endif

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

