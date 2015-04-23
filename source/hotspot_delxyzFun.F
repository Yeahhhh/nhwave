!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D

     subroutine delxFun_2D(Din,Dout)
!-------------------------------------------
!    Second-order derivative in x
!    Called by
!       delxyzFun
!    Last update: 04/01/2011, Gangfeng Ma
!------------------------------------------
     use global !, only: SP,Small,Zero,dx,Mloc,Nloc,Mask,Brks
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc) :: Dout
     real(SP) :: TMP1,TMP2,LIMITER
     integer :: i,j
#include "limiter.inc"

     do j = 1,Nloc
     do i = 2,Mloc-1
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
#include "limiter.inc"

     do j = 2,Nloc-1
     do i = 1,Mloc
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





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D





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

#include "limiter.inc"


! !$acc kernel

     do k = 1,Kloc
     do j = 1,Nloc
     do i = 2,Mloc-1
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

     do k = 1,Kloc
     do j = 1,Nloc
       Dout(1,j,k) = (Din(2,j,k)-Din(1,j,k))/dx
       Dout(Mloc,j,k) = (Din(Mloc,j,k)-Din(Mloc-1,j,k))/dx
     enddo
     enddo

! !$acc end kernel

     return
     end subroutine delxFun_3D







! hackathon profiling self 4.36%
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
#include "limiter.inc"

     do k = 1,Kloc
     do j = 2,Nloc-1
     do i = 1,Mloc
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

     do k = 1,Kloc
     do i = 1,Mloc
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
     real(SP) :: TMP1,TMP2 ,LIMITER
     integer :: i,j,k
#include "limiter.inc"

     do k = 2,Kloc-1
     do j = 1,Nloc
     do i = 1,Mloc
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

     do j = 1,Nloc
     do i = 1,Mloc
       Dout(i,j,1) = (Din(i,j,2)-Din(i,j,1))/(0.5*(dsig(1)+dsig(2)))
       Dout(i,j,Kloc) = (Din(i,j,Kloc)-Din(i,j,Kloc-1))/(0.5*(dsig(Kloc-1)+dsig(Kloc)))
     enddo
     enddo

     return
     end subroutine delzFun_3D






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



