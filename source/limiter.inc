!     real(SP) :: LIMITER, A, B
     real(SP) :: A, B
     LIMITER(A, B)=(A*ABS(B)+ABS(A)*B)/(ABS(A)+ABS(B))


!     FUNCTION LIMITER(A,B)
!     use global, only: SP,Zero,One,Small
!     IMPLICIT NONE
!     REAL(SP),INTENT(IN) :: A,B
!     REAL(SP) :: LIMITER
!
!!     ! minmod limiter
!!     LIMITER=max(Zero,min(A,B))
!
!!     ! van Leer limiter
!     LIMITER=(A*ABS(B)+ABS(A)*B)/(ABS(A)+ABS(B))
!
!!     ! superbee limiter
!!     LIMITER=SIGN(One,B)*MAX(Zero,MIN(2.0*ABS(B),SIGN(One,B)*A),  &
!!          MIN(ABS(B),2.0*SIGN(One,B)*A))
!
!     RETURN
!     END FUNCTION LIMITER
