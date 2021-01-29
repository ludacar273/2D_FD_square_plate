
subroutine metodoMKL(B, N, VRO,WO)

      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: B(N,N)  
      INTEGER, INTENT(IN) :: N
      DOUBLE PRECISION, INTENT(OUT) :: VRO( N, N ), WO(N)
      !REAL :: LDA, LDVL, LDVR
!     .. Local Scalars ..
      INTEGER          INFO, LWORK, LWMAX, i, j
      PARAMETER (LWMAX = 100000000)

!     .. Local Arrays ..
      DOUBLE PRECISION A( N, N ), VL( N, N ),WR( N ), WI( N ), W(N), WR1(N), VR1(N,n), VR(N,N), aux1(N)
      DOUBLE PRECISION, ALLOCATABLE :: WORK( : )
      DOUBLE PRECISION aux 
      
      
      allocate( WORK(LWMAX))

      A(:,:) = B(:,:)
print*, '----------------------'
print*,'running MKL Method'
! Computes optimal value of LWORK
LWORK = -1
      CALL DGEEV( 'Vectors', 'Vectors', N, A, N, WR, WI, VL, N, VR, N, WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      !CALL SGEEV( 'Vectors', 'Vectors', N, A, N, WR, WI, VL, N, VR, N, WORK, LWORK, INFO )
      !LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
!print*,LWORK
!print *,A
       !CALL SGEEV( 'V', 'V', N, A, N, WR, WI, VL, N, VR, N, WORK, LWORK, INFO )
      CALL DGEEV( 'V', 'V', N, A, N, WR, WI, VL, N, VR, N, WORK, LWORK, INFO )

print*, 'MKL method have finished'
print*, '---------------------------'

!                        **************ATENCION***************

VR1(:,:)=VR(:,:)
WR1(:)=WR(:)

!ordenado del los autovalores solo para la comparacion de este metodo con algun otro 
!esto solo es para la comprobacion visual, tener en cuenta que quiza ambos metodos no tengan
!sus autovalores ordenados de la isma forma y con eso sus autovectores sufren el mismo desorden.
do j = 1, N-1
  do i = 1, N-j
    if(WR1(i) > WR1(i+1) )then      !ordenamiento de menor a mayor por metodo de la burbuja

      aux=WR1(i)
      aux1(:)=VR1(:,i)
                  
      WR1(i)=WR1(i+1)             
      VR1(:,i)=VR1(:,i+1)             

      WR1(i+1)=aux
      VR1(:,i+1)=aux1(:)

    end if
  end do               
end do

do i = 1, N    !inversion de WR1 y VR1
  WO(i)=WR1(N+1-i)
  VRO(:,i)=VR1(:,N+1-i)
end do
!***************************************************************************************




9998 FORMAT( 11(:,1X,F6.2) )

deallocate( WORK)

END SUBROUTINE metodoMKL
