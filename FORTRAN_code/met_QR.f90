subroutine metodoQR(A,cant, N, Qf,Af)

   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(IN) :: A(N,N)  
   INTEGER, INTENT(IN) :: N
   INTEGER, INTENT(IN) :: cant
   DOUBLE PRECISION, INTENT(OUT) :: Qf(N,N), Af(N)
   INTEGER :: i,j,x, var1
   
   DOUBLE PRECISION, ALLOCATABLE :: U(:,:), Q(:,:), Qt(:,:), Ak(:,:), Qaux(:,:)	

   allocate(Ak(N,N))
   allocate(Qaux(N,N))
 
   Ak(:,:)=A(:,:)!la matriz A no puede sufrir cambios ya que es un OUT, por eso guardamos los datos en otra matriz
   
   !valores iniciales
   x = 0!numero de iteraciones
   var1 = 1!condición de Loop
   
   call inicializarmatriz(2, N, N, Qaux)!inicializamos una matriz como la identidad, para posterior multiplicacion con ella.
  
   DO WHILE (var1 == 1)
      
      allocate(U(N,N))

      call gran_schmidt( Ak, N, N, U ) !ortonormalización de gran_schmidt
       
      allocate(Q(N,N))
      !calculo de matriz ortonormal!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      do i = 1, N
         Q(:,i) = U(:,i) / sqrt( sum( U(:,i)*U(:,i) ) ) !Q matriz ortonormal
      end do
      deallocate(U)

      allocate(Qt(N,N))
      Qt=transpose(Q)
 
      !la nueva matriz Ak = Qt*A'*Q
      Ak=matmul(Qt,matmul(Ak,Q))
      deallocate(Qt)      

      Qaux=matmul(Qaux,Q)
      deallocate(Q)      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      x = x + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (cant == x)  then
         var1 = 0

         Qf(:,:)=Qaux(:,:)
         do i = 1, N
            Af(i)=Ak(i,i) 
         end do
        
         deallocate(Ak)
         deallocate(Qaux)        
      end if 
     
  END DO

9998 FORMAT( 11(:,1X,F6.2) )
END SUBROUTINE metodoQR
