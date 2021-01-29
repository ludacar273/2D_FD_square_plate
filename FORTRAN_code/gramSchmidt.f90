SUBROUTINE gran_schmidt( A, m, n, U )
   INTEGER, INTENT(IN) :: m,n
   DOUBLE PRECISION, INTENT(IN) :: A(m,n)
   DOUBLE PRECISION, INTENT(OUT) :: U(m,n) 

   INTEGER :: i,j
   DOUBLE PRECISION :: S(m)


   U(:,1) = A(:,1)
   do i = 2, n
 
      S = 0.0
      do j = 1, i-1
         S = S + ( sum( A(:,i)*U(:,j) )/sum( U(:,j)*U(:,j) )  )*U(:,j)        
      end do

      U(:,i) = A(:,i) - S  
      
   end do  
   
END SUBROUTINE gran_schmidt 
