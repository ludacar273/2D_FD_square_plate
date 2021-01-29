SUBROUTINE inicializarmatriz(cond, n, m, A) 
   INTEGER, INTENT(IN) :: n, m, cond
   DOUBLE PRECISION, INTENT(OUT) :: A(n,m) 

   INTEGER :: i,j
 
   if(cond==1)then ! si cond=1, se inicializa una matriz de ceros. 
     do i = 1, n
        do j =1, m
           A(i,j)=0.d0
        end do
     end do
   else if (cond==2) then! si cond=2, se inicializa una matriz identidad. 
    
     do i =1, n
       do j =1, m
         if (i == j)then
           A(i,j) = 1.d0
         else
           A(i,j) = 0.d0
         end if
       end do 
     end do
   end if
   
END SUBROUTINE inicializarmatriz
