SUBROUTINE matrix_repre_final(n,nangular,k, A )!k ultimo valor de los k de fourier, desde 0 hasta k
   INTEGER, INTENT(IN) :: n, nangular, k
   DOUBLE PRECISION, INTENT(OUT) :: A(n*(k+1),n*(k+1)) !(k+1) por que incluye al cero

   DOUBLE PRECISION :: B(n,n)
   INTEGER :: i


   call inicializarmatriz(1, n*(k+1), n*(k+1), A)!llenado de la matriz A, dimension (n(k+1)^2*n(k+1)^2), de puros ceros. revisar subrutina 'inicializarmatriz'.
   
   
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Tensor interno (completa) 
 

do l = 0, k

  call matrix_repre_aux(n,nangular,l, B)

  do i =1, n 
    do j =1, n 
       A(n*l + i, n*l + j)=B(i,j)
    end do
  end do 

end do
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!  verificacion simetria !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !comprobacion si la matriz es simetrica
  var = 1
  do i = 1, n*(k+1) - 1
    do j = i+1, n*(k+1)
      if ( A(i,j).NE.A(j,i) ) then
        var = 0
      end if
    end do
  end do
 
  if (var == 1) then
    print*, 'la matriz es simetrica'
  else if (var==0) then
    print*, 'la matriz no es simetrica'
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(unit=24,file='A.txt',STATUS="UNKNOWN",ACTION="WRITE")
      do  i= 1, n*(k+1)
        write(24,*) A(i,:)
      end do
  close(24)

 
     
END SUBROUTINE matrix_repre_final
