SUBROUTINE data_saving(n, filename0, MAVe, MAVa)
 
   INTEGER, INTENT(IN) :: n
   DOUBLE PRECISION, INTENT(IN) :: MAVe(n,n), MAVa(n)
   CHARACTER (len=4), INTENT(IN) :: filename0
   INTEGER :: i
   DOUBLE PRECISION, ALLOCATABLE :: auxMAVe(:,:)
   CHARACTER(len=30):: filename1


   allocate( auxMAVe(1,n**2) )
   auxMAVe = reshape(MAVe,(/1, n**2/) )
   !auxMAVe=reshape(MAVe,(/n, n/) )


   filename1 = filename0 // '_Mat_AutoVec.dat'

   open(unit=32,file=filename1,STATUS="UNKNOWN",ACTION="WRITE")
   !La matriz autovectores
   do i = 1, n**2
     
     write(32,9998)  auxMAVe(1,i)  !si quieres quitarle el formato por que los numeros son pequeños solo borra el 9998
    
   end do 
   close(32)
   
   filename1 = filename0 // '_Mat_AutoValo.dat'

   open(unit=34,file=filename1,STATUS="UNKNOWN",ACTION="WRITE")      
   !matriz autovalores
   do i = 1, n
     write(34,9998) MAVa(i)!si quieres quitarle el formato por que los numeros son pequeños solo borra el 9998
   end do
   close(34)


   deallocate (auxMAVe)

   9998 FORMAT( 1000(:,1X,F20.4) )

END SUBROUTINE data_saving
