program main
  IMPLICIT NONE
  DOUBLE PRECISION :: h
  INTEGER :: cantk, n1, nang, cant1, sub, j, i, autovector
  DOUBLE PRECISION, ALLOCATABLE :: A(:,:),MAVe(:,:), MAVa(:)
  CHARACTER (len=4) :: filename

!!!!!!!!!!!CREDITOS!!!!!!!!!!!!!!
print*,''
print*,''
print*,''
print*, '                **********FIGURAS DE CHLADNI*************'
print*, 'Este programa resuelve el problema de autovalores del operador biarmonico'
print*, 'mediante diferencias finitas. Esta solucion corresponde  a encontrar  los'
print*, 'nodos de vibración de una onda estacionaria sobre una placa circular, con'
print*, 'condiciones de Dirichlet. Esto es, encontrar las figuras de Chladni para'
print*, 'alguna frecuencia de vibracion.'
print*, ''
print*, 'AUTORES'
print*, '-Cardenas Andrade Luis Daniel.'
print*, '-Pachas Yeren Valeria Sofia.'
print*, '-Dr. Paredes Cabrel Alejandro.'
print*, '-Dr. Beltran Ramirez Jhosep.'
print*,''
print*,''
print*,''
print*,''
print*,'---------------------------'
print*,'Que autovector desea ver:'
read(*,*) autovector
print*,'---------------------------'



!lectura de los parametros!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 OPEN(UNIT = 20, FILE='main_parametros.txt',STATUS="UNKNOWN",ACTION="read")
 read(20,*) sub   !que metodo usar, 1  es para el QR y 0 es para el MKL
 read(20,*) n1    !numero de puntos discretizados del radio.
 read(20,*) cant1 !cantidad de iteraciones para el calculo QR
 read(20,*) cantk !ultimo valor de K, Los valores de K de la serie de fourier van de K =0,1,2,3,...,cantk | total de valores(cantk+1)
 read(20,*) nang  !numero de puntos discretizados del angulo (valor par).
 close(20) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 

  allocate(A(n1*(cantk+1),n1*(cantk+1)))!matriz de factores de la discretizacion del biarmonico
  allocate(MAVe(n1*(cantk+1),n1*(cantk+1)))!matriz de autovectores
  allocate(MAVa(n1*(cantk+1)))!array de autovalores
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print*, '-----------------------------------------'
print*, 'entrando a la matriz de representacion'
!!!!!!!!LLAMADO A a la subrutina para la obtencion de la matriz de la representación del biarmonico!!!!!!!!  
  
  call matrix_repre_final(n1,nang,cantk,A)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print*, 'paso sin problema la matriz de representación'
print*, '-----------------------------------------------'
print*, ' '
print*, ' '
print*, ' '
!!!!!!!!LLAMADO A LOS METODOS PARA LA OBTENCION DE AUTOVECTORES Y AUTOVALORES!!!!!!!!!  
  if (sub == 1) then
     call metodoQR(A, cant1, n1*(cantk+1), MAVe, MAVa)
     !MAVe: matriz autovectores
     !MAVa: matriz autovalores
     filename ='QR'         
  else if (sub == 0) then 
     call metodoMKL(A, n1*(cantk+1), MAVe, MAVa)
     filename='MKL' 
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print*, ' '
print*, ' '
print*, ' '

!-----------------------------------------------------------------------------------------------------------------------------            
!-----------------------------------------------------------------------------------------------------------------------------
  open(unit=40,file='test_Matrix_Autovecto.txt',STATUS="UNKNOWN",ACTION="WRITE")!solo para ver la forma de los autovalores      
   !matriz autovalores
   do i = 1, n1*(cantk+1)
     write(40,*) MAVe(i,autovector)  !write(40,*) MAVe(i,j) observar el autovector ubicado en la posicion j
   end do
   close(40)
!-----------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------

call data_saving(n1*(cantk+1), filename, MAVe, MAVa)!subrrutina para guardar los datos obtenidos en un archivo .dat


  deallocate(A)!matriz de factores de la discretizacion del biarmonico
  deallocate(MAVe)!matriz de autovectores del metodo QR
  deallocate(MAVa)!matriz de autovalores del metodo QR

  print*, 'finished :)'
end program








