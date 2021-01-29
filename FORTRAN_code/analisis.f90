program analisis

CHARACTER (len=4) :: filename
CHARACTER (len=30):: filename1
INTEGER :: i, j, n1, k, r, sub, cant1, cantk, nang, nav
DOUBLE PRECISION, ALLOCATABLE :: MAVe(:,:), aux1MAVe(:,:), aux2MAVe(:,:), aux3MAVe(:,:)
DOUBLE PRECISION, ALLOCATABLE :: aux4MAVe(:,:), O(:,:), Ob(:,:)
DOUBLE PRECISION :: hr, ho, ir, jr, tetap, pi, c

 OPEN(UNIT = 20, FILE='main_parametros.txt',STATUS="UNKNOWN",ACTION="read")
 read(20,*) sub   !que metodo usar, 1 es para el QR y 0 es para el MKL
 read(20,*) n1    !numero de puntos de discretizacion del radio. 
 read(20,*) cant1 !cantidad de iteraciones para el calculo QR
 read(20,*) cantk !cantidad de valores de K
 read(20,*) nang  !numero de puntos discretizados del angulo.
 close(20)

if (sub==1) then
  filename = 'QR'
else if (sub == 0)then 
  filename = 'MKL'
end if

!lectura de los parametros para el analisis de resultados!!!!!!!!!!!!!!!!!!!!
! OPEN(UNIT = 18, FILE= 'analisis_parametros.txt',STATUS="UNKNOWN",ACTION="read")
! read(18,*) k2!numero de autovalor 
! CLOSE(18)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!lectura de los datos!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!allocate(MAVa(n**2))
allocate(aux1MAVe(1,(n1*(cantk+1))**2 ))

!filename1 = filename//'_Mat_AutoValo.dat'
!OPEN(UNIT = 10, FILE= filename1,STATUS="UNKNOWN",ACTION="read")

!do i = 1, n**2
!    read(10,*) MAVa(i)
!end do

 filename1 = filename//'_Mat_AutoVec.dat'
 OPEN(UNIT = 12, FILE= filename1,STATUS="UNKNOWN",ACTION="read")

 do i = 1, (n1*(cantk+1))**2 !n1*cantk dimension del tensor.
    read(12,*) aux1MAVe(1,i) !LEEMOS LOS TODOS LOS PUNTOS QUE FUERON GUARDADOS EN UNA FILA
 end do
 close (12)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!reasignacion de dimensiones!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(aux2MAVe(n1*(cantk+1),n1*(cantk+1)))

aux2MAVe=reshape(aux1MAVe,(/n1*(cantk+1),n1*(cantk+1)/))!REDIMENSIONAMOS TODOS LOS PUNTOS DE UNA FILA A UNA MATRIZ NUEVAMENTE
deallocate(aux1MAVe)

do i =1, n1*(cantk+1)
   do j = 1, n1*(cantk+1)
      !print*, aux2MAVe(i,j)
   end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(aux3MAVe(1,n1*(cantk+1)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print*, 'ingrese el autovalor deseado:'
read(*,*) nav !nav numero de autovector que se quiere observar(numero de columna de la matriz aux2MAVe)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i = 1, n1*(cantk+1)
   !print*, aux2MAVe(i,nav)
end do

aux3MAVe(1,:)=aux2MAVe(:,nav)!ASIGNAMOS EL VECTOR A UTILIZAR GRACIAS AL AUTOVECTOR DESEADO

deallocate(aux2MAVe)

do i = 1, n1*(cantk+1)
!   aux3MAVe(1,i)
end do

aux4MAVe=reshape(aux3MAVe,(/(cantk+1),n1/), order= (/2,1/) )!el orden importa en esta parde de la asignacion
                                                            ! de las constantes C.
deallocate(aux3MAVe)

!!!!!!!!!!!!! CALCULO Y LLENADO DE LA MATRIZ SOLUCION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!! Datos para la gráfica !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(unit=32,file='dat_splot.txt',STATUS="UNKNOWN",ACTION="WRITE")

     allocate(O(n1,nang))
     pi = 3.14159265358979323846
     ho = 2.0*pi/dble(nang)
     hr = 1.0/(dble(n1) + 1.0)
     
     do r = 1, n1

        do i = 1, nang
        
              ir=dble(i) - 1.0 
              c=0.0
              tetap= ho*ir
              do k = 2, (cantk+1) 
                 c = c + aux4MAVe(k,r)*COS(tetap*(k-1))   !r: radio
              end do
              O(r,i) = aux4MAVe(1,r) + 2*c
        end do

     end do   

     deallocate(aux4MAVe)

     allocate(Ob(n1+1,nang))  
  
     call inicializarmatriz(1, n1+1, nang, Ob)

     do i = 1, n1
        do j =1, nang
           Ob(i,j)=O(i,j)
        end do
     end do 


     do i=1, n1+1
        !do j = 1, nang

           ir=dble(dble(i)*hr )
           jr=dble(j) - 1.0 
           tetap= ho*jr
           !write(32,*) ir*cos(tetap), ir*sin(tetap), Ob(i,j)      
           write(32,*)  Ob(i,:)       
           !write(32,*) Ob(i,j)
           
           
         !end do
      end do   

   close(32)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



end program analisis
