program test

  implicit none

  integer :: i,j

  do i=0, 10
    write(*,*) i
  end do

  write(*,*) "i is : ",i
  write(*,*)
   
  do j=0,10
    write(*,*) j
    if (j .GE. 10) EXIT
  end do 
  
  write(*,*) "j is : ", j

end program test
