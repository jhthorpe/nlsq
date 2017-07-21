!Program to fit multiple datasets with shared parameters

PROGRAM main

  USE nlsq_input     !gets the input
  USE optimize   !runs the optimization

  IMPLICIT NONE

  ! x		:	2D DP array of x values to fit
  ! y		:	2D DP array of y values to fit
  ! ik		:	1D DP array of initial parameters
  ! k		:	1D DP array of new parameters
  ! tol 	:	DP convergence tolerance 
  ! max_it	:	int max iterations
  ! stat	:	int status
  ! der_type	:	char - type of derivatives chosen, 0 - analytic, 1 - forward, 2 - central
  ! hscal	:	sp scale factor for h with numerical derivatives

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: x,y
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: beta0, beta
  REAL(KIND=8) :: tol,hscal
  INTEGER :: max_it,stat,i,j,der_type

  !testing lapack
  real(KIND=8), dimension(0:1,0:1) :: A,Ainv
  real(KIND=8), dimension(0:1) :: work
  integer(KIND=8), dimension(0:1) :: ipiv
  integer :: q,info

  stat = 0

  WRITE(*,*) 
  WRITE(*,*) "%%%%%%%%%%%%%%%%%%%%" 
  WRITE(*,*) "Non-linear, least squares solver started"
  WRITE(*,*) 

  CALL get_input(x,y,beta0,tol,max_it,der_type,hscal,stat) 
  beta = beta0
  CALL opt_MRQT(beta,beta0,x,y,tol,max_it,der_type,hscal,stat)
  WRITE(*,*)   
  WRITE(*,*) "Completed with status: ", stat
  WRITE(*,*) "%%%%%%%%%%%%%%%%%%%%" 
  WRITE(*,*) 

  !testing lapack
  WRITE(*,*) "Testing lapack"
  !external DGETRF
  !external DGETRI
  A(0,0) = 1.0D0 
  A(0,1) = 1.0D0
  A(1,0) = 0.0D0 
  A(1,1) = 2.0D0
  WRITE(*,*) A(:,0)
  WRITE(*,*) A(:,1)
  Ainv = A
  q = 2
  info = 0
  call DGETRF(q,q,Ainv,q,ipiv,info)

  if (info .NE. 0) then
    WRITE(*,*) "matrix is singular"
    STOP
  else 
    write(*,*)"here"
  end if

  WRITE(*,*) ipiv(:)

  call DGETRI(q,Ainv,q,ipiv,work,q,info)

  write(*,*) "here1"

  if (info .NE. 0) THEN
    write(*,*) "matrix inversion failed"
    stop
  else
    write(*,*) Ainv(0,0),Ainv(1,0)
    write(*,*) Ainv(0,1),Ainv(1,1)
  end if

  
END PROGRAM main
