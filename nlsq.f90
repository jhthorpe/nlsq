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

  stat = 0

  WRITE(*,*) 
  WRITE(*,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" 
  WRITE(*,*) "Non-linear, least squares solver started"
  WRITE(*,*) 

  CALL get_input(x,y,beta0,tol,max_it,der_type,hscal,stat) 
  beta = beta0
  CALL opt_MRQT(beta,beta0,x,y,tol,max_it,der_type,hscal,stat)
  WRITE(*,*)   
  WRITE(*,*) "Completed with status: ", stat
  WRITE(*,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" 
  WRITE(*,*) 
  
END PROGRAM main
