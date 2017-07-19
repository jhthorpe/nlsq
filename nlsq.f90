!Program to perform global, nonlinear least squares fit using Newton-Gauss algorithm

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

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: x,y
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: beta0, beta
  DOUBLE PRECISION :: tol,foo
  INTEGER :: max_it,stat,i,j

  stat = 0

  WRITE(*,*) 
  WRITE(*,*) "%%%%%%%%%%%%%%%%%%%%" 
  WRITE(*,*) "Non-linear, least squares solver started"
  WRITE(*,*) 

  CALL get_input(x,y,beta0,tol,max_it,stat) 
  beta = beta0
  CALL opt_NG(beta,beta0,x,y,tol,max_it,stat)
  WRITE(*,*)   
  WRITE(*,*) "Completed with status: ", stat
  WRITE(*,*) "%%%%%%%%%%%%%%%%%%%%" 
  WRITE(*,*) 

END PROGRAM main
