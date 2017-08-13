!Program to fit multiple datasets with shared parameters

PROGRAM main

  USE nlsq_input     !gets the input
  USE optimize   !runs the optimization

  IMPLICIT NONE

  ! x		:	2D DP array of x values to fit
  ! y		:	2D DP array of y values to fit
  ! ik		:	1D DP array of initial parameters
  ! ec		:       int number of equality constraints	
  ! nc		:	int number of inequality constraints
  ! k		:	1D DP array of new parameters
  ! tol 	:	DP convergence tolerance 
  ! max_com	:	int max number of lagrange iterations
  ! max_opt	:	int max number of optimizater iterations
  ! c		:	DP penalty constant
  ! cscal	:	DP scale factor for penalty constant
  ! lm		:	initial Lagrange multiplier for constraints
  ! stat	:	int status
  ! der_type	:	char - type of derivatives chosen, 0 - analytic, 1 - forward, 2 - central
  ! hscal	:	sp scale factor for h with numerical derivatives

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: x,y
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: beta0, beta
  REAL(KIND=8) :: tol,conv,hscal,c,cscal,lm
  INTEGER :: max_con,max_opt,stat,i,j,der_type,ec,nc

  stat = 0

  WRITE(*,*) 
  WRITE(*,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" 
  WRITE(*,*) "Non-linear, least squares solver started"
  WRITE(*,*) 

  CALL get_input(x,y,beta0,ec,nc,tol,conv,max_con,max_opt,der_type,hscal,c,cscal,lm,stat) 
  beta = beta0

  !check input
  
  IF (SIZE(x(:,0)) .LT. SIZE(beta0(:))) WRITE(*,*) "WARNING, your system is overdetermined." 

  IF (nc .NE. 0) THEN
    WRITE(*,*) "Sorry, inequality constraints not yet implimented"
    STOP
  ELSE IF (ec .EQ. 0 .AND. nc .EQ. 0) THEN
    CALL opt_LM(beta,beta0,x,y,tol,max_opt,der_type,hscal,stat)
  ELSE IF (ec .NE. 0 .AND. nc .EQ. 0) THEN
    CALL opt_eq(beta,beta0,ec,nc,x,y,tol,conv,max_con,max_opt,der_type,hscal,c,cscal,lm,stat) 
  END IF
  
  WRITE(*,*)   
  WRITE(*,*) "Completed with status: ", stat
  WRITE(*,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" 
  WRITE(*,*) 
  
END PROGRAM main
