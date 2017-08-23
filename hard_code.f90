!module for hard coding residuals and jacobians, input your functions here

MODULE hard_code

  IMPLICIT NONE
 
  CONTAINS 

!--------------------------------------------------------
!individual residuals (required)
 REAL(KIND=8) FUNCTION residual(r,i,y,x,beta)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i,r
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: beta
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: x,y

    ! r		:	integer control value, indicates which residual you want to calculate
    ! i		:	integer, idicates which value of the y vector you want to get
    ! beta	:	double prescision, 1D vector, contains parameters
    ! x,y	:	double precision, 1D vectors, contain x and y values

    !Construct your residuals here 
    IF (r .EQ. 0) THEN
      residual = y(i) - (EXP(-(beta(0)+beta(1))*x(i)))  
    ELSE IF (r .EQ. 1) THEN
      residual = y(i) - (beta(0)/(beta(0)+beta(1)))*(1.0D0 - EXP(-(beta(0)+beta(1))*x(i))) 
    ELSE IF (r .EQ. 2) THEN
      residual = y(i) - (beta(1)/(beta(0)+beta(1)))*(1.0D0 - EXP(-(beta(0)+beta(1))*x(i))) 
    ELSE
      WRITE(*,*) "WARNING - attempted to reference non-defined residual"
      STOP
    END IF


  END FUNCTION residual
!--------------------------------------------------------
!sum of individual constraints at a given time (required if constrained optimization)
  REAL(KIND=8) FUNCTION eq_con(i,y,x,beta)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: i
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: beta
    REAL(KIND=8),DIMENSION(0:,0:), INTENT(IN) :: x,y

    REAL (KIND=8) :: temp
    
    ! n		:	int selects which constraint to use
    ! i		:	int selects y or x values
    ! beta	:	1D DP list of parameters
    ! x,y	:	2D dp list of x,y values - just in case you need some strange constraints
    ! temp	:	running sum
    
    temp = 0.0D0

    !construct and sum up your constraints here, must be of form h(x) = 0
    
    !EXAMPLES
    !two parameters are equal : temp = temp + (beta(0) - beta(1))
    !one parameter has constant value constant: temp = tmep + (beta(1) - 0.33D0)

    temp = temp + (beta(1) - 0.25D0)

    eq_con = temp 

  END FUNCTION eq_con 
!--------------------------------------------------------
! hard coded jaccobian (optional)
  REAL(KIND=8) FUNCTION jacobian(r,p,i,y,x,beta)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i,p,r
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: beta
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: x,y
    

    ! r		:	integer control value, indicates which residual you want to use
    ! p		:	integer control value, indicates which parameter you want to use
    ! i		:	integer, idicates which index of the y vector you want to get
    ! beta	:	single prescision, 1D vector, contains parameters
    ! x,y	:	single precision, 1D vectors, contain x and y values

    IF (r .EQ. 0) THEN
      jacobian = -1.0E0*x(i)*EXP(-(beta(0)+beta(1))*x(i))
    ELSE IF (r .EQ. 1) THEN
      IF (p .EQ. 0) THEN
        jacobian = EXP(-(beta(0)+beta(1))*x(i)) &
                  *(beta(0)**2.0E0*x(i) + beta(1)*(EXP((beta(0)+beta(1))*x(i)) + beta(0)*x(i) - 1.0E0)) &
                  /(beta(0)+beta(1))**2.0E0
      ELSE
        jacobian = (beta(0)/(beta(0)+beta(1))**2.0E0)*EXP(-1.0E0*(beta(0)+beta(1))*x(i))&
                  * (-1.0E0*EXP((beta(0)+beta(1))*x(i)) + beta(0)*x(i)+beta(1)*x(i) + 1.0E0)
      END IF
    ELSE IF (r .EQ. 2) THEN
      IF (p .EQ. 0) THEN
        jacobian = (beta(1)/(beta(0)+beta(1))**2.0E0)*EXP(-1.0E0*(beta(0)+beta(1))*x(i))&
                  * (-1.0E0*EXP((beta(0)+beta(1))*x(i)) + beta(0)*x(i)+beta(1)*x(i) + 1.0E0)
      ELSE
        jacobian = EXP(-(beta(0)+beta(1))*x(i)) &
                  *(beta(1)**2.0E0*x(i) + beta(0)*(EXP((beta(0)+beta(1))*x(i)) + beta(1)*x(i) - 1.0E0)) &
                  /(beta(0)+beta(1))**2.0E0
      END IF
    ELSE
      WRITE(*,*) "WARNING - attempted to reference non-defined jacobian"
      STOP
    END IF
    
  END FUNCTION jacobian

END MODULE hard_code
