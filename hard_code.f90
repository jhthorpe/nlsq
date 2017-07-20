!module for hard coding residuals and jacobians, input your functions here

MODULE hard_code

  IMPLICIT NONE
 
  CONTAINS 

!--------------------------------------------------------
!individual residuals
 REAL FUNCTION residual(p,i,y,x,beta)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i,p
    REAL, DIMENSION(0:), INTENT(IN) :: beta
    REAL, DIMENSION(0:), INTENT(IN) :: x,y

    !Construct your residuals here 
    IF (p .EQ. 0) THEN
      residual = y(i) - (EXP(-(beta(0)+beta(1))*x(i)))  
    ELSE IF (p .EQ. 1) THEN
      residual = y(i) - (beta(0)/(beta(0)+beta(1)))*(1.0D0 - EXP(-(beta(0)+beta(1))*x(i))) 
    ELSE IF (p .EQ. 2) THEN
      residual = y(i) - (beta(1)/(beta(0)+beta(1)))*(1.0D0 - EXP(-(beta(0)+beta(1))*x(i))) 
    ELSE
      WRITE(*,*) "WARNING - attempted to reference non-defined residual"
      STOP
    END IF


  END FUNCTION residual
!--------------------------------------------------------
!analytical jacobians indicies
  REAL FUNCTION jacobian(i,j,beta)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i,j
    REAL, DIMENSION(0:), INTENT(IN) :: beta

    WRITE(*,*) i,j
    WRITE(*,*) beta(i)

  END FUNCTION jacobian
!--------------------------------------------------------

END MODULE hard_code
