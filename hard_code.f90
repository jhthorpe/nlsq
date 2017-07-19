!module for hard coding residuals and jacobians, input your functions here

MODULE hard_code

  IMPLICIT NONE
 
  CONTAINS 

!--------------------------------------------------------
!individual residuals
   SUBROUTINE residual(i,y,x,beta)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    DOUBLE PRECISION, DIMENSION(0:), INTENT(IN) :: beta,x,y
     

  END SUBROUTINE residual
!--------------------------------------------------------
!analytical jacobians indicies
  DOUBLE PRECISION FUNCTION jacobian(i,j,beta)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i,j
    DOUBLE PRECISION, DIMENSION(0:), INTENT(IN) :: beta

    WRITE(*,*) i,j
    WRITE(*,*) beta(i)

  END FUNCTION jacobian
!--------------------------------------------------------

END MODULE hard_code
