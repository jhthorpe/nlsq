MODULE optimize
  USE hard_code
  IMPLICIT NONE

  CONTAINS

!--------------------------------------------------------
! optimization runner

  SUBROUTINE opt_NG(beta,beta0,x,y,tol,max_it,stat)
    IMPLICIT NONE
    ! r         :       1D dp residuals (yi - f(xi,beta))
    ! Jf        :       2D dp function jacobians 
    ! Jft	:	2D dp transpose of Jf
    ! x         :       2D dp x values
    ! y         :       2D dp y values
    ! beta      :       1D dp parameters
    ! tol	:	dp convergence tolerance
    ! max_it	:	int maximum iterations
    ! stat	:	integer status

    !INOUT
    DOUBLE PRECISION, DIMENSION(0:,0:), INTENT(IN) :: x,y
    DOUBLE PRECISION, DIMENSION(0:), INTENT(INOUT) :: beta, beta0
    DOUBLE PRECISION, INTENT(IN) :: tol
    INTEGER, INTENT(INOUT) :: max_it,stat
    
    !internal
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: r
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Jf, Jft
    INTEGER :: np,nr,i,j

    stat = 2
    np = SIZE(beta)
    nr = SIZE(x(:,0))

    WRITE(*,*) "Starting Newton-Gauss optimization"
    WRITE(*,*)
  
    !Initialize values    
    ALLOCATE(r(0:nr-1))
    ALLOCATE(Jf(0:np-1,0:nr-1)) !check these dimensions
    stat = 0
    
  END SUBROUTINE opt_NG


!--------------------------------------------------------
  !update residuals and jacobians
  SUBROUTINE update_values(r,Jf,x,y,beta)
    IMPLICIT NONE
    ! r         :       1D dp residuals (yi - f(xi,beta))
    ! Jf        :       2D dp function jacobians 
    ! x         :       2D dp x values
    ! y         :       2D dp y values
    ! beta      :       1D dp parameters

    DOUBLE PRECISION, DIMENSION(0:,0:), INTENT(INOUT) :: Jf
    DOUBLE PRECISION, DIMENSION(0:), INTENT(INOUT) :: r,x,y,beta



  END SUBROUTINE update_values
!--------------------------------------------------------
 

END MODULE optimize
