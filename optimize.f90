MODULE optimize
  USE hard_code
  IMPLICIT NONE

  CONTAINS

!--------------------------------------------------------
! Marquardt optimization, from Marquardt-1963
  SUBROUTINE opt_mrqt(beta,beta0,x,y,tol,max_it,der_type,hscal,stat)
    IMPLICIT NONE
    ! r         :       1D sp residuals (yi - f(xi,beta))
    ! Jf        :       2D sp function jacobians 
    ! Jft	:	2D sp transpose of Jf
    ! x         :       2D sp x values
    ! y         :       2D sp y values
    ! beta      :       1D sp parameters
    ! tol	:	sp convergence tolerance
    ! max_it	:	int maximum iterations
    ! stat	:	integer status
    ! m		: 	int number of curves to fit
    ! nb	:	int number of parameters
    ! n		:	1D int number of datapoints
    ! l		: 	sp lambda, controls algorithm
    ! der_type	:	int, type of derivatives chosen: 0- analytic, 1-forward, 2-central
    ! hscal	:	sp, scale factor for numerical derivatives

    !INOUT    
    REAL, DIMENSION(0:,0:), INTENT(IN) :: x,y
    REAL, DIMENSION(0:), INTENT(INOUT) :: beta, beta0
    REAL, INTENT(IN) :: tol,hscal
    INTEGER, INTENT(IN) :: max_it,der_type
    INTEGER, INTENT(INOUT) :: stat
  
    !internal
    REAL :: old_scs,new_csc,l
    REAL, ALLOCATABLE, DIMENSION(:) :: chisq,delta 
    INTEGER, DIMENSION(:), ALLOCATABLE :: n
    INTEGER :: i,j,nb,m, iter

    stat = 2

    !get important numbers
    m = SIZE(x(:,1))
    nb = SIZE(beta(:))
    ALLOCATE(chisq(0:m-1))
    ALLOCATE(delta(0:nb-1))
    ALLOCATE(n(0:m-1))
    DO i=0,m-1
      n(i) = SIZE(x(i,:))
    END DO

    WRITE(*,*) 
    WRITE(*,*) "~~~~~~~~~~~~~~~~~~~~"
    WRITE(*,*) "Starting Marquardt optimization"
    WRITE(*,*) 

    !Initial Setup 

    ! 1) calcuate initial sum chi^2
    old_scs = 0.0E0
    DO i=0,m-1
      chisq(i) = 0.0E0
      DO j=0,n(i)-1
        chisq(i) = chisq(i) + (residual(i,j,y(i,:),x(i,:),beta0(:))**2.0E0) 
      END DO
      chisq(i) = chisq(i)/(1.0E0*n(i))
      old_scs = old_scs + chisq(i)
    END DO

    WRITE(*,*) "Initial sum of squares vector :"
    WRITE(*,*) chisq(:) 
    WRITE(*,*) "Initial sum of sum of squares is...", old_scs
    WRITE(*,*) 

    !2) set lambda
    l = 10.0E-03

    !3) Perform algorithm
    DO iter=0,max_it

      !3) find delta  


    END DO
     
   stat = 0
 
  END SUBROUTINE opt_mrqt

!--------------------------------------------------------
! Gauss-Newton Least Squares optimization 

  SUBROUTINE opt_GN(beta,beta0,x,y,tol,max_it,stat)
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
    INTEGER :: nb,nr,i,j

    stat = 2
    nb = SIZE(beta)
    nr = SIZE(x(:,0))

    WRITE(*,*) "Starting Newton-Gauss optimization"
    WRITE(*,*)
  
    !Initialize values    
    ALLOCATE(r(0:nr-1))
    ALLOCATE(Jf(0:nb-1,0:nr-1)) !check these dimensions
    stat = 0
    
  END SUBROUTINE opt_GN


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
