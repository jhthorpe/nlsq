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
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: x,y
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: beta, beta0
    REAL(KIND=8), INTENT(IN) :: tol,hscal
    INTEGER, INTENT(IN) :: max_it,der_type
    INTEGER, INTENT(INOUT) :: stat
  
    !internal
    REAL(KIND=8) :: old_scs,new_csc,l
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: chisq,delta,fdiv 
    INTEGER, DIMENSION(:), ALLOCATABLE :: n
    INTEGER :: i,j,nb,m, iter

    stat = 2

    !get important numbers
    m = SIZE(x(:,1))
    nb = SIZE(beta(:))
    ALLOCATE(chisq(0:m-1))
    ALLOCATE(delta(0:nb-1))
    ALLOCATE(fdiv(0:nb-1))
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
    DO iter=0,max_it-1

      !3) find delta  
      !3a) get first partial deriviative vector
      CALL get_fdiv(fdiv,der_type,hscal,y,x,beta0,n,stat) 
      WRITE(*,*) "first partial derivative vector"
      WRITE(*,*) fdiv(:)
      IF (stat .NE. 2) STOP

      !3b) get second partial derivative vector
      !CALL get_sdiv(sdiv,der_type,hscal,y,x,beta0,n,stat)


    END DO
     
   stat = 0
 
  END SUBROUTINE opt_mrqt

!--------------------------------------------------------
! Get the first partial derivative vector for Marquardt optimization 
  SUBROUTINE get_fdiv(fdiv,der_type,hscal,y,x,beta0,n,stat)
    IMPLICIT NONE

    ! fdiv	:	1D sp first partial deriviatives vector
    ! beta0	:	1D sp vector of parameters 
    ! hscal	:	sp scaling for numerical derivatives
    ! der_type	:	int derivative type
    ! m		:	number of datasets
    ! nb	:	number of parameters
    ! n		:	number of values in each vector in x 
    ! rs	:	running sum

    !INOUT
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: fdiv
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: x,y
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: beta0
    REAL(KIND=8), INTENT(IN) :: hscal
    INTEGER, DIMENSION(0:), INTENT(IN) :: n
    INTEGER, INTENT(IN) :: der_type
    INTEGER, INTENT(INOUT) :: stat

    !Internal
    INTEGER :: i,j,k,nb,m
    REAL(KIND=8) :: rs
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: Jr

    stat = stat + 1

    nb = SIZE(beta0(:))
    m = SIZE(x(:,1))

    !if analytical derivatives
    IF (der_type .EQ. 0) THEN
      DO i=0,nb-1 !loop over parameters
        rs = 0.0E0
        DO j=0,m-1 !loop over residuals
          DO k=0,n(j)-1 !loop over indicies
            rs = rs + residual(j,k,y(j,:),x(j,:),beta0(:))*jacobian(j,i,k,y(j,:),x(j,:),beta0(:))/(1.0E0*n(j))
          END DO
        END DO
        fdiv(i) = rs
      END DO
  
    ELSE IF(der_type .EQ. 2) THEN
      !calculate jacobian     
      ALLOCATE(Jr(0:m-1,0:nb-1,0:n(0)-1)) !Jr (residual,parameter,index)
      CALL get_Jac(y(:,:),x(:,:),beta0(:),der_type,hscal,Jr(:,:,:))

      DO i=0,nb-1 !loop over parameters
        rs=0.0E0
        DO j=0,m-1 !loop over residuals
          DO k=0,n(j)-1 !loop over indicies
          rs = rs + residual(j,k,y(j,:),x(j,:),beta0(:))*Jr(j,i,k)/(1.0E0*n(j)) 
          END DO
        END DO
        fdiv(i) = rs
      END DO
       
    END IF
 
    stat = stat - 1

  END SUBROUTINE get_fdiv
!--------------------------------------------------------
! calculate the jacobian
  SUBROUTINE get_Jac(y,x,beta0,der_type,hscal,Jr)
    IMPLICIT NONE
    
    !INOUT
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: Jr
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: x,y
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: beta0
    REAL(KIND=8), INTENT(IN) :: hscal
    INTEGER, INTENT(IN) :: der_type

    !internal
    INTEGER :: i,j,k,nr,nb,n
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: betaf,betar

    nr = SIZE(x(:,0))     ! number of residuals
    nb = SIZE(beta0(:))   ! number of parameters
    n = SIZE(x(0,:))      ! number of datapoints

    ALLOCATE(betaf(0:nb-1))
    ALLOCATE(betar(0:nb-1))

    !central derivatives
    IF (der_type .EQ. 2) THEN
      DO i=0,nb-1 !loop over parameters
        betaf(:) = beta0(:)
        betar(:) = beta0(:)
        !change betaf and r
        betaf(i) = beta0(i) + hscal 
        betar(i) = beta0(i) - hscal
        DO j=0,nr-1 !loop over residuals
          DO k=0,n-1 !look over indicies
           Jr(j,i,k) = (residual(j,k,y(j,:),x(j,:),betaf(:)) - residual(j,k,y(j,:),x(j,:),betar(:)))/(-2.0E0*hscal) 
           !negative sign, because we are using the y(i)-f(i,b) residuals, which adds in a negative sign
          END DO
        END DO
      END DO

    ELSE
      WRITE(*,*) "Sorry, that jacobian has not be implimented yet"
      STOP
    END IF

  END SUBROUTINE get_Jac
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
