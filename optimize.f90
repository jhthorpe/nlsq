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
    ! nr	: 	int number of curves to fit
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
    REAL(KIND=8) :: old_scs,new_scs,l
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: chisq,new_chisq,delta,fdiv,track_scs 
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: sdiv, track_param 
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: Jr 
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: Hr 
    INTEGER, DIMENSION(:), ALLOCATABLE :: n
    INTEGER :: i,j,nb,nr,iter,flag

    flag = 0
    stat = 2

    !get important numbers
    nr = SIZE(x(:,0))
    nb = SIZE(beta(:))
    ALLOCATE(chisq(0:nr-1))
    ALLOCATE(new_chisq(0:nr-1))
    ALLOCATE(delta(0:nb-1))
    ALLOCATE(fdiv(0:nb-1))
    ALLOCATE(sdiv(0:nb-1,0:nb-1))
    ALLOCATE(track_scs(0:max_it-1))
    ALLOCATE(track_param(0:nb-1,0:max_it-1))
    ALLOCATE(n(0:nr-1))
    DO i=0,nr-1
      n(i) = SIZE(x(i,:))
    END DO
    ALLOCATE(Jr(0:nr-1,0:nb-1,0:n(0)-1)) !Jr (residual,parameter,index)
    ALLOCATE(Hr(0:nr-1,0:nb-1,0:nb-1,0:n(0)-1)) !Hr (residual, parameter, parameter, index)

    WRITE(*,*) 
    WRITE(*,*) "~~~~~~~~~~~~~~~~~~~~"
    WRITE(*,*) "Starting Marquardt optimization"
    WRITE(*,*) 

    !Initial Setup 

    ! 1) calcuate initial sum chi^2
    old_scs = 0.0D0
    DO i=0,nr-1
      chisq(i) = 0.0D0
      DO j=0,n(i)-1
        chisq(i) = chisq(i) + (residual(i,j,y(i,:),x(i,:),beta0(:))**2.0D0) 
      END DO
      !chisq(i) = chisq(i)/(1.0D0*n(i))
      old_scs = old_scs + chisq(i)
    END DO

    WRITE(*,*) "Initial sum of squares vector :"
    WRITE(*,*) chisq(:) 
    WRITE(*,*) "Initial sum of sum of squares is...", old_scs
    WRITE(*,*) 

    !2) set lambda
    l = 10.0D-03

    !3) Perform algorithm
    DO iter=0,max_it-1


      IF (iter .GT. 0) THEN
        WRITE(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        WRITE(*,*) "Iteration : ", iter
        WRITE(*,*)
        WRITE(*,*) "Starting parameters" 
        WRITE(*,*) beta0(:) 
        WRITE(*,*) "Lambda : ", l
        WRITE(*,*) 
      END IF


      !3) find delta  
      !3a) get first partial deriviative vector
      CALL get_fdiv(fdiv,der_type,hscal,y,x,beta0,n,Jr,stat) 
      IF (stat .NE. 2) STOP
      WRITE(*,*) "first partial derivative vector"
      WRITE(*,*) fdiv(:)
      WRITE(*,*)

      !3b) get second partial derivative vector
      CALL get_sdiv(sdiv,der_type,hscal,y,x,beta0,n,Jr,Hr,stat)
      IF (stat .NE. 2) STOP
      WRITE(*,*) "second partial derivative matrix"
      DO j=0,nb-1
        WRITE(*,*) sdiv(:,j)
      END DO
      WRITE(*,*)

      !3c) multiply by identiy*(1+lambda)
      DO i=0,nb-1
        sdiv(i,i) = sdiv(i,i) * (1.0D0 + l)
      END DO
      WRITE(*,*) "matrix alpha"
      DO j=0,nb-1
        WRITE(*,*) sdiv(:,j) 
      END DO

      !3d) Invert alpha matrix
      CALL invert_2ddp(sdiv,stat)
      IF (stat .NE. 2) STOP
      WRITE(*,*)
      WRITE(*,*) "inverted alpha"
      DO j=0,nb-1
        WRITE(*,*) sdiv(:,j)
      END DO

      !3e) Get delta
      DO i=0,nb-1
        delta(i) = 0.0D0
        DO j=0,nb-1
          delta(i) = delta(i) + fdiv(j)*sdiv(j,i)
        END DO
      END DO

      !4) get new alpha
      DO i=0,nb-1
        beta(i) = beta0(i) + delta(i)
      END DO
      WRITE(*,*) 
      WRITE(*,*) "trial parameters"
      WRITE(*,*) beta(:)
     
     !5) get trail fit
      new_scs = 0.0D0
      DO i=0,nr-1
        new_chisq(i) = 0.0D0
        DO j=0,n(i)-1
          new_chisq(i) = new_chisq(i) + (residual(i,j,y(i,:),x(i,:),beta(:))**2.0D0) 
        END DO
        !new_chisq(i) = new_chisq(i)/(1.0D0*n(i))
        new_scs = new_scs + new_chisq(i)
      END DO
      WRITE(*,*)
      WRITE(*,*) "Trial sum of squares vector :"
      WRITE(*,*) new_chisq(:) 
      WRITE(*,*) "Trial sum of sum of squares is...", new_scs
      WRITE(*,*) 
  
      !6) if sum of sum of squares is worse...
      IF (new_scs .GE. old_scs) THEN
        WRITE(*,*) "No improvement, adjusting lambda" 
        CALL summary(beta0,beta0,old_scs,old_scs,chisq,chisq,iter)
        l = l * 10.0D0
        WRITE(*,*) "New lambda : ", l 
      ELSE
        IF (new_scs .LE. tol) flag = 1 !we have reached the cutoff we care about
        WRITE(*,*) "Adopting new parameters"
        CALL summary(beta0,beta,old_scs,new_scs,chisq,new_chisq,iter)
        l = l * 0.10D0
        beta0 = beta(:)
        old_scs = new_scs
        chisq = new_chisq(:)
        WRITE(*,*) "New lambda : ", l
      END IF

      !track output
      track_scs(iter) = old_scs
      track_param(:,iter) = beta0(:) 

      IF (flag .EQ. 1 .OR. flag .EQ. 2) EXIT

         
    END DO
     
    IF (flag .EQ. 0) THEN
      WRITE(*,*) "Max iterations reached"
    ELSE IF (flag .EQ. 1) THEN
      WRITE(*,*) "Sum of Sum of Square Difference below cutoff" 
    ELSE IF (flag .EQ. 2) THEN
      WRITE(*,*) "Optimization has converged"
    END IF

    CALL write_output(track_scs,track_param,iter)

   stat = 0
 
  END SUBROUTINE opt_mrqt

!--------------------------------------------------------
! print output of optimization
  SUBROUTINE write_output(goal,param,iter)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: goal 
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: param
    INTEGER, INTENT(IN) :: iter

    INTEGER :: i

    WRITE(*,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"  
    WRITE(*,*)  
    WRITE(*,*) "                          SUMMARY OF CALCULATION" 
    WRITE(*,*) "iteration      X^2     parameters" 
    DO i=0,iter-1
      WRITE(*,*) i, goal(i), param(:,i)
    END DO
    WRITE(*,*) 
        

  END SUBROUTINE write_output
!--------------------------------------------------------
!  print sumary of iteration
  SUBROUTINE summary(old_b,new_b,old_scs,new_scs,old_cs,new_cs,iter)
    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: old_scs,new_scs
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: old_cs,new_cs,old_b,new_b
    INTEGER, INTENT(IN) :: iter

    INTEGER :: nb,nr,i,j
    
    nr = SIZE(new_cs(:))
    nb = SIZE(new_b(:))
      WRITE(*,*)
      WRITE(*,*) "================================================================================"
      WRITE(*,*) "  SUMMARY OF ITERATION", iter 
      WRITE(*,*) "--------------------------------------------------------------------------------" 
      WRITE(*,*) "  PARAMETERS" 
      WRITE(*,*) "  old      new      change      %  "
      DO i=0,nb-1
        WRITE(*,*) old_b(i), new_b(i), new_b(i)-old_b(i), 100.0D0*(new_b(i)-old_b(i))/old_b(i)
      END DO
      WRITE(*,*) "-----------"
      WRITE(*,*) "SUM OF SQUARE DIFFERENCES"  
      WRITE(*,*) "  old      new      change      %  "
      DO i=0,nr-1
        WRITE(*,*) old_cs(i), new_cs(i), new_cs(i)-old_cs(i), 100.0E0*(new_cs(i)-old_cs(i))/old_cs(i)
      END DO
      WRITE(*,*) "-----------"
      WRITE(*,*) "SUM OF SUM OF SQUARE DIFFERENCES"
      WRITE(*,*) "  old      new      change      %  "
      WRITE(*,*) old_scs, new_scs, old_scs-new_scs, 100.0D0*(new_scs-old_scs)/old_scs
      WRITE(*,*) "================================================================================"
      WRITE(*,*)
   END SUBROUTINE summary
!--------------------------------------------------------
! Get the first partial derivative vector for Marquardt optimization 
  SUBROUTINE get_fdiv(fdiv,der_type,hscal,y,x,beta0,n,Jr,stat)
    IMPLICIT NONE

    ! fdiv	:	1D sp first partial deriviatives vector
    ! beta0	:	1D sp vector of parameters 
    ! hscal	:	sp scaling for numerical derivatives
    ! der_type	:	int derivative type
    ! Jr	:	3D sp jacobian 
    ! nr	:	number of datasets
    ! nb	:	number of parameters
    ! n		:	number of values in each vector in x 
    ! rs	:	running sum

    !INOUT
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: fdiv
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: x,y
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: Jr
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: beta0
    REAL(KIND=8), INTENT(IN) :: hscal
    INTEGER, DIMENSION(0:), INTENT(IN) :: n
    INTEGER, INTENT(IN) :: der_type
    INTEGER, INTENT(INOUT) :: stat

    !Internal
    INTEGER :: i,j,k,nb,nr
    REAL(KIND=8) :: rs

    stat = stat + 1

    nb = SIZE(beta0(:))
    nr = SIZE(x(:,0))

    !if analytical derivatives
    IF (der_type .EQ. 0) THEN
      DO i=0,nb-1 !loop over parameters
        rs = 0.0E0
        DO j=0,nr-1 !loop over residuals
          DO k=0,n(j)-1 !loop over indicies
            rs = rs + residual(j,k,y(j,:),x(j,:),beta0(:))*jacobian(j,i,k,y(j,:),x(j,:),beta0(:))!/(1.0E0*n(j))
          END DO
        END DO
        fdiv(i) = rs
      END DO
  
    ELSE IF(der_type .EQ. 2) THEN
      !calcualat the jacobian
      CALL get_Jac(y(:,:),x(:,:),beta0(:),der_type,hscal,Jr(:,:,:))
      DO i=0,nb-1 !loop over parameters
        rs=0.0E0
        DO j=0,nr-1 !loop over residuals
          DO k=0,n(j)-1 !loop over indicies
          rs = rs + residual(j,k,y(j,:),x(j,:),beta0(:))*Jr(j,i,k)!/(1.0E0*n(j)) 
          END DO
        END DO
        fdiv(i) = rs
      END DO
       
    END IF
 
    stat = stat - 1

  END SUBROUTINE get_fdiv
!--------------------------------------------------------
  SUBROUTINE get_sdiv(sdiv,der_type,hscal,y,x,beta0,n,Jr,Hr,stat)
    IMPLICIT NONE

    ! sdiv	:	2D sp second partial deriviatives vector
    ! Jr	:	3D sp jacobian
    ! Hr	:	4D sp hessian
    ! beta0	:	1D sp vector of parameters 
    ! hscal	:	sp scaling for numerical derivatives
    ! der_type	:	int derivative type
    ! nr	:	number of datasets
    ! nb	:	number of parameters
    ! n		:	number of values in each vector in x 
    ! rs	:	running sum

    !INOUT
    REAL(KIND=8), DIMENSION(0:,0:,0:,0:), INTENT(INOUT) :: Hr
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: Jr
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: sdiv
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: x,y
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: beta0
    REAL(KIND=8), INTENT(IN) :: hscal
    INTEGER, DIMENSION(0:), INTENT(IN) :: n
    INTEGER, INTENT(IN) :: der_type
    INTEGER, INTENT(INOUT) :: stat

    !Internal
    INTEGER :: i,j,k,l,nb,nr,ni
    REAL(KIND=8) :: rs

    stat = stat + 1

    nb = SIZE(beta0(:))
    nr = SIZE(x(:,0))
    ni = SIZE(x(0,:))

    DO i=0,nb-1      !loop over parameter i
      DO j=0,nb-1    !loop over parameter j
        rs = 0.0E0
        DO k=0,nr-1  !loop over residual 
          DO l=0,ni-1 ! loop over index
            rs = rs + Jr(k,i,l)*Jr(k,j,l)!/ni !this is a strange way to get the second derivative? 
          END DO
        END DO 
        sdiv(i,j) = rs
      END DO
    END DO
    

    stat = stat - 1

  END SUBROUTINE get_sdiv
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
! use lapack and blas to invert a 2d, dp matrix
  SUBROUTINE invert_2ddp(A,stat)
    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(0:,0:),INTENT(INOUT) :: A
    INTEGER, INTENT(INOUT) :: stat

    !LAPACK specific values
    REAL(KIND=8), DIMENSION(SIZE(A(:,0))) :: work !lapack working array
    INTEGER(KIND=8), DIMENSION(SIZE(A(:,0))) :: ipiv !lapack index storage
    INTEGER :: n,info

    n = SIZE(A(:,0))
    stat = stat + 1

    CALL DGETRF(n,n,A,n,ipiv,info)
    
    IF (info .NE. 0) THEN
      WRITE(*,*) "Not a square matrix"
      STOP
    END IF

    CALL DGETRI(n,A,n,ipiv,work,n,info)
 
    IF (info .NE. 0) THEN
      WRITE(*,*) "matrix inversion failed"
      STOP
    END IF

    stat = stat - 1

  END SUBROUTINE invert_2ddp 
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
