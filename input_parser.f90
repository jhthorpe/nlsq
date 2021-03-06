MODULE nlsq_input

  IMPLICIT NONE

  CONTAINS

!-------------------------------------------------------
  SUBROUTINE get_input(x,y,ik,ec,nc,tol,conv,max_con,max_opt,der_type,hscal,c,cscal,lm,stat)
    IMPLICIT NONE
    ! x           :       2D DP array of x values to fit
    ! y           :       2D DP array of y values to fit
    ! ik          :       1D DP array of initial parameters
    ! ec          :       int number of equality constraints    
    ! nc          :       int number of inequality constraints
    ! tol         :       DP convergence tolerance 
    ! max_con     :       int max iterations for each optimization
    ! max_opt     :       int max number of lagrange iterations
    ! stat	  :	  int status
    
    !INOUT variables
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: x,y
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ik 
    REAL(KIND=8), INTENT(INOUT) :: tol,hscal,c,cscal,lm,conv
    INTEGER, INTENT(INOUT) :: max_opt,max_con,stat,der_type,ec,nc
    CHARACTER(LEN=32) :: fname,line

    !internal varaibles
    INTEGER :: nd,np,nk,i,j,io

    stat = 1 !reset at end if nothing goes wrong
 
    !defaults
    hscal = 1.0E-4

    OPEN (unit=1, file="input.dat", status="old", access="sequential")
   
    !initial parameters
    READ(1,*) nd    !number of datasets
    READ(1,*) nk    !number of independent parameters
    READ(1,*) ec
    READ(1,*) nc
    ALLOCATE(ik(0:nk-1))
    READ(1,*) 
    READ(1,*)
    READ(1,*) tol
    READ(1,*) conv
    READ(1,*) max_opt
    READ(1,*) line 
    IF (line .EQ. "analytic") THEN
      der_type = 0
    ELSE IF (line .EQ. "forward") THEN
      der_type = 1
    ELSE IF (line .EQ. "central") THEN
      der_type = 2
    ELSE
      WRITE(*,*) "Bad derivative type input, using central differences"
      der_type = 2
    END IF
    IF (der_type .EQ. 0) THEN
      READ(1,*)
    ELSE
      READ(1,*) hscal
    END IF 
    READ(1,*)
    READ(1,*)
    READ(1,*) max_con
    READ(1,*) c
    READ(1,*) cscal
    READ(1,*) lm
    READ(1,*)
    READ(1,*) 
    READ(1,*) ik
    READ(1,*)
    READ(1,*)

    !get data
    np = 0
    DO i=0,nd-1

      READ(1,*) fname 

      !get number of datapoints
      IF (i .EQ. 0) THEN
        OPEN(unit=2,file=fname,status="old",access="sequential")
        DO
          READ(2,*,iostat=io)
          IF (io .NE. 0) EXIT
          np = np + 1
        END DO
        CLOSE(unit=2,status="keep")

        ALLOCATE(x(0:nd-1,0:np-1))
        ALLOCATE(y(0:nd-1,0:np-1))
      END IF

      !read file
      OPEN(unit=2,file=fname,status="old",access="sequential")
      DO j=0,np-1
        READ(2,*) x(i,j), y(i,j)  
      END DO 
      CLOSE(unit=2,status="keep")
    END DO     

    CLOSE (unit=1, status="keep")
    
    CALL print_input(nd,nk,ik,ec,nc,tol,conv,der_type,hscal,c,cscal,lm,max_con,max_opt)

    stat = 0

  END SUBROUTINE get_input

!-------------------------------------------------------

  SUBROUTINE print_input(nd,nk,ik,ec,nc,tol,conv,der_type,hscal,c,cscal,lm,max_con,max_opt)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: ik 
    REAL(KIND=8), INTENT(IN) :: tol,hscal,c,cscal,lm,conv
    INTEGER, INTENT(IN) :: nd,nk,max_con,max_opt,der_type,ec,nc
    INTEGER :: i

    WRITE(*,*) "Number of Datasets : ", nd
    WRITE(*,*) "Number of residuals : ", nk
    WRITE(*,*) "Number of equality constraints : ", ec
    WRITE(*,*) "Number of inequality constraints : ", nc
    WRITE(*,*) "Initial parameters :"
    DO i=0,SIZE(ik)-1
      WRITE(*,*) ik(i)
    END DO 
    IF (der_type .EQ. 0) THEN
      WRITE(*,*) "Derivative type : analytical" 
      WRITE(*,*) "Remember to edit 'hard_code.f90'"
    ELSE IF (der_type .EQ. 1) THEN
      WRITE(*,*) "Derivative type : forward"
      WRITE(*,*) "h scale factor of : ", hscal
    ELSE
      WRITE(*,*) "Derivative type : central"
      WRITE(*,*) "h scale factor of : ", hscal
    END IF 
    IF (ec .NE. 0 .OR. nc .NE. 0) THEN
      WRITE(*,*) "Initial penalty constant : ", c
      WRITE(*,*) "Penalty constant scale factor : ", cscal
      WRITE(*,*) "Initial Lagrange Multiplier : ", lm
    END IF

  END SUBROUTINE print_input

END MODULE nlsq_input
