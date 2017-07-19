MODULE nlsq_input

  IMPLICIT NONE

  CONTAINS

!-------------------------------------------------------
  SUBROUTINE get_input(x,y,ik,tol,max_it,stat)
    IMPLICIT NONE
    ! x           :       2D DP array of x values to fit
    ! y           :       2D DP array of y values to fit
    ! ik          :       1D DP array of initial parameters
    ! tol         :       DP convergence tolerance 
    ! max_it      :       int max iterations
    ! stat	  :	  int status
    
    !INOUT variables
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: x,y
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ik 
    DOUBLE PRECISION, INTENT(INOUT) :: tol
    INTEGER, INTENT(INOUT) :: max_it,stat
    CHARACTER(LEN=32) :: fname

    !internal varaibles
    INTEGER :: nd,np,nk,i,j,io

    stat = 1 !reset at end if nothing goes wrong

    OPEN (unit=1, file="input.dat", status="old", access="sequential")
   
    !initial parameters
    READ(1,*) nd    !number of datasets
    READ(1,*) nk    !number of independent parameters
    ALLOCATE(ik(0:nk-1))
    READ(1,*) tol
    READ(1,*) max_it
    READ(1,*)
    READ(1,*) 
    READ(1,*) ik
    READ(1,*)
    READ(1,*)


    !get number of datapoints
    

    !get data
    DO i=0,nd-1

      READ(1,*) fname 
      np = 0

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
        READ(2,*) y(0,j), x(0,j)  
      END DO 
      CLOSE(unit=2,status="keep")
    END DO     

    CLOSE (unit=1, status="keep")
    
    CALL print_input(nd,nk,ik,tol,max_it)

    stat = 0

  END SUBROUTINE get_input

!-------------------------------------------------------

  SUBROUTINE print_input(nd,nk,ik,tol,max_it)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(0:), INTENT(IN) :: ik 
    DOUBLE PRECISION, INTENT(IN) :: tol
    INTEGER, INTENT(IN) :: nd,nk,max_it
    INTEGER :: i

    WRITE(*,*) "Number of Datasets : ", nd
    WRITE(*,*) "Number of independent parameters : ", nk
    WRITE(*,*) "Initial parameters :"
    DO i=0,SIZE(ik)-1
      WRITE(*,*) ik(i)
    END DO 
    

  END SUBROUTINE print_input

END MODULE nlsq_input
