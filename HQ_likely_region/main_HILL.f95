! This file (main_QNSTOP.f95) contains a main program that calls
! QNSTOPP to minimize the KL value associated with the budding
! yeast cell cycle simulation.

PROGRAM MAIN_HILL
USE QNSTOPS_MOD
USE MODEL_HILL !, ONLY: HILL_OBJ  !module for mother simulation
IMPLICIT NONE

INTEGER:: NumInitValues, Q_EVALlim, Q_ITERlim, Q_STATUS, Q_SWITCH, Q_NSTART, Q_N
INTEGER:: Q_OMP, Q_TRACE, i
CHARACTER(LEN=1):: Q_MODE
REAL(KIND=R8):: Q_MIN_TAU, Q_GAIN, Q_TAU, Q_GAMMAV, Q_GAMMAW, Q_ETA, Q_FMIN


REAL(KIND=R8), DIMENSION(4):: Q_LB       ! Lower bounds.
REAL(KIND=R8), DIMENSION(4):: Q_UB       ! Upper bounds.
REAL(KIND=R8), DIMENSION(4):: Q_XI       ! Initial start point.
REAL(KIND=R8), DIMENSION(4):: hillParameters   !Parameters for hill optimization

 
! Name lists for problem configuration 'PROBLEM', optimization parameters
! 'OPTPARM'.
NAMELIST /PROBLEM/ NumInitValues, Q_LB, Q_UB, hillParameters, Q_MODE
NAMELIST /OPTPARM/ Q_ITERlim, Q_EVALlim, Q_SWITCH, Q_NSTART, Q_N,&
         Q_MIN_TAU, Q_GAIN, Q_TAU, Q_GAMMAV, Q_GAMMAW, Q_ETA, Q_OMP, Q_TRACE


OPEN(UNIT=111, FILE='parameter.nml', DELIM='APOSTROPHE')
READ(UNIT=111, NML=PROBLEM)
READ(UNIT=111, NML=OPTPARM)
CLOSE(111)

IF(Q_TRACE .NE. 0) Q_TRACE = 55

OPEN(UNIT=55, FILE = 'sOutput.txt', STATUS = 'replace', ACTION = 'write')
!WRITE(UNIT=55, NML=PROBLEM)
!WRITE(UNIT=55, NML=OPTPARM)

!Remove fixed and replicate values from parameter array 
DO i = 1,4
  Q_XI(i) = hillParameters(i)
END DO


CALL QNSTOPS(NumInitValues, Q_LB, Q_UB, HILL_OBJ, Q_XI, Q_MODE, Q_FMIN, Q_STATUS,&
  SWITCH=Q_SWITCH, NSTART=Q_NSTART, N=Q_N, MAX_ITER=Q_ITERlim, MAX_EVAL=Q_EVALlim, &
  MIN_TAU=Q_MIN_TAU, TAU=Q_TAU, GAIN=Q_GAIN, ETA=Q_ETA, TRACE=Q_TRACE)


IF ((Q_STATUS .GE. 10) .OR. (Q_FMIN > 400.0_R8)) THEN
  WRITE (*,113) Q_STATUS
   113 FORMAT('Stochastic optimization failed with status',I3.2,'.' )
  STOP
ELSE
   !Fill the array containing all the parameters for the hill model
   !using the values output by QNSTOP
   DO i = 1,4
     hillParameters(i) = Q_XI(i)
   END DO
   WRITE(*,115) Q_FMIN, hillParameters
   115 FORMAT('Stochastic optimization completed successfully with a &
      & minimum error of ',F15.8,'.',/ &
      'Optimized hill parameters ='/(8X,F8.4))
END IF

CLOSE(55)
END PROGRAM MAIN_HILL
