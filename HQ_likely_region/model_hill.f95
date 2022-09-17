MODULE MODEL_HILL
USE REAL_PRECISION, ONLY: R8
!USE DUNI_OMP
!$ USE OMP_LIB, ONLY: OMP_GET_THREAD_NUM

IMPLICIT NONE

! PRIVATE




PUBLIC:: HILL_OBJ


CONTAINS


FUNCTION HILL_OBJ(PARAMS, IFLAG) RESULT(Error)
! This is an objective function that fits experimental statistics from
! Table 1 of [1] from the budding yeast cell cycle model.

  REAL(KIND=R8), INTENT(IN):: PARAMS(:) ! Must be 4-dimensional,
! but is assumed shape for optimization algorithms.
  INTEGER, INTENT(OUT):: IFLAG
  REAL(KIND=R8):: Error ! Fitting error from experimental statistics.

! LIMITED DATA.
  INTEGER, PARAMETER:: NUMD = 50
  REAL(KIND=R8), PARAMETER:: t = 0.2_R8


! INTEGER, PARAMETER:: DataY(NUMD) = (/ &
!   15, 42, 60, 57, 49, 66, 70, 66, 60, 60 &
! /) ! e 10


! INTEGER, PARAMETER:: DataY(NUMD) = (/ &
!  0,   1,   7,  11,  15,  19,  29,  33,  38,  42, &
! 52,  54,  56,  58,  60,  71,  66,  61,  55,  57, &
! 54,  57,  60,  55,  49,  61,  63,  69,  74,  66, &
! 63,  60,  63,  67,  70,  68,  66,  68,  69,  66, &
! 61,  62,  61,  63,  60,  63,  65,  63,  61,  60 &
! /) ! e 50

! INTEGER, PARAMETER:: DataY(NUMD) = (/ &
! 0,   0,   0,   0,   0,   1,   1,   1,   1,   4, &
! 5,   7,   7,   6,   8,  11,  12,  14,  14,  15, &
! 15,  16,  17,  19,  21,  22,  23,  29,  31,  28, &
! 33,  33,  32,  32,  35,  38,  40,  39,  40,  42, &
! 42,  43,  46,  52,  51,  50,  53,  54,  53,  55, &
! 56,  56,  58,  57,  59,  58,  58,  57,  56,  60, &
! 62,  63,  63,  71,  70,  69,  66,  66,  64,  64, &
! 63,  61,  60,  58,  55,  55,  55,  57,  60,  57, &
! 58,  54,  52,  54,  55,  56,  55,  57,  59,  59, &
! 60,  60,  60,  59,  56,  55,  53,  51,  53,  49, &
! 49,  55,  54,  61,  61,  63,  63,  63,  63,  62, &
! 61,  69,  75,  75,  76,  74,  73,  69,  66,  66, &
! 66,  67,  63,  63,  63,  59,  58,  60,  61,  63, &
! 63,  63,  64,  60,  62,  67,  64,  65,  66,  70, &
! 69,  69,  68,  68,  68,  67,  67,  66,  66,  66, &
! 68,  68,  65,  68,  66,  69,  66,  67,  64,  66, &
! 69,  64,  62,  61,  63,  60,  59,  62,  63,  64, &
! 64,  61,  59,  60,  58,  63,  61,  59,  59,  60, &
! 59,  58,  62,  63,  63,  63,  67,  65,  66,  68, &
! 67,  63,  63,  62,  62,  61,  57,  62,  61,  60 &
! /) ! e 200

  DOUBLE PRECISION, PARAMETER:: DataT(NUMD) = (/ &
     0.2D0, 0.4D0, 0.6D0, 0.8D0, 1.0D0, 1.2D0, 1.4D0, 1.6D0, 1.8D0, 2.0D0, &
     2.2D0, 2.4D0, 2.6D0, 2.8D0, 3.0D0, 3.2D0, 3.4D0, 3.6D0, 3.8D0, 4.0D0, & 
     4.2D0, 4.4D0, 4.6D0, 4.8D0, 5.0D0, 5.2D0, 5.4D0, 5.6D0, 5.8D0, 6.0D0, &
     6.2D0, 6.4D0, 6.6D0, 6.8D0, 7.0D0, 7.2D0, 7.4D0, 7.6D0, 7.8D0, 8.0D0, &
     8.2D0, 8.4D0, 8.6D0, 8.8D0, 9.0D0, 9.2D0, 9.4D0, 9.6D0, 9.8D0, 10.0D0 /)

  INTEGER, PARAMETER:: DataY(NUMD) = (/ & 
     0,  1,   5,    7,  12,  19,  32,  38,  43,  50,  50,  54, &
    60,  59,  60,  54,  61,  61,  55,  60,  58,  56,  54,  55, &
    58,  63,  65,  58,  64,  63,  61,  63,  61,  58,  54,  47, &
    53,  55,  59,  61,  51,  55,  52,  52,  46,  57,  59,  54,  59,  54/) !a

! INTEGER, PARAMETER:: DataY(NUMD) = (/ & 
!   1 ,  2 ,   5,  12,  18,  24,  21,  33,  34,  38,  48, &
!   42,  46,  46,  51,  49,  52,  56,  57,  55,  55,  60, &
!   66,  66,  69,  66,  62,  62,  64,  63,  62,  62,  65, &  
!   63,  68,  71,  68,  62,  58,  64,  65,  61,  63,  62, &  
!   59,  61,  62,  64,  60,  54 /) !b

! INTEGER, PARAMETER:: DataY(NUMD) = (/ & 
!    1,   1,   6,   9,  15,  21,  32,  37,  37,  39,  41,  40,  44,  52, &
!   48,  44,  44,  49,  53,  54,  48,  49,  42,  40,  52,  59,  58, &
!   58,  54,  61,  61,  62,  59,  61,  66,  67,  66,  59,  60,  54, &
!   57,  65,  65,  61,  54,  54,  50,  56,  60,  62 /) !c


!   DOUBLE PRECISION, PARAMETER:: DataT(NUMD) = (/ &
! 0.05D0,  0.1D0, 0.15D0,  0.2D0, 0.25D0,  0.3D0, 0.35D0,  0.4D0, 0.45D0,  0.5D0, &
! 0.55D0,  0.6D0, 0.65D0,  0.7D0, 0.75D0,  0.8D0, 0.85D0,  0.9D0, 0.95D0,  1.0D0, &
! 1.05D0,  1.1D0, 1.15D0,  1.2D0, 1.25D0,  1.3D0, 1.35D0,  1.4D0, 1.45D0,  1.5D0, & 
! 1.55D0,  1.6D0, 1.65D0,  1.7D0, 1.75D0,  1.8D0, 1.85D0,  1.9D0, 1.95D0,  2.0D0, &
! 2.05D0,  2.1D0, 2.15D0,  2.2D0, 2.25D0,  2.3D0, 2.35D0,  2.4D0, 2.45D0,  2.5D0, & 
! 2.55D0,  2.6D0, 2.65D0,  2.7D0, 2.75D0,  2.8D0, 2.85D0,  2.9D0, 2.95D0,  3.0D0, &
! 3.05D0,  3.1D0, 3.15D0,  3.2D0, 3.25D0,  3.3D0, 3.35D0,  3.4D0, 3.45D0,  3.5D0, & 
! 3.55D0,  3.6D0, 3.65D0,  3.7D0, 3.75D0,  3.8D0, 3.85D0,  3.9D0, 3.95D0,  4.0D0, & 
! 4.05D0,  4.1D0, 4.15D0,  4.2D0, 4.25D0,  4.3D0, 4.35D0,  4.4D0, 4.45D0,  4.5D0, &  
! 4.55D0,  4.6D0, 4.65D0,  4.7D0, 4.75D0,  4.8D0, 4.85D0,  4.9D0, 4.95D0,  5.0D0, & 
! 5.05D0,  5.1D0, 5.15D0,  5.2D0, 5.25D0,  5.3D0, 5.35D0,  5.4D0, 5.45D0,  5.5D0, &  
! 5.55D0,  5.6D0, 5.65D0,  5.7D0, 5.75D0,  5.8D0, 5.85D0,  5.9D0, 5.95D0,  6.0D0, & 
! 6.05D0,  6.1D0, 6.15D0,  6.2D0, 6.25D0,  6.3D0, 6.35D0,  6.4D0, 6.45D0,  6.5D0, &  
! 6.55D0,  6.6D0, 6.65D0,  6.7D0, 6.75D0,  6.8D0, 6.85D0,  6.9D0, 6.95D0,  7.0D0, & 
! 7.05D0,  7.1D0, 7.15D0,  7.2D0, 7.25D0,  7.3D0, 7.35D0,  7.4D0, 7.45D0,  7.5D0, &  
! 7.55D0,  7.6D0, 7.65D0,  7.7D0, 7.75D0,  7.8D0, 7.85D0,  7.9D0, 7.95D0,  8.0D0, & 
! 8.05D0,  8.1D0, 8.15D0,  8.2D0, 8.25D0,  8.3D0, 8.35D0,  8.4D0, 8.45D0,  8.5D0, &  
! 8.55D0,  8.6D0, 8.65D0,  8.7D0, 8.75D0,  8.8D0, 8.85D0,  8.9D0, 8.95D0,  9.0D0, & 
! 9.05D0,  9.1D0, 9.15D0,  9.2D0, 9.25D0,  9.3D0, 9.35D0,  9.4D0, 9.45D0,  9.5D0, & 
! 9.55D0,  9.6D0, 9.65D0,  9.7D0, 9.75D0,  9.8D0, 9.85D0,  9.9D0, 9.95D0,  10.0D0 /)

! INTEGER, PARAMETER:: DataY(NUMD) = (/ & 
!    0,   0,   0,   0,   0,   0,   0,   0,   1, & 
!    1,   1,   3,   5,   5,   6,   6,   7,   9, &   
!   11,  11,  12,  14,  15,  18,  19,  26,  27, &  
!   29,  32,  31,  33,  34,  34,  35,  34,  33, &  
!   32,  35,  39,  40,  40,  42,  42,  41,  43, &  
!   42,  44,  41,  46,  50,  47,  51,  49,  48, &  
!   48,  49,  46,  48,  53,  52,  48,  49,  47, &  
!   50,  50,  54,  54,  53,  56,  57,  56,  58, &  
!   55,  59,  58,  56,  58,  57,  59,  52,  54, &  
!   57,  56,  57,  57,  54,  54,  57,  52,  51, &  
!   51,  51,  51,  51,  51,  54,  56,  57,  58, &  
!   59,  61,  61,  63,  65,  64,  65,  63,  60, &  
!   65,  64,  63,  59,  61,  59,  58,  55,  54, &  
!   57,  60,  60,  64,  65,  61,  61,  62,  61, &  
!   62,  60,  60,  60,  62,  60,  57,  59,  60, &  
!   59,  60,  61,  59,  61,  60,  61,  63,  65, &  
!   62,  59,  59,  64,  65,  60,  60,  65,  64, &  
!   62,  57,  57,  54,  54,  57,  59,  60,  60, &  
!   67,  65,  62,  60,  64,  66,  65,  66,  67, &  
!   65,  58,  63,  61,  62,  63,  65,  64,  66, &  
!   64,  66,  61,  64,  67,  63,  66,  66,  66, &  
!   65,  66,  68,  72,  72,  72,  73,  71,  70,  63,  67/) !n=200


! DOUBLE PRECISION, PARAMETER:: DataT(NUMD) = (/ &
! 1.0D0, 2.0D0, 3.0D0, 4.0D0, 5.0D0, 6.0D0, 7.0D0, 8.0D0, 9.0D0, 10.0D0/)
! INTEGER, PARAMETER:: DataY(NUMD) = (/ & 
!   ! 18,  40,  44,  55,  49,  61,  61,  65,  54,  50 &
! 10,  38,  50,  50,  61,  46,  54,  63,  56,  56 &
! /) ! n = 10


  INTEGER :: i
  INTEGER, PARAMETER:: NP = 100
  REAL(KIND=R8)::A(NP,NP)=0.0_R8, E(NP,NP)
  REAL(KIND=R8):: KN, KM, KS, KD
  REAL(KIND=R8):: k, g, L, TEMP
  REAL(KIND=R8), DIMENSION(3):: p
  !REAL(KIND=R8)::ISEED,USEED
  !CALL RANDOM_NUMBER(ISEED)
  !WRITE(*,*)'SEED',ISEED
  !USEED = DUSTAR(INT(2.0_R8**24*ISEED))


  KS = 10.0**(PARAMS(1))
  KD = 10.0**(PARAMS(2))
  KM = 10.0**(PARAMS(3))
  KN = 10.0**(PARAMS(4))

  p(1) = 100.0_R8
  p(2) = 1000.0_R8
  p(3) = 0.0_R8

  DO i=0,NP-2
     
    IF (KN*log10(KM/p(2))>13) THEN
      k = KS*(p(1)-i)/(1.0e+13+1.0)
    ELSE
      k = KS*(p(1)-i)/((KM/p(2))**KN+1.0)
    END IF

    g = KD*i
    A(i+1,i+1) =  -k-g   ! Flow out.
    A(i+2,i+1) = k       ! Flow in due to production
    IF (i .GT. 0) THEN
        A(i,i+1) = g     ! Flow in due to degradation
    END IF
  END DO


  E = expm(t,A)
  L = 0.0
  DO i=1,NUMD-1
    IF (E(DataY(i+1)+1,DataY(i)+1) .LT. 1.0e-10) THEN
      TEMP = -100.0_R8
      !print *, DataY(i+1)+1, DataY(i)+1
    ELSE
      TEMP = log( E(DataY(i+1)+1,DataY(i)+1) )
    END IF
    L = L +  TEMP
  END DO
  ERROR = - L
!print *, 'Param :', PARAMS
!print *, 'E :', A

  IFLAG = 0

  
  ! WRITE(*,*) ERROR
  !  WRITE(*,789) SIMY
  !  WRITE(*,789) KN, KM, KS, KD
  !789 format('***: ', 4ES12.4)
 

  RETURN


 END FUNCTION HILL_OBJ


FUNCTION expm(t, H) result(expH)
  real(KIND=R8), intent(in) :: t
  real(KIND=R8), dimension(:,:), intent(in) :: H
  real(KIND=R8), dimension(size(H,1),size(H,2)) :: expH
  
  ! Expokit variables
  !external :: DGPADM
  integer, parameter :: ideg = 6
  real(KIND=R8), dimension(4*size(H,1)*size(H,2) + ideg + 1) :: wsp
  integer, dimension(size(H,1))  :: iwsp
  integer :: iexp, ns, iflag, n
  
  if (size(H,1) /= size(H,2)) then
  stop 'expm: matrix must be square'
  end if
  
  n = size(H,1)
  call DGPADM(ideg, n, t, H, n, wsp, size(wsp,1), iwsp, iexp, ns, iflag)
  expH = reshape(wsp(iexp:iexp+n*n-1), shape(expH))

END FUNCTION expm




END MODULE MODEL_HILL
