MODULE TIMEINFO ! MODULE FOR TIME INFORMATION SUBROUTINES & FUNCTIONS
      USE OMP_LIB; USE MPI
      USE MISC
      PRIVATE
! ==========================================================================================================
!  PARAMS/TYPES ============================================================================================
! ==========================================================================================================
      TYPE :: TIME_INFO
        REAL(P8) :: TI, T, DT, DTF, TOTT
        INTEGER :: NI, N, TOTN
      END TYPE
      TYPE(TIME_INFO), PUBLIC :: TINFO
! ==========================================================================================================
!  DECLARATIONS ============================================================================================
! ==========================================================================================================
      PUBLIC :: TIME_INIT
! ==========================================================================================================
!  INTERFACES ==============================================================================================
! ==========================================================================================================
      INTERFACE TIME_INIT
        MODULE PROCEDURE TIME_INIT
      END INTERFACE

CONTAINS
! ==========================================================================================================
!  PUBLIC PROCEDURES =======================================================================================
! ==========================================================================================================
      SUBROUTINE TIME_INIT(DTIN, TIIN, TOTTIN, NIIN, TOTNIN)
! ==========================================================================================================
! [USAGE]:
! INITIALIZE THE TIME INFORMATION FOR IVP.
! [VARIABLES]:
! DTIN >> DELTA TIME INPUT.
! TIIN >> (OPTIONAL) INITIAL TIME INPUT.
! NIIN >> (OPTIONAL) INITIAL STEP #
! TOTTIN >> (OPTIONAL) TERMINATION CONDITION IN TERMS OF TOTAL TIME IN IVP
! TOTNIN >> (OPTIONAL) TERMINATION CONDITION IN TERMS OF TOTAL NO. OF STEPS RUN IN IVP
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), INTENT(IN) :: DTIN
      REAL(P8), OPTIONAL   :: TIIN, TOTTIN
      INTEGER, OPTIONAL    :: NIIN, TOTNIN

      INTEGER              :: NIIN_ = 0
      INTEGER              :: TOTNIN_ = 1
      REAL(P8)             :: TIIN_ = 0.D0
      REAL(P8)             :: TOTTIN_ = 1.D-10

      IF (DTIN .LE. 0.D0) STOP 'TIME_INIT: DELTA T MUST BE POSITIVE.'

      IF (PRESENT(TIIN)) TIIN_ = TIIN
      IF (PRESENT(TOTTIN)) TOTTIN_ = TOTTIN
      IF (PRESENT(NIIN)) NIIN_ = NIIN
      IF (PRESENT(TOTNIN)) TOTNIN_ = TOTNIN

      TINFO%DT = DTIN
      TINFO%TI = TIIN_
      TINFO%NI = NIIN_

      TINFO%TOTT = MAX(TOTTIN, TINFO%DT*TOTNIN)

      IF (TINFO%TOTT .EQ. TOTTIN) THEN
        TINFO%TOTN = FLOOR(TOTTIN/TINFO%DT)
      ELSE
        TINFO%TOTN = TOTNIN
      ENDIF

      TINFO%DTF = TINFO%TOTT - (TINFO%TOTN-1)*TINFO%DT

      TINFO%T = TINFO%TI
      TINFO%N = TINFO%NI

      RETURN
      END SUBROUTINE
! ==========================================================================================================
END MODULE