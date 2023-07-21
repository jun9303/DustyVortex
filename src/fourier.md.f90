MODULE FOURIER ! MODULE FOR FOURIER EXPANSION (IN Z DIRECTION) SUBROUTINES & FUNCTIONS
      USE OMP_LIB; USE MPI
      USE MISC; USE MATOPS
      USE LEGENDRE
      IMPLICIT NONE
      PRIVATE
! ==========================================================================================================
!  PARAMS/TYPES ============================================================================================
! ==========================================================================================================
      INTEGER, PUBLIC :: NZ                             ! NO. OF COLLOC. PTS IN THE PHI (AZIMUTH.) DIRECTION
      INTEGER, PUBLIC :: NZCHOP, NZCHOPU             ! NO. OF TOTAL & UPPER-EFFECTIVE SPECTRAL ELEMENTS IN Z
      INTEGER, PUBLIC :: NZDIM                                                  ! NO. OF DIMENSION IN Z DIR.
      INTEGER, PUBLIC :: NZCHOPDIM                                ! ORIGINAL (GLOBAL) CHOPING DIM. IN Z DIR.
      INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE :: AK                    ! ACCESSIBLE AZIMUTHAL WAVENUMBERS

      TYPE, PUBLIC :: FOU_TRANSFORM
        REAL(P8) :: ZLEN          ! COPMUTATIONAL DOMAIN LENGTH (OR THE LONGEST WAVELENGTH IN USE) IN Z DIR.
        REAL(P8), DIMENSION(:), ALLOCATABLE :: Z                        ! EQUISPACED COLLOCATION POINTS IN Z
      END TYPE
! ==========================================================================================================
!  DECLARATIONS ============================================================================================
! ==========================================================================================================
      PUBLIC :: FOU_INIT    ! INITIALIZE THE FOURIER-LEGENDRE SPECTRAL METHOD. LOAD LEG. TFM. KIT ON LEG_TFM
      PUBLIC :: FOU_CHOPSET                                  ! SET UP THE (GLOBAL) CHOPPING LOCATIONS WRT AK
      PUBLIC :: FOU_DZ                                                 ! MATRIX FOURIER SPACE OPERATOR: D/DZ
      PUBLIC :: FOU_DDZ                               ! MATRIX FOURIER SPACE OPERATOR: D^2/DZ^2 = D(D/DZ)/DZ
! ==========================================================================================================
!  INTERFACES ==============================================================================================
! ==========================================================================================================
      INTERFACE FOU_INIT
        MODULE PROCEDURE FOU_INIT
      END INTERFACE

      INTERFACE FOU_CHOPSET
        MODULE PROCEDURE FOU_CHOPSET
      END INTERFACE

      INTERFACE FOU_DZ
        MODULE PROCEDURE FOU_DZ_GEN
        MODULE PROCEDURE FOU_DZ_BND
      END INTERFACE

      INTERFACE FOU_DDZ
        MODULE PROCEDURE FOU_DDZ_GEN
        MODULE PROCEDURE FOU_DDZ_BND
      END INTERFACE

CONTAINS
! ==========================================================================================================
!  PUBLIC PROCEDURES =======================================================================================
! ==========================================================================================================
      FUNCTION FOU_INIT(NZIN, NZCHOPIN, ZLENIN) RESULT(TFM)
! ==========================================================================================================
! [USAGE]:
! INITIALIZE ALL PARAMETERS REQUIRED FOR FOURIER EXPANSIONS.
! [VARIABLES]:
! NZIN >> INTEGER VALUE USED FOR NZ (# OF COLLOCATION PTS IN THE Z DIRECTION)
! NZCHOPIN >> INTEGER VALUE USED FOR NZCHOP (# OF SPECTRAL ELEMENTS IN THE Z DIR. #CHOP. .LE. #COLLOC.)
! ZLENIN >> REAL VALUE USED FOR ZLEN (DOMAIN LENGTH IN Z)
! FOU_TFM >> FOURIER TRANSFOMRATION KIT. SEE WHAT VARIABLES ARE COLLECTED IN THIS TYPE
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: NZIN, NZCHOPIN
      REAL(P8), INTENT(IN) :: ZLENIN
      TYPE(FOU_TRANSFORM)  :: TFM

      INTEGER              :: I

      IF (MOD(NZIN, 2) .NE. 0) STOP 'FOU_INIT: NZ MUST BE EVEN FOR FASTER COMPUTATION.'
      IF (IAND(NZIN, NZIN-1) .NE. 0) WRITE(*,*) 'FOU_INIT: [WARN] 2^N PREFERABLE FOR NZ FOR EFFICIENCY.'
      IF (NZCHOPIN*2 - 1 .GT. NZIN)  STOP 'FOU_INIT: NZCHOP <= NZ/2 MUST BE SATISFIED.'

      NZ = NZIN
      NZCHOP = NZCHOPIN
      TFM%ZLEN = ZLENIN

      NZDIM = NZ

      NZCHOPDIM = MIN(NZCHOP, (NZ + 3)/2)
      NZCHOPU = NZ - NZCHOP + 2
      IF (NZCHOP .EQ. 1) NZCHOPU = 1

      CALL FOU_CHOPSET(0)

      ALLOCATE( AK(NZDIM) )
      DO I = 1, NZDIM
        IF (I .LE. NZCHOP) THEN
          AK(I) = 2.D0*PI/TFM%ZLEN * (I-1)
        ELSEIF (I .GE. NZCHOPU) THEN
          AK(I) = 2.D0*PI/TFM%ZLEN * (I-1-NZ)
        ELSE
          AK(I) = 0.D0
        ENDIF
      ENDDO

      ALLOCATE( TFM%Z(NZ+1) )
      TFM%Z = (/ (TFM%ZLEN*(I-1.D0)/NZ, I=1,NZ+1) /)

      RETURN
      END FUNCTION
! ==========================================================================================================
      SUBROUTINE FOU_CHOPSET(IOF)
! ==========================================================================================================
! [USAGE]:
! SET UP THE CHOP LOCATIONS IN Z
! [VARIABLES]:
! IOF >> OFFSET OF NZCHOP
! [NOTES]:
! SEE BELOW TO COMPREHEND HOW CHOPPING IS DONE.
!
! A PENCIL OF +---------------------+-----------------+-------------------+
! (R, P)   -> + 0 | 1 | .... | NC-1 +   INACTIVE(0)   + -NC+1 | .... | -1 +   
! (NXC X NPC) +---------------------+-----------------+-------------------+
!                                ↑                       ↑              ↑
!                           NZCHOP(NZCHOPDIM)          NZCHOPU        NZ(NZDIM)
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IOF

      NZCHOP = NZCHOP + IOF
      NZCHOPU = NZCHOPU - IOF

      IF (2 * IOF .GT. NZCHOPU - NZCHOP) STOP 'LEG_CHOPSET: NZCHOP TOO LARGE.'

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      FUNCTION FOU_DZ_GEN(LI, LO, DLEN, CHOP) RESULT(FOU_DZ)
! ==========================================================================================================
! [USAGE]:
! MATRIX OPERATOR TO MULTIPLY D/DZ IN FOURIER FUNCTION SPACE
! [VARIABLES]:
! LI >> LENGTH OF (A^M) (M = 0, 1, ... , NC-1, ..., -NC+1, ... -1)
! LO >> LENGTH OF (B^M) (M = 0, 1, ... , NC-1, ..., -NC+1, ... -1)
! DLEN >> DOMAIN LENGTH FOR FOURIER EXPANSION.
! CHOP >> (OPTIONAL) # OF FOURIER ELEMENTS OF POSITIVE WNS. IN USE. IF NOT PRESENT, SET AS MAXIMUM
! FOU_DZ >> MATRIX EXPRESSION OF THE OPERATOR MULTIPLYING D/DZ IN FOURIER FUNCTION SPACE
! [NOTES]:
! 1. FOR A FUNCTION F(Z) EXPENDED WITH COEFFCIEINT A^M AS
! F(Z) = SUM_N=-INFTY^INFTY A^M EXP(IU*2*PI*N*Z/ZLEN)
! COEFFICIENTS OF THE FUNCTION D/DZ(F(X)) B^M AS
! D/DZ(F(Z)) = SUM_N=-INFTY^INFTY B^M EXP(IU*2*PI*N*Z/ZLEN)
! CAN BE EVAULATED BY "B^M = FOU_DZ * A^M"
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)                      :: LI, LO
      REAL(P8), INTENT(IN)                     :: DLEN
      INTEGER, OPTIONAL                        :: CHOP
      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: FOU_DZ

      INTEGER                                  :: I, M, NN
      INTEGER                                  :: CHOP_
      REAL(P8), DIMENSION(:), ALLOCATABLE      :: K

      CHOP_ = (MIN(LI, LO) + 1) / 2

      IF (PRESENT(CHOP)) CHOP_ = CHOP
      IF (2 * CHOP_ - 1 .GT. MIN(LI, LO)) STOP 'FOU_DDZ: INVALID PARAM SET (2*CHOP-1 > MIN(LI,LO)).'


      ALLOCATE( FOU_DZ(LO, LI) )
      FOU_DZ = 0.D0

      M = MIN(LI, LO)
      ALLOCATE( K(M) )
      DO I = 1, M
        IF (I .LE. CHOP_) THEN
          K(I) = 2.D0*PI/DLEN * (I-1)
        ELSEIF (I .GE. M - CHOP_ + 2) THEN
          K(I) = 2.D0*PI/DLEN * (I-1-NX)
        ELSE
          K(I) = 0.D0
        ENDIF
      ENDDO

      DO NN = 1, M
        IF (NN .LE. LI) FOU_DZ(NN, NN) = IU * K(NN) ! MAIN DIAGONAL
      ENDDO

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION FOU_DZ_BND(LI, DLEN, CHOP) RESULT(FOU_DZ)
! ==========================================================================================================
! [USAGE]:
! BAND MATRIX OPERATOR TO MULTIPLY D/DZ IN FOURIER FUNCTION SPACE
! [VARIABLES]:
! LI >> LENGTH OF (A^M) (M = 0, 1, ... , NC-1, ..., -NC+1, ... -1)
! DLEN >> DOMAIN LENGTH FOR FOURIER EXPANSION.
! CHOP >> (OPTIONAL) # OF FOURIER ELEMENTS OF POSITIVE WNS. IN USE. IF NOT PRESENT, SET AS MAXIMUM
! FOU_DDZ >> MATRIX EXPRESSION OF THE OPERATOR MULTIPLYING D/DZ IN FOURIER FUNCTION SPACE
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: LI
      REAL(P8), INTENT(IN) :: DLEN
      INTEGER, OPTIONAL    :: CHOP
      TYPE(CMPX_BNDM)      :: FOU_DZ

      IF (PRESENT(CHOP)) THEN
        FOU_DZ = BANDMAT(FOU_DZ_GEN(LI, LI, DLEN, CHOP), KL=0, KU=0)
      ELSE
        FOU_DZ = BANDMAT(FOU_DZ_GEN(LI, LI, DLEN), KL=0, KU=0)
      ENDIF

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION FOU_DDZ_GEN(LI, LO, DLEN, CHOP) RESULT(FOU_DDZ)
! ==========================================================================================================
! [USAGE]:
! MATRIX OPERATOR TO MULTIPLY D(D/DZ)/DZ IN FOURIER FUNCTION SPACE
! [VARIABLES]:
! LI >> LENGTH OF (A^M) (M = 0, 1, ... , NC-1, ..., -NC+1, ... -1)
! LO >> LENGTH OF (B^M) (M = 0, 1, ... , NC-1, ..., -NC+1, ... -1)
! DLEN >> DOMAIN LENGTH FOR FOURIER EXPANSION.
! CHOP >> (OPTIONAL) # OF FOURIER ELEMENTS OF POSITIVE WNS. IN USE. IF NOT PRESENT, SET AS MAXIMUM
! FOU_DDZ >> MATRIX EXPRESSION OF THE OPERATOR MULTIPLYING D(D/DZ)/DZ IN FOURIER FUNCTION SPACE
! [NOTES]:
! 1. FOR A FUNCTION F(Z) EXPENDED WITH COEFFCIEINT A^M AS
! F(Z) = SUM_N=-INFTY^INFTY A^M EXP(IU*2*PI*N*Z/ZLEN)
! COEFFICIENTS OF THE FUNCTION D/DZ(F(X)) B^M AS
! D/DZ(F(Z)) = SUM_N=-INFTY^INFTY B^M EXP(IU*2*PI*N*Z/ZLEN)
! CAN BE EVAULATED BY "B^M = FOU_DZ * A^M"
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)                   :: LI, LO
      REAL(P8), INTENT(IN)                  :: DLEN
      INTEGER, OPTIONAL                     :: CHOP
      REAL(P8), DIMENSION(:,:), ALLOCATABLE :: FOU_DDZ

      INTEGER                               :: I, M, NN
      INTEGER                               :: CHOP_
      REAL(P8), DIMENSION(:), ALLOCATABLE   :: K

      CHOP_ = (MIN(LI, LO) + 1) / 2

      IF (PRESENT(CHOP)) CHOP_ = CHOP
      IF (2 * CHOP_ - 1 .GT. MIN(LI, LO)) STOP 'FOU_DDZ: INVALID PARAM SET (2*CHOP-1 > MIN(LI,LO)).'


      ALLOCATE( FOU_DDZ(LO, LI) )
      FOU_DDZ = 0.D0

      M = MIN(LI, LO)
      ALLOCATE( K(M) )
      DO I = 1, M
        IF (I .LE. CHOP_) THEN
          K(I) = 2.D0*PI/DLEN * (I-1)
        ELSEIF (I .GE. M - CHOP_ + 2) THEN
          K(I) = 2.D0*PI/DLEN * (I-1-NX)
        ELSE
          K(I) = 0.D0
        ENDIF
      ENDDO

      DO NN = 1, M
        IF (NN .LE. LI) FOU_DDZ(NN, NN) = - K(NN)**2. ! MAIN DIAGONAL
      ENDDO

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION FOU_DDZ_BND(LI, DLEN, CHOP) RESULT(FOU_DDZ)
! ==========================================================================================================
! [USAGE]:
! BAND MATRIX OPERATOR TO MULTIPLY D(D/DZ)/DZ IN FOURIER FUNCTION SPACE
! [VARIABLES]:
! LI >> LENGTH OF (A^M) (M = 0, 1, ... , NC-1, ..., -NC+1, ... -1)
! DLEN >> DOMAIN LENGTH FOR FOURIER EXPANSION.
! CHOP >> (OPTIONAL) # OF FOURIER ELEMENTS OF POSITIVE WNS. IN USE. IF NOT PRESENT, SET AS MAXIMUM
! FOU_DDZ >> MATRIX EXPRESSION OF THE OPERATOR MULTIPLYING D(D/DZ)/DZ IN FOURIER FUNCTION SPACE
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: LI
      REAL(P8), INTENT(IN) :: DLEN
      INTEGER, OPTIONAL    :: CHOP
      TYPE(REAL_BNDM)      :: FOU_DDZ

      IF (PRESENT(CHOP)) THEN
        FOU_DDZ = BANDMAT(FOU_DDZ_GEN(LI, LI, DLEN, CHOP), KL=0, KU=0)
      ELSE
        FOU_DDZ = BANDMAT(FOU_DDZ_GEN(LI, LI, DLEN), KL=0, KU=0)
      ENDIF

      RETURN
      END FUNCTION
! ==========================================================================================================
END MODULE