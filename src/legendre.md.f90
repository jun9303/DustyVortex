MODULE LEGENDRE ! MODULE FOR LEGENDRE EXPANSION (IN R-P POLAR COORDS) SUBROUTINES & FUNCTIONS
      USE OMP_LIB
      USE MISC; USE MATOPS
      IMPLICIT NONE
      PRIVATE
! ==========================================================================================================
!  PARAMS/TYPES ============================================================================================
! ==========================================================================================================
      INTEGER, PUBLIC :: NX, NXH                      ! NO. OF COLLOC. PTS. IN THE R/X DIRECTION; NXH = NX/2
      INTEGER, PUBLIC :: NP, NPH            ! NO. OF COLLOC. PTS IN THE PHI (AZIMUTH.) DIRECTION; NPH = NP/2
      INTEGER, PUBLIC :: NXCHOP, NPCHOP     ! NO. OF SPECTRAL ELEMENTS IN EACH DIR. LESS THAN # COLLOC. PTS.
      INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE :: NXCHOPS, NPCHOPS       ! CHOP INDICES IN EACH DIR. WRT M
      INTEGER, PUBLIC :: NXDIM, NPDIM                                       ! NO. OF DIMENSIONS IN EACH DIR.
      INTEGER, PUBLIC :: NXCHOPDIM, NPCHOPDIM              ! ORIGINAL DIMENSIONS IN EACH DIR. AFTER CHOPPING
      INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE :: M                     ! ACCESSIBLE AZIMUTHAL WAVENUMBERS

      TYPE, PUBLIC :: LEG_TRANSFORM
        REAL(P8) :: L                                                         ! MAP PARAMETER (0 < L < INF.)
        REAL(P8), DIMENSION(:), ALLOCATABLE :: X      ! UN-MAPPED, ORIGINAL GAUSS-LEGENDRE ABSCISSA (-1<X<1)
        REAL(P8), DIMENSION(:), ALLOCATABLE :: R             ! MAPPED COLLOCATION POINTS IN R (0 < R < INF.)
        REAL(P8), DIMENSION(:), ALLOCATABLE :: P, PR, PI  ! EQUISPACED COLLOCATAION POINTS IN PHI (AZIMUTH.)
        REAL(P8), DIMENSION(:), ALLOCATABLE :: W                                    ! GAUSS-LEGENDRE WEIGHTS
        REAL(P8), DIMENSION(:,:), ALLOCATABLE :: NORM, LOGNORM              ! NORMALIZATION FACTORS OF P^M_N
        REAL(P8), DIMENSION(:,:,:), ALLOCATABLE :: PF ! LEGENDRE VALUE TABLE FOR P^M_N(X_I) GIVEN M, N AND I
        REAL(P8), DIMENSION(:), ALLOCATABLE :: AT0, AT1     ! LEGENDRE VALUES AT X=-1 (R=0) & AT X=1 (R=INF)
        REAL(P8), DIMENSION(:), ALLOCATABLE :: LN                                    ! MINUS LOG VALUES OF X
      END TYPE
! ==========================================================================================================
!  DECLARATIONS ============================================================================================
! ==========================================================================================================
      PUBLIC :: LEG_INIT                                        ! INITIALIZE THE LEGENDRE TRANSFORMATION KIT
      PUBLIC :: LEG_CHOPSET                                   ! SET UP THE (GLOBAL) CHOPPING LOCATIONS WRT M
      PUBLIC :: LEG_LOGNORM             ! CALCULATE AND STORE THE LOG NORMS OF THE LEGENDRE FUNCTIONS IN USE
      PUBLIC :: LEG_TBL            ! CALCULATE AND STORE THE LEGENDRE FUNCTION VALUES AT SPECIFIED LOCATIONS
      PUBLIC :: LEG_XM                          ! MATRIX LEGENDRE SPACE OPERATOR: (1-X)* <-> 2R^2/(R^2+L^2)*
      PUBLIC :: LEG_XP                          ! MATRIX LEGENDRE SPACE OPERATOR: (1+X)* <-> 2L^2/(R^2+L^2)*
      PUBLIC :: LEG_X                          ! MATRIX LEGENDRE SPACE OPERATOR: X* <-> (R^2-L^2)/(R^2+L^2)*
      PUBLIC :: LEG_XXDX                      ! MATRIX LEGENDRE SPACE OPERATOR: (1-X**2)*D/DX() <-> R*D/DR()
      PUBLIC :: LEG_DEL2H     ! MATRIX LEGENDRE SPACE OPERATOR: DEL^2_PERP(M) (1/R*D/DR(R*D/DR()) - M^2/R^2)
      PUBLIC :: LEG_DEL2  ! MATRIX LEGENDRE SPACE OPERATOR: DEL^2(M, AK) (1/R*D/DR(R*D/DR()) - M^2/R^2-AK^2)
! ==========================================================================================================
!  INTERFACES ==============================================================================================
! ==========================================================================================================
      INTERFACE LEG_INIT
        MODULE PROCEDURE LEG_INIT
      END INTERFACE

      INTERFACE LEG_CHOPSET
        MODULE PROCEDURE LEG_CHOPSET
      END INTERFACE

      INTERFACE LEG_LOGNORM                                 ! TABULATE THE LOG OF THE NORMALIZATION FACTORS
        MODULE PROCEDURE LEG_LOGNORM_MULTI                            ! LOG(NORM) OF P^M_N FOR MULTIPLE M'S
        MODULE PROCEDURE LEG_LOGNORM_SINGLE                                  ! LOG(NORM) OF P^M_N FOR ONE M
      END INTERFACE

      INTERFACE LEG_TBL                                             ! TABULATE THE LEGENDRE FUNCTION VALUES
        MODULE PROCEDURE LEG_TBL_1D                              ! FOR FIXED M, COMPUTE P^M_N VAL @ POINT X
        MODULE PROCEDURE LEG_TBL_2D                                ! FOR FIXED M, COMPUTE P^M_N VAL AT X(:)
        MODULE PROCEDURE LEG_TBL_3D                                      ! MAKE FULL TABLE FOR MULTIPLE M'S
      END INTERFACE

      INTERFACE LEG_XM
        MODULE PROCEDURE LEG_XM_GEN
        MODULE PROCEDURE LEG_XM_BND
      END INTERFACE

      INTERFACE LEG_XP
        MODULE PROCEDURE LEG_XP_GEN
        MODULE PROCEDURE LEG_XP_BND
      END INTERFACE

      INTERFACE LEG_X
        MODULE PROCEDURE LEG_X_GEN
        MODULE PROCEDURE LEG_X_BND
      END INTERFACE

      INTERFACE LEG_XXDX
        MODULE PROCEDURE LEG_XXDX_GEN
        MODULE PROCEDURE LEG_XXDX_BND
      END INTERFACE

      INTERFACE LEG_DEL2H
        MODULE PROCEDURE LEG_DEL2H_GEN
        MODULE PROCEDURE LEG_DEL2H_BND
      END INTERFACE

      INTERFACE LEG_DEL2
        MODULE PROCEDURE LEG_DEL2_GEN
        MODULE PROCEDURE LEG_DEL2_BND
      END INTERFACE

CONTAINS
! ==========================================================================================================
!  PUBLIC PROCEDURES =======================================================================================
! ==========================================================================================================
      FUNCTION LEG_INIT(NXIN, NPIN, NXCHOPIN, NPCHOPIN, ELLIN, VISCPOW) RESULT(TFM)
! ==========================================================================================================
! [USAGE]:
! INITIALIZE ALL PARAMETERS REQUIRED FOR THE MAPPED LEGENDRE SPECTRAL METHOD, ALL OF WHICH ARE STORED UNDER
! THE LEG_TRANSFORM TYPE LEG_TFM.
! [VARIABLES]:
! NXIN >> INTEGER VALUE USED FOR NX (# OF COLLOCATION PTS IN THE X/R DIRECTION)
! NPIN >> INTEGER VALUE USED FOR NPH (# OF COLLOCATION PTS IN THE PHI DIRECTION)
! NXCHOPIN >> INTEGER VALUE USED FOR NXCHOP (# OF SPECTRAL ELEMENTS IN THE X/R DIR. #CHOP. .LE. #COLLOC.)
! NPCHOPIN >> INTEGER VALUE USED FOR NPCHOP (# OF SPECTRAL ELEMENTS IN THE PHI DIR. #CHOP. .LE. #COLLOC.)
! ELLIN >> REAL VALUE USED FOR ELL (LEGENDRE FUNCTIONS` MAP PARAMETER)
! VISCPOW >> (OPTIONAL) THE HIGHEST ORDER OF (HYPER)VISCOUS TERMS TO BE CONSIDERED. DEAFULT IS 2.
! LEG_TFM >> LEGENDRE TRANSFOMRATION KIT. SEE WHAT VARIABLES ARE COLLECTED IN THIS TYPE
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: NXIN, NPIN, NXCHOPIN, NPCHOPIN
      REAL(P8), INTENT(IN) :: ELLIN
      INTEGER, OPTIONAL    :: VISCPOW
      TYPE(LEG_TRANSFORM)  :: TFM

      INTEGER              :: I, J
      INTEGER              :: VISCPOW_

      VISCPOW_ = 2

      IF (PRESENT(VISCPOW)) VISCPOW_ = VISCPOW
      IF (MOD(VISCPOW_, 2) .NE. 0) STOP 'LEGINIT: VISCOUS TERM ORDER (VISCPOW) MUST BE EVEN.'
      IF (VISCPOW_ .GT. 8) STOP 'LEGINIT: VISCOUS TERM ORDER (VISCPOW) LESS THAN 8 ACCEPTED.'

      IF (MOD(NXIN, 2) .NE. 0) STOP 'LEGINIT: N(X/R) MUST BE EVEN FOR FASTER COMPUTATION.'
      IF (MOD(NPIN, 2) .NE. 0) STOP 'LEGINIT: NP MUST BE EVEN.'
      IF (IAND(NPIN, NPIN-1) .NE. 0) THEN
        IF (WARN_ON) WRITE(*,*) 'LEGINIT: [WARN] 2^N PREFERABLE FOR NP FOR EFFICIENCY.'
      ENDIF
      IF (NXCHOPIN .GT. NXIN) STOP 'LEGINIT: N(X/R)CHOP MUST BE SMALLER THAN OR EQUAL TO N(X/R).'
      IF (NPCHOPIN*2 - 1 .GT. NPIN)  STOP 'LEG_INIT: NPCHOP <= NP/2 MUST BE SATISFIED.'

      TFM%L = ELLIN
      ALLOCATE( TFM%X(NXIN) )
      ALLOCATE( TFM%W(NXIN) )
      ALLOCATE( TFM%R(NXIN) )
      CALL FIND_GAULEG_ZEROS(TFM%X, TFM%W, TFM%R, TFM%L)

      NX = NXIN
      NXH = NX/2
      NP = NPIN
      NPH = NP/2
      NXCHOP = NXCHOPIN
      NPCHOP = NPCHOPIN

      NXDIM = NX + MAX(3, VISCPOW_)
      NPDIM = NPH + 1
      IF (NPDIM .EQ. 2) NPDIM = 1

      NXCHOPDIM = MIN(NXCHOP, NX)
      NPCHOPDIM = MIN(NPCHOP, NPH+1)

      ALLOCATE( M(NPCHOP) )
      M = (/ (I, I = 0, NPCHOP-1) /)

      ALLOCATE( NXCHOPS(NPH+1) )
      ALLOCATE( NPCHOPS(NX) )
      CALL LEG_CHOPSET(0)

      ALLOCATE( TFM%LN(NX) )
      TFM%LN = -LOG(1.D0 - TFM%X)

      ALLOCATE( TFM%LOGNORM(NXCHOP+14, NPCHOP) )
      TFM%LOGNORM = LEG_LOGNORM(NXCHOP+14, M)

      ALLOCATE( TFM%NORM(NXCHOP+14, NPCHOP) )
      DO I = 1, SIZE(TFM%NORM, 1)
        DO J = 1, SIZE(TFM%NORM, 2)
          TFM%NORM(I,J) = EXP(TFM%LOGNORM(I,J))
        ENDDO
      ENDDO

      ALLOCATE( TFM%PF(NXH, NXCHOPDIM+14, NPCHOPDIM) )
      TFM%PF = LEG_TBL(TFM%X(:NXH),NXCHOPDIM+14, M ,TFM%LOGNORM)

      ALLOCATE( TFM%AT0(NXCHOPDIM) )
      TFM%AT0 = LEG_TBL( -1.D0, NXCHOPDIM, 0, TFM%LOGNORM(:,1))
      
      ALLOCATE( TFM%AT1(NXCHOPDIM) )
      TFM%AT1 = LEG_TBL(  1.D0, NXCHOPDIM, 0, TFM%LOGNORM(:,1))

      ALLOCATE( TFM%PR(NPH+1) )
      TFM%PR = (/ (2*PI/NPH*(I-1.D0), I=1, NPH+1) /)
      
      ALLOCATE( TFM%PI(NPH+1) )
      TFM%PI = (/ (2*PI/NPH*(I-.5D0), I=1, NPH+1) /)

      ALLOCATE( TFM%P(2*NPH+1) )
      TFM%P = (/ (PI/NPH*(I-1.D0), I=1, 2*NPH+1) /)

      RETURN
      END FUNCTION
! ==========================================================================================================
      SUBROUTINE LEG_CHOPSET(IOF, IOF2)
! ==========================================================================================================
! [USAGE]:
! SET UP THE CHOP LOCATIONS WITH RESPECT TO AZIMUTHAL WNS (M)
! [VARIABLES]:
! IOF >> OFFSET OF NZCHOP
! IOF2 >> OFFSET OF NPCHOP (OPTIONAL)
! [NOTES]:
! SEE BELOW TO COMPREHEND HOW CHOPPING IS DONE.
!     NXDIM-> +----------------------------------------+---------+----+
!             +     RESERVED FOR HIGHER CHOPSET        +   INAC  +    +
!        NX-> +----------------------------------------+---------+----+
!             +             INACTIVE(0)                +         +    +
!             +                                        +    I    +    +
!             +  ↑                                     +    N    +    +
!    NXCHOP-> +--+ . . ↑(IOF) . . . . . . . . . . . . .+    A    +    +
! (NXCHOPDIM) +  +-----+     ↑                         +    C    +    +
!             +        +-----+     ↑                   +    T    +    +
!             +              +-----+     ↑             +    I    +    +
!             +                    +-----+ GRAD=-1     +    V    +    +
!             +                          +-----+     ↑ +    E    +    +
!             +         ACTIVE COEFFICIENTS    +-----+ +   (0)   +    +
!             +                                      +-+         +    +
!             +----------------------------------------+---------+----+
!                                                      ↑         ↑    ↑
!                                           (NPCHOPDIM)NPCHOP    NPH  NPDIM (RESERVED FOR REAL-2-CPLX DFT)
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IOF
      INTEGER, OPTIONAL, INTENT(IN) :: IOF2

      INTEGER             :: MM, NN

      NXCHOP = NXCHOP + IOF
      IF (PRESENT(IOF2)) NPCHOP = NPCHOP + IOF2

      IF (NXCHOP .GT. NXDIM) STOP 'LEG_CHOPSET: NXCHOP TOO LARGE.'
      IF (NPCHOP .GT. NPDIM) STOP 'LEG_CHOPSET: NPCHOP TOO LARGE.'

      DO MM = 1, NPCHOP
        NXCHOPS(MM) = MAX(MIN(NXCHOP, NXCHOP-M(MM)), 0)
      ENDDO
      NXCHOPS(NPCHOP+1:) = 0

      DO NN = 1, NX
        NPCHOPS(NN) = MAX(MIN(NPCHOP, NXCHOP-NN), 0)
      ENDDO

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      FUNCTION LEG_LOGNORM_MULTI(NDIM, MS) RESULT(LEG_LOG_NORM)
! ==========================================================================================================
! [USAGE]:
! FIND THE LOG OF THE NORMALIZATION FACTORS OF P_N^M FOR MULTIPLE AZIMUTHAL WNS IN MS.
! EXAMPLE. F = LEG_LOG_NORM (64, /0, MS/)) <-- ARRAY INPUT (/0, M/) = MS
! [VARIABLES]:
! NDIM >> THE (MAXIMUM) DEGREE OF THE ASSOCIATED LEGENDRE FUNCTION IN USE
! MS >> 1-D ARRAY OF THE ORDER OF THE ASSOCIATED LEGENDRE FUNCTIONS IN USE
! LEG_LOG_NORM >> LOG NORMALIZATION FACTOR TABLE. FOR A COLUMN CORRESPONDING TO M, IT CALCULATES
!                 (/ LOG(NORM^M_|M|), LOG(NORM^M_|M|+1), ... ,LOG(NORM^M_|M|+NDIM-1)) ) /)^T
! [NOTES]:
! 1. CHECK THE FOLLOWING RECURSIVE FORMULAS
! P^L_L  (X) = (-1)**L*(2*L-1)!!*(1-X**2)**(L/2)
! P^L_L+1(X) = X*(2*L+1)*P^L_L(X) 
! P^K_L(X)   = 1/(L-K)*[(2*L-1)*X*P^K_L-1(X) - (L+K-1)*P^K_L-2(X)]
! TO OBTAIN THE NORMALIZATION FACTORS AS FOLLOWS:
! N^M_N = SQRT(K(2N+1)*(N-|M|!)/(N+|M|)!) WHERE K = 1 FOR M=0; 2 OTHER.
! 2. IN THIS FUNCTION RECIPROCAL IS OUTPUT OF N^M_N.
! GIVEN N = NN-1 & M = MS(MM), LEG_NORMT(NN,MM) * P^M_N(R)
! IS AN NORMALIZED FUNCTION FOR THE ORDER M AND DEGREE N!
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)                     :: NDIM
      INTEGER, DIMENSION(:), INTENT(IN)       :: MS
      REAL(P8), DIMENSION(1:NDIM, 1:SIZE(MS)) :: LEG_LOG_NORM

      INTEGER                                 :: NN, MM, N, M, ME
      REAL(P8), DIMENSION(:), ALLOCATABLE     :: WK_LOG

      ME = MAXVAL(ABS(MS))
      ALLOCATE(WK_LOG(0:ME))

      WK_LOG(0) = LOG(0.5D0)

      DO M = 1, ME
        ! EQN. (13) (SEE MATSUSHIMA & MARCUS, 1997) WITHOUT SQRT
        WK_LOG(M) = WK_LOG(M-1)+LOG(2*M+1.D0)-LOG(2*M*(2*M-1.D0)**2)  
      ENDDO

      DO MM = 1, SIZE(MS)
        M = ABS(MS(MM))
        LEG_LOG_NORM(1,MM) = WK_LOG(M)
      ENDDO

      DO MM = 1, SIZE(MS)
        M = ABS(MS(MM))
        DO NN = 2, NDIM
          N = M+(NN-1)
          LEG_LOG_NORM(NN,MM) = LEG_LOG_NORM(NN-1,MM)&
                                +LOG(2*N+1.D0)-LOG(2*N-1.D0)&
                                +LOG(1.D0*(N-M))-LOG(1.D0*(N+M))
        ENDDO
      ENDDO

      DO MM = 1, SIZE(MS)
        LEG_LOG_NORM(:,MM) = .5D0*LEG_LOG_NORM(:,MM)
      ENDDO

      DEALLOCATE(WK_LOG)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION LEG_LOGNORM_SINGLE(NDIM, M) RESULT(LEG_LOG_NORM)
! ==========================================================================================================
! [USAGE]:
! FIND THE LOG OF THE NORMALIZATION FACTORS OF P_N^M FOR A SINGLE AZIMUTHAL WN M.
! EXAMPLE. F = LEG_LOG_NORM (64, M) <-- SCALAR INPUT M
! [VARIABLES]:
! NDIM >> THE (MAXIMUM) DEGREE OF THE ASSOCIATED LEGENDRE FUNCTION IN USE
! M >> THE ORDER OF THE ASSOCIATED LEGENDRE FUNCTIONS
! LEG_LOG_NORM >> LOG NORMALIZATION FACTOR TABLE. IT CALCULATES
!                 (/ LOG(NORM^M_|M|), LOG(NORM^M_|M|+1), ... ,LOG(NORM^M_|M|+NDIM-1)) ) /)^T
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)         :: NDIM
      INTEGER, INTENT(IN)         :: M
      REAL(P8), DIMENSION(1:NDIM) :: LEG_LOG_NORM

      LEG_LOG_NORM = RESHAPE( LEG_LOGNORM_MULTI(NDIM, (/ ABS(M) /)), (/ NDIM /) )

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION LEG_TBL_3D(X, NE, MS, LOG_NORM) RESULT(LEG_TBL)
! ==========================================================================================================
! [USAGE]:
! MAKE A FULL TABLE OF THE ASSOCIATED LEGENDRE FUNCTIONS
! [VARIABLES]:
! X >> MULTIPLE POINTS (EXAMPLE. COLLOCATION POINTS) TO EVAULATE LEG. FUN. VALS.
! NE >> MAX. DEGREE OF THE ASSOCIATED LEGENDRE FUNCTION
! MS >> SET OF ORDERS OF THE ASSOCIATED LEGENDRE FUNCTIONS
! LOG_NORM >> TABLE OF THE LOG NORMALIZATION FACTORS, FROM LEG_LOGNORM()
! LEG_TBL >> THE LEGENDRE FUNCTION COEFFICIENTS AT X & M FROM N = 1 TO NE. FOR GIVEN X AND M, THE COLUMN IS
!            LEG_TBL(X,:,M) = (/ P^M_|M|(X), ..., P^M_|M|+NE-1(X) /)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      USE FMVALS ! EXTERNAL LIB. REQUIRED -- MULTIPLE PRECISION (FM PACKAGE)
      USE FMZM   ! EXTERNAL LIB. REQUIRED -- MULTIPLE PRECISION (FM PACKAGE)

      IMPLICIT NONE

      REAL(P8), DIMENSION(:), INTENT(IN)         :: X
      INTEGER, INTENT(IN)                        :: NE
      INTEGER, DIMENSION(:), INTENT(IN)          :: MS
      REAL(P8), DIMENSION(:,:), INTENT(IN)       :: LOG_NORM
      REAL(P8), DIMENSION(SIZE(X), NE, SIZE(MS)) :: LEG_TBL

      INTEGER                                    :: M, ME, MM, NN, N, XX
      TYPE(FM), ALLOCATABLE, DIMENSION(:,:)      :: FM_LEG_TBL

      CALL FM_SET(50)

      LEG_TBL = 0.D0
      ME = SIZE(MS)

      ALLOCATE( FM_LEG_TBL(SIZE(X), NE) )

      DO MM = 1, ME
        M = ABS(MS(MM))
        IF (M .EQ. 0) THEN
          FM_LEG_TBL(:,1) = TO_FM('1.0')
        ELSE
          DO XX = 1, SIZE(X)
            FM_LEG_TBL(XX,1) = (TO_FM('-1.0')**M)*SQRT(TO_FM('1.0')-TO_FM(X(XX))**2)**M
          ENDDO
        ENDIF

        DO XX = 1, SIZE(X)
          FM_LEG_TBL(XX,2) = X(XX)*FM_LEG_TBL(XX,1) 
        ENDDO

        DO NN = 3, NE
          N = M + NN - 1
          DO XX = 1, SIZE(X)
            FM_LEG_TBL(XX,NN) = TO_FM('1.0')/(N-M)*(       & 
            FM_LEG_TBL(XX,NN-1) * TO_FM(X(XX)) -           & 
            FM_LEG_TBL(XX,NN-2) * (N+M-1)/(2*N-1)/(2*N-3))
          ENDDO
        ENDDO

        IF ((SIZE(LOG_NORM, 1) .LT. NE) .OR. (SIZE(LOG_NORM, 2) .LT. ME)) THEN
          WRITE(*,*) 'LEG_TBL: NORM TABLE SIZE INCOMPATIBLE -- TOO SMALL.'
          STOP
        ENDIF

        DO NN = 1, NE
          FM_LEG_TBL(:,NN) = FM_LEG_TBL(:,NN)*EXP(TO_FM(LOG_FACT(M+NN-1))+TO_FM(LOG_NORM(NN,MM)))
        ENDDO

        DO NN = 1, NE
          DO XX = 1, SIZE(X)
            LEG_TBL(XX,NN,MM) = TO_DP(FM_LEG_TBL(XX,NN)) 
          ENDDO
        ENDDO

        FM_LEG_TBL(:,:) = TO_FM(0.D0)
      ENDDO

      CALL FM_DEALLOCATE( FM_LEG_TBL )
      DEALLOCATE( FM_LEG_TBL )

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION LEG_TBL_2D(X, NE, M, LOG_NORM) RESULT(LEG_TBL)
! ==========================================================================================================
! [USAGE]:
! MAKE A 2-DIMENSIONAL TABLE OF THE ASSOCIATED LEGENDRE FUNCTIONS WITH A FIXED AZIMUTH. WAVE NO. M
! [VARIABLES]:
! X >> MULTIPLE POINTS (EXAMPLE. COLLOCATION POINTS) TO EVAULATE LEG. FUN. VALS.
! NE >> MAX. DEGREE OF THE ASSOCIATED LEGENDRE FUNCTION
! M >> THE ORDER OF THE ASSOCIATED LEGENDRE FUNCTION IN USE
! LOG_NORM >> TABLE OF THE LOG NORMALIZATION FACTORS, FROM LEG_LOGNORM()
! LEG_TBL >> THE LEGENDRE FUNCTION COEFFICIENTS AT X FROM N = 1 TO NE. FOR GIVEN X, THE COLUMN IS
!            LEG_TBL(X,:) = (/ P^M_|M|(X), ..., P^M_|M|+NE-1(X) /)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:), INTENT(IN) :: X
      INTEGER, INTENT(IN)                :: NE
      INTEGER, INTENT(IN)                :: M
      REAL(P8), DIMENSION(:), INTENT(IN) :: LOG_NORM
      REAL(P8), DIMENSION(SIZE(X), NE)   :: LEG_TBL

      REAL(P8), DIMENSION(NE, 1) :: LNT

      LNT(:,1) = LOG_NORM(:NE)
      LEG_TBL = RESHAPE(LEG_TBL_3D(X, NE, (/ ABS(M) /), LNT), (/ SIZE(X), NE /))

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION LEG_TBL_1D(X, NE, M, LOG_NORM) RESULT(LEG_TBL)
! ==========================================================================================================
! [USAGE]:
! MAKE A 1-DIMENSIONAL TABLE OF THE ASSOCIATED LEGENDRE FUNCTIONS WITH FIXED M AND X
! [VARIABLES]:
! X >> A POINT TO EVAULATE LEG. FUN. VAL.
! NE >> MAX. DEGREE OF THE ASSOCIATED LEGENDRE FUNCTION
! M >> THE ORDER OF THE ASSOCIATED LEGENDRE FUNCTION IN USE
! LOG_NORM >> TABLE OF THE LOG NORMALIZATION FACTORS, FROM LEG_LOGNORM()
! LEG_TBL >> THE LEGENDRE FUNCTION COEFFICIENTS AT X FROM N = 1 TO NE. THE OUTPUT IS
!            LEG_TBL(:) = (/ P^M_|M|(X), ..., P^M_|M|+NE-1(X) /)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), INTENT(IN)               :: X
      INTEGER, INTENT(IN)                :: NE
      INTEGER, INTENT(IN)                :: M
      REAL(P8), DIMENSION(:), INTENT(IN) :: LOG_NORM
      REAL(P8), DIMENSION(NE)            :: LEG_TBL

      REAL(P8), DIMENSION(NE, 1) :: LNT

      LNT(:,1) = LOG_NORM(:NE)
      LEG_TBL = RESHAPE(LEG_TBL_3D( (/ X /), NE, (/ ABS(M) /), LNT), (/ NE /))

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION LEG_XM_GEN(M, LI, LO, NORM) RESULT(LEG_XM)
! ==========================================================================================================
! [USAGE]:
! MATRIX OPERATOR TO MULTIPLY (1-X) IN FUNCTION SPACE
! [VARIABLES]:
! M >> THE ORDER OF LEGENDRE FUNCTION TO OPERATOR ON
! LI >> LENGTH OF (A^M_N) (N = |M|, |M|+1, ... , |M|+LI-1)
! LO >> LENGTH OF (B^M_N) (N = |M|, |M|+1, ... , |B|+LO-1)
! NORM >> OPTIONAL (FROM LEG_NORM). IF EXIST, NORMALIZED OPERATOR IS OBTAINED
! LEG_XM >> MATRIX EXPRESSION OF THE OPERATOR MULTIPLYING (1-X) IN FUNCTION SPACE
! [NOTES]:
! 1. FOR A FUNCTION F(X) EXPENDED WITH COEFFCIEINT A^M_N AS
! F(X) = SUM_N=|M|^INFTY A^M_N P^M_N(X)
! COEFFICIENTS OF THE FUNCTION (1-X)F(X) B^M_N AS
! (1-X)F(X) = SUM_N=|M|^INFTY B^M_N P^M_N(X)
! CAN BE EVAULATED BY "B^M_N = LEG_XM * A^M_N"
! 2. (MATSUSHIMA & MARCUS, 1997)
! GIVEN THE ORDER M,
! B^M_N = -(N-ABS(M))/(2*N-1)A^M_N-1+A^M_N-(N+ABS(M)+1)/(2*N+3)A^M_N+1
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)                   :: M, LI, LO
      REAL(P8), DIMENSION(:), OPTIONAL      :: NORM
      REAL(P8), DIMENSION(:,:), ALLOCATABLE :: LEG_XM

      INTEGER                               :: NN, N, I, J, AM

      AM = ABS(M)

      ALLOCATE( LEG_XM(LO, LI) )
      LEG_XM = 0.D0

      DO NN = 1, LO
        N = AM + (NN - 1)
        IF (NN-1 .GE. 1)  LEG_XM(NN,NN-1) = -(N-AM)/(2*N-1.D0) ! SUBDIAGONAL
        IF (NN .LE. LI)   LEG_XM(NN,NN)   = 1.0D0 ! MAIN DIAGONAL
        IF (NN+1 .LE. LI) LEG_XM(NN,NN+1) = -(N+AM+1.D0)/(2*N+3.D0) ! SUPERDIAGONAL
      ENDDO

      IF (PRESENT(NORM)) THEN
        IF (SIZE(NORM) .LT. MAX(LO, LI)) THEN
          WRITE(*,*) 'LEG_XM: SIZE(NORM) TOO SMALL'
          WRITE(*,*) 'SIZE(NORM)=',SIZE(NORM)
          WRITE(*,*) 'NI=',LO,'  NJ=',LI
        STOP
        ENDIF

        DO I = 1, LO
          DO J = 1, LI
            IF (LEG_XM(I,J) .NE. 0.D0) LEG_XM(I,J) = LEG_XM(I,J)*(NORM(J)/NORM(I))
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION LEG_XM_BND(M, LI, NORM) RESULT(LEG_XM)
! ==========================================================================================================
! [USAGE]:
! BAND MATRIX OPERATOR TO MULTIPLY (1-X) IN FUNCTION SPACE
! [VARIABLES]:
! M >> THE ORDER OF LEGENDRE FUNCTION TO OPERATOR ON
! LI >> LENGTH OF (A^M_N) (N = |M|, |M|+1, ... , |M|+LI-1)
! LO >> LENGTH OF (B^M_N) (N = |M|, |M|+1, ... , |B|+LO-1)
! NORM >> OPTIONAL (FROM LEG_NORM). IF EXIST, NORMALIZED OPERATOR IS OBTAINED
! LEG_XM >> BAND MATRIX EXPRESSION OF THE OPERATOR MULTIPLYING (1-X) IN FUNCTION SPACE
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)              :: M, LI
      REAL(P8), DIMENSION(:), OPTIONAL :: NORM
      TYPE(REAL_BNDM)                  :: LEG_XM

      IF (PRESENT(NORM)) THEN
        LEG_XM = BANDMAT(LEG_XM_GEN(M, LI, LI, NORM), KL=1, KU=1)
      ELSE
        LEG_XM = BANDMAT(LEG_XM_GEN(M, LI, LI), KL=1, KU=1)
      ENDIF

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION LEG_XP_GEN(M, LI, LO, NORM) RESULT(LEG_XP)
! ==========================================================================================================
! [USAGE]:
! MATRIX OPERATOR TO MULTIPLY (1+X) IN FUNCTION SPACE
! [VARIABLES]:
! M >> THE ORDER OF LEGENDRE FUNCTION TO OPERATOR ON
! LI >> LENGTH OF (A^M_N) (N = |M|, |M|+1, ... , |M|+LI-1)
! LO >> LENGTH OF (B^M_N) (N = |M|, |M|+1, ... , |B|+LO-1)
! NORM >> OPTIONAL (FROM LEG_NORM). IF EXIST, NORMALIZED OPERATOR IS OBTAINED
! LEG_XP >> MATRIX EXPRESSION OF THE OPERATOR MULTIPLYING (1+X) IN FUNCTION SPACE
! [NOTES]:
! 1. FOR A FUNCTION F(X) EXPENDED WITH COEFFCIEINT A^M_N AS
! F(X) = SUM_N=|M|^INFTY A^M_N P^M_N(X)
! COEFFICIENTS OF THE FUNCTION (1-X)F(X) B^M_N AS
! (1-X)F(X) = SUM_N=|M|^INFTY B^M_N P^M_N(X)
! CAN BE EVAULATED BY "B^M_N = LEG_XP * A^M_N"
! 2. (MATSUSHIMA & MARCUS, 1997)
! GIVEN THE ORDER M,
! B^M_N = (N-ABS(M))/(2*N-1)A^M_N-1+A^M_N+(N+ABS(M)+1)/(2*N+3)A^M_N+1
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)                   :: M, LI, LO
      REAL(P8), DIMENSION(:), OPTIONAL      :: NORM
      REAL(P8), DIMENSION(:,:), ALLOCATABLE :: LEG_XP

      INTEGER                               :: NN, N, I, J, AM

      AM = ABS(M)

      ALLOCATE( LEG_XP(LO, LI) )
      LEG_XP = 0.D0

      DO NN = 1, LO
        N = AM + (NN - 1)
        IF (NN-1 .GE. 1)  LEG_XP(NN,NN-1) = (N-AM)/(2*N-1.D0) ! SUBDIAGONAL
        IF (NN .LE. LI)   LEG_XP(NN,NN)   = 1.0D0 ! MAIN DIAGONAL
        IF (NN+1 .LE. LI) LEG_XP(NN,NN+1) = (N+AM+1.D0)/(2*N+3.D0) ! SUPERDIAGONAL
      ENDDO

      IF (PRESENT(NORM)) THEN
        IF (SIZE(NORM) .LT. MAX(LO, LI)) THEN
          WRITE(*,*) 'LEG_XP: SIZE(NORM) TOO SMALL'
          WRITE(*,*) 'SIZE(NORM)=',SIZE(NORM)
          WRITE(*,*) 'NI=',LO,'  NJ=',LI
        STOP
        ENDIF

        DO I = 1, LO
          DO J = 1, LI
            IF (LEG_XP(I,J) .NE. 0.D0) LEG_XP(I,J) = LEG_XP(I,J)*(NORM(J)/NORM(I))
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION LEG_XP_BND(M, LI, NORM) RESULT(LEG_XP)
! ==========================================================================================================
! [USAGE]:
! BAND MATRIX OPERATOR TO MULTIPLY (1+X) IN FUNCTION SPACE
! [VARIABLES]:
! M >> THE ORDER OF LEGENDRE FUNCTION TO OPERATOR ON
! LI >> LENGTH OF (A^M_N) (N = |M|, |M|+1, ... , |M|+LI-1)
! LO >> LENGTH OF (B^M_N) (N = |M|, |M|+1, ... , |B|+LO-1)
! NORM >> OPTIONAL (FROM LEG_NORM). IF EXIST, NORMALIZED OPERATOR IS OBTAINED
! LEG_XP >> BAND MATRIX EXPRESSION OF THE OPERATOR MULTIPLYING (1+X) IN FUNCTION SPACE
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)              :: M, LI
      REAL(P8), DIMENSION(:), OPTIONAL :: NORM
      TYPE(REAL_BNDM)                  :: LEG_XP

      IF (PRESENT(NORM)) THEN
        LEG_XP = BANDMAT(LEG_XP_GEN(M, LI, LI, NORM), KL=1, KU=1)
      ELSE
        LEG_XP = BANDMAT(LEG_XP_GEN(M, LI, LI), KL=1, KU=1)
      ENDIF

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION LEG_X_GEN(M, LI, LO, NORM) RESULT(LEG_X)
! ==========================================================================================================
! [USAGE]:
! MATRIX OPERATOR TO MULTIPLY X IN FUNCTION SPACE
! [VARIABLES]:
! M >> THE ORDER OF LEGENDRE FUNCTION TO OPERATOR ON
! LI >> LENGTH OF (A^M_N) (N = |M|, |M|+1, ... , |M|+LI-1)
! LO >> LENGTH OF (B^M_N) (N = |M|, |M|+1, ... , |B|+LO-1)
! NORM >> OPTIONAL (FROM LEG_NORM). IF EXIST, NORMALIZED OPERATOR IS OBTAINED
! LEG_X >> MATRIX EXPRESSION OF THE OPERATOR MULTIPLYING X IN FUNCTION SPACE
! [NOTES]:
! 1. FOR A FUNCTION F(X) EXPENDED WITH COEFFCIEINT A^M_N AS
! F(X) = SUM_N=|M|^INFTY A^M_N P^M_N(X)
! COEFFICIENTS OF THE FUNCTION (1-X)F(X) B^M_N AS
! X*F(X) = SUM_N=|M|^INFTY B^M_N P^M_N(X)
! CAN BE EVAULATED BY "B^M_N = LEG_X * A^M_N"
! 2. (MATSUSHIMA & MARCUS, 1997)
! GIVEN THE ORDER M,
! B^M_N = (N-ABS(M))/(2*N-1)A^M_N-1+(N+ABS(M)+1)/(2*N+3)A^M_N+1
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)                   :: M, LI, LO
      REAL(P8), DIMENSION(:), OPTIONAL      :: NORM
      REAL(P8), DIMENSION(:,:), ALLOCATABLE :: LEG_X

      INTEGER                               :: NN, N, I, J, AM

      AM = ABS(M)

      ALLOCATE( LEG_X(LO, LI) )
      LEG_X = 0.D0

      DO NN = 1, LO
        N = AM + (NN - 1)
        IF (NN-1 .GE. 1)  LEG_X(NN,NN-1) = (N-AM)/(2*N-1.D0) ! SUBDIAGONAL
        IF (NN+1 .LE. LI) LEG_X(NN,NN+1) = (N+AM+1.D0)/(2*N+3.D0) ! SUPERDIAGONAL
      ENDDO

      IF (PRESENT(NORM)) THEN
        IF (SIZE(NORM) .LT. MAX(LO, LI)) THEN
          WRITE(*,*) 'LEG_X: SIZE(NORM) TOO SMALL'
          WRITE(*,*) 'SIZE(NORM)=',SIZE(NORM)
          WRITE(*,*) 'NI=',LO,'  NJ=',LI
        STOP
        ENDIF

        DO I = 1, LO
          DO J = 1, LI
            IF (LEG_X(I,J) .NE. 0.D0) LEG_X(I,J) = LEG_X(I,J)*(NORM(J)/NORM(I))
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION LEG_X_BND(M, LI, NORM) RESULT(LEG_X)
! ==========================================================================================================
! [USAGE]:
! BAND MATRIX OPERATOR TO MULTIPLY X IN FUNCTION SPACE
! [VARIABLES]:
! M >> THE ORDER OF LEGENDRE FUNCTION TO OPERATOR ON
! LI >> LENGTH OF (A^M_N) (N = |M|, |M|+1, ... , |M|+LI-1)
! LO >> LENGTH OF (B^M_N) (N = |M|, |M|+1, ... , |B|+LO-1)
! NORM >> OPTIONAL (FROM LEG_NORM). IF EXIST, NORMALIZED OPERATOR IS OBTAINED
! LEG_X >> BAND MATRIX EXPRESSION OF THE OPERATOR MULTIPLYING X IN FUNCTION SPACE
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)              :: M, LI
      REAL(P8), DIMENSION(:), OPTIONAL :: NORM
      TYPE(REAL_BNDM)                  :: LEG_X

      IF (PRESENT(NORM)) THEN
        LEG_X = BANDMAT(LEG_X_GEN(M, LI, LI, NORM), KL=1, KU=1)
      ELSE
        LEG_X = BANDMAT(LEG_X_GEN(M, LI, LI), KL=1, KU=1)
      ENDIF

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION LEG_XXDX_GEN(M, LI, LO, NORM) RESULT(LEG_XXDX)
! ==========================================================================================================
! [USAGE]:
! MATRIX OPERATOR TO MULTIPLY (1-X**2)*D/DX() IN FUNCTION SPACE
! [VARIABLES]:
! M >> THE ORDER OF LEGENDRE FUNCTION TO OPERATOR ON
! LI >> LENGTH OF (A^M_N) (N = |M|, |M|+1, ... , |M|+LI-1)
! LO >> LENGTH OF (B^M_N) (N = |M|, |M|+1, ... , |B|+LO-1)
! NORM >> OPTIONAL (FROM LEG_NORM). IF EXIST, NORMALIZED OPERATOR IS OBTAINED
! LEG_XXDX >> MATRIX EXPRESSION OF THE OPERATOR MULTIPLYING (1-X**2)*D/DX() IN FUNCTION SPACE
! [NOTES]:
! 1. FOR A FUNCTION F(X) EXPENDED WITH COEFFCIEINT A^M_N AS
! F(X) = SUM_N=|M|^INFTY A^M_N P^M_N(X)
! COEFFICIENTS OF THE FUNCTION (1-X)F(X) B^M_N AS
! X*F(X) = SUM_N=|M|^INFTY B^M_N P^M_N(X)
! CAN BE EVAULATED BY "B^M_N = LEG_XXDX * A^M_N"
! 2. (MATSUSHIMA & MARCUS, 1997)
! GIVEN THE ORDER M,
! B^M_N =(N-1)(N-ABS(M))/(2*N-1)A^M_N-1+(N+2)(N+ABS(M)+1)/(2*N+3)A^M_N+1
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)                   :: M, LI, LO
      REAL(P8), DIMENSION(:), OPTIONAL      :: NORM
      REAL(P8), DIMENSION(:,:), ALLOCATABLE :: LEG_XXDX

      INTEGER                               :: NN, N, I, J, AM

      AM = ABS(M)

      ALLOCATE( LEG_XXDX(LO, LI) )
      LEG_XXDX = 0.D0

      DO NN = 1, LO
        N = AM + (NN - 1)
        IF (NN-1 .GE. 1)  LEG_XXDX(NN,NN-1) = -(N-1.D0)*(N-AM)/(2*N-1.D0) ! SUBDIAGONAL
        IF (NN+1 .LE. LI) LEG_XXDX(NN,NN+1) = (N+2.D0)*(N+AM+1.D0)/(2*N+3.D0) ! SUPERDIAGONAL
      ENDDO

      IF (PRESENT(NORM)) THEN
        IF (SIZE(NORM) .LT. MAX(LO, LI)) THEN
          WRITE(*,*) 'LEG_XXDX: SIZE(NORM) TOO SMALL'
          WRITE(*,*) 'SIZE(NORM)=',SIZE(NORM)
          WRITE(*,*) 'NI=',LO,'  NJ=',LI
        STOP
        ENDIF

        DO I = 1, LO
          DO J = 1, LI
            IF (LEG_XXDX(I,J) .NE. 0.D0) LEG_XXDX(I,J) = LEG_XXDX(I,J)*(NORM(J)/NORM(I))
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION LEG_XXDX_BND(M, LI, NORM) RESULT(LEG_XXDX)
! ==========================================================================================================
! [USAGE]:
! BAND MATRIX OPERATOR TO MULTIPLY (1-X**2)*D/DX() IN FUNCTION SPACE
! [VARIABLES]:
! M >> THE ORDER OF LEGENDRE FUNCTION TO OPERATOR ON
! LI >> LENGTH OF (A^M_N) (N = |M|, |M|+1, ... , |M|+LI-1)
! LO >> LENGTH OF (B^M_N) (N = |M|, |M|+1, ... , |B|+LO-1)
! NORM >> OPTIONAL (FROM LEG_NORM). IF EXIST, NORMALIZED OPERATOR IS OBTAINED
! LEG_XXDX >> BAND MATRIX EXPRESSION OF THE OPERATOR MULTIPLYING (1-X**2)*D/DX() IN FUNCTION SPACE
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)              :: M, LI
      REAL(P8), DIMENSION(:), OPTIONAL :: NORM
      TYPE(REAL_BNDM)                  :: LEG_XXDX

      IF (PRESENT(NORM)) THEN
        LEG_XXDX = BANDMAT(LEG_XXDX_GEN(M, LI, LI, NORM), KL=1, KU=1)
      ELSE
        LEG_XXDX = BANDMAT(LEG_XXDX_GEN(M, LI, LI), KL=1, KU=1)
      ENDIF

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION LEG_DEL2H_GEN(M, LI, LO, ELL, NORM) RESULT(LEG_DEL2H)
! ==========================================================================================================
! [USAGE]:
! MATRIX OPERATOR TO MULTIPLY DEL^2_PERP (HOR. LAPLACIAN) IN FUNCTION SPACE
! [VARIABLES]:
! M >> THE ORDER OF LEGENDRE FUNCTION TO OPERATOR ON
! LI >> LENGTH OF (A^M_N) (N = |M|, |M|+1, ... , |M|+LI-1)
! LO >> LENGTH OF (B^M_N) (N = |M|, |M|+1, ... , |B|+LO-1)
! NORM >> OPTIONAL (FROM LEG_NORM). IF EXIST, NORMALIZED OPERATOR IS OBTAINED
! LEG_DEL2H >> MATRIX EXPRESSION OF THE OPERATOR MULTIPLYING DEL^2_PERP IN FUNCTION SPACE
! [NOTES]:
! 1. FOR A FUNCTION F(X) EXPENDED WITH COEFFCIEINT A^M_N AS
! F(X(R),PHI,Z) = SUM_N=|M|^INFTY A^M_N P^M_N(X(R))*EXP(I*M*PHI+I*AK*Z)
! COEFFICIENTS OF THE FUNCTION (1-X)F(X) B^M_N AS
! DEL^2_PERP(F(X(R))) = SUM_N=|M|^INFTY B^M_N P^M_N(X(R), PHI, Z)
!                                                  *EXP(I*M*PHI+I*AK*Z)
! CAN BE EVAULATED BY "B^M_N = LEG_DEL2H * A^M_N"
! 2. TO FIND THE RECURRENCE RELATIONSHIP,
! SEE EQUATION (107) IN THE APPENDIX OF MATSUSHIMA & MARCUS (1997)
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)                   :: M, LI, LO
      REAL(P8), INTENT(IN)                  :: ELL
      REAL(P8), DIMENSION(:), OPTIONAL      :: NORM
      REAL(P8), DIMENSION(:,:), ALLOCATABLE :: LEG_DEL2H

      INTEGER                               :: NN, N, I, J, AM
      REAL(P8)                              :: ELL2
      REAL(P8), PARAMETER                   :: S = 0.D0

      AM = ABS(M)
      ELL2 = ELL**2.D0

      ALLOCATE( LEG_DEL2H(LO, LI) )
      LEG_DEL2H = 0.D0

      DO NN = 1, LO
        N = AM + (NN - 1)
        IF (NN-2 .GE. 1)  LEG_DEL2H(NN,NN-2) = -(N-AM-1.D0)*(N-AM)*(N-2+S)*(N-1+S)/(2*N-3.)/(2*N-1.D0)/ELL2
        IF (NN-1 .GE. 1)  LEG_DEL2H(NN,NN-1) = 2.D0*N*(N-AM)*(N-1+S)/(2*N-1.D0) /ELL2
        IF (NN .LE. LI)   LEG_DEL2H(NN,NN)   = (-2.D0*N*(N+1.D0)*(3.D0*N*N+3*N-AM**2-2)&
                                               +2*S*(S-2)**(N*N+N+M*M-1.D0))/(2*N-1.D0)/(2*N+3.D0)/ELL2  
        IF (NN+1 .LE. LI) LEG_DEL2H(NN,NN+1) = 2.D0*(N+1)*(N+AM+1.D0)*(N+2-S)/(2*N+3.D0)/ELL2
        IF (NN+2 .LE. LI) LEG_DEL2H(NN,NN+2) = -(N+AM+1.D0)*(N+AM+2.D0)*(N+3-S)*(N+2-S)&
                                               /(2*N+3.D0)/(2*N+5.D0)/ELL2
      ENDDO

      IF (PRESENT(NORM)) THEN
        IF (SIZE(NORM) .LT. MAX(LO, LI)) THEN
          WRITE(*,*) 'LEG_DEL2H: SIZE(NORM) TOO SMALL'
          WRITE(*,*) 'SIZE(NORM)=',SIZE(NORM)
          WRITE(*,*) 'NI=',LO,'  NJ=',LI
        STOP
        ENDIF

        DO I = 1, LO
          DO J = 1, LI
            IF (LEG_DEL2H(I,J) .NE. 0.D0) LEG_DEL2H(I,J) = LEG_DEL2H(I,J)*(NORM(J)/NORM(I))
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION LEG_DEL2H_BND(M, LI, ELL, NORM) RESULT(LEG_DEL2H)
! ==========================================================================================================
! [USAGE]:
! BAND MATRIX OPERATOR TO MULTIPLY DEL^2_PERP (HOR. LAPLACIAN) IN FUNCTION SPACE
! [VARIABLES]:
! M >> THE ORDER OF LEGENDRE FUNCTION TO OPERATOR ON
! LI >> LENGTH OF (A^M_N) (N = |M|, |M|+1, ... , |M|+LI-1)
! LO >> LENGTH OF (B^M_N) (N = |M|, |M|+1, ... , |B|+LO-1)
! NORM >> OPTIONAL (FROM LEG_NORM). IF EXIST, NORMALIZED OPERATOR IS OBTAINED
! LEG_DEL2H >> BAND MATRIX EXPRESSION OF THE OPERATOR MULTIPLYING DEL^2_PERP IN FUNCTION SPACE
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)              :: M, LI
      REAL(P8), INTENT(IN)             :: ELL
      REAL(P8), DIMENSION(:), OPTIONAL :: NORM
      TYPE(REAL_BNDM)                  :: LEG_DEL2H

      IF (PRESENT(NORM)) THEN
        LEG_DEL2H = BANDMAT(LEG_DEL2H_GEN(M, LI, LI, ELL, NORM), KL=2, KU=2)
      ELSE
        LEG_DEL2H = BANDMAT(LEG_DEL2H_GEN(M, LI, LI, ELL), KL=2, KU=2)
      ENDIF

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION LEG_DEL2_GEN(M, AK, LI, LO, ELL, NORM) RESULT(LEG_DEL2)
! ==========================================================================================================
! [USAGE]:
! MATRIX OPERATOR TO MULTIPLY (DEL^2_PERP - AK**2) IN FUNCTION SPACE
! [VARIABLES]:
! M >> THE ORDER OF LEGENDRE FUNCTION TO OPERATOR ON
! AK >> REAL VALUE (INTENDED TO REPRESENT THE AXIAL WAVELENGTH, AK, IN USE. SEE NOTES)
! LI >> LENGTH OF (A^M_N) (N = |M|, |M|+1, ... , |M|+LI-1)
! LO >> LENGTH OF (B^M_N) (N = |M|, |M|+1, ... , |B|+LO-1)
! NORM >> OPTIONAL (FROM LEG_NORM). IF EXIST, NORMALIZED OPERATOR IS OBTAINED
! LEG_DEL2 >> MATRIX EXPRESSION OF THE OPERATOR MULTIPLYING (DEL^2_PERP - AK**2) IN FUNCTION SPACE
! [NOTES]:
! 1. FOR A FUNCTION F(X) EXPENDED WITH COEFFCIEINT A^M_N AS
! F(X(R),PHI,Z) = SUM_N=|M|^INFTY A^M_N P^M_N(X(R))*EXP(I*M*PHI+I*AK*Z)
! COEFFICIENTS OF THE FUNCTION (1-X)F(X) B^M_N AS
! DEL^2 (F(X(R))) = SUM_N=|M|^INFTY B^M_N P^M_N(X(R), PHI, Z)
!                                                  *EXP(I*M*PHI+I*AK*Z)
! CAN BE EVAULATED BY "B^M_N = LEG_DEL2 * A^M_N"
! 2. TO FIND THE RECURRENCE RELATIONSHIP,
! SEE EQUATION (107) IN THE APPENDIX OF MATSUSHIMA & MARCUS (1997)
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)                   :: M, LI, LO
      REAL(P8), INTENT(IN)                  :: AK, ELL
      REAL(P8), DIMENSION(:), OPTIONAL      :: NORM
      REAL(P8), DIMENSION(:,:), ALLOCATABLE :: LEG_DEL2

      INTEGER                               :: NN, N, I, J, AM
      REAL(P8)                              :: ELL2
      REAL(P8), PARAMETER                   :: S = 0.D0

      AM = ABS(M)
      ELL2 = ELL**2.D0

      ALLOCATE( LEG_DEL2(LO, LI) )
      LEG_DEL2 = 0.D0

      DO NN = 1, LO
        N = AM + (NN - 1)
        IF (NN-2 .GE. 1)  LEG_DEL2(NN,NN-2) = -(N-AM-1.D0)*(N-AM)*(N-2+S)*(N-1+S)/(2*N-3.)/(2*N-1.D0)/ELL2
        IF (NN-1 .GE. 1)  LEG_DEL2(NN,NN-1) = 2.D0*N*(N-AM)*(N-1+S)/(2*N-1.D0) /ELL2
        IF (NN .LE. LI)   LEG_DEL2(NN,NN)   = (-2.D0*N*(N+1.D0)*(3.D0*N*N+3*N-AM**2-2)&
                                               +2*S*(S-2)**(N*N+N+M*M-1.D0))/(2*N-1.D0)/(2*N+3.D0)/ELL2  
        IF (NN+1 .LE. LI) LEG_DEL2(NN,NN+1) = 2.D0*(N+1)*(N+AM+1.D0)*(N+2-S)/(2*N+3.D0)/ELL2
        IF (NN+2 .LE. LI) LEG_DEL2(NN,NN+2) = -(N+AM+1.D0)*(N+AM+2.D0)*(N+3-S)*(N+2-S)&
                                               /(2*N+3.D0)/(2*N+5.D0)/ELL2
      ENDDO

      DO NN = 1, LO
        IF (NN .LE. LI)   LEG_DEL2(NN,NN)   = LEG_DEL2(NN,NN) - AK**2.
      ENDDO

      IF (PRESENT(NORM)) THEN
        IF (SIZE(NORM) .LT. MAX(LO, LI)) THEN
          WRITE(*,*) 'LEG_DEL2: SIZE(NORM) TOO SMALL'
          WRITE(*,*) 'SIZE(NORM)=',SIZE(NORM)
          WRITE(*,*) 'NI=',LO,'  NJ=',LI
        STOP
        ENDIF

        DO I = 1, LO
          DO J = 1, LI
            IF (LEG_DEL2(I,J) .NE. 0.D0) LEG_DEL2(I,J) = LEG_DEL2(I,J)*(NORM(J)/NORM(I))
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION LEG_DEL2_BND(M, AK, LI, ELL, NORM) RESULT(LEG_DEL2)
! ==========================================================================================================
! [USAGE]:
! BAND MATRIX OPERATOR TO MULTIPLY (DEL^2_PERP - AK**2) IN FUNCTION SPACE
! [VARIABLES]:
! M >> THE ORDER OF LEGENDRE FUNCTION TO OPERATOR ON
! AK >> REAL VALUE REPRESENTING THE AXIAL WAVELENGTH IN USE
! LI >> LENGTH OF (A^M_N) (N = |M|, |M|+1, ... , |M|+LI-1)
! LO >> LENGTH OF (B^M_N) (N = |M|, |M|+1, ... , |B|+LO-1)
! NORM >> OPTIONAL (FROM LEG_NORM). IF EXIST, NORMALIZED OPERATOR IS OBTAINED
! LEG_DEL2H >> BAND MATRIX EXPRESSION OF THE OPERATOR MULTIPLYING (DEL^2_PERP - AK**2) IN FUNCTION SPACE
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)              :: M, LI
      REAL(P8), INTENT(IN)             :: AK, ELL
      REAL(P8), DIMENSION(:), OPTIONAL :: NORM
      TYPE(REAL_BNDM)                  :: LEG_DEL2

      IF (PRESENT(NORM)) THEN
        LEG_DEL2 = BANDMAT(LEG_DEL2_GEN(M, AK, LI, LI, ELL, NORM), KL=2, KU=2)
      ELSE
        LEG_DEL2 = BANDMAT(LEG_DEL2_GEN(M, AK, LI, LI, ELL), KL=2, KU=2)
      ENDIF

      RETURN
      END FUNCTION
! ==========================================================================================================
!  PRIVATE PROCEDURES ======================================================================================
! ==========================================================================================================
      SUBROUTINE FIND_GAULEG_ZEROS(X, W, R, ELL)
! ==========================================================================================================
! [USAGE]:
! FIND THE ZEROS AND ABSCISSA OF LEGENDRE FUNCTIONS P^M_N(X), OR P^M_L_N(R)
! HERE P^M_L_N(R) == P^M_N([R^2 - ELL^2]/[R^2 + ELL^2])
! [VARIABLES]:
! X >> ABSCISSA (OR COLLOCATION POINTS) OF P^M_N(X) WHERE -1 < X < 1.
! R >> COLLOCATION POINTS AFTER MAPPING X = (R^2 - ELL^2) / (R^2 + ELL^2).
! ELL >> MAP PARAMETER
! W >> GAUSS-LEGENDRE WEIGHTS FOR P^M_N(X) = P^M_L_N(R) FOR EACH ZERO
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:), INTENT(INOUT) :: X
      REAL(P8), DIMENSION(:), INTENT(INOUT) :: W
      REAL(P8), DIMENSION(:), INTENT(INOUT) :: R
      REAL(P8), INTENT(IN)                  :: ELL

      CALL GAUSS_LEGENDRE(-1.D0, 1.D0, X, W, SIZE(X))

      R = ELL * SQRT((1+X)/(1-X))

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE GAUSS_LEGENDRE(X1, X2, X, W, N)
! ==========================================================================================================
! [USAGE]: 
! GET THE ABSCISSAS AND WEIGHTS OF GAUSS-LEGENDRE N-POINT QUADRATURES
! FROM NUMERICAL RECIPES IN FORTRAN90 (2ND ED.) (REF. Press et al., 1996)
! [VARIABLES]:
! X1,X2 >> LOWER AND UPPER BOUNDS
! X >> ABSCISSAS (OR COLLOCATION POINTS) FOR P^M_N
! W >> WEIGHTS FOR P^M_N FOR EACH ZERO (* P^M_N: GENERAL LEGENDRE POLY.)
! N >> # OF QUADRATURES TO BE CALCULATED
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8) , INTENT(IN)   :: X1,X2 ! LOWER & UPPER BOUNDS
      REAL(P8) , DIMENSION(:) :: X,W   ! VALUE & WEIGHT (1:N)
      INTEGER                 :: N     ! TOTAL NUMBER OF QUATRATURES

      REAL(P8) , PARAMETER    :: EPS=1.D-15 ! ERROR LIMIT
      REAL(P8)                :: XM, XL, P1, P2, P3, Z, Z1, PP
      INTEGER                 :: M, I, J

      M = CEILING((N+1)/2.D0)
      XM = (X2+X1)/2.D0
      XL = (X2-X1)/2.D0

      DO I = 1, M ! REFINEMENT USING NEWTON'S METHOD
        Z = COS(PI*(I-.25D0)/(N+.5D0))
        DO WHILE (ABS(Z-Z1) .GT. EPS)
          P1 = 1.D0
          P2 = 0.D0
          DO J = 1,N
            P3 = P2
            P2 = P1
            P1 = ((2*J-1)*Z*P2-(J-1)*P3)/J
          ENDDO
          PP = N*(Z*P1-P2)/(Z*Z-1)
          Z1 = Z
          Z  = Z1-P1/PP
        ENDDO
        X(I)     = XM-XL*Z
        X(N+1-I) = XM+XL*Z
        W(I)     = 2*XL/((1-Z*Z)*PP*PP)
        W(N+1-I) = W(I)
      ENDDO

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      FUNCTION LOG_FACT(M) RESULT(LF)
! ==========================================================================================================
! [USAGE]:
! COMPUTE LOG((2M)! / (2^M * M!)) OR EQUIVALENTLY LOG((2M-1)!!)
! [VARIABLES]:
! M >> FACTORIAL ARGUMENT
! LF >> LOG((2M)! / (2^M * M!)) OR EQUIVALENTLY LOG((2M-1)!!)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: M
      REAL(P8)            :: LF
      
      LF = LOG_GAMMA(2*M+1.D0) - M*LOG(2.D0) - LOG_GAMMA(M+1.D0)

      RETURN
      END FUNCTION LOG_FACT
! ==========================================================================================================
END MODULE