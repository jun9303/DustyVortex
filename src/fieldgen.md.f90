MODULE FIELDGEN ! MODULE FOR 3D SCALAR/VECTOR FIELD GENERATION SUBROUTINES & FUNCTIONS
      USE OMP_LIB; USE MPI
      USE MISC; USE MATOPS
      USE LEGENDRE; USE FOURIER
      IMPLICIT NONE
      PRIVATE
! ==========================================================================================================
!  PARAMS/TYPES ============================================================================================
! ==========================================================================================================
      TYPE(LEG_TRANSFORM), PUBLIC :: TFM

      TYPE :: FIELD_INFO
        REAL(P8), DIMENSION(:), ALLOCATABLE :: R, P, PR, PI, Z
        REAL(P8) :: ZLEN, ELL, ELL2, VISC, VISCP
        INTEGER :: NR, NRH, NP, NPH, NZ, VISCPOW
        INTEGER :: NRDIM, NPDIM, NZDIM
        INTEGER :: NRCHOP, NPCHOP, NZCHOP, NZCHOPU
        INTEGER, DIMENSION(:), ALLOCATABLE :: NRCHOPS, NPCHOPS
        INTEGER, DIMENSION(:), ALLOCATABLE :: M
        REAL(P8), DIMENSION(:), ALLOCATABLE :: AK
      END TYPE
      TYPE(FIELD_INFO), PUBLIC :: FINFO

      TYPE, PUBLIC :: SCALAR_FIELD
        COMPLEX(P8), DIMENSION(:,:,:), ALLOCATABLE :: E                           ! SCALAR QTY. IN (R, P, Z)
        INTEGER :: NRC, NPC, NZC, NZCU                            ! NO. OF SPECTRAL ELMTS. IN EACH DIRECTION
        INTEGER, DIMENSION(:), ALLOCATABLE :: NRCS, NPCS             ! CHOP INDICES IN R/P WITH RESPECT TO M
        REAL(P8) :: LN                    ! LOG TERM EXPRESSION. ONLY USED FOR THE TOROIDAL-POLOIDAL DECOMP.
        CHARACTER(LEN=3) :: SPACE  ! IN R-P-Z, RESPECTIVELY, INDICATES WHETHER THE QTY. IS IN PHYS. OR FUNC.
      END TYPE

      TYPE, PUBLIC :: VECTOR_FIELD
        COMPLEX(P8), DIMENSION(:,:,:), ALLOCATABLE :: ER                          ! R-COMP QTY. IN (R, P, Z)
        COMPLEX(P8), DIMENSION(:,:,:), ALLOCATABLE :: EP                          ! P-COMP QTY. IN (R, P, Z)
        COMPLEX(P8), DIMENSION(:,:,:), ALLOCATABLE :: EZ                          ! Z-COMP QTY. IN (R, P, Z)
        INTEGER :: NRC, NPC, NZC, NZCU                            ! NO. OF SPECTRAL ELMTS. IN EACH DIRECTION
        INTEGER, DIMENSION(:), ALLOCATABLE :: NRCS, NPCS             ! CHOP INDICES IN R/P WITH RESPECT TO M
        CHARACTER(LEN=3) :: SPACE  ! IN R-P-Z, RESPECTIVELY, INDICATES WHETHER THE QTY. IS IN PHYS. OR FUNC.
      END TYPE
! ==========================================================================================================
!  DECLARATIONS ============================================================================================
! ==========================================================================================================
      PUBLIC :: FIELD_INIT                          ! GLOBAL (EULERIAN) FIELD SETUP FOR FLUID VELOCITY, ETC.
      PUBLIC :: FIELD_SET                                                 ! SET UP A Q-VORTEX VELOCITY FIELD
      PUBLIC :: ALLOC, DEALLOC                                                        ! ALLOCATE A NEW FIELD
      PUBLIC :: CHOP_SET                                            ! UPDATE THE CHOP LOCATIONS IN THE FIELD
      PUBLIC :: CHOP_RST                                           ! RESET THE CHOP LOCATIONS TO THE DEFAULT
      PUBLIC :: CHOP_DO                               ! PERFORM THE CHOPPING ACCORDING TO THE CHOP VARIABLES
      PUBLIC :: TRANS                                  ! SPACE TRANSFORM 'PPP' <-> 'PFP' <-> 'FFP' <-> 'FFF' 
      PUBLIC :: CALCAT_EXACT, CALCAT_LINTP             ! CALCULATES THE FIELD AT AN ARBITRARY LOC. (R, P, Z)
      PUBLIC :: CALCAT0, CALCAT1                     ! CALCULATES THE FIELD VALUE AT INF.(AT1) OR ORIG.(AT0)
      PUBLIC :: ZEROAT1                   ! ZERO THE FIELD VALUE AT INF.(AT1) BY ADJUSTING (0,0,K)-TH COEFF.
      PUBLIC :: INTEG                                           ! INTEGRATE THE FIELD OVER THE ENTIRE DOMAIN
      PUBLIC :: MSAVE, MLOAD                                               ! SAVE/LOAD THE FIELD INFORMATION
! ==========================================================================================================
!  INTERFACES ==============================================================================================
! ==========================================================================================================
      INTERFACE FIELD_INIT
        MODULE PROCEDURE FIELD_INIT
      END INTERFACE

      INTERFACE FIELD_SET
        MODULE PROCEDURE FIELD_SET_VECTOR
      END INTERFACE

      INTERFACE ALLOC
        MODULE PROCEDURE ALLOC_SCALAR
        MODULE PROCEDURE ALLOC_VECTOR
      END INTERFACE

      INTERFACE DEALLOC
        MODULE PROCEDURE DEALLOC_SCALAR
        MODULE PROCEDURE DEALLOC_VECTOR
      END INTERFACE

      INTERFACE CHOP_SET
        MODULE PROCEDURE CHOPSET_SCALAR
        MODULE PROCEDURE CHOPSET_VECTOR
      END INTERFACE

      INTERFACE CHOP_RST
        MODULE PROCEDURE CHOPRST_SCALAR
        MODULE PROCEDURE CHOPRST_VECTOR
      END INTERFACE

      INTERFACE CHOP_DO
        MODULE PROCEDURE CHOPDO_SCALAR
        MODULE PROCEDURE CHOPDO_VECTOR
      END INTERFACE

      INTERFACE TRANS
        MODULE PROCEDURE TRANS_SCALAR
        MODULE PROCEDURE TRANS_VECTOR
      END INTERFACE

      INTERFACE CALCAT_EXACT
        MODULE PROCEDURE CALCAT_EXACT_SCALAR
        MODULE PROCEDURE CALCAT_EXACT_VECTOR
      END INTERFACE

      INTERFACE CALCAT_LINTP
        MODULE PROCEDURE CALCAT_LINTP_SCALAR_MULTI
        MODULE PROCEDURE CALCAT_LINTP_SCALAR_SINGLE
        MODULE PROCEDURE CALCAT_LINTP_VECTOR_MULTI
        MODULE PROCEDURE CALCAT_LINTP_VECTOR_SINGLE
      END INTERFACE

      INTERFACE CALCAT0
        MODULE PROCEDURE CALCAT0_SCALAR
        MODULE PROCEDURE CALCAT0_VECTOR
      END INTERFACE

      INTERFACE CALCAT1
        MODULE PROCEDURE CALCAT1_SCALAR
        MODULE PROCEDURE CALCAT1_VECTOR
      END INTERFACE

      INTERFACE ZEROAT1
        MODULE PROCEDURE ZEROAT1_SCALAR
        MODULE PROCEDURE ZEROAT1_VECTOR
      END INTERFACE

      INTERFACE INTEG
        MODULE PROCEDURE INTEG_SCALAR
        MODULE PROCEDURE INTEG_VECTOR
      END INTERFACE

      INTERFACE MSAVE
        MODULE PROCEDURE MSAVE_SCALAR
        MODULE PROCEDURE MSAVE_VECTOR
      END INTERFACE  

      INTERFACE MLOAD
        MODULE PROCEDURE MLOAD_SCALAR
        MODULE PROCEDURE MLOAD_VECTOR
      END INTERFACE

CONTAINS
! ==========================================================================================================
!  PUBLIC PROCEDURES =======================================================================================
! ==========================================================================================================
      SUBROUTINE FIELD_INIT(NRIN, NPIN, NZIN, ELLIN, ZLENIN, VISCIN,  &
                                   NRCHOPIN, NPCHOPIN, NZCHOPIN, VISCPOW, VISCPIN)
! ==========================================================================================================
! [USAGE]:
! INITIALIZE THE GLOBAL FIELD (DOMAIN) SETUP, PROXYING FOU_INIT (IN FOURIER.MD) & LEG_INIT (IN LEGENDRE.MD)
! [PARAMETERS]:
! NRIN >> INTEGER VALUE USED FOR NR (# OF COLLOCATION PTS IN THE X/R DIRECTION)
! NPIN >> INTEGER VALUE USED FOR NP (# OF COLLOCATION PTS IN THE PHI DIRECTION)
! NZIN >> INTEGER VALUE USED FOR NZ (# OF COLLOCATION PTS IN THE Z DIRECTION)
! NRCHOPIN >> INTEGER VALUE USED FOR NRCHOP (# OF SPECTRAL ELEMENTS IN THE X/R DIR. #CHOP. .LE. #COLLOC.)
! NPCHOPIN >> INTEGER VALUE USED FOR NPCHOP (# OF SPECTRAL ELEMENTS IN THE PHI DIR. #CHOP. .LE. #COLLOC.)
! NZCHOPIN >> INTEGER VALUE USED FOR NZCHOP (# OF SPECTRAL ELEMENTS IN THE Z DIR. #CHOP. .LE. #COLLOC.)
! ELLIN >> REAL VALUE USED FOR ELL (LEGENDRE FUNCTIONS` MAP PARAMETER)
! ZLENIN >> REAL VALUE USED FOR ZLEN (DOMAIN LENGTH IN Z)
! VISCPOW >> THE HIGHEST ORDER OF (HYPER)VISCOUS TERMS TO BE CONSIDERED. DEAFULT IS 2.
! VISCIN >> GENERAL VISCOSITY NU.
! VISCPIN >> HYPERVISCOSITY (IF VISCPOW .LT. 2) NUP.
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: NRIN, NPIN, NZIN
      REAL(P8), INTENT(IN) :: ELLIN, ZLENIN
      INTEGER, OPTIONAL    :: NRCHOPIN, NPCHOPIN, NZCHOPIN
      INTEGER, OPTIONAL    :: VISCPOW
      REAL(P8), OPTIONAL   :: VISCIN, VISCPIN
      TYPE(FOU_TRANSFORM)  :: FTF

      ! INTEGER              :: I, J
      INTEGER              :: NRCHOPIN_, NPCHOPIN_, NZCHOPIN_, VISCPOW_
      REAL(P8)             :: VISCIN_, VISCPIN_

      NRCHOPIN_ = NRIN   ! DEFAULT: CHOP. INDEX IN R = NR
      NPCHOPIN_ = NPIN/2 ! DEFAULT: CHOP. INDEX IN PHI = NPH
      NZCHOPIN_ = NZIN/2 ! DEFAULT: CHOP. INDEX IN Z = NZ / 2
      VISCPOW_ = 2
      VISCIN_  = 0.D0    ! DEFAULT: INVISCID
      VISCPIN_ = 0.D0    ! DEFAULT: NO HYPERVISCOSITY

      IF (PRESENT(VISCPOW)) VISCPOW_ = VISCPOW
      IF (PRESENT(NRCHOPIN)) NRCHOPIN_ = NRCHOPIN
      IF (PRESENT(NPCHOPIN)) NPCHOPIN_ = NPCHOPIN
      IF (PRESENT(NZCHOPIN)) NZCHOPIN_ = NZCHOPIN
      IF (PRESENT(VISCIN)) VISCIN_ = VISCIN
      IF (PRESENT(VISCPIN)) VISCPIN_ = VISCPIN

      IF (VISCIN_ .LT. 0.D0) STOP 'FIELD_INIT: VISCOSITY CANNOT BE NEGATIVE.'
      IF (VISCPIN_ .LT. 0.D0) STOP 'FIELD_INIT: HYPERVISCOSITY CANNOT BE NEGATIVE.'
      IF ((VISCPOW_ .GT. 2) .AND. (VISCPIN_ .LT. 1.D-15)) THEN
        STOP 'FIELD_INIT: HYPERVISCOSITY CAN`T BE ZERO WHEN HYPERVISOCISTY ORDER >=4 IS SET.'
      ENDIF

      FTF = FOU_INIT(NZIN, NZCHOPIN_, ZLENIN)
      TFM = LEG_INIT(NRIN, NPIN, NRCHOPIN_, NPCHOPIN_, ELLIN, VISCPOW_)

      FINFO%R = TFM%R   ! COPYING INFO FROM TFM TO FINFO FOR CONVENIENCE. CAN REFER TO EITHER OF TWO.
      FINFO%P = TFM%P
      FINFO%PR = TFM%PR
      FINFO%PI = TFM%PI
      FINFO%ELL = TFM%L
      FINFO%ELL2 = FINFO%ELL**2.D0

      FINFO%Z = FTF%Z
      FINFO%ZLEN = FTF%ZLEN

      FINFO%NR = NX
      FINFO%NRH = NXH
      FINFO%NP = NP
      FINFO%NPH = NPH
      FINFO%NZ = NZ

      FINFO%NRDIM = NXDIM
      FINFO%NPDIM = NPDIM
      FINFO%NZDIM = NZDIM

      FINFO%NRCHOP = NXCHOP
      FINFO%NPCHOP = NPCHOP
      FINFO%NZCHOP = NZCHOP
      FINFO%NZCHOPU = NZCHOPU

      FINFO%NRCHOPS = NXCHOPS
      FINFO%NPCHOPS = NPCHOPS

      FINFO%M = M
      FINFO%AK = AK

      FINFO%VISCPOW = VISCPOW_
      FINFO%VISC = VISCIN_
      FINFO%VISCP = VISCPIN_

      RETURN
      END SUBROUTINE
! =========================================================================================================
      SUBROUTINE FIELD_SET_VECTOR(VFLD, Q, XO, YO)
! =========================================================================================================
! [USAGE]:
! GENERATE THE Q-VORTEX VELOCITY PROFILE (IN PPP SPACE)
! [VARAIABLES]:
! VFLD >> VECTOR FIELD WHERE (U_R, U_P, U_Z) = (0, (1-EXP(-R^2))/R, EXP(-R^2)/Q)
! Q >> Q VORTEX PARAMETER
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! =========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(INOUT) :: VFLD
      REAL(P8), INTENT(IN)              :: Q, XO, YO

      INTEGER                           :: I, J
      REAL(P8)                          :: RR, RI

      CALL ALLOC(VFLD, 'PPP')

      DO I = 1, FINFO%NR
        DO J = 1, FINFO%NPH
          RR = SQRT((FINFO%R(I)*COS(FINFO%PR(J)-XO))**2.D0+(FINFO%R(I)*SIN(FINFO%PR(J)-YO))**2.D0)
          RI = SQRT((FINFO%R(I)*COS(FINFO%PI(J)-XO))**2.D0+(FINFO%R(I)*SIN(FINFO%PI(J)-YO))**2.D0)
          VFLD%EP(I,J,:) = CMPLX((1.D0-EXP(-RR**2.D0))/RR, (1.D0-EXP(-RI**2.D0))/RI, P8)
          VFLD%EZ(I,J,:) = CMPLX(EXP(-RR**2.D0)/Q, EXP(-RI**2.D0)/Q, P8)
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE ALLOC_SCALAR(A, SP)
! ==========================================================================================================
! [USAGE]:
! ALLOCATE A SCALAR FIELD IN A SPECIFIED SPACE TAGGED BY 'PPP', 'PFP', 'FFP' OR 'FFF'.
! [VARIABLES]:
! A >> SCALAR FIELD DEFINED IN (R, P, Z)
! SP >> (OPTIONAL) TAG INDICATING WHICH SPACE A IS LOCATED. IF NOT PRESENT, 'FFF' IS DEFAULT.
! [NOTES]:
! SANGJOON LEE @ JUNE 2023 
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: A
      CHARACTER(LEN=3), INTENT(IN), OPTIONAL :: SP

      CHARACTER(LEN=3) :: SP_ = 'FFF'

      IF (PRESENT(SP)) SP_ = SP

      IF (SP_ .EQ. 'FFF') THEN
        A%SPACE = 'FFF'
      ELSEIF (SP_ .EQ. 'FFP') THEN
        A%SPACE = 'FFP'
      ELSEIF (SP_ .EQ. 'PFP') THEN
        A%SPACE = 'PFP'
      ELSEIF (SP_ .EQ. 'PPP') THEN
        A%SPACE = 'PPP'
      ELSE
        STOP 'ALLOC_SCALAR: INVALID SPACE TAG. ONLY ALLOWS FFF, FFP, PFP AND PPP.'
      ENDIF

      IF (ALLOCATED(A%E)) THEN
        A%E = 0.D0
      ELSE
        ALLOCATE( A%E(FINFO%NRDIM, FINFO%NPDIM, FINFO%NZDIM) )
        A%E = 0.D0
      ENDIF

      A%NRC = FINFO%NRCHOP
      A%NPC = FINFO%NPCHOP
      A%NZC = FINFO%NZCHOP
      A%NZCU = FINFO%NZCHOPU

      A%NRCS = FINFO%NRCHOPS
      A%NPCS = FINFO%NPCHOPS

      A%LN = 0.D0

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE ALLOC_VECTOR(A, SP)
! ==========================================================================================================
! [USAGE]:
! ALLOCATE A VECTOR FIELD IN A SPECIFIED SPACE TAGGED BY 'PPP', 'PFP', 'FFP' OR 'FFF'.
! [VARIABLES]:
! A >> VECTOR FIELD DEFINED IN (R, P, Z)
! SP >> (OPTIONAL) TAG INDICATING WHICH SPACE A IS LOCATED. IF NOT PRESENT, 'FFF' IS DEFAULT.
! [NOTES]:
! SANGJOON LEE @ JUNE 2023 
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(INOUT) :: A
      CHARACTER(LEN=3), INTENT(IN), OPTIONAL :: SP

      CHARACTER(LEN=3) :: SP_ = 'FFF'

      IF (PRESENT(SP)) SP_ = SP

      IF (SP_ .EQ. 'FFF') THEN
        A%SPACE = 'FFF'
      ELSEIF (SP_ .EQ. 'FFP') THEN
        A%SPACE = 'FFP'
      ELSEIF (SP_ .EQ. 'PFP') THEN
        A%SPACE = 'PFP'
      ELSEIF (SP_ .EQ. 'PPP') THEN
        A%SPACE = 'PPP'
      ELSE
        STOP 'ALLOC_SCALAR: INVALID SPACE TAG. ONLY ALLOWS FFF, FFP, PFP AND PPP.'
      ENDIF

      IF (ALLOCATED(A%ER)) THEN
        A%ER = 0.D0
      ELSE
        ALLOCATE( A%ER(FINFO%NRDIM, FINFO%NPDIM, FINFO%NZDIM) )
        A%ER = 0.D0
      ENDIF

      IF (ALLOCATED(A%EP)) THEN
        A%EP = 0.D0
      ELSE
        ALLOCATE( A%EP(FINFO%NRDIM, FINFO%NPDIM, FINFO%NZDIM) )
        A%EP = 0.D0
      ENDIF

      IF (ALLOCATED(A%EZ)) THEN
        A%EZ = 0.D0
      ELSE
        ALLOCATE( A%EZ(FINFO%NRDIM, FINFO%NPDIM, FINFO%NZDIM) )
        A%EZ = 0.D0
      ENDIF
      
      A%NRC = FINFO%NRCHOP
      A%NPC = FINFO%NPCHOP
      A%NZC = FINFO%NZCHOP
      A%NZCU = FINFO%NZCHOPU

      A%NRCS = FINFO%NRCHOPS
      A%NPCS = FINFO%NPCHOPS

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE DEALLOC_SCALAR(A)
! ==========================================================================================================
! [USAGE]:
! DEALLOCATE A SCALAR FIELD
! [VARIABLES]:
! A >> SCALAR FIELD TO BE DEALLOCATED
! [NOTES]:
! SANGJOON LEE @ JUNE 2023 
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: A

      IF (ALLOCATED(A%E)) DEALLOCATE(A%E)

      A%NRC = 0
      A%NPC = 0
      A%NZC = 0
      A%NZCU = 0

      IF (ALLOCATED(A%NRCS)) DEALLOCATE(A%NRCS)
      IF (ALLOCATED(A%NPCS)) DEALLOCATE(A%NPCS)

      A%LN = 0.D0

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE DEALLOC_VECTOR(A)
! ==========================================================================================================
! [USAGE]:
! DEALLOCATE A VECTOR FIELD
! [VARIABLES]:
! A >> VECTOR FIELD TO BE DEALLOCATED
! [NOTES]:
! SANGJOON LEE @ JUNE 2023 
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(INOUT) :: A

      IF (ALLOCATED(A%ER)) DEALLOCATE(A%ER)
      IF (ALLOCATED(A%EP)) DEALLOCATE(A%EP)
      IF (ALLOCATED(A%EZ)) DEALLOCATE(A%EZ)

      A%NRC = 0
      A%NPC = 0
      A%NZC = 0
      A%NZCU = 0

      IF (ALLOCATED(A%NRCS)) DEALLOCATE(A%NRCS)
      IF (ALLOCATED(A%NPCS)) DEALLOCATE(A%NPCS)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE CHOPSET_SCALAR(A, IOF1, IOF2, IOF3)
! ==========================================================================================================
! [USAGE]:
! SET UP NEW CHOPPING LOCATIONS FOR A SCALAR FIELD A.
! [VARIABLES]:
! A >> SCALAR FIELD
! IOF1 >> OFFSET OF NRCHOP
! IOF2 >> (OPTIONAL) OFFSET OF NPCHOP. DEFAULT IS 0.
! IOF3 >> (OPTIONAL) OFFSET OF NZCHOP. DEFAULT IS 0.
! [NOTES]:
! SANGJOON LEE @ JUNE 2023 
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: A
      INTEGER, INTENT(IN)               :: IOF1
      INTEGER, OPTIONAL                 :: IOF2, IOF3

      INTEGER                           :: MM, NN
      INTEGER                           :: IOF2_, IOF3_

      IOF2_ = 0
      IOF3_ = 0

      IF (PRESENT(IOF2)) IOF2_ = IOF2
      IF (PRESENT(IOF3)) IOF3_ = IOF3

      A%NRC = A%NRC + IOF1
      A%NPC = A%NPC + IOF2_
      A%NZC = A%NZC + IOF3_

      A%NZCU = A%NZCU - IOF3_

      IF (A%NRC .GT. FINFO%NRDIM) STOP 'CHOPSET: NRCHOP TOO LARGE.'
      IF (A%NPC .GT. FINFO%NPDIM) STOP 'CHOPSET: NPCHOP TOO LARGE.'
      IF (2*IOF3_ .GT. A%NZCU - A%NZC) STOP 'CHOPSET: NZCHOP TOO LARGE.'

      DO MM = 1, A%NPC
        A%NRCS(MM) = MAX(MIN(A%NRC, A%NRC-FINFO%M(MM)), 0)
      ENDDO
      A%NRCS(A%NPC+1:) = 0

      DO NN = 1, FINFO%NR
        A%NPCS(NN) = MAX(MIN(A%NPC, A%NRC-NN), 0)
      ENDDO

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE CHOPSET_VECTOR(A, IOF1, IOF2, IOF3)
! ==========================================================================================================
! [USAGE]:
! SET UP NEW CHOPPING LOCATIONS FOR A VECTOR FIELD A.
! [VARIABLES]:
! A >> VECTOR FIELD
! IOF1 >> OFFSET OF NRCHOP
! IOF2 >> (OPTIONAL) OFFSET OF NPCHOP. DEFAULT IS 0.
! IOF3 >> (OPTIONAL) OFFSET OF NZCHOP. DEFAULT IS 0.
! [NOTES]:
! SANGJOON LEE @ JUNE 2023 
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(INOUT) :: A
      INTEGER, INTENT(IN)               :: IOF1
      INTEGER, OPTIONAL                 :: IOF2, IOF3

      INTEGER                           :: MM, NN
      INTEGER                           :: IOF2_, IOF3_

      IOF2_ = 0
      IOF3_ = 0

      IF (PRESENT(IOF2)) IOF2_ = IOF2
      IF (PRESENT(IOF3)) IOF3_ = IOF3

      A%NRC = A%NRC + IOF1
      A%NPC = A%NPC + IOF2_
      A%NZC = A%NZC + IOF3_

      A%NZCU = A%NZCU - IOF3_

      IF (A%NRC .GT. FINFO%NRDIM) STOP 'CHOPSET: NRCHOP TOO LARGE.'
      IF (A%NPC .GT. FINFO%NPDIM) STOP 'CHOPSET: NPCHOP TOO LARGE.'
      IF (2*IOF3_ .GT. A%NZCU - A%NZC) STOP 'CHOPSET: NZCHOP TOO LARGE.'

      DO MM = 1, A%NPC
        A%NRCS(MM) = MAX(MIN(A%NRC, A%NRC-FINFO%M(MM)), 0)
      ENDDO
      A%NRCS(A%NPC+1:) = 0

      DO NN = 1, FINFO%NR
        A%NPCS(NN) = MAX(MIN(A%NPC, A%NRC-NN), 0)
      ENDDO

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE CHOPRST_SCALAR(A)
! ==========================================================================================================
! [USAGE]:
! RESET THE CHOPPING LOCATION OF A GIVEN SCALAR FIELD A BASED ON THE GLOBAL FIELD INFORMATION IN FINFO
! [VARIABLES]:
! A >> SCALAR FIELD TO RESET THE CHOPPING LOCATIONS
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: A

      INTEGER                           :: I1, I2, I3

      I1 = FINFO%NRCHOP - A%NRC
      I2 = FINFO%NPCHOP - A%NPC
      I3 = FINFO%NZCHOP - A%NZC

      CALL CHOPSET_SCALAR(A, I1, I2, I3)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE CHOPRST_VECTOR(A)
! ==========================================================================================================
! [USAGE]:
! RESET THE CHOPPING LOCATION OF A GIVEN VECTOR FIELD A BASED ON THE GLOBAL FIELD INFORMATION IN FINFO
! [VARIABLES]:
! A >> VECTOR FIELD TO RESET THE CHOPPING LOCATIONS
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(INOUT) :: A

      INTEGER                           :: I1, I2, I3

      I1 = FINFO%NRCHOP - A%NRC
      I2 = FINFO%NPCHOP - A%NPC
      I3 = FINFO%NZCHOP - A%NZC

      CALL CHOPSET_VECTOR(A, I1, I2, I3)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE CHOPDO_SCALAR(A)
! ==========================================================================================================
! [USAGE]:
! PERFORM THE CHOPPING FOLLOWING THE CHOPPING SETUP BY CHOPSET (OR BY DEFAULT)
! [VARIABLES]:
! A >> SCALAR FIELD TO BE CHOPPED
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: A

      INTEGER                           :: MM

      IF (A%SPACE(1:1) .EQ. 'F') THEN ! RADIAL (R)
        DO MM = 1, A%NPC
          A%E(A%NRCS(MM)+1:, MM, :) = 0.D0
        ENDDO
      ENDIF 

      IF (A%SPACE(2:2) .EQ. 'F') THEN ! AZIMUTHAL (P)
        A%E(:, A%NPC+1:, :) = 0.D0
      ENDIF

      IF (A%SPACE(3:3) .EQ. 'F') THEN ! AXIAL (Z)
        IF (A%NZC .LT. A%NZCU) A%E(:, :, A%NZC+1:A%NZCU-1) = 0.D0
      ENDIF

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE CHOPDO_VECTOR(A)
! ==========================================================================================================
! [USAGE]:
! PERFORM THE CHOPPING FOLLOWING THE CHOPPING SETUP BY CHOPSET (OR BY DEFAULT)
! [VARIABLES]:
! A >> VECTOR FIELD TO BE CHOPPED
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(INOUT) :: A

      INTEGER                           :: MM

      IF (A%SPACE(1:1) .EQ. 'F') THEN ! RADIAL (R)
        DO MM = 1, A%NPC
          A%ER(A%NRCS(MM)+1:, MM, :) = 0.D0
          A%EP(A%NRCS(MM)+1:, MM, :) = 0.D0
          A%EZ(A%NRCS(MM)+1:, MM, :) = 0.D0
        ENDDO
      ENDIF 

      IF (A%SPACE(2:2) .EQ. 'F') THEN ! AZIMUTHAL (P)
        A%ER(:, A%NPC+1:, :) = 0.D0
        A%EP(:, A%NPC+1:, :) = 0.D0
        A%EZ(:, A%NPC+1:, :) = 0.D0
      ENDIF

      IF (A%SPACE(3:3) .EQ. 'F') THEN ! AXIAL (Z)
        IF (A%NZC .LT. A%NZCU) A%ER(:, :, A%NZC+1:A%NZCU-1) = 0.D0
        IF (A%NZC .LT. A%NZCU) A%EP(:, :, A%NZC+1:A%NZCU-1) = 0.D0
        IF (A%NZC .LT. A%NZCU) A%EZ(:, :, A%NZC+1:A%NZCU-1) = 0.D0
      ENDIF

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE TRANS_SCALAR(A, SP)
! ==========================================================================================================
! [USAGE]:
! TRANSFORM THE SPACE WHERE A IS IN.
! [VARIABLES]:
! A >> SCALAR FIELD TO TRANSFORM
! SP >> SPACE TAG TO WHICH WE WANT TO TRANSFORM THE SPACE OF A
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: A
      CHARACTER(LEN=3), INTENT(IN)      :: SP

      INTEGER                           :: KK

      IF (.NOT. ((SP .EQ. 'PPP') .OR. (SP .EQ. 'PFP') .OR. (SP .EQ. 'FFP') .OR. (SP .EQ. 'FFF'))) THEN
        STOP 'TRANS: ONLY PPP, PFP, FFP AND FFF SPACES ALLOWED FOR WHERE A FIELD IS.'
      ENDIF

      IF ((A%SPACE(2:2) .EQ. 'P') .AND. (SP(2:2) .EQ. 'F')) THEN
        CALL HORFFT(A%E, A%NPC, 0)
        A%SPACE(2:2) = SP(2:2)
        CALL CHOPDO_SCALAR(A)
      ENDIF
      IF ((A%SPACE(1:1) .EQ. 'P') .AND. (SP(1:1) .EQ. 'F')) THEN
        IF (ABS(A%LN) .GT. 1.D-14) THEN
          DO KK = 1, FINFO%NZ
            A%E(:FINFO%NR, 1, KK) = A%E(:FINFO%NR, 1, KK) - A%LN*TFM%LN(:FINFO%NR)
          ENDDO
        ENDIF
        CALL RTRANS(A%E, A%NRCS, A%NPC, 0)
        A%SPACE(1:1) = SP(1:1)
        CALL CHOPDO_SCALAR(A)
      ENDIF
      IF ((A%SPACE(3:3) .EQ. 'P') .AND. (SP(3:3) .EQ. 'F')) THEN
        CALL VERFFT(A%E, A%NZC, A%NRCS, A%NPC, 0)
        A%SPACE(3:3) = SP(3:3)
        CALL CHOPDO_SCALAR(A)
      ENDIF

      IF ((A%SPACE(3:3) .EQ. 'F') .AND. (SP(3:3) .EQ. 'P')) THEN
        CALL VERFFT(A%E, A%NZC, A%NRCS, A%NPC, 1)
        A%SPACE(3:3) = SP(3:3)
        CALL CHOPDO_SCALAR(A)
      ENDIF
      IF ((A%SPACE(1:1) .EQ. 'F') .AND. (SP(1:1) .EQ. 'P')) THEN
        CALL RTRANS(A%E, A%NRCS, A%NPC, 1)
        IF (ABS(A%LN) .GT. 1.D-14) THEN
            DO KK = 1, FINFO%NZ
              A%E(:FINFO%NR, 1, KK) = A%E(:FINFO%NR, 1, KK) + A%LN*TFM%LN(:FINFO%NR)
            ENDDO
        ENDIF
        A%SPACE(1:1) = SP(1:1)
        CALL CHOPDO_SCALAR(A)
      ENDIF
      IF ((A%SPACE(2:2) .EQ. 'F') .AND. (SP(2:2) .EQ. 'P')) THEN
        CALL HORFFT(A%E, A%NPC, 1)
        A%SPACE(2:2) = SP(2:2)
        CALL CHOPDO_SCALAR(A)
      ENDIF

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE TRANS_VECTOR(A, SP)
! ==========================================================================================================
! [USAGE]:
! TRANSFORM THE SPACE WHERE A IS IN.
! [VARIABLES]:
! A >> SCALAR FIELD TO TRANSFORM
! SP >> SPACE TAG TO WHICH WE WANT TO TRANSFORM THE SPACE OF A
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(INOUT) :: A
      CHARACTER(LEN=3), INTENT(IN)      :: SP

      IF (.NOT. ((SP .EQ. 'PPP') .OR. (SP .EQ. 'PFP') .OR. (SP .EQ. 'FFP') .OR. (SP .EQ. 'FFF'))) THEN
        STOP 'TRANS: ONLY PPP, PFP, FFP AND FFF SPACES ALLOWED FOR WHERE A FIELD IS.'
      ENDIF

      IF ((A%SPACE(2:2) .EQ. 'P') .AND. (SP(2:2) .EQ. 'F')) THEN
        CALL HORFFT(A%ER, A%NPC, 0)
        CALL HORFFT(A%EP, A%NPC, 0)
        CALL HORFFT(A%EZ, A%NPC, 0)
        A%SPACE(2:2) = SP(2:2)
        CALL CHOPDO_VECTOR(A)
      ENDIF
      IF ((A%SPACE(1:1) .EQ. 'P') .AND. (SP(1:1) .EQ. 'F')) THEN
        CALL RTRANS(A%ER, A%NRCS, A%NPC, 0)
        CALL RTRANS(A%EP, A%NRCS, A%NPC, 0)
        CALL RTRANS(A%EZ, A%NRCS, A%NPC, 0)
        A%SPACE(1:1) = SP(1:1)
        CALL CHOPDO_VECTOR(A)
      ENDIF
      IF ((A%SPACE(3:3) .EQ. 'P') .AND. (SP(3:3) .EQ. 'F')) THEN
        CALL VERFFT(A%ER, A%NZC, A%NRCS, A%NPC, 0)
        CALL VERFFT(A%EP, A%NZC, A%NRCS, A%NPC, 0)
        CALL VERFFT(A%EZ, A%NZC, A%NRCS, A%NPC, 0)
        A%SPACE(3:3) = SP(3:3)
        CALL CHOPDO_VECTOR(A)
      ENDIF

      IF ((A%SPACE(3:3) .EQ. 'F') .AND. (SP(3:3) .EQ. 'P')) THEN
        CALL VERFFT(A%ER, A%NZC, A%NRCS, A%NPC, 1)
        CALL VERFFT(A%EP, A%NZC, A%NRCS, A%NPC, 1)
        CALL VERFFT(A%EZ, A%NZC, A%NRCS, A%NPC, 1)
        A%SPACE(3:3) = SP(3:3)
        CALL CHOPDO_VECTOR(A)
      ENDIF
      IF ((A%SPACE(1:1) .EQ. 'F') .AND. (SP(1:1) .EQ. 'P')) THEN
        CALL RTRANS(A%ER, A%NRCS, A%NPC, 1)
        CALL RTRANS(A%EP, A%NRCS, A%NPC, 1)
        CALL RTRANS(A%EZ, A%NRCS, A%NPC, 1)
        A%SPACE(1:1) = SP(1:1)
        CALL CHOPDO_VECTOR(A)
      ENDIF
      IF ((A%SPACE(2:2) .EQ. 'F') .AND. (SP(2:2) .EQ. 'P')) THEN
        CALL HORFFT(A%ER, A%NPC, 1)
        CALL HORFFT(A%EP, A%NPC, 1)
        CALL HORFFT(A%EZ, A%NPC, 1)
        A%SPACE(2:2) = SP(2:2)
        CALL CHOPDO_VECTOR(A)
      ENDIF

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      FUNCTION CALCAT_EXACT_SCALAR(A, R, P, Z) RESULT(CALC)
! ==========================================================================================================
! [USAGE]:
! CALCULATES THE FUNCTION VALUES AT AN ARBITRARY LOCATION (R, P, Z)
! [VARIABLES]:
! A >> SCALAR FIELD
! R >> RADIAL COORDINATE R
! P >> AZIMUTHAL COORDINIATE PHI
! Z >> AXIAL COORDINATE Z
! CALC >> FUNCTION VALUE OF A(R, P, Z)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(IN)           :: A
      REAL(P8), INTENT(IN)                     :: R, P, Z
      REAL(P8)                                 :: CALC

      INTEGER                                  :: NN, MM, KK, MV
      REAL(P8)                                 :: KV, X
      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: AZ
      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: FT
      TYPE(SCALAR_FIELD)                       :: AT

      AT = A
      IF (AT%SPACE .NE. 'FFF') THEN
        IF (WARN_ON) WRITE(*,*) 'CALCAT: [WARNING] THE FIELD TO BE EVALUATED IS NOT IN FFF SPACE.'
        CALL TRANS_SCALAR(AT, 'FFF')
      ENDIF

      ALLOCATE( AZ(AT%NRC,AT%NPC) )
      AZ = 0.D0
      DO KK = 1, FINFO%NZ
        IF ((KK .LE. AT%NZC) .OR. (KK .GE. AT%NZCU)) THEN
          KV = FINFO%AK(KK)
          DO MM = 1, AT%NPC
            DO NN = 1, AT%NRCS(MM)
              AZ(NN,MM) = AZ(NN,MM) + AT%E(NN,MM,KK)*EXP(IU*KV*Z)
            ENDDO
          ENDDO
        ENDIF
      ENDDO

      X = -(R**2.D0 - FINFO%ELL2)/(R**2.D0 + FINFO%ELL2)
      FT = RESHAPE(LEG_TBL((/X/), AT%NRC, FINFO%M(1:AT%NPC), TFM%LOGNORM), (/ AT%NRC, AT%NPC /))

      DO MM = 1, AT%NPC
        MV = FINFO%M(MM)
        FT(:,MM) = FT(:,MM) * EXP(IU*MV*P)
      ENDDO

      CALC = SUM(AZ*FT)*2.D0

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION CALCAT_EXACT_VECTOR(A, R, P, Z) RESULT(CALC)
! ==========================================================================================================
! [USAGE]:
! CALCULATES THE FUNCTION VALUES AT AN ARBITRARY LOCATION (R, P, Z)
! [VARIABLES]:
! A >> VECTOR FIELD
! R >> RADIAL COORDINATE R
! P >> AZIMUTHAL COORDINIATE PHI
! Z >> AXIAL COORDINATE Z
! CALC >> FUNCTION VALUE OF A(R, P, Z). 1ST: A_R, 2ND: A_P, 3RD: A_Z
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(IN)             :: A
      REAL(P8), INTENT(IN)                       :: R, P, Z
      REAL(P8), DIMENSION(3)                     :: CALC

      INTEGER                                    :: NN, MM, KK, MV
      REAL(P8)                                   :: KV, X
      COMPLEX(P8), DIMENSION(:,:,:), ALLOCATABLE :: AZ
      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE   :: FT
      TYPE(VECTOR_FIELD)                         :: AT

      AT = A
      IF (AT%SPACE .NE. 'FFF') THEN
        IF (WARN_ON) WRITE(*,*) 'CALCAT: [WARNING] THE FIELD TO BE EVALUATED IS NOT IN FFF SPACE.'
        CALL TRANS_VECTOR(AT, 'FFF')
      ENDIF

      ALLOCATE( AZ(AT%NRC,AT%NPC,3) )
      AZ = 0.D0
      DO KK = 1, FINFO%NZ
        IF ((KK .LE. AT%NZC) .OR. (KK .GE. AT%NZCU)) THEN
          KV = FINFO%AK(KK)
          DO MM = 1, AT%NPC
            DO NN = 1, AT%NRCS(MM)
              AZ(NN,MM,1) = AZ(NN,MM,1) + AT%ER(NN,MM,KK)*EXP(IU*KV*Z)
              AZ(NN,MM,2) = AZ(NN,MM,2) + AT%EP(NN,MM,KK)*EXP(IU*KV*Z)
              AZ(NN,MM,3) = AZ(NN,MM,3) + AT%EZ(NN,MM,KK)*EXP(IU*KV*Z)
            ENDDO
          ENDDO
        ENDIF
      ENDDO

      X = -(R**2.D0 - FINFO%ELL2)/(R**2.D0 + FINFO%ELL2)
      FT = RESHAPE(LEG_TBL((/X/), AT%NRC, FINFO%M(1:AT%NPC), TFM%LOGNORM), (/ AT%NRC, AT%NPC /))

      DO MM = 1, AT%NPC
        MV = FINFO%M(MM)
        FT(:,MM) = FT(:,MM) * EXP(IU*MV*P)
      ENDDO

      CALC(1) = SUM(AZ(:,:,1)*FT)*2.D0
      CALC(2) = SUM(AZ(:,:,2)*FT)*2.D0
      CALC(3) = SUM(AZ(:,:,3)*FT)*2.D0

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION CALCAT_LINTP_SCALAR_MULTI(A, RIN, PIN, ZIN) RESULT(CALC)
! ==========================================================================================================
! [USAGE]:
! CALCULATES AN APPROXIMATE VALUE OF THE FUNCTION VALUE AT AN ARBITRARY LOCATION (R, P, Z)
! [VARIABLES]:
! A >> SCALAR FIELD
! RIN >> RADIAL COORDINATE R
! PIN >> AZIMUTHAL COORDINIATE PHI
! ZIN >> AXIAL COORDINATE Z
! CALC >> APPROX. FUNCTION VALUE OF A(R, P, Z)
! [NOTES]:
! 1. THIS SUROUTINE UES TRILINEAR INTERPOLATION. IN A UNIT BOX [0,1] X [0,1] X [0,1], IF ALL 8 VALUES AT
! THE VERTICES ARE KNOWN AS V000, ... V111, THEN A FUNCTION VALUE AT (X,Y,Z) IS INTERPOLATIVELY KNOWN AS
! F(X,Y,Z) = V000 (1-X)(1-Y)(1-Z) + V100 X(1-Y)(1-Z) + V010 (1-X)Y(1-Z) + V001 (1-X)(1-Y)Z +
!            V111 XYZ             + V110 XY(1-Z)     + V101 X(1-Y)Z     + V011 (1-X)YZ 
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(IN)          :: A
      REAL(P8), DIMENSION(:), INTENT(IN)      :: RIN, PIN, ZIN
      REAL(P8), DIMENSION(:), ALLOCATABLE     :: CALC

      INTEGER                                 :: I, J, K, N
      REAL(P8)                                :: X, Y, Z
      REAL(P8), DIMENSION(0:1,0:1,0:1)        :: V
      REAL(P8), DIMENSION(:,:,:), ALLOCATABLE :: AF
      TYPE(SCALAR_FIELD)                      :: AT

      IF (.NOT. (SIZE(RIN) .EQ. SIZE(PIN))) STOP 'CALCAT: INPUT SIZE INCOMPATIBLE.'
      IF (.NOT. (SIZE(RIN) .EQ. SIZE(ZIN))) STOP 'CALCAT: INPUT SIZE INCOMPATIBLE.'

      ALLOCATE( CALC(SIZE(RIN)) )

      AT = A
      IF (AT%SPACE .NE. 'PPP') THEN
        IF (WARN_ON) WRITE(*,*) 'CALCAT: [WARNING] THE FIELD TO BE EVALUATED IS NOT IN PPP SPACE.'
        CALL TRANS_SCALAR(AT, 'PPP')
      ENDIF

      ALLOCATE( AF(0:FINFO%NR+1, 1:FINFO%NP+1, 1:FINFO%NZ+1) )

      AF(1:FINFO%NR,1:FINFO%NP:2,1:FINFO%NZ) = REAL(AT%E)
      AF(1:FINFO%NR,2:FINFO%NP:2,1:FINFO%NZ) = AIMAG(AT%E)
      AF(0,FINFO%NP+1,:FINFO%NZ) = CALCAT0_SCALAR(AT)
      AF(FINFO%NR+1,FINFO%NP+1,:FINFO%NZ) = CALCAT1_SCALAR(AT)
      DO J = 1, FINFO%NP
        AF(0,J,:FINFO%NZ) = AF(0,FINFO%NP+1,:FINFO%NZ)
        AF(FINFO%NR+1,J,:FINFO%NZ) = AF(FINFO%NR+1,FINFO%NP+1,:FINFO%NZ)
      ENDDO
      AF(:,FINFO%NP+1,:) = AF(:,1,:)
      AF(:,:,FINFO%NZ+1) = AF(:,:,1)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(N,I,J,K,X,Y,Z,V)
      DO N = 1, SIZE(RIN)
        IF (RIN(N) .LT. FINFO%R(1)) THEN
          I = 0
        ELSE
          DO I = 1, FINFO%NR
            IF ((RIN(N) - FINFO%R(I))*(RIN(N) - FINFO%R(I+1)) .LE. 0.D0) EXIT
          ENDDO
        ENDIF
        DO J = 1, FINFO%NP
          IF ((PIN(N) - FINFO%P(J))*(PIN(N) - FINFO%P(J+1)) .LE. 0.D0) EXIT
        ENDDO
        DO K = 1, FINFO%NZ
          IF ((ZIN(N) - FINFO%Z(K))*(ZIN(N) - FINFO%Z(K+1)) .LE. 0.D0) EXIT
        ENDDO

        IF (I .EQ. 0) THEN
          X = RIN(N) / FINFO%R(I+1)
        ELSEIF (I .EQ. FINFO%NR) THEN
          X = MAX((RIN(N) - FINFO%R(I))/ FINFO%R(I), 1.D0)
        ELSE
          X = (RIN(N) - FINFO%R(I))/(FINFO%R(I+1) - FINFO%R(I))
        ENDIF
        Y = (PIN(N) - FINFO%P(J))/(FINFO%P(J+1) - FINFO%P(J))
        Z = (ZIN(N) - FINFO%Z(K))/(FINFO%Z(K+1) - FINFO%Z(K))

        V(0,0,0) = AF(I  , J  , K  )
        V(0,1,0) = AF(I  , J+1, K  )
        V(1,0,0) = AF(I+1, J  , K  )
        V(1,1,0) = AF(I+1, J+1, K  )
        V(0,0,1) = AF(I  , J  , K+1)
        V(0,1,1) = AF(I  , J+1, K+1)
        V(1,0,1) = AF(I+1, J  , K+1)
        V(1,1,1) = AF(I+1, J+1, K+1)

        CALC(N) = V(0,0,0)*(1.D0-X)*(1.D0-Y)*(1.D0-Z) + V(1,0,0)*      X *(1.D0-Y)*(1.D0-Z) + &
                  V(0,1,0)*(1.D0-X)*      Y *(1.D0-Z) + V(0,0,1)*(1.D0-X)*(1.D0-Y)*      Z  + &
                  V(1,1,1)*      X *      Y *      Z  + V(1,1,0)*      X *      Y *(1.D0-Z) + &
                  V(1,0,1)*      X *(1.D0-Y)*      Z  + V(0,1,1)*(1.D0-X)*      Y *      Z
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION CALCAT_LINTP_SCALAR_SINGLE(A, RIN, PIN, ZIN) RESULT(CALC)
! ==========================================================================================================
! [USAGE]:
! CALCULATES AN APPROXIMATE VALUE OF THE FUNCTION VALUE AT AN ARBITRARY LOCATION (R, P, Z)
! [VARIABLES]:
! A >> SCALAR FIELD
! RIN >> RADIAL COORDINATE R
! PIN >> AZIMUTHAL COORDINIATE PHI
! ZIN >> AXIAL COORDINATE Z
! CALC >> APPROX. FUNCTION VALUE OF A(R, P, Z)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(IN) :: A
      REAL(P8), INTENT(IN)           :: RIN, PIN, ZIN
      REAL(P8)                       :: CALC

      REAL(P8), DIMENSION(1)         :: T

      T = CALCAT_LINTP_SCALAR_MULTI(A, (/RIN/), (/PIN/), (/ZIN/))
      CALC = T(1)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION CALCAT_LINTP_VECTOR_MULTI(A, RIN, PIN, ZIN) RESULT(CALC)
! ==========================================================================================================
! [USAGE]:
! CALCULATES AN APPROXIMATE VALUE OF THE FUNCTION VALUE AT AN ARBITRARY LOCATION (R, P, Z)
! [VARIABLES]:
! A >> VECTOR FIELD
! RIN >> RADIAL COORDINATE R
! PIN >> AZIMUTHAL COORDINIATE PHI
! ZIN >> AXIAL COORDINATE Z
! CALC >> APPROX. FUNCTION VALUE OF A(R, P, Z). 1ST: A_R, 2ND: A_P, 3RD: A_Z
! [NOTES]:
! 1. THIS SUROUTINE UES TRILINEAR INTERPOLATION. IN A UNIT BOX [0,1] X [0,1] X [0,1], IF ALL 8 VALUES AT
! THE VERTICES ARE KNOWN AS V000, ... V111, THEN A FUNCTION VALUE AT (X,Y,Z) IS INTERPOLATIVELY KNOWN AS
! F(X,Y,Z) = V000 (1-X)(1-Y)(1-Z) + V100 X(1-Y)(1-Z) + V010 (1-X)Y(1-Z) + V001 (1-X)(1-Y)Z +
!            V111 XYZ             + V110 XY(1-Z)     + V101 X(1-Y)Z     + V011 (1-X)YZ 
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(IN)             :: A
      REAL(P8), DIMENSION(:), INTENT(IN)         :: RIN, PIN, ZIN
      REAL(P8), DIMENSION(:,:), ALLOCATABLE      :: CALC

      INTEGER                                    :: I, J, K, N
      REAL(P8)                                   :: X, Y, Z
      REAL(P8), DIMENSION(0:1,0:1,0:1,3)         :: V
      REAL(P8), DIMENSION(:,:,:,:), ALLOCATABLE  :: AF
      TYPE(VECTOR_FIELD)                         :: AT

      IF (.NOT. (SIZE(RIN) .EQ. SIZE(PIN))) STOP 'CALCAT: INPUT SIZE INCOMPATIBLE.'
      IF (.NOT. (SIZE(RIN) .EQ. SIZE(ZIN))) STOP 'CALCAT: INPUT SIZE INCOMPATIBLE.'

      ALLOCATE( CALC(3,SIZE(RIN)) )
      CALC = 0.D0

      AT = A
      IF (AT%SPACE .NE. 'PPP') THEN
        IF (WARN_ON) WRITE(*,*) 'CALCAT: [WARNING] THE FIELD TO BE EVALUATED IS NOT IN PPP SPACE.'
        CALL TRANS_VECTOR(AT, 'PPP')
      ENDIF

      ALLOCATE( AF(0:FINFO%NR+1, 1:FINFO%NP+1, 1:FINFO%NZ+1,3) )

      AF(1:FINFO%NR,1:FINFO%NP:2,1:FINFO%NZ,1) = REAL(AT%ER)
      AF(1:FINFO%NR,2:FINFO%NP:2,1:FINFO%NZ,1) = AIMAG(AT%ER)
      AF(1:FINFO%NR,1:FINFO%NP:2,1:FINFO%NZ,2) = REAL(AT%EP)
      AF(1:FINFO%NR,2:FINFO%NP:2,1:FINFO%NZ,2) = AIMAG(AT%EP)
      AF(1:FINFO%NR,1:FINFO%NP:2,1:FINFO%NZ,3) = REAL(AT%EZ)
      AF(1:FINFO%NR,2:FINFO%NP:2,1:FINFO%NZ,3) = AIMAG(AT%EZ)
      AF(0,FINFO%NP+1,:FINFO%NZ,:) = CALCAT0_VECTOR(AT)
      AF(FINFO%NR+1,FINFO%NP+1,:FINFO%NZ,:) = CALCAT1_VECTOR(AT)
      DO J = 1, FINFO%NP
        AF(0,J,:FINFO%NZ,:) = AF(0,FINFO%NP+1,:FINFO%NZ,:)
        AF(FINFO%NR+1,J,:FINFO%NZ,:) = AF(FINFO%NR+1,FINFO%NP+1,:FINFO%NZ,:)
      ENDDO
      AF(:,FINFO%NP+1,:,:) = AF(:,1,:,:)
      AF(:,:,FINFO%NZ+1,:) = AF(:,:,1,:)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(N,I,J,K,X,Y,Z,V)
      DO N = 1, SIZE(RIN)
        IF (RIN(N) .LT. FINFO%R(1)) THEN
          I = 0
        ELSE
          DO I = 1, FINFO%NR
            IF ((RIN(N) - FINFO%R(I))*(RIN(N) - FINFO%R(I+1)) .LE. 0.D0) EXIT
          ENDDO
        ENDIF
        DO J = 1, FINFO%NP
          IF ((PIN(N) - FINFO%P(J))*(PIN(N) - FINFO%P(J+1)) .LE. 0.D0) EXIT
        ENDDO
        DO K = 1, FINFO%NZ
          IF ((ZIN(N) - FINFO%Z(K))*(ZIN(N) - FINFO%Z(K+1)) .LE. 0.D0) EXIT
        ENDDO

        IF (I .EQ. 0) THEN
          X = RIN(N) / FINFO%R(I+1)
        ELSEIF (I .EQ. FINFO%NR) THEN
          X = MAX((RIN(N) - FINFO%R(I))/ FINFO%R(I), 1.D0)
        ELSE
          X = (RIN(N) - FINFO%R(I))/(FINFO%R(I+1) - FINFO%R(I))
        ENDIF
        Y = (PIN(N) - FINFO%P(J))/(FINFO%P(J+1) - FINFO%P(J))
        Z = (ZIN(N) - FINFO%Z(K))/(FINFO%Z(K+1) - FINFO%Z(K))

        V(0,0,0,:) = AF(I  , J  , K  , :)
        V(0,1,0,:) = AF(I  , J+1, K  , :)
        V(1,0,0,:) = AF(I+1, J  , K  , :)
        V(1,1,0,:) = AF(I+1, J+1, K  , :)
        V(0,0,1,:) = AF(I  , J  , K+1, :)
        V(0,1,1,:) = AF(I  , J+1, K+1, :)
        V(1,0,1,:) = AF(I+1, J  , K+1, :)
        V(1,1,1,:) = AF(I+1, J+1, K+1, :)

        CALC(:,N) = V(0,0,0,:)*(1.D0-X)*(1.D0-Y)*(1.D0-Z) + V(1,0,0,:)*      X *(1.D0-Y)*(1.D0-Z) + &
                    V(0,1,0,:)*(1.D0-X)*      Y *(1.D0-Z) + V(0,0,1,:)*(1.D0-X)*(1.D0-Y)*      Z  + &
                    V(1,1,1,:)*      X *      Y *      Z  + V(1,1,0,:)*      X *      Y *(1.D0-Z) + &
                    V(1,0,1,:)*      X *(1.D0-Y)*      Z  + V(0,1,1,:)*(1.D0-X)*      Y *      Z
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION CALCAT_LINTP_VECTOR_SINGLE(A, RIN, PIN, ZIN) RESULT(CALC)
! ==========================================================================================================
! [USAGE]:
! CALCULATES AN APPROXIMATE VALUE OF THE FUNCTION VALUE AT AN ARBITRARY LOCATION (R, P, Z)
! [VARIABLES]:
! A >> VECTOR FIELD
! RIN >> RADIAL COORDINATE R
! PIN >> AZIMUTHAL COORDINIATE PHI
! ZIN >> AXIAL COORDINATE Z
! CALC >> APPROX. FUNCTION VALUE OF A(R, P, Z). 1ST: A_R, 2ND: A_P, 3RD: A_Z
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(IN)             :: A
      REAL(P8), INTENT(IN)                       :: RIN, PIN, ZIN
      REAL(P8), DIMENSION(3)                     :: CALC

      CALC = RESHAPE(CALCAT_LINTP_VECTOR_MULTI(A, (/RIN/), (/PIN/), (/ZIN/)), (/ 3 /))

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION CALCAT0_SCALAR(A) RESULT(CALC)
! ==========================================================================================================
! [USAGE]:
! CALCULATES THE FUNCTION VALUES AT R = 0 (AND PHI = 0)
! [VARIABLES]:
! A >> SCALAR FIELD
! CALC >> FUNCTION VALUES OF A AT R = 0 (AND PHI = 0)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(IN)    :: A
      COMPLEX(P8), DIMENSION(FINFO%NZDIM) :: CALC

      TYPE(SCALAR_FIELD)                :: AT

      AT = A
      CALL TRANS_SCALAR(AT, 'FFF')

      CALC = TRANSPOSE(AT%E(:AT%NRC, 1, :)) .MUL. TFM%AT0(:AT%NRC)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION CALCAT0_VECTOR(A) RESULT(CALC)
! ==========================================================================================================
! [USAGE]:
! CALCULATES THE FUNCTION VALUES AT R = 0 (AND PHI = 0)
! [VARIABLES]:
! A >> VECTOR FIELD
! CALC >> FUNCTION VALUES OF A AT R = 0 (AND PHI = 0)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(IN)      :: A
      COMPLEX(P8), DIMENSION(FINFO%NZDIM,3) :: CALC

      TYPE(VECTOR_FIELD)                  :: AT

      AT = A
      CALL TRANS_VECTOR(AT, 'FFF')

      CALC(:,1) = TRANSPOSE(AT%ER(:AT%NRC, 1, :)) .MUL. TFM%AT0(:AT%NRC)
      CALC(:,2) = TRANSPOSE(AT%EP(:AT%NRC, 1, :)) .MUL. TFM%AT0(:AT%NRC)
      CALC(:,3) = TRANSPOSE(AT%EZ(:AT%NRC, 1, :)) .MUL. TFM%AT0(:AT%NRC)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION CALCAT1_SCALAR(A) RESULT(CALC)
! ==========================================================================================================
! [USAGE]:
! CALCULATES THE FUNCTION VALUES AT R = INF (AND PHI = 0)
! [VARIABLES]:
! A >> SCALAR FIELD
! CALC >> FUNCTION VALUES OF A AT R = INF (AND PHI = 0)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(IN)    :: A
      COMPLEX(P8), DIMENSION(FINFO%NZDIM) :: CALC

      TYPE(SCALAR_FIELD)                :: AT

      AT = A
      CALL TRANS_SCALAR(AT, 'FFF')

      CALC = TRANSPOSE(AT%E(:AT%NRC, 1, :)) .MUL. TFM%AT1(:AT%NRC)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION CALCAT1_VECTOR(A) RESULT(CALC)
! ==========================================================================================================
! [USAGE]:
! CALCULATES THE FUNCTION VALUES AT R = INF (AND PHI = 0)
! [VARIABLES]:
! A >> VECTOR FIELD
! CALC >> FUNCTION VALUES OF A AT R = INF (AND PHI = 0)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(IN)      :: A
      COMPLEX(P8), DIMENSION(FINFO%NZDIM,3) :: CALC

      TYPE(VECTOR_FIELD)                  :: AT

      AT = A
      CALL TRANS_VECTOR(AT, 'FFF')

      CALC(:,1) = TRANSPOSE(AT%ER(:AT%NRC, 1, :)) .MUL. TFM%AT1(:AT%NRC)
      CALC(:,2) = TRANSPOSE(AT%EP(:AT%NRC, 1, :)) .MUL. TFM%AT1(:AT%NRC)
      CALC(:,3) = TRANSPOSE(AT%EZ(:AT%NRC, 1, :)) .MUL. TFM%AT1(:AT%NRC)

      RETURN
      END FUNCTION
! ==========================================================================================================
      SUBROUTINE ZEROAT1_SCALAR(A)
! ==========================================================================================================
! [USAGE]:
! MAKE THE FIELD VALUE AT INFITITY ZERO BY ADJUSTING THE P_L_0^0(R)'S COEFFICIENTS AT EACH K.
! [VARIABLES]:
! A >> SCALAR FIELD
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT)        :: A

      COMPLEX(P8), DIMENSION(:), ALLOCATABLE   :: CALC
      CHARACTER(LEN=3)                         :: ORIG_SP

      ORIG_SP = A%SPACE
      CALL TRANS_SCALAR(A, 'FFF')

      CALC = CALCAT1_SCALAR(A)

      A%E(1,1,:) = A%E(1,1,:) - CALC/TFM%AT1(1)

      CALL TRANS_SCALAR(A, ORIG_SP)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE ZEROAT1_VECTOR(A)
! ==========================================================================================================
! [USAGE]:
! MAKE THE FIELD VALUE AT INFITITY ZERO BY ADJUSTING THE P_L_0^0(R)'S COEFFICIENTS AT EACH K.
! [VARIABLES]:
! A >> VECTOR FIELD
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(INOUT)           :: A

      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: CALC
      CHARACTER(LEN=3)                         :: ORIG_SP

      ORIG_SP = A%SPACE
      CALL TRANS_VECTOR(A, 'FFF')

      CALC = CALCAT1_VECTOR(A)

      A%ER(1,1,:) = A%ER(1,1,:) - CALC(:,1)/TFM%AT1(1)
      A%EP(1,1,:) = A%EP(1,1,:) - CALC(:,2)/TFM%AT1(1)
      A%EZ(1,1,:) = A%EZ(1,1,:) - CALC(:,3)/TFM%AT1(1)

      CALL TRANS_VECTOR(A, ORIG_SP)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      FUNCTION INTEG_SCALAR(A) RESULT(INTEG)
! ==========================================================================================================
! [USAGE]:
! INTEGRATION OF A SCALAR FIELD A OVER A DOMAIN 
!     0  <  R   < INFTY
!     0  <  PHI < 2*PI
!     0  <  Z   < ZLEN
! [VARIABLES]:
! A >> A SCALAR FIELD INPUT FOR INTEGRATION, WHICH MUST DECAY AS FAST AS O(1/R^4) AS R -> INF.
! INTEG >> INTEGRATE(A, RDRDPDZ)_R=0~INF, P=0~2PI, Z=0~ZLEN
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(IN) :: A
      REAL(P8)                       :: INTEG

      INTEGER                        :: MM, KK
      TYPE(SCALAR_FIELD)             :: AT

      AT = A
      CALL TRANS_SCALAR(AT, 'PPP')


!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, KK)
      DO MM = 1, FINFO%NPH
        DO KK = 1, FINFO%NZ
          AT%E(:FINFO%NR,MM,KK) = AT%E(:FINFO%NR,MM,KK)/(1.D0-TFM%X)**2.D0
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL TRANS_SCALAR(AT, 'FFF')

      INTEG = 4.D0*PI*FINFO%ZLEN*FINFO%ELL2*AT%E(1,1,1)*TFM%NORM(1,1)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION INTEG_VECTOR(A) RESULT(INTEG)
! ==========================================================================================================
! [USAGE]:
! INTEGRATION OF A SCALAR FIELD A OVER A DOMAIN 
!     0  <  R   < INFTY
!     0  <  PHI < 2*PI
!     0  <  Z   < ZLEN
! [VARIABLES]:
! A >> A VECTOR FIELD INPUT FOR INTEGRATION, WHOSE COMPONENTS MUST DECAY AS FAST AS O(1/R^4) AS R -> INF.
! INTEG >> INTEGRATE(A, RDRDPDZ)_R=0~INF, P=0~2PI, Z=0~ZLEN. 1ST: R, 2ND: P, 3RD: Z
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(IN) :: A
      REAL(P8), DIMENSION(3)         :: INTEG

      INTEGER                        :: MM, KK
      TYPE(VECTOR_FIELD)             :: AT

      AT = A
      CALL TRANS_VECTOR(AT, 'PPP')

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, KK)
      DO MM = 1, FINFO%NPH
        DO KK = 1, FINFO%NZ
          AT%ER(:FINFO%NR,MM,KK) = AT%ER(:FINFO%NR,MM,KK)/(1.D0-TFM%X)**2.D0
          AT%EP(:FINFO%NR,MM,KK) = AT%EP(:FINFO%NR,MM,KK)/(1.D0-TFM%X)**2.D0
          AT%EZ(:FINFO%NR,MM,KK) = AT%EZ(:FINFO%NR,MM,KK)/(1.D0-TFM%X)**2.D0
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL TRANS_VECTOR(AT, 'FFF')

      INTEG(1) = 4.D0*PI*FINFO%ZLEN*FINFO%ELL2*AT%ER(1,1,1)*TFM%NORM(1,1)
      INTEG(2) = 4.D0*PI*FINFO%ZLEN*FINFO%ELL2*AT%EP(1,1,1)*TFM%NORM(1,1)
      INTEG(3) = 4.D0*PI*FINFO%ZLEN*FINFO%ELL2*AT%EZ(1,1,1)*TFM%NORM(1,1)

      RETURN
      END FUNCTION
! ==========================================================================================================
      SUBROUTINE MSAVE_SCALAR(A, FN)
! ==========================================================================================================
! [USAGE]:
! SAVE THE SCALAR-TYPE VARIABLE A INTO FN
! [VARIABLES]:
! A >> SCALAR-TYPE VARIABLE TO SAVE
! FN >> FILENAME
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(IN) :: A
      CHARACTER(LEN=*)               :: FN

      INTEGER                        :: STATUS, I, J, K, FU

      OPEN(UNIT=FU, FILE=TRIM(ADJUSTL(FN)), STATUS='UNKNOWN', FORM='UNFORMATTED', IOSTAT=STATUS)
      IF (STATUS .NE.0) STOP 'MSAVE: FAILED TO CREATE/OPEN A FILE.'

      WRITE(FU) 'S'                             ! SAVING SCALAR ... 
      WRITE(FU) A%SPACE                            ! FIELD CURRENT SPACE
      WRITE(FU) FINFO%NRDIM, FINFO%NPDIM, FINFO%NZDIM ! FIELD STORAGE SIZE
      WRITE(FU) FINFO%NR, FINFO%NPH, FINFO%NZ         ! COLLOC PTS. SIZE
      WRITE(FU) FINFO%NRCHOP, FINFO%NPCHOP, FINFO%NZCHOP ! GLOBAL CHOPPING LOCATIONS
      WRITE(FU) A%NRC, A%NPC, A%NZC             ! SPECTRAL ELMTS. SIZE (CHOPPING LOCS.)
      WRITE(FU) FINFO%ZLEN, FINFO%ELL               ! Z-DOMAIN LEN. & MAP PARAMETER

      WRITE(FU) (((A%E(I,J,K), I=1, FINFO%NRDIM), J=1, FINFO%NPDIM), K=1, FINFO%NZDIM)
      WRITE(FU) A%LN

      CLOSE(FU)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE MSAVE_VECTOR(A, FN)
! ==========================================================================================================
! [USAGE]:
! SAVE THE SCALAR-TYPE VARIABLE A INTO FN
! [VARIABLES]:
! A >> SCALAR-TYPE VARIABLE TO SAVE
! FN >> FILENAME
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(IN) :: A
      CHARACTER(LEN=*)               :: FN

      INTEGER                        :: STATUS, FU
      INTEGER                        :: I, J, K

      OPEN(UNIT=FU, FILE=TRIM(ADJUSTL(FN)), STATUS='UNKNOWN', FORM='UNFORMATTED', IOSTAT=STATUS)
      IF (STATUS .NE.0) STOP 'MSAVE: FAILED TO CREATE/OPEN A FILE.'

      WRITE(FU) 'V'                                ! SAVING VECTOR ... 
      WRITE(FU) A%SPACE                            ! FIELD CURRENT SPACE
      WRITE(FU) FINFO%NRDIM, FINFO%NPDIM, FINFO%NZDIM    ! FIELD STORAGE SIZE
      WRITE(FU) FINFO%NR, FINFO%NPH, FINFO%NZ            ! COLLOC PTS. SIZE
      WRITE(FU) FINFO%NRCHOP, FINFO%NPCHOP, FINFO%NZCHOP ! GLOBAL CHOPPING LOCATIONS
      WRITE(FU) A%NRC, A%NPC, A%NZC                ! SPECTRAL ELMTS. SIZE (CHOPPING LOCS.)
      WRITE(FU) FINFO%ZLEN, FINFO%ELL                  ! Z-DOMAIN LEN. & MAP PARAMETER

      WRITE(FU) (((A%ER(I,J,K), I=1, FINFO%NRDIM), J=1, FINFO%NPDIM), K=1, FINFO%NZDIM)
      WRITE(FU) (((A%EP(I,J,K), I=1, FINFO%NRDIM), J=1, FINFO%NPDIM), K=1, FINFO%NZDIM)
      WRITE(FU) (((A%EZ(I,J,K), I=1, FINFO%NRDIM), J=1, FINFO%NPDIM), K=1, FINFO%NZDIM)

      CLOSE(FU)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE MLOAD_SCALAR(FN, A)
! ==========================================================================================================
! [USAGE]:
! LOAD THE SCALAR-TYPE VARIABLE A FROM FN
! [VARIABLES]:
! FN >> FILENAME WHERE THE INFORMATION IS STORED
! A >> SCALAR-TYPE VARIABLE TO BE LOADED
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      CHARACTER(LEN=*)                  :: FN
      TYPE(SCALAR_FIELD), INTENT(INOUT) :: A

      INTEGER                           :: STATUS, FU
      INTEGER                           :: I, J, K, NRDIM
      REAL(P8)                          :: Z, L
      CHARACTER(LEN=1)                  :: SORV
      CHARACTER(LEN=3)                  :: SP

      OPEN(UNIT=FU, FILE=TRIM(ADJUSTL(FN)), STATUS='OLD', FORM='UNFORMATTED', IOSTAT=STATUS)
      IF (STATUS .NE.0) STOP 'MLOAD: FAILED TO CREATE/OPEN A FILE.'

      READ(FU) SORV
      IF (SORV .NE. 'S') STOP 'MLOAD: SCALAR EXPECTED.'

      READ(FU) SP
      CALL ALLOC_SCALAR(A, SP)
      
      READ(FU) I, J, K
      IF (.NOT. ((I .EQ. FINFO%NRDIM) .AND. (J .EQ. FINFO%NPDIM) .AND. (K .EQ. FINFO%NZDIM))) THEN
        IF (WARN_ON) WRITE(*,*) 'MLOAD: [WARNING] FIELD SIZE INCOMPATIBLE.'
        IF (WARN_ON) WRITE(*,*) 'DIM1: FILE --> ', I, ', PROGRAM -->', FINFO%NRDIM
        IF (WARN_ON) WRITE(*,*) 'DIM2: FILE --> ', J, ', PROGRAM -->', FINFO%NPDIM
        IF (WARN_ON) WRITE(*,*) 'DIM3: FILE --> ', K, ', PROGRAM -->', FINFO%NZDIM
      ENDIF
      NRDIM = MIN(I, FINFO%NRDIM)
      IF (.NOT. (J .EQ. FINFO%NPDIM)) STOP
      IF (.NOT. (K .EQ. FINFO%NZDIM)) STOP

      READ(FU) I, J, K
      IF (.NOT. ((I .EQ. FINFO%NR) .AND. (J .EQ. FINFO%NPH) .AND. (K .EQ. FINFO%NZ))) THEN
        WRITE(*,*) 'MLOAD: COLLOCATION NOS. INCOMPATIBLE.'
        WRITE(*,*) 'R.COLLOC: FILE --> ', I, ', PROGRAM -->', FINFO%NR
        WRITE(*,*) 'P.COLLOC: FILE --> ', J, ', PROGRAM -->', FINFO%NPH
        WRITE(*,*) 'Z.COLLOC: FILE --> ', K, ', PROGRAM -->', FINFO%NZ
        STOP
      ENDIF
      
      READ(FU) I, J, K
      IF (.NOT. ((I .EQ. FINFO%NRCHOP) .AND. (J .EQ. FINFO%NPCHOP) .AND. (K .EQ. FINFO%NZCHOP))) THEN
        IF (WARN_ON) WRITE(*,*) 'MLOAD: [WARNING] CHOPPING LOCATION(S) INCONSISTENT.'
        IF (WARN_ON) WRITE(*,*) 'R.CHOP: FILE --> ', I, ', PROGRAM -->', FINFO%NRCHOP
        IF (WARN_ON) WRITE(*,*) 'P.CHOP: FILE --> ', J, ', PROGRAM -->', FINFO%NPCHOP
        IF (WARN_ON) WRITE(*,*) 'Z.CHOP: FILE --> ', K, ', PROGRAM -->', FINFO%NZCHOP
      ENDIF

      READ(FU) I, J, K
      A%NRC = I
      A%NPC = J
      A%NZC = K
      CALL CHOPSET_SCALAR(A, I-FINFO%NRCHOP, J-FINFO%NPCHOP, K-FINFO%NZCHOP)

      READ(FU) Z, L
      IF (.NOT. ((Z .EQ. FINFO%ZLEN) .AND. (L .EQ. FINFO%ELL))) THEN
        IF (WARN_ON) WRITE(*,*) 'MLOAD: [WARNING] ZLEN & ELL INCONSISTENT.'
        IF (WARN_ON) WRITE(*,*) 'ZLEN: FILE --> ', Z, ', PROGRAM -->', FINFO%ZLEN
        IF (WARN_ON) WRITE(*,*) 'ELL: FILE --> ', L, ', PROGRAM -->', FINFO%ELL
      ENDIF

      READ(FU) (((A%E(I,J,K), I=1, NRDIM), J=1, FINFO%NPDIM), K=1, FINFO%NZDIM)
      READ(FU) A%LN

      CLOSE(FU)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE MLOAD_VECTOR(FN, A)
! ==========================================================================================================
! [USAGE]:
! LOAD THE VECTOR-TYPE VARIABLE A FROM FN
! [VARIABLES]:
! FN >> FILENAME WHERE THE INFORMATION IS STORED
! A >> VECTOR-TYPE VARIABLE TO BE LOADED
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      CHARACTER(LEN=*)                  :: FN
      TYPE(VECTOR_FIELD), INTENT(INOUT) :: A

      INTEGER                           :: STATUS, FU
      INTEGER                           :: I, J, K, NRDIM
      REAL(P8)                          :: Z, L
      CHARACTER(LEN=1)                  :: SORV
      CHARACTER(LEN=3)                  :: SP

      OPEN(UNIT=FU, FILE=TRIM(ADJUSTL(FN)), STATUS='OLD', FORM='UNFORMATTED', IOSTAT=STATUS)
      IF (STATUS .NE.0) STOP 'MLOAD: FAILED TO CREATE/OPEN A FILE.'

      READ(FU) SORV
      IF (SORV .NE. 'V') STOP 'MLOAD: VECTOR EXPECTED.'

      READ(FU) SP
      CALL ALLOC_VECTOR(A, SP)
      
      READ(FU) I, J, K
      IF (.NOT. ((I .EQ. FINFO%NRDIM) .AND. (J .EQ. FINFO%NPDIM) .AND. (K .EQ. FINFO%NZDIM))) THEN
        IF (WARN_ON) WRITE(*,*) 'MLOAD: [WARNING] FIELD SIZE INCOMPATIBLE.'
        IF (WARN_ON) WRITE(*,*) 'DIM1: FILE --> ', I, ', PROGRAM -->', FINFO%NRDIM
        IF (WARN_ON) WRITE(*,*) 'DIM2: FILE --> ', J, ', PROGRAM -->', FINFO%NPDIM
        IF (WARN_ON) WRITE(*,*) 'DIM3: FILE --> ', K, ', PROGRAM -->', FINFO%NZDIM
      ENDIF
      NRDIM = MIN(I, FINFO%NRDIM)
      IF (.NOT. (J .EQ. FINFO%NPDIM)) STOP
      IF (.NOT. (K .EQ. FINFO%NZDIM)) STOP

      READ(FU) I, J, K
      IF (.NOT. ((I .EQ. FINFO%NR) .AND. (J .EQ. FINFO%NPH) .AND. (K .EQ. FINFO%NZ))) THEN
        WRITE(*,*) 'MLOAD: COLLOCATION NOS. INCOMPATIBLE.'
        WRITE(*,*) 'R.COLLOC: FILE --> ', I, ', PROGRAM -->', FINFO%NR
        WRITE(*,*) 'P.COLLOC: FILE --> ', J, ', PROGRAM -->', FINFO%NPH
        WRITE(*,*) 'Z.COLLOC: FILE --> ', K, ', PROGRAM -->', FINFO%NZ
        STOP
      ENDIF
      
      READ(FU) I, J, K
      IF (.NOT. ((I .EQ. FINFO%NRCHOP) .AND. (J .EQ. FINFO%NPCHOP) .AND. (K .EQ. FINFO%NZCHOP))) THEN
        IF (WARN_ON) WRITE(*,*) 'MLOAD: [WARNING] CHOPPING LOCATION(S) INCONSISTENT.'
        IF (WARN_ON) WRITE(*,*) 'R.CHOP: FILE --> ', I, ', PROGRAM -->', FINFO%NRCHOP
        IF (WARN_ON) WRITE(*,*) 'P.CHOP: FILE --> ', J, ', PROGRAM -->', FINFO%NPCHOP
        IF (WARN_ON) WRITE(*,*) 'Z.CHOP: FILE --> ', K, ', PROGRAM -->', FINFO%NZCHOP
      ENDIF

      READ(FU) I, J, K
      A%NRC = I
      A%NPC = J
      A%NZC = K
      CALL CHOPSET_VECTOR(A, I-FINFO%NRCHOP, J-FINFO%NPCHOP, K-FINFO%NZCHOP)

      READ(FU) Z, L
      IF (.NOT. ((Z .EQ. FINFO%ZLEN) .AND. (L .EQ. FINFO%ELL))) THEN
        IF (WARN_ON) WRITE(*,*) 'MLOAD: [WARNING] ZLEN & ELL INCONSISTENT.'
        IF (WARN_ON) WRITE(*,*) 'ZLEN: FILE --> ', Z, ', PROGRAM -->', FINFO%ZLEN
        IF (WARN_ON) WRITE(*,*) 'ELL: FILE --> ', L, ', PROGRAM -->', FINFO%ELL
      ENDIF

      READ(FU) (((A%ER(I,J,K), I=1, NRDIM), J=1, FINFO%NPDIM), K=1, FINFO%NZDIM)
      READ(FU) (((A%EP(I,J,K), I=1, NRDIM), J=1, FINFO%NPDIM), K=1, FINFO%NZDIM)
      READ(FU) (((A%EZ(I,J,K), I=1, NRDIM), J=1, FINFO%NPDIM), K=1, FINFO%NZDIM)

      CLOSE(FU)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
!  PRIVATE PROCEDURES ======================================================================================
! ==========================================================================================================
      SUBROUTINE HORFFT(AE, NPC, IS)
! ==========================================================================================================
! [USAGE]:
! PERFORM HORIZONTAL (AZIMUTHAL DIR.) FFT WITH RESPECT TO THE 2ND AXE OF A 3D ARRAY AE
! [VARIABLES]:
! A >> 3D ARRAY CONTAINING EITHER SCALAR OR ONE COMPONENT OF A VECTOR FIELD. (DIM1, DIM2, DIM3)=(R, PHI, Z)
! NPC >> CHOP INDEX IN P (2ND AXE)
! IS >> FORWARD OR BACKWARD (1: BACKWARD (F TO P), 0: FORWARD (P TO F))
! [NOTES]:
! 1. WHEN GOING FROM FFF SPACE TO PPP SPACE, EXPECTED CALLING SEQUENCE =
!    CALL VERFFT(A,1) (FFF -> FFP)
!    CALL RTRAN(A,1) (FFP -> PFP)
!    CALL HORFFT(A,1) (PFP -> PPP)
! 2. WHEN GOING FROM PPP SPACE TO FFF SPACE, EXPECTED CALLING SEQUENCE =
!    CALL HORFFT(A,0) (PPP -> PFP)
!    CALL RTRAN(A,0) (PFP -> FFP)
!    CALL VERFFT(A,0) (FFP -> FFF)
! FOR FFT, THE EXTERNAL FFTW3 LIBRARY IS UTILIZED.
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      INCLUDE 'fftw3.f03'

      COMPLEX(P8), DIMENSION(:,:,:), ALLOCATABLE           :: AE
      INTEGER, INTENT(IN)                                  :: NPC
      INTEGER, INTENT(IN)                                  :: IS

      COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), ALLOCATABLE :: B
      REAL(C_DOUBLE), DIMENSION(:), ALLOCATABLE            :: C
      INTEGER                                              :: II, KK, MM
      TYPE(C_PTR)                                          :: PLAN

      IF (NPC .NE. FINFO%NPCHOP) THEN
        IF (WARN_ON) WRITE(*,*) 'HORFFT: [WARNING] GLOBAL CHOP .NE. LOCAL CHOP IN PHI'
      ENDIF

      IF (IS .EQ. 0) THEN
        ALLOCATE( B(FINFO%NPDIM) )
        ALLOCATE( C(FINFO%NPH*2) )
        B = AE(1, :FINFO%NPDIM, 1)
        DO II = 1, FINFO%NPH
          C(2*II-1) = REAL(B(II))
          C(2*II) = AIMAG(B(II))
        ENDDO
        PLAN = FFTW_PLAN_DFT_R2C_1D(2*FINFO%NPH, C, B, FFTW_ESTIMATE)

        B = 0.D0
        C = 0.D0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(II, KK, MM, B, C)
        DO II = 1, FINFO%NR
          DO KK = 1, FINFO%NZ
            B = AE(II, :FINFO%NPDIM, KK)
            C = 0.D0
            DO MM = 1, FINFO%NPH
              C(2*MM-1) = REAL(B(MM))
              C(2*MM) = AIMAG(B(MM))
            ENDDO
            CALL FFTW_EXECUTE_DFT_R2C(PLAN, C, B)
            AE(II, :FINFO%NPDIM, KK) = B/(2*FINFO%NPH)
          ENDDO
        ENDDO
!$OMP END PARALLEL DO

        IF (.NOT. ALL(ABS(AE(:,NPC+1:,:)) .LT. 1.D-14)) THEN
          IF (WARN_ON) WRITE(*,*) 'HORFFT/F: [WARNING] RESIDUES IN THE REGION TO BE CHOPPED.'
        ENDIF

        DEALLOCATE( B )
        DEALLOCATE( C )

        CALL FFTW_DESTROY_PLAN(PLAN)
      ELSE
        IF ((.NOT. (ALL(AIMAG(AE(:,1,:)) .LT. 1.D-14) .AND. ALL(AIMAG(AE(:,NPC,:)) .LT. 1.D-14))) .OR. &
            (.NOT. ALL(ABS(AE(:,NPC+1:,:)) .LT. 1.D-14))) THEN
          IF (WARN_ON) WRITE(*,*) 'HORFFT/B: [WARNING] UNNECESSARY NONZERO ENTRIES WILL LEAD TO RESIDUES.'
        ENDIF

        ALLOCATE( B(FINFO%NPDIM) )
        ALLOCATE( C(FINFO%NPH*2) )
        B = AE(1, :FINFO%NPDIM, 1)
        C = 0.D0
        PLAN = FFTW_PLAN_DFT_C2R_1D(2*FINFO%NPH, B, C, FFTW_ESTIMATE)

        B = 0.D0
        C = 0.D0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(II, KK, MM, B, C)
        DO II = 1, FINFO%NR
          DO KK = 1, FINFO%NZ
            B = AE(II, :FINFO%NPDIM, KK)
            C = 0.D0
            CALL FFTW_EXECUTE_DFT_C2R(PLAN, B, C)
            DO MM = 1, FINFO%NPH
              AE(II, MM, KK) = CMPLX(C(2*MM-1), C(2*MM), P8)
            ENDDO
            AE(II, FINFO%NPH+1:, KK) = 0.D0
          ENDDO
        ENDDO
!$OMP END PARALLEL DO

        DEALLOCATE( B )
        DEALLOCATE( C )

        CALL FFTW_DESTROY_PLAN(PLAN)
      ENDIF

      WHERE (ABS(REAL(AE)) .LT. 1.D-14) AE = CMPLX(0.D0, AIMAG(AE), P8)
      WHERE (ABS(AIMAG(AE)) .LT. 1.D-14) AE = CMPLX(REAL(AE), 0.D0, P8)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE RTRANS(AE, NRCS, NPC, IS)
! ==========================================================================================================
! [USAGE]:
! PERFORM RADIAL LEGENDRE TRANSFORM WITH RESPECT TO THE 1ST AXE OF A 3D ARRAY AE
! [VARIABLES]:
! A >> 3D ARRAY CONTAINING EITHER SCALAR OR ONE COMPONENT OF A VECTOR FIELD. (DIM1, DIM2, DIM3)=(R, PHI, Z)
! NRCS >> CHOP INDEX IN R (1ST AXE)
! NPC >> CHOP INDEX IN P (2ND AXE)
! IS >> FORWARD OR BACKWARD (1: BACKWARD (F TO P), 0: FORWARD (P TO F))
! [NOTES]:
! 1. WHEN GOING FROM FFF SPACE TO PPP SPACE, EXPECTED CALLING SEQUENCE =
!    CALL VERFFT(A,1) (FFF -> FFP)
!    CALL RTRAN(A,1) (FFP -> PFP)
!    CALL HORFFT(A,1) (PFP -> PPP)
! 2. WHEN GOING FROM PPP SPACE TO FFF SPACE, EXPECTED CALLING SEQUENCE =
!    CALL HORFFT(A,0) (PPP -> PFP)
!    CALL RTRAN(A,0) (PFP -> FFP)
!    CALL VERFFT(A,0) (FFP -> FFF)
! FOR FFT, THE EXTERNAL FFTW3 LIBRARY IS UTILIZED.
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:,:,:), ALLOCATABLE     :: AE
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: NRCS 
      INTEGER, INTENT(IN)                            :: NPC
      INTEGER, INTENT(IN)                            :: IS

      COMPLEX(P8), DIMENSION(:,:,:), ALLOCATABLE     :: A
      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE       :: BE, BO
      INTEGER                                        :: I, NN, MM

      IF (NRCS(1) .NE. FINFO%NRCHOP) THEN
        IF (WARN_ON) WRITE(*,*) 'RTRANS: [WARNING] GLOBAL CHOP .NE. LOCAL CHOP IN R'
      ENDIF
      IF (NPC .NE. FINFO%NPCHOP) THEN
        IF (WARN_ON) WRITE(*,*) 'RTRANS: [WARNING] GLOBAL CHOP .NE. LOCAL CHOP IN PHI'
      ENDIF

      ALLOCATE( A(SIZE(AE,1), SIZE(AE,2), SIZE(AE,3)) )
      A = 0.D0

      IF( IS .EQ. 0 ) THEN
        ALLOCATE( BE(FINFO%NRH, FINFO%NZDIM) )
        ALLOCATE( BO(FINFO%NRH, FINFO%NZDIM) )

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, BE, BO)
        DO MM = 1, NPC
          BE = 0.D0
          BO = 0.D0

          DO I = 1, FINFO%NRH
            BE(I,:) = (AE(I,MM,:)+AE(FINFO%NR-I+1,MM,:))*TFM%W(I)
            BO(I,:) = (AE(I,MM,:)-AE(FINFO%NR-I+1,MM,:))*TFM%W(I)
          ENDDO

          NN = NRCS(MM)
          IF (NN .GE. 1) THEN
            A(1:NN:2, MM, :) = TRANSPOSE(TFM%PF(:FINFO%NRH, 1:NN:2, MM)) .MUL. BE
          ENDIF
          IF (NN .GE. 2) THEN
            A(2:NN:2, MM, :) = TRANSPOSE(TFM%PF(:FINFO%NRH, 2:NN:2, MM)) .MUL. BO
          ENDIF
        ENDDO
!$OMP END PARALLEL DO

        AE = A

        DEALLOCATE( BE )
        DEALLOCATE( BO )
      ELSE
        ALLOCATE( BE(FINFO%NRH, FINFO%NZDIM) )
        ALLOCATE( BO(FINFO%NRH, FINFO%NZDIM) )

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, BE, BO)
        DO MM = 1, NPC
          BE = 0.D0
          BO = 0.D0

          NN = NRCS(MM)
          IF (NN .GE. 1) THEN
            BE = TFM%PF(:FINFO%NRH, 1:NN:2, MM) .MUL. AE(1:NN:2, MM, :)
          ELSE
            BE = 0.D0
          ENDIF
          IF (NN .GE. 2) THEN
            BO = TFM%PF(:FINFO%NRH, 2:NN:2, MM) .MUL. AE(2:NN:2, MM, :)
          ELSE
            BO = 0.D0
          ENDIF

          A(1:FINFO%NRH, MM, :) = BE+BO
          A(FINFO%NR:FINFO%NRH+1:-1, MM, :) = BE-BO
        ENDDO
!$OMP END PARALLEL DO

        AE = A

        DEALLOCATE( BE )
        DEALLOCATE( BO )
      ENDIF

      DEALLOCATE( A )

      WHERE (ABS(REAL(AE)) .LT. 1.D-14) AE = CMPLX(0.D0, AIMAG(AE), P8)
      WHERE (ABS(AIMAG(AE)) .LT. 1.D-14) AE = CMPLX(REAL(AE), 0.D0, P8)

      RETURN
      END SUBROUTINE 
! ==========================================================================================================
      SUBROUTINE VERFFT(AE, NZC, NRCS, NPC, IS)
! ==========================================================================================================
! [USAGE]:
! PERFORM VERTICAL (AXIAL DIR.) FFT WITH RESPECT TO THE 3RD AXE OF A 3D ARRAY AE
! [VARIABLES]:
! A >> 3D ARRAY CONTAINING EITHER SCALAR OR ONE COMPONENT OF A VECTOR FIELD. (DIM1, DIM2, DIM3)=(R, PHI, Z)
! NZC >> CHOP INDEX IN Z (3RD AXE)
! NRCS >> CHOP INDEX IN R (1ST AXE)
! NPC >> CHOP INDEX IN P (2ND AXE)
! IS >> FORWARD OR BACKWARD (1: BACKWARD (F TO P), 0: FORWARD (P TO F))
! [NOTES]:
! 1. WHEN GOING FROM FFF SPACE TO PPP SPACE, EXPECTED CALLING SEQUENCE =
!    CALL VERFFT(A,1) (FFF -> FFP)
!    CALL RTRAN(A,1) (FFP -> PFP)
!    CALL HORFFT(A,1) (PFP -> PPP)
! 2. WHEN GOING FROM PPP SPACE TO FFF SPACE, EXPECTED CALLING SEQUENCE =
!    CALL HORFFT(A,0) (PPP -> PFP)
!    CALL RTRAN(A,0) (PFP -> FFP)
!    CALL VERFFT(A,0) (FFP -> FFF)
! FOR FFT, THE EXTERNAL FFTW3 LIBRARY IS UTILIZED.
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      INCLUDE 'fftw3.f03'

      COMPLEX(P8), DIMENSION(:,:,:), ALLOCATABLE           :: AE
      INTEGER, INTENT(IN)                                  :: NZC
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN)       :: NRCS 
      INTEGER, INTENT(IN)                                  :: NPC
      INTEGER, INTENT(IN)                                  :: IS

      COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), ALLOCATABLE :: B
      INTEGER                                              :: II, JJ
      TYPE(C_PTR)                                          :: PLAN

      IF (NRCS(1) .NE. FINFO%NRCHOP) THEN
        IF (WARN_ON) WRITE(*,*) 'VERFFT: [WARNING] GLOBAL CHOP .NE. LOCAL CHOP IN R'
      ENDIF
      IF (NPC .NE. FINFO%NPCHOP) THEN
        IF (WARN_ON) WRITE(*,*) 'VERFFT: [WARNING] GLOBAL CHOP .NE. LOCAL CHOP IN PHI'
      ENDIF

      IF (IS .EQ. 0) THEN
        ALLOCATE( B(FINFO%NZ) )
        B = AE(1, 1, :FINFO%NZ)
        PLAN = FFTW_PLAN_DFT_1D(FINFO%NZ, B, B, FFTW_FORWARD, FFTW_ESTIMATE)

        B = 0.D0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(JJ, II, B)
        DO JJ = 1, NPC
          DO II = 1, NRCS(JJ)
            B = AE(II, JJ, :FINFO%NZ)
            CALL FFTW_EXECUTE_DFT(PLAN, B, B)
            AE(II, JJ, :FINFO%NZ) = B/FINFO%NZ
          ENDDO
        ENDDO
!$OMP END PARALLEL DO

        IF (NZC .LE. FINFO%NZ-NZC+2) THEN
          IF (.NOT. ALL(ABS(AE(:,:,NZC+1:FINFO%NZ-NZC+1)) .LT. 1.D-14)) THEN
            IF (WARN_ON) WRITE(*,*) 'VERFFT/F: [WARNING] RESIDUES IN THE REGION TO BE CHOPPED.'
          ENDIF
        ENDIF

        DEALLOCATE( B )

        CALL FFTW_DESTROY_PLAN(PLAN)
      ELSE
        IF (NZC .LE. FINFO%NZ-NZC+2) THEN
          IF (.NOT. ALL(ABS(AE(:,:,NZC+1:FINFO%NZ-NZC+1)) .LT. 1.D-14)) THEN
            IF (WARN_ON) WRITE(*,*) 'VERFFT/B: [WARNING] UNNECESSARY NONZERO ENTRIES WILL LEAD TO RESIDUES.'
          ENDIF
        ENDIF

        ALLOCATE( B(FINFO%NZ) )
        B = AE(1, 1, :FINFO%NZ)
        PLAN = FFTW_PLAN_DFT_1D(FINFO%NZ, B, B, FFTW_BACKWARD, FFTW_ESTIMATE)

        B = 0.D0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(JJ, II, B)
        DO JJ = 1, NPC
          DO II = 1, NRCS(JJ)
            B = AE(II, JJ, :FINFO%NZ)
            CALL FFTW_EXECUTE_DFT(PLAN, B, B)
            AE(II, JJ, :FINFO%NZ) = B
          ENDDO
        ENDDO
!$OMP END PARALLEL DO

        DEALLOCATE( B )

        CALL FFTW_DESTROY_PLAN(PLAN)
      ENDIF

      WHERE (ABS(REAL(AE)) .LT. 1.D-14) AE = CMPLX(0.D0, AIMAG(AE), P8)
      WHERE (ABS(AIMAG(AE)) .LT. 1.D-14) AE = CMPLX(REAL(AE), 0.D0, P8)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
END MODULE