MODULE PARTICLEGEN ! MODULE FOR PARTICLE ARRAY GENERATION SUBROUTINES & FUNCTIONS
      USE OMP_LIB; USE MPI
      USE MISC; USE MATOPS
      IMPLICIT NONE
      PRIVATE
! ==========================================================================================================
!  PARAMS/TYPES ============================================================================================
! ==========================================================================================================
      TYPE :: PARTICLE_INFO
        REAL(P8) :: D, VOL, DENS, RLXT
        INTEGER :: NPARTS
      END TYPE
      TYPE(PARTICLE_INFO), PUBLIC :: PINFO

      TYPE, PUBLIC :: PARTICLE
        REAL(P8), DIMENSION(:,:), ALLOCATABLE :: POS
        REAL(P8), DIMENSION(:,:), ALLOCATABLE :: VEL
        REAL(P8), DIMENSION(:,:), ALLOCATABLE :: ACC
      END TYPE
! ==========================================================================================================
!  DECLARATIONS ============================================================================================
! ==========================================================================================================
      PUBLIC :: PARTICLE_INIT                                            ! GLOBAL PARTICLE INFORMATION SETUP
      PUBLIC :: PARTICLE_ALLOC, PARTICLE_DEALLOC                         ! ALLOCATE A PARTICLE-TYPE VARIABLE
      PUBLIC :: PARTICLE_SET                                   ! SET UP A PARTICLE'S KINETIC INFO IN (R,P,Z)
      PUBLIC :: PSAVE, PLOAD                                            ! SAVE/LOAD THE PARTICLE INFORMATION
! ==========================================================================================================
!  INTERFACES ==============================================================================================
! ==========================================================================================================
      INTERFACE PARTICLE_INIT
        MODULE PROCEDURE PARTICLE_INIT
      END INTERFACE

      INTERFACE PARTICLE_ALLOC
        MODULE PROCEDURE PARTICLE_ALLOC
      END INTERFACE

      INTERFACE PARTICLE_DEALLOC
        MODULE PROCEDURE PARTICLE_DEALLOC
      END INTERFACE

      INTERFACE PARTICLE_SET
        MODULE PROCEDURE PARTICLE_SET_INDIVIDUAL
        MODULE PROCEDURE PARTICLE_SET_BULK
      END INTERFACE

      INTERFACE PSAVE
        MODULE PROCEDURE PSAVE_MULTI
      END INTERFACE  

      INTERFACE PLOAD
        MODULE PROCEDURE PLOAD_MULTI
      END INTERFACE

CONTAINS
! ==========================================================================================================
!  PUBLIC PROCEDURES =======================================================================================
! ==========================================================================================================
      SUBROUTINE PARTICLE_INIT(DIAMIN, DENSIN, NPARTSIN, VISCIN)
! ==========================================================================================================
! [USAGE]:
! INITIALIZE THE GLOBAL FIELD (DOMAIN) SETUP, PROXYING FOU_INIT (IN FOURIER.MD) & LEG_INIT (IN LEGENDRE.MD)
! [VARIABLES]:
! DIAMIN >> REAL VALUE USED FOR D (EFF. DIAMETER OF A PARTICLE)
! DESNIN >> REAL VALUE USED FOR DENS (DENSITY OF A PARTICLE)
! NPARTSIN >> INTEGER VALUE USED FOR NP (# OF PARTICLES IN CONSIDERATION)
! VISCIN >> GENERAL VISCOSITY NU.
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)           :: NPARTSIN
      REAL(P8), INTENT(IN)          :: DIAMIN, DENSIN
      REAL(P8), OPTIONAL            :: VISCIN

      ! INTEGER                       :: I, J, N
      REAL(P8)                      :: VISCIN_

      VISCIN_  = 0.D0    ! DEFAULT: INVISCID

      IF (PRESENT(VISCIN)) VISCIN_ = VISCIN

      IF (VISCIN_ .LT. 0.D0) STOP 'PARTICLE_INIT: VISCOSITY CANNOT BE NEGATIVE.'
      IF (VISCIN_ .EQ. 0.D0) VISCIN_ = NEAREST(0.D0, 1.D0)

      PINFO%D = DIAMIN 
      PINFO%NPARTS = NPARTSIN
      PINFO%DENS = DENSIN
      PINFO%VOL = PI/6.D0*(DIAMIN**3.D0)
      PINFO%RLXT = DENSIN*(DIAMIN**2.D0)/VISCIN_/18.D0

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE PARTICLE_ALLOC(PTCL)
! ==========================================================================================================
! [USAGE]:
! ALLOCATE A PARTICLE-TYPE VARIABLE PTCL BASED ON THE GIVEN INITIALIZATION DATA
! [VARIABLES]:
! PTCL >> A PARTICLE-TYPE VARAIBLE TO BE ALLOCATED
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(PARTICLE), INTENT(INOUT)  :: PTCL

      IF (ALLOCATED(PTCL%POS)) THEN
        PTCL%POS = 0.D0
      ELSE
        ALLOCATE( PTCL%POS(3, PINFO%NPARTS) )
        PTCL%POS = 0.D0
      ENDIF

      IF (ALLOCATED(PTCL%VEL)) THEN
        PTCL%VEL = 0.D0
      ELSE
        ALLOCATE( PTCL%VEL(3, PINFO%NPARTS) )
        PTCL%VEL = 0.D0
      ENDIF

      IF (ALLOCATED(PTCL%ACC)) THEN
        PTCL%ACC = 0.D0
      ELSE
        ALLOCATE( PTCL%ACC(3, PINFO%NPARTS) )
        PTCL%ACC = 0.D0
      ENDIF

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE PARTICLE_DEALLOC(PTCL)
! ==========================================================================================================
! [USAGE]:
! DEALLOCATE A PARTICLE-TYPE VARIABLE PTCL 
! [VARIABLES]:
! PTCL >> A PARTICLE-TYPE VARAIBLE TO BE DEALLOCATED
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(PARTICLE), INTENT(INOUT)  :: PTCL

      IF (ALLOCATED(PTCL%POS)) DEALLOCATE(PTCL%POS)
      IF (ALLOCATED(PTCL%VEL)) DEALLOCATE(PTCL%VEL)
      IF (ALLOCATED(PTCL%ACC)) DEALLOCATE(PTCL%ACC)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE PARTICLE_SET_INDIVIDUAL(PTCL, N, POS, VEL)
! ==========================================================================================================
! [USAGE]:
! SET UP (INITIALIZE) A PARTICLE'S POSATION AND, OPTIONALLY, VELOCITY INDIVIDUALLY
! [VARIABLES]:
! N >> PARTICLE #
! POS >> LOCAION INFORMATION IN (R, P, Z)
! VEL >> (OPTIONAL) VELOCITY INFORMATION IN (R, P, Z). DEFAULT IS ZERO.
! [NOTES]:
! SANGJOON LEE @ JUNE 2023 
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(PARTICLE), INTENT(INOUT)    :: PTCL

      INTEGER, INTENT(IN)              :: N
      REAL(P8), DIMENSION(3)           :: POS
      REAL(P8), DIMENSION(3), OPTIONAL :: VEL

      REAL(P8), DIMENSION(3)           :: VEL_ = (/ 0.D0, 0.D0, 0.D0 /)

      IF (.NOT. ALLOCATED(PTCL%POS)) STOP 'PARTICLE_SET: AN INPUT PARTICLE VARIABLE NOT ALLOCATED.'

      IF (PRESENT(VEL)) VEL_ = VEL

      PTCL%POS(:,N) = POS
      PTCL%VEL(:,N) = VEL_

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE PARTICLE_SET_BULK(PTCL, TYPE, CENTER, CL, VEL)
! ==========================================================================================================
! [USAGE]:
! SET UP (INITIALIZE) A PARTICLE'S POSATION AND, OPTIONALLY, VELOCITY IN A BULK MANNER. IN A SPCEIFIED
! REGION, PARTICLES ARE RANDOMMLY SCATTERED.
! [VARIABLES]:
! TYPE >> PARTICLE DISTRIBUTION TYPE. 'B': BOX, 'C': CYLINDRICAL & 'S': SPHERICAL
! CENTER >> CENTER OF THE PARTICLE DISTRIBUTION (X, Y, Z)
! CL >> SCALING FACTORS. *(X_REGION,Y_REGION,Z_REGION)
! VEL >> (OPTIONAL) INITIAL BULK VELOCITY. DEFAULT IS ZERO.
! [NOTES]:
! SANGJOON LEE @ JUNE 2023 
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(PARTICLE), INTENT(INOUT)        :: PTCL

      CHARACTER(LEN=1), INTENT(IN)         :: TYPE
      REAL(P8), DIMENSION(3), INTENT(IN)   :: CENTER
      REAL(P8), DIMENSION(3), INTENT(IN)   :: CL
      REAL(P8), DIMENSION(3), OPTIONAL     :: VEL

      INTEGER                              :: N
      REAL(P8)                             :: X, Y, Z
      REAL(P8), DIMENSION(3)               :: VEL_ = (/ 0.D0, 0.D0, 0.D0 /)
      REAL(P8), DIMENSION(3, PINFO%NPARTS) :: PO

      IF (PRESENT(VEL)) VEL_ = VEL

      CALL PARTICLE_ALLOC(PTCL)

      CALL RANDOM_NUMBER( PO )

      IF (TYPE(1:1) .EQ. 'B') THEN
        PO = PO - 0.5D0
      ELSEIF (TYPE(1:1) .EQ. 'C') THEN
        PO(1,:) = PO(1,:)**(1.D0/2.D0) * 0.5D0
        PO(2,:) = PO(2,:) * 2.D0*PI
        PO(3,:) = PO(3,:) - 0.5D0
      ELSEIF (TYPE(1:1) .EQ. 'S') THEN
        PO(1,:) = PO(1,:)**(1.D0/3.D0) * 0.5D0
        PO(2,:) = PO(2,:) * PI
        PO(3,:) = PO(3,:) * 2.D0*PI
      ENDIF

      DO N = 1, PINFO%NPARTS
        IF (TYPE(1:1) .EQ. 'B') THEN
          X = PO(1,N)
          Y = PO(2,N)
          Z = PO(3,N)
        ELSEIF (TYPE(1:1) .EQ. 'C') THEN
          X = PO(1,N) * COS(PO(2,N))
          Y = PO(1,N) * SIN(PO(2,N)) 
          Z = PO(3,N)
        ELSEIF (TYPE(1:1) .EQ. 'S') THEN
          X = PO(1,N) * SIN(PO(2,N)) * COS(PO(3,N))
          Y = PO(1,N) * SIN(PO(2,N)) * SIN(PO(3,N))
          Z = PO(1,N) * COS(PO(2,N))
        ELSE
          STOP 'PARTICLE_SET: INVALID TYPE SETUP. TYPE MUST ONLY BE `B`OX, `C`YLINDER, OR `S`PHERE.'
        ENDIF

        PO(1,N) = X * CL(1) + CENTER(1)
        PO(2,N) = Y * CL(2) + CENTER(2)
        PO(3,N) = Z * CL(3) + CENTER(3)

        PTCL%POS(:,N) = PO(:,N)
        PTCL%VEL(:,N) = VEL_
      ENDDO

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE PSAVE_MULTI(PTCL, FN)
! ==========================================================================================================
! [USAGE]:
! SAVE MULTIPLE PARTICLES PS' INFORMATION INTO FN
! [VARIABLES]:
! PS >> ARRAY OF PARTICLES TO SAVE
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(PARTICLE), INTENT(INOUT) :: PTCL
      CHARACTER(LEN=*)              :: FN

      INTEGER                       :: STATUS, FU
      INTEGER                       :: I, N, NP, ND

      NP = PINFO%NPARTS
      ND = 3

      OPEN(UNIT=FU, FILE=TRIM(ADJUSTL(FN)), STATUS='UNKNOWN', FORM='UNFORMATTED', IOSTAT=STATUS)
      IF (STATUS .NE.0) STOP 'PSAVE: FAILED TO CREATE/OPEN A FILE.'

      WRITE(FU) 'P'                                ! SAVING PARTICLE ... 
      WRITE(FU) PINFO%D, PINFO%VOL                 ! PARTICLE DIAM. & VOL.
      WRITE(FU) PINFO%DENS, PINFO%RLXT            ! PARTICLE DENS. & RELAX. TIME
      WRITE(FU) NP, ND                             ! # OF PARTICLES (1 BC IT STORES A SINGLE PTCL)
      
      WRITE(FU) ((PTCL%POS(I,N), I=1, ND), N=1, NP)
      WRITE(FU) ((PTCL%VEL(I,N), I=1, ND), N=1, NP)
      WRITE(FU) ((PTCL%ACC(I,N), I=1, ND), N=1, NP)

      CLOSE(FU)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE PLOAD_MULTI(FN, PTCL)
! ==========================================================================================================
! [USAGE]:
! LOAD MULTIPLE PARTICLES' INFORMATION IN FN TO AN ARRAY OF PARTICLE-TYPE VARIABLE PS
! [VARIABLES]:
! FN >> FILENAME WHERE THE INFORMATION IS STORED
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      CHARACTER(LEN=*)              :: FN
      TYPE(PARTICLE), INTENT(INOUT) :: PTCL

      INTEGER                       :: STATUS, FU
      INTEGER                       :: I, N, NP, ND
      REAL(P8)                      :: D, V, DE, RT
      CHARACTER(LEN=1)              :: ISP

      OPEN(UNIT=FU, FILE=TRIM(ADJUSTL(FN)), STATUS='OLD', FORM='UNFORMATTED', IOSTAT=STATUS)
      IF (STATUS .NE.0) STOP 'PLOAD: FAILED TO CREATE/OPEN A FILE.'

      READ(FU) ISP
      IF (ISP .NE. 'P') STOP 'PLOAD: PARTICLE EXPECTED.'

      READ(FU) D, V
      IF (.NOT. ((D .EQ. PINFO%D) .AND. (V .EQ. PINFO%VOL))) THEN
        WRITE(*,*) 'PLOAD: PARTICLE INFO. INCOMPATIBLE.'
        WRITE(*,*) 'DIA.: FILE --> ', D, ', PROGRAM -->', PINFO%D
        WRITE(*,*) 'VOL.: FILE --> ', V, ', PROGRAM -->', PINFO%VOL
        STOP
      ENDIF
      
      READ(FU) DE, RT
      IF (.NOT. ((D .EQ. PINFO%D) .AND. (V .EQ. PINFO%VOL))) THEN
        WRITE(*,*) 'PLOAD: PARTICLE INFO. INCOMPATIBLE.'
        WRITE(*,*) 'DENSITY: FILE --> ', DE, ', PROGRAM -->', PINFO%DENS
        WRITE(*,*) 'RELAX.T: FILE --> ', RT, ', PROGRAM -->', PINFO%RLXT
        STOP
      ENDIF
      
      READ(FU) NP, ND
      IF (.NOT. (NP .EQ. PINFO%NPARTS)) THEN
        WRITE(*,*) 'PLOAD: PARTICLE INFO. INCOMPATIBLE.'
        WRITE(*,*) '# PARTICLES: FILE --> ', NP, ', PROGRAM -->', PINFO%NPARTS
        WRITE(*,*) 'DIMENSIONS : FILE --> ', ND, ', PROGRAM -->', 3
        STOP
      ENDIF
      
      IF (.NOT. ALLOCATED(PTCL%POS)) ALLOCATE(PTCL%POS(3, PINFO%NPARTS))
      IF (.NOT. ALLOCATED(PTCL%VEL)) ALLOCATE(PTCL%VEL(3, PINFO%NPARTS))
      IF (.NOT. ALLOCATED(PTCL%ACC)) ALLOCATE(PTCL%ACC(3, PINFO%NPARTS))
      READ(FU) ((PTCL%POS(I,N), I=1, ND), N=1, NP)
      READ(FU) ((PTCL%VEL(I,N), I=1, ND), N=1, NP)
      READ(FU) ((PTCL%ACC(I,N), I=1, ND), N=1, NP)

      CLOSE(FU)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
!  PRIVATE PROCEDURES ======================================================================================
! ==========================================================================================================

! ==========================================================================================================
END MODULE