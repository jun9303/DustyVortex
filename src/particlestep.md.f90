MODULE PARTICLESTEP ! MODULE FOR PARTICLE TIME STEPPING SUBROUTINES & FUNCTIONS
      USE MISC; USE MATOPS; USE TIMEINFO
      USE LEGENDRE; USE FOURIER; USE FIELDGEN; USE FIELDOPS; USE FIELDSTEP
      USE PARTICLEGEN; USE PARTICLEINTF
      PRIVATE
! ==========================================================================================================
!  DECLARATIONS ============================================================================================
! ==========================================================================================================
      PUBLIC :: PFFLDCOMP
      PUBLIC :: PTCSTEP_EXPLICIT
      PUBLIC :: PTCSTEP_BOOTSTR
! ==========================================================================================================
!  INTERFACES ==============================================================================================
! ==========================================================================================================
      INTERFACE PFFLDCOMP
        MODULE PROCEDURE PFFLDCOMP
      END INTERFACE

      INTERFACE PTCSTEP_EXPLICIT
        MODULE PROCEDURE PFE
        MODULE PROCEDURE PAB
      END INTERFACE

      INTERFACE PTCSTEP_BOOTSTR
        MODULE PROCEDURE PRICH
      END INTERFACE

CONTAINS
! ==========================================================================================================
!  PUBLIC PROCEDURES =======================================================================================
! ==========================================================================================================
      SUBROUTINE PFFLDCOMP(PTCL, VFLD, PFFLD, ITERM, EPS, PFFLD_EST)
! =========================================================================================================
! [USAGE]:
! COMPUTE THE PARTICLE-FLUID INTERACTION FORCE VECTOR FIELD.
! [VARIABLES]:
! VFLD >> VELOCITY VECTOR FIELD AT THE CURRENT STEP, KNOWN IN ADVANCE (PASSED FROM FIELD TIME STEPPING)
! ITERM >> MAX. NO OF ITERATION FOR THE SECANT METHOD
! EPS >> CONVERGENCE EPSILON FOR THE SECANT METHOD
! PFFLD >> COMPUTED PARTICLE-FLUID INTERACTION FORCE VECTOR FIELD
! PFFLD_EST >> (OPTIONAL) PARTICLE-FLUID INTERACTION FORCE VECTOR FIELD ESTIMATE (OFTEN FROM THE PREV. ONE)
! [NOTES]:
! 1. F_P = - VOLFRAC_P * [ -GRAD(P) + RE^(-1)*DEL2(U) ] + DENS_P * VOLFRAC_P / RLXT_P * [ U_P - U ]
!    ALSO, IT IS KNOWN THAT -GRAD(P) = GRAD(U^2/2) - POTENTIAL(U X W) - POTENTIAL(F_P)
! 2. COMBINING THE ABOVE EQUATIONS, WE OBTAIN F_P = G(F_P). IF WE LET H(F) == G(F) - F, WE CAN FIND A SOLN
!    ITERATIVELY USING THE SECANT METHOD. THAT IS,
!    GIVEN: F_(0), F_(1) == G(F_(0)), I = 1
!    ITERATION: F_(I+1) = F_(I) - H(F_(I))*( F_(I) - F(I-1) )/( H(F_(I)) - H(F_(I-1)) ); I++
!    STOP: IF ||F_(I+1)-F(I)|| < EPS, OUTPUT F_P = F_(I+1). OTHERWISE, IF I REACHES ITERM, F_P = F(ITERM+1)
! SANGJOON LEE @ JUNE 2023
! =========================================================================================================
      IMPLICIT NONE

      TYPE(PARTICLE), INTENT(INOUT)     :: PTCL
      TYPE(VECTOR_FIELD), INTENT(INOUT) :: VFLD
      TYPE(VECTOR_FIELD), INTENT(INOUT) :: PFFLD
      INTEGER                           :: ITERM
      REAL(P8)                          :: EPS
      TYPE(VECTOR_FIELD), OPTIONAL      :: PFFLD_EST
     
      INTEGER                           :: I
      REAL(P8)                          :: EPSR, EPSP, EPSZ
      TYPE(VECTOR_FIELD)                :: PFFLDP
      TYPE(VECTOR_FIELD)                :: X, XX, XXX

      I = 1

      CALL ALLOC( X, 'PPP' )
      CALL ALLOC( XX, 'PPP' )
      CALL ALLOC( XXX, 'PPP' )

      IF (PRESENT(PFFLD_EST)) THEN
        CALL PACCCOMP(PTCL, VFLD, PFFLD_EST)
      ELSE
        CALL PACCCOMP(PTCL, VFLD)
      ENDIF

      CALL PFORCEFLD(PTCL, X)

      CALL PACCCOMP(PTCL, VFLD, X)
      CALL PFORCEFLD(PTCL, XX)
      
      EPSR = MAX(MAXVAL(ABS(REAL((XX%ER-X%ER)))), MAXVAL(ABS(AIMAG((XX%ER-X%ER)))))
      EPSP = MAX(MAXVAL(ABS(REAL((XX%EP-X%EP)))), MAXVAL(ABS(AIMAG((XX%EP-X%EP)))))
      EPSZ = MAX(MAXVAL(ABS(REAL((XX%EZ-X%EZ)))), MAXVAL(ABS(AIMAG((XX%EZ-X%EZ)))))
      ! WRITE(*,101) I, EPSR, EPSP, EPSZ

      PFFLDP = XX
      IF ( (EPSR .LE. EPS) .AND. (EPSP .LE. EPS) .AND. (EPSZ .LE. EPS) ) THEN
        PFFLD = XX
      ELSE
        DO WHILE ( (I .LT. ITERM) .AND. ( (EPSR .GT. EPS) .OR. (EPSP .GT. EPS) .OR. (EPSZ .GT. EPS) ) )
          CALL PACCCOMP(PTCL, VFLD, XX)
          CALL PFORCEFLD(PTCL, PFFLD)

          XXX%ER = XX%ER
          WHERE ( ABS((PFFLD%ER-XX%ER)-(PFFLDP%ER-X%ER)) .GT. EPS ) &
          XXX%ER=XX%ER-(PFFLD%ER-XX%ER)*(XX%ER-X%ER)/((PFFLD%ER-XX%ER)-(PFFLDP%ER-X%ER))

          XXX%EP = XX%EP
          WHERE ( ABS((PFFLD%EP-XX%EP)-(PFFLDP%EP-X%EP)) .GT. EPS ) &
          XXX%EP=XX%EP-(PFFLD%EP-XX%EP)*(XX%EP-X%EP)/((PFFLD%EP-XX%EP)-(PFFLDP%EP-X%EP))

          XXX%EZ = XX%EZ
          WHERE ( ABS((PFFLD%EZ-XX%EZ)-(PFFLDP%EZ-X%EZ)) .GT. EPS ) &
          XXX%EZ=XX%EZ-(PFFLD%EZ-XX%EZ)*(XX%EZ-X%EZ)/((PFFLD%EZ-XX%EZ)-(PFFLDP%EZ-X%EZ))
          I = I+1

          EPSR = MAX(MAXVAL(ABS(REAL((XXX%ER-XX%ER)))), MAXVAL(ABS(AIMAG((XXX%ER-XX%ER)))))
          EPSP = MAX(MAXVAL(ABS(REAL((XXX%EP-XX%EP)))), MAXVAL(ABS(AIMAG((XXX%EP-XX%EP)))))
          EPSZ = MAX(MAXVAL(ABS(REAL((XXX%EZ-XX%EZ)))), MAXVAL(ABS(AIMAG((XXX%EZ-XX%EZ)))))
          ! WRITE(*,101) I, EPSR, EPSP, EPSZ

          PFFLDP = PFFLD

          X = XX
          XX = XXX

        ENDDO
        PFFLD = XXX
      ENDIF

 101  FORMAT('PARTICLE-FLUID FORCE COMP. ITER. #'I03,': EPSR=',E10.3,', EPSP=',E10.3,', EPSZ=',E10.3)

      CALL DEALLOC( PFFLDP ) 
      CALL DEALLOC( X )
      CALL DEALLOC( XX )
      CALL DEALLOC( XXX )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE PFE(PTCL)
! ==========================================================================================================
! [USAGE]:
! EXPLICIT TIME STEPPING USING FORWARD EULER. GLOBALLY FIRST-ORDER IN TIME
! [VARIABLES]:
! PTCLS >> A SET OF PARTICLES. ON ENTRY, AT THE CURRENT STEP. ON EXIT. AT THE NEXT STEP.
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(PARTICLE), INTENT(INOUT) :: PTCL

      REAL(P8)                      :: DT
      
      IF (TINFO%N - TINFO%NI .EQ. TINFO%TOTN) THEN
        DT = TINFO%DTF
      ELSE
        DT = TINFO%DT
      ENDIF

      PTCL%POS = PTCL%POS + DT*PTCL%VEL
      PTCL%VEL = PTCL%VEL + DT*PTCL%ACC

      PTCL%POS(3,:) = MODULO(PTCL%POS(3,:),FINFO%ZLEN)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE PAB(PTCL, PTCLP)
! ==========================================================================================================
! [USAGE]:
! EXPLICIT TIME STEPPING USING FORWARD EULER. GLOBALLY FIRST-ORDER IN TIME
! [VARIABLES]:
! PTCLS >> A SET OF PARTICLES. ON ENTRY, AT THE CURRENT STEP. ON EXIT. AT THE NEXT STEP.
! PTCLSP >> A SET OF PARTICLES IN THE PREVIOUS STEP. ON EXIT. AT THE CURRENT STEP.
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(PARTICLE), INTENT(INOUT) :: PTCL
      TYPE(PARTICLE), INTENT(INOUT) :: PTCLP

      REAL(P8)                      :: DT
      TYPE(PARTICLE)                :: PTCLO

      IF (TINFO%N - TINFO%NI .EQ. TINFO%TOTN) THEN
        DT = TINFO%DTF
      ELSE
        DT = TINFO%DT
      ENDIF

      PTCLO%POS = PTCL%POS
      PTCLO%VEL = PTCL%VEL
      PTCLO%ACC = PTCL%ACC

      PTCL%POS = PTCL%POS + DT*(1.5D0*PTCL%VEL - 0.5D0*PTCLP%VEL)
      PTCL%VEL = PTCL%VEL + DT*(1.5D0*PTCL%ACC - 0.5D0*PTCLP%ACC)

      PTCL%POS(3,:) = MODULO(PTCL%POS(3,:),FINFO%ZLEN)

      PTCLP%POS = PTCLO%POS
      PTCLP%VEL = PTCLO%VEL
      PTCLP%ACC = PTCLO%ACC

      CALL PARTICLE_DEALLOC( PTCLO )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE PRICH(PSI, CHI, PTCL, F, N, ORDER)
! ==========================================================================================================
! [USAGE]:
! AT THE INITIAL STEP WHERE THE PREVIOUS STEP'S INFORMATION DOESN'T EXIST, DO A SINGLE-STEP EXTRAPOLATION
! [VARIABLES]:
! PTCL >> ON ENTRY, THE INITIAL PTCL. ON EXIT. ONE-STEP EXTRAPOLATED PTCL AT T = TI+DT. ERR ~ O(DT^2)
! N >> (OPTIONAL) SUBSTEPS TO BE TAKEN. DEFAULT IS 2.
! ORDER >> (OPTIONAL). IF 1, USE 1ST ORDER RICH (== FORWARD EULER). IF 2, USE 2ND ORDER RICH (DEFAULT)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT)   :: PSI, CHI
      TYPE(PARTICLE), INTENT(INOUT)       :: PTCL
      TYPE(VECTOR_FIELD), INTENT(INOUT)   :: F
      INTEGER, OPTIONAL                   :: N, ORDER

      INTEGER                             :: N_ = 2
      INTEGER                             :: ORDER_ = 2
      INTEGER                             :: I
      REAL(P8)                            :: DT, DDT, DDTH

      TYPE(SCALAR_FIELD)                  :: PSI1, CHI1, PSI2, CHI2
      TYPE(PARTICLE)                      :: PTCL1, PTCL2
      TYPE(VECTOR_FIELD)                  :: U, F1, F2, FF

      IF (PRESENT(N)) N_ = N
      IF (PRESENT(ORDER)) ORDER_ = ORDER

      IF (N_ .LE. 0) STOP 'PTCSTEP_BOOTSTR: WRONG SUBSTEP INPUT (< 1).'
      IF (ORDER_ .LE. 0) STOP 'PTCSTEP_BOOTSTR: TARGET ACCURACY ORDER MUST BE POSITIVE INTEGER.'
      IF (TINFO%N .NE. TINFO%NI) STOP 'PTCSTEP_BOOTSTR: BOOTSTRAPPING UNNECESSARY EXCEPT THE INITIAL STEP.'

      DT = TINFO%DT
      DDT = DT/N_
      DDTH = 0.5D0*DDT

      TINFO%DT = DDT ! SUBTIME STEPPING
      PSI1 = PSI; CHI1 = CHI; PTCL1 = PTCL; F1 = F
      DO I = 1, N_
        CALL FLDSTEP_SEMIIMP(PSI1, CHI1, F1)
        U = PT2VEL(PSI1, CHI1)
        CALL PFE(PTCL1)
        FF = F1
        CALL PFFLDCOMP(PTCL1, U, F1, 30, 1.D-12, FF)
      ENDDO

      IF (ORDER_ .EQ. 1) THEN
        PSI%E = PSI1%E
        CHI%E = CHI1%E
        PTCL%POS = PTCL1%POS
        PTCL%VEL = PTCL1%VEL
        PTCL%ACC = PTCL1%ACC
        F%ER = F1%ER; F%EP = F1%EP; F%EZ = F1%EZ
        GOTO 999
      ENDIF

      TINFO%DT = DDTH ! HALF-SUBTIME STEPPING
      PSI2 = PSI; CHI2 = CHI; PTCL2 = PTCL; F2 = F
      DO I = 1, 2*N_
        CALL FLDSTEP_SEMIIMP(PSI2, CHI2, F2)
        U = PT2VEL(PSI2, CHI2)
        CALL PFE(PTCL2)
        FF = F2
        CALL PFFLDCOMP(PTCL2, U, F2, 30, 1.D-12, FF)
      ENDDO

      IF (ORDER_ .EQ. 2) THEN
        PSI%E = 2.D0*PSI2%E - PSI1%E
        CHI%E = 2.D0*CHI2%E - CHI1%E
        PTCL%POS = 2.D0*PTCL2%POS - PTCL1%POS
        PTCL%VEL = 2.D0*PTCL2%VEL - PTCL1%VEL
        PTCL%ACC = 2.D0*PTCL2%ACC - PTCL1%ACC
        F%ER = 2.D0*F2%ER - F1%ER; F%EP = 2.D0*F2%EP - F1%EP; F%EZ = 2.D0*F2%EZ - F1%EZ
        GOTO 999
      ENDIF

      IF (ORDER_ .GT. 2) STOP 'PTCSTEP_BOOTSTR: CURRENTLY ORDER = 1 OR 2 SUPPORTED ONLY.'

999   TINFO%DT = DT ! RECOVER TO ORIGINAL TIME STEPPING
      PTCL%POS(3,:) = MODULO(PTCL%POS(3,:),FINFO%ZLEN)

      CALL DEALLOC( PSI1 )
      CALL DEALLOC( CHI1 )
      CALL DEALLOC( PSI2 )
      CALL DEALLOC( CHI2 )
      CALL PARTICLE_DEALLOC( PTCL1 )
      CALL PARTICLE_DEALLOC( PTCL2 )
      CALL DEALLOC( U )
      CALL DEALLOC( F1 )
      CALL DEALLOC( F2 )
      CALL DEALLOC( FF )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
END MODULE
