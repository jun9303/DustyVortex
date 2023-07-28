MODULE FIELDSTEP ! MODULE FOR 3D SCALAR/VECTOR FIELD TIME STEPPING SUBROUTINES & FUNCTIONS
      USE MISC; USE MATOPS; USE TIMEINFO
      USE LEGENDRE; USE FOURIER; USE FIELDGEN; USE FIELDOPS
      PRIVATE
! ==========================================================================================================
!  DECLARATIONS ============================================================================================
! ==========================================================================================================
      PUBLIC :: NLTERM
      PUBLIC :: LTERM
      PUBLIC :: FLDSTEP_EXPLICIT
      PUBLIC :: FLDSTEP_SEMIIMP
      PUBLIC :: FLDSTEP_BOOTSTR
! ==========================================================================================================
!  INTERFACES ==============================================================================================
! ==========================================================================================================
      INTERFACE NLTERM
        MODULE PROCEDURE NLTERM_PT
        MODULE PROCEDURE NLTERM_U
      END INTERFACE

      INTERFACE LTERM 
        MODULE PROCEDURE LTERM_PT
        MODULE PROCEDURE LTERM_U
      END INTERFACE

      INTERFACE FLDSTEP_EXPLICIT
        MODULE PROCEDURE FEFE
        MODULE PROCEDURE ABAB
      END INTERFACE

      INTERFACE FLDSTEP_SEMIIMP
        MODULE PROCEDURE FEBE
        MODULE PROCEDURE ABBE_OR_ABCN
      END INTERFACE

      INTERFACE FLDSTEP_BOOTSTR
        MODULE PROCEDURE RICH
      END INTERFACE

CONTAINS
! ==========================================================================================================
!  PUBLIC PROCEDURES =======================================================================================
! ==========================================================================================================
      SUBROUTINE NLTERM_PT(PSI, CHI, PSIUXW, CHIUXW, F)
! ==========================================================================================================
! [USAGE]:
! COMPUTE THE NONLINEAR TERM PROJ(U X W) == (PSIUXW, CHIUXW)
! [VARIABLES]:
! PSI >> TOROIDAL COMPONENT OF U (VELCOTIY FIELD) (IN FFF)
! CHI >> POLOIDAL COMPONENT OF U (VELOCITY FIELD) (IN FFF)
! PSIUXW >> TOROIDAL COMPONENT OF THE NONLINEAR TERM U X CURL(U) (OR U X W) (IN FFF)
! CHIUXW >> POLOIDAL COMPONENT OF THE NONLINEAR TERM U X CURL(U) (OR U X W) (IN FFF)
! F >> (OPTIONAL) AN EXTERNAL (NONLINEAR) FORCE FIELD.
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: PSI, CHI
      TYPE(SCALAR_FIELD), INTENT(OUT)   :: PSIUXW, CHIUXW
      TYPE(VECTOR_FIELD), OPTIONAL      :: F

      TYPE(SCALAR_FIELD)                :: PSIF, CHIF
      TYPE(VECTOR_FIELD)                :: U

      CALL ALLOC(U, 'PPP')

      U = PT2VEL(PSI, CHI)
      CALL UXWPROJ(U, PSIUXW, CHIUXW)

      IF (PRESENT(F)) THEN
        CALL VELPROJ(F, PSIF, CHIF, PSILN=0.D0)
        PSIUXW%E = PSIUXW%E + PSIF%E
        CHIUXW%E = CHIUXW%E + CHIF%E

        CALL DEALLOC( PSIF )
        CALL DEALLOC( CHIF )
      ENDIF

      CALL DEALLOC( U )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE NLTERM_U(U, PSIUXW, CHIUXW, F)
! ==========================================================================================================
! [USAGE]:
! COMPUTE THE NONLINEAR TERM PROJ(U X W) == (PSIUXW, CHIUXW)
! [VARIABLES]:
! U >> VELOCITY VECTOR FIELD, WHICH MUST BE SOLENOIDAL. (IN PPP SPACE)
! PSIUXW >> TOROIDAL COMPONENT OF THE NONLINEAR TERM U X CURL(U) (OR U X W) (IN FFF)
! CHIUXW >> POLOIDAL COMPONENT OF THE NONLINEAR TERM U X CURL(U) (OR U X W) (IN FFF)
! F >> (OPTIONAL) AN EXTERNAL (NONLINEAR) FORCE FIELD.
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(INOUT) :: U
      TYPE(SCALAR_FIELD), INTENT(OUT)   :: PSIUXW, CHIUXW
      TYPE(VECTOR_FIELD), OPTIONAL      :: F

      TYPE(SCALAR_FIELD)                :: PSIF, CHIF

      CALL UXWPROJ(U, PSIUXW, CHIUXW)

      IF (PRESENT(F)) THEN
        CALL VELPROJ(F, PSIF, CHIF, PSILN=0.D0)
        PSIUXW%E = PSIUXW%E + PSIF%E
        CHIUXW%E = CHIUXW%E + CHIF%E

        CALL DEALLOC( PSIF )
        CALL DEALLOC( CHIF )        
      ENDIF

      CALL DEALLOC( U )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE LTERM_PT(PSI, CHI, PSIVIS, CHIVIS)
! ==========================================================================================================
! [USAGE]:
! COMPUTE THE INEAR TERM VISC(U) == (PSIVIS, CHIVIS)
! [VARIABLES]:
! PSI >> TOROIDAL COMPONENT OF U (VELCOTIY FIELD) (IN FFF)
! CHI >> POLOIDAL COMPONENT OF U (VELOCITY FIELD) (IN FFF)
! PSIUXW >> TOROIDAL COMPONENT OF THE NONLINEAR TERM U X CURL(U) (OR U X W) (IN FFF)
! CHIUXW >> POLOIDAL COMPONENT OF THE NONLINEAR TERM U X CURL(U) (OR U X W) (IN FFF)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: PSI, CHI
      TYPE(SCALAR_FIELD), INTENT(OUT)   :: PSIVIS, CHIVIS

      CALL ALLOC(PSIVIS, 'FFF')
      CALL ALLOC(CHIVIS, 'FFF')

      IF (FINFO%VISCPOW .EQ. 2) THEN ! ZERO HYPERVISCOSITY
        IF (FINFO%VISC .LT. 1.D-15) THEN ! INVISCID (ZERO VISCOSITY)
          PSIVIS%E = 0.D0
          CHIVIS%E = 0.D0
        ELSE ! VISCOUS (NONZERO VISCOSITY)
          PSIVIS = DEL2(PSI)
          PSIVIS%E = PSIVIS%E * FINFO%VISC
          CHIVIS = DEL2(CHI) 
          CHIVIS%E = CHIVIS%E * FINFO%VISC
        ENDIF
      ELSE ! NONZERO HYPERVISCOSITY
        PSIVIS = HELMP(PSI, FINFO%VISCPOW, 0.D0, FINFO%VISC/FINFO%VISCP)
        PSIVIS%E = PSIVIS%E * FINFO%VISCP
        CHIVIS = HELMP(CHI, FINFO%VISCPOW, 0.D0, FINFO%VISC/FINFO%VISCP)
        CHIVIS%E = CHIVIS%E * FINFO%VISCP
      ENDIF

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE LTERM_U(U, PSIVIS, CHIVIS, RUP00INF)
! ==========================================================================================================
! [USAGE]:
! COMPUTE THE INEAR TERM VISC(U) == (PSIVIS, CHIVIS)
! [VARIABLES]:
! U >> VELOCITY VECTOR FIELD, WHICH MUST BE SOLENOIDAL. (IN PPP SPACE)
! PSIVIS >> TOROIDAL COMPONENT OF THE LINEAR TERM VISC(U) (IN FFF)
! CHIVIS >> POLOIDAL COMPONENT OF THE LINEAR TERM VISC(U) (IN FFF)
! RUP00INF >> (OPTIONAL) LIMIT OF R*U_P AS R -> INF. IF NOT GIVEN, ESTIMATED AT THE DOMAIN END. 
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(INOUT) :: U
      TYPE(SCALAR_FIELD), INTENT(OUT)   :: PSIVIS, CHIVIS
      REAL(P8), OPTIONAL                :: RUP00INF

      TYPE(SCALAR_FIELD)                :: PSI, CHI

      IF (PRESENT(RUP00INF)) THEN
        CALL VELPROJ(U,PSI,CHI,PSILN=-.5D0*RUP00INF)
      ELSE
        CALL VELPROJ(U,PSI,CHI)
      ENDIF

      CALL LTERM_PT(PSI,CHI,PSIVIS,CHIVIS)

      CALL DEALLOC( PSI )
      CALL DEALLOC( CHI )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE FEFE(PSI, CHI, F)
! ==========================================================================================================
! [USAGE]:
! EXPLICIT TIME STEPPING USING FORWARD EULER. GLOBALLY FIRST-ORDER IN TIME
! [VARIABLES]:
! PSI >> ON ENTRY, PSI AT THE CURRENT STEP. ON EXIT, PSI AT THE NEXT STEP. (IN FFF)
! CHI >> ON ENTRY, CHI AT THE CURRENT STEP. ON EXIT, CHI AT THE NEXT STEP. (IN FFF)
! F >> (OPTIONAL) AN EXTERNAL (NONLINEAR) FORCE FIELD AT THE CURRENT STEP.
! [NOTES]:
! FOWRARD EULER: F(K+1) = F(K) + DT * (NLTERM(K) + LTERM(K))
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT)   :: PSI, CHI
      TYPE(VECTOR_FIELD), OPTIONAL        :: F

      REAL(P8)                            :: DT
      TYPE(SCALAR_FIELD)                  :: PSIVIS, CHIVIS, PSIUXW, CHIUXW

      IF (TINFO%N - TINFO%NI .EQ. TINFO%TOTN) THEN
        DT = TINFO%DTF
      ELSE
        DT = TINFO%DT
      ENDIF

      IF (PRESENT(F)) THEN
        CALL NLTERM_PT(PSI, CHI, PSIUXW, CHIUXW, F)
      ELSE
        CALL NLTERM_PT(PSI, CHI, PSIUXW, CHIUXW)
      ENDIF
      CALL LTERM_PT(PSI, CHI, PSIVIS, CHIVIS)

      PSI%E = PSI%E + DT * (PSIUXW%E + PSIVIS%E)
      PSI%LN = PSI%LN + DT * (PSIUXW%LN + PSIVIS%LN)
      CHI%E = CHI%E + DT * (CHIUXW%E + CHIVIS%E)
      CHI%LN = CHI%LN + DT * (CHIUXW%LN + CHIVIS%LN)

      CALL DEALLOC( PSIVIS )
      CALL DEALLOC( CHIVIS )
      CALL DEALLOC( PSIUXW )
      CALL DEALLOC( CHIUXW )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE ABAB(PSI, CHI, PSIP, CHIP, F, FP)
! ==========================================================================================================
! [USAGE]:
! IMPLICIT TIME STEPPING USING ADAMS-BASHFORTH. GLOBALLY SECOND-ORDER IN TIME
! [VARIABLES]:
! PSI >> ON ENTRY, PSI AT THE CURRENT STEP. ON EXIT, PSI AT THE NEXT STEP. (IN FFF)
! CHI >> ON ENTRY, CHI AT THE CURRENT STEP. ON EXIT, CHI AT THE NEXT STEP. (IN FFF)
! PSIP >> ON ENTRY, PSI AT THE PREVIOUS STEP. ON EXIT, PSI AT THE CURRENT STEP. (IN FFF)
! CHIP >> ON ENTRY, CHI AT THE PREVIOUS STEP. ON EXIT, CHI AT THE CURRENT STEP. (IN FFF)
! F >> (OPTIONAL) AN EXTERNAL (NONLINEAR) FORCE FIELD AT THE CURRENT STEP.
! FP >> (OPTIONAL) AN EXTERNAL (NONLINEAR) FORCE FIELD AT THE PREVIOUS STEP.
! [NOTES]:
! ADAMS-BASHFORTH: F(K+1) = F(K) + DT * [1.5D0*(NLTERM(K)+LTERM(K)) - 0.5D0*(NLTERM(K-1)+LTERM(K-1))]
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT)   :: PSI, CHI
      TYPE(SCALAR_FIELD), INTENT(INOUT)   :: PSIP, CHIP
      TYPE(VECTOR_FIELD), OPTIONAL        :: F, FP

      REAL(P8)                            :: DT
      TYPE(SCALAR_FIELD)                  :: PSIVIS, CHIVIS, PSIUXW, CHIUXW
      TYPE(SCALAR_FIELD)                  :: PSIVISP, CHIVISP, PSIUXWP, CHIUXWP

      IF (XOR(PRESENT(F),PRESENT(FP))) STOP 'FLDSTEP_EXPLICIT: INSUFFICIENT FORCE FIELD INPUTS.'

      IF (TINFO%N - TINFO%NI .EQ. TINFO%TOTN) THEN
        DT = TINFO%DTF
      ELSE
        DT = TINFO%DT
      ENDIF

      IF (PRESENT(F)) THEN
        CALL NLTERM_PT(PSI, CHI, PSIUXW, CHIUXW, F)
      ELSE
        CALL NLTERM_PT(PSI, CHI, PSIUXW, CHIUXW)
      ENDIF
      CALL LTERM_PT(PSI, CHI, PSIVIS, CHIVIS)

      IF (PRESENT(FP)) THEN
        CALL NLTERM_PT(PSIP, CHIP, PSIUXWP, CHIUXWP, FP)
      ELSE
        CALL NLTERM_PT(PSIP, CHIP, PSIUXWP, CHIUXWP)
      ENDIF
      CALL LTERM_PT(PSIP, CHIP, PSIVISP, CHIVISP)

      PSIP = PSI
      CHIP = CHI

      PSI%E = PSI%E + DT * (1.5D0*(PSIUXW%E + PSIVIS%E) - 0.5D0*(PSIUXWP%E + PSIVISP%E))
      PSI%LN = PSI%LN + DT * (1.5D0*(PSIUXW%LN + PSIVIS%LN) - 0.5D0*(PSIUXWP%LN + PSIVISP%LN))
      CHI%E = CHI%E + DT * (1.5D0*(CHIUXW%E + CHIVIS%E) - 0.5D0*(CHIUXWP%E + CHIVISP%E))
      CHI%LN = CHI%LN + DT * (1.5D0*(CHIUXW%LN + CHIVIS%LN) - 0.5D0*(CHIUXWP%LN + CHIVISP%LN))


      CALL DEALLOC( PSIVIS )
      CALL DEALLOC( CHIVIS )
      CALL DEALLOC( PSIUXW )
      CALL DEALLOC( CHIUXW )
      CALL DEALLOC( PSIVISP )
      CALL DEALLOC( CHIVISP )
      CALL DEALLOC( PSIUXWP )
      CALL DEALLOC( CHIUXWP )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE FEBE(PSI, CHI, F)
! ==========================================================================================================
! [USAGE]:
! SEMI-IMPLICIT TIME STEPPING USING FORWARD EULER FOR THE NONLINEAR TERMS AND BACKWARD EULER FOR THE LINEAR 
! TERMS.
! [VARAIBLES]:
! PSI >> ON ENTRY, PSI AT THE CURRENT STEP. ON EXIT, PSI AT THE NEXT STEP. (IN FFF)
! CHI >> ON ENTRY, CHI AT THE CURRENT STEP. ON EXIT, CHI AT THE NEXT STEP. (IN FFF)
! F >> (OPTIONAL) AN EXTERNAL (NONLINEAR) FORCE FIELD AT THE CURRENT STEP.
! [NOTES]:
! FORWARD EULER: F(K+1/2) = F(K) + DT * NLTERM(K)
! BACKWARD EULER: F(K+1) = F(K+1/2) + DT * LTERM(K+1)
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT)   :: PSI, CHI
      TYPE(VECTOR_FIELD), OPTIONAL        :: F

      REAL(P8)                            :: DT, A, B
      TYPE(SCALAR_FIELD)                  :: PSIH, CHIH, PSIUXW, CHIUXW

      IF (TINFO%N - TINFO%NI .EQ. TINFO%TOTN) THEN
        DT = TINFO%DTF
      ELSE
        DT = TINFO%DT
      ENDIF

      IF (PRESENT(F)) THEN
        CALL NLTERM_PT(PSI, CHI, PSIUXW, CHIUXW, F)
      ELSE
        CALL NLTERM_PT(PSI, CHI, PSIUXW, CHIUXW)
      ENDIF

      CALL ALLOC(PSIH, 'FFF')
      CALL ALLOC(CHIH, 'FFF')

      PSIH%E = PSI%E + DT * PSIUXW%E
      PSIH%LN = PSI%LN + DT * PSIUXW%LN
      CHIH%E = CHI%E + DT * CHIUXW%E
      CHIH%LN = CHI%LN + DT * CHIUXW%LN

      IF (FINFO%VISCPOW .EQ. 2) THEN ! ZERO HYPERVISCOSITY
        IF (FINFO%VISC .LT. 1.D-15) THEN ! F(K+1) = F(K+1/2)
          STOP 'FLDSTEP_SEMIIMP: INVISCID CASE, NO LINEAR TERMS IN RHS. USE EXPLICIT METHODS.'
        ELSE ! F(K+1) = (1-DT*VISC*DEL2)^(-1) F(K+1/2) 
             ! = IHELM( (-1/(DT*VISC))*F(K+1/2)), ALP=(-1/(DT*VISC)) )
          A = -1.D0/(DT*FINFO%VISC)
          PSIH%LN = A*PSIH%LN
          PSIH%E = A*PSIH%E
          PSI = IHELM(PSIH, ALP=A, ISTP=.TRUE.)
          CHIH%LN = A*CHIH%LN
          CHIH%E = A*CHIH%E
          CHI = IHELM(CHIH, ALP=A, ISTP=.TRUE.)
        ENDIF
      ELSE ! F(K+1) = (1-DT*VISC*DEL2-DT*VISCP*DELP)^(-1) F(K+1/2)
           ! = [(-DT*VISCP)*(DELP+(VISC/VISCP)*DEL2+(-1/(DT*VISCP)))]^(-1) F(K+1/2)
           ! = IHELMP( (-1/(DT*VISCP))*F(K+1/2)), ALP=(-1/(DT*VISCP)), BET=(VISC/VISCP) )
        A = -1.D0/(DT*FINFO%VISCP)
        B = FINFO%VISC/FINFO%VISCP
        PSIH%LN = A*PSIH%LN
        PSIH%E = A*PSIH%E
        PSI = IHELMP(PSIH, FINFO%VISCPOW, ALP=A, BET=B, ISTP=.TRUE.)
        CHIH%LN = A*CHIH%LN
        CHIH%E = A*CHIH%E
        CHI = IHELMP(CHIH, FINFO%VISCPOW, ALP=A, BET=B, ISTP=.TRUE.)
      ENDIF

      CALL DEALLOC( PSIH )
      CALL DEALLOC( CHIH )
      CALL DEALLOC( PSIUXW )
      CALL DEALLOC( CHIUXW )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE ABBE_OR_ABCN(PSI, CHI, PSIP, CHIP, ORDER, F, FP)
! ==========================================================================================================
! [USAGE]:
! SEMI-IMPLICIT TIME STEPPING USING ADAMS-BASHFORTH FOR THE NONLINEAR TERMS AND BACKWARD EULER OR
! CRANK-NICOLSON FOR THE LINEAR TERMS. 
! [VARIABLES]:
! PSI >> ON ENTRY, PSI AT THE CURRENT STEP. ON EXIT, PSI AT THE NEXT STEP. (IN FFF)
! CHI >> ON ENTRY, CHI AT THE CURRENT STEP. ON EXIT, CHI AT THE NEXT STEP. (IN FFF)
! PSIP >> ON ENTRY, PSI AT THE PREVIOUS STEP. ON EXIT, PSI AT THE CURRENT STEP. (IN FFF)
! CHIP >> ON ENTRY, CHI AT THE PREVIOUS STEP. ON EXIT, CHI AT THE CURRENT STEP. (IN FFF)
! ORDER >> (OPTIONAL). IF 1, USE BE1 (BACKWRAD EULER). IF 2, USE CN2 (CRANK-NICOLSON, DEFAULT)
! F >> (OPTIONAL) AN EXTERNAL (NONLINEAR) FORCE FIELD AT THE CURRENT STEP.
! FP >> (OPTIONAL) AN EXTERNAL (NONLINEAR) FORCE FIELD AT THE PREVIOUS STEP.
! [NOTES]:
! ADAMS-BASHFORTH: F(K+1/2) = F(K) + DT * [1.5D0*NLTERM(K) - 0.5D0*NLTERM(K-1)]
! BACKWARD EULER: F(K+1) = F(K+1/2) + DT * LTERM(K+1)
! CRANK-NICOLSON: F(K+1) = F(K+1/2) + DT/2 * [LTERM(K+1) + LTERM(K)]
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT)   :: PSI, CHI
      TYPE(SCALAR_FIELD), INTENT(INOUT)   :: PSIP, CHIP
      INTEGER, OPTIONAL                   :: ORDER
      TYPE(VECTOR_FIELD), OPTIONAL        :: F, FP

      INTEGER                             :: ORDER_ = 2
      REAL(P8)                            :: DT, A, B
      TYPE(SCALAR_FIELD)                  :: PSIH, CHIH, PSIUXW, CHIUXW
      TYPE(SCALAR_FIELD)                  :: PSIUXWP, CHIUXWP, PSIVISP, CHIVISP

      IF (PRESENT(ORDER)) ORDER_ = ORDER

      IF (ORDER_ .LE. 0) STOP 'FLDSTEP_SEMIIMP: GLOBAL ACCURACY ORDER MUST BE POSITIVE INTEGER.'

      IF (XOR(PRESENT(F),PRESENT(FP))) STOP 'FLDSTEP_SEMIIMP: INSUFFICIENT FORCE FIELD INPUTS.'

      IF (TINFO%N - TINFO%NI .EQ. TINFO%TOTN) THEN
        DT = TINFO%DTF
      ELSE
        DT = TINFO%DT
      ENDIF

      IF (PRESENT(F)) THEN
        CALL NLTERM_PT(PSI, CHI, PSIUXW, CHIUXW, F)
      ELSE
        CALL NLTERM_PT(PSI, CHI, PSIUXW, CHIUXW)
      ENDIF

      IF (PRESENT(FP)) THEN
        CALL NLTERM_PT(PSIP, CHIP, PSIUXWP, CHIUXWP, FP)
      ELSE
        CALL NLTERM_PT(PSIP, CHIP, PSIUXWP, CHIUXWP)
      ENDIF

      IF (ORDER_ .GT. 1) CALL LTERM_PT(PSIP, CHIP, PSIVISP, CHIVISP)

      PSIP = PSI
      CHIP = CHI

      CALL ALLOC(PSIH, 'FFF')
      CALL ALLOC(CHIH, 'FFF')

      PSIH%E = PSI%E + DT * (1.5D0*PSIUXW%E - 0.5D0*PSIUXWP%E)
      PSIH%LN = PSI%LN + DT * (1.5D0*PSIUXW%LN - 0.5D0*PSIUXWP%LN)
      CHIH%E = CHI%E + DT * (1.5D0*CHIUXW%E - 0.5D0*CHIUXWP%E)
      CHIH%LN = CHI%LN + DT * (1.5D0*CHIUXW%LN - 0.5D0*CHIUXWP%LN)

      IF (ORDER_ .EQ. 1) THEN

        IF (FINFO%VISCPOW .EQ. 2) THEN ! ZERO HYPERVISCOSITY
          IF (FINFO%VISC .LT. 1.D-15) THEN ! F(K+1) = F(K+1/2)
            STOP 'FLDSTEP_SEMIIMP: INVISCID CASE, NO LINEAR TERMS IN RHS. USE EXPLICIT METHODS.'
          ELSE ! F(K+1) = (1-DT*VISC*DEL2)^(-1) F(K+1/2) 
               ! = IHELM( (-1/(DT*VISC))*F(K+1/2)), ALP=(-1/(DT*VISC)) )
            A = -1.D0/(DT*FINFO%VISC)
            PSIH%LN = A*PSIH%LN
            PSIH%E = A*PSIH%E
            PSI = IHELM(PSIH, ALP=A, ISTP=.TRUE.)
            CHIH%LN = A*CHIH%LN
            CHIH%E = A*CHIH%E
            CHI = IHELM(CHIH, ALP=A, ISTP=.TRUE.)
          ENDIF
        ELSE ! F(K+1) = (1-DT*VISC*DEL2-DT*VISCP*DELP)^(-1) F(K+1/2)
             ! = [(-DT*VISCP)*(DELP+(VISC/VISCP)*DEL2+(-1/(DT*VISCP)))]^(-1) F(K+1/2)
             ! = IHELMP( (-1/(DT*VISCP))*F(K+1/2)), ALP=(-1/(DT*VISCP)), BET=(VISC/VISCP) )
          A = -1.D0/(DT*FINFO%VISCP)
          B = FINFO%VISC/FINFO%VISCP
          PSIH%LN = A*PSIH%LN
          PSIH%E = A*PSIH%E
          PSI = IHELMP(PSIH, FINFO%VISCPOW, ALP=A, BET=B, ISTP=.TRUE.)
          CHIH%LN = A*CHIH%LN
          CHIH%E = A*CHIH%E
          CHI = IHELMP(CHIH, FINFO%VISCPOW, ALP=A, BET=B, ISTP=.TRUE.)
        ENDIF

      ELSE IF (ORDER_ .EQ. 2) THEN

        IF (FINFO%VISCPOW .EQ. 2) THEN ! ZERO HYPERVISCOSITY
          IF (FINFO%VISC .LT. 1.D-15) THEN ! F(K+1) = F(K+1/2)
            STOP 'FLDSTEP_SEMIIMP: INVISCID CASE, NO LINEAR TERMS IN RHS. USE EXPLICIT METHODS.'
          ELSE ! F(K+1) = (1-DT/2*VISC*DEL2)^(-1) [F(K+1/2) + DT/2*LTERM(K)] 
               ! = IHELM( (-2/(DT*VISC))*[F(K+1/2)+DT/2*LTERM(K)], ALP=(-2/(DT*VISC)) )
            A = -2.D0/(DT*FINFO%VISC)
            PSIH%LN = A*(PSIH%LN + DT/2.D0*PSIVISP%LN)
            PSIH%E = A*(PSIH%E + DT/2.D0*PSIVISP%E)
            PSI = IHELM(PSIH, ALP=A, ISTP=.TRUE.)
            CHIH%LN = A*(CHIH%LN + DT/2.D0*CHIVISP%LN)
            CHIH%E = A*(CHIH%E + DT/2.D0*CHIVISP%E)
            CHI = IHELM(CHIH, ALP=A, ISTP=.TRUE.)
          ENDIF
        ELSE ! F(K+1) = (1-DT/2*VISC*DEL2-DT/2*VISCP*DELP)^(-1) [F(K+1/2) + DT/2*LTERM(K)]
             ! = [(-DT/2*VISCP)*(DELP+(VISC/VISCP)*DEL2+(-2/(DT*VISCP)))]^(-1) [F(K+1/2) + DT/2*LTERM(K)]
             ! = IHELMP((-2/(DT*VISCP))*[F(K+1/2) + DT/2*LTERM(K)]), ALP=(-2/(DT*VISCP)), BET=(VISC/VISCP))
          A = -2.D0/(DT*FINFO%VISCP)
          B = FINFO%VISC/FINFO%VISCP
          PSIH%LN = A*(PSIH%LN + DT/2.D0*PSIVISP%LN)
          PSIH%E = A*(PSIH%E + DT/2.D0*PSIVISP%E)
          PSI = IHELMP(PSIH, FINFO%VISCPOW, ALP=A, BET=B, ISTP=.TRUE.)
          CHIH%LN = A*(CHIH%LN + DT/2.D0*CHIVISP%LN)
          CHIH%E = A*(CHIH%E + DT/2.D0*CHIVISP%E)
          CHI = IHELMP(CHIH, FINFO%VISCPOW, ALP=A, BET=B, ISTP=.TRUE.)
        ENDIF

      ENDIF

      IF (ORDER_ .GT. 2) STOP 'FLDSTEP_SEMIIMP: CURRENTLY ORDER = 1 (AB_BE) OR 2 (AB_CN) SUPPORTED ONLY.'

      CALL DEALLOC( PSIH )
      CALL DEALLOC( CHIH )
      CALL DEALLOC( PSIUXWP )
      CALL DEALLOC( CHIUXWP )
      CALL DEALLOC( PSIUXW )
      CALL DEALLOC( CHIUXW )
      IF (ORDER_ .GT. 1) CALL DEALLOC( PSIVISP )
      IF (ORDER_ .GT. 1) CALL DEALLOC( CHIVISP )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE RICH(PSI, CHI, N, ORDER)
! ==========================================================================================================
! [USAGE]:
! AT THE INITIAL STEP WHERE THE PREVIOUS STEP'S INFORMATION DOESN'T EXIST, DO A SINGLE-STEP EXTRAPOLATION
! [VARIABLES]:
! PSI >> ON ENTRY, THE INITIAL PSI. ON EXIT, ONE-STEP EXTRAPOLATED PSI AT T = TI+DT. ERR ~ O(DT^2)
! CHI >> ON ENTRY, THE INITIAL CHI. ON EXIT. ONE-STEP EXTRAPOLATED CHI AT T = TI+DT. ERR ~ O(DT^2)
! N >> (OPTIONAL) SUBSTEPS TO BE TAKEN. DEFAULT IS 2.
! ORDER >> (OPTIONAL). IF 1, USE 1ST ORDER RICH (== FORWARD EULER). IF 2, USE 2ND ORDER RICH (DEFAULT)
! [NOTES]:
! DESIRED VALUE: PSI_EXACT(TI+DT), CHI_EXACT(TI+DT)
! KNOWN VALUE: PSI(TI), CHI(TI), DT. ESTIMATION BY THE SEMI-IMPLICIT EULER STEPPING (FE1_BE1).
! USING N SUBSTEPS, I.E. DDT = DT/N, PSI_EXACT(TI+DT) - PSI_EULER_DDT(TI+DT) = C*DDT*DT + O(DT^2)        [1]
! IF USING 2N SUBSTEPS, PSI_EXACT(TI+DT) - PSI_EULER_DDT/2(TI+DT) = C*DDT*2*DT + O(DT^2)                 [2]
! RICHARDSON EXTP.: 2*[2]-[1] => PSI_EXACT(TI+DT)-(2*PSI_EULER_DDT/2 - PSI_EULER_DDT) = O(DT^2)
! * USED FOR PURE VORTEX (NON-DUSTY) RUN ONLY. IN CASE OF DUSTY SIMULATIONS, USE PRICH IN PARTICLESTEP.MD
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT)   :: PSI, CHI
      INTEGER, OPTIONAL                   :: N, ORDER

      INTEGER                             :: N_ = 2
      INTEGER                             :: ORDER_ = 2
      INTEGER                             :: I
      REAL(P8)                            :: DT, DDT, DDTH

      TYPE(SCALAR_FIELD)                  :: PSI1, CHI1, PSI2, CHI2

      IF (PRESENT(N)) N_ = N
      IF (PRESENT(ORDER)) ORDER_ = ORDER

      IF (N_ .LE. 0) STOP 'FLDSTEP_BOOTSTR: WRONG SUBSTEP INPUT (< 1).'
      IF (ORDER_ .LE. 0) STOP 'FLDSTEP_BOOTSTR: TARGET ACCURACY ORDER MUST BE POSITIVE INTEGER.'
      IF (TINFO%N .NE. TINFO%NI) STOP 'FLDSTEP_BOOTSTR: BOOTSTRAPPING UNNECESSARY EXCEPT THE INITIAL STEP.'

      DT = TINFO%DT
      DDT = DT/N_
      DDTH = 0.5D0*DDT

      TINFO%DT = DDT ! SUBTIME STEPPING
      PSI1 = PSI; CHI1 = CHI
      DO I = 1, N_
        CALL FEBE(PSI1, CHI1)
      ENDDO

      IF (ORDER_ .EQ. 1) THEN
        ! PSI_EULER_DDT(TI+DT) ~ O(DT))
        PSI%E = PSI1%E
        CHI%E = CHI1%E
        GOTO 999
      ENDIF

      TINFO%DT = DDTH ! HALF-SUBTIME STEPPING
      PSI2 = PSI; CHI2 = CHI
      DO I = 1, 2*N_
        CALL FEBE(PSI2, CHI2)
      ENDDO

      IF (ORDER_ .EQ. 2) THEN
        ! 2*PSI_EULER_DDT/2(TI+DT) - PSI_EULER_DDT(TI+DT) ~ O(DT^2))
        PSI%E = 2.D0*PSI2%E - PSI1%E
        CHI%E = 1.D0*CHI2%E - CHI1%E
        GOTO 999
      ENDIF

      IF (ORDER_ .GT. 2) STOP 'STEP_BOOT_STR: CURRENTLY ORDER = 1 OR 2 SUPPORTED ONLY.'

999   TINFO%DT = DT ! RECOVER TO ORIGINAL TIME STEPPING

      CALL DEALLOC( PSI1 )
      CALL DEALLOC( CHI1 )
      CALL DEALLOC( PSI2 )
      CALL DEALLOC( CHI2 )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
END MODULE