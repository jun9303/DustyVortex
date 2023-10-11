MODULE PARTICLESTEP ! MODULE FOR 3D PARTICLE CONCENTRATION SCALAR FIELD TIME STEPPING SUBROUTINES & FUNCTIONS
      USE OMP_LIB; USE MPI
      USE MISC; USE MATOPS; USE TIMEINFO
      USE LEGENDRE; USE FOURIER; USE FIELDGEN; USE FIELDOPS; USE FIELDSTEP
      PRIVATE
! ==========================================================================================================
!  DECLARATIONS ============================================================================================
! ==========================================================================================================
      PUBLIC :: PNLTERM
      PUBLIC :: PFLDSTEP_EXPLICIT
      PUBLIC :: PFLDSTEP_BOOTSTR
      PUBLIC :: PFTERM
! ==========================================================================================================
!  INTERFACES ==============================================================================================
! ==========================================================================================================
      INTERFACE PNLTERM
        MODULE PROCEDURE PNLTERM_C
      END INTERFACE

      INTERFACE PFLDSTEP_EXPLICIT
        MODULE PROCEDURE FE
        MODULE PROCEDURE AB
      END INTERFACE

      INTERFACE PFLDSTEP_BOOTSTR
        MODULE PROCEDURE PRICH
      END INTERFACE

      INTERFACE PFTERM
        MODULE PROCEDURE PFTERM_C_U_UP
        MODULE PROCEDURE PFTERM_C_U_UP_UPP
      END INTERFACE
CONTAINS
! ==========================================================================================================
!  PUBLIC PROCEDURES =======================================================================================
! ==========================================================================================================
      SUBROUTINE PNLTERM_C(PSI, CHI, C, CNL, F)
! ==========================================================================================================
! [USAGE]:
! COMPUTE THE NONLINEAR TERM - U DOT (GRAD C) + ( 1 - BETA ) TAU (DIV F)
! [VARIABLES]:
! PSI >> TOROIDAL COMPONENT OF U (FLUID VELCOTIY FIELD) (IN FFF)
! CHI >> POLOIDAL COMPONENT OF U (FLUID VELOCITY FIELD) (IN FFF)
! C >> PARTICLE CONCENTRATION SCALAR FIELD (IN PPP)
! F >> (OPTIONAL) AN EXTERNAL FORCE VECTOR FIELD (IN PPP)
! CNL >> COMPUTED NONLINEAR TERM -U DOT (GRAD C)
! [NOTES]:
! SANGJOON LEE @ OCTOBER 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: PSI, CHI, C
      TYPE(SCALAR_FIELD), INTENT(OUT)   :: CNL
      TYPE(VECTOR_FIELD), OPTIONAL      :: F

      INTEGER                           :: MM, NN, KK
      REAL(P8)                          :: BETA
      TYPE(VECTOR_FIELD)                :: GRADC, U
      TYPE(SCALAR_FIELD)                :: UGRADC, DIVF

      CALL ALLOC(U, 'PPP')

      U = PT2VEL(PSI, CHI)

      BETA = 3.D0/(2.D0*FINFO%PTCDENS+1.D0)

      GRADC = GRAD(C)

      CALL ALLOC( UGRADC, 'PPP' )
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, KK)
          DO NN = 1, FINFO%NR
            DO MM = 1, FINFO%NPH
              DO KK = 1, FINFO%NZ
                UGRADC%E(NN,MM,KK) = REAL(U%ER(NN,MM,KK))*REAL(GRADC%ER(NN,MM,KK)) &
                                     + REAL(U%EP(NN,MM,KK))*REAL(GRADC%EP(NN,MM,KK)) &
                                     + REAL(U%EZ(NN,MM,KK))*REAL(GRADC%EZ(NN,MM,KK)) &
                                     + AIMAG(U%ER(NN,MM,KK))*AIMAG(GRADC%ER(NN,MM,KK))*IU &
                                     + AIMAG(U%EP(NN,MM,KK))*AIMAG(GRADC%EP(NN,MM,KK))*IU &
                                     + AIMAG(U%EZ(NN,MM,KK))*AIMAG(GRADC%EZ(NN,MM,KK))*IU
              ENDDO
            ENDDO
          ENDDO 
!$OMP END PARALLEL DO

      CALL ALLOC( DIVF, 'PPP' )
      IF (PRESENT(F)) DIVF = DIVERG(F)

      CALL ALLOC( CNL, 'PPP' )
      CNL%E = -UGRADC%E + (1.D0 - BETA) * FINFO%PTCRESP * DIVF%E
      CNL%LN = 0.D0

      CALL DEALLOC( DIVF )
      CALL DEALLOC( UGRADC )
      CALL DEALLOC( GRADC )
      CALL DEALLOC( U )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE FE(C, PSI, CHI, F)
! ==========================================================================================================
! [USAGE]:
! EXPLICIT TIME STEPPING USING FORWARD EULER. GLOBALLY FIRST-ORDER IN TIME
! [VARIABLES]:
! C >> ON ENTRY, C AT THE CURRENT STEP. ON EXIT, C AT THE NEXT STEP. (IN PPP)
! PSI >> PSI AT THE CURRENT STEP. (IN FFF)
! CHI >> CHI AT THE CURRENT STEP. (IN FFF)
! F >> (OPTIONAL) AN EXTERNAL (NONLINEAR) FORCE FIELD AT THE CURRENT STEP.
! [NOTES]:
! FOWRARD EULER: F(K+1) = F(K) + DT * (NLTERM(K) + LTERM(K))
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT)   :: C, PSI, CHI
      TYPE(VECTOR_FIELD), OPTIONAL        :: F

      REAL(P8)                            :: DT
      TYPE(SCALAR_FIELD)                  :: CNL

      IF (TINFO%N - TINFO%NI .EQ. TINFO%TOTN) THEN
        DT = TINFO%DTF
      ELSE
        DT = TINFO%DT
      ENDIF

      IF (PRESENT(F)) THEN
        CALL PNLTERM_C(PSI, CHI, C, CNL, F)
      ELSE
        CALL PNLTERM_C(PSI, CHI, C, CNL)
      ENDIF

      C%E = C%E + DT * CNL%E
      C%LN = C%LN + DT * CNL%LN

      CALL DEALLOC( CNL )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE AB(C, CP, PSI, CHI, PSIP, CHIP, F, FP)
! ==========================================================================================================
! [USAGE]:
! EXPLICIT TIME STEPPING USING FORWARD EULER. GLOBALLY FIRST-ORDER IN TIME
! [VARIABLES]:
! C >> ON ENTRY, C AT THE CURRENT STEP. ON EXIT, C AT THE NEXT STEP. (IN PPP)
! CP >> ON ENTRY, C AT THE PREVIOUS STEP. ON EXIT, C AT THE CURENT STEP. (IN PPP)
! PSI >> PSI AT THE CURRENT STEP. (IN FFF)
! CHI >> CHI AT THE CURRENT STEP. (IN FFF)
! PSIP >> PSI AT THE PREVIOUS STEP. (IN FFF)
! CHIP >> CHI AT THE PREVIOUS STEP. (IN FFF)
! F >> (OPTIONAL) AN EXTERNAL (NONLINEAR) FORCE FIELD AT THE CURRENT STEP.
! FP >> (OPTIONAL) AN EXTERNAL (NONLINEAR) FORCE FIELD AT THE PREVIOUS STEP.
! [NOTES]:
! ADAMS-BASHFORTH: F(K+1) = F(K) + DT * [1.5D0*(NLTERM(K)+LTERM(K)) - 0.5D0*(NLTERM(K-1)+LTERM(K-1))]
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT)   :: C, CP, PSI, CHI, PSIP, CHIP
      TYPE(VECTOR_FIELD), OPTIONAL        :: F, FP

      REAL(P8)                            :: DT
      TYPE(SCALAR_FIELD)                  :: CNL, CNLP

      IF (XOR(PRESENT(F),PRESENT(FP))) STOP 'FLDSTEP_EXPLICIT: INSUFFICIENT FORCE FIELD INPUTS.'

      IF (TINFO%N - TINFO%NI .EQ. TINFO%TOTN) THEN
        DT = TINFO%DTF
      ELSE
        DT = TINFO%DT
      ENDIF

      IF (PRESENT(F)) THEN
        CALL PNLTERM_C(PSI, CHI, C, CNL, F)
      ELSE
        CALL PNLTERM_C(PSI, CHI, C, CNL)
      ENDIF

      IF (PRESENT(FP)) THEN
        CALL PNLTERM_C(PSIP, CHIP, CP, CNLP, FP)
      ELSE
        CALL PNLTERM_C(PSIP, CHIP, CP, CNLP)
      ENDIF

      CP = C

      C%E = C%E + DT * (1.5D0*CNL%E - 0.5D0*CNLP%E)
      C%LN = C%LN + DT * (1.5D0*CNL%LN - 0.5D0*CNLP%LN)

      CALL DEALLOC( CNL )
      CALL DEALLOC( CNLP )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE PRICH(PSI, CHI, C, N, ORDER)
! ==========================================================================================================
! [USAGE]:
! AT THE INITIAL STEP WHERE THE PREVIOUS STEP'S INFORMATION DOESN'T EXIST, DO A SINGLE-STEP EXTRAPOLATION
! [VARIABLES]:
! PSI >> ON ENTRY, THE INITIAL PSI. ON EXIT, ONE-STEP EXTRAPOLATED PSI AT T = TI+DT. ERR ~ O(DT^2)
! CHI >> ON ENTRY, THE INITIAL CHI. ON EXIT, ONE-STEP EXTRAPOLATED CHI AT T = TI+DT. ERR ~ O(DT^2)
! C >> ON ENTRY, THE INITIAL C. ON EXIT, ONE-STEP EXTRAPOLATED C AT T = TI+DT. ERR ~ O(DT^2)
! N >> (OPTIONAL) SUBSTEPS TO BE TAKEN. DEFAULT IS 2.
! ORDER >> (OPTIONAL). IF 1, USE 1ST ORDER RICH (== FORWARD EULER). IF 2, USE 2ND ORDER RICH (DEFAULT)
! [NOTES]:
! DESIRED VALUE: PSI_EXACT(TI+DT), CHI_EXACT(TI+DT), C_EXACT(TI+DT)
! KNOWN VALUE: PSI(TI), CHI(TI), C(TI), DT. ESTIMATION BY THE SEMI-IMPLICIT EULER STEPPING (FE1_BE1).
! USING N SUBSTEPS, I.E. DDT = DT/N, PSI_EXACT(TI+DT) - PSI_EULER_DDT(TI+DT) = C*DDT*DT + O(DT^2)        [1]
! IF USING 2N SUBSTEPS, PSI_EXACT(TI+DT) - PSI_EULER_DDT/2(TI+DT) = C*DDT*2*DT + O(DT^2)                 [2]
! RICHARDSON EXTP.: 2*[2]-[1] => PSI_EXACT(TI+DT)-(2*PSI_EULER_DDT/2 - PSI_EULER_DDT) = O(DT^2)
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT)   :: PSI, CHI, C
      INTEGER, OPTIONAL                   :: N, ORDER

      INTEGER                             :: N_ = 2
      INTEGER                             :: ORDER_ = 2
      INTEGER                             :: I, MM, NN, KK
      REAL(P8)                            :: DT, DDT, DDTH

      TYPE(SCALAR_FIELD)                  :: PSI1, CHI1, PSI2, CHI2, C1, C2, PSI1P, CHI1P, PSI2P, CHI2P
      TYPE(SCALAR_FIELD)                  :: PSIVIS, CHIVIS, PSIUXW, CHIUXW
      TYPE(VECTOR_FIELD)                  :: U, UXW, UXWSOL, VIS, F

      IF (PRESENT(N)) N_ = N
      IF (PRESENT(ORDER)) ORDER_ = ORDER

      IF (N_ .LE. 0) STOP 'PFLDSTEP_BOOTSTR: WRONG SUBSTEP INPUT (< 1).'
      IF (ORDER_ .LE. 0) STOP 'PFLDSTEP_BOOTSTR: TARGET ACCURACY ORDER MUST BE POSITIVE INTEGER.'
      IF (TINFO%N .NE. TINFO%NI) STOP 'PFLDSTEP_BOOTSTR: BOOTSTRAPPING UNNECESSARY EXCEPT THE INITIAL STEP.'

      DT = TINFO%DT
      DDT = DT/N_
      DDTH = 0.5D0*DDT

      TINFO%DT = DDT ! SUBTIME STEPPING
      PSI1 = PSI; CHI1 = CHI; C1 = C
      DO I = 1, N_
        IF (I .EQ. 1) THEN
          U = PT2VEL(PSI1, CHI1)
          UXW = PT2VOR(PSI1, CHI1) ! UXW NOW CONTAINS CURL(U)
          UXW = VPROD(U, UXW)      ! UXW NOW CONTAINS U X CURL(U) IN PPP
          CALL LTERM(PSI1, CHI1, PSIVIS, CHIVIS)
          VIS = PT2VEL(PSIVIS, CHIVIS) ! VIS NOW CONTAINS VISC*DEL2(U) IN PPP
          CALL VELPROJ(UXW,PSIUXW,CHIUXW,PSILN=0.D0,CHILN=0.D0)
          UXWSOL = PT2VEL(PSIUXW, CHIUXW) ! UXWSOL IS THE SOLENOIDAL PORTION OF UXW IN PPP
          F = GRADSQH(U) ! F NOW CONTAINS GRAD(U^2/2)

          F%ER = F%ER - (UXW%ER - UXWSOL%ER) + VIS%ER
          F%EP = F%EP - (UXW%EP - UXWSOL%EP) + VIS%EP
          F%EZ = F%EZ - (UXW%EZ - UXWSOL%EZ) + VIS%EZ

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, KK)
          DO NN = 1, FINFO%NR
            DO MM = 1, FINFO%NPH
              DO KK = 1, FINFO%NZ
                F%ER(NN,MM,KK) = REAL(C%E(NN,MM,KK))*REAL(F%ER(NN,MM,KK)) &
                                 + AIMAG(C%E(NN,MM,KK))*AIMAG(F%ER(NN,MM,KK))*IU
                F%EP(NN,MM,KK) = REAL(C%E(NN,MM,KK))*REAL(F%EP(NN,MM,KK)) &
                                 + AIMAG(C%E(NN,MM,KK))*AIMAG(F%EP(NN,MM,KK))*IU
                F%EZ(NN,MM,KK) = REAL(C%E(NN,MM,KK))*REAL(F%EZ(NN,MM,KK)) &
                                 + AIMAG(C%E(NN,MM,KK))*AIMAG(F%EZ(NN,MM,KK))*IU
                                 
              ENDDO
            ENDDO
          ENDDO 
!$OMP END PARALLEL DO

          CALL DEALLOC( UXWSOL ); CALL DEALLOC( UXW ); CALL DEALLOC( U ); CALL DEALLOC( VIS )
          CALL DEALLOC( PSIVIS ); CALL DEALLOC( CHIVIS ); CALL DEALLOC( PSIUXW ); CALL DEALLOC( CHIUXW )
        ELSE
          F = PFTERM_C_U_UP(C, PSI1, CHI1, PSI1P, CHI1P)
        ENDIF
        PSI1P = PSI1; CHI1P = CHI1
        CALL FE(C1, PSI1, CHI1, F) ! C1 ADVANCES
        F%ER = F%ER*(FINFO%PTCDENS-1.D0)
        F%EP = F%EP*(FINFO%PTCDENS-1.D0)
        F%EZ = F%EZ*(FINFO%PTCDENS-1.D0)
        CALL FLDSTEP_SEMIIMP(PSI1, CHI1, F) ! PSI1, CHI1 ADVANCES
      ENDDO

      IF (ORDER_ .EQ. 1) THEN
        ! PSI_EULER_DDT(TI+DT) ~ O(DT))
        PSI%E = PSI1%E
        PSI%LN = PSI1%LN
        CHI%E = CHI1%E
        CHI%LN = CHI1%LN
        C%E = C1%E
        C%LN = C1%LN
        GOTO 999
      ENDIF

      TINFO%DT = DDTH ! HALF-SUBTIME STEPPING
      PSI2 = PSI; CHI2 = CHI; C2 = C
      DO I = 1, 2*N_
        IF (I .EQ. 1) THEN
          U = PT2VEL(PSI2, CHI2)
          UXW = PT2VOR(PSI2, CHI2) ! UXW NOW CONTAINS CURL(U)
          UXW = VPROD(U, UXW)      ! UXW NOW CONTAINS U X CURL(U) IN PPP
          CALL LTERM(PSI2, CHI2, PSIVIS, CHIVIS)
          VIS = PT2VEL(PSIVIS, CHIVIS) ! VIS NOW CONTAINS VISC*DEL2(U) IN PPP
          CALL VELPROJ(UXW,PSIUXW,CHIUXW,PSILN=0.D0,CHILN=0.D0)
          UXWSOL = PT2VEL(PSIUXW, CHIUXW) ! UXWSOL IS THE SOLENOIDAL PORTION OF UXW IN PPP
          F = GRADSQH(U) ! F NOW CONTAINS GRAD(U^2/2)

          F%ER = F%ER - (UXW%ER - UXWSOL%ER) + VIS%ER
          F%EP = F%EP - (UXW%EP - UXWSOL%EP) + VIS%EP
          F%EZ = F%EZ - (UXW%EZ - UXWSOL%EZ) + VIS%EZ

          CALL DEALLOC( UXWSOL ); CALL DEALLOC( UXW ); CALL DEALLOC( U ); CALL DEALLOC( VIS )
          CALL DEALLOC( PSIVIS ); CALL DEALLOC( CHIVIS ); CALL DEALLOC( PSIUXW ); CALL DEALLOC( CHIUXW )
        ELSE
          F = PFTERM_C_U_UP(C, PSI2, CHI2, PSI2P, CHI2P)
        ENDIF
        PSI2P = PSI2; CHI2P = CHI2
        CALL FE(C2, PSI2, CHI2, F) ! C2 ADVANCES
        F%ER = F%ER*(FINFO%PTCDENS-1.D0)
        F%EP = F%EP*(FINFO%PTCDENS-1.D0)
        F%EZ = F%EZ*(FINFO%PTCDENS-1.D0)
        CALL FLDSTEP_SEMIIMP(PSI2, CHI2, F) ! PSI2, CHI2 ADVANCES
      ENDDO

      IF (ORDER_ .EQ. 2) THEN
        ! 2*PSI_EULER_DDT/2(TI+DT) - PSI_EULER_DDT(TI+DT) ~ O(DT^2))
        PSI%E = 2.D0*PSI2%E - PSI1%E
        PSI%LN = 2.D0*PSI2%LN - PSI1%LN
        CHI%E = 2.D0*CHI2%E - CHI1%E
        CHI%LN = 2.D0*CHI2%LN - CHI1%LN
        C%E = 2.D0*C2%E - C1%E
        C%LN = 2.D0*C2%LN - C1%LN
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
      FUNCTION PFTERM_C_U_UP(C, PSI, CHI, PSIP, CHIP) RESULT(F)
! ==========================================================================================================
! [USAGE]:
! COMPUTE THE FORCE COUPLING TERM F = C DOT (DU/DT + (U DOT DEL)U)
! [VARIABLES]:
! C >> PARTICLE CONCENTRATION SCALAR FIELD AT THE CURRENT STEP (IN PPP)
! PSI >> TOROIDAL COMPONENT OF U AT THE CURRENT STEP (IN FFF)
! CHI >> POLOIDAL COMPONENT OF U AT THE CURRENT STEP (IN FFF)
! PSIP >> TOROIDAL COMPONENT OF U AT THE PREVIOUS STEP (IN FFF)
! PSIP >> POLOIDAL COMPONENT OF U AT THE PREVIOUS STEP (IN FFF)
! F >> OUTPUT VECTOR FIELD AS C DOT (DU/DT + (U DOT DEL)U) (IN PPP)
! [NOTES]:
! SANGJOON LEE @ OCTOBER 2023
! ==========================================================================================================
      IMPLICIT NONE
      TYPE(SCALAR_FIELD), INTENT(INOUT) :: C, PSI, CHI, PSIP, CHIP
      TYPE(VECTOR_FIELD)                :: F

      INTEGER                           :: MM, NN, KK
      REAL(P8)                          :: DT
      TYPE(VECTOR_FIELD)                :: U, UP
      TYPE(VECTOR_FIELD)                :: DUDT, UXW, GRADUSQH

      IF (TINFO%N - TINFO%NI .EQ. TINFO%TOTN) THEN
        DT = TINFO%DTF
      ELSE
        DT = TINFO%DT
      ENDIF

      U = PT2VEL(PSI, CHI)
      UP = PT2VEL(PSIP, CHIP)

      CALL ALLOC(DUDT, 'PPP')
      DUDT%ER = (U%ER - UP%ER) / DT
      DUDT%EP = (U%EP - UP%EP) / DT
      DUDT%EZ = (U%EZ - UP%EZ) / DT

      GRADUSQH = GRADSQH(U)
      UXW = PT2VOR(PSI, CHI)
      UXW = VPROD(U, UXW)

      DUDT%ER = DUDT%ER + GRADUSQH%ER - UXW%ER
      DUDT%EP = DUDT%EP + GRADUSQH%EP - UXW%EP
      DUDT%EZ = DUDT%EZ + GRADUSQH%EZ - UXW%EZ

      CALL ALLOC(F, 'PPP')
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, KK)
        DO NN = 1, FINFO%NR
          DO MM = 1, FINFO%NPH
            DO KK = 1, FINFO%NZ
              F%ER(NN,MM,KK) = REAL(C%E(NN,MM,KK))*REAL(DUDT%ER(NN,MM,KK)) &
                               + AIMAG(C%E(NN,MM,KK))*AIMAG(DUDT%ER(NN,MM,KK))*IU
              F%EP(NN,MM,KK) = REAL(C%E(NN,MM,KK))*REAL(DUDT%EP(NN,MM,KK)) &
                               + AIMAG(C%E(NN,MM,KK))*AIMAG(DUDT%EP(NN,MM,KK))*IU
              F%EZ(NN,MM,KK) = REAL(C%E(NN,MM,KK))*REAL(DUDT%EZ(NN,MM,KK)) &
                               + AIMAG(C%E(NN,MM,KK))*AIMAG(DUDT%EZ(NN,MM,KK))*IU
                               
            ENDDO
          ENDDO
        ENDDO 
!$OMP END PARALLEL DO

      CALL DEALLOC( DUDT )
      CALL DEALLOC( UXW )
      CALL DEALLOC( GRADUSQH )
      CALL DEALLOC( U )
      CALL DEALLOC( UP )

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION PFTERM_C_U_UP_UPP(C, PSI, CHI, PSIP, CHIP, PSIPP, CHIPP) RESULT(F)
! ==========================================================================================================
! [USAGE]:
! COMPUTE THE FORCE COUPLING TERM F = C DOT (DU/DT + (U DOT DEL)U)
! [VARIABLES]:
! C >> PARTICLE CONCENTRATION SCALAR FIELD AT THE CURRENT STEP (IN PPP)
! PSI >> TOROIDAL COMPONENT OF U AT THE CURRENT STEP (IN FFF)
! CHI >> POLOIDAL COMPONENT OF U AT THE CURRENT STEP (IN FFF)
! PSIP >> TOROIDAL COMPONENT OF U AT THE PREVIOUS STEP (IN FFF)
! PSIP >> POLOIDAL COMPONENT OF U AT THE PREVIOUS STEP (IN FFF)
! PSIPP >> TOROIDAL COMPONENT OF U AT THE SECOND PREVIOUS STEP (IN FFF)
! PSIPP >> POLOIDAL COMPONENT OF U AT THE SECOND PREVIOUS STEP (IN FFF)
! F >> OUTPUT VECTOR FIELD AS C DOT (DU/DT + (U DOT DEL)U) (IN PPP)
! [NOTES]:
! SANGJOON LEE @ OCTOBER 2023
! ==========================================================================================================
      IMPLICIT NONE
      TYPE(SCALAR_FIELD), INTENT(INOUT) :: C, PSI, CHI, PSIP, CHIP, PSIPP, CHIPP
      TYPE(VECTOR_FIELD)                :: F

      INTEGER                           :: MM, NN, KK
      REAL(P8)                          :: DT
      TYPE(VECTOR_FIELD)                :: U, UP, UPP
      TYPE(VECTOR_FIELD)                :: DUDT, UXW, GRADUSQH

      IF (TINFO%N - TINFO%NI .EQ. TINFO%TOTN) THEN
        DT = TINFO%DTF
      ELSE
        DT = TINFO%DT
      ENDIF

      U = PT2VEL(PSI, CHI)
      UP = PT2VEL(PSIP, CHIP)
      UPP = PT2VEL(PSIPP, CHIPP)

      CALL ALLOC(DUDT, 'PPP')
      DUDT%ER = (1.5D0*U%ER - 2.D0*UP%ER + .5D0*UPP%ER) / DT
      DUDT%EP = (1.5D0*U%EP - 2.D0*UP%EP + .5D0*UPP%EP) / DT
      DUDT%EZ = (1.5D0*U%EZ - 2.D0*UP%EZ + .5D0*UPP%EZ) / DT

      GRADUSQH = GRADSQH(U)
      UXW = PT2VOR(PSI, CHI)
      UXW = VPROD(U, UXW)

      DUDT%ER = DUDT%ER + GRADUSQH%ER - UXW%ER
      DUDT%EP = DUDT%EP + GRADUSQH%EP - UXW%EP
      DUDT%EZ = DUDT%EZ + GRADUSQH%EZ - UXW%EZ

      CALL ALLOC(F, 'PPP')
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, KK)
        DO NN = 1, FINFO%NR
          DO MM = 1, FINFO%NPH
            DO KK = 1, FINFO%NZ
              F%ER(NN,MM,KK) = REAL(C%E(NN,MM,KK))*REAL(DUDT%ER(NN,MM,KK)) &
                               + AIMAG(C%E(NN,MM,KK))*AIMAG(DUDT%ER(NN,MM,KK))*IU
              F%EP(NN,MM,KK) = REAL(C%E(NN,MM,KK))*REAL(DUDT%EP(NN,MM,KK)) &
                               + AIMAG(C%E(NN,MM,KK))*AIMAG(DUDT%EP(NN,MM,KK))*IU
              F%EZ(NN,MM,KK) = REAL(C%E(NN,MM,KK))*REAL(DUDT%EZ(NN,MM,KK)) &
                               + AIMAG(C%E(NN,MM,KK))*AIMAG(DUDT%EZ(NN,MM,KK))*IU
                               
            ENDDO
          ENDDO
        ENDDO 
!$OMP END PARALLEL DO

      CALL DEALLOC( DUDT )
      CALL DEALLOC( UXW )
      CALL DEALLOC( GRADUSQH )
      CALL DEALLOC( U )
      CALL DEALLOC( UP )
      CALL DEALLOC( UPP )

      RETURN
      END FUNCTION
! ==========================================================================================================
END MODULE