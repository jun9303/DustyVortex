MODULE FIELDOPS ! MODULE FOR 3D SCALAR/VECTOR FIELD OPERATION SUBROUTINES & FUNCTIONS
      USE OMP_LIB
      USE MISC; USE MATOPS
      USE LEGENDRE; USE FOURIER; USE FIELDGEN
      IMPLICIT NONE
      PRIVATE
! ==========================================================================================================
!  DECLARATIONS ============================================================================================
! ==========================================================================================================
      PUBLIC :: DELSQP ! B = (1-X)^(-2)*DEL2_PERP(A)
      PUBLIC :: IDELSQP ! B = DEL2_PERP^(-1)((1-X)^2*A)
      PUBLIC :: MULXP ! B = (1+X)*A <==> 2L^2 / (R^2 + L^2) *
      PUBLIC :: DIVXP ! B = (1+X)^(-1)*A <==> (R^2 + L^2) / 2L^2 * 
      PUBLIC :: MULXM ! B = (1-X)*A <==> 2R^2 / (R^2 + L^2) *
      PUBLIC :: DIVXM ! B = (1-X)^(-1)*A <==> (R^2 + L^2) / 2R^2 *
      PUBLIC :: XXDX  ! B = (1-X)^2*D(A)/DX <==> R*D()/DR
      PUBLIC :: DEL2 ! B = DEL2(A)
      PUBLIC :: IDEL2 ! B = DEL2^(-1)(A)
      PUBLIC :: HELM ! B = HELM(A)
      PUBLIC :: IHELM ! B = HELM^(-1)(A)
      PUBLIC :: HELMP ! B = (DEL^P + BET*DEL2 + ALP)*A
      PUBLIC :: IHELMP ! B = (DEL^P + BET*DEL2 + ALP)^(-1)*A
      PUBLIC :: VPROD ! C = A X B
      PUBLIC :: PT2VEL ! (PSI, CHI) ---> (U_R, U_P, U_Z)
      PUBLIC :: PT2VOR ! (PSI, CHI) ---> (O_R, O_P, O_Z) == CURL OF U
      PUBLIC :: PTFFRMV ! (PSI, CHI) ---> (PSI_FFRMV, CHIFFRMV) == FAR-FIELD VORTICITY IS REMOVED
      PUBLIC :: VELPROJ ! (U_R, U_P, U_Z) ---> (PSI, CHI)
      PUBLIC :: UXWPROJ ! (U_R, U_P, U_Z)  X (O_R, O_P, O_Z) ---> (PSI_UXW, CHI_UXW) == (PSIN, CHIN)
      PUBLIC :: VEL2GRADP ! (U_R, U_P, U_Z) ---> GRAD(PRESSURE) VEC. FIELD
      PUBLIC :: VEL2FORCE ! (U_R, U_P, U_Z) ---> -GRAD(P) + NU*VECDEL^2(U) (FLUID FORCE DENSITY VEC. FIELD)
! ==========================================================================================================
!  INTERFACES ==============================================================================================
! ==========================================================================================================
      INTERFACE DELSQP
        MODULE PROCEDURE DELSQP_SCALAR
      END INTERFACE

      INTERFACE IDELSQP
        MODULE PROCEDURE IDELSQP_SCALAR
      END INTERFACE

      INTERFACE MULXP
        MODULE PROCEDURE MULXP_SCALAR
      END INTERFACE

      INTERFACE DIVXP
        MODULE PROCEDURE DIVXP_SCALAR
      END INTERFACE

      INTERFACE MULXM
        MODULE PROCEDURE MULXM_SCALAR
      END INTERFACE

      INTERFACE DIVXM
        MODULE PROCEDURE DIVXM_SCALAR
      END INTERFACE

      INTERFACE XXDX
        MODULE PROCEDURE XXDX_SCALAR
      END INTERFACE

      INTERFACE DEL2
        MODULE PROCEDURE DEL2_SCALAR
      END INTERFACE

      INTERFACE IDEL2
        MODULE PROCEDURE IDEL2_PRELN_SCALAR
        MODULE PROCEDURE IDEL2_PROLN_SCALAR
      END INTERFACE

      INTERFACE HELM
        MODULE PROCEDURE HELM_SCALAR
      END INTERFACE

      INTERFACE IHELM
        MODULE PROCEDURE IHELM_SCALAR
      END INTERFACE

      INTERFACE HELMP
        MODULE PROCEDURE HELMP_SCALAR
      END INTERFACE

      INTERFACE IHELMP
        MODULE PROCEDURE IHELMP_SCALAR
      END INTERFACE

      INTERFACE VPROD
        MODULE PROCEDURE VPROD_VECTOR_VECTOR
      END INTERFACE

      INTERFACE PT2VEL
        MODULE PROCEDURE PT2VEL_SCALAR_SCALAR
      END INTERFACE

      INTERFACE PT2VOR
        MODULE PROCEDURE PT2VOR_SCALAR_SCALAR
      END INTERFACE

      INTERFACE PTFFRMV
        MODULE PROCEDURE PTFFRMV_SCALAR_SCALAR
      END INTERFACE

      INTERFACE VELPROJ
        MODULE PROCEDURE VELPROJ_VECTOR
      END INTERFACE

      INTERFACE UXWPROJ
        MODULE PROCEDURE UXWPROJ_VECTOR
      END INTERFACE

      INTERFACE VEL2GRADP
        MODULE PROCEDURE VEL2GRADP_VECTOR
      END INTERFACE

      INTERFACE VEL2FORCE
        MODULE PROCEDURE VEL2FORCE_VECTOR
      END INTERFACE

CONTAINS
! ==========================================================================================================
!  PUBLIC PROCEDURES =======================================================================================
! ==========================================================================================================
      FUNCTION DELSQP_SCALAR(A) RESULT(B)
! ==========================================================================================================
! [USAGE]:
! MODIFIED HORIZONTAL LAPLACIAN (1-X)^(-2)*DEL^2_PERP() OPERATOR.
! NOTE THAT A(1,1,1) (P_L^0^(R) (= 1)'S COEFF.) WILL BE LOST AFTER THE OPERATION, AS B(1,1,1) = 0 (CONST.).
! [VARAIBLES]:
! A >> SCALAR FIELD
! B >> (1-X)^(-2)*DEL^2_PERP(A)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: A
      TYPE(SCALAR_FIELD)                :: B

      INTEGER                           :: MM, NN, N

      IF (A%SPACE .NE. 'FFF') THEN  
        IF (WARN_ON) WRITE(*,*) 'DELSQP: [WARNING] INPUT SCALAR FIELD NOT IN FFF. TO FFF BEFORE OPERATION...'
        CALL TRANS(A, 'FFF')
      ENDIF

      CALL ALLOC(B, 'FFF')
      B = A

      B%LN = 0.D0

      IF (FINFO%M(1) .NE. 0) STOP 'DELSQP: THE FIRST AZIMUTHAL WAVENUMBER MUST BE 0.'

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, N)
      DO MM = 1, A%NPC
        DO NN = 1, A%NRCS(MM)
          N = FINFO%M(MM) + NN - 1
          B%E(NN,MM,:) = -A%E(NN,MM,:)*N*(N+1)/FINFO%ELL2
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      B%E(1,1,1) = 1.D0/FINFO%ELL2/TFM%NORM(1,1)*A%LN

      RETURN

      END FUNCTION
! ==========================================================================================================
      FUNCTION IDELSQP_SCALAR(A,ISTP) RESULT(B)
! ==========================================================================================================
! [USAGE]:
! INVERSE MODIFIED HORIZONTAL LAPLACIAN [(1-X)^(-2)*DEL^2_PERP()]^(-1) OPERATOR.
! NOTE THAT B(1,1,1) (P_L_0^0(R) (= 1)'S COEFF.) IS NOT RECONSTRUCTUABLE (SEE DELSQP), SO EXPECT THAT
! IDELSQP(DELSQP(A)) = A + CONSTANT, UNLESS ORIGINALLY A(1,1,1) = 0.D0.
! [VARAIBLES]:
! A >> SCALAR FIELD
! ISTP >> (OPTIONAL) IF THE FIELD IS A TOROIDAL OR POLOIDAL FIELD, PERFORM ZEROAT1().
! B >> [(1-X)^(-2)*DEL^2_PERP]^(-1)(A)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: A
      LOGICAL, OPTIONAL, INTENT(IN)     :: ISTP
      TYPE(SCALAR_FIELD)                :: B

      INTEGER                           :: MM, NN, N
      LOGICAL                           :: ISTP_ = .FALSE.

      IF (PRESENT(ISTP)) ISTP_ = ISTP

      IF (A%SPACE .NE. 'FFF') THEN  
        IF (WARN_ON) WRITE(*,*) 'DELSQP: [WARNING] INPUT SCALAR FIELD NOT IN FFF. TO FFF BEFORE OPERATION...'
        CALL TRANS(A, 'FFF')
      ENDIF

      CALL ALLOC(B, 'FFF')
      B = A

      IF (FINFO%M(1) .NE. 0) STOP 'DELSQP: THE FIRST AZIMUTHAL WAVENUMBER MUST BE 0.'

      B%LN = A%E(1,1,1)*FINFO%ELL2*TFM%NORM(1,1)

      B%E(1,1,:) = 0.D0
      DO NN = 2, A%NRCS(1)
        N = NN - 1
        B%E(NN,1,:) = -A%E(NN,1,:)/N/(N+1)*FINFO%ELL2
      ENDDO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, N)
      DO MM = 2, A%NPC
        DO NN = 1, A%NRCS(MM)
          N = FINFO%M(MM) + NN - 1
          B%E(NN,MM,:) = -A%E(NN,MM,:)/N/(N+1)*FINFO%ELL2
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      IF (ISTP_) CALL ZEROAT1(B)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MULXP_SCALAR(A) RESULT(B)
! ==========================================================================================================
! [USAGE]:
! (1+X)*, OR 2R^2/(R^2+L^2)* OPERATOR.
! [VARAIBLES]:
! A >> SCALAR FIELD
! B >> 2R^2/(R^2+L^2)*A
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: A
      TYPE(SCALAR_FIELD)                :: B

      INTEGER                           :: MM, NN
      TYPE(REAL_BNDM)                   :: XP

      IF (A%SPACE .NE. 'FFF') THEN  
        IF (WARN_ON) WRITE(*,*) 'MULXP: [WARNING] INPUT SCALAR FIELD NOT IN FFF. TO FFF BEFORE OPERATION...'
        CALL TRANS(A, 'FFF')
      ENDIF

      CALL ALLOC(B, 'FFF')
      B = A

      IF (ABS(A%LN) .GT. 1.D-14) THEN
        IF (WARN_ON) WRITE(*,*) 'MULXP: [WARNING] INPUT HAS NONZERO LOG TERM. OPERATION MAY BE INACCURATE.'
      ENDIF

      B%LN = 0.D0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, XP)
      DO MM = 1, A%NPC
        NN = A%NRCS(MM)
        XP = LEG_XP(FINFO%M(MM), NN, TFM%NORM(:, MM))
        B%E(:NN, MM, :B%NZC) = XP .MUL. A%E(:NN, MM, :A%NZC)
        B%E(:NN, MM, B%NZCU:) = XP .MUL. A%E(:NN, MM, A%NZCU:)
        DEALLOCATE( XP%E )
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION DIVXP_SCALAR(A) RESULT(B)
! ==========================================================================================================
! [USAGE]:
! 1/(1+X)*, OR (R^2+L^2)/2R^2* OPERATOR.
! [VARAIBLES]:
! A >> SCALAR FIELD
! B >> (R^2+L^2)/2R^2*A
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: A
      TYPE(SCALAR_FIELD)                :: B

      INTEGER                           :: MM, NN
      TYPE(REAL_BNDM)                   :: XP

      IF (A%SPACE .NE. 'FFF') THEN  
        IF (WARN_ON) WRITE(*,*) 'DIVXP: [WARNING] INPUT SCALAR FIELD NOT IN FFF. TO FFF BEFORE OPERATION...'
        CALL TRANS(A, 'FFF')
      ENDIF

      CALL ALLOC(B, 'FFF')
      B = A

      IF (ABS(A%LN) .GT. 1.D-14) THEN
        IF (WARN_ON) WRITE(*,*) 'DIVXP: [WARNING] INPUT HAS NONZERO LOG TERM. OPERATION MAY BE INACCURATE.'
      ENDIF

      B%LN = 0.D0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, XP)
      DO MM = 1, A%NPC
        NN = A%NRCS(MM)
        XP = LEG_XP(FINFO%M(MM), NN, TFM%NORM(:,MM))
        CALL SOLVE(XP, B%E(:NN,MM,:B%NZC))
        CALL SOLVE(XP, B%E(:NN,MM,B%NZCU:))
        DEALLOCATE( XP%E )
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MULXM_SCALAR(A) RESULT(B)
! ==========================================================================================================
! [USAGE]:
! (1-X)*, OR 2L^2/(R^2+L^2)* OPERATOR.
! [VARAIBLES]:
! A >> SCALAR FIELD
! B >> 2L^2/(R^2+L^2)*A
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: A
      TYPE(SCALAR_FIELD)                :: B

      INTEGER                           :: MM, NN
      TYPE(REAL_BNDM)                   :: XM

      IF (A%SPACE .NE. 'FFF') THEN  
        IF (WARN_ON) WRITE(*,*) 'MULXM: [WARNING] INPUT SCALAR FIELD NOT IN FFF. TO FFF BEFORE OPERATION...'
        CALL TRANS(A, 'FFF')
      ENDIF

      CALL ALLOC(B, 'FFF')
      B = A

      IF (ABS(A%LN) .GT. 1.D-14) THEN
        IF (WARN_ON) WRITE(*,*) 'MULXM: [WARNING] INPUT HAS NONZERO LOG TERM. OPERATION MAY BE INACCURATE.'
      ENDIF

      B%LN = 0.D0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, XM)
      DO MM = 1, A%NPC
        NN = A%NRCS(MM)
        XM = LEG_XM(FINFO%M(MM), NN, TFM%NORM(:, MM))
        B%E(:NN, MM, :B%NZC) = XM .MUL. A%E(:NN, MM, :A%NZC)
        B%E(:NN, MM, B%NZCU:) = XM .MUL. A%E(:NN, MM, A%NZCU:)
        DEALLOCATE( XM%E )
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION DIVXM_SCALAR(A) RESULT(B)
! ==========================================================================================================
! [USAGE]:
! 1/(1-X)*, OR (R^2+L^2)/2L^2* OPERATOR.
! [VARAIBLES]:
! A >> SCALAR FIELD
! B >> (R^2+L^2)/2L^2*A
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: A
      TYPE(SCALAR_FIELD)                :: B

      INTEGER                           :: MM, NN
      TYPE(REAL_BNDM)                   :: XM

      IF (A%SPACE .NE. 'FFF') THEN  
        IF (WARN_ON) WRITE(*,*) 'DIVXM: [WARNING] INPUT SCALAR FIELD NOT IN FFF. TO FFF BEFORE OPERATION...'
        CALL TRANS(A, 'FFF')
      ENDIF

      CALL ALLOC(B, 'FFF')
      B = A

      IF (ABS(A%LN) .GT. 1.D-14) THEN
        IF (WARN_ON) WRITE(*,*) 'DIVXM: [WARNING] INPUT HAS NONZERO LOG TERM. OPERATION MAY BE INACCURATE.'
      ENDIF

      B%LN = 0.D0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, XM)
      DO MM = 1, A%NPC
        NN = A%NRCS(MM)
        XM = LEG_XM(FINFO%M(MM), NN, TFM%NORM(:,MM))
        CALL SOLVE(XM, B%E(:NN,MM,:B%NZC))
        CALL SOLVE(XM, B%E(:NN,MM,B%NZCU:))
        DEALLOCATE( XM%E )
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION XXDX_SCALAR(A) RESULT(B)
! ==========================================================================================================
! [USAGE]:
! (1-X)^2*D()/DX, OR R*D()/DR OPERATOR.
! [VARAIBLES]:
! A >> SCALAR FIELD
! B >> R*D(A)/DR
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: A
      TYPE(SCALAR_FIELD)                :: B

      INTEGER                           :: MM, NN
      TYPE(REAL_BNDM)                   :: XXDX

      IF (A%SPACE .NE. 'FFF') THEN  
        IF (WARN_ON) WRITE(*,*) 'MULXM: [WARNING] INPUT SCALAR FIELD NOT IN FFF. TO FFF BEFORE OPERATION...'
        CALL TRANS(A, 'FFF')
      ENDIF

      CALL ALLOC(B, 'FFF')
      B = A

      B%LN = 0.D0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, XXDX)
      DO MM = 1, A%NPC
        NN = A%NRCS(MM)
        XXDX = LEG_XXDX(FINFO%M(MM), NN, TFM%NORM(:, MM))
        B%E(:NN, MM, :B%NZC) = XXDX .MUL. A%E(:NN, MM, :A%NZC)
        B%E(:NN, MM, B%NZCU:) = XXDX .MUL. A%E(:NN, MM, A%NZCU:)
        DEALLOCATE( XXDX%E )
      ENDDO
!$OMP END PARALLEL DO

      ! R*D(P_L(R))/DR = P_L_1^0(R) + P_L_0^0(R)
      B%E(1,1,1) = B%E(1,1,1) + 1.D0/TFM%NORM(1,1)*A%LN
      B%E(2,1,1) = B%E(2,1,1) + 1.D0/TFM%NORM(2,1)*A%LN

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION DEL2_SCALAR(A) RESULT(B)
! ==========================================================================================================
! [USAGE]:
! LAPLACIAN OPERATOR.
! [VARAIBLES]:
! A >> SCALAR FIELD
! B >> DEL^2(A)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: A
      TYPE(SCALAR_FIELD)                :: B

      INTEGER                           :: MM, NN, KK
      TYPE(REAL_BNDM)                   :: DEL2

      IF (A%SPACE .NE. 'FFF') THEN  
        IF (WARN_ON) WRITE(*,*) 'MULXM: [WARNING] INPUT SCALAR FIELD NOT IN FFF. TO FFF BEFORE OPERATION...'
        CALL TRANS(A, 'FFF')
      ENDIF

      CALL ALLOC(B, 'FFF')
      B = A

      B%LN = 0.D0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, KK, DEL2)
      DO MM = 1, A%NPC
        NN = A%NRCS(MM)
        DO KK = 1, FINFO%NZ
          IF ((KK .LE. A%NZC) .OR. (KK .GE. A%NZCU)) THEN
            DEL2 = LEG_DEL2(FINFO%M(MM), FINFO%AK(KK), NN, FINFO%ELL, TFM%NORM(:, MM))
            B%E(:NN, MM, KK) = DEL2 .MUL. A%E(:NN, MM, KK)
            DEALLOCATE( DEL2%E )
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! 1/R*D/DR(R*D(P_L(R))/DR) = 4/3*P_L_0^0(R) - 2*P_L_1^0(R) + 2/3*P_L_2^0(R)
      B%E(1,1,1) = B%E(1,1,1) + 4.D0/3.D0/FINFO%ELL2/TFM%NORM(1,1)*A%LN
      B%E(2,1,1) = B%E(2,1,1) - 2.D0     /FINFO%ELL2/TFM%NORM(2,1)*A%LN
      B%E(3,1,1) = B%E(3,1,1) + 2.D0/3.D0/FINFO%ELL2/TFM%NORM(3,1)*A%LN

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION IDEL2_PRELN_SCALAR(A,LN,ISTP) RESULT(B)
! ==========================================================================================================
! [USAGE]:
! INVERSE LAPLACIAN OPERATOR. AN ADDITIONAL CONDITION FOR THE INVERSION IS LOGTERM = LN, GIVEN AS INPUT.
! [VARAIBLES]:
! A >> SCALAR FIELD
! LN >> FIELD'S LOGTERM VALUE AFTER TAKING THE INVERSE
! ISTP >> (OPTIONAL) IF THE FIELD IS A TOROIDAL OR POLOIDAL FIELD, PERFORM ZEROAT1().
! B >> [DEL^(-2)](A)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT)     :: A
      LOGICAL, OPTIONAL, INTENT(IN)         :: ISTP
      REAL(P8), INTENT(IN)                  :: LN

      TYPE(SCALAR_FIELD)                    :: B

      INTEGER                               :: MM, NN, KK
      TYPE(REAL_BNDM)                       :: DEL2
      REAL(P8), DIMENSION(:,:), ALLOCATABLE :: DEL2PRE
      LOGICAL                               :: ISTP_ = .FALSE.

      IF (PRESENT(ISTP)) ISTP_ = ISTP

      IF (A%SPACE .NE. 'FFF') THEN  
        IF (WARN_ON) WRITE(*,*) 'IDEL2: [WARNING] INPUT SCALAR FIELD NOT IN FFF. TO FFF BEFORE OPERATION...'
        CALL TRANS(A, 'FFF')
      ENDIF

      CALL ALLOC(B, 'FFF')
      B = A

      IF (ABS(A%LN) .GT. 1.D-14) THEN
        IF (WARN_ON) WRITE(*,*) 'IDEL2: [WARNING] INPUT HAS NONZERO LOG TERM. OPERATION MAY BE INACCURATE.'
      ENDIF

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, KK, DEL2PRE, DEL2)
      DO MM = 1, A%NPC
        NN = A%NRCS(MM)
        DO KK = 1, FINFO%NZ
          IF ((KK .LE. A%NZC) .OR. (KK .GE. A%NZCU)) THEN
            IF ((MM .EQ. 1) .AND. (KK .EQ. 1)) THEN
              DEL2PRE = LEG_DEL2(FINFO%M(MM), FINFO%AK(KK), NN, NN, FINFO%ELL, TFM%NORM(:,MM))
              DEL2PRE(2:NN, :) = DEL2PRE(1:NN-1, :)
              DEL2PRE(1,:) = 0.D0

              DEL2PRE(1,1) = 1.D0
              DEL2PRE(2,1) = DEL2PRE(2,1) + 4.D0/3.D0/FINFO%ELL2
              DEL2PRE(3,1) = DEL2PRE(3,1) - 2.D0     /FINFO%ELL2*TFM%NORM(1,1)/TFM%NORM(2,1)
              DEL2PRE(4,1) = DEL2PRE(4,1) + 2.D0/3.D0/FINFO%ELL2*TFM%NORM(1,1)/TFM%NORM(3,1)

              DEL2 = BANDMAT(DEL2PRE, KU = 2, KL = 3)
              DEALLOCATE( DEL2PRE )

              B%E(2:NN,MM,KK) = B%E(1:NN-1,MM,KK)
              B%E(1,MM,KK) = LN/TFM%NORM(1,1)
              CALL SOLVE(DEL2, B%E(:NN, MM, KK))

              B%LN = B%E(1,MM,KK)*TFM%NORM(1,1)
            ELSE
              DEL2 = LEG_DEL2(FINFO%M(MM), FINFO%AK(KK), NN, FINFO%ELL, TFM%NORM(:,MM))
              CALL SOLVE(DEL2, B%E(:NN, MM, KK))
              DEALLOCATE( DEL2%E )
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      IF (ISTP_) CALL ZEROAT1(B)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION IDEL2_PROLN_SCALAR(A,ISTP) RESULT(B)
! ==========================================================================================================
! [USAGE]:
! INVERSE LAPLACIAN OPERATOR. AN ADDITIONAL CONDITION FOR THE INVERSION IS LN = A(1,1,1), SO LN IS COMPUTED.
! [VARAIBLES]:
! A >> SCALAR FIELD
! ISTP >> (OPTIONAL) IF THE FIELD IS A TOROIDAL OR POLOIDAL FIELD, PERFORM ZEROAT1().
! LN >> FIELD'S LOGTERM VALUE AFTER TAKING THE INVERSE
! B >> [DEL^(-2)](A)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT)     :: A
      LOGICAL, OPTIONAL, INTENT(IN)         :: ISTP
      TYPE(SCALAR_FIELD)                    :: B

      INTEGER                               :: MM, NN, KK
      TYPE(REAL_BNDM)                       :: DEL2
      REAL(P8), DIMENSION(:,:), ALLOCATABLE :: DEL2PRE
      LOGICAL                               :: ISTP_ = .FALSE.

      IF (PRESENT(ISTP)) ISTP_ = ISTP

      IF (A%SPACE .NE. 'FFF') THEN  
        IF (WARN_ON) WRITE(*,*) 'IDEL2: [WARNING] INPUT SCALAR FIELD NOT IN FFF. TO FFF BEFORE OPERATION...'
        CALL TRANS(A, 'FFF')
      ENDIF

      CALL ALLOC(B, 'FFF')
      B = A

      IF (ABS(A%LN) .GT. 1.D-14) THEN
        IF (WARN_ON) WRITE(*,*) 'IDEL2: [WARNING] INPUT HAS NONZERO LOG TERM. OPERATION MAY BE INACCURATE.'
      ENDIF

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, KK, DEL2PRE, DEL2)
      DO MM = 1, A%NPC
        NN = A%NRCS(MM)
        DO KK = 1, FINFO%NZ
          IF ((KK .LE. A%NZC) .OR. (KK .GE. A%NZCU)) THEN
            IF ((MM .EQ. 1) .AND. (KK .EQ. 1)) THEN
              DEL2PRE = LEG_DEL2(FINFO%M(MM), FINFO%AK(KK), NN, NN, FINFO%ELL, TFM%NORM(:,MM))
              DEL2PRE(1,1) = DEL2PRE(1,1) + 4.D0/3.D0/FINFO%ELL2
              DEL2PRE(2,1) = DEL2PRE(2,1) - 2.D0     /FINFO%ELL2*TFM%NORM(1,1)/TFM%NORM(2,1)
              DEL2PRE(3,1) = DEL2PRE(3,1) + 2.D0/3.D0/FINFO%ELL2*TFM%NORM(1,1)/TFM%NORM(3,1)

              DEL2 = BANDMAT(DEL2PRE, KU = 2, KL = 2)
              DEALLOCATE( DEL2PRE )

              CALL SOLVE(DEL2, B%E(:NN, MM, KK))

              B%LN = B%E(1,MM,KK)*TFM%NORM(1,1)
            ELSE
              DEL2 = LEG_DEL2(FINFO%M(MM), FINFO%AK(KK), NN, FINFO%ELL, TFM%NORM(:,MM))
              CALL SOLVE(DEL2, B%E(:NN, MM, KK))
              DEALLOCATE( DEL2%E )
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      IF (ISTP_) CALL ZEROAT1(B)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION HELM_SCALAR(A, ALP) RESULT(B)
! ==========================================================================================================
! [USAGE]:
! HELMHOLTZ OPERATOR.
! [VARIABLES]:
! A >> SCALAR FIELD
! ALP >> HELMHOLTZ CONSTANT. SEE BELOW
! B >> HELM(A) == (DEL2 + ALP*I) A
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: A
      REAL(P8), INTENT(IN)              :: ALP
      TYPE(SCALAR_FIELD)                :: B

      IF (ABS(ALP) .LT. 1.D-14) THEN
        IF (WARN_ON) WRITE(*,*) 'HELM: [WARNING] ALP=0. WHY NOT USING DEL2(A) INSTEAD?'
      ENDIF

      IF (A%SPACE .NE. 'FFF') THEN  
        IF (WARN_ON) WRITE(*,*) 'HELM: [WARNING] INPUT SCALAR FIELD NOT IN FFF. TO FFF BEFORE OPERATION...'
        CALL TRANS(A, 'FFF')
      ENDIF

      CALL ALLOC(B, 'FFF')
      B = DEL2(A)

      B%E = B%E + ALP*A%E
      B%LN = ALP*A%LN

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION IHELM_SCALAR(A, ALP, ISTP) RESULT(B)
! ==========================================================================================================
! [USAGE]:
! INVERSE HELMHOLTZ OPERATOR.
! [VARIABLES]:
! A >> SCALAR FIELD
! ALP >> HELMHOLTZ CONSTANT. SEE BELOW
! ISTP >> (OPTIONAL) IF THE FIELD IS A TOROIDAL OR POLOIDAL FIELD, PERFORM ZEROAT1().
! B >> HELM(A) == (DEL2 + ALP*I)^(-1) A
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: A
      REAL(P8), INTENT(IN)              :: ALP
      LOGICAL, OPTIONAL, INTENT(IN)     :: ISTP
      TYPE(SCALAR_FIELD)                :: B

      INTEGER                           :: MM, NN, KK
      TYPE(REAL_BNDM)                   :: HELM
      LOGICAL                           :: ISTP_ = .FALSE.

      IF (PRESENT(ISTP)) ISTP_ = ISTP

      IF (ABS(ALP) .LT. 1.D-14) STOP 'IHELM: ALPHA CANNOT BE ZERO. CONSIDER USING IDEL2(A) INSTEAD.'

      IF (A%SPACE .NE. 'FFF') THEN
        IF (WARN_ON) WRITE(*,*) 'IHELM: [WARNING] INPUT SCALAR FIELD NOT IN FFF. TO FFF BEFORE OPERATION...'
        CALL TRANS(A, 'FFF')
      ENDIF

      CALL ALLOC(B, 'FFF')
      B = A

      B%LN = A%LN/ALP

      B%E(1,1,1) = B%E(1,1,1) - 4.D0/3.D0/FINFO%ELL2/TFM%NORM(1,1)*B%LN
      B%E(2,1,1) = B%E(2,1,1) + 2.D0     /FINFO%ELL2/TFM%NORM(2,1)*B%LN
      B%E(3,1,1) = B%E(3,1,1) - 2.D0/3.D0/FINFO%ELL2/TFM%NORM(3,1)*B%LN

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, KK, HELM)
      DO MM = 1, A%NPC
        NN = A%NRCS(MM)
        DO KK = 1, FINFO%NZ
          IF ((KK .LE. A%NZC) .OR. (KK .GE. A%NZCU)) THEN
            HELM = LEG_DEL2(FINFO%M(MM), FINFO%AK(KK), NN, FINFO%ELL, TFM%NORM(:,MM))
            HELM%E(HELM%SUPLEN+1,:) = HELM%E(HELM%SUPLEN+1,:) + ALP
            CALL SOLVE(HELM, B%E(:NN,MM,KK))
            DEALLOCATE( HELM%E )
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      IF (ISTP_) CALL ZEROAT1(B)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION HELMP_SCALAR(A, P, ALP, BET) RESULT(B)
! ==========================================================================================================
! [USAGE]:
! HELMHOLTZ OPERATOR WITH THE POWER OF LAPLACIAN ADDED.
! [VARIABLES]:
! A >> SCALAR FIELD
! P >> POWER OF THE LAPLACIAN. EVEN & LARGER THAN OR EQUAL TO 4
! ALP >> HELMHOLTZ CONSTANT. SEE BELOW
! BETA >> VISCOSITY CONSTANT. SEE BELOW
! B >> HELMP(A) == (DEL^P + BET*DEL2 + ALP*I) A
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: A
      INTEGER, INTENT(IN)               :: P
      REAL(P8), INTENT(IN)              :: ALP, BET
      TYPE(SCALAR_FIELD)                :: B

      TYPE(SCALAR_FIELD)                :: D2, DP
      INTEGER                           :: I

      IF (.NOT. (MOD(P, 2) .EQ. 0) .AND. (P .GE. 4)) STOP 'HELMP: EVEN POWER >= 4 REQUIRED.'
      IF (P .GT. 8) STOP 'HELMP: CURRENTLY P UP TO 8 IS SUPPORTED ONLY.'

      IF (A%SPACE .NE. 'FFF') THEN  
        IF (WARN_ON) WRITE(*,*) 'HELMP: [WARNING] INPUT SCALAR FIELD NOT IN FFF. TO FFF BEFORE OPERATION...'
        CALL TRANS(A, 'FFF')
      ENDIF

      CALL ALLOC(B, 'FFF')
      B = DEL2(A)
      D2 = B

      DO I = 1, P/2-1
        DP = DEL2(B)
        B = DP
      ENDDO

      B%E = B%E + BET*D2%E + ALP*A%E
      B%LN = ALP*A%LN

      CALL DEALLOC( D2 )
      CALL DEALLOC( DP )

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION IHELMP_SCALAR(A, P, ALP, BET, ISTP) RESULT(B)
! ==========================================================================================================
! [USAGE]:
! INVERSE HELMHOLTZ OPERATOR WITH THE POWER OF LAPLACIAN ADDED.
! [VARIABLES]:
! A >> SCALAR FIELD
! P >> POWER OF THE LAPLACIAN. EVEN & LARGER THAN OR EQUAL TO 4
! ALP >> HELMHOLTZ CONSTANT. SEE BELOW
! BETA >> VISCOSITY CONSTANT. SEE BELOW
! ISTP >> (OPTIONAL) IF THE FIELD IS A TOROIDAL OR POLOIDAL FIELD, PERFORM ZEROAT1().
! B >> HELMP(A) == (DEL^P + BET*DEL2 + ALP*I)^(-1) A
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT)     :: A
      INTEGER, INTENT(IN)                   :: P
      REAL(P8), INTENT(IN)                  :: ALP, BET
      LOGICAL, OPTIONAL, INTENT(IN)         :: ISTP
      TYPE(SCALAR_FIELD)                    :: B

      INTEGER                               :: MM, NN, KK, I, NP
      TYPE(REAL_BNDM)                       :: HELMP, DEL2
      REAL(P8), DIMENSION(FINFO%NRCHOP)     :: BL, BL2
      REAL(P8), DIMENSION(:,:), ALLOCATABLE :: W
      LOGICAL                               :: ISTP_ = .FALSE.

      IF (PRESENT(ISTP)) ISTP_ = ISTP

      IF (.NOT. (MOD(P, 2) .EQ. 0) .AND. (P .GE. 4)) STOP 'IHELMP: EVEN POWER >= 4 REQUIRED.'
      IF (P .GT. 8) STOP 'IHELMP: CURRENTLY P UP TO 8 IS SUPPORTED ONLY.'
      IF (ABS(ALP) .LT. 1.D-14) STOP 'IHELMP: ALPHA CANNOT BE ZERO.'

      IF (A%SPACE .NE. 'FFF') THEN
        IF (WARN_ON) WRITE(*,*) 'IHELMP: [WARNING] INPUT SCALAR FIELD NOT IN FFF. TO FFF BEFORE OPERATION...'
        CALL TRANS(A, 'FFF')
      ENDIF

      CALL ALLOC(B, 'FFF')
      B = A

      B%LN = A%LN/ALP

      BL2 = 0.D0
      BL2(1) = 4.D0/3.D0/FINFO%ELL2/TFM%NORM(1,1)*B%LN
      BL2(2) = -2.D0    /FINFO%ELL2/TFM%NORM(2,1)*B%LN
      BL2(3) = 2.D0/3.D0/FINFO%ELL2/TFM%NORM(3,1)*B%LN
      BL = BL2

      DO I = 1, P/2-1
        BL = LEG_DEL2(FINFO%M(1), FINFO%AK(1), FINFO%NRCHOP, FINFO%ELL, TFM%NORM(:,1)) .MUL. BL
      ENDDO

      B%E(:FINFO%NRCHOP,1,1) = B%E(:FINFO%NRCHOP,1,1) - BL - BET*BL2

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, NP, KK, HELMP, DEL2, W)
      DO MM = 1, A%NPC
        NN = A%NRCS(MM)
        NP = NN+P
        DO KK = 1, FINFO%NZ
          IF ((KK .LE. A%NZC) .OR. (KK .GE. A%NZCU)) THEN
            HELMP = LEG_DEL2(FINFO%M(MM), FINFO%AK(KK), NP, FINFO%ELL, TFM%NORM(:,MM))
            DEL2 = HELMP
            DO I = 1, P/2-1
              HELMP = LEG_DEL2(FINFO%M(MM), FINFO%AK(KK), NP, FINFO%ELL, TFM%NORM(:,MM)) .MUL. HELMP
            ENDDO
            W = FULLMAT(HELMP) + BET*FULLMAT(DEL2)
            DO I = 1, NP
              W(I,I) = W(I,I) + ALP
            ENDDO
            HELMP = BANDMAT(W(:NN, :NN), KU=P, KL=P)
            DEALLOCATE( W )
            CALL SOLVE(HELMP, B%E(:NN,MM,KK))
            DEALLOCATE( HELMP%E )
            DEALLOCATE( DEL2%E )
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      IF (ISTP_) CALL ZEROAT1(B)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION VPROD_VECTOR_VECTOR(A,B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! CALCULATE THE CROSS PRODUCT OF TWO VECTOR FIELDS A AND B AND STORE THE RESULTS IN C.
! [VARIABLES]:
! A >> THE FIRST VECTOR FIELD (IN PPP)
! B >> THE SECOND VECTOR FIELD (IN PPP)
! C >> A X B (IN PPP SPACE)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(INOUT) :: A
      TYPE(VECTOR_FIELD), INTENT(INOUT) :: B
      TYPE(VECTOR_FIELD)                :: C

      INTEGER                           :: MM, NN, KK
      REAL(P8)                          :: A1, A2, A3, C1, C2, C3, B1, B2, B3, D1, D2, D3
      TYPE(VECTOR_FIELD)                :: AT, BT

      AT = A
      BT = B

      IF (A%SPACE .NE. 'PPP') CALL TRANS(AT, 'PPP')
      IF (B%SPACE .NE. 'PPP') CALL TRANS(BT, 'PPP')

      CALL ALLOC(C, 'PPP')

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(KK, MM, NN, A1, A2, A3, C1, C2, C3, B1, B2, B3, D1, D2, D3)
      DO KK = 1, FINFO%NZ
        DO MM = 1, FINFO%NPH
          DO NN = 1, FINFO%NR
            A1 = REAL(AT%ER(NN,MM,KK))
            A2 = REAL(AT%EP(NN,MM,KK))
            A3 = REAL(AT%EZ(NN,MM,KK))
            B1 = REAL(BT%ER(NN,MM,KK))
            B2 = REAL(BT%EP(NN,MM,KK))
            B3 = REAL(BT%EZ(NN,MM,KK))
            C1 = AIMAG(AT%ER(NN,MM,KK))
            C2 = AIMAG(AT%EP(NN,MM,KK))
            C3 = AIMAG(AT%EZ(NN,MM,KK))
            D1 = AIMAG(BT%ER(NN,MM,KK))
            D2 = AIMAG(BT%EP(NN,MM,KK))
            D3 = AIMAG(BT%EZ(NN,MM,KK))
            C%ER(NN,MM,KK) = CMPLX(A2*B3-A3*B2, C2*D3-C3*D2, P8)
            C%EP(NN,MM,KK) = CMPLX(A3*B1-A1*B3, C3*D1-C1*D3, P8)
            C%EZ(NN,MM,KK) = CMPLX(A1*B2-A2*B1, C1*D2-C2*D1, P8)
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL DEALLOC( AT )
      CALL DEALLOC( BT )

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION PT2VEL_SCALAR_SCALAR(PSI, CHI) RESULT(U)
! ==========================================================================================================
! [USAGE]:
! OBTAIN A SOLENOIDAL VELOCITY FIELD U FROM THE EQUIVALENT POLOIDAL AND TOROIDAL SCALAR FIELDS P AND T.
! [VARIABLES]:
! PSI >> SCALAR TOROIDAL FIELD
! CHI >> SCALAR POLOIDAL FIELD
! U >> RESULTING VELCOTIY FIELD (IN PPP SPACE). TO BE EQUIVALENT, SOLENOIDAL ASSUMPTION MUST BE HELD.
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: PSI
      TYPE(SCALAR_FIELD), INTENT(INOUT) :: CHI
      TYPE(VECTOR_FIELD)                :: U

      INTEGER                           :: MM, NN, KK, MV
      REAL(P8)                          :: KV
      TYPE(SCALAR_FIELD)                :: RUR, RUP, UZ, W

      IF (PSI%NRC .NE. CHI%NRC) STOP 'PT2VEL: PSI AND CHI HAVE DIFF. CHOP. INDEX IN R.'
      IF (PSI%NPC .NE. CHI%NPC) STOP 'PT2VEL: PSI AND CHI HAVE DIFF. CHOP. INDEX IN PHI.'
      IF (PSI%NZC .NE. CHI%NZC) STOP 'PT2VEL: PSI AND CHI HAVE DIFF. CHOP. INDEX IN Z.'

      RUR = XXDX_SCALAR(CHI) ! RUR = R * D(CHI)/DR
      RUP = XXDX_SCALAR(PSI) ! RUP = R * D(PSI)/DR

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, MV, KK, KV)
      DO MM = 1, PSI%NPC
        NN = PSI%NRCS(MM)
        MV = FINFO%M(MM)
        DO KK = 1, FINFO%NZ
          IF ((KK .LE. PSI%NZC) .OR. (KK .GE. PSI%NZCU)) THEN
            KV = FINFO%AK(KK)
            RUR%E(:NN,MM,KK) = IU*MV*PSI%E(:NN,MM,KK)+IU*KV*RUR%E(:NN,MM,KK) ! RUR = D/DP(PSI) + R*D/DR(D/DZ(CHI))
            RUP%E(:NN,MM,KK) = -RUP%E(:NN,MM,KK)-MV*KV*CHI%E(:NN,MM,KK) ! RUP =-R*D/DR(PSI) + D/DP(D/DZ(CHI))
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL TRANS(RUR, 'PPP')
      CALL TRANS(RUP, 'PPP')

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, KK)
      DO MM = 1, FINFO%NPH
        DO KK = 1, FINFO%NZ
          RUR%E(:FINFO%NR,MM,KK) = RUR%E(:FINFO%NR,MM,KK) / FINFO%R
          RUP%E(:FINFO%NR,MM,KK) = RUP%E(:FINFO%NR,MM,KK) / FINFO%R
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      UZ = DELSQP_SCALAR(CHI) ! UZ = (1-X)^(-2)*DEL2H(CHI)
      W = MULXM_SCALAR(UZ)    ! UZ = (1-X)^(-1)*DEL2H(CHI)
      UZ = MULXM_SCALAR(W)    ! UZ = DEL2H(CHI)
      UZ%E = -UZ%E            ! UZ = -DEL2H(CHI)

      CALL TRANS(UZ, 'PPP')

      CALL ALLOC(U, 'PPP')
      U%ER = RUR%E
      U%EP = RUP%E
      U%EZ = UZ%E

      CALL DEALLOC( UZ )
      CALL DEALLOC( W )
      CALL DEALLOC( RUP )
      CALL DEALLOC( RUR )

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION PT2VOR_SCALAR_SCALAR(PSI, CHI) RESULT(OMEGA)
! ==========================================================================================================
! [USAGE]:
! OBTAIN A SOLENOIDAL VORTICITY FIELD OMEGA FROM THE EQUIVALENT POLOIDAL AND TOROIDAL SCALAR FIELDS P AND T.
! [VARIABLES]:
! PSI >> SCALAR TOROIDAL FIELD
! CHI >> SCALAR POLOIDAL FIELD
! OMEGA >> RESULTING VORTIICTY FIELD (IN PPP SPACE). TO BE EQUIVALENT, SOLENOIDAL ASSUMPTION MUST BE HELD.
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: PSI
      TYPE(SCALAR_FIELD), INTENT(INOUT) :: CHI
      TYPE(VECTOR_FIELD)                :: OMEGA

      INTEGER                           :: MM, NN, KK, MV
      REAL(P8)                          :: KV
      TYPE(SCALAR_FIELD)                :: ROR, ROP, OZ, W

      IF (PSI%NRC .NE. CHI%NRC) STOP 'PT2VOR: PSI AND CHI HAVE DIFF. CHOP. INDEX IN R.'
      IF (PSI%NPC .NE. CHI%NPC) STOP 'PT2VOR: PSI AND CHI HAVE DIFF. CHOP. INDEX IN PHI.'
      IF (PSI%NZC .NE. CHI%NZC) STOP 'PT2VOR: PSI AND CHI HAVE DIFF. CHOP. INDEX IN Z.'

      OZ = DEL2_SCALAR(CHI) ! OZ = DEL2(CHI)
      ROR = XXDX_SCALAR(PSI) ! RUR = R * D(PSI)/DR
      ROP = XXDX_SCALAR(OZ) ! RUP = R * D(DEL2(CHI))/DR

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, MV, KK, KV)
      DO MM = 1, PSI%NPC
        NN = PSI%NRCS(MM)
        MV = FINFO%M(MM)
        DO KK = 1, FINFO%NZ
          IF ((KK .LE. PSI%NZC) .OR. (KK .GE. PSI%NZCU)) THEN
            KV = FINFO%AK(KK)
            ROR%E(:NN,MM,KK) = -IU*MV*OZ%E(:NN,MM,KK)+IU*KV*ROR%E(:NN,MM,KK) ! ROR = -D/DP(DEL^2(CHI)) + R*D/DR(D/DZ(PSI))
            ROP%E(:NN,MM,KK) = ROP%E(:NN,MM,KK)-MV*KV*PSI%E(:NN,MM,KK) ! ROP = R*D/DR(DEL^2(CHI)) + D/DP(D/DZ(PSI))
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL TRANS(ROR, 'PPP')
      CALL TRANS(ROP, 'PPP')

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, KK)
      DO MM = 1, FINFO%NPH
        DO KK = 1, FINFO%NZ
          ROR%E(:FINFO%NR,MM,KK) = ROR%E(:FINFO%NR,MM,KK) / FINFO%R
          ROP%E(:FINFO%NR,MM,KK) = ROP%E(:FINFO%NR,MM,KK) / FINFO%R
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      OZ = DELSQP_SCALAR(PSI) ! OZ = (1-X)^(-2)*DEL2H(PSI)
      W = MULXM_SCALAR(OZ)    ! OZ = (1-X)^(-1)*DEL2H(PSI)
      OZ = MULXM_SCALAR(W)    ! OZ = DEL2H(PSI)
      OZ%E = -OZ%E            ! OZ = -DEL2H(PSI)

      CALL TRANS(OZ, 'PPP')

      CALL ALLOC(OMEGA, 'PPP')
      OMEGA%ER = ROR%E
      OMEGA%EP = ROP%E
      OMEGA%EZ = OZ%E

      CALL DEALLOC( W )
      CALL DEALLOC( ROR )
      CALL DEALLOC( ROP )
      CALL DEALLOC( OZ )

      RETURN
      END FUNCTION
! ==========================================================================================================
      SUBROUTINE PTFFRMV_SCALAR_SCALAR(PSI, CHI, U)
! ==========================================================================================================
! [USAGE]:
! REMOVE THE FAR-FIELD COMPONENTS (GENERALLY FOR STABLE COMPUTATION)
! [VARIABLES]:
! PSI >> TOROIDAL COMPONENT OF U (VELOCITY FIELD) (IN FFF)
! CHI >> POLOIDAL COMPONENT OF U (VELOCITY FIELD) (IN FFF)
! U >> (OPTIONAL) VELOCITY FIELD U (IN PPP)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(SCALAR_FIELD), INTENT(INOUT) :: PSI
      TYPE(SCALAR_FIELD), INTENT(INOUT) :: CHI
      TYPE(VECTOR_FIELD), OPTIONAL      :: U

      INTEGER                           :: NS, NS0, MM, KK, I
      REAL(P8)                          :: LN
      TYPE(SCALAR_FIELD)                :: WK
  
      NS = FINFO%NR * 3/4
      NS0 = MIN(NS+4, FINFO%NR)

      WK = DELSQP_SCALAR(PSI)
      CALL TRANS(WK, 'PFP')
      CALL VERFFT_RPHYS(WK%E, WK%NZC, WK%NPC, 0)
      WK%SPACE = 'PFF' ! PFF IS NON-STANDARD YET NECESSARY FOR CALCULATION IN THIS SUBROUTINE 

      WK%E(NS0:FINFO%NR, 2:WK%NPC, 1) = 0.D0
      IF (FINFO%NZ .GT. 1) WK%E(NS0:FINFO%NR, 1, 2:FINFO%NZ) = 0.D0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, KK, I)
      DO MM = 1, WK%NPC
        DO KK = 1, FINFO%NZ
          IF ((KK .LE. WK%NZC) .OR. (KK .GE. WK%NZCU)) THEN
            WK%E(NS:FINFO%NR, MM, KK) = WK%E(NS:FINFO%NR, MM, KK)*TFM%W(NS:FINFO%NR)
            DO I = 1,5
              WK%E(NS:FINFO%NR, MM, KK) = SMOOTH(WK%E(NS:FINFO%NR, MM, KK))
            ENDDO
            WK%E(NS:FINFO%NR, MM, KK) = WK%E(NS:FINFO%NR, MM, KK)/TFM%W(NS:FINFO%NR)
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      CALL VERFFT_RPHYS(WK%E, WK%NZC, WK%NPC, 1)
      WK%SPACE = 'PFP'
      CALL TRANS(WK, 'FFF')
      PSI = IDELSQP_SCALAR(WK, ISTP=.TRUE.)

      CALL CHOP_SET(CHI, 2); CALL CHOP_SET(WK, 2)
      WK = DEL2_SCALAR(CHI)
      CALL TRANS(WK, 'PFP')
      CALL VERFFT_RPHYS(WK%E, WK%NZC, WK%NPC, 0)
      WK%SPACE = 'PFF' ! PFF IS NON-STANDARD YET NECESSARY FOR CALCULATION IN THIS SUBROUTINE 

      WK%E(NS0:FINFO%NR, 2:WK%NPC, 1) = 0.D0
      IF (FINFO%NZ .GT. 1) WK%E(NS0:FINFO%NR, 1, 2:FINFO%NZ) = 0.D0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, KK, I)
      DO MM = 1, WK%NPC
        DO KK = 1, FINFO%NZ
          IF ((KK .LE. WK%NZC) .OR. (KK .GE. WK%NZCU)) THEN
            WK%E(NS:FINFO%NR, MM, KK) = WK%E(NS:FINFO%NR, MM, KK)*TFM%W(NS:FINFO%NR) &
                                                                 /(1.D0-TFM%X(NS:FINFO%NR))**2.D0
            DO I = 1,5
              WK%E(NS:FINFO%NR, MM, KK) = SMOOTH(WK%E(NS:FINFO%NR, MM, KK))
            ENDDO
            WK%E(NS:FINFO%NR, MM, KK) = WK%E(NS:FINFO%NR, MM, KK)/TFM%W(NS:FINFO%NR) &
                                                                 *(1.D0-TFM%X(NS:FINFO%NR))**2.D0
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      CALL CHOP_RST(WK); CALL CHOP_DO(WK)
      CALL VERFFT_RPHYS(WK%E, WK%NZC, WK%NPC, 1)
      WK%SPACE = 'PFP'
      CALL TRANS(WK, 'FFF')
      LN = CHI%LN
      CHI = IDEL2_PRELN_SCALAR(WK, LN, ISTP=.TRUE.)

      IF (PRESENT(U)) U = PT2VEL(PSI, CHI)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE VELPROJ_VECTOR(U,PSI,CHI,PSILN,CHILN)
! ==========================================================================================================
! [USAGE]:
! PROJECT A (SOLENOIDAL) VECTOR FIELD (UR, UP, UZ) TO POLOIDAL-TOROIDAL SCALAR FIELDS
! [VARIABLES]:
! U >> A (SOLENOIDAL) VECTOR FIELD (UR, UP, UZ) (IN PPP SPACE)
! PSI >> TOROIDAL STREAMFUNCTION WITH RESPECT TO U (IN FFF SPACE)
! CHI >> POLOIDAL STREAMFUNCTION WITH RESPECT TO U (IN FFF SPACE). CHI%LN = 0.D0 IS ASSUMED.
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(INOUT)      :: U
      TYPE(SCALAR_FIELD), INTENT(OUT)        :: PSI
      TYPE(SCALAR_FIELD), INTENT(OUT)        :: CHI
      REAL(P8), OPTIONAL                     :: PSILN, CHILN

      TYPE(VECTOR_FIELD)                     :: UT
      INTEGER                                :: MM, NN, KK, I, N, MV
      REAL(P8)                               :: KV
      REAL(P8), DIMENSION(:), ALLOCATABLE    :: FFF
      REAL(P8), DIMENSION(:,:), ALLOCATABLE  :: V, D, T, PFD
      COMPLEX(P8), DIMENSION(:), ALLOCATABLE :: W1, W2

      IF (U%NRC .NE. FINFO%NRCHOP) STOP 'VELPROJ: U`S CHOP .NE. GLOBAL CHOP IN R.'
      IF (U%NPC .NE. FINFO%NPCHOP) STOP 'VELPROJ: U`S CHOP .NE. GLOBAL CHOP IN PHI.'
      IF (U%NZC .NE. FINFO%NZCHOP) STOP 'VELPROJ: U`S CHOP .NE. GLOBAL CHOP IN Z.'

      CALL ALLOC(PSI, 'FFF')
      CALL ALLOC(CHI, 'FFF')

      IF (.NOT. (U%SPACE .EQ. 'PPP')) STOP 'VELPROJ: INPUT MUST BE IN PPP SPACE.'
      UT = U

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, KK)
      DO MM = 1, FINFO%NPH
        DO KK = 1, FINFO%NZ
          UT%ER(:FINFO%NR,MM,KK) = UT%ER(:FINFO%NR,MM,KK) * FINFO%R ! NOW UR -> RUR
          UT%EP(:FINFO%NR,MM,KK) = UT%EP(:FINFO%NR,MM,KK) * FINFO%R ! NOW UP -> RUP
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL TRANS(UT, 'PFP')
      CALL VERFFT_RPHYS(UT%ER, UT%NZC, UT%NPC, 0)
      CALL VERFFT_RPHYS(UT%EP, UT%NZC, UT%NPC, 0)
      CALL VERFFT_RPHYS(UT%EZ, UT%NZC, UT%NPC, 0)
      UT%SPACE = 'PFF' ! PFF IS NON-STANDARD YET NECESSARY FOR CALCULATION IN THIS SUBROUTINE 

      IF (PRESENT(PSILN)) THEN
        PSI%LN = PSILN
      ELSE
        PSI%LN = -UT%EP(FINFO%NR,1,1)/(1/TFM%NORM(1,1)+1/TFM%NORM(2,1))
      ENDIF
      UT%EP(:FINFO%NR,1,1) = UT%EP(:FINFO%NR,1,1) + PSI%LN*(1.D0+TFM%X(:FINFO%NR))

      CALL CHOP_SET(UT, 3)

      ALLOCATE( FFF(FINFO%NRH) )
      ALLOCATE( V(UT%NRC, FINFO%NRH) )
      ALLOCATE( D(UT%NRC, FINFO%NRH) )
      ALLOCATE( T(UT%NRC, FINFO%NRH) )
      ALLOCATE( PFD(FINFO%NRH, UT%NRC+1) )
      ALLOCATE( W1(UT%NRC) )
      ALLOCATE( W2(UT%NRC) )

      FFF = (1.D0 - TFM%X(:FINFO%NRH)**2.D0) / TFM%W(:FINFO%NRH)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, MV, PFD, I, N, V, D, T, KK, W1, W2, KV)
      DO MM = 1, UT%NPC
        NN = UT%NRCS(MM)
        MV = FINFO%M(MM)
        PFD(:FINFO%NRH, :NN+1) = TFM%PF(:FINFO%NRH, :NN+1, MM) .MUL. LEG_XXDX(MV, NN+1, TFM%NORM(:, MM))
        DO I = 1, NN
          N = MAX(1, MV+I-1)
          V(I, :FINFO%NRH) = TFM%PF(:FINFO%NRH,I,MM)/N/(N+1)/FFF
          D(I, :FINFO%NRH) = PFD(:FINFO%NRH,I)/N/(N+1)/FFF
          T(I, :FINFO%NRH) = TFM%PF(:FINFO%NRH,I,MM)*TFM%W(:FINFO%NRH)
        ENDDO

        DO KK = 1, FINFO%NZ
          IF ((KK .LE. UT%NZC) .OR. (KK .GE. UT%NZCU)) THEN
            W1 = 0.D0; W2 = 0.D0
            IF (MV .NE. 0) THEN
              CALL EOMUL(V(:NN,:), UT%ER(:FINFO%NR,MM,KK), W1(:NN))
            ENDIF
            CALL OEMUL(D(:NN,:), UT%EP(:FINFO%NR,MM,KK), W2(:NN))
            PSI%E(:NN,MM,KK) = -IU*MV*W1(:NN) - W2(:NN)

            W1 = 0.D0; W2 = 0.D0
            CALL OEMUL(D(:NN,:), UT%ER(:FINFO%NR,MM,KK), W1(:NN))
            IF (MV .NE. 0) THEN
              CALL EOMUL(V(:NN,:), UT%EP(:FINFO%NR,MM,KK), W2(:NN))
            ENDIF
            CALL EOMUL(T(:NN,:), UT%EZ(:FINFO%NR,MM,KK), CHI%E(:NN,MM,KK))
            KV = FINFO%AK(KK)
            CHI%E(:NN,MM,KK) = IU*KV*W1(:NN)+MV*KV*W2(:NN)-CHI%E(:NN,MM,KK)
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      DEALLOCATE( FFF, V, D, T, PFD, W1, W2 )

      IF (PRESENT(CHILN)) THEN
        CHI = IDEL2_PRELN_SCALAR(CHI, LN=CHILN, ISTP=.TRUE.)
      ELSE
        CHI = IDEL2_PROLN_SCALAR(CHI, ISTP=.TRUE.)
      ENDIF

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, KK)
      DO MM = 1, FINFO%NRDIM
        DO NN = 1, FINFO%NPDIM
          DO KK = 1, FINFO%NZDIM
            PSI%E(MM,NN,KK) = CMPLX(DNINT(REAL(PSI%E(MM,NN,KK))*1.D14)*1.D-14, &
                                    DNINT(AIMAG(PSI%E(MM,NN,KK))*1.D14)*1.D-14, P8)
            CHI%E(MM,NN,KK) = CMPLX(DNINT(REAL(CHI%E(MM,NN,KK))*1.D14)*1.D-14, &
                                    DNINT(AIMAG(CHI%E(MM,NN,KK))*1.D14)*1.D-14, P8)
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL CHOP_DO(PSI); CALL ZEROAT1(PSI)
      CALL CHOP_DO(CHI); CALL ZEROAT1(CHI)

      CALL DEALLOC( UT )

      RETURN
      END SUBROUTINE 
! ==========================================================================================================
      SUBROUTINE UXWPROJ_VECTOR(U,PSIUXW,CHIUXW,PSILN,CHILN)
! ==========================================================================================================
! [USAGE]:
! PROJECT THE CROSS PRODUCT OF A VECTOR FIELD U AND CURL(U) TO POLOIDAL-TOROIDAL SCALAR FIELDS
! [VARIABLES]:
! U >> A (SOLENOIDAL) VECTOR FIELD (UR, UP, UZ) (IN PPP SPACE)
! PSI >> TOROIDAL STREAMFUNCTION WITH RESPECT TO U X CURL(U) (IN FFF SPACE)
! CHI >> POLOIDAL STREAMFUNCTION WITH RESPECT TO U X CURL(U) (IN FFF SPACE). CHI%LN = 0.D0 IS ASSUMED.
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(INOUT) :: U
      TYPE(VECTOR_FIELD)                :: W, UXW
      TYPE(SCALAR_FIELD)                :: PSI, CHI
      TYPE(SCALAR_FIELD), INTENT(OUT)   :: PSIUXW, CHIUXW
      REAL(P8), OPTIONAL                :: PSILN, CHILN

      IF (U%NRC .NE. FINFO%NRCHOP) STOP 'VELPROJ: U`S CHOP .NE. GLOBAL CHOP IN R.'
      IF (U%NPC .NE. FINFO%NPCHOP) STOP 'VELPROJ: U`S CHOP .NE. GLOBAL CHOP IN PHI.'
      IF (U%NZC .NE. FINFO%NZCHOP) STOP 'VELPROJ: U`S CHOP .NE. GLOBAL CHOP IN Z.'

      CALL ALLOC(PSI, 'FFF')
      CALL ALLOC(CHI, 'FFF')

      IF (.NOT. (U%SPACE .EQ. 'PPP')) STOP 'VELPROJ: INPUT MUST BE IN PPP SPACE.'

      IF (PRESENT(PSILN)) THEN
        IF (PRESENT(CHILN)) THEN
          CALL VELPROJ_VECTOR(U,PSI,CHI,PSILN=PSILN,CHILN=CHILN)
        ELSE
          CALL VELPROJ_VECTOR(U,PSI,CHI,PSILN=PSILN)
        ENDIF
      ELSE
        IF (PRESENT(CHILN)) THEN
          CALL VELPROJ_VECTOR(U,PSI,CHI,CHILN=CHILN)
        ELSE
          CALL VELPROJ_VECTOR(U,PSI,CHI)
        ENDIF
      ENDIF
      
      W = PT2VOR_SCALAR_SCALAR(PSI,CHI)
      UXW = VPROD_VECTOR_VECTOR(U, W)

      CALL ALLOC(PSIUXW, 'FFF')
      CALL ALLOC(CHIUXW, 'FFF')

      CALL VELPROJ_VECTOR(UXW,PSIUXW,CHIUXW,PSILN=0.D0,CHILN=0.D0)

      CALL DEALLOC( W )
      CALL DEALLOC( UXW )
      CALL DEALLOC( PSI )
      CALL DEALLOC( CHI )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      FUNCTION VEL2GRADP_VECTOR(U,PSILN,CHILN,F) RESULT(GRADP)
! ==========================================================================================================
! [USAGE]:
! COMPUTE THE GRADIENT OF PRESSURE FIELD FROM A GIVEN VELOCITY FIELD U
! [VARIABLES]:
! U >> A VECTOR FIELD (UR, UP, UZ) ASSUMED TO BE SOLENOIDAL (DIV(U) = 0) (IN PPP SPACE)
! F >> (OPTIONAL) AN EXTERNAL (NONLINEAR) FORCE FIELD (IN PPP SPACE)
! GRADP >> GRAD(P) (IN PPP SPACE)
! [NOTES]:
! 1. GRAD(P) IS A POTENTIAL VECTOR FIELD AND THEREFORE TAKING ITS T-P PROJECTION RESULTS IN ZERO. IN OTHER
! WORDS, VECLPROJ(GRAD(P)) WOULD BE ZERO. VEL2PT & VELPROJ ARE INVERSE TO EACH OTHER ONLY WHEN U HAS ONLY
! SOLENOIDAL PART (AND NO POTENTIAL PART).
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(INOUT) :: U
      TYPE(VECTOR_FIELD)                :: GRADP
      REAL(P8), OPTIONAL                :: PSILN, CHILN
      TYPE(VECTOR_FIELD), OPTIONAL      :: F

      TYPE(SCALAR_FIELD)                :: PSI, CHI
      TYPE(VECTOR_FIELD)                :: GRADUSQH, UXW, UXWSOL

      CALL ALLOC(GRADP, 'PPP')

      GRADUSQH = GRADUSQH_VECTOR(U)

      IF (PRESENT(PSILN)) THEN
        IF (PRESENT(CHILN)) THEN
          CALL VELPROJ_VECTOR(U,PSI,CHI,PSILN=PSILN,CHILN=CHILN)
        ELSE
          CALL VELPROJ_VECTOR(U,PSI,CHI,PSILN=PSILN)
        ENDIF
      ELSE
        IF (PRESENT(CHILN)) THEN
          CALL VELPROJ_VECTOR(U,PSI,CHI,CHILN=CHILN)
        ELSE
          CALL VELPROJ_VECTOR(U,PSI,CHI)
        ENDIF
      ENDIF
      UXW = PT2VOR_SCALAR_SCALAR(PSI, CHI) ! UXW NOW CONTAINS CURL(U)
      UXW = VPROD_VECTOR_VECTOR(U, UXW)    ! UXW NOW CONTAINS U X CURL(U)

      CALL VELPROJ_VECTOR(UXW,PSI,CHI,PSILN=0.D0,CHILN=0.D0)
      UXWSOL = PT2VEL_SCALAR_SCALAR(PSI, CHI) ! UXWSOL IS THE SOLENOIDAL PORTION OF U X W. (IN PPP SPACE)

      GRADP%ER =  - GRADUSQH%ER + UXW%ER - UXWSOL%ER
      GRADP%EP =  - GRADUSQH%EP + UXW%EP - UXWSOL%EP
      GRADP%EZ =  - GRADUSQH%EZ + UXW%EZ - UXWSOL%EZ

      IF (PRESENT(F)) THEN
        CALL VELPROJ_VECTOR(F,PSI,CHI,PSILN=0.D0,CHILN=0.D0)
        UXWSOL = PT2VEL_SCALAR_SCALAR(PSI, CHI) ! UXWSOL IS THE SOLENOIDAL PORTION OF F. (IN PPP SPACE)

        GRADP%ER = GRADP%ER + F%ER - UXWSOL%ER
        GRADP%EP = GRADP%EP + F%EP - UXWSOL%EP
        GRADP%EZ = GRADP%EZ + F%EZ - UXWSOL%EZ
      ENDIF

      CALL DEALLOC( PSI )
      CALL DEALLOC( CHI )
      CALL DEALLOC( GRADUSQH )
      CALL DEALLOC( UXW )
      CALL DEALLOC( UXWSOL )

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION VEL2FORCE_VECTOR(U,PSILN,CHILN,F) RESULT(FORCE)
! ==========================================================================================================
! [USAGE]:
! COMPUTE THE FORCE DENSITY FIELD FROM A GIVEN VELOCITY FIELD U
! [VARIABLES]:
! U >> A VECTOR FIELD (UR, UP, UZ) ASSUMED TO BE SOLENOIDAL (DIV(U) = 0) (IN PPP SPACE)
! F >> (OPTIONAL) AN EXTERNAL (NONLINEAR) FORCE FIELD (IN PPP SPACE)
! FORCE >> DIV(STRESS_TENSOR) = -GRAD(P) + NU*DEL2(U) (IN PPP SPACE)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(INOUT) :: U
      TYPE(VECTOR_FIELD)                :: FORCE
      REAL(P8), OPTIONAL                :: PSILN, CHILN
      TYPE(VECTOR_FIELD), OPTIONAL      :: F

      TYPE(SCALAR_FIELD)                :: PSI, CHI
      TYPE(VECTOR_FIELD)                :: GRADP, DEL2U

      CALL ALLOC(FORCE, 'PPP')

      IF (PRESENT(PSILN)) THEN
        IF (PRESENT(CHILN)) THEN
          IF (PRESENT(F)) THEN
            GRADP = VEL2GRADP_VECTOR(U, PSILN=PSILN, CHILN=CHILN, F=F)
          ELSE
            GRADP = VEL2GRADP_VECTOR(U, PSILN=PSILN, CHILN=CHILN)
          ENDIF
        ELSE
          IF (PRESENT(F)) THEN
            GRADP = VEL2GRADP_VECTOR(U, PSILN=PSILN, F=F)
          ELSE
            GRADP = VEL2GRADP_VECTOR(U, PSILN=PSILN)
          ENDIF
        ENDIF
      ELSE
        IF (PRESENT(CHILN)) THEN
          IF (PRESENT(F)) THEN
            GRADP = VEL2GRADP_VECTOR(U, CHILN=CHILN, F=F)
          ELSE
            GRADP = VEL2GRADP_VECTOR(U, CHILN=CHILN)
          ENDIF
        ELSE
          IF (PRESENT(F)) THEN
            GRADP = VEL2GRADP_VECTOR(U, F=F)
          ELSE
            GRADP = VEL2GRADP_VECTOR(U)
          ENDIF
        ENDIF
      ENDIF

      IF (PRESENT(PSILN)) THEN
        IF (PRESENT(CHILN)) THEN
          CALL VELPROJ_VECTOR(U,PSI,CHI,PSILN=PSILN,CHILN=CHILN)
        ELSE
          CALL VELPROJ_VECTOR(U,PSI,CHI,PSILN=PSILN)
        ENDIF
      ELSE
        IF (PRESENT(CHILN)) THEN
          CALL VELPROJ_VECTOR(U,PSI,CHI,CHILN=CHILN)
        ELSE
          CALL VELPROJ_VECTOR(U,PSI,CHI)
        ENDIF
      ENDIF

      PSI = DEL2_SCALAR(PSI) ! PSI NOW CONTAINS DEL2(PSI_U)
      CHI = DEL2_SCALAR(CHI) ! CHI NOW CONTAINS DEL2(CHI_U)
      DEL2U = PT2VEL_SCALAR_SCALAR(PSI, CHI) ! DEL2U=PROJ(DEL2(PSI_U),DEL2(CHI_U))=DEL2(PROJ(PSI_U,CHI_U))

      FORCE%ER = FINFO%VISC*DEL2U%ER - GRADP%ER
      FORCE%EP = FINFO%VISC*DEL2U%EP - GRADP%EP
      FORCE%EZ = FINFO%VISC*DEL2U%EZ - GRADP%EZ

      CALL DEALLOC( PSI )
      CALL DEALLOC( CHI )
      CALL DEALLOC( GRADP )
      CALL DEALLOC( DEL2U )

      RETURN
      END FUNCTION
! ==========================================================================================================
!  PRIVATE PROCEDURES ======================================================================================
! ==========================================================================================================
      SUBROUTINE EOMUL(A,B,C)
! ==========================================================================================================
! [USAGE]:
! COMPUTE C(:NI) = A(:NI, :NJ) .MUL. B(:NJ) WHERE A HAS PATTERN AS FOLLOWS:
! A(2*I  :J) =  A(2*I  :NJ-J+1) AND
! A(2*I-1:J) = -A(2*I-1:NJ-J+1)
! ONLY HALF OF A SHOULD BE GIVEN ON INPUT.
! [VARIABLES]:
! A >> A SPECIAL REAL-VALUED MATRIX WITH THE PATTERN ABOVE. DIM: NI X NJ/2
! B >> A MATRIX (VECTOR) OF THE DIMENSION NJ X 1
! C >> A .MUL. B
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:,:), INTENT(IN)   :: A
      COMPLEX(P8), DIMENSION(:), INTENT(IN)  :: B
      COMPLEX(P8), DIMENSION(:), INTENT(OUT) :: C

      COMPLEX(P8), DIMENSION(:), ALLOCATABLE :: BE, BO
      INTEGER                                :: NI, NJ, NJH

      NI = SIZE(A,1)
      NJH = SIZE(A,2)
      NJ = SIZE(B,1)

      IF (NJH*2 .NE. NJ) STOP 'EOMUL: SIZE MISMATCH.'

      ALLOCATE( BE(NJH) )
      ALLOCATE( BO(NJH) )

      BE = B(1:NJH) + B(NJ:NJH+1:-1)
      BO = B(1:NJH) - B(NJ:NJH+1:-1)

      C(1::2) = A(1::2,:) .MUL. BE
      IF (NI .GT. 1) C(2::2) = A(2::2,:) .MUL. BO

      DEALLOCATE( BE, BO )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE OEMUL(A,B,C)
! ==========================================================================================================
! [USAGE]:
! COMPUTE C(:NI) = A(:NI, :NJ) .MUL. B(:NJ) WHERE A HAS PATTERN AS FOLLOWS:
! A(2*I  :J) = -A(2*I  :NJ-J+1) AND
! A(2*I-1:J) =  A(2*I-1:NJ-J+1)
! ONLY HALF OF A SHOULD BE GIVEN ON INPUT.
! [VARIABLES]:
! A >> A SPECIAL REAL-VALUED MATRIX WITH THE PATTERN ABOVE. DIM: NI X NJ/2
! B >> A MATRIX (VECTOR) OF THE DIMENSION NJ X 1
! C >> A .MUL. B
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:,:), INTENT(IN)   :: A
      COMPLEX(P8), DIMENSION(:), INTENT(IN)  :: B
      COMPLEX(P8), DIMENSION(:), INTENT(OUT) :: C

      COMPLEX(P8), DIMENSION(:), ALLOCATABLE :: BE, BO
      INTEGER                                :: NI, NJ, NJH

      NI = SIZE(A,1)
      NJH = SIZE(A,2)
      NJ = SIZE(B,1)

      IF (NJH*2 .NE. NJ) STOP 'OEMUL: SIZE MISMATCH.'

      ALLOCATE( BE(NJH) )
      ALLOCATE( BO(NJH) )

      BE = B(1:NJH) + B(NJ:NJH+1:-1)
      BO = B(1:NJH) - B(NJ:NJH+1:-1)

      C(1::2) = A(1::2,:) .MUL. BO
      IF (NI .GT. 1) C(2::2) = A(2::2,:) .MUL. BE

      DEALLOCATE( BE, BO )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE VERFFT_RPHYS(AE, NZC, NPC, IS)
! ==========================================================================================================
! [USAGE]:
! PERFORM VERTICAL (AXIAL DIR.) FFT WITH RESPECT TO THE 3RD AXE OF A 3D ARRAY AE
! [VARIABLES]:
! AE >> 3D ARRAY CONTAINING EITHER SCALAR OR ONE COMPONENT OF A VECTOR FIELD. (DIM1, DIM2, DIM3)=(R, PHI, Z)
! NZC >> CHOP INDEX IN Z (3RD AXE)
! NPC >> CHOP INDEX IN P (2ND AXE)
! IS >> FORWARD OR BACKWARD (1: BACKWARD (F TO P), 0: FORWARD (P TO F))
! [NOTES]:
! THIS SUBROUTINE IS ALMOST SIMILAR TO VERFFT, EXCEPT THAT IT ASSUMES THE FIELD TO BE IN THE
! PHYSICAL(R)-FUNCTION(PHI)-X(Z) SPACE.
! FOR FFT, THE EXTERNAL FFTW3 LIBRARY IS UTILIZED.
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      INCLUDE 'fftw3.f03'

      COMPLEX(P8), DIMENSION(:,:,:), ALLOCATABLE           :: AE
      INTEGER, INTENT(IN)                                  :: NZC
      INTEGER, INTENT(IN)                                  :: NPC
      INTEGER, INTENT(IN)                                  :: IS

      COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), ALLOCATABLE :: B
      INTEGER                                              :: II, JJ
      TYPE(C_PTR)                                          :: PLAN

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
          DO II = 1, FINFO%NR
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
          DO II = 1, FINFO%NR
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
      FUNCTION GRADUSQH_VECTOR(U) RESULT(GRADUSQH)
! ==========================================================================================================
! [USAGE]:
! COMPUTE GRAD(U^2/2) WHERE U IS AN INPUT VECTOR FIELD
! [VARIABLES]:
! U >> A VECTOR FIELD (UR, UP, UZ) (IN PPP SPACE)
! GRADUSQH >> GRAD(U^2/2) (IN PPP SPACE)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(INOUT) :: U
      TYPE(VECTOR_FIELD)                :: GRADUSQH

      TYPE(SCALAR_FIELD)                :: USQH, W
      INTEGER                           :: NN, MM, KK, MV
      REAL(P8)                          :: R, I, KV

      CALL ALLOC(USQH, 'PPP')

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, NN, KK, R, I)
      DO NN = 1, FINFO%NR
        DO MM = 1, FINFO%NPH
          DO KK = 1, FINFO%NZ
            R = REAL(U%ER(NN,MM,KK))**2.D0+REAL(U%EP(NN,MM,KK))**2.D0+REAL(U%EZ(NN,MM,KK))**2.D0
            I = AIMAG(U%ER(NN,MM,KK))**2.D0+AIMAG(U%EP(NN,MM,KK))**2.D0+AIMAG(U%EZ(NN,MM,KK))**2.D0
            USQH%E(NN,MM,KK) = CMPLX(R, I, P8) / 2.D0
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL TRANS(USQH, 'FFF')
      CALL ALLOC(GRADUSQH, 'FFF')

      W = XXDX_SCALAR(USQH) ! R*D/DR(USQH)
      GRADUSQH%ER = W%E

      W%E = 0.D0
      DO MM = 1, USQH%NPC
        MV = FINFO%M(MM)
        W%E(:,MM,:) = USQH%E(:,MM,:)*IU*MV ! D/DP(USQH)
      ENDDO
      GRADUSQH%EP = W%E

      W%E = 0.D0
      DO KK = 1, FINFO%NZ
        IF ((KK .LE. USQH%NZC) .OR. (KK .GE. USQH%NZCU)) THEN
          KV = FINFO%AK(KK)
          W%E(:,:,KK) = USQH%E(:,:,KK)*IU*KV ! D/DZ(USQH)
        ENDIF
      ENDDO
      GRADUSQH%EZ = W%E

      CALL TRANS(GRADUSQH, 'PPP')

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(MM, KK)
      DO MM = 1, FINFO%NPH
        DO KK = 1, FINFO%NZ
          GRADUSQH%ER(:FINFO%NR,MM,KK) = GRADUSQH%ER(:FINFO%NR,MM,KK) / FINFO%R
          GRADUSQH%EP(:FINFO%NR,MM,KK) = GRADUSQH%EP(:FINFO%NR,MM,KK) / FINFO%R
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION SMOOTH(A) RESULT(AS)
! ==========================================================================================================
! [USAGE]:
! SMOOTH A 1D (COMPLEX) SEQUENCE USING A CONVOLUTION FILTERING.
! [VARIABLES]:
! A >> A SEQUENCE TO BE SMOOTHEN
! AS >> SMOOTHENED SEQUENCE CORRESPONDING TO SMOOTH(A)
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:), INTENT(IN) :: A
      COMPLEX(P8), DIMENSION(SIZE(A))       :: AS

      INTEGER                               :: I, NI
      REAL(P8), DIMENSION(3), PARAMETER     :: FAE = (/ .2D0, .3D0, .5D0 /)
      REAL(P8)                              :: F

      NI = SIZE(A)

      IF (NI .LT. 3) THEN
        STOP 'SMOOTH: THE INPUT LENGTH MUST BE LONGER THAN OR EQUAL TO 3.'
      ENDIF

      AS = 0.D0
      AS(1) = A(1)
      AS(NI-2:NI) = AS(NI-2:NI) + A(NI)*FAE

      DO I = 2, NI-1
        F = (1.D0 + (NI - I)/(NI - 1.D0)) * 0.5D0
        AS(I-1:I+1) = AS(I-1:I+1) + A(I) * (/ (1.D0-F)/2.D0, F, (1.D0-F)/2.D0 /)
      ENDDO

      RETURN
      END FUNCTION
! ==========================================================================================================
END MODULE