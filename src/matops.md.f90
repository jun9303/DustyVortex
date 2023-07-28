MODULE MATOPS ! MODULE FOR ESSENTIAL MATRIX OPERATION SUBROUTINES & FUNCTIONS
      USE OMP_LIB
      USE MISC
      IMPLICIT NONE
      PRIVATE
! ==========================================================================================================
!  PARAMS/TYPES ============================================================================================
! ==========================================================================================================
      TYPE, PUBLIC :: REAL_BNDM
        REAL(P8), DIMENSION(:,:), ALLOCATABLE :: E ! MATRIX ENTRY ELEMENTS
        INTEGER :: SUBLEN = 0 ! DEFAULT IS DIAGONAL
        INTEGER :: SUPLEN = 0 ! DEFAULT IS DIAGONAL
      END TYPE

      TYPE, PUBLIC :: CMPX_BNDM
        COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: E ! MATRIX ENTRY ELEMENTS
        INTEGER :: SUBLEN = 0 ! SUBDIAGONALS. DEFAULT IS DIAGONAL
        INTEGER :: SUPLEN = 0 ! SUPERRDIAGONALS. DEFAULT IS DIAGONAL
      END TYPE
! ==========================================================================================================
!  DECLARATIONS ============================================================================================
! ==========================================================================================================
      PUBLIC :: INV                                                            ! INVERSTION OF A MATRIX
      PUBLIC :: BANDMAT                                                        ! GENERAL MAT. TO BANDED MAT.
      PUBLIC :: BAND_TRANSPOSE                                                 ! TRANPOSE OF A BAND MATRIX
      PUBLIC :: FULLMAT                                                        ! BANDED MAT. TO GENERAL MAT.
      PUBLIC :: EIG                                                            ! EIGENVALUE PROBLEM SOLVER
      PUBLIC :: OPERATOR(.MUL.)                                                ! MATRIX MULTIPLICATION OP.
      PUBLIC :: LU                                                             ! PERFORM LU DECOMPOSITION
      PUBLIC :: SOLVE                                                         ! LINEAR PROBLEM AX=B SOLVER
! ==========================================================================================================
!  INTERFACES ==============================================================================================
! ==========================================================================================================
      INTERFACE INV
        MODULE PROCEDURE INV_C
        MODULE PROCEDURE INV_R
      END INTERFACE

      INTERFACE BANDMAT
        MODULE PROCEDURE BANDMAT_C
        MODULE PROCEDURE BANDMAT_R
      END INTERFACE

      INTERFACE BAND_TRANSPOSE
        MODULE PROCEDURE BAND_TRANSPOSE_C
        MODULE PROCEDURE BAND_TRANSPOSE_R
      END INTERFACE

      INTERFACE FULLMAT
        MODULE PROCEDURE FULLMAT_C
        MODULE PROCEDURE FULLMAT_R
      END INTERFACE

      INTERFACE EIG
        MODULE PROCEDURE GEIG_CC
        MODULE PROCEDURE GEIG_CR
        MODULE PROCEDURE GEIG_RC
        MODULE PROCEDURE GEIG_RR
        MODULE PROCEDURE EIG_C
        MODULE PROCEDURE EIG_R
      END INTERFACE

      INTERFACE OPERATOR(.MUL.)
        MODULE PROCEDURE MUL_C_C
        MODULE PROCEDURE MUL_C_C1
        MODULE PROCEDURE MUL_C1_C
        MODULE PROCEDURE MUL_C1_C1
        MODULE PROCEDURE MUL_C_CB
        MODULE PROCEDURE MUL_C1_CB
        MODULE PROCEDURE MUL_CB_C
        MODULE PROCEDURE MUL_CB_C1
        MODULE PROCEDURE MUL_CB_CB
        MODULE PROCEDURE MUL_C_R
        MODULE PROCEDURE MUL_C_R1
        MODULE PROCEDURE MUL_C1_R
        MODULE PROCEDURE MUL_C1_R1
        MODULE PROCEDURE MUL_C_RB
        MODULE PROCEDURE MUL_C1_RB
        MODULE PROCEDURE MUL_CB_R
        MODULE PROCEDURE MUL_CB_R1
        MODULE PROCEDURE MUL_CB_RB
        MODULE PROCEDURE MUL_R_C
        MODULE PROCEDURE MUL_R_C1
        MODULE PROCEDURE MUL_R1_C
        MODULE PROCEDURE MUL_R1_C1
        MODULE PROCEDURE MUL_R_CB
        MODULE PROCEDURE MUL_R1_CB
        MODULE PROCEDURE MUL_RB_C
        MODULE PROCEDURE MUL_RB_C1
        MODULE PROCEDURE MUL_RB_CB
        MODULE PROCEDURE MUL_R_R
        MODULE PROCEDURE MUL_R_R1
        MODULE PROCEDURE MUL_R1_R
        MODULE PROCEDURE MUL_R1_R1
        MODULE PROCEDURE MUL_R_RB
        MODULE PROCEDURE MUL_R1_RB
        MODULE PROCEDURE MUL_RB_R
        MODULE PROCEDURE MUL_RB_R1
        MODULE PROCEDURE MUL_RB_RB
      END INTERFACE

      INTERFACE LU
        MODULE PROCEDURE LU_C
        MODULE PROCEDURE LU_CB
        MODULE PROCEDURE LU_R
        MODULE PROCEDURE LU_RB
      END INTERFACE

      INTERFACE SOLVE
        MODULE PROCEDURE SOLVE_C_C
        MODULE PROCEDURE SOLVE_C_C_MULTI
        MODULE PROCEDURE SOLVE_CB_C
        MODULE PROCEDURE SOLVE_CB_C_MULTI
        MODULE PROCEDURE SOLVE_R_C
        MODULE PROCEDURE SOLVE_R_C_MULTI
        MODULE PROCEDURE SOLVE_RB_C
        MODULE PROCEDURE SOLVE_RB_C_MULTI
        MODULE PROCEDURE SOLVE_R_R
        MODULE PROCEDURE SOLVE_R_R_MULTI
        MODULE PROCEDURE SOLVE_RB_R
        MODULE PROCEDURE SOLVE_RB_R_MULTI
      END INTERFACE
CONTAINS
! ==========================================================================================================
!  PUBLIC PROCEDURES =======================================================================================
! ==========================================================================================================
      FUNCTION INV_R(A) RESULT(AINV)
! ==========================================================================================================
! [USAGE]: 
! INVERSE OF A GIVEN MATRIX A. GET A^(-1) USING LAPACK
! [VARIABLES]:
! A >> NI X NJ REAL MATRIX
! AINV >> INVERSE OF A (A^(-1)) CORRESPONDING TO INV_R(A)
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:,:)                    :: A
      REAL(P8), DIMENSION(SIZE(A, 1), SIZE(A, 2)) :: AINV
      REAL(P8), DIMENSION(SIZE(A, 1))             :: WORK ! WORK ARRAY FOR LAPACK
      INTEGER, DIMENSION(SIZE(A, 1))              :: IPIV ! PIVOT INDICES
      INTEGER                                     :: N, INFO

      ! LAPACK LIBRARY
      EXTERNAL DGETRF
      EXTERNAL DGETRI

      AINV = A
      N = SIZE(A, 1)

      IF (SIZE(A, 2) .NE. N) THEN
        WRITE(*,*) 'INV: INPUT MATRIX NOT SQUARE.'
        STOP
      ENDIF

      ! DGETRF COMPUTES AN LU FACTORIZATION OF A GENERAL NI X NJ MATRIX A
      CALL DGETRF(N, N, AINV, N, IPIV, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'INV: MATRIX IS NUMERICALLY SINGULAR.'
        STOP
      ENDIF

      ! DGETRI COMPUTES THE INVERSE OF A MATRIX USING THE LU FACTORIZATION COMPUTED BY DGETRF.
      CALL DGETRI(N, AINV, N, IPIV, WORK, N, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'INV: MATRIX INVERSION FAILED.'
        STOP
      ENDIF

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION INV_C(A) RESULT(AINV)
! ==========================================================================================================
! [USAGE]: 
! INVERSE OF A GIVEN MATRIX A. GET A^(-1) USING LAPACK
! [VARIABLES]:
! A >> NI X NJ COMPLEX MATRIX
! AINV >> INVERSE OF A (A^(-1)) CORRESPONDING TO INV_R(A)
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:,:)                    :: A
      COMPLEX(P8), DIMENSION(SIZE(A, 1), SIZE(A, 2)) :: AINV
      COMPLEX(P8), DIMENSION(SIZE(A, 1))             :: WORK ! WORK ARRAY FOR LAPACK
      INTEGER, DIMENSION(SIZE(A, 1))                 :: IPIV ! PIVOT INDICES
      INTEGER                                        :: N, INFO

      ! LAPACK LIBRARY
      EXTERNAL ZGETRF
      EXTERNAL ZGETRI

      AINV = A
      N = SIZE(A, 1)

      IF (SIZE(A, 2) .NE. N) THEN
        WRITE(*,*) 'INV: INPUT MATRIX NOT SQUARE.'
        STOP
      ENDIF

      ! ZGETRF COMPUTES AN LU FACTORIZATION OF A GENERAL NI X NJ MATRIX A
      CALL ZGETRF(N, N, AINV, N, IPIV, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'INV: MATRIX IS NUMERICALLY SINGULAR.'
        STOP
      ENDIF

      ! ZGETRI COMPUTES THE INVERSE OF A MATRIX USING THE LU FACTORIZATION COMPUTED BY DGETRF.
      CALL ZGETRI(N, AINV, N, IPIV, WORK, N, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'INV: MATRIX INVERSION FAILED.'
        STOP
      ENDIF

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION BANDMAT_C(A, KL, KU) RESULT(ABND)
! ==========================================================================================================
! [USAGE]: 
! GET A BAND MATRIX ABND FROM A BY TRUNCATING ITS KL SUBDIAGONALS AND KU SUPERDIAGONALS
! [VARIABLES]:
! A >> NI X NJ COMPLEX MATRIX
! KL >> NUMBER OF SUBDIAGONALS TO REMAIN. ALWAYS LARGER THAN OR EQUAL TO 0
! KU >> NUMBER OF SUPERDIAGONALS TO REMAIN. ALWAYS LARGER THAN OR EQUAL TO 0
! ABND >> RESULTING BAND MATRIX CORRESPONDING TO BANDMAT_C(A, KL, KU)
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER                     :: KL, KU
      COMPLEX(P8), DIMENSION(:,:) :: A
      TYPE(CMPX_BNDM)             :: ABND

      INTEGER                     :: NI, NJ
      INTEGER                     :: I, J

      IF (KL .LT. 0) KL = 0
      IF (KU .LT. 0) KU = 0

      NI = SIZE(A, 1)
      NJ = SIZE(A, 2)

      IF (NI .LT. KL) THEN
        WRITE(*,*) 'BANDMAT: INCOMPATIBLE SUBDIAGONAL NUMBER.'
        STOP
      ENDIF

      IF (NJ .LT. KU) THEN
        WRITE(*,*) 'BANDMAT: INCOMPATIBLE SUBDIAGONAL NUMBER.'
        STOP
      ENDIF

      ABND%SUBLEN = KL
      ABND%SUPLEN = KU
      ALLOCATE( ABND%E(KL+KU+1, NJ) )
      ABND%E = HUGE(1.D0) ! INITIALIZATION AS MACHINE INFINITY IN ORDER TO RECOGNIZE UNUSED ENTRIES

      DO J = 1, NJ
        DO I = MAX(1, J-KU), MIN(NI, J+KL)
          ABND%E(KU+1+I-J, J) = A(I, J)
        ENDDO
      ENDDO

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION BANDMAT_R(A, KL, KU) RESULT(ABND)
! ==========================================================================================================
! [USAGE]: 
! GET A BAND MATRIX ABND FROM A BY TRUNCATING ITS KL SUBDIAGONALS AND KU SUPERDIAGONALS
! [VARIABLES]:
! A >> NI X NJ REAL MATRIX
! KL >> NUMBER OF SUBDIAGONALS TO REMAIN. ALWAYS LARGER THAN OR EQUAL TO 0
! KU >> NUMBER OF SUPERDIAGONALS TO REMAIN. ALWAYS LARGER THAN OR EQUAL TO 0
! ABND >> RESULTING BAND MATRIX CORRESPONDING TO BANDMAT_R(A, KL, KU)
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER                  :: KL, KU
      REAL(P8), DIMENSION(:,:) :: A
      TYPE(REAL_BNDM)          :: ABND

      INTEGER                  :: NI, NJ
      INTEGER                  :: I, J

      IF (KL .LT. 0) KL = 0
      IF (KU .LT. 0) KU = 0

      NI = SIZE(A, 1)
      NJ = SIZE(A, 2)

      IF (NI .LT. KL) THEN
        WRITE(*,*) 'BANDMAT: INCOMPATIBLE SUBDIAGONAL NUMBER.'
        STOP
      ENDIF

      IF (NJ .LT. KU) THEN
        WRITE(*,*) 'BANDMAT: INCOMPATIBLE SUBDIAGONAL NUMBER.'
        STOP
      ENDIF

      ABND%SUBLEN = KL
      ABND%SUPLEN = KU
      ALLOCATE( ABND%E(KL+KU+1, NJ) )
      ABND%E = HUGE(1.D0) ! INITIALIZATION AS MACHINE INFINITY IN ORDER TO RECOGNIZE UNUSED ENTRIES

      DO J = 1, NJ
        DO I = MAX(1, J-KU), MIN(NI, J+KL)
          ABND%E(KU+1+I-J, J) = A(I, J)
        ENDDO
      ENDDO

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION BAND_TRANSPOSE_C(B) RESULT (BT)
! ==========================================================================================================
! [USAGE]:                                   
! CREATE A TRANPOSED COMPLEX BAND MATRIX
! [VARIABLES]:                                  
! B >> AN ARBITRARY COMPLEX BAND MATRIX
! BT >> B'S TRANSPOSE IN A BAND FORM, CORRESPONDING TO BAND_TRANSPOSE_C(B)
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(CMPX_BNDM) :: B
      TYPE(CMPX_BNDM) :: BT

      INTEGER         :: KU, KL
      INTEGER         :: NI, NJ
      INTEGER         :: I, J

      NI = SIZE(B%E, 1)
      NJ = SIZE(B%E, 2)

      KU = B%SUPLEN
      KL = B%SUBLEN

      I = COUNT(B%E(NI,:) .EQ. HUGE(1.D0))
      J = COUNT(B%E(:,1) .EQ. HUGE(1.D0))

      ALLOCATE( BT%E(NI, (NI-I)+(NJ-J)-1) )
      BT%E = HUGE(1.D0)

      BT%SUPLEN = KL
      BT%SUBLEN = KU

      DO I = 1, NI
        J = MIN(SIZE(B%E, 2), SIZE(BT%E, 2))
        BT%E(I, 1:J) = B%E(NI+1-I,:)
        BT%E(I,:) = CSHIFT(BT%E(I,:), -NI+I+KU )
      ENDDO

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION BAND_TRANSPOSE_R(B) RESULT (BT)
! ==========================================================================================================
! [USAGE]:                                   
! CREATE A TRANPOSED REAL BAND MATRIX
! [VARIABLES]:                                  
! B >> AN ARBITRARY REAL BAND MATRIX
! BT >> B'S TRANSPOSE IN A BAND FORM, CORRESPONDING TO BAND_TRANSPOSE_R(B)
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(REAL_BNDM) :: B
      TYPE(REAL_BNDM) :: BT

      INTEGER         :: KU, KL
      INTEGER         :: NI, NJ
      INTEGER         :: I, J

      NI = SIZE(B%E, 1)
      NJ = SIZE(B%E, 2)

      KU = B%SUPLEN
      KL = B%SUBLEN

      I = COUNT(B%E(NI,:) .EQ. HUGE(1.D0))
      J = COUNT(B%E(:,1) .EQ. HUGE(1.D0))

      ALLOCATE( BT%E(NI, (NI-I)+(NJ-J)-1) )
      BT%E = HUGE(1.D0)

      BT%SUPLEN = KL
      BT%SUBLEN = KU

      DO I = 1, NI
        J = MIN(SIZE(B%E, 2), SIZE(BT%E, 2))
        BT%E(I, 1:J) = B%E(NI+1-I,:)
        BT%E(I,:) = CSHIFT(BT%E(I,:), -NI+I+KU )
      ENDDO


      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION FULLMAT_C(ABND) RESULT(A)
! ==========================================================================================================
! [USAGE]:
! GET A FULL (GENERAL) MATRIX CONSTRUCTED FROM A GIVEN BAND MATRIX ABND
! [VARIABLES]:
! ABND >> BAND MATRIX (COMPLEX)
! A >> RESULTING FULL MATRIX CORRESPONDING TO FULLMAT_C(ABND) 
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(CMPX_BNDM)                          :: ABND
      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: A

      INTEGER                               :: KU, KL
      INTEGER                               :: NI, NJ
      INTEGER                               :: I, J

      NI = SIZE(ABND%E, 1)
      NJ = SIZE(ABND%E, 2)

      KU = ABND%SUPLEN
      KL = ABND%SUBLEN

      I = COUNT(ABND%E(NI,:) .EQ. HUGE(1.D0))
      J = COUNT(ABND%E(:,1) .EQ. HUGE(1.D0))

      ALLOCATE(A((NI-I)+(NJ-J)-1,NJ))
      A = 0.D0

      NI = SIZE(A, 1)
      NJ = SIZE(A, 2)

      DO J = 1, NJ
        DO I = MAX(1, J-KU), MIN(NI, J+KL)
          A(I, J) = ABND%E(KU+1+I-J, J)
        ENDDO
      ENDDO

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION FULLMAT_R(ABND) RESULT(A)
! ==========================================================================================================
! [USAGE]:
! GET A FULL (GENERAL) MATRIX CONSTRUCTED FROM A GIVEN BAND MATRIX ABND
! [VARIABLES]:
! ABND >> BAND MATRIX (REAL)
! A >> RESULTING FULL MATRIX CORRESPONDING TO FULLMAT_R(ABND) 
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(REAL_BNDM)                       :: ABND
      REAL(P8), DIMENSION(:,:), ALLOCATABLE :: A

      INTEGER                               :: KU, KL
      INTEGER                               :: NI, NJ
      INTEGER                               :: I, J

      NI = SIZE(ABND%E, 1)
      NJ = SIZE(ABND%E, 2)

      KU = ABND%SUPLEN
      KL = ABND%SUBLEN

      I = COUNT(ABND%E(NI,:) .EQ. HUGE(1.D0))
      J = COUNT(ABND%E(:,1) .EQ. HUGE(1.D0))

      ALLOCATE(A((NI-I)+(NJ-J)-1,NJ))
      A = 0.D0

      NI = SIZE(A, 1)
      NJ = SIZE(A, 2)

      DO J = 1, NJ
        DO I = MAX(1, J-KU), MIN(NI, J+KL)
          A(I, J) = ABND%E(KU+1+I-J, J)
        ENDDO
      ENDDO

      RETURN
      END FUNCTION
! ==========================================================================================================
      SUBROUTINE GEIG_CC(A, B, EV, ER, EL)
! ==========================================================================================================
! [USAGE]:                                   
! SOLVE THE EIGENVALUE PROBLEM AV = BVD FOR GIVEN SQUARE MATRICES A, B IN ADDITION TO THE CALCULATION OF
! THE RIGHT EIGENVECTOR MATRIX R. THE LEFT EIGENVECTOR MATRIX CAN BE COMPUTED OPTINALLY
! [VARIABLES]:                
! A >> N X N COMPLEX NON-HERMITIAN (GENERAL) MATRIX
! B >> N X N COMPLEX NON-HERMITIAN (GENERAL) MATRIX
! EV >> N COMPLEX EIGENVALUES IN N-VECTOR
! ER >> RIGHT EIGENVECTOR MATRIX SATISFYING AR = R X DIAG(EV)
! EL >> (OPTIONAL) LEFT EIGENVECTOR MATRIX SATISFYING L*A = DIAG(EV)X L* 
!                  NORMALIZATION IS PERFORMED TO SATISFY EL* X ER = I_N 
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:,:), INTENT(IN)              :: A, B
      COMPLEX(P8), DIMENSION(:), INTENT(INOUT)             :: EV
      COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT)           :: ER
      COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: EL

      COMPLEX(P8), DIMENSION(SIZE(A,1),SIZE(A,2))          :: W1, W2

      COMPLEX(P8), DIMENSION(SIZE(A,1))                    :: AEIG, BEIG
      COMPLEX(P8), DIMENSION(2*SIZE(A,1))                  :: WK
      REAL(P8)   , DIMENSION(8*SIZE(A,1))                  :: RWK

      INTEGER                                              :: I, NI, LW, INFO
      REAL(P8), DIMENSION(SIZE(EV))                        :: IMEV

      ! LAPACK LIBRARY
      EXTERNAL ZGGEV

      IF (SIZE(A,1) .NE. SIZE(A,2)) THEN
        WRITE (*,*) 'EIG: FIRST MATRIX NOT SQUARE.'
        WRITE (*,*) SIZE(A,1),'X',SIZE(A,2)
        STOP
      ENDIF

      IF (SIZE(B,1) .NE. SIZE(B,2)) THEN
        WRITE (*,*) 'EIG: SECOND MATRIX NOT SQUARE.'
        WRITE (*,*) SIZE(B,1),'X',SIZE(B,2)
        STOP
      ENDIF

      IF (SIZE(A,1) .NE. SIZE(B,1)) THEN
        WRITE (*,*) 'EIG: INPUT MATRICES INCOMPATIBLE IN SIZE.'
        STOP
      ENDIF

      IF (SIZE(A,1) .GT. SIZE(EV)) THEN
        WRITE (*,*) 'EIG: EV VECTOR SIZE TOO SMALL.'
        WRITE (*,*) SIZE(A,1),', BUT SIZE(EV) = ',SIZE(EV)
        STOP
      ENDIF

      NI = SIZE(A, 1)

      W1 = A
      W2 = B ! PRIVATE FUNCTION

      LW = SIZE(WK)
      ER = 0.D0*IU

      IF (PRESENT(EL)) EL = 0.D0*IU

      IF (PRESENT(EL)) THEN
        CALL ZGGEV('V', 'V', NI, W1, NI, W2, NI, AEIG, BEIG, EL, NI, ER, NI, WK, LW, RWK, INFO)
        ! NORMALIZE LEFT EIGENVECTORS IN ORDER TO KEEP L*R = ID
        DO I = 1, NI
          EL(:,I) = EL(:,I) / CONJG(DOT_PRODUCT(EL(:,I),ER(:,I)))
        ENDDO
      ELSE
        CALL ZGGEV('N', 'V', NI, W1, NI, W2, NI, AEIG, BEIG, EL, NI, ER, NI, WK, LW, RWK, INFO)
      ENDIF

      IF (INFO .NE. 0) THEN
        WRITE (*,*) 'EIG: EVP SOLVER FAILED.'
        STOP
      ENDIF

      DO I = 1, NI
        IF(ABS(BEIG(I)).EQ.0.D0) THEN
          ! ASSIGN MACHINE INFINITY DIRECTLY TO AVOID THE DIVISION-BY-ZERO ERROR FROM AEIG/BEIG
          EV(I) = HUGE(1.D0)
        ELSE
          EV(I) = AEIG(I) / BEIG(I)
        ENDIF
      ENDDO

      IMEV = AIMAG(EV)
      CALL QUICKSORT_NORECUR(IMEV, EV, ER, EL) ! PRIVATE SUBROUTINE

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE GEIG_CR(A, B, EV, ER, EL)
! ==========================================================================================================
! [USAGE]:                                   
! SOLVE THE EIGENVALUE PROBLEM AV = BVD FOR GIVEN SQUARE MATRICES A, B IN ADDITION TO THE CALCULATION OF
! THE RIGHT EIGENVECTOR MATRIX R. THE LEFT EIGENVECTOR MATRIX CAN BE COMPUTED OPTINALLY
! [VARIABLES]:                
! A >> N X N COMPLEX NON-HERMITIAN (GENERAL) MATRIX
! B >> N X N REAL NON-HERMITIAN (GENERAL) MATRIX
! EV >> N COMPLEX EIGENVALUES IN N-VECTOR
! ER >> RIGHT EIGENVECTOR MATRIX SATISFYING AR = R X DIAG(EV)
! EL >> (OPTIONAL) LEFT EIGENVECTOR MATRIX SATISFYING L*A = DIAG(EV)X L* 
!                  NORMALIZATION IS PERFORMED TO SATISFY EL* X ER = I_N 
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:,:), INTENT(IN)              :: A
      REAL(P8), DIMENSION(:,:), INTENT(IN)                 :: B
      COMPLEX(P8), DIMENSION(:), INTENT(INOUT)             :: EV
      COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT)           :: ER
      COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: EL

      CALL GEIG_CC(A, CMPLX(B), EV, ER, EL)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE GEIG_RC(A, B, EV, ER, EL)
! ==========================================================================================================
! [USAGE]:                                   
! SOLVE THE EIGENVALUE PROBLEM AV = BVD FOR GIVEN SQUARE MATRICES A, B IN ADDITION TO THE CALCULATION OF
! THE RIGHT EIGENVECTOR MATRIX R. THE LEFT EIGENVECTOR MATRIX CAN BE COMPUTED OPTINALLY
! [VARIABLES]:                
! A >> N X N REAL NON-HERMITIAN (GENERAL) MATRIX
! B >> N X N COMPLEX NON-HERMITIAN (GENERAL) MATRIX
! EV >> N COMPLEX EIGENVALUES IN N-VECTOR
! ER >> RIGHT EIGENVECTOR MATRIX SATISFYING AR = R X DIAG(EV)
! EL >> (OPTIONAL) LEFT EIGENVECTOR MATRIX SATISFYING L*A = DIAG(EV)X L* 
!                  NORMALIZATION IS PERFORMED TO SATISFY EL* X ER = I_N 
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:,:), INTENT(IN)                 :: A
      COMPLEX(P8), DIMENSION(:,:), INTENT(IN)              :: B
      COMPLEX(P8), DIMENSION(:), INTENT(INOUT)             :: EV
      COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT)           :: ER
      COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: EL

      CALL GEIG_CC(CMPLX(A), B, EV, ER, EL)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE GEIG_RR(A, B, EV, ER, EL)
! ==========================================================================================================
! [USAGE]:                                   
! SOLVE THE EIGENVALUE PROBLEM AV = BVD FOR GIVEN SQUARE MATRICES A, B IN ADDITION TO THE CALCULATION OF
! THE RIGHT EIGENVECTOR MATRIX R. THE LEFT EIGENVECTOR MATRIX CAN BE COMPUTED OPTINALLY
! [VARIABLES]:                
! A >> N X N REAL NON-HERMITIAN (GENERAL) MATRIX
! B >> N X N REAL NON-HERMITIAN (GENERAL) MATRIX
! EV >> N COMPLEX EIGENVALUES IN N-VECTOR
! ER >> RIGHT EIGENVECTOR MATRIX SATISFYING AR = R X DIAG(EV)
! EL >> (OPTIONAL) LEFT EIGENVECTOR MATRIX SATISFYING L*A = DIAG(EV)X L* 
!                  NORMALIZATION IS PERFORMED TO SATISFY EL* X ER = I_N 
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:,:), INTENT(IN)                 :: A
      REAL(P8), DIMENSION(:,:), INTENT(IN)                 :: B
      COMPLEX(P8), DIMENSION(:), INTENT(INOUT)             :: EV
      COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT)           :: ER
      COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: EL

      CALL GEIG_CC(CMPLX(A), CMPLX(B), EV, ER, EL)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE EIG_C(A, EV, ER, EL)
! ==========================================================================================================
! [USAGE]:                                   
! SOLVE THE EIGENVALUE PROBLEM AV = VD FOR GIVEN SQUARE MATRICES A IN ADDITION TO THE CALCULATION OF
! THE RIGHT EIGENVECTOR MATRIX R. THE LEFT EIGENVECTOR MATRIX CAN BE COMPUTED OPTINALLY
! [VARIABLES]:                
! A >> N X N COMPLEX NON-HERMITIAN (GENERAL) MATRIX
! EV >> N COMPLEX EIGENVALUES IN N-VECTOR
! ER >> RIGHT EIGENVECTOR MATRIX SATISFYING AR = R X DIAG(EV)
! EL >> (OPTIONAL) LEFT EIGENVECTOR MATRIX SATISFYING L*A = DIAG(EV)X L* 
!                  NORMALIZATION IS PERFORMED TO SATISFY EL* X ER = I_N 
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:,:), INTENT(IN)              :: A
      COMPLEX(P8), DIMENSION(:), INTENT(INOUT)             :: EV
      COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT)           :: ER
      COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: EL

      CALL GEIG_CR(A, EYE(SIZE(A, 1)), EV, ER, EL) ! EYE: CREATES IDENTITY MAT. PRIVATE PROCEDURE

      RETURN
      END SUBROUTINE
! ==========================================================================================================      
      SUBROUTINE EIG_R(A, EV, ER, EL)
! ==========================================================================================================
! [USAGE]:                                   
! SOLVE THE EIGENVALUE PROBLEM AV = VD FOR GIVEN SQUARE MATRICES A IN ADDITION TO THE CALCULATION OF
! THE RIGHT EIGENVECTOR MATRIX R. THE LEFT EIGENVECTOR MATRIX CAN BE COMPUTED OPTINALLY
! [VARIABLES]:                
! A >> N X N REAL NON-HERMITIAN (GENERAL) MATRIX
! EV >> N COMPLEX EIGENVALUES IN N-VECTOR
! ER >> RIGHT EIGENVECTOR MATRIX SATISFYING AR = R X DIAG(EV)
! EL >> (OPTIONAL) LEFT EIGENVECTOR MATRIX SATISFYING L*A = DIAG(EV)X L* 
!                  NORMALIZATION IS PERFORMED TO SATISFY EL* X ER = I_N 
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:,:), INTENT(IN)                 :: A
      COMPLEX(P8), DIMENSION(:), INTENT(INOUT)             :: EV
      COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT)           :: ER
      COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: EL

      CALL GEIG_RR(A, EYE(SIZE(A, 1)), EV, ER, EL) ! EYE: CREATES IDENTITY MAT. PRIVATE PROCEDURE

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      FUNCTION MUL_C_C(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X N COMPLEX MATRIX
! B >> N X K COMPLEX MATRIX
! C >> M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:,:), INTENT(IN)     :: A
      COMPLEX(P8), DIMENSION(:,:), INTENT(IN)     :: B
      COMPLEX(P8), DIMENSION(SIZE(A,1),SIZE(B,2)) :: C

      IF (SIZE(A,2) .NE. SIZE(B,1)) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ENDIF

      C = MATMUL(A, B)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_C_C1(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X N COMPLEX MATRIX
! B >> N X 1 COMPLEX 1-D MATRIX (VECTOR)
! C >> M X 1 COMPLEX 1-D MATRIX (VECTOR) CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:,:), INTENT(IN) :: A
      COMPLEX(P8), DIMENSION(:), INTENT(IN)   :: B
      COMPLEX(P8), DIMENSION(SIZE(A,1))       :: C

      IF (SIZE(A,2) .NE. SIZE(B,1)) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ENDIF

      C = RESHAPE(MATMUL(A, RESHAPE(B, (/SIZE(B), 1/))), (/SIZE(A, 1)/))

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_C1_C(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X 1 COMPLEX 1-D MATRIX (VECTOR)
! B >> N X K COMPLEX MATRIX
! C >> M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:), INTENT(IN)        :: A
      COMPLEX(P8), DIMENSION(:,:), INTENT(IN)      :: B
      COMPLEX(P8), DIMENSION(SIZE(A,1), SIZE(B,2)) :: C

      IF (1 .NE. SIZE(B,1)) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ENDIF

      C = MATMUL(RESHAPE(A, (/SIZE(A), 1/)), B)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_C1_C1(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X 1 COMPLEX 1-D MATRIX (VECTOR)
! B >> N X 1 COMPLEX 1-D MATRIX (VECTOR)
! C >> M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:), INTENT(IN)   :: A
      COMPLEX(P8), DIMENSION(:), INTENT(IN)   :: B
      COMPLEX(P8), DIMENSION(SIZE(A,1))       :: C

      IF (1 .NE. SIZE(B,1)) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ENDIF

      C = RESHAPE(MATMUL(RESHAPE(A, (/SIZE(A), 1/)), RESHAPE(B, (/SIZE(B), 1/))), (/SIZE(A, 1)/))

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_C_CB(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X N COMPLEX MATRIX
! B >> BAND MATRIX THAT IS ORIGINALLY N X K COMPLEX MATRIX IN A FULL FORM
! C >> M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:,:), INTENT(IN)  :: A
      TYPE(CMPX_BNDM), INTENT(IN)              :: B
      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: C

      INTEGER                                  :: N, NI, NJ, I, J
      INTEGER                                  :: A0, A1, B0, B1

      N = SIZE(A, 2)

      NI = SIZE(B%E, 1)
      NJ = SIZE(B%E, 2)

      I = COUNT(B%E(NI,:) .EQ. HUGE(1.D0))
      J = COUNT(B%E(:,1) .EQ. HUGE(1.D0))

      IF (N .LT. (NI-I)+(NJ-J)-1) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ELSE
        N = (NI-I)+(NJ-J)-1
      ENDIF

      ALLOCATE( C(SIZE(A, 1), SIZE(B%E, 2)) )
      C = 0.D0

      DO I = 1, SIZE(C, 1)
        DO J = 1, SIZE(C, 2)
          A0 = MAX(1, J-B%SUPLEN)
          A1 = MIN(N, J+B%SUBLEN)
          B0 = MAX(  2+B%SUPLEN-J, 1)
          B1 = MIN(N+1+B%SUPLEN-J, N)
          C(I,J) = SUM(A(I,A0:A1)*B%E(B0:B1,J))
        ENDDO
      ENDDO

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_C1_CB(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X 1 COMPLEX MATRIX (VECTOR)
! B >> BAND MATRIX THAT IS ORIGINALLY N X K COMPLEX MATRIX IN A FULL FORM
! C >> M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:), INTENT(IN)    :: A
      TYPE(CMPX_BNDM), INTENT(IN)              :: B
      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: C

      C = MUL_C_CB(RESHAPE(A, (/SIZE(A), 1/)), B)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_CB_C(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> BAND MATRIX THAT IS ORIGINALLY M X N COMPLEX MATRIX IN A FULL FORM
! B >> N X K COMPLEX MATRIX
! C >> M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(CMPX_BNDM), INTENT(IN)              :: A
      COMPLEX(P8), DIMENSION(:,:), INTENT(IN)  :: B
      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: C

      C = TRANSPOSE(MUL_C_CB(TRANSPOSE(B),BAND_TRANSPOSE(A)))

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_CB_C1(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> BAND MATRIX THAT IS ORIGINALLY M X N COMPLEX MATRIX IN A FULL FORM
! B >> N X 1 COMPLEX MATRIX (VECTOR)
! C >> M X 1 COMPLEX MATRIX (VECTOR) CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(CMPX_BNDM), INTENT(IN)            :: A
      COMPLEX(P8), DIMENSION(:), INTENT(IN)  :: B
      COMPLEX(P8), DIMENSION(:), ALLOCATABLE :: C

      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: C2

      C2 = MUL_CB_C(A, RESHAPE(B, (/SIZE(B), 1/)))
      C = RESHAPE(C2, (/SIZE(C2)/))

      DEALLOCATE( C2 )

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_CB_CB(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> BAND MATRIX THAT IS ORIGINALLY M X N COMPLEX MATRIX IN A FULL FORM
! B >> BAND MATRIX THAT IS ORIGINALLY N X K COMPLEX MATRIX IN A FULL FORM
! C >> BAND MATRIX OF THE M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(CMPX_BNDM), INTENT(IN)              :: A
      TYPE(CMPX_BNDM), INTENT(IN)              :: B
      TYPE(CMPX_BNDM)                          :: C

      INTEGER                                  :: N, NI, NJ, I, J, K
      INTEGER                                  :: A0, A1, CT
      COMPLEX(P8)                              :: C0
      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: CE ! TEMPORARY STORAGE W/ REDUNDANT DIAGONALS

      C%SUPLEN = A%SUPLEN + B%SUPLEN
      C%SUBLEN = A%SUBLEN + B%SUBLEN

      N = SIZE(A%E, 2)

      NI = SIZE(B%E, 1)
      NJ = SIZE(B%E, 2)

      I = COUNT(B%E(NI,:) .EQ. HUGE(1.D0))
      J = COUNT(B%E(:,1) .EQ. HUGE(1.D0))

      IF (N .LT. (NI-I)+(NJ-J)-1) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ELSE
        N = (NI-I)+(NJ-J)-1
      ENDIF

      ALLOCATE( CE(C%SUPLEN+C%SUBLEN+1, NJ) )

      CE = HUGE(1.D0)

      DO I = 1, SIZE(CE, 1)
        DO J = 1, SIZE(CE, 2)
          A0 = MAX(1, J-B%SUPLEN)
          A1 = MIN(N, J+B%SUBLEN)
          C0 = 0.D0
          CT = 0
          DO K = A0, A1
            IF ((I+J-B%SUPLEN-K .GE. 1) .AND. (I+J-B%SUPLEN-K .LE. SIZE(A%E, 1)) .AND. &
                (B%SUPLEN+1+K-J .GE. 1) .AND. (B%SUPLEN+1+K-J .LE. SIZE(B%E, 1))) THEN
              IF ((A%E(I+J-B%SUPLEN-K,K) .NE. HUGE(1.D0)) .AND. &
                  (B%E(B%SUPLEN+1+K-J,J) .NE. HUGE(1.D0))) THEN
                C0 = C0 + A%E(I+J-B%SUPLEN-K,K) * B%E(B%SUPLEN+1+K-J,J)
                CT = CT + 1
              ENDIF
            ENDIF
          ENDDO
          IF (CT .NE. 0) THEN
            CE(I,J) = C0
          ENDIF
        ENDDO
      ENDDO

      I = 1 ! TRUNCATE UNNECESSARY SUPERDIAGONALS IN CE
      DO WHILE (.TRUE.)
        IF (ALL(CE(I,:) .EQ. HUGE(1.D0))) THEN
          IF (C%SUPLEN .EQ. 0) EXIT
          C%SUPLEN = C%SUPLEN - 1
          I = I+1
        ELSE
          EXIT
        ENDIF
      ENDDO

      J = SIZE(CE, 1) ! TRUNCATE UNNECESSARY SUBDIAGONALS IN CE
      DO WHILE (.TRUE.)
        IF (ALL(CE(J,:) .EQ. HUGE(1.D0))) THEN
          IF (C%SUBLEN .EQ. 0) EXIT
          C%SUBLEN = C%SUBLEN - 1
          J = J-1
        ELSE
          EXIT
        ENDIF
      ENDDO

      ALLOCATE( C%E(J-I+1, NJ) )
      C%E = CE(I:J,:)

      DEALLOCATE( CE )

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_C_R(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X N COMPLEX MATRIX
! B >> N X K REAL MATRIX
! C >> M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:,:), INTENT(IN)     :: A
      REAL(P8), DIMENSION(:,:), INTENT(IN)        :: B
      COMPLEX(P8), DIMENSION(SIZE(A,1),SIZE(B,2)) :: C

      IF (SIZE(A,2) .NE. SIZE(B,1)) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ENDIF

      C = MATMUL(A, B)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_C_R1(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X N COMPLEX MATRIX
! B >> N X 1 REAL 1-D MATRIX (VECTOR)
! C >> M X 1 COMPLEX 1-D MATRIX (VECTOR) CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:,:), INTENT(IN) :: A
      REAL(P8), DIMENSION(:), INTENT(IN)      :: B
      COMPLEX(P8), DIMENSION(SIZE(A,1))       :: C

      IF (SIZE(A,2) .NE. SIZE(B,1)) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ENDIF

      C = RESHAPE(MATMUL(A, RESHAPE(B, (/SIZE(B), 1/))), (/SIZE(A, 1)/))

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_C1_R(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X 1 COMPLEX 1-D MATRIX (VECTOR)
! B >> N X K REAL MATRIX
! C >> M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:), INTENT(IN)        :: A
      REAL(P8), DIMENSION(:,:), INTENT(IN)         :: B
      COMPLEX(P8), DIMENSION(SIZE(A,1), SIZE(B,2)) :: C

      IF (1 .NE. SIZE(B,1)) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ENDIF

      C = MATMUL(RESHAPE(A, (/SIZE(A), 1/)), B)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_C1_R1(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X 1 COMPLEX 1-D MATRIX (VECTOR)
! R >> N X 1 REAL 1-D MATRIX (VECTOR)
! C >> M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:), INTENT(IN)   :: A
      REAL(P8), DIMENSION(:), INTENT(IN)      :: B
      COMPLEX(P8), DIMENSION(SIZE(A,1))       :: C

      IF (1 .NE. SIZE(B,1)) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ENDIF

      C = RESHAPE(MATMUL(RESHAPE(A, (/SIZE(A), 1/)), RESHAPE(B, (/SIZE(B), 1/))), (/SIZE(A, 1)/))

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_C_RB(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X N COMPLEX MATRIX
! B >> BAND MATRIX THAT IS ORIGINALLY N X K REAL MATRIX IN A FULL FORM
! C >> M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:,:), INTENT(IN)  :: A
      TYPE(REAL_BNDM), INTENT(IN)              :: B
      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: C

      INTEGER                                  :: N, NI, NJ, I, J
      INTEGER                                  :: A0, A1, B0, B1

      N = SIZE(A, 2)

      NI = SIZE(B%E, 1)
      NJ = SIZE(B%E, 2)

      I = COUNT(B%E(NI,:) .EQ. HUGE(1.D0))
      J = COUNT(B%E(:,1) .EQ. HUGE(1.D0))

      IF (N .LT. (NI-I)+(NJ-J)-1) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ELSE
        N = (NI-I)+(NJ-J)-1
      ENDIF

      ALLOCATE( C(SIZE(A, 1), SIZE(B%E, 2)) )
      C = 0.D0

      DO I = 1, SIZE(C, 1)
        DO J = 1, SIZE(C, 2)
          A0 = MAX(1, J-B%SUPLEN)
          A1 = MIN(N, J+B%SUBLEN)
          B0 = MAX(  2+B%SUPLEN-J, 1)
          B1 = MIN(N+1+B%SUPLEN-J, N)
          C(I,J) = SUM(A(I,A0:A1)*B%E(B0:B1,J))
        ENDDO
      ENDDO

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_C1_RB(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X 1 COMPLEX MATRIX (VECTOR)
! B >> BAND MATRIX THAT IS ORIGINALLY N X K REAL MATRIX IN A FULL FORM
! C >> M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:), INTENT(IN)    :: A
      TYPE(REAL_BNDM), INTENT(IN)              :: B
      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: C

      C = MUL_C_RB(RESHAPE(A, (/SIZE(A), 1/)), B)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_CB_R(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> BAND MATRIX THAT IS ORIGINALLY M X N COMPLEX MATRIX IN A FULL FORM
! B >> N X K REAL MATRIX
! C >> M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(CMPX_BNDM), INTENT(IN)              :: A
      REAL(P8), DIMENSION(:,:), INTENT(IN)     :: B
      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: C

      C = TRANSPOSE(MUL_R_CB(TRANSPOSE(B),BAND_TRANSPOSE(A)))

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_CB_R1(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> BAND MATRIX THAT IS ORIGINALLY M X N COMPLEX MATRIX IN A FULL FORM
! B >> N X 1 REAL MATRIX (VECTOR)
! C >> M X 1 COMPLEX MATRIX (VECTOR) CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(CMPX_BNDM), INTENT(IN)            :: A
      REAL(P8), DIMENSION(:), INTENT(IN)     :: B
      COMPLEX(P8), DIMENSION(:), ALLOCATABLE :: C

      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: C2

      C2 = MUL_CB_R(A, RESHAPE(B, (/SIZE(B), 1/)))
      C = RESHAPE(C2, (/SIZE(C2)/))

      DEALLOCATE( C2 )

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_CB_RB(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> BAND MATRIX THAT IS ORIGINALLY M X N COMPLEX MATRIX IN A FULL FORM
! B >> BAND MATRIX THAT IS ORIGINALLY N X K REAL MATRIX IN A FULL FORM
! C >> BAND MATRIX OF THE M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(CMPX_BNDM), INTENT(IN)              :: A
      TYPE(REAL_BNDM), INTENT(IN)              :: B
      TYPE(CMPX_BNDM)                          :: C

      INTEGER                                  :: N, NI, NJ, I, J, K
      INTEGER                                  :: A0, A1, CT
      COMPLEX(P8)                              :: C0
      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: CE ! TEMPORARY STORAGE W/ REDUNDANT DIAGONALS

      C%SUPLEN = A%SUPLEN + B%SUPLEN
      C%SUBLEN = A%SUBLEN + B%SUBLEN

      N = SIZE(A%E, 2)

      NI = SIZE(B%E, 1)
      NJ = SIZE(B%E, 2)

      I = COUNT(B%E(NI,:) .EQ. HUGE(1.D0))
      J = COUNT(B%E(:,1) .EQ. HUGE(1.D0))

      IF (N .LT. (NI-I)+(NJ-J)-1) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ELSE
        N = (NI-I)+(NJ-J)-1
      ENDIF

      ALLOCATE( CE(C%SUPLEN+C%SUBLEN+1, NJ) )

      CE = HUGE(1.D0)

      DO I = 1, SIZE(CE, 1)
        DO J = 1, SIZE(CE, 2)
          A0 = MAX(1, J-B%SUPLEN)
          A1 = MIN(N, J+B%SUBLEN)
          C0 = 0.D0
          CT = 0
          DO K = A0, A1
            IF ((I+J-B%SUPLEN-K .GE. 1) .AND. (I+J-B%SUPLEN-K .LE. SIZE(A%E, 1)) .AND. &
                (B%SUPLEN+1+K-J .GE. 1) .AND. (B%SUPLEN+1+K-J .LE. SIZE(B%E, 1))) THEN
              IF ((A%E(I+J-B%SUPLEN-K,K) .NE. HUGE(1.D0)) .AND. &
                  (B%E(B%SUPLEN+1+K-J,J) .NE. HUGE(1.D0))) THEN
                C0 = C0 + A%E(I+J-B%SUPLEN-K,K) * B%E(B%SUPLEN+1+K-J,J)
                CT = CT + 1
              ENDIF
            ENDIF
          ENDDO
          IF (CT .NE. 0) THEN
            CE(I,J) = C0
          ENDIF
        ENDDO
      ENDDO

      I = 1 ! TRUNCATE UNNECESSARY SUPERDIAGONALS IN CE
      DO WHILE (.TRUE.)
        IF (ALL(CE(I,:) .EQ. HUGE(1.D0))) THEN
          IF (C%SUPLEN .EQ. 0) EXIT
          C%SUPLEN = C%SUPLEN - 1
          I = I+1
        ELSE
          EXIT
        ENDIF
      ENDDO

      J = SIZE(CE, 1) ! TRUNCATE UNNECESSARY SUBDIAGONALS IN CE
      DO WHILE (.TRUE.)
        IF (ALL(CE(J,:) .EQ. HUGE(1.D0))) THEN
          IF (C%SUBLEN .EQ. 0) EXIT
          C%SUBLEN = C%SUBLEN - 1
          J = J-1
        ELSE
          EXIT
        ENDIF
      ENDDO

      ALLOCATE( C%E(J-I+1, NJ) )
      C%E = CE(I:J,:)

      DEALLOCATE( CE )

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_R_C(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X N REAL MATRIX
! B >> N X K COMPLEX MATRIX
! C >> M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:,:), INTENT(IN)        :: A
      COMPLEX(P8), DIMENSION(:,:), INTENT(IN)     :: B
      COMPLEX(P8), DIMENSION(SIZE(A,1),SIZE(B,2)) :: C

      IF (SIZE(A,2) .NE. SIZE(B,1)) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ENDIF

      C = MATMUL(A, B)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_R_C1(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X N REAL MATRIX
! B >> N X 1 COMPLEX 1-D MATRIX (VECTOR)
! C >> M X 1 COMPLEX 1-D MATRIX (VECTOR) CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:,:), INTENT(IN)    :: A
      COMPLEX(P8), DIMENSION(:), INTENT(IN)   :: B
      COMPLEX(P8), DIMENSION(SIZE(A,1))       :: C

      IF (SIZE(A,2) .NE. SIZE(B,1)) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ENDIF

      C = RESHAPE(MATMUL(A, RESHAPE(B, (/SIZE(B), 1/))), (/SIZE(A, 1)/))

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_R1_C(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X 1 REAL 1-D MATRIX (VECTOR)
! B >> N X K COMPLEX MATRIX
! C >> M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:), INTENT(IN)           :: A
      COMPLEX(P8), DIMENSION(:,:), INTENT(IN)      :: B
      COMPLEX(P8), DIMENSION(SIZE(A,1), SIZE(B,2)) :: C

      IF (1 .NE. SIZE(B,1)) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ENDIF

      C = MATMUL(RESHAPE(A, (/SIZE(A), 1/)), B)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_R1_C1(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X 1 REAL 1-D MATRIX (VECTOR)
! R >> N X 1 COMPLEX 1-D MATRIX (VECTOR)
! C >> M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:), INTENT(IN)      :: A
      COMPLEX(P8), DIMENSION(:), INTENT(IN)   :: B
      COMPLEX(P8), DIMENSION(SIZE(A,1))       :: C

      IF (1 .NE. SIZE(B,1)) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ENDIF

      C = RESHAPE(MATMUL(RESHAPE(A, (/SIZE(A), 1/)), RESHAPE(B, (/SIZE(B), 1/))), (/SIZE(A, 1)/))

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_R_CB(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X N REAL MATRIX
! B >> BAND MATRIX THAT IS ORIGINALLY N X K COMPLEX MATRIX IN A FULL FORM
! C >> M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:,:), INTENT(IN)     :: A
      TYPE(CMPX_BNDM), INTENT(IN)              :: B
      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: C

      INTEGER                                  :: N, NI, NJ, I, J
      INTEGER                                  :: A0, A1, B0, B1

      N = SIZE(A, 2)

      NI = SIZE(B%E, 1)
      NJ = SIZE(B%E, 2)

      I = COUNT(B%E(NI,:) .EQ. HUGE(1.D0))
      J = COUNT(B%E(:,1) .EQ. HUGE(1.D0))

      IF (N .LT. (NI-I)+(NJ-J)-1) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ELSE
        N = (NI-I)+(NJ-J)-1
      ENDIF

      ALLOCATE( C(SIZE(A, 1), SIZE(B%E, 2)) )
      C = 0.D0

      DO I = 1, SIZE(C, 1)
        DO J = 1, SIZE(C, 2)
          A0 = MAX(1, J-B%SUPLEN)
          A1 = MIN(N, J+B%SUBLEN)
          B0 = MAX(  2+B%SUPLEN-J, 1)
          B1 = MIN(N+1+B%SUPLEN-J, N)
          C(I,J) = SUM(A(I,A0:A1)*B%E(B0:B1,J))
        ENDDO
      ENDDO

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_R1_CB(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X 1 REAL MATRIX (VECTOR)
! B >> BAND MATRIX THAT IS ORIGINALLY N X K COMPLEX MATRIX IN A FULL FORM
! C >> M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:), INTENT(IN)       :: A
      TYPE(CMPX_BNDM), INTENT(IN)              :: B
      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: C

      C = MUL_R_CB(RESHAPE(A, (/SIZE(A), 1/)), B)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_RB_C(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> BAND MATRIX THAT IS ORIGINALLY M X N REAL MATRIX IN A FULL FORM
! B >> N X K COMPLEX MATRIX
! C >> M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(REAL_BNDM), INTENT(IN)              :: A
      COMPLEX(P8), DIMENSION(:,:), INTENT(IN)  :: B
      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: C

      C = TRANSPOSE(MUL_C_RB(TRANSPOSE(B),BAND_TRANSPOSE(A)))

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_RB_C1(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> BAND MATRIX THAT IS ORIGINALLY M X N REAL MATRIX IN A FULL FORM
! B >> N X 1 COMPLEX MATRIX (VECTOR)
! C >> M X 1 COMPLEX MATRIX (VECTOR) CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(REAL_BNDM), INTENT(IN)            :: A
      COMPLEX(P8), DIMENSION(:), INTENT(IN)     :: B
      COMPLEX(P8), DIMENSION(:), ALLOCATABLE :: C

      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: C2

      C2 = MUL_RB_C(A, RESHAPE(B, (/SIZE(B), 1/)))
      C = RESHAPE(C2, (/SIZE(C2)/))

      DEALLOCATE( C2 )

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_RB_CB(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> BAND MATRIX THAT IS ORIGINALLY M X N REAL MATRIX IN A FULL FORM
! B >> BAND MATRIX THAT IS ORIGINALLY N X K COMPLEX MATRIX IN A FULL FORM
! C >> BAND MATRIX OF THE M X K COMPLEX MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(REAL_BNDM), INTENT(IN)              :: A
      TYPE(CMPX_BNDM), INTENT(IN)              :: B
      TYPE(CMPX_BNDM)                          :: C

      INTEGER                                  :: N, NI, NJ, I, J, K
      INTEGER                                  :: A0, A1, CT
      COMPLEX(P8)                              :: C0
      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: CE ! TEMPORARY STORAGE W/ REDUNDANT DIAGONALS

      C%SUPLEN = A%SUPLEN + B%SUPLEN
      C%SUBLEN = A%SUBLEN + B%SUBLEN

      N = SIZE(A%E, 2)

      NI = SIZE(B%E, 1)
      NJ = SIZE(B%E, 2)

      I = COUNT(B%E(NI,:) .EQ. HUGE(1.D0))
      J = COUNT(B%E(:,1) .EQ. HUGE(1.D0))

      IF (N .LT. (NI-I)+(NJ-J)-1) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ELSE
        N = (NI-I)+(NJ-J)-1
      ENDIF

      ALLOCATE( CE(C%SUPLEN+C%SUBLEN+1, NJ) )

      CE = HUGE(1.D0)

      DO I = 1, SIZE(CE, 1)
        DO J = 1, SIZE(CE, 2)
          A0 = MAX(1, J-B%SUPLEN)
          A1 = MIN(N, J+B%SUBLEN)
          C0 = 0.D0
          CT = 0
          DO K = A0, A1
            IF ((I+J-B%SUPLEN-K .GE. 1) .AND. (I+J-B%SUPLEN-K .LE. SIZE(A%E, 1)) .AND. &
                (B%SUPLEN+1+K-J .GE. 1) .AND. (B%SUPLEN+1+K-J .LE. SIZE(B%E, 1))) THEN
              IF ((A%E(I+J-B%SUPLEN-K,K) .NE. HUGE(1.D0)) .AND. &
                  (B%E(B%SUPLEN+1+K-J,J) .NE. HUGE(1.D0))) THEN
                C0 = C0 + A%E(I+J-B%SUPLEN-K,K) * B%E(B%SUPLEN+1+K-J,J)
                CT = CT + 1
              ENDIF
            ENDIF
          ENDDO
          IF (CT .NE. 0) THEN
            CE(I,J) = C0
          ENDIF
        ENDDO
      ENDDO

      I = 1 ! TRUNCATE UNNECESSARY SUPERDIAGONALS IN CE
      DO WHILE (.TRUE.)
        IF (ALL(CE(I,:) .EQ. HUGE(1.D0))) THEN
          IF (C%SUPLEN .EQ. 0) EXIT
          C%SUPLEN = C%SUPLEN - 1
          I = I+1
        ELSE
          EXIT
        ENDIF
      ENDDO

      J = SIZE(CE, 1) ! TRUNCATE UNNECESSARY SUBDIAGONALS IN CE
      DO WHILE (.TRUE.)
        IF (ALL(CE(J,:) .EQ. HUGE(1.D0))) THEN
          IF (C%SUBLEN .EQ. 0) EXIT
          C%SUBLEN = C%SUBLEN - 1
          J = J-1
        ELSE
          EXIT
        ENDIF
      ENDDO

      ALLOCATE( C%E(J-I+1, NJ) )
      C%E = CE(I:J,:)

      DEALLOCATE( CE )

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_R_R(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X N REAL MATRIX
! B >> N X K REAL MATRIX
! C >> M X K REAL MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:,:), INTENT(IN)        :: A
      REAL(P8), DIMENSION(:,:), INTENT(IN)        :: B
      REAL(P8), DIMENSION(SIZE(A,1),SIZE(B,2))    :: C

      IF (SIZE(A,2) .NE. SIZE(B,1)) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ENDIF

      C = MATMUL(A, B)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_R_R1(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X N REAL MATRIX
! B >> N X 1 REAL 1-D MATRIX (VECTOR)
! C >> M X 1 REAL 1-D MATRIX (VECTOR) CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:,:), INTENT(IN)    :: A
      REAL(P8), DIMENSION(:), INTENT(IN)      :: B
      REAL(P8), DIMENSION(SIZE(A,1))          :: C

      IF (SIZE(A,2) .NE. SIZE(B,1)) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ENDIF

      C = RESHAPE(MATMUL(A, RESHAPE(B, (/SIZE(B), 1/))), (/SIZE(A, 1)/))

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_R1_R(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X 1 REAL 1-D MATRIX (VECTOR)
! B >> N X K REAL MATRIX
! C >> M X K REAL MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:), INTENT(IN)           :: A
      REAL(P8), DIMENSION(:,:), INTENT(IN)         :: B
      REAL(P8), DIMENSION(SIZE(A,1), SIZE(B,2))    :: C

      IF (1 .NE. SIZE(B,1)) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ENDIF

      C = MATMUL(RESHAPE(A, (/SIZE(A), 1/)), B)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_R1_R1(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X 1 REAL 1-D MATRIX (VECTOR)
! R >> N X 1 REAL 1-D MATRIX (VECTOR)
! C >> M X K REAL MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:), INTENT(IN)      :: A
      REAL(P8), DIMENSION(:), INTENT(IN)      :: B
      REAL(P8), DIMENSION(SIZE(A,1))          :: C

      IF (1 .NE. SIZE(B,1)) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ENDIF

      C = RESHAPE(MATMUL(RESHAPE(A, (/SIZE(A), 1/)), RESHAPE(B, (/SIZE(B), 1/))), (/SIZE(A, 1)/))

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_R_RB(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X N REAL MATRIX
! B >> BAND MATRIX THAT IS ORIGINALLY N X K REAL MATRIX IN A FULL FORM
! C >> M X K REAL MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:,:), INTENT(IN)     :: A
      TYPE(REAL_BNDM), INTENT(IN)              :: B
      REAL(P8), DIMENSION(:,:), ALLOCATABLE    :: C

      INTEGER                                  :: N, NI, NJ, I, J
      INTEGER                                  :: A0, A1, B0, B1

      N = SIZE(A, 2)

      NI = SIZE(B%E, 1)
      NJ = SIZE(B%E, 2)

      I = COUNT(B%E(NI,:) .EQ. HUGE(1.D0))
      J = COUNT(B%E(:,1) .EQ. HUGE(1.D0))

      IF (N .LT. (NI-I)+(NJ-J)-1) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ELSE
        N = (NI-I)+(NJ-J)-1
      ENDIF

      ALLOCATE( C(SIZE(A, 1), SIZE(B%E, 2)) )
      C = 0.D0

      DO I = 1, SIZE(C, 1)
        DO J = 1, SIZE(C, 2)
          A0 = MAX(1, J-B%SUPLEN)
          A1 = MIN(N, J+B%SUBLEN)
          B0 = MAX(  2+B%SUPLEN-J, 1)
          B1 = MIN(N+1+B%SUPLEN-J, N)
          C(I,J) = SUM(A(I,A0:A1)*B%E(B0:B1,J))
        ENDDO
      ENDDO

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_R1_RB(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> M X 1 REAL MATRIX (VECTOR)
! B >> BAND MATRIX THAT IS ORIGINALLY N X K REAL MATRIX IN A FULL FORM
! C >> M X K REAL MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:), INTENT(IN)       :: A
      TYPE(REAL_BNDM), INTENT(IN)              :: B
      REAL(P8), DIMENSION(:,:), ALLOCATABLE    :: C

      C = MUL_R_RB(RESHAPE(A, (/SIZE(A), 1/)), B)

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_RB_R(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> BAND MATRIX THAT IS ORIGINALLY M X N REAL MATRIX IN A FULL FORM
! B >> N X K REAL MATRIX
! C >> M X K REAL MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(REAL_BNDM), INTENT(IN)              :: A
      REAL(P8), DIMENSION(:,:), INTENT(IN)     :: B
      REAL(P8), DIMENSION(:,:), ALLOCATABLE    :: C

      C = TRANSPOSE(MUL_R_RB(TRANSPOSE(B),BAND_TRANSPOSE(A)))

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_RB_R1(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> BAND MATRIX THAT IS ORIGINALLY M X N REAL MATRIX IN A FULL FORM
! B >> N X 1 REAL MATRIX (VECTOR)
! C >> M X 1 REAL MATRIX (VECTOR) CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(REAL_BNDM), INTENT(IN)            :: A
      REAL(P8), DIMENSION(:), INTENT(IN)     :: B
      REAL(P8), DIMENSION(:), ALLOCATABLE    :: C

      REAL(P8), DIMENSION(:,:), ALLOCATABLE  :: C2

      C2 = MUL_RB_R(A, RESHAPE(B, (/SIZE(B), 1/)))
      C = RESHAPE(C2, (/SIZE(C2)/))

      DEALLOCATE( C2 )

      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION MUL_RB_RB(A, B) RESULT(C)
! ==========================================================================================================
! [USAGE]:
! UNIFYING MATRIX MULTIPLICATION OPERATOR EXPRESSION AS C = A .MUL. B
! [VARIABLES]:
! A >> BAND MATRIX THAT IS ORIGINALLY M X N REAL MATRIX IN A FULL FORM
! B >> BAND MATRIX THAT IS ORIGINALLY N X K REAL MATRIX IN A FULL FORM
! C >> BAND MATRIX OF THE M X K REAL MATRIX CORRESPONDING TO A .MUL. B
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(REAL_BNDM), INTENT(IN)              :: A
      TYPE(REAL_BNDM), INTENT(IN)              :: B
      TYPE(REAL_BNDM)                          :: C

      INTEGER                                  :: N, NI, NJ, I, J, K
      INTEGER                                  :: A0, A1, CT
      REAL(P8)                                 :: C0
      REAL(P8), DIMENSION(:,:), ALLOCATABLE    :: CE ! TEMPORARY STORAGE W/ REDUNDANT DIAGONALS

      C%SUPLEN = A%SUPLEN + B%SUPLEN
      C%SUBLEN = A%SUBLEN + B%SUBLEN

      N = SIZE(A%E, 2)

      NI = SIZE(B%E, 1)
      NJ = SIZE(B%E, 2)

      I = COUNT(B%E(NI,:) .EQ. HUGE(1.D0))
      J = COUNT(B%E(:,1) .EQ. HUGE(1.D0))

      IF (N .LT. (NI-I)+(NJ-J)-1) THEN
        WRITE(*,*) '.MUL.: MATRIX SIZE INCOMPATIBLE.'
        STOP
      ELSE
        N = (NI-I)+(NJ-J)-1
      ENDIF

      ALLOCATE( CE(C%SUPLEN+C%SUBLEN+1, NJ) )

      CE = HUGE(1.D0)

      DO I = 1, SIZE(CE, 1)
        DO J = 1, SIZE(CE, 2)
          A0 = MAX(1, J-B%SUPLEN)
          A1 = MIN(N, J+B%SUBLEN)
          C0 = 0.D0
          CT = 0
          DO K = A0, A1
            IF ((I+J-B%SUPLEN-K .GE. 1) .AND. (I+J-B%SUPLEN-K .LE. SIZE(A%E, 1)) .AND. &
                (B%SUPLEN+1+K-J .GE. 1) .AND. (B%SUPLEN+1+K-J .LE. SIZE(B%E, 1))) THEN
              IF ((A%E(I+J-B%SUPLEN-K,K) .NE. HUGE(1.D0)) .AND. &
                  (B%E(B%SUPLEN+1+K-J,J) .NE. HUGE(1.D0))) THEN
                C0 = C0 + A%E(I+J-B%SUPLEN-K,K) * B%E(B%SUPLEN+1+K-J,J)
                CT = CT + 1
              ENDIF
            ENDIF
          ENDDO
          IF (CT .NE. 0) THEN
            CE(I,J) = C0
          ENDIF
        ENDDO
      ENDDO

      I = 1 ! TRUNCATE UNNECESSARY SUPERDIAGONALS IN CE
      DO WHILE (.TRUE.)
        IF (ALL(CE(I,:) .EQ. HUGE(1.D0))) THEN
          IF (C%SUPLEN .EQ. 0) EXIT
          C%SUPLEN = C%SUPLEN - 1
          I = I+1
        ELSE
          EXIT
        ENDIF
      ENDDO

      J = SIZE(CE, 1) ! TRUNCATE UNNECESSARY SUBDIAGONALS IN CE
      DO WHILE (.TRUE.)
        IF (ALL(CE(J,:) .EQ. HUGE(1.D0))) THEN
          IF (C%SUBLEN .EQ. 0) EXIT
          C%SUBLEN = C%SUBLEN - 1
          J = J-1
        ELSE
          EXIT
        ENDIF
      ENDDO

      ALLOCATE( C%E(J-I+1, NJ) )
      C%E = CE(I:J,:)

      DEALLOCATE( CE )

      RETURN
      END FUNCTION
! ==========================================================================================================
      SUBROUTINE LU_C(A, IPIV)
! ==========================================================================================================
! [USAGE]:
! PERFORM THE LU FACTORIZATION WITH RESPECT TO AN INPUT M X N COMPLEX MATRIX A. 
! ON EXIT, A IS REPLACED WITH L+U-I.
! [VARIABLES]:
! A >> ON ENTRY, THE M X N MATRIX TO BE FACTORED. ON EXIT, THE FACTORS L AND U VIA FACTORIAZATION (L+U-I)
! IPIV >> WORKSPACE ARRAY, SIZE AT LEAST MAX(1, MIN(#ROWS, #COLS))
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT) :: A
      INTEGER, DIMENSION(:), INTENT(OUT)         :: IPIV

      INTEGER                                    :: M, N, INFO

      EXTERNAL ZGETRF

      M = SIZE(A,1)
      N = SIZE(A,2)

      IF (SIZE(IPIV) .NE. MIN(M, N)) THEN
        WRITE(*,*) 'LU: INCONSISTENCY IN IPIV`S DIMENSION.'
        STOP
      ENDIF

      CALL ZGETRF(M, N, A, M, IPIV, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'LU: LU FACTORIZATION FAILED.'
        STOP
      ENDIF

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE LU_CB(A, IPIV)
! ==========================================================================================================
! [USAGE]:
! PERFORM THE LU FACTORIZATION WITH RESPECT TO AN INPUT M X N COMPLEX MATRIX A IN A BAND FORM. 
! ON EXIT, A IS REPLACED WITH L+U-I.
! [VARIABLES]:
! A >> ON ENTRY, THE M X N MATRIX TO BE FACTORED. ON EXIT, THE FACTORS L AND U VIA FACTORIAZATION (L+U-I)
! IPIV >> WORKSPACE ARRAY, SIZE AT LEAST MAX(1, MIN(#ROWS, #COLS))
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(CMPX_BNDM), INTENT(INOUT)           :: A
      INTEGER, DIMENSION(:), INTENT(OUT)       :: IPIV

      INTEGER                                  :: M, N, KU, KL, INFO, NI, NJ, I, J
      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: AE

      EXTERNAL ZGBTRF

      KU = A%SUPLEN
      KL = A%SUBLEN

      NI = SIZE(A%E, 1)
      NJ = SIZE(A%E, 2)

      I = COUNT(A%E(NI,:) .EQ. HUGE(1.D0))
      J = COUNT(A%E(:,1) .EQ. HUGE(1.D0))

      M = (NI-I)+(NJ-J)-1
      N = NJ

      ALLOCATE( AE(2*KL+KU+1, N) )
      WHERE (A%E .NE. HUGE(1.D0)) AE(KL+1:,:) = A%E

      IF (SIZE(IPIV) .NE. MIN(M, N)) THEN
        WRITE(*,*) 'LU: INCONSISTENCY IN IPIV`S DIMENSION.'
        STOP
      ENDIF

      CALL ZGBTRF(M, N, KL, KU, AE, 2*KL+KU+1, IPIV, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'LU: LU FACTORIZATION FAILED.'
        STOP
      ENDIF

      CALL MOVE_ALLOC(AE, A%E)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE LU_R(A, IPIV)
! ==========================================================================================================
! [USAGE]:
! PERFORM THE LU FACTORIZATION WITH RESPECT TO AN INPUT M X N REAL MATRIX A. 
! ON EXIT, A IS REPLACED WITH L+U-I.
! [VARIABLES]:
! A >> ON ENTRY, THE M X N MATRIX TO BE FACTORED. ON EXIT, THE FACTORS L AND U VIA FACTORIAZATION (L+U-I)
! IPIV >> WORKSPACE ARRAY, SIZE AT LEAST MAX(1, MIN(#ROWS, #COLS))
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:,:), INTENT(INOUT)    :: A
      INTEGER, DIMENSION(:), INTENT(OUT)         :: IPIV

      INTEGER                                    :: M, N, INFO

      EXTERNAL DGETRF

      M = SIZE(A,1)
      N = SIZE(A,2)

      IF (SIZE(IPIV) .NE. MIN(M, N)) THEN
        WRITE(*,*) 'LU: INCONSISTENCY IN IPIV`S DIMENSION.'
        STOP
      ENDIF

      CALL DGETRF(M, N, A, M, IPIV, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'LU: LU FACTORIZATION FAILED.'
        STOP
      ENDIF

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE LU_RB(A, IPIV)
! ==========================================================================================================
! [USAGE]:
! PERFORM THE LU FACTORIZATION WITH RESPECT TO AN INPUT M X N REAL MATRIX A IN A BAND FORM. 
! ON EXIT, A IS REPLACED WITH L+U-I.
! [VARIABLES]:
! A >> ON ENTRY, THE M X N MATRIX TO BE FACTORED. ON EXIT, THE FACTORS L AND U VIA FACTORIAZATION (L+U-I)
! IPIV >> WORKSPACE ARRAY, SIZE AT LEAST MAX(1, MIN(#ROWS, #COLS))
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(REAL_BNDM), INTENT(INOUT)        :: A
      INTEGER, DIMENSION(:), INTENT(OUT)    :: IPIV

      INTEGER                               :: M, N, KU, KL, INFO, NI, NJ, I, J
      REAL(P8), DIMENSION(:,:), ALLOCATABLE :: AE

      EXTERNAL DGBTRF

      KU = A%SUPLEN
      KL = A%SUBLEN

      NI = SIZE(A%E, 1)
      NJ = SIZE(A%E, 2)

      I = COUNT(A%E(NI,:) .EQ. HUGE(1.D0))
      J = COUNT(A%E(:,1) .EQ. HUGE(1.D0))

      M = (NI-I)+(NJ-J)-1
      N = NJ

      ALLOCATE( AE(2*KL+KU+1, N) )
      WHERE (A%E .NE. HUGE(1.D0)) AE(KL+1:,:) = A%E

      IF (SIZE(IPIV) .NE. N) THEN
        WRITE(*,*) 'LU: INCONSISTENCY IN IPIV`S DIMENSION.'
        STOP
      ENDIF

      CALL DGBTRF(M, N, KL, KU, AE, 2*KL+KU+1, IPIV, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'LU: LU FACTORIZATION FAILED.'
        STOP
      ENDIF

      CALL MOVE_ALLOC(AE, A%E)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE SOLVE_C_C(A,B)
! ==========================================================================================================
! [USAGE]:
! SOLVE THE SYSTEM(S) OF LINEAR EQUATIONS AX = B, I.E., X = A^(-1)B
! WHERE A NEEDS TO BE LU-FACTORIZED IN ADVANCE AS L+U-I
! [VARIABLES]:
! A >> N X N SQUARE COMPLEX MATRIX AFTER LU FACTORIZATION. IF NOT FACTORIZED, IT WILL BE FACTORIZED ON EXIT
! B >> ON ENTRY, COMPLEX N X 1 RHS VECTOR OF THE SYSTEM. ON EXIT, THE SOLUTION WILL BE RECORDED.
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT)     :: A
      COMPLEX(P8), DIMENSION(:), INTENT(INOUT)       :: B

      INTEGER, DIMENSION(MIN(SIZE(A, 1),SIZE(A, 2))) :: IPIV
      INTEGER                                        :: NI, NJ, INFO

      EXTERNAL ZGETRS

      IF (SIZE(A, 1) .NE. SIZE(A, 2)) THEN
        WRITE (*,*) 'SOLVE: MATRIX NOT SQUARE'
        WRITE (*,*) SIZE(A,1),'X',SIZE(A,2)
        STOP
      ENDIF

      IF ( SIZE(A, 2) .NE. SIZE(B) ) THEN
        WRITE (*,*) 'SOLVE: MATRIX SIZE INCONSISTENT'
        STOP
      ENDIF

      NI = SIZE(A, 1)
      NJ = 1

      CALL LU_C(A, IPIV)
      CALL ZGETRS('N', NI, NJ, A, NI, IPIV, B, NI, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'SOLVE: LINEAR SYSTEM UNABLE TO BE SOLVED.'
        STOP
      ENDIF

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE SOLVE_C_C_MULTI(A,B)
! ==========================================================================================================
! [USAGE]:
! SOLVE THE SYSTEM(S) OF LINEAR EQUATIONS AX = B, I.E., X = A^(-1)B
! WHERE A NEEDS TO BE LU-FACTORIZED IN ADVANCE AS L+U-I
! [VARIABLES]:
! A >> N X N SQUARE COMPLEX MATRIX AFTER LU FACTORIZATION. IF NOT FACTORIZED, IT WILL BE FACTORIZED ON EXIT
! B >> ON ENTRY, COMPLEX N X K RHS MATRIX OF THE SYSTEM. ON EXIT, THE SOLUTION WILL BE RECORDED.
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT)     :: A
      COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT)     :: B

      INTEGER, DIMENSION(MIN(SIZE(A, 1),SIZE(A, 2))) :: IPIV
      INTEGER                                        :: NI, NJ, INFO

      EXTERNAL ZGETRS

      IF (SIZE(A, 1) .NE. SIZE(A, 2)) THEN
        WRITE (*,*) 'SOLVE: MATRIX NOT SQUARE'
        WRITE (*,*) SIZE(A,1),'X',SIZE(A,2)
        STOP
      ENDIF

      IF ( SIZE(A, 2) .NE. SIZE(B, 1) ) THEN
        WRITE (*,*) 'SOLVE: MATRIX SIZE INCONSISTENT'
        STOP
      ENDIF

      NI = SIZE(A, 1)
      NJ = SIZE(B, 2)

      CALL LU_C(A, IPIV)
      CALL ZGETRS('N', NI, NJ, A, NI, IPIV, B, NI, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'SOLVE: LINEAR SYSTEM UNABLE TO BE SOLVED.'
        STOP
      ENDIF

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE SOLVE_CB_C(A,B)
! ==========================================================================================================
! [USAGE]:
! SOLVE THE SYSTEM(S) OF LINEAR EQUATIONS AX = B, I.E., X = A^(-1)B
! WHERE A NEEDS TO BE LU-FACTORIZED IN ADVANCE AS L+U-I
! [VARIABLES]:
! A >> N X N SQUARE COMPLEX MATRIX AFTER LU FACTORIZATION IN A BAND FORMAT.
! B >> ON ENTRY, COMPLEX N X 1 RHS VECTOR OF THE SYSTEM. ON EXIT, THE SOLUTION WILL BE RECORDED.
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(CMPX_BNDM), INTENT(INOUT)           :: A
      COMPLEX(P8), DIMENSION(:), INTENT(INOUT) :: B

      INTEGER, DIMENSION(:), ALLOCATABLE       :: IPIV
      INTEGER                                  :: N, NRHS, KL, KU, INFO, NI, NJ, I, J

      EXTERNAL ZGBTRS

      IF ( SIZE(A%E, 2) .NE. SIZE(B) ) THEN
        WRITE (*,*) 'SOLVE: MATRIX SIZE INCONSISTENT'
        STOP
      ENDIF

      N = SIZE(B)
      NRHS = 1
      KL = A%SUBLEN
      KU = A%SUPLEN

      NI = SIZE(A%E, 1)
      NJ = SIZE(A%E, 2)

      I = COUNT(A%E(NI,:) .EQ. HUGE(1.D0))
      J = COUNT(A%E(:,1) .EQ. HUGE(1.D0))

      ALLOCATE( IPIV(MIN((NI-I)+(NJ-J)-1, NJ)) )

      CALL LU_CB(A, IPIV)
      CALL ZGBTRS('N', N, KL, KU, NRHS, A%E, 2*KL+KU+1, IPIV, B, N, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'SOLVE: LINEAR SYSTEM UNABLE TO BE SOLVED.'
        STOP
      ENDIF

      ! DO I = 1, SIZE(B)
      !   B(I) = CMPLX(DNINT(REAL(B(I))*1.D30)*1.D-30, DNINT(AIMAG(B(I))*1.D30)*1.D-30, P8)
      ! ENDDO

      DEALLOCATE( IPIV )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE SOLVE_CB_C_MULTI(A,B)
! ==========================================================================================================
! [USAGE]:
! SOLVE THE SYSTEM(S) OF LINEAR EQUATIONS AX = B, I.E., X = A^(-1)B
! WHERE A NEEDS TO BE LU-FACTORIZED IN ADVANCE AS L+U-I
! [VARIABLES]:
! A >> N X N SQUARE COMPLEX MATRIX AFTER LU FACTORIZATION IN A BAND FORMAT.
! B >> ON ENTRY, COMPLEX N X K RHS MATRIX OF THE SYSTEM. ON EXIT, THE SOLUTION WILL BE RECORDED.
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(CMPX_BNDM), INTENT(INOUT)             :: A
      COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT) :: B

      INTEGER, DIMENSION(:), ALLOCATABLE         :: IPIV
      INTEGER                                    :: N, NRHS, KL, KU, INFO, NI, NJ, I, J

      EXTERNAL ZGBTRS

      IF ( SIZE(A%E, 2) .NE. SIZE(B, 1) ) THEN
        WRITE (*,*) 'SOLVE: MATRIX SIZE INCONSISTENT'
        STOP
      ENDIF

      N = SIZE(B, 1)
      NRHS = SIZE(B, 2)
      KL = A%SUBLEN
      KU = A%SUPLEN

      NI = SIZE(A%E, 1)
      NJ = SIZE(A%E, 2)

      I = COUNT(A%E(NI,:) .EQ. HUGE(1.D0))
      J = COUNT(A%E(:,1) .EQ. HUGE(1.D0))

      ALLOCATE( IPIV(MIN((NI-I)+(NJ-J)-1, NJ)) )

      CALL LU_CB(A, IPIV)
      CALL ZGBTRS('N', N, KL, KU, NRHS, A%E, 2*KL+KU+1, IPIV, B, N, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'SOLVE: LINEAR SYSTEM UNABLE TO BE SOLVED.'
        STOP
      ENDIF

      ! DO I = 1, SIZE(B,1)
      !   DO J = 1, SIZE(B,2)
      !     B(I,J) = CMPLX(DNINT(REAL(B(I,J))*1.D30)*1.D-30, DNINT(AIMAG(B(I,J))*1.D30)*1.D-30, P8)
      !   ENDDO
      ! ENDDO

      DEALLOCATE( IPIV )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE SOLVE_R_C(A,B)
! ==========================================================================================================
! [USAGE]:
! SOLVE THE SYSTEM(S) OF LINEAR EQUATIONS AX = B, I.E., X = A^(-1)B
! WHERE A NEEDS TO BE LU-FACTORIZED IN ADVANCE AS L+U-I
! [VARIABLES]:
! A >> N X N SQUARE REAL MATRIX AFTER LU FACTORIZATION. IF NOT FACTORIZED, IT WILL BE FACTORIZED ON EXIT
! B >> ON ENTRY, COMPLEX N X 1 RHS VECTOR OF THE SYSTEM. ON EXIT, THE SOLUTION WILL BE RECORDED.
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:,:), INTENT(INOUT)        :: A
      COMPLEX(P8), DIMENSION(:), INTENT(INOUT)       :: B

      REAL(P8), DIMENSION(SIZE(B))                   :: BR, BI
      INTEGER, DIMENSION(MIN(SIZE(A, 1),SIZE(A, 2))) :: IPIV, IPIVC
      INTEGER                                        :: NI, NJ, INFO

      EXTERNAL DGETRS

      IF (SIZE(A, 1) .NE. SIZE(A, 2)) THEN
        WRITE (*,*) 'SOLVE: MATRIX NOT SQUARE'
        WRITE (*,*) SIZE(A,1),'X',SIZE(A,2)
        STOP
      ENDIF

      IF ( SIZE(A, 2) .NE. SIZE(B) ) THEN
        WRITE (*,*) 'SOLVE: MATRIX SIZE INCONSISTENT'
        STOP
      ENDIF

      NI = SIZE(A, 1)
      NJ = 1

      CALL LU_R(A, IPIV)
      IPIVC = IPIV

      BR = REAL(B)
      CALL DGETRS('N', NI, NJ, A, NI, IPIV, BR, NI, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'SOLVE: LINEAR SYSTEM UNABLE TO BE SOLVED.'
        STOP
      ENDIF

      BI = AIMAG(B)
      CALL DGETRS('N', NI, NJ, A, NI, IPIV, BI, NI, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'SOLVE: LINEAR SYSTEM UNABLE TO BE SOLVED.'
        STOP
      ENDIF

      B = BR + BI * IU

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE SOLVE_R_C_MULTI(A,B)
! ==========================================================================================================
! [USAGE]:
! SOLVE THE SYSTEM(S) OF LINEAR EQUATIONS AX = B, I.E., X = A^(-1)B
! WHERE A NEEDS TO BE LU-FACTORIZED IN ADVANCE AS L+U-I
! [VARIABLES]:
! A >> N X N SQUARE REAL MATRIX AFTER LU FACTORIZATION. IF NOT FACTORIZED, IT WILL BE FACTORIZED ON EXIT
! B >> ON ENTRY, COMPLEX N X K RHS MATRIX OF THE SYSTEM. ON EXIT, THE SOLUTION WILL BE RECORDED.
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:,:), INTENT(INOUT)        :: A
      COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT)     :: B

      REAL(P8), DIMENSION(SIZE(B,1), SIZE(B,2))      :: BR, BI
      INTEGER, DIMENSION(MIN(SIZE(A, 1),SIZE(A, 2))) :: IPIV, IPIVC
      INTEGER                                        :: NI, NJ, INFO

      EXTERNAL DGETRS

      IF (SIZE(A, 1) .NE. SIZE(A, 2)) THEN
        WRITE (*,*) 'SOLVE: MATRIX NOT SQUARE'
        WRITE (*,*) SIZE(A,1),'X',SIZE(A,2)
        STOP
      ENDIF

      IF ( SIZE(A, 2) .NE. SIZE(B, 1) ) THEN
        WRITE (*,*) 'SOLVE: MATRIX SIZE INCONSISTENT'
        STOP
      ENDIF

      NI = SIZE(A, 1)
      NJ = SIZE(B, 2)

      CALL LU_R(A, IPIV)
      IPIVC = IPIV

      BR = REAL(B)
      CALL DGETRS('N', NI, NJ, A, NI, IPIV, BR, NI, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'SOLVE: LINEAR SYSTEM UNABLE TO BE SOLVED.'
        STOP
      ENDIF

      BI = AIMAG(B)
      CALL DGETRS('N', NI, NJ, A, NI, IPIVC, BI, NI, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'SOLVE: LINEAR SYSTEM UNABLE TO BE SOLVED.'
        STOP
      ENDIF

      B = BR + BI * IU

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE SOLVE_RB_C(A,B)
! ==========================================================================================================
! [USAGE]:
! SOLVE THE SYSTEM(S) OF LINEAR EQUATIONS AX = B, I.E., X = A^(-1)B
! WHERE A NEEDS TO BE LU-FACTORIZED IN ADVANCE AS L+U-I
! [VARIABLES]:
! A >> N X N SQUARE REAL MATRIX AFTER LU FACTORIZATION IN A BAND FORMAT.
! B >> ON ENTRY, COMPLEX N X 1 RHS VECTOR OF THE SYSTEM. ON EXIT, THE SOLUTION WILL BE RECORDED.
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(REAL_BNDM), INTENT(INOUT)           :: A
      COMPLEX(P8), DIMENSION(:), INTENT(INOUT) :: B

      REAL(P8), DIMENSION(:,:), ALLOCATABLE    :: AE
      REAL(P8), DIMENSION(SIZE(B))             :: BR, BI
      INTEGER, DIMENSION(:), ALLOCATABLE       :: IPIV, IPIVC
      INTEGER                                  :: N, NRHS, KL, KU, INFO, NI, NJ, I, J

      EXTERNAL DGBTRS

      IF ( SIZE(A%E, 2) .NE. SIZE(B) ) THEN
        WRITE (*,*) 'SOLVE: MATRIX SIZE INCONSISTENT'
        STOP
      ENDIF

      N = SIZE(B)
      NRHS = 1
      KL = A%SUBLEN
      KU = A%SUPLEN

      NI = SIZE(A%E, 1)
      NJ = SIZE(A%E, 2)

      I = COUNT(A%E(NI,:) .EQ. HUGE(1.D0))
      J = COUNT(A%E(:,1) .EQ. HUGE(1.D0))

      ALLOCATE( IPIV(MIN((NI-I)+(NJ-J)-1, NJ)) )

      CALL LU_RB(A, IPIV)
      IPIVC = IPIV
      AE = A%E

      BR = REAL(B)
      CALL DGBTRS('N', N, KL, KU, NRHS, AE, 2*KL+KU+1, IPIV, BR, N, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'SOLVE: LINEAR SYSTEM UNABLE TO BE SOLVED.'
        STOP
      ENDIF

      BI = AIMAG(B)
      CALL DGBTRS('N', N, KL, KU, NRHS, AE, 2*KL+KU+1, IPIVC, BI, N, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'SOLVE: LINEAR SYSTEM UNABLE TO BE SOLVED.'
        STOP
      ENDIF

      B = BR + BI * IU

      ! DO I = 1, SIZE(B)
      !   B(I) = CMPLX(DNINT(REAL(B(I))*1.D30)*1.D-30, DNINT(AIMAG(B(I))*1.D30)*1.D-30, P8)
      ! ENDDO

      DEALLOCATE( AE, IPIV, IPIVC )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE SOLVE_RB_C_MULTI(A,B)
! ==========================================================================================================
! [USAGE]:
! SOLVE THE SYSTEM(S) OF LINEAR EQUATIONS AX = B, I.E., X = A^(-1)B
! WHERE A NEEDS TO BE LU-FACTORIZED IN ADVANCE AS L+U-I
! [VARIABLES]:
! A >> N X N SQUARE COMPLEX MATRIX AFTER LU FACTORIZATION IN A BAND FORMAT.
! B >> ON ENTRY, COMPLEX N X K RHS MATRIX OF THE SYSTEM. ON EXIT, THE SOLUTION WILL BE RECORDED.
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(REAL_BNDM), INTENT(INOUT)              :: A
      COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT)  :: B

      REAL(P8), DIMENSION(:,:), ALLOCATABLE       :: AE
      REAL(P8), DIMENSION(SIZE(B, 1), SIZE(B, 2)) :: BR, BI
      INTEGER, DIMENSION(:), ALLOCATABLE          :: IPIV, IPIVC
      INTEGER                                     :: N, NRHS, KL, KU, INFO, NI, NJ, I, J

      EXTERNAL DGBTRS

      IF ( SIZE(A%E, 2) .NE. SIZE(B, 1) ) THEN
        WRITE (*,*) 'SOLVE: MATRIX SIZE INCONSISTENT'
        STOP
      ENDIF

      N = SIZE(B, 1)
      NRHS = SIZE(B, 2)
      KL = A%SUBLEN
      KU = A%SUPLEN

      NI = SIZE(A%E, 1)
      NJ = SIZE(A%E, 2)

      I = COUNT(A%E(NI,:) .EQ. HUGE(1.D0))
      J = COUNT(A%E(:,1) .EQ. HUGE(1.D0))

      ALLOCATE( IPIV(MIN((NI-I)+(NJ-J)-1, NJ)) )

      CALL LU_RB(A, IPIV)
      IPIVC = IPIV
      AE = A%E

      BR = REAL(B)
      CALL DGBTRS('N', N, KL, KU, NRHS, AE, 2*KL+KU+1, IPIV, BR, N, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'SOLVE: LINEAR SYSTEM UNABLE TO BE SOLVED.'
        STOP
      ENDIF

      BI = AIMAG(B)
      CALL DGBTRS('N', N, KL, KU, NRHS, AE, 2*KL+KU+1, IPIVC, BI, N, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'SOLVE: LINEAR SYSTEM UNABLE TO BE SOLVED.'
        STOP
      ENDIF

      B = BR + BI * IU

      ! DO I = 1, SIZE(B,1)
      !   DO J = 1, SIZE(B,2)
      !     B(I,J) = CMPLX(DNINT(REAL(B(I,J))*1.D30)*1.D-30, DNINT(AIMAG(B(I,J))*1.D30)*1.D-30, P8)
      !   ENDDO
      ! ENDDO

      DEALLOCATE( AE, IPIV, IPIVC )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE SOLVE_R_R(A,B)
! ==========================================================================================================
! [USAGE]:
! SOLVE THE SYSTEM(S) OF LINEAR EQUATIONS AX = B, I.E., X = A^(-1)B
! WHERE A NEEDS TO BE LU-FACTORIZED IN ADVANCE AS L+U-I
! [VARIABLES]:
! A >> N X N SQUARE REAL MATRIX AFTER LU FACTORIZATION. IF NOT FACTORIZED, IT WILL BE FACTORIZED ON EXIT
! B >> ON ENTRY, REAL N X 1 RHS VECTOR OF THE SYSTEM. ON EXIT, THE SOLUTION WILL BE RECORDED.
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:,:), INTENT(INOUT)         :: A
      REAL(P8), DIMENSION(:), INTENT(INOUT)           :: B

      INTEGER, DIMENSION(MIN(SIZE(A, 1),SIZE(A, 2)))  :: IPIV
      INTEGER                                         :: NI, NJ, INFO

      EXTERNAL DGETRS

      IF (SIZE(A, 1) .NE. SIZE(A, 2)) THEN
        WRITE (*,*) 'SOLVE: MATRIX NOT SQUARE'
        WRITE (*,*) SIZE(A,1),'X',SIZE(A,2)
        STOP
      ENDIF

      IF ( SIZE(A, 2) .NE. SIZE(B) ) THEN
        WRITE (*,*) 'SOLVE: MATRIX SIZE INCONSISTENT'
        STOP
      ENDIF

      NI = SIZE(A, 1)
      NJ = 1

      CALL LU_R(A, IPIV)
      CALL DGETRS('N', NI, NJ, A, NI, IPIV, B, NI, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'SOLVE: LINEAR SYSTEM UNABLE TO BE SOLVED.'
        STOP
      ENDIF

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE SOLVE_R_R_MULTI(A,B)
! ==========================================================================================================
! [USAGE]:
! SOLVE THE SYSTEM(S) OF LINEAR EQUATIONS AX = B, I.E., X = A^(-1)B
! WHERE A NEEDS TO BE LU-FACTORIZED IN ADVANCE AS L+U-I
! [VARIABLES]:
! A >> N X N SQUARE REAL MATRIX AFTER LU FACTORIZATION. IF NOT FACTORIZED, IT WILL BE FACTORIZED ON EXIT
! B >> ON ENTRY, REAL N X K RHS MATRIX OF THE SYSTEM. ON EXIT, THE SOLUTION WILL BE RECORDED.
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:,:), INTENT(INOUT)        :: A
      REAL(P8), DIMENSION(:,:), INTENT(INOUT)        :: B

      INTEGER, DIMENSION(MIN(SIZE(A, 1),SIZE(A, 2))) :: IPIV
      INTEGER                                        :: NI, NJ, INFO

      EXTERNAL DGETRS

      IF (SIZE(A, 1) .NE. SIZE(A, 2)) THEN
        WRITE (*,*) 'SOLVE: MATRIX NOT SQUARE'
        WRITE (*,*) SIZE(A,1),'X',SIZE(A,2)
        STOP
      ENDIF

      IF ( SIZE(A, 2) .NE. SIZE(B, 1) ) THEN
        WRITE (*,*) 'SOLVE: MATRIX SIZE INCONSISTENT'
        STOP
      ENDIF

      NI = SIZE(A, 1)
      NJ = SIZE(B, 2)

      CALL LU_R(A, IPIV)
      CALL DGETRS('N', NI, NJ, A, NI, IPIV, B, NI, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'SOLVE: LINEAR SYSTEM UNABLE TO BE SOLVED.'
        STOP
      ENDIF

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE SOLVE_RB_R(A,B)
! ==========================================================================================================
! [USAGE]:
! SOLVE THE SYSTEM(S) OF LINEAR EQUATIONS AX = B, I.E., X = A^(-1)B
! WHERE A NEEDS TO BE LU-FACTORIZED IN ADVANCE AS L+U-I
! [VARIABLES]:
! A >> N X N SQUARE REAL MATRIX AFTER LU FACTORIZATION IN A BAND FORMAT.
! B >> ON ENTRY, REAL N X 1 RHS VECTOR OF THE SYSTEM. ON EXIT, THE SOLUTION WILL BE RECORDED.
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(REAL_BNDM), INTENT(INOUT)        :: A
      REAL(P8), DIMENSION(:), INTENT(INOUT) :: B

      INTEGER, DIMENSION(:), ALLOCATABLE    :: IPIV
      INTEGER                               :: N, NRHS, KL, KU, INFO, NI, NJ, I, J

      EXTERNAL DGBTRS

      IF ( SIZE(A%E, 2) .NE. SIZE(B) ) THEN
        WRITE (*,*) 'SOLVE: MATRIX SIZE INCONSISTENT'
        STOP
      ENDIF

      N = SIZE(B)
      NRHS = 1
      KL = A%SUBLEN
      KU = A%SUPLEN

      NI = SIZE(A%E, 1)
      NJ = SIZE(A%E, 2)

      I = COUNT(A%E(NI,:) .EQ. HUGE(1.D0))
      J = COUNT(A%E(:,1) .EQ. HUGE(1.D0))

      ALLOCATE( IPIV(MIN((NI-I)+(NJ-J)-1, NJ)) )

      CALL LU_RB(A, IPIV)
      CALL DGBTRS('N', N, KL, KU, NRHS, A%E, 2*KL+KU+1, IPIV, B, N, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'SOLVE: LINEAR SYSTEM UNABLE TO BE SOLVED.'
        STOP
      ENDIF

      ! DO I = 1, SIZE(B)
      !   B(I) = DNINT(B(I)*1.D30)*1.D-30
      ! ENDDO

      DEALLOCATE( IPIV )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE SOLVE_RB_R_MULTI(A,B)
! ==========================================================================================================
! [USAGE]:
! SOLVE THE SYSTEM(S) OF LINEAR EQUATIONS AX = B, I.E., X = A^(-1)B
! WHERE A NEEDS TO BE LU-FACTORIZED IN ADVANCE AS L+U-I
! [VARIABLES]:
! A >> N X N SQUARE REAL MATRIX AFTER LU FACTORIZATION IN A BAND FORMAT.
! B >> ON ENTRY, REAL N X K RHS MATRIX OF THE SYSTEM. ON EXIT, THE SOLUTION WILL BE RECORDED.
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(REAL_BNDM), INTENT(INOUT)          :: A
      REAL(P8), DIMENSION(:,:), INTENT(INOUT) :: B

      INTEGER, DIMENSION(:), ALLOCATABLE      :: IPIV
      INTEGER                                 :: N, NRHS, KL, KU, INFO, NI, NJ, I, J

      EXTERNAL DGBTRS

      IF ( SIZE(A%E, 2) .NE. SIZE(B, 1) ) THEN
        WRITE (*,*) 'SOLVE: MATRIX SIZE INCONSISTENT'
        STOP
      ENDIF

      N = SIZE(B, 1)
      NRHS = SIZE(B, 2)
      KL = A%SUBLEN
      KU = A%SUPLEN

      NI = SIZE(A%E, 1)
      NJ = SIZE(A%E, 2)

      I = COUNT(A%E(NI,:) .EQ. HUGE(1.D0))
      J = COUNT(A%E(:,1) .EQ. HUGE(1.D0))

      ALLOCATE( IPIV(MIN((NI-I)+(NJ-J)-1, NJ)) )
      
      CALL LU_RB(A, IPIV)
      CALL DGBTRS('N', N, KL, KU, NRHS, A%E, 2*KL+KU+1, IPIV, B, N, INFO)

      IF (INFO .NE. 0) THEN
        WRITE(*,*) 'SOLVE: LINEAR SYSTEM UNABLE TO BE SOLVED.'
        STOP
      ENDIF

      ! DO I = 1, SIZE(B,1)
      !   DO J = 1, SIZE(B,2)
      !     B(I,J) = DNINT(B(I,J)*1.D30)*1.D-30
      !   ENDDO
      ! ENDDO

      DEALLOCATE( IPIV )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
!  PRIVATE PROCEDURES ======================================================================================
! ==========================================================================================================
      FUNCTION EYE(NI) RESULT (E)
! ==========================================================================================================
! [USAGE]:                                   
! CREATE AN IDENTITY MATRIX OF ORDER NI
! [VARIABLES]:                                  
! NI >> THE ORDER OF THE IDENTITY MATRIX TO BE CREATED
! E >> NI X NI IDENTITY MATRIX
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER,INTENT(IN)        :: NI
      REAL(P8),DIMENSION(NI,NI) :: E

      INTEGER                   :: I

      E = 0.D0

      DO I = 1, NI
        E(I,I) = 1.D0
      ENDDO

      RETURN
      END FUNCTION
! ==========================================================================================================
      SUBROUTINE QUICKSORT_NORECUR(ARRAY,ARRAY1,ARRAY2,ARRAY3)
! ==========================================================================================================
! [USAGE]:                                   
! EXECUTE THE QUICK SORT IN ARRAY, THEN SORT ARRAY1-3 BASED ON THE RESULTING PERMUTATION FOR ARRAY
! [VARIABLES]:                                  
! ARRAY >> REAL 1-D ARRAY WHERE THE QUICK SORT TAKES PLACE
! ARRAY1-3 >> 1-D OR 2-D COMPLEX ARRAYS TO BE PERMUTATED AS A RESULT OF ARRAY'S SORTING RESULTS.
! [NOTE]:
! REF. https://www.mjr19.org.uk/IT/sorts/sorts.f90
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8),DIMENSION(:), INTENT(INOUT)      :: ARRAY
      COMPLEX(P8),DIMENSION(:), INTENT(INOUT)   :: ARRAY1
      COMPLEX(P8),DIMENSION(:,:), INTENT(INOUT) :: ARRAY2, ARRAY3
      
      REAL(P8)                                  :: TEMP,PIVOT
      COMPLEX(P8)                               :: TEMP1
      COMPLEX(P8),DIMENSION(SIZE(ARRAY))        :: TEMP2, TEMP3
      INTEGER                                   :: I, J, LEFT, RIGHT, LOW, HIGH
      INTEGER                                   :: STACK(2, STORAGE_SIZE(I)), STACK_PTR

      LOW = 1
      HIGH = SIZE(ARRAY)
      STACK_PTR = 1

      DO
        IF (HIGH-LOW .LT. 50) THEN ! USE INSERTION SORT ON SMALL ARRAYS
          DO I = LOW+1, HIGH
            TEMP = ARRAY(I)
            TEMP1 = ARRAY1(I)
            TEMP2 = ARRAY2(:,I)
            TEMP3 = ARRAY3(:,I)
            DO J = I-1, LOW, -1
              IF (ARRAY(J) .LE. TEMP) EXIT
              ARRAY(J+1) = ARRAY(J)
              ARRAY1(J+1) = ARRAY1(J)
              ARRAY2(:,J+1) = ARRAY2(:,J)
              ARRAY3(:,J+1) = ARRAY3(:,J)
            ENDDO
            ARRAY(J+1) = TEMP
            ARRAY1(J+1) = TEMP1
            ARRAY2(:,J+1) = TEMP2
            ARRAY3(:,J+1) = TEMP3
          ENDDO
          ! NOW POP FROM STACK
          IF (STACK_PTR.EQ.1) RETURN
          STACK_PTR = STACK_PTR - 1
          LOW = STACK(1,STACK_PTR)
          HIGH = STACK(2,STACK_PTR)
          CYCLE
        ENDIF

        ! FIND MEDIAN OF THREE PIVOT
        ! AND PLACE SENTINELS AT FIRST AND LAST ELEMENTS
        TEMP = ARRAY((LOW+HIGH)/2)
        TEMP1 = ARRAY1((LOW+HIGH)/2)
        TEMP2 = ARRAY2(:,(LOW+HIGH)/2)
        TEMP3 = ARRAY3(:,(LOW+HIGH)/2)
        ARRAY((LOW+HIGH)/2) = ARRAY(LOW+1)
        ARRAY1((LOW+HIGH)/2) = ARRAY1(LOW+1)
        ARRAY2(:,(LOW+HIGH)/2) = ARRAY2(:,LOW+1)
        ARRAY3(:,(LOW+HIGH)/2) = ARRAY3(:,LOW+1)
        IF (TEMP .GT. ARRAY(HIGH)) THEN
          ARRAY(LOW+1) = ARRAY(HIGH)
          ARRAY1(LOW+1) = ARRAY1(HIGH)
          ARRAY2(:,LOW+1) = ARRAY2(:,HIGH)
          ARRAY3(:,LOW+1) = ARRAY3(:,HIGH)
          ARRAY(HIGH) = TEMP
          ARRAY1(HIGH) = TEMP1
          ARRAY2(:,HIGH) = TEMP2
          ARRAY3(:,HIGH) = TEMP3
        ELSE
          ARRAY(LOW+1) = TEMP
          ARRAY1(LOW+1) = TEMP1
          ARRAY2(:,LOW+1) = TEMP2
          ARRAY3(:,LOW+1) = TEMP3
        ENDIF
        IF (ARRAY(LOW).GT.ARRAY(HIGH)) THEN
          TEMP = ARRAY(LOW)
          TEMP1 = ARRAY1(LOW)
          TEMP2 = ARRAY2(:,LOW)
          TEMP3 = ARRAY3(:,LOW)
          ARRAY(LOW) = ARRAY(HIGH)
          ARRAY1(LOW) = ARRAY1(HIGH)
          ARRAY2(:,LOW) = ARRAY2(:,HIGH)
          ARRAY3(:,LOW) = ARRAY3(:,HIGH)
          ARRAY(HIGH) = TEMP
          ARRAY1(HIGH) = TEMP1
          ARRAY2(:,HIGH) = TEMP2
          ARRAY3(:,HIGH) = TEMP3
        ENDIF
        IF (ARRAY(LOW) .GT. ARRAY(LOW+1)) THEN
          TEMP = ARRAY(LOW)
          TEMP1 = ARRAY1(LOW)
          TEMP2 = ARRAY2(:,LOW)
          TEMP3 = ARRAY3(:,LOW)
          ARRAY(LOW) = ARRAY(LOW+1)
          ARRAY1(LOW) = ARRAY1(LOW+1)
          ARRAY2(:,LOW) = ARRAY2(:,LOW+1)
          ARRAY3(:,LOW) = ARRAY3(:,LOW+1)
          ARRAY(LOW+1) = TEMP
          ARRAY1(LOW+1) = TEMP1
          ARRAY2(:,LOW+1) = TEMP2
          ARRAY3(:,LOW+1) = TEMP3
        ENDIF
        PIVOT = ARRAY(LOW+1)

        LEFT = LOW + 2
        RIGHT = HIGH - 1
        DO
          DO WHILE(ARRAY(LEFT) .LT. PIVOT)
            LEFT = LEFT + 1
          ENDDO
          DO WHILE(ARRAY(RIGHT) .GT. PIVOT)
            RIGHT = RIGHT - 1
          ENDDO
          IF (LEFT .GE. RIGHT) EXIT
          TEMP = ARRAY(LEFT)
          TEMP1 = ARRAY1(LEFT)
          TEMP2 = ARRAY2(:,LEFT)
          TEMP3 = ARRAY3(:,LEFT)
          ARRAY(LEFT) = ARRAY(RIGHT)
          ARRAY1(LEFT) = ARRAY1(RIGHT)
          ARRAY2(:,LEFT) = ARRAY2(:,RIGHT)
          ARRAY3(:,LEFT) = ARRAY3(:,RIGHT)
          ARRAY(RIGHT) = TEMP
          ARRAY1(RIGHT) = TEMP1
          ARRAY2(:,RIGHT) = TEMP2
          ARRAY3(:,RIGHT) = TEMP3
          LEFT = LEFT + 1
          RIGHT = RIGHT - 1
        ENDDO
        IF (LEFT .EQ. RIGHT) LEFT = LEFT + 1
        IF (LEFT .LT. (LOW+HIGH)/2) THEN
          STACK(1,STACK_PTR) = LEFT
          STACK(2,STACK_PTR) = HIGH
          STACK_PTR = STACK_PTR+1
          HIGH = LEFT-1
        ELSE
          STACK(1,STACK_PTR) = LOW
          STACK(2,STACK_PTR) = LEFT-1
          STACK_PTR = STACK_PTR+1
          LOW = LEFT
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE INSERTION_SORT(ARRAY,ARRAY1,ARRAY2,ARRAY3)
! ==========================================================================================================
! [USAGE]:                                   
! EXECUTE THE INSERTION SORT IN ARRAY, THEN SORT ARRAY1-3 BASED ON THE RESULTING PERMUTATION FOR ARRAY
! [VARIABLES]:                                  
! ARRAY >> REAL 1-D ARRAY WHERE THE INSERTION SORT TAKES PLACE
! ARRAY1-3 >> 1-D OR 2-D COMPLEX ARRAYS TO BE PERMUTATED AS A RESULT OF ARRAY'S SORTING RESULTS.
! [NOTE]:
! REF. https://www.mjr19.org.uk/IT/sorts/sorts.f90
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8),DIMENSION(:), INTENT(INOUT)      ::ARRAY
      COMPLEX(P8),DIMENSION(:), INTENT(INOUT)   :: ARRAY1
      COMPLEX(P8),DIMENSION(:,:), INTENT(INOUT) :: ARRAY2,ARRAY3
      
      INTEGER                                   :: I,J
      REAL(P8)                                  :: TEMP
      COMPLEX(P8)                               :: TEMP1
      COMPLEX(P8),DIMENSION(SIZE(ARRAY))        :: TEMP2, TEMP3

      DO I = 2, SIZE(ARRAY)
        TEMP = ARRAY(I)
        TEMP1 = ARRAY1(I)
        TEMP2 = ARRAY2(:,I)
        TEMP3 = ARRAY3(:,I)
        DO J = I-1, 1, -1
          IF (ARRAY(J) .LE. TEMP) EXIT
          ARRAY(J+1) = ARRAY(J)
          ARRAY1(J+1) = ARRAY1(J)
          ARRAY2(:,J+1) = ARRAY2(:,J)
          ARRAY3(:,J+1) = ARRAY3(:,J)
        ENDDO
        ARRAY(J+1) = TEMP
        ARRAY1(J+1) = TEMP1
        ARRAY2(:,J+1) = TEMP2
        ARRAY3(:,J+1) = TEMP3
      ENDDO
      
      RETURN
      END SUBROUTINE
! ==========================================================================================================
END MODULE