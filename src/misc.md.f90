MODULE MISC ! MODULE FOR MISCCELANEOUS SUBROUTINES & FUNCTIONS
      USE OMP_LIB
      IMPLICIT NONE
      PRIVATE
! ==========================================================================================================
!  PARAMS/TYPES ============================================================================================
! ==========================================================================================================
      LOGICAL, PARAMETER, PUBLIC     :: WARN_ON = .FALSE. ! WARNING INDICATOR ON/OFF (DEBUGGING PURPOSES)
      INTEGER, PARAMETER, PUBLIC     :: P4 = SELECTED_REAL_KIND(P=6, R=37)   ! S.PRECISION DATATYPE (32-BIT)
      INTEGER, PARAMETER, PUBLIC     :: P8 = SELECTED_REAL_KIND(P=15, R=307) ! D.PRECISION DATATYPE (64-BIT)
      REAL(P8), PARAMETER, PUBLIC    :: PI = ACOS(-1.D0)                     ! 3.141592653589793 ...
      COMPLEX(P8), PARAMETER, PUBLIC :: IU = (0.D0, 1.D0)                    ! IMAGINARY UNIT
      TYPE, PUBLIC :: INPUT_PARAMS                                           ! INPUT.PARAMS
        INTEGER           :: NR, NP, NZ, VISCPOW 
        REAL(P8)          :: ELL, ZLEN, VISC, VISCP

        INTEGER           :: NI, TOTN
        REAL(P8)          :: DT, TI, TOTT

        INTEGER           :: NPARTS
        REAL(P8)          :: DIAM, DENS

        CHARACTER(LEN=1)  :: PDISTTYPE
        REAL(P8)          :: XPC, YPC, ZPC, XSCAL, YSCAL, ZSCAL
        REAL(P8)          :: VXPC, VYPC, VZPC

        REAL(P8)          :: Q, XQO, YQO

        INTEGER           :: PFFITER
        REAL(P8)          :: PFFEPS

        LOGICAL           :: ISFFRMV
        INTEGER           :: FFRMVINTVL

        CHARACTER(LEN=72) :: DATDIR
        LOGICAL           :: ISDATSAV
        INTEGER           :: DATSAVINTVL

        CHARACTER(LEN=72) :: LOGDIR
        LOGICAL           :: ISPTCLTRCK
        INTEGER           :: PTCLTRCKN
      END TYPE
! ==========================================================================================================
!  DECLARATIONS ============================================================================================
! ==========================================================================================================
      PUBLIC :: PRINT_REAL_TIME                                                    ! PRINT A REAL-TIME CLOCK
      PUBLIC :: MCAT                                                                   ! ENTRIES OF A MATRIX
      PUBLIC :: MSAVE                                                            ! SAVE A MATRIX INTO A FILE
      PUBLIC :: MLOAD                                                            ! LOAD A MATRIX FROM A FILE
      PUBLIC :: ATON                                                                ! ASCII STRING TO NUMBER
      PUBLIC :: NTOA                                                             ! NUMBER TO FIXED-LEN ASCII
      PUBLIC :: READ_INPUTPARAMS                                  ! READ THE INPUT VARIABLES IN input.params
! ==========================================================================================================
!  INTERFACES ==============================================================================================
! ==========================================================================================================
      INTERFACE PRINT_REAL_TIME
        MODULE PROCEDURE PRINT_REAL_TIME
      END INTERFACE

      INTERFACE MCAT
        MODULE PROCEDURE MCAT_2DC
        MODULE PROCEDURE MCAT_2DR
        MODULE PROCEDURE MCAT_1DC
        MODULE PROCEDURE MCAT_1DR
      END INTERFACE

      INTERFACE MSAVE
        MODULE PROCEDURE MSAVE_3DC 
        MODULE PROCEDURE MSAVE_3DR
        MODULE PROCEDURE MSAVE_2DC
        MODULE PROCEDURE MSAVE_2DR
        MODULE PROCEDURE MSAVE_1DC
        MODULE PROCEDURE MSAVE_1DR
      END INTERFACE

      INTERFACE MLOAD
        MODULE PROCEDURE MLOAD_3DC
        MODULE PROCEDURE MLOAD_3DR
        MODULE PROCEDURE MLOAD_2DC
        MODULE PROCEDURE MLOAD_2DR
        MODULE PROCEDURE MLOAD_1DC
        MODULE PROCEDURE MLOAD_1DR
      END INTERFACE

      INTERFACE ATON
        MODULE PROCEDURE ATOF
        MODULE PROCEDURE ATOI
      END INTERFACE

      INTERFACE NTOA
        MODULE PROCEDURE ITOA
        MODULE PROCEDURE FTOA
      END INTERFACE

      INTERFACE READ_INPUTPARAMS
        MODULE PROCEDURE READ_INPUTPARAMS
      END INTERFACE

CONTAINS
! ==========================================================================================================
!  PUBLIC PROCEDURES =======================================================================================
! ==========================================================================================================
      SUBROUTINE PRINT_REAL_TIME()
! ==========================================================================================================
! [USAGE]: 
! PRINT THE CURRENT DATE AND TIME ON CLI
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      CHARACTER(8)  ::  DATE
      CHARACTER(10) ::  TIME
      CHARACTER(5)  ::  ZONE
      INTEGER       ::  VALS(1:8)

      CALL DATE_AND_TIME(DATE,TIME,ZONE,VALS)                                ! INTRINSIC SUBROUTINE

      WRITE(*,101) VALS(1),VALS(2),VALS(3),VALS(5),VALS(6),VALS(7)
 101  FORMAT(' @ 'I0.4,'-',I0.2,'-',I0.2,' ',I0.2,':',I0.2,':',I0.2)
      WRITE(*,*) ''

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE MCAT_2DC(A)
! ==========================================================================================================
! [USAGE]: 
! PRINT A GIVEN 2D COMPLEX MATRIX ON CLI
! [VARIABLES]:
! A >> NI X NJ COMPLEX MATRIX TO BE PRINTED ON CLI
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER                                 :: NI, NJ, PRINT_ROW
      INTEGER                                 :: I, J, JB, JE, JJ
      COMPLEX(P8), DIMENSION(:,:), INTENT(IN) :: A

      PRINT_ROW = 4

      NI = SIZE(A, 1)
      NJ = SIZE(A, 2)

      DO J = 1, (NJ+PRINT_ROW-1)/PRINT_ROW
        JB = PRINT_ROW*J-(PRINT_ROW-1)
        JE = MIN(PRINT_ROW*J, NJ)
        WRITE(*,101) (JJ, JJ = JB, JE)
 101    FORMAT(5X,8I24)
        DO I = 1, NI
          WRITE(*, 102) I, (A(I,JJ), JJ = JB, JE)
 102      FORMAT(S,I4.4, ':', 8(S,ES12.3E3,SP,ES11.3E3,'i'))
        ENDDO
      ENDDO
      WRITE(*,*) ''

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE MCAT_2DR(A)
! ==========================================================================================================
! [USAGE]: 
! PRINT A GIVEN 2D REAL MATRIX ON CLI
! [VARIABLES]:
! A >> NI X NJ REAL MATRIX TO BE PRINTED ON CLI
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER                              :: NI, NJ, PRINT_ROW
      INTEGER                              :: I, J, JB, JE, JJ
      REAL(P8), DIMENSION(:,:), INTENT(IN) :: A

      PRINT_ROW = 8

      NI = SIZE(A, 1)
      NJ = SIZE(A, 2)

      DO J = 1, (NJ+PRINT_ROW-1)/PRINT_ROW
        JB = PRINT_ROW*J-(PRINT_ROW-1)
        JE = MIN(PRINT_ROW*J, NJ)
        WRITE(*,101) (JJ, JJ = JB, JE)
 101    FORMAT(5X,8I12)
        DO I = 1, NI
          WRITE(*, 102) I, (A(I,JJ), JJ = JB, JE)
 102      FORMAT(I4.4, ':', 8ES12.3E3)
        ENDDO
      ENDDO
      WRITE(*,*) ''

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE MCAT_1DC(A)
! ==========================================================================================================
! [USAGE]: 
! PRINT A GIVEN 1D COMPLEX MATRIX (VECTOR) ON CLI. MORE PRECISION THAN THE COUNTERPART FOR 2D MATRICES
! [VARIABLES]:
! A >> N X 1 COMPLEX MATRIX (VECTOR) TO BE PRINTED ON CLI
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER                               :: I
      COMPLEX(P8), DIMENSION(:), INTENT(IN) :: A

      WRITE(*,101) (I, A(I), I = 1, SIZE(A))
 101  FORMAT(S,I4.4, ':', S,ES24.15E3, SP,ES23.15E3,'i')
      WRITE(*,*) ''

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE MCAT_1DR(A)
! ==========================================================================================================
! [USAGE]: 
! PRINT A GIVEN 1D REAL MATRIX (VECTOR) ON CLI. MORE PRECISION THAN THE COUNTERPART FOR 2D MATRICES
! [VARIABLES]:
! A >> N X 1 REAL MATRIX (VECTOR) TO BE PRINTED ON CLI
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER                            :: I
      REAL(P8), DIMENSION(:), INTENT(IN) :: A

      WRITE(*,101) (I, A(I), I = 1, SIZE(A))
 101  FORMAT(I4, ':', ES24.15E3)
      WRITE(*,*) ''

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE MSAVE_3DC(A, FN)
! ==========================================================================================================
! [USAGE]: 
! SAVE A GIVEN 3D COMPLEX MATRIX (TENSOR) INTO FILE LOCATED AT FN
! [VARIABLES]:
! A >> NI X NJ X NK COMPLEX MATRIX TO BE SAVED
! FN >> FILENAME (PATH) TO SAVE THE INPUT MATRIX
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER                                   :: NI, NJ, NK
      INTEGER                                   :: I, J, K, FU
      COMPLEX(P8), DIMENSION(:,:,:), INTENT(IN) :: A
      CHARACTER(LEN=*), INTENT(IN)              :: FN

      NI = SIZE(A, 1)
      NJ = SIZE(A, 2)
      NK = SIZE(A, 3)

      CALL LABOPEN_(FN, NJ*2, FU) ! PRIVATE SUBROUTINE
      WRITE(FU, 101) NI, NJ, NK
 101  FORMAT('/* C ', I4.4 ' ', I4.4' ', I4.4, ' */')

      DO K = 1, NK
        DO I = 1, NI
          DO J = 1, NJ
            CALL LABPR_(REAL(A(I,J,K)), FU) ! PRIVATE SUBROUTINE
            CALL LABPR_(AIMAG(A(I,J,K)), FU) ! PRIVATE SUBROUTINE
          ENDDO
          CALL LABCR_(FU) ! PRIVATE SUBROUTINE
        ENDDO
      ENDDO

      CALL LABCLOSE_(FU) ! PRIVATE SUBROUTINE

      WRITE(*,102) NI, NJ, NK, TRIM(ADJUSTL(FN))
 102  FORMAT('CMPX MATRIX OF DIM. ',I4.4,'X',I4.4,'X',I4.4,' WRITTEN TO ',A)
      WRITE(*,*) ''

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE MSAVE_3DR(A, FN)
! ==========================================================================================================
! [USAGE]: 
! SAVE A GIVEN 3D REAL MATRIX (TENSOR) INTO FILE LOCATED AT FN
! [VARIABLES]:
! A >> NI X NJ X NK REAL MATRIX TO BE SAVED
! FN >> FILENAME (PATH) TO SAVE THE INPUT MATRIX
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER                                :: NI, NJ, NK
      INTEGER                                :: I, J, K, FU
      REAL(P8), DIMENSION(:,:,:), INTENT(IN) :: A
      CHARACTER(LEN=*), INTENT(IN)           :: FN

      NI = SIZE(A, 1)
      NJ = SIZE(A, 2)
      NK = SIZE(A, 3)

      CALL LABOPEN_(FN, NJ, FU) ! PRIVATE SUBROUTINE
      WRITE(FU, 101) NI, NJ, NK
 101  FORMAT('/* R ', I4.4 ' ', I4.4' ', I4.4, ' */')

      DO K = 1, NK
        DO I = 1, NI
          DO J = 1, NJ
            CALL LABPR_(A(I,J,K), FU) ! PRIVATE SUBROUTINE
          ENDDO
          CALL LABCR_(FU) ! PRIVATE SUBROUTINE
        ENDDO
      ENDDO

      CALL LABCLOSE_(FU) ! PRIVATE SUBROUTINE

      WRITE(*,102) NI, NJ, NK, TRIM(ADJUSTL(FN))
 102  FORMAT('REAL MATRIX OF DIM. ',I4.4,'X',I4.4,'X',I4.4,' WRITTEN TO ',A)
      WRITE(*,*) ''

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE MSAVE_2DC(A, FN)
! ==========================================================================================================
! [USAGE]: 
! SAVE A GIVEN 2D COMPLEX MATRIX INTO FILE LOCATED AT FN
! [VARIABLES]:
! A >> NI X NJ COMPLEX MATRIX TO BE SAVED
! FN >> FILENAME (PATH) TO SAVE THE INPUT MATRIX
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:,:), INTENT(IN) :: A
      CHARACTER(LEN=*), INTENT(IN)            :: FN

      CALL MSAVE_3DC( RESHAPE(A, (/ SIZE(A,1 ), SIZE(A, 2), 1 /)), FN)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE MSAVE_2DR(A, FN)
! ==========================================================================================================
! [USAGE]: 
! SAVE A GIVEN 2D REAL MATRIX INTO FILE LOCATED AT FN
! [VARIABLES]:
! A >> NI X NJ REAL MATRIX TO BE SAVED
! FN >> FILENAME (PATH) TO SAVE THE INPUT MATRIX
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:,:), INTENT(IN) :: A
      CHARACTER(LEN=*), INTENT(IN)         :: FN

      CALL MSAVE_3DR( RESHAPE(A, (/ SIZE(A,1 ), SIZE(A, 2), 1 /)), FN)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE MSAVE_1DC(A, FN)
! ==========================================================================================================
! [USAGE]: 
! SAVE A GIVEN 1D COMPLEX MATRIX (VECTOR) INTO FILE LOCATED AT FN
! [VARIABLES]:
! A >> NI X 1 COMPLEX MATRIX (VECTOR) TO BE SAVED
! FN >> FILENAME (PATH) TO SAVE THE INPUT MATRIX
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      COMPLEX(P8), DIMENSION(:), INTENT(IN) :: A
      CHARACTER(LEN=*), INTENT(IN)          :: FN

      CALL MSAVE_3DC( RESHAPE(A, (/ SIZE(A, 1), 1, 1 /)), FN)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE MSAVE_1DR(A, FN)
! ==========================================================================================================
! [USAGE]: 
! SAVE A GIVEN 1D REAL MATRIX (VECTOR) INTO FILE LOCATED AT FN
! [VARIABLES]:
! A >> NI X 1 REAL MATRIX TO BE SAVED
! FN >> FILENAME (PATH) TO SAVE THE INPUT MATRIX
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8), DIMENSION(:), INTENT(IN) :: A
      CHARACTER(LEN=*), INTENT(IN)       :: FN

      CALL MSAVE_3DR( RESHAPE(A, (/ SIZE(A, 1), 1, 1 /)), FN)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE MLOAD_3DC(FN, A)
! ==========================================================================================================
! [USAGE]: 
! LOAD A 3D COMPLEX MATRIX (TENSOR) FROM A MATRIX FILE LOCATED AT FN
! [VARIABLES]:
! FN >> FILENAME (PATH) WHERE THE MATRIX TO BE LOADED IS STORED
! A >> NI X NJ X NK COMPLEX MATRIX TO BE LOADED
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER                                      :: NI, NJ, NK
      INTEGER                                      :: I, J, K, FU
      REAL(P8), DIMENSION(:,:,:), ALLOCATABLE      :: AL
      COMPLEX(P8), DIMENSION(:,:,:), INTENT(INOUT) :: A
      CHARACTER(LEN=*), INTENT(IN)                 :: FN
      CHARACTER(LEN=22)                            :: HD 

      NI = SIZE(A, 1)
      NJ = SIZE(A, 2)
      NK = SIZE(A, 3)

      ALLOCATE( AL(NI, NJ*2, NK) )
      AL = 0.D0

      CALL LABOPEN_(FN, NJ*2, FU)
      READ(FU, 101) HD
 101  FORMAT(A22)

      IF (HD(4:4) .EQ. 'R') THEN
        WRITE(*,*) 'MLOAD: FILE (REAL) INCOMPATIBLE WITH MATRIX TYPE (COMPLEX).'
        STOP
      ELSEIF (HD(4:4) .NE. 'C') THEN
        WRITE(*,*) 'MLOAD: HEADER CORRUPTED. MATRIX TYPE IS NEITHER `R` NOR `C`.'
        STOP
      ENDIF 

      CALL ATON(HD( 6: 9),I)
      CALL ATON(HD(11:14),J)
      CALL ATON(HD(16:19),K)

      IF (NI .NE. I) THEN
        WRITE(*,*) 'MLOAD: [WARNING] I-DIMENSION INCOMPATIBLE. MATRIX LOADING MAY BE INCOMPLETE.'
        WRITE(*,*) ''
      ENDIF

      IF (NJ .NE. J) THEN
        WRITE(*,*) 'MLOAD: [WARNING] J-DIMENSION INCOMPATIBLE. MATRIX LOADING MAY BE INCOMPLETE.'
        WRITE(*,*) ''
      ENDIF

      IF (NK .NE. K) THEN
        WRITE(*,*) 'MLOAD: [WARNING] K-DIMENSION INCOMPATIBLE. MATRIX LOADING MAY BE INCOMPLETE.'
        WRITE(*,*) ''
      ENDIF

      DO K = 1, NK
        DO I = 1, NI
          CALL LABRD_(AL(I,:,K), FU)
        ENDDO
      ENDDO

      A = AL(:,1:2*NJ-1:2,:) + AL(:,2:NJ*2:2,:) * IU

      DEALLOCATE( AL )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE MLOAD_3DR(FN, A)
! ==========================================================================================================
! [USAGE]: 
! LOAD A 3D REAL MATRIX (TENSOR) FROM A MATRIX FILE LOCATED AT FN
! [VARIABLES]:
! FN >> FILENAME (PATH) WHERE THE MATRIX TO BE LOADED IS STORED
! A >> NI X NJ X NK REAL MATRIX TO BE LOADED
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER                                   :: NI, NJ, NK
      INTEGER                                   :: I, J, K, FU
      REAL(P8), DIMENSION(:,:,:), ALLOCATABLE   :: AL
      REAL(P8), DIMENSION(:,:,:), INTENT(INOUT) :: A
      CHARACTER(LEN=*), INTENT(IN)              :: FN
      CHARACTER(LEN=22)                         :: HD

      NI = SIZE(A, 1)
      NJ = SIZE(A, 2)
      NK = SIZE(A, 3)

      ALLOCATE( AL(NI, NJ, NK) )
      AL = 0.D0

      CALL LABOPEN_(FN, NJ, FU)
      READ(FU, 101) HD
 101  FORMAT(A22)

      IF (HD(4:4) .EQ. 'C') THEN
        WRITE(*,*) 'MLOAD: FILE (COMPLEX) INCOMPATIBLE WITH MATRIX TYPE (REAL).'
        STOP
      ELSEIF (HD(4:4) .NE. 'R') THEN
        WRITE(*,*) 'MLOAD: HEADER CORRUPTED. MATRIX TYPE IS NEITHER `R` NOR `C`.'
        STOP
      ENDIF 

      CALL ATON(HD( 6: 9),I)
      CALL ATON(HD(11:14),J)
      CALL ATON(HD(16:19),K)

      IF (NI .NE. I) THEN
        WRITE(*,*) 'MLOAD: [WARNING] I-DIMENSION INCOMPATIBLE. MATRIX LOADING MAY BE INCOMPLETE.'
        WRITE(*,*) ''
      ENDIF

      IF (NJ .NE. J) THEN
        WRITE(*,*) 'MLOAD: [WARNING] J-DIMENSION INCOMPATIBLE. MATRIX LOADING MAY BE INCOMPLETE.'
        WRITE(*,*) ''
      ENDIF

      IF (NK .NE. K) THEN
        WRITE(*,*) 'MLOAD: [WARNING] K-DIMENSION INCOMPATIBLE. MATRIX LOADING MAY BE INCOMPLETE.'
        WRITE(*,*) ''
      ENDIF

      DO K = 1, NK
        DO I = 1, NI
          CALL LABRD_(AL(I,:,K), FU)
        ENDDO
      ENDDO

      A = AL(:,:,:)

      DEALLOCATE( AL )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE MLOAD_2DC(FN, A)
! ==========================================================================================================
! [USAGE]: 
! LOAD A 2D COMPLEX MATRIX FROM A MATRIX FILE LOCATED AT FN
! [VARIABLES]:
! FN >> FILENAME (PATH) WHERE THE MATRIX TO BE LOADED IS STORED
! A >> NI X NJ COMPLEX MATRIX TO BE LOADED
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER                                    :: NI, NJ
      COMPLEX(P8), DIMENSION(:,:,:), ALLOCATABLE :: AL
      COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT) :: A
      CHARACTER(LEN=*), INTENT(IN)               :: FN

      NI = SIZE(A, 1)
      NJ = SIZE(A, 2)

      ALLOCATE( AL(NI,NJ,1) )

      CALL MLOAD_3DC(FN, AL)

      A = RESHAPE(AL, (/ NI, NJ /))

      DEALLOCATE( AL )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE MLOAD_2DR(FN, A)
! ==========================================================================================================
! [USAGE]: 
! LOAD A 2D REAL MATRIX FROM A MATRIX FILE LOCATED AT FN
! [VARIABLES]:
! FN >> FILENAME (PATH) WHERE THE MATRIX TO BE LOADED IS STORED
! A >> NI X NJ REAL MATRIX TO BE LOADED
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER                                 :: NI, NJ
      REAL(P8), DIMENSION(:,:,:), ALLOCATABLE :: AL
      REAL(P8), DIMENSION(:,:), INTENT(INOUT) :: A
      CHARACTER(LEN=*), INTENT(IN)            :: FN

      NI = SIZE(A, 1)
      NJ = SIZE(A, 2)

      ALLOCATE( AL(NI,NJ,1) )

      CALL MLOAD_3DR(FN, AL)

      A = RESHAPE(AL, (/ NI, NJ /))

      DEALLOCATE( AL )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE MLOAD_1DC(FN, A)
! ==========================================================================================================
! [USAGE]: 
! LOAD A 1D COMPLEX MATRIX (VECTOR) FROM A MATRIX FILE LOCATED AT FN
! [VARIABLES]:
! FN >> FILENAME (PATH) WHERE THE MATRIX TO BE LOADED IS STORED
! A >> NI X 1 COMPLEX MATRIX (VECTOR) TO BE LOADED
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER                                    :: NI
      COMPLEX(P8), DIMENSION(:,:,:), ALLOCATABLE :: AL
      COMPLEX(P8), DIMENSION(:), INTENT(INOUT)   :: A
      CHARACTER(LEN=*), INTENT(IN)               :: FN

      NI = SIZE(A, 1)

      ALLOCATE( AL(NI,1,1) )

      CALL MLOAD_3DC(FN, AL)

      A = RESHAPE(AL, (/ NI /))

      DEALLOCATE( AL )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE MLOAD_1DR(FN, A)
! ==========================================================================================================
! [USAGE]: 
! LOAD A 1D REAL MATRIX (VECTOR) FROM A MATRIX FILE LOCATED AT FN
! [VARIABLES]:
! FN >> FILENAME (PATH) WHERE THE MATRIX TO BE LOADED IS STORED
! A >> NI X 1 REAL MATRIX (VECTOR) TO BE LOADED
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER                                 :: NI
      REAL(P8), DIMENSION(:,:,:), ALLOCATABLE :: AL
      REAL(P8), DIMENSION(:), INTENT(INOUT)   :: A
      CHARACTER(LEN=*), INTENT(IN)            :: FN

      NI = SIZE(A, 1)

      ALLOCATE( AL(NI,1,1) )

      CALL MLOAD_3DR(FN, AL)

      A = RESHAPE(AL, (/ NI /))

      DEALLOCATE( AL )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE ATOI(COMM, I)
! ==========================================================================================================
! [USAGE]:
! CONVERT ASCII TO INTEGER
! [VARIABLES]:
! COMM >> ASCII STRING [0-9]
! I    >> CONVERTED INTEGER
! NOTE:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE
      
      INTEGER, INTENT(OUT)         :: I
      REAL(P8)                     :: VA
      CHARACTER(LEN=*), INTENT(IN) :: COMM
      CHARACTER(LEN=256)           :: ST

      IF (LEN(COMM)>256) THEN
        WRITE(*,*) 'ATON: ASCII TOO LONG (>256).'
        STOP
      ENDIF

      ST = COMM(1:LEN(COMM)) // REPEAT(" ", 256)
      READ(UNIT=ST, FMT=*, ERR=901) VA

      I = VA

      RETURN

 901  CONTINUE
      WRITE(*,*) 'ATON: ERROR IN CONVESION. ASCII MUST REPRESENT INTEGER.'
      WRITE(*,*) 'ASCII= ', COMM
      STOP

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE ATOF(COMM, F)
! ==========================================================================================================
! [USAGE]:
! CONVERT ASCII TO FLOATING NUMBER
! [VARIABLES]:
! COMM >> ASCII STRING [0-9]
! F    >> CONVERTED FLOATING NUMBER
! NOTE:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE
      
      REAL(P8), INTENT(OUT)        :: F
      REAL(P8)                     :: VA
      CHARACTER(LEN=*), INTENT(IN) :: COMM
      CHARACTER(LEN=256)           :: ST

      IF (LEN(COMM)>256) THEN
        WRITE(*,*) 'ATON: ASCII TOO LONG (>256).'
        STOP
      ENDIF

      ST = COMM(1:LEN(COMM)) // REPEAT(" ", 256)
      READ(UNIT=ST, FMT='(F256.0)', ERR=901) VA

      F = VA

      RETURN

 901  CONTINUE
      WRITE(*,*) 'ATON: ERROR IN CONVESION. ASCII MUST REPRESENT FLOATING NUMBER.'
      WRITE(*,*) 'ASCII= ', COMM
      STOP

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      FUNCTION ITOA(I, FMT) RESULT(RES)
! ==========================================================================================================
! [USAGE]:
! CONVERT INTEGER TO STRING IN COMPLIANCE WITH THE PROVIDED FORMATTING FMT
! [VARIABLES]:
! I >> INTEGER
! FMT >> PRINT FORMAT E.G. I4.4 -> 0001, 0002, ...
! RES >> CONVERTED STRING CORRESPONDING TO ITOA(I, FMT)
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER                       :: I
      CHARACTER(LEN=*)              :: FMT
      CHARACTER(LEN=72)             :: TMP
      CHARACTER(LEN=:), ALLOCATABLE :: RES
    
      WRITE(TMP,FMT) I
      RES = TRIM(TMP)
  
      RETURN
      END FUNCTION
! ==========================================================================================================
      FUNCTION FTOA(F, FMT) RESULT(RES)
! ==========================================================================================================
! [USAGE]:
! CONVERT FLOATING NUMBER TO STRING IN COMPLIANCE WITH THE PROVIDED FORMATTING FMT
! [VARIABLES]:
! F >> FLOATING NUMBER
! FMT >> PRINT FORMAT E.G. I4.4 -> 0001, 0002, ...
! RES >> CONVERTED STRING CORRESPONDING TO ITOA(I, FMT)
! ==========================================================================================================
      IMPLICIT NONE

      REAL(P8)                      :: F
      CHARACTER(LEN=*)              :: FMT
      CHARACTER(LEN=72)             :: TMP
      CHARACTER(LEN=:), ALLOCATABLE :: RES
    
      WRITE(TMP,FMT) F
      RES = TRIM(TMP)
  
      RETURN
      END FUNCTION
! ==========================================================================================================
      SUBROUTINE READ_INPUTPARAMS(PARAMS)
! ==========================================================================================================
! [USAGE]:
! LOAD THE VARIABLES AS WRITTEN IN input.params
! [VARIABLES]:
! PARAMS >> INPUT_PARAMS-TYPE VARIABLE
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(INPUT_PARAMS), INTENT(INOUT) :: PARAMS

      INTEGER                       :: FU
      CHARACTER(LEN=1)              :: DUM


      OPEN(UNIT=FU, FILE='./input.params', STATUS='UNKNOWN', FORM='FORMATTED')

      READ(FU,*) DUM ! COMPUTATIONAL DOMAIN INFO.
      READ(FU,*) DUM; READ(FU,*) PARAMS%NR, PARAMS%NP, PARAMS%NZ
      READ(FU,*) DUM; READ(FU,*) PARAMS%ELL, PARAMS%ZLEN
      READ(FU,*) DUM

      READ(FU,*) DUM ! TIME STEPPING INFO.
      READ(FU,*) DUM; READ(FU,*) PARAMS%DT, PARAMS%TI, PARAMS%TOTT, PARAMS%NI, PARAMS%TOTN
      READ(FU,*) DUM

      READ(FU,*) DUM ! FLOW FIELD PROPERTY INFO.
      READ(FU,*) DUM; READ(FU,*) PARAMS%VISC, PARAMS%VISCPOW, PARAMS%VISCP
      READ(FU,*) DUM

      READ(FU,*) DUM ! INITIAL Q-VORTEX INFO.
      READ(FU,*) DUM; READ(FU,*) PARAMS%Q, PARAMS%XQO, PARAMS%YQO
      READ(FU,*) DUM

      READ(FU,*) DUM ! PARTICLE PROPERTY INFO.
      READ(FU,*) DUM; READ(FU,*) PARAMS%DIAM, PARAMS%DENS, PARAMS%NPARTS
      READ(FU,*) DUM
      
      READ(FU,*) DUM ! INITIAL PARTICLE INFO.
      READ(FU,*) DUM; READ(FU,*) PARAMS%PDISTTYPE
      READ(FU,*) DUM; READ(FU,*) PARAMS%XPC, PARAMS%YPC, PARAMS%ZPC
      READ(FU,*) DUM; READ(FU,*) PARAMS%XSCAL, PARAMS%YSCAL, PARAMS%ZSCAL
      READ(FU,*) DUM; READ(FU,*) PARAMS%VXPC, PARAMS%VYPC, PARAMS%VZPC
      READ(FU,*) DUM

      READ(FU,*) DUM ! FLOW-PARTICLE COUPLING INFO.
      READ(FU,*) DUM; READ(FU,*) PARAMS%PFFITER, PARAMS%PFFEPS
      READ(FU,*) DUM; READ(FU,*) PARAMS%ISFFRMV, PARAMS%FFRMVINTVL
      READ(FU,*) DUM

      READ(FU,*) DUM ! DATA STORAGE INFO.
      READ(FU,*) DUM; READ(FU,'(A72)') PARAMS%DATDIR; PARAMS%DATDIR = TRIM(ADJUSTL(PARAMS%DATDIR))
      IF (PARAMS%DATDIR(LEN(PARAMS%DATDIR):LEN(PARAMS%DATDIR)) .NE. '/') THEN
        PARAMS%DATDIR = TRIM(ADJUSTL(PARAMS%DATDIR))//'/'
      ENDIF
      READ(FU,*) DUM; READ(FU,*) PARAMS%ISDATSAV,PARAMS%DATSAVINTVL
      READ(FU,*) DUM

      READ(FU,*) DUM ! DATA LOGGING INFO.
      READ(FU,*) DUM; READ(FU,'(A72)') PARAMS%LOGDIR; PARAMS%LOGDIR = TRIM(ADJUSTL(PARAMS%LOGDIR))
      IF (PARAMS%LOGDIR(LEN(PARAMS%LOGDIR):LEN(PARAMS%LOGDIR)) .NE. '/') THEN
        PARAMS%LOGDIR = TRIM(ADJUSTL(PARAMS%LOGDIR))//'/'
      ENDIF
      READ(FU,*) DUM; READ(FU,*) PARAMS%ISPTCLTRCK,PARAMS%PTCLTRCKN
      READ(FU,*) DUM

      CLOSE(FU)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
!  PRIVATE PROCEDURES ======================================================================================
! ==========================================================================================================
      SUBROUTINE LABOPEN_(FN, NC, FU)
! ==========================================================================================================
! [USAGE]:
! OPEN A FILE LOCATED AT FN
! [VARIABLES]:
! FN >> FILENAME (PATH) WHERE THE VALUES TO BE LOADED ARE STORED
! NC >> NUMBER OF VALUES TO BE CALLED IN EACH LINE
! FU >> OPENED FILE'S ASSIGNED UNIT NUMBER
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER                      :: RECLEN
      INTEGER                      :: IS
      INTEGER, INTENT(IN)          :: NC, FU
      CHARACTER(LEN=*), INTENT(IN) :: FN

      RECLEN = 24*(NC+1)
      OPEN(UNIT=FU, FILE=FN, STATUS='UNKNOWN', IOSTAT=IS, RECL=RECLEN)

      IF (IS .NE. 0) THEN
        WRITE(*,*) '_LABOPEN: CANNOT OPEN', TRIM(ADJUSTL(FN))
      ENDIF

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE LABCLOSE_(FU)
! ==========================================================================================================
! [USAGE]:
! CLOSE A FILE ASSIGNED TO THE FILE UNIT NUMBER FU
! [VARIABLES]:
! FU >> OPENED FILE'S ASSIGNED UNIT NUMBER
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)      :: FU

      CLOSE(FU)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE LABPR_(A, FU)
! ==========================================================================================================
! [USAGE]:
! WRITE A SINGLE REAL NUMERIC VALUE INTO FILE ASSIGNED TO UNIT NUMBER FU
! [VARIABLES]:
! A >> A REAL NUMERIC VALUE TO BE WRITTEN
! FU >> FILE UNIT NUMBER
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)     :: FU
      REAL(P8), INTENT(IN)    :: A

      WRITE(UNIT=FU, FMT=110, ADVANCE='NO') A
110   FORMAT(ES24.15E3)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE LABCR_(FU)
! ==========================================================================================================
! [USAGE]:
! WRITE A CARRAGE RETURN INTO FILE ASSIGNED TO UNIT NUMBER FU
! [VARIABLES]:
! FU >> FILE UNIT NUMBER
! [NOTE]:
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: FU

      WRITE(UNIT=FU, FMT=111) ''
111   FORMAT(A1)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE LABRD_(A, FU)
! ==========================================================================================================
! [USAGE]:
! READ A SET OF REAL NUMERIC VALUES TO A FROM THE FILE ASSIGNED TO UNIT NUMBER FU
! [VARIABLES]:
! A >> A REAL NUMERIC VALUE WHERE THE READ VALUE TO BE STORED
! FU >> UNIT NUMBER OF THE MATRIX FILE
! ==========================================================================================================
      IMPLICIT NONE

      INTEGER, INTENT(IN)                   :: FU
      REAL(P8), DIMENSION(:), INTENT(INOUT) :: A

      READ(UNIT=FU, FMT=*) A(:)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
END MODULE