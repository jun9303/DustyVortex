! =========================================================================================================
!  WRITTEN BY SANGJOON (JOON) LEE, DEPT. OF MECHANICAL ENGINEERING, UNIVERSITY OF CALIFORNIA, BERKELEY    =
!  CONTACT: SANGJOONLEE@BERKELEY.EDU                                                                      =
!  UC BERKELEY COMPUTATIONAL FLUID DYNAMICS (CFD) LAB, HTTPS://CFD.ME.BERKELEY.EDU/                       =
!  SANGJOON LEE (C) 2023                                                                                  =
!  NONCOMMERCIAL USE. FOR RESEARCH PURPOSES ONLY                                                          =
! =========================================================================================================
      PROGRAM QVORTEXGEN
! =========================================================================================================
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! =========================================================================================================
      USE OMP_LIB; USE MPI
      USE MISC; USE MATOPS; USE TIMEINFO 
      USE FIELDGEN; USE FIELDOPS; USE FIELDSTEP

      IMPLICIT NONE

      INTEGER(P8)                   :: I, J, K
      REAL(P8)                      :: START, FINISH
      TYPE(SCALAR_FIELD)            :: PSI, CHI
      TYPE(VECTOR_FIELD)            :: U

      INTEGER                       :: NR, NP, NZ, VISCPOW 
      REAL(P8)                      :: ELL, ZLEN, VISC, VISCP

      INTEGER                       :: NI, TOTN
      REAL(P8)                      :: DT, TI, TOTT

      REAL(P8)                      :: Q, XQO, YQO

      CHARACTER(LEN=:), ALLOCATABLE :: DATDIR

      IF (MPI_RANK .EQ. 0) THEN
        WRITE(*,*) ''
        WRITE(*,*) '*** PROGRAM STARTED ***'
        CALL PRINT_REAL_TIME() ! MISC
        START = OMP_GET_WTIME() ! OMP_LIB
      ENDIF

! ....................................................................................................... !

      ! ------------------------------------------------------------------------------------------- !
        NR = 192; NP = 128; NZ = 128; ELL = 3.D0; ZLEN = 2*PI
        VISC = 1.D-5; VISCPOW = 2; VISCP = 0.D0
        DT = 1.D-3; TI = 0.D0; TOTT = 1.D0; NI = 0; TOTN = 1000
      ! ------------------------------------------------------------------------------------------- !      

      CALL FIELD_INIT(NRIN=NR, NPIN=NP, NZIN=NZ, ELLIN=ELL, ZLENIN=ZLEN, &
                      VISCIN=VISC, VISCPOW=VISCPOW, VISCPIN=VISCP)
      CALL TIME_INIT(DTIN=DT, TIIN=TI, TOTTIN=TOTT, NIIN=NI, TOTNIN=TOTN)

      IF (MPI_RANK .EQ. 0) THEN
        FINISH = OMP_GET_WTIME() ! OMP_LIB
        WRITE(*,101) 'FIELD & TIME INITALIZATION'//' '//REPEAT('.',48), FINISH-START
        START = OMP_GET_WTIME() ! OMP_LIB
      ENDIF

! ....................................................................................................... !

      ! ------------------------------------------------------------------------------------------- !      
        Q = 4.D0; XQO = 0.D0; YQO = 0.D0
      ! ------------------------------------------------------------------------------------------- !      

      CALL QVORTEX(VFLD=U, Q=Q, XO=XQO, YO=YQO)

      CALL VELPROJ(U,PSI,CHI,PSILN=-.5)
      U = PT2VEL(PSI,CHI)

      IF (MPI_RANK .EQ. 0) THEN
        FINISH = OMP_GET_WTIME() ! OMP_LIB
        WRITE(*,101) 'Q-VORTEX VELOCITY FIELD GENERATION'//' '//REPEAT('.',48), FINISH-START
        START = OMP_GET_WTIME() ! OMP_LIB
      ENDIF

! ....................................................................................................... !

      ! ------------------------------------------------------------------------------------------- !      
        DATDIR = './dat/'
      ! ------------------------------------------------------------------------------------------- !            

      IF (DATDIR(LEN(DATDIR):LEN(DATDIR)) .NE. '/') DATDIR = DATDIR // '/'
      CALL EXECUTE_COMMAND_LINE('mkdir -p '//TRIM(ADJUSTL(DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))))

      CALL MSAVE(U,   TRIM(ADJUSTL(DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/u.vfld')
      CALL MSAVE(PSI, TRIM(ADJUSTL(DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/psi.sfld')
      CALL MSAVE(CHI, TRIM(ADJUSTL(DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/chi.sfld')

      IF (MPI_RANK .EQ. 0) THEN
        FINISH = OMP_GET_WTIME() ! OMP_LIB
        WRITE(*,101) 'INITIAL FIELDS STORED'//' '//REPEAT('.',48), FINISH-START
        START = OMP_GET_WTIME() ! OMP_LIB
      ENDIF

! ....................................................................................................... !

      IF (MPI_RANK .EQ. 0) THEN
        WRITE(*,*) ''
        WRITE(*,*) '*** PROGRAM FINISHED ***'
        CALL PRINT_REAL_TIME() ! MISC
      ENDIF

 101  FORMAT(A48,1X,F13.6,' SECONDS.')

CONTAINS
! =========================================================================================================
!  PROGRAM-DEPENDENT PROCEDRUES ===========================================================================
! =========================================================================================================
      SUBROUTINE QVORTEX(VFLD, Q, XO, YO)
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

      INTEGER                           :: I, J, K
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
! =========================================================================================================
END PROGRAM
