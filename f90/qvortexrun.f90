! =========================================================================================================
!  WRITTEN BY SANGJOON (JOON) LEE, DEPT. OF MECHANICAL ENGINEERING, UNIVERSITY OF CALIFORNIA, BERKELEY    =
!  CONTACT: SANGJOONLEE@BERKELEY.EDU                                                                      =
!  UC BERKELEY COMPUTATIONAL FLUID DYNAMICS (CFD) LAB, HTTPS://CFD.ME.BERKELEY.EDU/                       =
!  SANGJOON LEE (C) 2023                                                                                  =
!  NONCOMMERCIAL USE. FOR RESEARCH PURPOSES ONLY                                                          =
! =========================================================================================================
      PROGRAM QVORTEXRUN
! =========================================================================================================
! [NOTES]:
! SANGJOON LEE @ JUNE 2023
! =========================================================================================================
      USE OMP_LIB; USE MPI
      USE MISC; USE MATOPS; USE TIMEINFO 
      USE FIELDGEN; USE FIELDOPS; USE FIELDSTEP
      USE PARTICLEGEN; USE PARTICLEINTF; USE PARTICLESTEP

      IMPLICIT NONE

      INTEGER                       :: I, J, K
      REAL(P8)                      :: START, FINISH
      TYPE(SCALAR_FIELD)            :: PSI, CHI, PSIP, CHIP
      TYPE(VECTOR_FIELD)            :: U, UP

      TYPE(INPUT_PARAMS)            :: P

      IF (MPI_RANK .EQ. 0) THEN
        WRITE(*,*) ''
        WRITE(*,*) '*** PROGRAM STARTED ***'
        CALL PRINT_REAL_TIME() ! MISC
        START = OMP_GET_WTIME() ! OMP_LIB
      ENDIF

! ....................................................................................................... !  

      CALL READ_INPUTPARAMS(PARAMS=P)
      CALL FIELD_INIT(NRIN=P%NR, NPIN=P%NP, NZIN=P%NZ, ELLIN=P%ELL, ZLENIN=P%ZLEN, &
                      VISCIN=P%VISC, VISCPOW=P%VISCPOW, VISCPIN=P%VISCP)
      CALL TIME_INIT(DTIN=P%DT, TIIN=P%TI, TOTTIN=P%TOTT, NIIN=P%NI, TOTNIN=P%TOTN)

      IF (MPI_RANK .EQ. 0) THEN
        FINISH = OMP_GET_WTIME() ! OMP_LIB
        WRITE(*,101) 'FIELD, TIME & PARTICLE INITALIZATION'//' '//REPEAT('.',48), FINISH-START
        START = OMP_GET_WTIME() ! OMP_LIB
      ENDIF

! ....................................................................................................... !

      CALL SYSTEM('mkdir -p '//TRIM(ADJUSTL(P%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))))

      CALL MLOAD(TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/u.vfld'  , U   )
      CALL MLOAD(TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/psi.sfld', PSI )
      CALL MLOAD(TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/chi.sfld', CHI )
      
      CALL OUTPUT3D(U, PSI, CHI, &
                    TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/tec3d.dat')

      IF (MPI_RANK .EQ. 0) THEN
        FINISH = OMP_GET_WTIME() ! OMP_LIB
        WRITE(*,101) 'INITIAL FIELDS LOADED'//' '//REPEAT('.',48), FINISH-START
        START = OMP_GET_WTIME() ! OMP_LIB
      ENDIF

! ....................................................................................................... !

      PSIP = PSI; CHIP = CHI; UP = U
      CALL FLDSTEP_BOOTSTR(PSI, CHI, N=1, ORDER=2)
      U = PT2VEL(PSI, CHI)

      TINFO%N = TINFO%N + 1
      IF (TINFO%N - TINFO%NI .EQ. TINFO%TOTN) THEN
        TINFO%T = TINFO%T + TINFO%DTF
      ELSE
        TINFO%T = TINFO%T + TINFO%DT
      ENDIF

      IF (MPI_RANK .EQ. 0) THEN
        FINISH = OMP_GET_WTIME() ! OMP_LIB
        WRITE(*,101) '1ST TIME STEP BOOTSTRAP - FLOW'//' '//REPEAT('.',48), FINISH-START
        START = OMP_GET_WTIME() ! OMP_LIB
      ENDIF

! ....................................................................................................... !
! ....................................................................................................... !
      DO WHILE (TINFO%N .LT. TINFO%NI+TINFO%TOTN)
! ....................................................................................................... !

      UP = U
      CALL FLDSTEP_SEMIIMP(PSI=PSI, CHI=CHI, PSIP=PSIP, CHIP=CHIP, ORDER=2)
      U = PT2VEL(PSI, CHI)

      TINFO%N = TINFO%N+1
      IF (TINFO%N - TINFO%NI .EQ. TINFO%TOTN) THEN
        TINFO%T = TINFO%T + TINFO%DTF
      ELSE
        TINFO%T = TINFO%T + TINFO%DT
      ENDIF

      IF (P%ISDATSAV) THEN
        IF (MOD(TINFO%N, P%DATSAVINTVL) .EQ. 0) THEN
          CALL SYSTEM('mkdir -p '//TRIM(ADJUSTL(P%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))))
          CALL MSAVE(U,   TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/u.vfld')
          CALL MSAVE(PSI, TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/psi.sfld')
          CALL MSAVE(CHI, TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/chi.sfld')

          CALL OUTPUT3D(U, PSI, CHI, &
                        TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/tec3d.dat')
        ENDIF
      ENDIF

      IF (MPI_RANK .EQ. 0) THEN
        FINISH = OMP_GET_WTIME() ! OMP_LIB
        WRITE(*,101) 'TIME STEPPING # '//NTOA(TINFO%N,'(I10)')//' '//'(T: '//NTOA(TINFO%T,'(F11.6)')//')'&
                     //' '//REPEAT('.',48), FINISH-START
        START = OMP_GET_WTIME() ! OMP_LIB
      ENDIF

! ....................................................................................................... !
      ENDDO
! ....................................................................................................... !
! ....................................................................................................... !

      IF (MPI_RANK .EQ. 0) THEN
        WRITE(*,*) ''
        WRITE(*,*) '*** PROGRAM FINISHED ***'
        CALL PRINT_REAL_TIME() ! MISC
      ENDIF

      CLOSE(11)

 101  FORMAT(A48,1X,F13.6,' SECONDS.')

CONTAINS
! =========================================================================================================
!  PROGRAM-DEPENDENT PROCEDRUES ===========================================================================
! =========================================================================================================
      SUBROUTINE OUTPUT3D(U, PSI, CHI, OUTPUTPATH)
! =========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(IN) :: U
      TYPE(SCALAR_FIELD), INTENT(IN) :: PSI, CHI
      CHARACTER(LEN=*), INTENT(IN)   :: OUTPUTPATH

      INTEGER                        :: I, J, JJ, K, KK, INUM, JNUM, KNUM, RSKIP, PSKIP, ZSKIP
      INTEGER                        :: FU

      INUM = 0; JNUM = 0; KNUM = 0
      RSKIP = 1; PSKIP = 1; ZSKIP = 1

      DO I=1,FINFO%NR,RSKIP
        INUM=INUM+1
      ENDDO
      DO J=1,2*FINFO%NPH+1,PSKIP
        JNUM=JNUM+1
      ENDDO
      DO K=1,FINFO%NZ+1,ZSKIP
        KNUM=KNUM+1
      ENDDO

      OPEN(FU, FILE=OUTPUTPATH)

      WRITE(FU,'(A110)')'variables= "x","y","z","ux","uy","uz","psi","chi"'
      WRITE(FU,*)'ZONE T="ZONE1" , I=',INUM,', J=',JNUM,', K=',KNUM,', ZONETYPE=Ordered'
      WRITE(FU,*)'DATAPACKING=POINT'
      DO K=1,FINFO%NZ+1,ZSKIP
        DO J=1,2*FINFO%NPH+1,PSKIP
          DO I=1,FINFO%NR,RSKIP
            KK = K
            JJ = J
            IF (K.EQ.FINFO%NZ+1) KK = 1
            IF (J.EQ.2*FINFO%NPH+1) JJ = 1
            IF (MOD(J,2).EQ.1) THEN
              WRITE(FU,103) FINFO%R(I)*COS(FINFO%P(JJ)),&
                            FINFO%R(I)*SIN(FINFO%P(JJ)),&
                            FINFO%Z(K),&
                            REAL(U%ER(I,(JJ+1)/2,KK))*COS(FINFO%P(JJ))&
                              -REAL(U%EP(I,(JJ+1)/2,KK))*SIN(FINFO%P(JJ)),&
                            REAL(U%ER(I,(JJ+1)/2,KK))*SIN(FINFO%P(JJ))&
                              +REAL(U%EP(I,(JJ+1)/2,KK))*COS(FINFO%P(JJ)),&
                            REAL(U%EZ(I,(JJ+1)/2,KK)),&
                            REAL(PSI%E(I,(JJ+1)/2,KK)),&
                            REAL(CHI%E(I,(JJ+1)/2,KK))
            ELSE
              WRITE(FU,103) FINFO%R(I)*COS(FINFO%P(JJ)),&
                            FINFO%R(I)*SIN(FINFO%P(JJ)),&
                            FINFO%Z(K),&
                            AIMAG(U%ER(I,(JJ+1)/2,KK))*COS(FINFO%P(JJ))&
                              -AIMAG(U%EP(I,(JJ+1)/2,KK))*SIN(FINFO%P(JJ)),&
                            AIMAG(U%ER(I,(JJ+1)/2,KK))*SIN(FINFO%P(JJ))&
                              +AIMAG(U%EP(I,(JJ+1)/2,KK))*COS(FINFO%P(JJ)),&
                            AIMAG(U%EZ(I,(JJ+1)/2,KK)),&
                            AIMAG(PSI%E(I,(JJ+1)/2,KK)),&
                            AIMAG(CHI%E(I,(JJ+1)/2,KK))
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      CLOSE(FU)

103   FORMAT(8E23.15)

      RETURN
      END SUBROUTINE
! =========================================================================================================
END PROGRAM
