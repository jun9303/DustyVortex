! =========================================================================================================
!  WRITTEN BY SANGJOON (JOON) LEE, DEPT. OF MECHANICAL ENGINEERING, UNIVERSITY OF CALIFORNIA, BERKELEY    =
!  CONTACT: SANGJOONLEE@BERKELEY.EDU                                                                      =
!  UC BERKELEY COMPUTATIONAL FLUID DYNAMICS (CFD) LAB, HTTPS://CFD.ME.BERKELEY.EDU/                       =
!  SANGJOON LEE (C) 2023                                                                                  =
!  NONCOMMERCIAL USE. FOR RESEARCH PURPOSES ONLY                                                          =
! =========================================================================================================
      PROGRAM DUSTYQVORTRUN
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
      TYPE(VECTOR_FIELD)            :: F, FP
      TYPE(VECTOR_FIELD)            :: U, UP
      TYPE(PARTICLE)                :: PTCL, PTCLP

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
      CALL PARTICLE_INIT(DIAMIN=P%DIAM, DENSIN=P%DENS, NPARTSIN=P%NPARTS, VISCIN=P%VISC)

      IF (MPI_RANK .EQ. 0) THEN
        FINISH = OMP_GET_WTIME() ! OMP_LIB
        WRITE(*,101) 'FIELD, TIME & PARTICLE INITALIZATION'//' '//REPEAT('.',48), FINISH-START
        START = OMP_GET_WTIME() ! OMP_LIB
      ENDIF

! ....................................................................................................... !
      IF (P%ISPTCLTRCK) OPEN(UNIT=11, FILE=TRIM(ADJUSTL(P%LOGDIR))//'ptcl.log')
! ....................................................................................................... !

      CALL SYSTEM('mkdir -p '//TRIM(ADJUSTL(P%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))))

      CALL MLOAD(TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/u.vfld'  , U   )
      CALL MLOAD(TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/psi.sfld', PSI )
      CALL MLOAD(TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/chi.sfld', CHI )
      
      CALL MLOAD(TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/fp.vfld' , F   )
      CALL PLOAD(TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/ptcl'    , PTCL)

      CALL OUTPUT3D(U, PSI, CHI, F, PTCL, &
                    TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/tec3d.dat')

      IF (MPI_RANK .EQ. 0) THEN
        FINISH = OMP_GET_WTIME() ! OMP_LIB
        WRITE(*,101) 'INITIAL FIELDS LOADED'//' '//REPEAT('.',48), FINISH-START
        START = OMP_GET_WTIME() ! OMP_LIB
      ENDIF

      IF (P%ISPTCLTRCK) THEN
        WRITE(11,'(1F15.6,<3*P%PTCLTRCKN>F15.8)') TINFO%T, (PTCL%POS(:3,I), I = 1, P%PTCLTRCKN)
      ENDIF

! ....................................................................................................... !

      PSIP = PSI; CHIP = CHI; UP = U
      PTCLP = PTCL; FP = F
      CALL PTCSTEP_BOOTSTR(PSI, CHI, PTCL, F, N=1, ORDER=2)
      U = PT2VEL(PSI, CHI)

      TINFO%N = TINFO%N + 1
      IF (TINFO%N - TINFO%NI .EQ. TINFO%TOTN) THEN
        TINFO%T = TINFO%T + TINFO%DTF
      ELSE
        TINFO%T = TINFO%T + TINFO%DT
      ENDIF

      IF (MPI_RANK .EQ. 0) THEN
        FINISH = OMP_GET_WTIME() ! OMP_LIB
        WRITE(*,101) '1ST TIME STEP BOOTSTRAP - FLOW & PARTICLES'//' '//REPEAT('.',48), FINISH-START
        START = OMP_GET_WTIME() ! OMP_LIB
      ENDIF

      IF (P%ISPTCLTRCK) THEN
        WRITE(11,'(1F15.6,<3*P%PTCLTRCKN>F15.8)') TINFO%T, (PTCL%POS(:3,I), I = 1, P%PTCLTRCKN)
      ENDIF

! ....................................................................................................... !
! ....................................................................................................... !
      DO WHILE (TINFO%N .LT. TINFO%NI+TINFO%TOTN)
! ....................................................................................................... !

      UP = U
      CALL FLDSTEP_SEMIIMP(PSI=PSI, CHI=CHI, PSIP=PSIP, CHIP=CHIP, ORDER=2, F=F, FP=FP)
      U = PT2VEL(PSI, CHI)

      CALL PTCSTEP_EXPLICIT(PTCL, PTCLP)

      FP = F
      CALL PFFLDCOMP(PTCL=PTCL, VFLD=U, PFFLD=F, EPS=P%PFFEPS, ITERM=P%PFFITER, PFFLD_EST=FP)

      TINFO%N = TINFO%N+1
      IF (TINFO%N - TINFO%NI .EQ. TINFO%TOTN) THEN
        TINFO%T = TINFO%T + TINFO%DTF
      ELSE
        TINFO%T = TINFO%T + TINFO%DT
      ENDIF
 
      IF (P%ISFFRMV) THEN
        IF (MOD(TINFO%N, P%FFRMVINTVL) .EQ. 0) THEN
          CALL PTFFRMV(PSI=PSI, CHI=CHI, U=U)
          CALL PTFFRMV(PSI=PSIP, CHI=CHIP, U=UP)
          IF (MPI_RANK .EQ. 0) THEN
            WRITE(*,102) 'FAR-FIELD REMOVAL @ TIMESTEP # '//NTOA(TINFO%N,'(I10)')&
                         //' '//REPEAT('.',72)
          ENDIF
        ENDIF
      ENDIF

      IF (P%ISPTCLTRCK) THEN
        WRITE(11,'(1F15.6,<3*P%PTCLTRCKN>F15.8)') TINFO%T, (PTCL%POS(:3,I), I = 1, P%PTCLTRCKN)
      ENDIF

      IF (P%ISDATSAV) THEN
        IF (MOD(TINFO%N, P%DATSAVINTVL) .EQ. 0) THEN
          CALL SYSTEM('mkdir -p '//TRIM(ADJUSTL(P%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))))
          CALL MSAVE(U,   TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/u.vfld')
          CALL MSAVE(PSI, TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/psi.sfld')
          CALL MSAVE(CHI, TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/chi.sfld')

          CALL MSAVE(F,   TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/fp.vfld')
          CALL PSAVE(PTCL,TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/ptcl')

          CALL OUTPUT3D(U, PSI, CHI, F, PTCL, &
                        TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/tec3d.dat')

          IF (MPI_RANK .EQ. 0) THEN
            WRITE(*,102) 'FIELD & PTCL DATA STORAGE @ TIMESTEP # '//NTOA(TINFO%N,'(I10)')&
                         //' '//REPEAT('.',72)
          ENDIF
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

! ....................................................................................................... !
      CLOSE(11)
! ....................................................................................................... !

 101  FORMAT(A48,2X,F13.6,' SECONDS.')
 102  FORMAT(A72)

CONTAINS
! =========================================================================================================
!  PROGRAM-DEPENDENT PROCEDRUES ===========================================================================
! =========================================================================================================
      SUBROUTINE OUTPUT3D(U, PSI, CHI, F, PTCL, OUTPUTPATH)
! =========================================================================================================
      IMPLICIT NONE

      TYPE(VECTOR_FIELD), INTENT(IN) :: U, F
      TYPE(SCALAR_FIELD), INTENT(IN) :: PSI, CHI
      TYPE(PARTICLE), INTENT(INOUT)  :: PTCL
      CHARACTER(LEN=*), INTENT(IN)   :: OUTPUTPATH

      TYPE(SCALAR_FIELD)             :: PTCLVOL, PSIPPP, CHIPPP
      TYPE(VECTOR_FIELD)             :: O

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

      CALL PVOLFRAC(PTCL, PTCLVOL)
      
      PSIPPP = PSI
      CALL TRANS(PSIPPP, 'PPP')

      CHIPPP = CHI
      CALL TRANS(CHIPPP, 'PPP')

      O = PT2VOR(PSI, CHI)

      OPEN(FU, FILE=OUTPUTPATH)

      WRITE(FU,'(A110)')'variables= "x","y","z","ux","uy","uz","wx","wy","wz","fpx","fpy","fpz","ptcl"'
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
                            REAL(O%ER(I,(JJ+1)/2,KK))*COS(FINFO%P(JJ))&
                              -REAL(O%EP(I,(JJ+1)/2,KK))*SIN(FINFO%P(JJ)),&
                            REAL(O%ER(I,(JJ+1)/2,KK))*SIN(FINFO%P(JJ))&
                              +REAL(O%EP(I,(JJ+1)/2,KK))*COS(FINFO%P(JJ)),&
                            REAL(O%EZ(I,(JJ+1)/2,KK)),&
                            REAL(F%ER(I,(JJ+1)/2,KK))*COS(FINFO%P(JJ))&
                              -REAL(F%EP(I,(JJ+1)/2,KK))*SIN(FINFO%P(JJ)),&
                            REAL(F%ER(I,(JJ+1)/2,KK))*SIN(FINFO%P(JJ))&
                              +REAL(F%EP(I,(JJ+1)/2,KK))*COS(FINFO%P(JJ)),&
                            REAL(F%EZ(I,(JJ+1)/2,KK)),&
                            REAL(PTCLVOL%E(I,(JJ+1)/2,KK))
            ELSE
              WRITE(FU,103) FINFO%R(I)*COS(FINFO%P(JJ)),&
                            FINFO%R(I)*SIN(FINFO%P(JJ)),&
                            FINFO%Z(K),&
                            AIMAG(U%ER(I,(JJ+1)/2,KK))*COS(FINFO%P(JJ))&
                              -AIMAG(U%EP(I,(JJ+1)/2,KK))*SIN(FINFO%P(JJ)),&
                            AIMAG(U%ER(I,(JJ+1)/2,KK))*SIN(FINFO%P(JJ))&
                              +AIMAG(U%EP(I,(JJ+1)/2,KK))*COS(FINFO%P(JJ)),&
                            AIMAG(U%EZ(I,(JJ+1)/2,KK)),&
                            AIMAG(O%ER(I,(JJ+1)/2,KK))*COS(FINFO%P(JJ))&
                              -AIMAG(O%EP(I,(JJ+1)/2,KK))*SIN(FINFO%P(JJ)),&
                            AIMAG(O%ER(I,(JJ+1)/2,KK))*SIN(FINFO%P(JJ))&
                              +AIMAG(O%EP(I,(JJ+1)/2,KK))*COS(FINFO%P(JJ)),&
                            AIMAG(O%EZ(I,(JJ+1)/2,KK)),&
                            AIMAG(F%ER(I,(JJ+1)/2,KK))*COS(FINFO%P(JJ))&
                              -AIMAG(F%EP(I,(JJ+1)/2,KK))*SIN(FINFO%P(JJ)),&
                            AIMAG(F%ER(I,(JJ+1)/2,KK))*SIN(FINFO%P(JJ))&
                              +AIMAG(F%EP(I,(JJ+1)/2,KK))*COS(FINFO%P(JJ)),&
                            AIMAG(F%EZ(I,(JJ+1)/2,KK)),&
                            AIMAG(PTCLVOL%E(I,(JJ+1)/2,KK))
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      CLOSE(FU)

      CALL DEALLOC( PTCLVOL )
      CALL DEALLOC( PSIPPP )
      CALL DEALLOC( CHIPPP )

103   FORMAT(13E23.15)

      RETURN
      END SUBROUTINE
! =========================================================================================================
END PROGRAM
