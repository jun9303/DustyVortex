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
      USE PARTICLESTEP

      IMPLICIT NONE

      INTEGER                       :: MM, NN, KK
      REAL(P8)                      :: START, FINISH
      TYPE(SCALAR_FIELD)            :: PSI, CHI, PSI_UNP, CHI_UNP, C, CP
      TYPE(SCALAR_FIELD)            :: PSIP, CHIP, PSIPP, CHIPP, PSIP_UNP, CHIP_UNP, PSIPP_UNP, CHIPP_UNP
      TYPE(VECTOR_FIELD)            :: F, FP, FPP, U, U_UNP, W, W_UNP

      TYPE(INPUT_PARAMS)            :: PRMS

      IF (MPI_RANK .EQ. 0) THEN
        WRITE(*,*) ''
        WRITE(*,*) '*** PROGRAM STARTED ***'
        CALL PRINT_REAL_TIME() ! MISC
        START = OMP_GET_WTIME() ! OMP_LIB
      ENDIF

! ....................................................................................................... !

      CALL READ_INPUTPARAMS(PARAMS=PRMS)
      CALL FIELD_INIT(NRIN=PRMS%NR, NPIN=PRMS%NP, NZIN=PRMS%NZ, ELLIN=PRMS%ELL, ZLENIN=PRMS%ZLEN, &
                      NRCHOPIN=PRMS%NRCHOP, NPCHOPIN=PRMS%NPCHOP, NZCHOPIN=PRMS%NZCHOP, &
                      VISCIN=PRMS%VISC, VISCPOW=PRMS%VISCPOW, VISCPIN=PRMS%VISCP, &
                      PTCRESPIN=PRMS%PTCRESP, PTCDENSIN=PRMS%PTCDENS)
      CALL TIME_INIT(DTIN=PRMS%DT, TIIN=PRMS%TI, TOTTIN=PRMS%TOTT, NIIN=PRMS%NI, TOTNIN=PRMS%TOTN)

      IF (MPI_RANK .EQ. 0) THEN
        FINISH = OMP_GET_WTIME() ! OMP_LIB
        WRITE(*,101) 'FIELD, TIME INITIALIZATION'//' '//REPEAT('.',48), FINISH-START
        START = OMP_GET_WTIME() ! OMP_LIB
      ENDIF

! ....................................................................................................... !

      ! TIMESTEP 0 (T==0)
      CALL MLOAD(TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/psi.sfld', PSIP)
      CALL MLOAD(TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/chi.sfld', CHIP)
      CALL MLOAD(TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/psi_unp.sfld', PSIP_UNP)
      CALL MLOAD(TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/chi_unp.sfld', CHIP_UNP)
      CALL MLOAD(TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/c.sfld', CP)
      CALL MLOAD(TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/f.vfld', FP)

      TINFO%N = TINFO%N + 1
      IF (TINFO%N - TINFO%NI .EQ. TINFO%TOTN) THEN
        TINFO%T = TINFO%T + TINFO%DTF
      ELSE
        TINFO%T = TINFO%T + TINFO%DT
      ENDIF

      ! TIMESTEP 1 (T==DT)
      CALL MLOAD(TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/psi.sfld', PSI)
      CALL MLOAD(TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/chi.sfld', CHI)
      CALL MLOAD(TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/psi_unp.sfld', PSI_UNP)
      CALL MLOAD(TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/chi_unp.sfld', CHI_UNP)
      CALL MLOAD(TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/c.sfld', C)
      CALL MLOAD(TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/f.vfld', F)

      TINFO%N = TINFO%N + 1
      IF (TINFO%N - TINFO%NI .EQ. TINFO%TOTN) THEN
        TINFO%T = TINFO%T + TINFO%DTF
      ELSE
        TINFO%T = TINFO%T + TINFO%DT
      ENDIF

      IF (MPI_RANK .EQ. 0) THEN
        FINISH = OMP_GET_WTIME() ! OMP_LIB
        WRITE(*,101) 'INITIAL FIELDS LOADED'//' '//REPEAT('.',48), FINISH-START
        START = OMP_GET_WTIME() ! OMP_LIB
      ENDIF

! ....................................................................................................... !
      DO WHILE (TINFO%N .LT. TINFO%NI+TINFO%TOTN)
! ....................................................................................................... !
      PSIPP = PSIP; CHIPP = CHIP ! PSIPP == PSI AT TIMESTEP N-1, CHIPP = CHI AT TIMESTEP N-1

      F%ER = F%ER*(FINFO%PTCDENS-1.D0); FP%ER = FP%ER*(FINFO%PTCDENS-1.D0)
      F%EP = F%EP*(FINFO%PTCDENS-1.D0); FP%EP = FP%EP*(FINFO%PTCDENS-1.D0)
      F%EZ = F%EZ*(FINFO%PTCDENS-1.D0); FP%EZ = FP%EZ*(FINFO%PTCDENS-1.D0)
      CALL FLDSTEP_SEMIIMP(PSI, CHI, PSIP, CHIP, ORDER=2, F=F, FP=FP) ! ON EXIT, PSI = PSI AT TIMESTEP N+1, CHI = CHI AT TIMESTEP N+1 & PSIP = PSI AT TIMESTEP N, CHIP = CHI AT TIMESTEP N

      F%ER = F%ER/(FINFO%PTCDENS-1.D0); FP%ER = FP%ER/(FINFO%PTCDENS-1.D0)
      F%EP = F%EP/(FINFO%PTCDENS-1.D0); FP%EP = FP%EP/(FINFO%PTCDENS-1.D0)
      F%EZ = F%EZ/(FINFO%PTCDENS-1.D0); FP%EZ = FP%EZ/(FINFO%PTCDENS-1.D0)
      CALL PFLDSTEP_EXPLICIT(C, CP, PSIP, CHIP, PSIPP, CHIPP, F, FP) ! ON EXIT, C = C AT TIMESTEP N+1, CP = C AT TIMESTEP N

      FP = F ! FP == F AT TIMESTEP N
      F = PFTERM(C, PSI, CHI, PSIP, CHIP, PSIPP, CHIPP) ! CALCULATE THE NONLINEAR COUPLING FORCING AT TIMESTEP N+1

      U = PT2VEL(PSI, CHI)
      W = PT2VOR(PSI, CHI)

      CALL FLDSTEP_SEMIIMP(PSI_UNP, CHI_UNP, PSIP_UNP, CHIP_UNP, ORDER=2)

      U_UNP = PT2VEL(PSI_UNP, CHI_UNP)
      W_UNP = PT2VOR(PSI_UNP, CHI_UNP)

      TINFO%N = TINFO%N+1
      IF (TINFO%N - TINFO%NI .EQ. TINFO%TOTN) THEN
        TINFO%T = TINFO%T + TINFO%DTF
      ELSE
        TINFO%T = TINFO%T + TINFO%DT
      ENDIF

      IF (PRMS%ISFFRMV) THEN
        IF (MOD(TINFO%N, PRMS%FFRMVINTVL) .EQ. 0) THEN
          CALL PTFFRMV(PSI=PSI, CHI=CHI)
          CALL PTFFRMV(PSI=PSIP, CHI=CHIP)
          CALL PTFFRMV(PSI=PSI_UNP, CHI=CHI_UNP)
          CALL PTFFRMV(PSI=PSIP_UNP, CHI=CHIP_UNP)
          IF (MPI_RANK .EQ. 0) THEN
            WRITE(*,102) 'FAR-FIELD REMOVAL @ TIMESTEP # '//NTOA(TINFO%N,'(I10)')&
                         //' '//REPEAT('.',72)
          ENDIF
        ENDIF
      ENDIF

      IF (PRMS%ISDATSAV) THEN
        IF (MOD(TINFO%N, PRMS%DATSAVINTVL) .EQ. 0) THEN
          CALL SYSTEM('mkdir -p '//TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)'))))
          CALL MSAVE(U,     TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/u.vfld')
          CALL MSAVE(U_UNP,     TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/u_unp.vfld')
          CALL MSAVE(W,     TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/w.vfld')
          CALL MSAVE(W_UNP,     TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/w_unp.vfld')
          CALL MSAVE(PSI,   TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/psi.sfld')
          CALL MSAVE(CHI,   TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/chi.sfld')
          CALL MSAVE(PSI_UNP,   TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/psi_unp.sfld')
          CALL MSAVE(CHI_UNP,   TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/chi_unp.sfld')
          CALL MSAVE(C,     TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/c.sfld')

          CALL MSAVE(F, TRIM(ADJUSTL(PRMS%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)')))//'/f.vfld')

          CALL OUTPUT3D(U, W, C, U_UNP, W_UNP, F, &
                        TRIM(ADJUSTL(PRMS%DATDIR)) & ! // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)'))) &
                        // '/tec3d_' // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.4)'))) // '.dat')
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

      IF (MPI_RANK .EQ. 0) THEN
        WRITE(*,*) ''
        WRITE(*,*) '*** PROGRAM FINISHED ***'
        CALL PRINT_REAL_TIME() ! MISC
      ENDIF

      CLOSE(11)

 101  FORMAT(A48,1X,F13.6,' SECONDS.')
 102  FORMAT(A72)

CONTAINS
! =========================================================================================================
!  PROGRAM-DEPENDENT PROCEDRUES ===========================================================================
! =========================================================================================================
END PROGRAM
