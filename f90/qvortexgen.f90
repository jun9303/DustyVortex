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
      USE PARTICLESTEP

      IMPLICIT NONE

      INTEGER                       :: MM, NN, KK
      REAL(P8)                      :: START, FINISH
      TYPE(SCALAR_FIELD)            :: PSI, CHI, PSI_UNP, CHI_UNP, C, PSIVIS, CHIVIS, PSIUXW, CHIUXW
      TYPE(SCALAR_FIELD)            :: PSIP, CHIP
      TYPE(VECTOR_FIELD)            :: U, U_UNP, W, W_UNP, F, UXW, UXWSOL, VIS

      TYPE(INPUT_PARAMS)            :: PRMS

      LOGICAL                       :: ISOPTIM
      COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE :: OPTIM   

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

      CALL FIELD_SET(VFLD=U_UNP, Q=PRMS%Q, XO=PRMS%XQO, YO=PRMS%YQO)
      CALL FIELD_SET(SFLD=C, TYPE=PRMS%PDISTTYPE, CENTER=(/PRMS%XPC, PRMS%YPC, PRMS%ZPC/), &
                     CL=PRMS%PSCALE, INTENSITY=PRMS%PINTENSITY)

      CALL VELPROJ(U_UNP,PSI_UNP,CHI_UNP,PSILN=-.5)
      U_UNP = PT2VEL(PSI_UNP,CHI_UNP)
      W_UNP = PT2VOR(PSI_UNP,CHI_UNP)

      PSI = PSI_UNP
      CHI = CHI_UNP

      ! INQUIRE( FILE='optpert.mtrx', EXIST=ISOPTIM)
      ! IF ( ISOPTIM ) THEN
      !   CALL ALLOC(PSI, 'FFF'); CALL ALLOC(CHI, 'FFF'); CALL ALLOC(C, 'PPP')
      !   ALLOCATE( OPTIM(FINFO%NRCHOP*2,1) )
      !   CALL MLOAD('optpert.mtrx', OPTIM)
      !   ! CALL MCAT(OPTIM)
      !   PSI%E(:FINFO%NRCHOP,2,1) = OPTIM(:FINFO%NRCHOP,1)*PRMS%PINTENSITY
      !   CHI%E(:FINFO%NRCHOP,2,1) = OPTIM(FINFO%NRCHOP+1:,1)*PRMS%PINTENSITY
      !   CHI = IDEL2(CHI, ISTP=.TRUE.)
      !   CHI%E = CHI%E + CHI_UNP%E
      !   CHI%LN = CHI%LN + CHI_UNP%LN
      !   PSI%E = PSI%E + PSI_UNP%E
      !   PSI%LN = PSI%LN + PSI_UNP%LN
      !   DEALLOCATE( OPTIM )
      ! ENDIF

      U = PT2VEL(PSI, CHI)
      W = PT2VOR(PSI, CHI)

      UXW = PT2VOR(PSI, CHI) ! UXW NOW CONTAINS CURL(U)
      UXW = VPROD(U, UXW)      ! UXW NOW CONTAINS U X CURL(U) IN PPP
      CALL LTERM(PSI, CHI, PSIVIS, CHIVIS)
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

      CALL DEALLOC( UXWSOL ); CALL DEALLOC( UXW ); CALL DEALLOC( VIS )
      CALL DEALLOC( PSIVIS ); CALL DEALLOC( CHIVIS ); CALL DEALLOC( PSIUXW ); CALL DEALLOC( CHIUXW )

      IF (MPI_RANK .EQ. 0) THEN
        FINISH = OMP_GET_WTIME() ! OMP_LIB
        WRITE(*,101) 'Q-VORTEX VELOCITY FIELD SETUP'//' '//REPEAT('.',48), FINISH-START
        START = OMP_GET_WTIME() ! OMP_LIB
      ENDIF

! ....................................................................................................... !

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

      IF (MPI_RANK .EQ. 0) THEN
        FINISH = OMP_GET_WTIME() ! OMP_LIB
        WRITE(*,101) 'INITIAL FIELDS STORED'//' '//REPEAT('.',48), FINISH-START
        START = OMP_GET_WTIME() ! OMP_LIB
      ENDIF

! ....................................................................................................... !

      PSIP = PSI; CHIP = CHI
      CALL PFLDSTEP_BOOTSTR(PSI, CHI, C=C, N=1, ORDER=2)
      U = PT2VEL(PSI, CHI)
      W = PT2VOR(PSI, CHI)

      CALL FLDSTEP_BOOTSTR(PSI_UNP, CHI_UNP, N=1, ORDER=2)
      U_UNP = PT2VEL(PSI_UNP, CHI_UNP)
      W_UNP = PT2VOR(PSI_UNP, CHI_UNP)

      F = PFTERM(C, PSI, CHI, PSIP, CHIP)

      TINFO%N = TINFO%N + 1
      IF (TINFO%N - TINFO%NI .EQ. TINFO%TOTN) THEN
        TINFO%T = TINFO%T + TINFO%DTF
      ELSE
        TINFO%T = TINFO%T + TINFO%DT
      ENDIF

      IF (MPI_RANK .EQ. 0) THEN
        FINISH = OMP_GET_WTIME() ! OMP_LIB
        WRITE(*,101) '1ST TIME STEP BOOTSTRAP'//' '//REPEAT('.',48), FINISH-START
        START = OMP_GET_WTIME() ! OMP_LIB
      ENDIF

! ....................................................................................................... !

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

      IF (MPI_RANK .EQ. 0) THEN
        FINISH = OMP_GET_WTIME() ! OMP_LIB
        WRITE(*,101) 'BOOTSTRAPPED FIELDS STORED'//' '//REPEAT('.',48), FINISH-START
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
END PROGRAM
