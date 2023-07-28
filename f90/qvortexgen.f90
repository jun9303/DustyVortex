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
      USE MISC; USE MATOPS; USE TIMEINFO 
      USE FIELDGEN; USE FIELDOPS; USE FIELDSTEP

      IMPLICIT NONE

      INTEGER                       :: I, J, K
      REAL(P8)                      :: START, FINISH
      TYPE(SCALAR_FIELD)            :: PSI, CHI
      TYPE(VECTOR_FIELD)            :: U

      TYPE(INPUT_PARAMS)            :: P

      WRITE(*,*) ''
      WRITE(*,*) '*** PROGRAM STARTED ***'
      CALL PRINT_REAL_TIME() ! MISC
      CALL CPU_TIME(START) ! OMP_LIB
      
! ....................................................................................................... !

      CALL READ_INPUTPARAMS(PARAMS=P)
      CALL FIELD_INIT(NRIN=P%NR, NPIN=P%NP, NZIN=P%NZ, ELLIN=P%ELL, ZLENIN=P%ZLEN, &
                      VISCIN=P%VISC, VISCPOW=P%VISCPOW, VISCPIN=P%VISCP)
      CALL TIME_INIT(DTIN=P%DT, TIIN=P%TI, TOTTIN=P%TOTT, NIIN=P%NI, TOTNIN=P%TOTN)

      CALL CPU_TIME(FINISH) ! OMP_LIB
      WRITE(*,101) 'FIELD, TIME INITALIZATION'//' '//REPEAT('.',48), FINISH-START
      CALL CPU_TIME(START) ! OMP_LIB

! ....................................................................................................... !

      CALL FIELD_SET(VFLD=U, Q=P%Q, XO=P%XQO, YO=P%YQO)

      CALL VELPROJ(U,PSI,CHI,PSILN=-.5)
      U = PT2VEL(PSI,CHI)

      CALL CPU_TIME(FINISH) ! OMP_LIB
      WRITE(*,101) 'Q-VORTEX VELOCITY FIELD SETUP'//' '//REPEAT('.',48), FINISH-START
      CALL CPU_TIME(START) ! OMP_LIB

! ....................................................................................................... !

      CALL EXECUTE_COMMAND_LINE('mkdir -p '//TRIM(ADJUSTL(P%DATDIR))//TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))))
      CALL MSAVE(U,   TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/u.vfld')
      CALL MSAVE(PSI, TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/psi.sfld')
      CALL MSAVE(CHI, TRIM(ADJUSTL(P%DATDIR)) // TRIM(ADJUSTL(NTOA(TINFO%T,'(F8.3)'))) // '/chi.sfld')

      CALL CPU_TIME(FINISH) ! OMP_LIB
      WRITE(*,101) 'INITIAL FIELDS STORED'//' '//REPEAT('.',48), FINISH-START
      CALL CPU_TIME(START) ! OMP_LIB

! ....................................................................................................... !

      WRITE(*,*) ''
      WRITE(*,*) '*** PROGRAM FINISHED ***'
      CALL PRINT_REAL_TIME() ! MISC

 101  FORMAT(A48,1X,F13.6,' SECONDS.')

CONTAINS
! =========================================================================================================
!  PROGRAM-DEPENDENT PROCEDRUES ===========================================================================
! =========================================================================================================

! =========================================================================================================
END PROGRAM
