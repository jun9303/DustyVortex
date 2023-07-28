MODULE PARTICLEINTF ! MODULE FOR PARTICLE INTERFACE (WRT FLOW) SUBROUTINES & FUNCTIONS
      USE OMP_LIB
      USE MISC; USE MATOPS
      USE LEGENDRE; USE FOURIER; USE FIELDGEN; USE FIELDOPS
      USE PARTICLEGEN
      IMPLICIT NONE
      PRIVATE
! ==========================================================================================================
!  DECLARATIONS ============================================================================================
! ==========================================================================================================
      PUBLIC :: PACCCOMP
      PUBLIC :: PVOLFRAC
      PUBLIC :: PVELFLD
      PUBLIC :: PFORCEFLD
! ==========================================================================================================
!  INTERFACES ==============================================================================================
! ==========================================================================================================
      INTERFACE PACCCOMP
        MODULE PROCEDURE PACCCOMP_PARTICLE
      END INTERFACE

      INTERFACE PVOLFRAC
        MODULE PROCEDURE PVOLFRAC_PARTICLE_SCALAR
      END INTERFACE

      INTERFACE PVELFLD
        MODULE PROCEDURE PVELFLD_PARTICLE_VECTOR
      END INTERFACE

      INTERFACE PFORCEFLD
        MODULE PROCEDURE PFORCEFLD_PARTICLE_VECTOR
      END INTERFACE

CONTAINS
! ==========================================================================================================
!  PUBLIC PROCEDURES =======================================================================================
! ==========================================================================================================
      SUBROUTINE PACCCOMP_PARTICLE(PTCL, VELFLD, FORCEFLD)
! ==========================================================================================================
! [USAGE]:
! COMPUTE THE ACCELERATION COMPONENTS OF EACH PARTICLE IN PTCL AS A RESULT OF FLUID-PARTICLE INTERACTIONS.
! [VARIABLES]:
! VELFLD >> SURROUNDING FLUID'S VELOCITY FIELD COMPUTED FROM THE FIELD MODULES, TO COMPUTE THE STOKES FORCE.
! FORCEFLD >> (OPTIONAL) FLUID'S INTERNAL FORCE FIELD. IF PRESENT, COMPUTES THE MOMENTUM EXCHANGE.
! [NOTES]:
! DU_P/DT = 1/DENS_P*DIV(STRESS_TENSOR)@X=PARTICLE + 1/RELXT_P*(U_F@X=PARTICLE - U_P)
!           <------ MOMENTUM EXCHANGE TERM ------>    <----- STOKES FORCE TERM ------>
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(PARTICLE), INTENT(INOUT)         :: PTCL
      TYPE(VECTOR_FIELD), INTENT(INOUT)     :: VELFLD
      TYPE(VECTOR_FIELD), OPTIONAL          :: FORCEFLD

      REAL(P8), DIMENSION(:), ALLOCATABLE   :: RPTS, PPTS, ZPTS
      REAL(P8), DIMENSION(:,:), ALLOCATABLE :: MOMEXC, STOKES
      TYPE(VECTOR_FIELD)                    :: FS

      RPTS = SQRT(PTCL%POS(1,:)**2.D0 + PTCL%POS(2,:)**2.D0)
      PPTS = ATAN2(PTCL%POS(2,:),PTCL%POS(1,:)); WHERE (PPTS .LT. 0.D0) PPTS = PPTS + 2.D0*PI
      ZPTS = MODULO(PTCL%POS(3,:),FINFO%ZLEN)

      IF (PRESENT(FORCEFLD)) THEN
        FS = VEL2FORCE(VELFLD, PSILN=-.5D0, CHILN=0.D0, F=FORCEFLD)
      ELSE
        FS = VEL2FORCE(VELFLD, PSILN=-.5D0, CHILN=0.D0)
      ENDIF
      MOMEXC = CALCAT_LINTP(FS, RPTS, PPTS, ZPTS) / MAX(PINFO%DENS, 1.D-30)
      STOKES = CALCAT_LINTP(VELFLD, RPTS, PPTS, ZPTS) / MAX(PINFO%RLXT, 1.D-30)

      MOMEXC = MOMEXC + STOKES

      PTCL%ACC(1,:) = MOMEXC(1,:)*COS(PPTS) - MOMEXC(2,:)*SIN(PPTS)
      PTCL%ACC(2,:) = MOMEXC(1,:)*SIN(PPTS) + MOMEXC(2,:)*COS(PPTS)
      PTCL%ACC(3,:) = MOMEXC(3,:)

      PTCL%ACC = PTCL%ACC - PTCL%VEL / MAX(PINFO%RLXT, 1.D-30)
      
      CALL DEALLOC(FS)
      DEALLOCATE( RPTS, PPTS, ZPTS )
      DEALLOCATE( MOMEXC, STOKES )

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE PVOLFRAC_PARTICLE_SCALAR(PTCL, PVOL)
! ==========================================================================================================
! [USAGE]:
! COMPUTE THE PARTICLE VOLUME FRACTION IN A EULERIAN MANNER.
! [PARAMETERS]:
! PVOL >> V_FRACTION_PARTICLE
! [NOTES]
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE

      TYPE(PARTICLE), INTENT(INOUT)       :: PTCL
      TYPE(SCALAR_FIELD), INTENT(INOUT)   :: PVOL

      INTEGER                             :: I, J, K, N, M
      REAL(P8), DIMENSION(3)              :: XO, X
      REAL(P8)                            :: DP, DZ, R2, SIGMA
      REAL(P8)                            :: INTE

      CALL ALLOC(PVOL, 'PPP')
      PVOL%E = 0.D0

      DP = 2.D0*PI / FINFO%NP
      DZ = FINFO%ZLEN / FINFO%NZ
      SIGMA = (FINFO%ELL/FINFO%NRH*DP*DZ)**(1.D0/3.D0)

      DO N = 1, PINFO%NPARTS
        XO(1) = PTCL%POS(1,N)
        XO(2) = PTCL%POS(2,N)
        XO(3) = PTCL%POS(3,N)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,K,X,M,R2)
        DO I = 1, FINFO%NR
          DO J = 1, FINFO%NP
            DO K = 1, FINFO%NZ
              X(1) = FINFO%R(I)*COS(FINFO%P(J))
              X(2) = FINFO%R(I)*SIN(FINFO%P(J))
              X(3) = FINFO%Z(K)
              DO M = -1, 1
                R2 = ((X(1)-XO(1))**2.D0+(X(2)-XO(2))**2.D0+(X(3)-XO(3)-M*FINFO%ZLEN)**2.D0)/SIGMA**2.D0
                IF (R2 .LT. 64.D0) THEN
                  IF (MOD(J,2) .EQ. 0) THEN
                    PVOL%E(I,(J+1)/2,K) = PVOL%E(I,(J+1)/2,K) + EXP(-.5D0*R2)/SIGMA**3.D0*IU
                  ELSE
                    PVOL%E(I,(J+1)/2,K) = PVOL%E(I,(J+1)/2,K) + EXP(-.5D0*R2)/SIGMA**3.D0
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDDO

      INTE = INTEG(PVOL); IF (INTE .EQ. 0.D0) INTE = 1.D-100
      PVOL%E = PVOL%E/INTE*PINFO%NPARTS*PINFO%VOL

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE PVELFLD_PARTICLE_VECTOR(PTCL, PVFLD)
! ==========================================================================================================
! [USAGE]:
! COMPUTE THE PARTICLE VELOCITY VECTOR FIELD IN A EULERIAN MANNER.
! [PARAMETERS]:
! PVFLD >> PARTICLE VELOCITY VECTOR FIELD
! [NOTES]
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE


      TYPE(PARTICLE), INTENT(INOUT)       :: PTCL
      TYPE(VECTOR_FIELD), INTENT(INOUT)   :: PVFLD

      INTEGER                             :: I, J, K, N, M
      REAL(P8), DIMENSION(3)              :: XO, X
      REAL(P8)                            :: DP, DZ, R2, SIGMA
      REAL(P8), DIMENSION(3)              :: INTE

      CALL ALLOC(PVFLD, 'PPP')
      PVFLD%ER = 0.D0
      PVFLD%EP = 0.D0
      PVFLD%EZ = 0.D0

      DP = 2.D0*PI / FINFO%NP
      DZ = FINFO%ZLEN / FINFO%NZ
      SIGMA = (FINFO%ELL/FINFO%NRH*DP*DZ)**(1.D0/3.D0)

      DO N = 1, PINFO%NPARTS
        XO(1) = PTCL%POS(1,N)
        XO(2) = PTCL%POS(2,N)
        XO(3) = PTCL%POS(3,N)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,K,X,M,R2)
        DO I = 1, FINFO%NR
          DO J = 1, FINFO%NP
            DO K = 1, FINFO%NZ
              X(1) = FINFO%R(I)*COS(FINFO%P(J))
              X(2) = FINFO%R(I)*SIN(FINFO%P(J))
              X(3) = FINFO%Z(K)
              DO M = -1, 1
                R2 = ((X(1)-XO(1))**2.D0+(X(2)-XO(2))**2.D0+(X(3)-XO(3)-M*FINFO%ZLEN)**2.D0)/SIGMA**2.D0
                IF (R2 .LT. 64.D0) THEN
                  IF (MOD(J,2) .EQ. 0) THEN
                    PVFLD%ER(I,(J+1)/2,K)=PVFLD%ER(I,(J+1)/2,K)+EXP(-.5D0*R2)/SIGMA**3.D0*IU*PTCL%VEL(1,N)
                    PVFLD%EP(I,(J+1)/2,K)=PVFLD%EP(I,(J+1)/2,K)+EXP(-.5D0*R2)/SIGMA**3.D0*IU*PTCL%VEL(2,N)
                    PVFLD%EZ(I,(J+1)/2,K)=PVFLD%EZ(I,(J+1)/2,K)+EXP(-.5D0*R2)/SIGMA**3.D0*IU*PTCL%VEL(3,N)
                  ELSE
                    PVFLD%ER(I,(J+1)/2,K)=PVFLD%ER(I,(J+1)/2,K)+EXP(-.5D0*R2)/SIGMA**3.D0*PTCL%VEL(1,N)
                    PVFLD%EP(I,(J+1)/2,K)=PVFLD%EP(I,(J+1)/2,K)+EXP(-.5D0*R2)/SIGMA**3.D0*PTCL%VEL(2,N)
                    PVFLD%EZ(I,(J+1)/2,K)=PVFLD%EZ(I,(J+1)/2,K)+EXP(-.5D0*R2)/SIGMA**3.D0*PTCL%VEL(3,N)
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDDO

      INTE = INTEG(PVFLD); WHERE (INTE .EQ. 0.D0) INTE = 1.D-100
      PVFLD%ER = PVFLD%ER/INTE(1)*SUM(PTCL%VEL(1,:))*PINFO%VOL
      PVFLD%EP = PVFLD%EP/INTE(2)*SUM(PTCL%VEL(2,:))*PINFO%VOL
      PVFLD%EZ = PVFLD%EZ/INTE(3)*SUM(PTCL%VEL(3,:))*PINFO%VOL

      RETURN
      END SUBROUTINE
! ==========================================================================================================
      SUBROUTINE PFORCEFLD_PARTICLE_VECTOR(PTCL, PFFLD)
! ==========================================================================================================
! [USAGE]:
! COMPUTE THE PARTICLE-FLUID INTERACTION FORCE VECTOR FIELD IN A EULERIAN MANNER.
! [PARAMETERS]:
! PFFLD >> FLUID_PARTICLE_INTERACT. FORCE VECTOR FIELD
! [NOTES]
! SANGJOON LEE @ JUNE 2023
! ==========================================================================================================
      IMPLICIT NONE


      TYPE(PARTICLE), INTENT(INOUT)       :: PTCL
      TYPE(VECTOR_FIELD), INTENT(INOUT)   :: PFFLD

      INTEGER                             :: I, J, K, N, M
      REAL(P8), DIMENSION(3)              :: XO, X
      REAL(P8)                            :: DP, DZ, R2, SIGMA
      REAL(P8), DIMENSION(3)              :: INTE

      CALL ALLOC(PFFLD, 'PPP')
      PFFLD%ER = 0.D0
      PFFLD%EP = 0.D0
      PFFLD%EZ = 0.D0

      DP = 2.D0*PI / FINFO%NP
      DZ = FINFO%ZLEN / FINFO%NZ
      SIGMA = (FINFO%ELL/FINFO%NRH*DP*DZ)**(1.D0/3.D0)

      DO N = 1, PINFO%NPARTS
        XO(1) = PTCL%POS(1,N)
        XO(2) = PTCL%POS(2,N)
        XO(3) = PTCL%POS(3,N)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,K,X,M,R2)
        DO I = 1, FINFO%NR
          DO J = 1, FINFO%NP
            DO K = 1, FINFO%NZ
              X(1) = FINFO%R(I)*COS(FINFO%P(J))
              X(2) = FINFO%R(I)*SIN(FINFO%P(J))
              X(3) = FINFO%Z(K)
              DO M = -1, 1
                R2 = ((X(1)-XO(1))**2.D0+(X(2)-XO(2))**2.D0+(X(3)-XO(3)-M*FINFO%ZLEN)**2.D0)/SIGMA**2.D0
                IF (R2 .LT. 64.D0) THEN
                  IF (MOD(J,2) .EQ. 0) THEN
                    PFFLD%ER(I,(J+1)/2,K)=PFFLD%ER(I,(J+1)/2,K)+EXP(-.5D0*R2)/SIGMA**3.D0*IU*PTCL%ACC(1,N)
                    PFFLD%EP(I,(J+1)/2,K)=PFFLD%EP(I,(J+1)/2,K)+EXP(-.5D0*R2)/SIGMA**3.D0*IU*PTCL%ACC(2,N)
                    PFFLD%EZ(I,(J+1)/2,K)=PFFLD%EZ(I,(J+1)/2,K)+EXP(-.5D0*R2)/SIGMA**3.D0*IU*PTCL%ACC(3,N)
                  ELSE
                    PFFLD%ER(I,(J+1)/2,K)=PFFLD%ER(I,(J+1)/2,K)+EXP(-.5D0*R2)/SIGMA**3.D0*PTCL%ACC(1,N)
                    PFFLD%EP(I,(J+1)/2,K)=PFFLD%EP(I,(J+1)/2,K)+EXP(-.5D0*R2)/SIGMA**3.D0*PTCL%ACC(2,N)
                    PFFLD%EZ(I,(J+1)/2,K)=PFFLD%EZ(I,(J+1)/2,K)+EXP(-.5D0*R2)/SIGMA**3.D0*PTCL%ACC(3,N)
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDDO

      INTE = INTEG(PFFLD); WHERE (INTE .EQ. 0.D0) INTE = 1.D-100
      PFFLD%ER = PFFLD%ER/INTE(1)*SUM(PTCL%ACC(1,:))*PINFO%VOL*PINFO%DENS*(-1.D0)
      PFFLD%EP = PFFLD%EP/INTE(2)*SUM(PTCL%ACC(2,:))*PINFO%VOL*PINFO%DENS*(-1.D0)
      PFFLD%EZ = PFFLD%EZ/INTE(3)*SUM(PTCL%ACC(3,:))*PINFO%VOL*PINFO%DENS*(-1.D0)

      RETURN
      END SUBROUTINE
! ==========================================================================================================
END MODULE