      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,SIX=6.D0,
     1 NINE=9.D0,PI=3.141592653589793D0)
C
      DIMENSION EELAS(6),EPLAS(6),FLOW(6),SDEV(6),
     1 HARD(3),STRESS_TRIAL(6)
C
C     Get material properties
      YOUNG = PROPS(1)    ! Young's modulus
      POISSON = PROPS(2)  ! Poisson's ratio
      PHI = PROPS(3)      ! Friction angle in degrees
      PSI = PROPS(4)      ! Dilation angle in degrees
      COHESION = PROPS(5) ! Cohesion
C
C     Convert angles to radians
      PHI = PHI*PI/180.0D0
      PSI = PSI*PI/180.0D0
C
C     Calculate elastic constants
      BULK = YOUNG/(THREE*(ONE-TWO*POISSON))
      G = YOUNG/(TWO*(ONE+POISSON))
      LAMBDA = BULK - TWO*G/THREE
C
C     Initialize elastic stiffness matrix
      DO I = 1,NTENS
        DO J = 1,NTENS
          DDSDDE(I,J) = ZERO
        END DO
      END DO
C
C     Build elastic stiffness matrix
      DO I = 1,NDI
        DO J = 1,NDI
          DDSDDE(I,J) = LAMBDA
        END DO
        DDSDDE(I,I) = LAMBDA + TWO*G
      END DO
C
      DO I = NDI+1,NTENS
        DDSDDE(I,I) = G
      END DO
C
C     Calculate trial stress
      DO I = 1,NTENS
        STRESS_TRIAL(I) = STRESS(I)
        DO J = 1,NTENS
          STRESS_TRIAL(I) = STRESS_TRIAL(I) + 
     1      DDSDDE(I,J)*DSTRAN(J)
        END DO
      END DO
C
C     Calculate principal stresses
      CALL SPRINC(STRESS_TRIAL,PRINCIPAL,DIRECTION,1,NDI,NSHR)
C
C     Mohr-Coulomb parameters
      SINPHI = SIN(PHI)
      COSPHI = COS(PHI)
      SINPSI = SIN(PSI)
      TANPHI = TAN(PHI)
      TANPSI = TAN(PSI)
C
C     Calculate Mohr-Coulomb yield function
      FK = (ONE+SINPHI)/(ONE-SINPHI)
      F = PRINCIPAL(1) - FK*PRINCIPAL(3) - 
     1    TWO*COHESION*SQRT(FK)
C
      IF (F .GT. ZERO) THEN
C       Plastic state
        DLAM = F/(TWO*G*(ONE+FK*FK) + TWO*LAMBDA*(FK-ONE))
C
C       Update stresses
        DP1 = DLAM
        DP3 = -DLAM*FK
C
        PRINCIPAL(1) = PRINCIPAL(1) - (TWO*G*DP1 + 
     1    LAMBDA*(DP1+DP3))
        PRINCIPAL(2) = PRINCIPAL(2) - LAMBDA*(DP1+DP3)
        PRINCIPAL(3) = PRINCIPAL(3) - (TWO*G*DP3 + 
     1    LAMBDA*(DP1+DP3))
C
C       Transform back to global coordinates
        CALL SCLA(PRINCIPAL,DIRECTION,STRESS,NDI,NSHR)
      ELSE
C       Elastic state
        DO I = 1,NTENS
          STRESS(I) = STRESS_TRIAL(I)
        END DO
      END IF
C
      RETURN
      END
