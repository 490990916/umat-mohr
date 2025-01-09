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
     1 DDSDDE(NTENS,NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),
     4 DFGRD0(3,3),DFGRD1(3,3),
     5 DDSDDT(NTENS),DRPLDE(NTENS)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,SIX=6.D0,
     1 TOLER=1.D-10)
C
      DIMENSION STRESS_TRIAL(6),PRINCIPAL(3),DIRECTION(3,3)
      DIMENSION DSTRESS(6), DE_PRIN(3), STRESS_OLD(6)
C
C     保存旧应力状态
      DO I = 1,NTENS
        STRESS_OLD(I) = STRESS(I)
      END DO
C
C     材料参数
      YOUNG = PROPS(1)    ! 弹性模量 200000MPa
      POISSON = PROPS(2)  ! 泊松比 0.3
      PHI = PROPS(3)      ! 摩擦角 30度
      PSI = PROPS(4)      ! 膨胀角 10度
      COHESION = PROPS(5) ! 黏聚力 1MPa
C
C     角度转换为弧度
      PI = 3.141592653589793D0
      PHI = PHI*PI/180.0D0
      PSI = PSI*PI/180.0D0
C
C     计算弹性常数
      BULK = YOUNG/(THREE*(ONE-TWO*POISSON))
      G = YOUNG/(TWO*(ONE+POISSON))
      LAMBDA = BULK - TWO*G/THREE
C
C     初始化弹性刚度矩阵
      DO I = 1,NTENS
        DO J = 1,NTENS
          DDSDDE(I,J) = ZERO
        END DO
      END DO
C
C     构建弹性刚度矩阵
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
C     计算弹性预测应力增量
      DO I = 1,NTENS
        DSTRESS(I) = ZERO
        DO J = 1,NTENS
          DSTRESS(I) = DSTRESS(I) + DDSDDE(I,J)*DSTRAN(J)
        END DO
      END DO
C
C     更新试探应力
      DO I = 1,NTENS
        STRESS_TRIAL(I) = STRESS(I) + DSTRESS(I)
      END DO
C
C     计算主应力
      CALL SPRIND(STRESS_TRIAL,PRINCIPAL,DIRECTION,1,NDI,NSHR)
C
C     计算摩尔库伦参数
      SINPHI = SIN(PHI)
      COSPHI = COS(PHI)
      SINPSI = SIN(PSI)
      COSPSI = COS(PSI)
      
C     摩尔库伦参数
      ETA = (ONE+SINPHI)/(ONE-SINPHI)      ! 摩擦系数
      ETAG = (ONE+SINPSI)/(ONE-SINPSI)     ! 膨胀系数
      COHES = TWO*COHESION*COSPHI/(ONE-SINPHI)
C
C     计算屈服函数
      F = PRINCIPAL(1) - ETA*PRINCIPAL(3) - COHES
C
      IF (F .GT. TOLER) THEN
C       塑性状态 - 使用迭代返回映射
        MAX_ITER = 50
        TOL_R = 1.D-10
        DLAM_TOT = ZERO
        
        DO ITER = 1, MAX_ITER
          H = TWO*G + TWO*LAMBDA*SINPSI + TWO*G*ETA*ETAG
          DLAM = F/H
          
C         限制每次迭代的塑性乘子增量
          DLAM = MIN(DLAM, 1.D-7)
          DLAM_TOT = DLAM_TOT + DLAM
          
C         更新主应力
          DE_PRIN(1) = -TWO*G*DLAM - LAMBDA*DLAM*(ONE+ETAG)
          DE_PRIN(2) = -LAMBDA*DLAM*(ONE+ETAG)
          DE_PRIN(3) = TWO*G*ETA*DLAM - LAMBDA*DLAM*(ONE+ETAG)
          
          DO I = 1,3
            PRINCIPAL(I) = PRINCIPAL(I) + DE_PRIN(I)
          END DO
          
C         重新计算屈服函数
          F = PRINCIPAL(1) - ETA*PRINCIPAL(3) - COHES
          
C         检查收敛性
          IF (ABS(F) .LT. TOL_R) THEN
            EXIT
          END IF
          
C         如果迭代不收敛或塑性乘子过大，减小时间增量
          IF (ITER .EQ. MAX_ITER .OR. DLAM_TOT .GT. 1.D-4) THEN
            PNEWDT = 0.25D0
            DO I = 1,NTENS
              STRESS(I) = STRESS_OLD(I)
            END DO
            RETURN
          END IF
        END DO
        
C       转换回全局坐标系
        CALL SPRIND(PRINCIPAL,DIRECTION,STRESS,-1,NDI,NSHR)
        
C       计算一致切线模量（使用亚切线刚度）
        ALPHA = 0.5D0  ! 亚切线因子
        DO I = 1,NTENS
          DO J = 1,NTENS
            DDSDDE(I,J) = ALPHA*DDSDDE(I,J)
          END DO
        END DO
        
C       检查应力状态的合理性
        IF (PRINCIPAL(1) .LT. PRINCIPAL(2) .OR. 
     1      PRINCIPAL(2) .LT. PRINCIPAL(3)) THEN
          PNEWDT = 0.25D0
          DO I = 1,NTENS
            STRESS(I) = STRESS_OLD(I)
          END DO
          RETURN
        END IF
        
      ELSE
C       弹性状态
        DO I = 1,NTENS
          STRESS(I) = STRESS_TRIAL(I)
        END DO
      END IF
C
C     存储状态变量
      IF (NSTATV .GE. 1) THEN
        STATEV(1) = F  ! 屈服函数值
      END IF
C
      RETURN
      END
