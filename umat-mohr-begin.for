      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
C     声明局部变量
      DOUBLE PRECISION ZERO,ONE,TWO,THREE,SIX,PI
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,SIX=6.D0)
      DOUBLE PRECISION E,XNU,COHESION,PHI,PSI
      DOUBLE PRECISION GMODU,BULK,ELAM
      DOUBLE PRECISION SINPHI,COSPHI,SINPSI,TANPHI
      DOUBLE PRECISION STRESS_TR,F,ALPHA
      DOUBLE PRECISION PRINC(3),DPRINC(3,3)
C
C     材料参数
      E = PROPS(1)        ! 杨氏模量
      XNU = PROPS(2)      ! 泊松比
      COHESION = PROPS(3) ! 粘聚力
      PHI = PROPS(4)      ! 内摩擦角(度)
      PSI = PROPS(5)      ! 膨胀角(度)
C
C     将角度转换为弧度
      PI = 4.D0*ATAN(ONE)
      PHI = PHI*PI/180.D0
      PSI = PSI*PI/180.D0
C
C     计算弹性常数
      GMODU = E/(TWO*(ONE+XNU))
      BULK = E/(THREE*(ONE-TWO*XNU))
      ELAM = BULK - TWO*GMODU/THREE
C
C     计算MC模型参数
      SINPHI = SIN(PHI)
      COSPHI = COS(PHI)
      SINPSI = SIN(PSI)
      TANPHI = TAN(PHI)
      
C     初始化弹性刚度矩阵
      DO I=1,NTENS
        DO J=1,NTENS
          DDSDDE(I,J) = ZERO
        END DO
      END DO
C
C     填充弹性刚度矩阵
      DO I=1,NDI
        DO J=1,NDI
          DDSDDE(I,J) = ELAM
        END DO
        DDSDDE(I,I) = ELAM + TWO*GMODU
      END DO
C
      DO I=NDI+1,NTENS
        DDSDDE(I,I) = GMODU
      END DO
C
C     计算试探应力
      DO I=1,NTENS
        STRESS_TR = STRESS(I)
        DO J=1,NTENS
          STRESS_TR = STRESS_TR + DDSDDE(I,J)*DSTRAN(J)
        END DO
        STATEV(I) = STRESS_TR    ! 保存试探应力
      END DO
C
C     计算主应力
      CALL SPRIND(STATEV,PRINC,DPRINC,3,3,3)
      
C     排序主应力（σ1 ≥ σ2 ≥ σ3）
      CALL SORT_STRESS(PRINC)
      
C     计算屈服函数值
      F = PRINC(1) - PRINC(3)*SINPHI - TWO*COHESION*COSPHI
      
C     判断是否屈服
      IF (F .GT. ZERO) THEN
C       塑性修正
        ALPHA = F/(TWO*GMODU*(ONE+SINPHI*SINPSI))
        
C       更新主应力
        PRINC(1) = PRINC(1) - TWO*GMODU*ALPHA
        PRINC(3) = PRINC(3) + TWO*GMODU*ALPHA*SINPSI
        PRINC(2) = PRINC(2) - GMODU*ALPHA*(ONE+SINPSI)
        
C       将主应力转换回应力分量
        CALL TRANSFORM_STRESS_BACK(PRINC,STRESS,DPRINC)
        
C       更新状态变量
        STATEV(NTENS+1) = ONE  ! 塑性标志
        STATEV(NTENS+2) = F    ! 屈服函数值
      ELSE
C       弹性状态
        DO I=1,NTENS
          STRESS(I) = STATEV(I)
        END DO
        STATEV(NTENS+1) = ZERO
        STATEV(NTENS+2) = F
      ENDIF
C
      RETURN
      END
C
C     辅助子程序：对主应力排序
      SUBROUTINE SORT_STRESS(PRINC)
      IMPLICIT NONE
      DOUBLE PRECISION PRINC(3), TEMP
      INTEGER I, J
      
      DO I=1,2
        DO J=I+1,3
          IF (PRINC(I) .LT. PRINC(J)) THEN
            TEMP = PRINC(I)
            PRINC(I) = PRINC(J)
            PRINC(J) = TEMP
          ENDIF
        END DO
      END DO
      RETURN
      END
C
C     辅助子程序：应力转换
      SUBROUTINE TRANSFORM_STRESS_BACK(PRINC,STRESS,DPRINC)
      IMPLICIT NONE
      DOUBLE PRECISION PRINC(3),STRESS(6),DPRINC(3,3)
      DOUBLE PRECISION ZERO
      PARAMETER(ZERO=0.D0)
      INTEGER I,J
      
      DO I=1,3
        STRESS(I) = ZERO
        DO J=1,3
          STRESS(I) = STRESS(I) + PRINC(J)*DPRINC(I,J)**2
        END DO
      END DO
      
      IF (SIZE(STRESS) .GT. 3) THEN
        DO I=4,6
          STRESS(I) = ZERO
          J = I-3
          STRESS(I) = PRINC(1)*DPRINC(1,J)*DPRINC(1,J+1) +
     1                PRINC(2)*DPRINC(2,J)*DPRINC(2,J+1) +
     2                PRINC(3)*DPRINC(3,J)*DPRINC(3,J+1)
        END DO
      ENDIF
      
      RETURN
      END
