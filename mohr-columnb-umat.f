! 辅助函数
      FUNCTION ISNAN(X) RESULT(RES)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: X
      LOGICAL :: RES
      RES = (X .NE. X)
      END FUNCTION ISNAN

      FUNCTION TRACE(TENSOR) RESULT(TR)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: TENSOR(3,3)
      DOUBLE PRECISION :: TR
      TR = TENSOR(1,1) + TENSOR(2,2) + TENSOR(3,3)
      END FUNCTION TRACE

      FUNCTION IDENTITY() RESULT(I_MATRIX)
      IMPLICIT NONE
      DOUBLE PRECISION :: I_MATRIX(3,3)
      INTEGER :: I
      I_MATRIX = 0.0D0
      DO I = 1, 3
        I_MATRIX(I,I) = 1.0D0
      END DO
      END FUNCTION IDENTITY

      FUNCTION MAP_TO_TENSOR(STRESS_VECTOR) RESULT(TENSOR_FORM)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: STRESS_VECTOR(6)
      DOUBLE PRECISION :: TENSOR_FORM(3,3)
    
      TENSOR_FORM = 0.0D0  ! 初始化
    
    ! 正应力分量
      TENSOR_FORM(1,1) = STRESS_VECTOR(1)
      TENSOR_FORM(2,2) = STRESS_VECTOR(2)
      TENSOR_FORM(3,3) = STRESS_VECTOR(3)
    
    ! 剪应力分量
      TENSOR_FORM(1,2) = STRESS_VECTOR(4)
      TENSOR_FORM(2,1) = STRESS_VECTOR(4)
      TENSOR_FORM(1,3) = STRESS_VECTOR(5)
      TENSOR_FORM(3,1) = STRESS_VECTOR(5)
      TENSOR_FORM(2,3) = STRESS_VECTOR(6)
      TENSOR_FORM(3,2) = STRESS_VECTOR(6)
      END FUNCTION MAP_TO_TENSOR

      FUNCTION MAP_TO_VECTOR(STRESS_TENSOR) RESULT(VECTOR_FORM)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: STRESS_TENSOR(3,3)
      DOUBLE PRECISION :: VECTOR_FORM(6)
    
    ! 正应力分量
      VECTOR_FORM(1) = STRESS_TENSOR(1,1)
      VECTOR_FORM(2) = STRESS_TENSOR(2,2)
      VECTOR_FORM(3) = STRESS_TENSOR(3,3)
    
    ! 剪应力分量
      VECTOR_FORM(4) = STRESS_TENSOR(1,2)
      VECTOR_FORM(5) = STRESS_TENSOR(1,3)
      VECTOR_FORM(6) = STRESS_TENSOR(2,3)
      END FUNCTION MAP_TO_VECTOR

      SUBROUTINE HANDLE_ERROR(MSG)
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: MSG
      WRITE(*,*) 'ERROR in UMAT: ', MSG
      CALL XIT
	END SUBROUTINE HANDLE_ERROR
	SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     & RPL,DDSDDT,DRPLDE,DRPLDT,
     & STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     & NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     & CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,STRESSOLD)

      INCLUDE 'ABA_PARAM.INC'

    ! 声明接口
      INTERFACE
        FUNCTION ISNAN(X)
            DOUBLE PRECISION, INTENT(IN) :: X
            LOGICAL :: ISNAN
        END FUNCTION ISNAN

        FUNCTION TRACE(TENSOR)
            DOUBLE PRECISION, INTENT(IN) :: TENSOR(3,3)
            DOUBLE PRECISION :: TRACE
        END FUNCTION TRACE

        FUNCTION IDENTITY()
            DOUBLE PRECISION :: IDENTITY(3,3)
        END FUNCTION IDENTITY

        FUNCTION MAP_TO_TENSOR(STRESS_VECTOR)
            DOUBLE PRECISION, INTENT(IN) :: STRESS_VECTOR(6)
            DOUBLE PRECISION :: MAP_TO_TENSOR(3,3)
        END FUNCTION MAP_TO_TENSOR

        FUNCTION MAP_TO_VECTOR(STRESS_TENSOR)
            DOUBLE PRECISION, INTENT(IN) :: STRESS_TENSOR(3,3)
            DOUBLE PRECISION :: MAP_TO_VECTOR(6)
        END FUNCTION MAP_TO_VECTOR
	END INTERFACE
    ! ABAQUS 参数声明
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     &  DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     &  STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     &  PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

      ! 添加调试输出
      WRITE(*,*) 'Debug: Starting UMAT, Element=', NOEL, 
     & ' IntPt=', NPT

      ! 声明局部变量
      INTEGER I, J, K
      DOUBLE PRECISION E, XNU, EBULK3, EG2, EG, ELAM
      DOUBLE PRECISION SIG1, SIG2, SIG3
      DOUBLE PRECISION SMEAN, PHI, PSI, COHESION
      DOUBLE PRECISION SYIELD, HARD, DYDP
      DOUBLE PRECISION SIGMEAN, DSIGMEAN
      DOUBLE PRECISION SDEVIA(6), DDEVIA(6)
      DOUBLE PRECISION DPSTRAN(6), DSSTRAN(6)
      DOUBLE PRECISION BETA1, BETA2, DGAMMA
      DOUBLE PRECISION F1, F2, DENOM
      DOUBLE PRECISION STRESS_NEW(6)
      DOUBLE PRECISION DSTRESS(6)
      DOUBLE PRECISION PRINCIPAL(3)
      DOUBLE PRECISION DIRECTION(3,3)
	 ! 添加调试输出
      WRITE(*,*) 'Debug: Starting UMAT, Element=', NOEL, 
     & ' IntPt=', NPT

      ! 安全检查
      IF (NTENS .NE. 6) THEN
          WRITE(*,*) 'Error: NTENS must be 6 for C3D8R element'
          CALL XIT
      ENDIF

      IF (NPROPS .LT. 2) THEN
          WRITE(*,*) 'Error: Not enough material properties'
          CALL XIT
      ENDIF
       ! 声明您的常量
      DOUBLE PRECISION SQR2, TOL
      DOUBLE PRECISION STRESS_TOL, YIELD_TOL, DEV_TOL
      DOUBLE PRECISION STIFF_SCALE, MAX_GAMMA, MAX_STRESS_RATIO    
	!定义参数值
	PARAMETER (SQR2 = 1.41421356237D0)
	PARAMETER (TOL = 1.0D-20)
	PARAMETER (STRESS_TOL = 1.0D-8)
	PARAMETER (YIELD_TOL = 1.0D-8)
	PARAMETER (DEV_TOL = 1.0D-8)
	PARAMETER (STIFF_SCALE = 1.0D-3)
	PARAMETER (MAX_GAMMA = 0.01D0)
	PARAMETER (MAX_STRESS_RATIO = 100.0D0)
	WRITE(*,*) 'Debug: Parameters defined successfully'

    ! 收敛控制变量
      DOUBLE PRECISION :: STRESS_NORM, PREV_STRESS_NORM
      INTEGER :: ITER_COUNT, MAX_ITER
      LOGICAL :: CONVERGED

    ! 初始化检查
      IF (NTENS .LT. 6) THEN
        CALL HANDLE_ERROR('NTENS must be at least 6')
      END IF

      IF (NSTATV .LT. 4) THEN
        CALL HANDLE_ERROR('NSTATV must be at least 4')
      END IF

    ! 初始化数组
      DDSDDE = 0.0D0
      DEPS_E = 0.0D0
      STR_TRIAL = 0.0D0
      SIG_TRIAL = 0.0D0
      SIG_CORR = 0.0D0
      DEVIATORIC_STRESS = 0.0D0
      TEMP_MATRIX = 0.0D0
      I_MATRIX = 0.0D0

    ! 设置单位矩阵
      I_MATRIX = IDENTITY()
	
    ! 读取材料参数
      K = PROPS(1)       ! 体积模量
      G = PROPS(2)       ! 剪切模量
      PHI = PROPS(3)     ! 摩尔-库伦摩擦角 (弧度)
      PSI = PROPS(4)     ! 脆性角 (弧度)
      C = PROPS(5)       ! 黏聚力

    ! 材料参数验证
      IF (K .LE. 0.0D0 .OR. G .LE. 0.0D0 .OR. C .LE. 0.0D0) THEN
        CALL HANDLE_ERROR('Invalid material parameters')
      END IF

    ! 初始化状态变量
      IF (KSTEP .EQ. 1 .AND. KINC .EQ. 1) THEN
        STATEV(1) = 0.0D0  ! 塑性乘子
        STATEV(2) = 0.0D0  ! 屈服函数值
        STATEV(3) = 0.0D0  ! 偏应力范数
        STATEV(4) = 0.0D0  ! 塑性状态标志
      END IF

    ! 初始化收敛控制变量
      ITER_COUNT = 0
      MAX_ITER = 50
      CONVERGED = .FALSE.

    ! 计算弹性预测应力
      DO I = 1, NTENS
        DEPS_E(I) = DSTRAN(I)
        STR_TRIAL(I) = STRESS(I)
        DO J = 1, NTENS
            STR_TRIAL(I) = STR_TRIAL(I) + DDSDDE(I,J) * DEPS_E(J)
        END DO
      END DO

    ! 将向量形式转换为张量形式
      SIG_TRIAL = MAP_TO_TENSOR(STR_TRIAL)

    ! 计算偏应力
      DEVIATORIC_STRESS = SIG_TRIAL - (TRACE(SIG_TRIAL)/3.0D0) * I_MATRIX
      NORM_DEVIATORIC = 0.0D0
      DO I = 1, 3
        DO J = 1, 3
            NORM_DEVIATORIC = NORM_DEVIATORIC + 
     &          DEVIATORIC_STRESS(I,J) * DEVIATORIC_STRESS(I,J)
        END DO
      END DO
      NORM_DEVIATORIC = SQRT(0.5D0 * NORM_DEVIATORIC)

    ! 计算摩尔库仑参数
	M_COULOMB = 6.0D0 * SIN(PHI)/(SQR2 * (3.0D0 - SIN(PHI)))
      M = 6.0D0 * SIN(PSI)/(SQR2 * (3.0D0 - SIN(PSI)))

    ! 计算屈服函数值
      F_TRIAL = NORM_DEVIATORIC + 
     &    (M_COULOMB/3.0D0) * TRACE(SIG_TRIAL) - C * COS(PHI)

    ! 检查屈服条件
      IF (F_TRIAL .LE. YIELD_TOL) THEN
        ! 弹性状态
        STATEV(4) = 0.0D0
        SIG_CORR = SIG_TRIAL
      ELSE
        ! 塑性状态
        STATEV(4) = 1.0D0

        ! 计算塑性乘子
        GAMMA = F_TRIAL / (2.0D0 * G + (K * M * M_COULOMB)/3.0D0)
        GAMMA = MIN(GAMMA, MAX_GAMMA)  ! 限制塑性修正量

        ! 应力修正
        SIG_CORR = SIG_TRIAL - 2.0D0 * G * GAMMA * 
     &    (DEVIATORIC_STRESS/(2.0D0 * NORM_DEVIATORIC)) - 
     &    K * M * GAMMA * I_MATRIX

        ! 检查应力修正的合理性
        DO I = 1, 3
            DO J = 1, 3
                IF (ABS(SIG_CORR(I,J)) .GT. 
     &              MAX_STRESS_RATIO * ABS(SIG_TRIAL(I,J))) THEN
                    WRITE(*,*) 'WARNING: Excessive stress correction'
                    PNEWDT = 0.25D0
                    RETURN
                END IF
            END DO
        END DO

        ! 更新状态变量
        STATEV(1) = GAMMA
        STATEV(2) = F_TRIAL
        STATEV(3) = NORM_DEVIATORIC
      END IF
    ! 计算切线刚度矩阵
      IF (STATEV(4) .LT. 0.5D0) THEN
        ! 弹性切线刚度
        DO I = 1, NTENS
            DO J = 1, NTENS
                IF (I .EQ. J) THEN
                    IF (I .LE. NDI) THEN
                        DDSDDE(I,J) = K + 4.0D0*G/3.0D0
                    ELSE
                        DDSDDE(I,J) = G
                    END IF
                ELSE IF (I .LE. NDI .AND. J .LE. NDI) THEN
                    DDSDDE(I,J) = K - 2.0D0*G/3.0D0
                ELSE
                    DDSDDE(I,J) = 0.0D0
                END IF
            END DO
        END DO
      ELSE
        ! 弹塑性切线刚度
        DO I = 1, NTENS
            DO J = 1, NTENS
                DDSDDE(I,J) = DDSDDE(I,J) * STIFF_SCALE
            END DO
        END DO
      END IF

    ! 将修正后的应力转换回向量形式
      STRESS = MAP_TO_VECTOR(SIG_CORR)

    ! 检查数值稳定性
      HAS_NAN = .FALSE.
      DO I = 1, NTENS
        IF (ISNAN(STRESS(I))) THEN
            HAS_NAN = .TRUE.
            EXIT
        END IF
      END DO

      IF (HAS_NAN) THEN
        WRITE(*,*) 'WARNING: NaN detected in stress'
        PNEWDT = 0.25D0
        RETURN
      END IF

    ! 计算应力变化率用于时间步长控制
      STRESS_DOT = 0.0D0
      DO I = 1, NTENS
        STRESS_DOT = STRESS_DOT + 
     &    (STRESS(I) - STRESSOLD(I))**2
      END DO
      STRESS_DOT = SQRT(STRESS_DOT)/DTIME

    ! 根据应力变化率调整时间步长
      IF (STRESS_DOT .GT. 1.0D3) THEN
        PNEWDT = 0.5D0
      ELSE IF (STRESS_DOT .LT. 1.0D1) THEN
        PNEWDT = 1.5D0
      END IF

      RETURN
      END SUBROUTINE UMAT

