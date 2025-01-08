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
C     弹性参数
      E = PROPS(1)    ! 杨氏模量
      XNU = PROPS(2)  ! 泊松比
C
C     计算弹性常数
      EBULK3 = E/(1.0D0-2.0D0*XNU)
      EG2 = E/(1.0D0+XNU)
      EG = EG2/2.0D0
      ELAM = (EBULK3-EG2)/3.0D0
C
C     初始化刚度矩阵
      DO K1=1,NTENS
        DO K2=1,NTENS
          DDSDDE(K2,K1)=0.0D0
        END DO
      END DO
C
C     填充刚度矩阵
      DO K1=1,NDI
        DO K2=1,NDI
          DDSDDE(K2,K1)=ELAM
        END DO
        DDSDDE(K1,K1)=EG2+ELAM
      END DO
      DO K1=NDI+1,NTENS
        DDSDDE(K1,K1)=EG
      END DO
C
C     计算应力增量
      DO K1=1,NTENS
        DO K2=1,NTENS
          STRESS(K1)=STRESS(K1)+DDSDDE(K1,K2)*DSTRAN(K2)
        END DO
      END DO
C
      RETURN
      END
