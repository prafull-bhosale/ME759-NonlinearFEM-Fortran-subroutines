!11111111122222222223333333333444444444455555555556666666666777777777777
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                      C
!     UEL9_ECLA.for    9-Noded Classical Elasticity                    C
!     GENERAL FINITE ELEMENT LIBRARY TO BE USED WITHIN ABAQUS          C
!     UNIVERSIDAD EAFIT                                                C
!     LABORATORIO DE MECANICA APLICADA                                 C
!     BLOQUE 14-PISO 2                                                 C
!     MEDELLIN, COLOMBIA                                               C
!                                                                      C
!     LAST UPDATED FEB 21/2022                                         C
!                                                                      C
!     CURRENTLY:                                                       C
!                                                                      C
!     4 NODED ISOPARAMETRIC PURE DISPLACEMENT ELEMENT                  C
!                                                                      C
!----------------------------------------------------------------------C
!                                                                      C
!----------------------------------------------------------------------C
!23456789012345678901234567890123456789012345678901234567890123456789012
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!11111111122222222223333333333444444444455555555556666666666777777777777
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                      C
!       --U S E R   E L E M E N T    S U B R O U T I N E S---          C
!                                                                      C
!23456789012345678901234567890123456789012345678901234567890123456789012
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                      C
!  USER ELEMENT SUBROUTINE- ELASTIC MATERIAL BEHAVIOR                  C
!                                                                      C
!  4 NODED ISOPARAMETRIC plain strain ELEMENT                                 C
!  2X2 GAUSS INTEGRATION                                               C
!  EDITED BY Prafull Bhosale for 4 noded quadrilateral plain strain element  C
!  CREATED BY JUAN GOMEZ                                               C
!                                                                      C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
!     
      INCLUDE 'ABA_PARAM.INC'
!
      CHARACTER*80 CMNAME
!  
      PARAMETER (ZERO=0.D0,HALF=0.5D0,ONE=1.D0,NTENS=4,TWO=2.D0,THREE=3.D0)
!
!     Parameters required by UMAT.f
!
      PARAMETER (NDI=3,NSTATV=14,SSE=0.D0,SCD=0.D0,RPL=0.D0,DRPLDT=0.D0,TEMP=0.D0,DTEMP=0.D0,NSHR=1,CELENT=2.D0,LAYER=1,KSPT=1)
!
!     Parameter arrays from UEL.f
!
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
!
!     User defined arrays
!
      DIMENSION B(NTENS,NDOFEL),BT(NDOFEL,NTENS),BLM(1,NDOFEL),BLMT(NDOFEL,1),FRST1(NDOFEL,NDOFEL),FRST2(NDOFEL,1),XX(2,NNODE),XW(9),XP(2,9),AUX1(NTENS,NDOFEL),STRESS(NTENS),STRAN(NTENS),DSTRAN(NTENS),EELAS(NTENS),EPLAS(NTENS),SR(NDOFEL)
!
!     Arrays to be used in UMAT.f
!
      DIMENSION STATEV(NSTATV),DDSDDT(NTENS),DRPLDE(NTENS),UPREDEF(1), DPRED(1),UCOORDS(3),DROT(3, 3), DFGRD0(3, 3),DFGRD1(3, 3),DDSDDE(NTENS, NTENS)
!
!     Initializes parameters required by UMAT.f
!
      CALL KCLEARV(SR,NDOFEL)
      CALL KCLEARV(STRESS,NTENS)
      CALL KCLEARV(STRAN,NTENS)
      CALL KCLEARV(DSTRAN,NTENS)
      CALL KCLEARV(DDSDDT,NTENS)
      CALL KCLEARV(DRPLDE,NTENS)
      CALL KCLEAR(DROT,3,3)
      CALL KCLEAR(DFGRD0,3,3)
      CALL KCLEAR(DFGRD1,3,3)
      CALL KCLEAR(DDSDDE,NTENS,NTENS)
      UPREDEF(1)=ZERO
      DPRED(1)=ZERO
!
!     Clears RHS vector and Stiffness matrix
!
      CALL KCLEARV(RHS,NDOFEL)
      CALL KCLEAR(AMATRX,NDOFEL,NDOFEL)
      CALL KCLEAR(AUX1,NTENS,NDOFEL)
!
      CALL KGPOINTS2X2(XP,XW)
      NGPT=4
!
!     Loops around all Gauss points
!
      DO NN=1,NGPT
!
        RII=XP(1,NN)
        SII=XP(2,NN)
        ALF=XW(NN)
!
!       Compute State variable index corresponding
!       to current Gauss point and load stress,total strain
!       elastic strain,plastic strain and equivalent plastic
!       strain from state variables as defined in USER ELEMENT. 
!       Different variables are required for different constitutive models.
!
        ISVINT=1+(NN-1)*NSVARS/NGPT
        JJ=1
        DO II=ISVINT,ISVINT+3
          STRESS(JJ)=SVARS(II)
          STRAN(JJ )=SVARS(II+4)
          JJ=JJ+1
        END DO
!
!       Starts STATEV array definition as required by UMAT.f
!       Different variables for different constitutive models.
!
!
!       Ends STATEV array definition as required by UMAT.f
!
!       Assembles strain-displacement and relative rotation-displacement
!       matrices B and BLM=(L-M)
!
        CALL KSTDM(JELEM,NNODE,NDOFEL,NTENS,COORDS,B,DDET,RII,SII,XBAR)
!
!       Computes strain increment and updates strain.
!
        CALL KMAVEC(B,NTENS,NDOFEL,DU,DSTRAN)
        CALL KUPDVEC(STRAN,NTENS,DSTRAN)
!
!       Assembles material matrix and updates state variables
!
        CALL UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP,DTEMP, UPREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV,PROPS, NPROPS, UCOORDS, DROT, PNEWDT, CELENT, DFGRD0,DFGRD1, JELEM, NN, LAYER, KSPT, KSTEP, KINC)
!
!       Assembles STIFFNESS matrix
!                                        T
!       Generalized stress component.===B M B (stress and couple stress)
!
        CALL KMMULT(DDSDDE,NTENS,NTENS,B,NTENS,NDOFEL,AUX1)
        CALL KMTRAN(B,NTENS,NDOFEL,BT)
        CALL KMMULT(BT,NDOFEL,NTENS,AUX1,NTENS,NDOFEL,FRST1)
!
!       Assembles RHS vector contribution
!
        CALL KMAVEC(BT,NDOFEL,NTENS,STRESS,FRST2)
!
!       Considers Gauss weight and Jacobian determinant representing
!       volume differential.
!
        CALL KSMULT(FRST1,NDOFEL,NDOFEL,ALF*DDET*XBAR)
        CALL KSMULT(FRST2,NDOFEL,1,-ALF*DDET*XBAR)
!
!       Updates Stiffness matrix and RHS vector
!       
        CALL KUPDMAT(AMATRX,NDOFEL,NDOFEL,FRST1)
        CALL KUPDVEC(RHS,NDOFEL,FRST2)
!
!       Clears material Jacobian and temporary stiffness matrix
!       array for new Gauss point 
!
        CALL KCLEAR(DDSDDE,NTENS,NTENS)
        CALL KCLEAR(FRST1,NDOFEL,NDOFEL)
        CALL KCLEAR(FRST2,NDOFEL,1)
!
!       Starts updating of state variables with updated values from UMAT.f
!
        JJ=1
        DO II=ISVINT,ISVINT+3
          SVARS(II)=STRESS(JJ)
          SVARS(II+4 )=STRAN(JJ)
          JJ=JJ+1
        END DO
!
!       Ends updating of state variables with updated values from UMAT.f
!
      END DO
!
!     Assembles surface tractions contribution.
!
      DO II=1,NDLOAD
        IDFACE=JDLTYP(II,1)
        PMAG=ADLMAG(II,1)
        CALL KSURTRAC(IDFACE,PMAG,NDOFEL,NNODE,COORDS,SR)
        CALL KUPDVEC(RHS,NDOFEL,SR)
      END DO

!
!     EXTRAPOLATE STRAIN TO THE NODES
!
!      RII=-1.0
!      SII=1.0
!      CALL KSTDM(JELEM,JTYPE,NNODE,NDOFEL,NTENS,COORDS,B,BD,BV,HP,
!     1           DDET,RII,SII,XBAR)
!      CALL KMAVEC(B,NTENS,NDOFEL,U,STRAN)
!      svars(109)=stress(2)
!      SVARS(110)=STRAN(2)
!cc
      RETURN
!      
      END
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!11111111122222222223333333333444444444455555555556666666666777777777777
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                      C
!--C O N S T I T U T I V E  M A T E R I A L S  S U B R O U T I N E S---C
!                                                                      C
!23456789012345678901234567890123456789012345678901234567890123456789012
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!11111111122222222223333333333444444444455555555556666666666777777777777
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                      C
!     SUBROUTINE UMAT                                                  C
!     ISOTROPIC ELASTICITY                                             C
!     CANNOT BE USED FOR PLANE STRESS                                  C
!     NTENS: LENGTH OF STRESS VECTOR                                   C
!     NDI:   NUMBER OF NORMAL STRESS COMPONENTS                        C
!     JANUARY 26/2004                                                  C
!                                                                      C
!     LOCAL ARRAYS                                                     C
!                                                                      C
!     PROPS(1) - E                                                     C
!     PROPS(2) - NU                                                    C
!                                                                      C
!     DDSDDE() - MATERIAL JACOBIAN                                     C
!     STRESS() - UPDATED STRESS VECTOR                                 C
!                                                                      C
!23456789012345678901234567890123456789012345678901234567890123456789012
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP,DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV,PROPS, NPROPS, COORDS, DROT, PNEWDT, CELENT, DFGRD0,DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)
!
      INCLUDE 'ABA_PARAM.INC'
!
      CHARACTER*80 CMNAME
! 
      DIMENSION STRESS(NTENS), STATEV(NSTATV),DDSDDE(NTENS, NTENS),DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1), DPRED(1), PROPS(NPROPS),COORDS(3),DROT(3, 3),DFGRD0(3, 3), DFGRD1(3, 3)
!
      DIMENSION DS(NTENS)
!
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.0D0,ENUMAX=0.4999D0, NEWTON=10, TOLER=1.0D-7,MAXITER=30,FOUR=4.D0)
!
!**********************************************************************
!               C O S S E R A T  E L A S T I C I T Y
!**********************************************************************
!
!     Elastic properties
!
      EMOD=PROPS(1)
      ENU=MIN(PROPS(2),ENUMAX)
      PLS=PROPS(3)
      EBULK3=EMOD/(ONE-TWO*ENU)
      EG2=EMOD/(ONE+ENU)
      EG=EG2/TWO
      EG3=THREE*EG
      ELAM=(EBULK3-EG2)/THREE
!
!     Elastic stiffness
!
      DO K1=1, 3
        DO K2=1, 3
          DDSDDE(K2, K1)=ELAM
        END DO
        DDSDDE(K1, K1)=EG2+ELAM
      END DO
      DO K1=4, 4
        DDSDDE(K1, K1)=EG
      END DO
!
!     Calculate the stress increment and
!     updates the stress vector.
!
      CALL KMAVEC(DDSDDE,NTENS,NTENS,DSTRAN,DS)
      CALL KUPDVEC(STRESS,NTENS,DS)
!
      RETURN
!
      END
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!11111111122222222223333333333444444444455555555556666666666777777777777
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                      C
!         --G E N E R A L  F E M  S U B R O U T I N E S--              C
!                                                                      C
!23456789012345678901234567890123456789012345678901234567890123456789012
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
!     SUBROUTINE KSTDM:GENERATES THE STRAIN-DISPLACEMENT MATRIX B      C
!     AND JACOBIAN DETERMINANT DDET AT THE POINT r ,s                  C
!                                                 i  j                 C
!     FOR AN 9-NODED 2D ELEMENT-PLANE STRAIN                           C
!     B    =STRAIN-DISPLACEMENT MATRIX                                 C
!     BLM  =RELATIVE STRAIN DISPLACEMENT MATRIX                        C
!     DDET =JACOBIAN DETERMINANT                                       C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!
      SUBROUTINE KSTDM(IDEL,NNE,NDOFEL,NTENS,XX,B,DDET,R,S,XBAR)
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=0.5D0)
!
      DIMENSION XX(2,NNE),B(NTENS,NDOFEL),P(2,NNE),AUX1(2,NNE),XJ(2,2),XJI(2,2)
!
!     Initialize arrays
!
      XBAR=ONE
      CALL KCLEAR(B,NTENS,NDOFEL)
      CALL KCLEAR(XJ,2,2)
      CALL KCLEAR(XJI,2,2)
      CALL KCLEAR(P,2,NNE)
!
!     Shape functions derivatives w.r.t natural coordinates
!
      CALL KSFDER(IELT,NDOFEL,NNE,R,S,P)
!
!     Computes the Jacobian operator and its determinant at point (r,s)
!
      CALL KJACOPER(NNE,XJ,XX,P)
      DDET=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1)
!
!     Computes the inverse of the Jacobiam operator
!
      CALL KJACINVE(XJ,DDET,XJI)
!
!     Jacobian Inverse times Shape Functions derivatives w.r.t natural coordinates
!     to produce shape function derivatives w.r.t x,y coordinates.
!          
      CALL KMMULT(XJI,2,2,P,2,NNE,AUX1)
!
!     Assembles B matrix for a
!     Cosserat element.
!
      DO I=1,NNE
        II=2*(I-1)+1
        B(1,II)=AUX1(1,I)
        B(2,II+1)=AUX1(2,I)
        B(4,II)=AUX1(2,I)
        B(4,II+1)=AUX1(1,I)
      END DO
!
      RETURN
!
      END
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
!     SUBROUTINE KSFDER:GENERATES THE SHAPE FUNCTION DERIVATIVES       C
!     ACCORDING TO ELEMENT TYPE AT THE POINT r ,s                      C
!                                             i  j                     C
!     B    =STRAIN-DISPLACEMENT MATRIX                                 C
!     DDET =JACOBIAN DETERMINANT                                       C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!
      SUBROUTINE KSFDER(IELT,NDOFEL,NNE,R,S,P)
!
      IMPLICIT REAL*8(A-H,O-Z)
!
      PARAMETER(ONE=1.D0,TWO=2.D0, HALF=0.5D0,QUARTER=0.25D0,FOUR=4.D0)
!
      DIMENSION P(2,NNE)
!
!
!     4-NODED ELEMENT
!     Derivatives w.r.t the natural coordinates
!     w.r.t.r
!
      P(1,1)=(1./4.)*(S-1)
      P(1,2)=(1./4.)*(1-S)
      P(1,3)=(1./4.)*(1+S)
      P(1,4)=(1./4.)*(-S-1)
!
!     w.r.t.s
!
      P(2,1)=(1./4.)*(R-1)
      P(2,2)=(1./4.)*(-R-1)
      P(2,3)=(1./4.)*(1+R)
      P(2,4)=(1./4.)*(1-R)
!
      RETURN
!
      END
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                      C
! 3X3 GAUSS POINTS GENERATION                                          C
!                                                                      C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KGPOINTS3X3(XP,XW)
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      PARAMETER (ZERO=0.D0,ONE=1.D0)
!
      DIMENSION XP(2,9),XW(9),RLGP(2,9)
!
      RLGP(1,1)=-ONE
      RLGP(1,2)=ZERO
      RLGP(1,3)= ONE
      RLGP(1,4)=-ONE
      RLGP(1,5)=ZERO
      RLGP(1,6)= ONE
      RLGP(1,7)=-ONE
      RLGP(1,8)=ZERO
      RLGP(1,9)= ONE
!
      RLGP(2,1)=-ONE
      RLGP(2,2)=-ONE
      RLGP(2,3)=-ONE
      RLGP(2,4)=ZERO
      RLGP(2,5)=ZERO
      RLGP(2,6)=ZERO
      RLGP(2,7)= ONE
      RLGP(2,8)= ONE
      RLGP(2,9)= ONE
!
      XW(1)=0.555555555555556D0**2
      XW(2)=0.555555555555556D0*0.888888888888889D0
      XW(3)=0.555555555555556D0**2
!
      XW(4)=0.555555555555556D0*0.888888888888889D0
      XW(5)=0.888888888888889D0**2
      XW(6)=0.555555555555556D0*0.888888888888889D0
!
      XW(7)=0.555555555555556D0**2
      XW(8)=0.555555555555556D0*0.888888888888889D0
      XW(9)=0.555555555555556D0**2
!
      G=DSQRT(0.60D0)
      DO I=1,2
        DO J=1,9 
          XP(I,J)=G*RLGP(I,J)
        END DO
      END DO
!     
      RETURN
!
      END
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
!   SUBROUTINE KGPOINTS2X2                                             C
!                                                                      C
!   2X2 GAUSS POINTS GENERATION                                        C
!                                                                      C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!
      SUBROUTINE KGPOINTS2X2(XP,XW)
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      PARAMETER (ZERO=0.D0,ONE=1.D0)
!
      DIMENSION XP(2,4),XW(4)
!
      CALL KCLEARV(XW,4)
      CALL KCLEAR(XP,2,4)
!
      DO K1=1,4
        XW(K1)=ONE
      END DO
!
      XP(1,1)=-0.577350269189626
      XP(2,1)=-0.577350269189626
      XP(1,2)= 0.577350269189626
      XP(2,2)=-0.577350269189626
      XP(1,3)=-0.577350269189626
      XP(2,3)= 0.577350269189626
      XP(1,4)= 0.577350269189626
      XP(2,4)= 0.577350269189626
!
      RETURN
!
      END
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
!   SUBROUTINE KJACOPER                                                C
!                                                                      C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!
      SUBROUTINE KJACOPER(NNE,XJA,XCORD,RDER)
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      PARAMETER(ZERO=0.0D0)
!
      DIMENSION XJA(2,2),XCORD(2,NNE),RDER(2,NNE)
!
      CALL KCLEAR(XJA,2,2)
!
      DUM=ZERO
      DO K=1,2
        DO J=1,2
          DO I=1,NNE
            DUM=DUM+RDER(K,I)*XCORD(J,I)
          END DO
          XJA(K,J)=DUM
          DUM=ZERO
        END DO
      END DO
!
      RETURN
!
      END
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                      C
!      2X2 Jacobian inverse                                            C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!
      SUBROUTINE KJACINVE(XJA,DD,XJAI)
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      DIMENSION XJA(2,2),XJAI(2,2),XADJ(2,2)
!
      XADJ(1,1)=XJA(2,2)
      XADJ(1,2)=XJA(1,2)
      XADJ(2,1)=XJA(2,1)
      XADJ(2,2)=XJA(1,1)
      DO J=1,2
        DO K=1,2
          COFA=((-1)**(J+K))*XADJ(J,K)
          XJAI(J,K)=COFA/DD
        END DO
      END DO
!
      RETURN
!
      END
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
!             M A T R I X   H A N D L I N G                            C
!-------------U T I L I T I E S   B L O C K--------------              C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                      C
!      SUBROUTINE KCLEAR(A,N,M)                                        C
!      Clear a real matrix                                             C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KCLEAR(A,N,M)
!
      IMPLICIT REAL*8(A-H,O-Z)
!
      PARAMETER(ZERO=0.0D0)
      DIMENSION A(N,M)
!
      DO I=1,N
        DO J=1,M
          A(I,J)=ZERO
        END DO
      END DO
!
      RETURN
!
      END
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                      C
!      SUBROUTINE KMMULT(A,NRA,NCA,B,NRB,NCB,C)                        C
!      Real matrix product                                             C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KMMULT(A,NRA,NCA,B,NRB,NCB,C)
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      PARAMETER(ZERO=0.D0)
      DIMENSION A(NRA,NCA),B(NRB,NCB),C(NRA,NCB)
!
      CALL KCLEAR(C,NRA,NCB)
      DUM=ZERO
      DO I=1,NRA
        DO J=1,NCB
         DO K=1,NCA
           DUM=DUM+A(I,K)*B(K,J)
          END DO
          C(I,J)=DUM
          DUM=ZERO
        END DO
      END DO
!
      RETURN
!
      END
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                      C
!      SUBROUTINE KSMULT(A,NR,NC,S)                                    C
!      Matrix times a scalar.                                          C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KSMULT(A,NR,NC,S)
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      DIMENSION A(NR,NC)
!
      DO I=1,NR
        DO J=1,NC
          DUM=A(I,J)
          A(I,J)=S*DUM
          DUM=0.D0
        END DO  
      END DO
!
      RETURN
!
      END
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                      C
!      SUBROUTINE KUPDMAT(A,NR,NC,B)                                   C
!      Updates an existing matrix with an incremental matrix.          C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KUPDMAT(A,NR,NC,B)
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      PARAMETER(ZERO=0.D0)
!
      DIMENSION A(NR,NC),B(NR,NC)
!
      DO I=1,NR
        DO J=1,NC
          DUM=A(I,J)
          A(I,J)=ZERO
          A(I,J)=DUM+B(I,J)
          DUM=ZERO
        END DO
      END DO
!
      RETURN
!
      END
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                      C
!      SUBROUTINE KMTRAN(A,NRA,NCA,B)                                  C      
!      Matrix transpose                                                C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KMTRAN(A,NRA,NCA,B)
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      DIMENSION A(NRA,NCA),B(NCA,NRA)
!
      CALL KCLEAR(B,NCA,NRA)
      DO I=1,NRA
       DO J=1,NCA
         B(J,I)=A(I,J)
        END DO
      END DO
!
      RETURN
!
      END
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                      C
!      SUBROUTINE KMAVEC(A,NRA,NCA,B,C)                                C
!      Real matrix times vector                                        C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KMAVEC(A,NRA,NCA,B,C)
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      PARAMETER(ZERO=0.D0)
      DIMENSION A(NRA,NCA),B(NCA),C(NRA)
!
      CALL KCLEARV(C,NRA)
!
      DO K1=1,NRA
        DO K2=1,NCA
          C(K1)=C(K1)+A(K1,K2)*B(K2)	    
        END DO
      END DO     
!
      RETURN
!
      END
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                      C
!      SUBROUTINE KCLEARV(A,N)                                         C
!      Clear a real vector                                             C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KCLEARV(A,N)
!
      IMPLICIT REAL*8(A-H,O-Z)
!
      PARAMETER(ZERO=0.0D0)
!
      DIMENSION A(N)
!
      DO I=1,N
        A(I)=ZERO
      END DO
!
      RETURN
!
      END
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                      C
!      SUBROUTINE KUPDVEC(A,NR,B)                                      C
!      Updates an existing vector with an incremental vector.          C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KUPDVEC(A,NR,B)
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      PARAMETER(ZERO=0.D0)
!
      DIMENSION A(NR),B(NR)
!
      DO I=1,NR
        DUM=A(I)
        A(I)=ZERO
        A(I)=DUM+B(I)
        DUM=ZERO
      END DO
!
      RETURN
!
      END
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                      C
!      SUBROUTINE KVECSUB(A,NRA,B,NRB,C)                               C
!      Substracts one column vector from another column vector         C
!      IFLAG=0 for substraction                                        C
!      IFLAG=1 for addition                                            C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KVECSUB(A,NRA,B,NRB,C,IFLAG)
!
      IMPLICIT REAL*8(A-H,O-Z)
!
      PARAMETER (ONE=1.0D0, ONENEG=-1.0D0)
!
      DIMENSION A(NRA,1),B(NRB,1),C(NRB,1)
!
      SCALAR=ONENEG
!
      IF (IFLAG.EQ.1) SCALAR=ONE
!
      DO I=1,NRA
        C(I,1)=A(I,1)+B(I,1)*SCALAR
      END DO
!
      RETURN
!
      END
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                      C
!      SUBROUTINE KMATSUB(A,NRA,NCA,B,C,IFLAG)                         C
!      Substracts one rectangular matrix from another rectangular      C
!      matrix                                                          C
!      IFLAG=0 for substraction                                        C
!      IFLAG=1 for addition                                            C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KMATSUB(A,NRA,NCA,B,C,IFLAG)
!
      IMPLICIT REAL*8(A-H,O-Z)
!
      PARAMETER (ONE=1.0D0, ONENEG=-1.0D0)
!
      DIMENSION A(NRA,NCA),B(NRA,NCA),C(NRA,NCA)
!
      CALL KCLEAR(C,NRA,NCA)
!
      SCALAR=ONENEG
!
      IF (IFLAG.EQ.1) SCALAR=ONE
!
      DO I=1,NRA
        DO J=1,NCA
          C(I,J)=A(I,J)+B(I,J)*SCALAR
        END DO
      END DO
!
      RETURN
!
      END
!
!
!11111111122222222223333333333444444444455555555556666666666777777777777
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                      C
!     SUBROUTINE IDENTITY                                              C
!     CREATES AN IDENTITY MATRIX OF DIMENSIONS NDIM,NDIM               C
!                                                                      C
!23456789012345678901234567890123456789012345678901234567890123456789012
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      SUBROUTINE KIDENTITY(DEL,NDIM)
!
      IMPLICIT REAL*8(A-H,O-Z)
!
      PARAMETER(ONE=1.D0)
!
      DIMENSION DEL(NDIM,NDIM)
!
      CALL KCLEAR(DEL,NDIM,NDIM)
!
      DO K1=1,NDIM
        DEL(K1,K1)=ONE
      END DO
!
      RETURN
!
      END
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
!   SUBROUTINE KINVERSE                                                 C
!                                                                      C
!   IVEERSE OF A MATRIX USING LU DECOMPOSITION                         C
!   TAKEN FROM NUMERICAL RECIPES By Press et al                        C
!                                                                      C
!   A   Matrix to be inverted.                                         C
!   Y   Inverse of A                                                   C
!   N   Dimension                                                      C
!                                                                      C
!                                                                      C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!
      SUBROUTINE KINVERSE(A,Y,NP,N)
!
      IMPLICIT REAL*8(A-H,O-Z)
!
      PARAMETER (ZERO=0.D0,ONE=1.D0)
!
      DIMENSION A(NP,NP),Y(NP,NP),INDX(NP),AUX(NP,NP)
!
      CALL KCLEAR(AUX,NP,NP)
      CALL KCOPYMAT(A,AUX,N)
!
      DO I=1,N
        DO J=1,N
          Y(I,J)=ZERO
        END DO
        Y(I,I)=ONE
      END DO
      CALL KLUDCMP(AUX,N,NP,INDX,D)
      DO J=1,N
        CALL KLUBKSB(AUX,N,NP,INDX,Y(1,J))
      END DO
!
      RETURN
!
      END
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
!   SUBROUTINE KLUDCMP                                                 C
!                                                                      C
!   LU MATRIX DECOMPOSITION                                            C
!   TAKEN FROM NUMERICAL RECIPES By Press et al                        C
!                                                                      C
!                                                                      C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!
      SUBROUTINE KLUDCMP(A,N,NP,INDX,D)
!
      IMPLICIT REAL*8(A-H,O-Z)
!
      PARAMETER(NMAX=500,TINY=1.0E-20,ZERO=0.D0,ONE=1.D0)
!
      DIMENSION INDX(N),A(NP,NP),VV(NMAX)
!
      D=ONE
      DO I=1,N
        AAMAX=ZERO
        DO J=1,N
          IF(ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
        END DO
        IF(AAMAX.EQ.0.) PAUSE 'SINGULAR MATRIX IN LUDCMP'
        VV(I)=ONE/AAMAX
      END DO
!
      DO J=1,N
        DO I=1,J-1
          SUM=A(I,J)
          DO K=1,I-1
            SUM=SUM-A(I,K)*A(K,J)
          END DO
          A(I,J)=SUM
        END DO
        AAMAX=ZERO
        DO I=J,N
          SUM=A(I,J)
          DO K=1,J-1
            SUM=SUM-A(I,K)*A(K,J)
          END DO
          A(I,J)=SUM
          DUM=VV(I)*ABS(SUM)
          IF(DUM.GE.AAMAX) THEN
             IMAX=I
             AAMAX=DUM
          END IF
        END DO
        IF(J.NE.IMAX) THEN
           DO K=1,N
             DUM=A(IMAX,K)
             A(IMAX,K)=A(J,K)
             A(J,K)=DUM
           END DO
           D=-D
           VV(IMAX)=-VV(J)
        END IF
        INDX(J)=IMAX
        IF(A(J,J).EQ.0.) A(J,J)=TINY
        IF(J.NE.N) THEN
           DUM=ONE/A(J,J)
           DO I=J+1,N
             A(I,J)=A(I,J)*DUM
           END DO
        END IF
      END DO
!
      RETURN
!
      END
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
!   SUBROUTINE KLUBKSB                                                 C
!                                                                      C
!   FORWARD SUBSTITUTION                                               C
!   TAKEN FROM NUMERICAL RECIPES By Press et al                        C
!                                                                      C
!                                                                      C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!
      SUBROUTINE KLUBKSB(A,N,NP,INDX,B)
!
      IMPLICIT REAL*8(A-H,O-Z)
!
      PARAMETER (ZERO=0.D0)
!
      DIMENSION INDX(N),A(NP,NP),B(NP)
!
      II=0
      DO I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF(II.NE.0) THEN
           DO J=II,I-1
             SUM=SUM-A(I,J)*B(J)
           END DO
        ELSE IF(SUM.NE.ZERO) THEN
           II=I
        END IF
        B(I)=SUM
      END DO
!
      DO I=N,1,-1
        SUM=B(I)
        DO J=I+1,N
          SUM=SUM-A(I,J)*B(J)
        END DO
        B(I)=SUM/A(I,I)
      END DO
!
      RETURN
!
      END
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                      C
!     SUBROUTINE KCOPYMAT                                              C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      SUBROUTINE KCOPYMAT(A,B,N)
!
      IMPLICIT REAL*8(A-H,O-Z)
!
      DIMENSION A(N,N),B(N,N)
!
      CALL KCLEAR(B,N,N)
!
      DO K1=1,N
        DO K2=1,N
          B(K1,K2)=A(K1,K2)
        END DO
      END DO
!
      RETURN
!
      END
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
! SUBROUTINE SURPRESS                                                  C
!                                                                      C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      SUBROUTINE KSURTRAC(IDFACE,PMAG,NDOFEL,NNODE,COORD,SR)
!
      IMPLICIT REAL*8(A-H,O-Z)
!
      DIMENSION COORD(2,NNODE),SR(NDOFEL),XGP(6),WGP(6),AUX(6,2),FF(6),F(2),VAUX(6),NODE(3),TH(6)
!
!     Reads pressure load data for each loaded element
!
      DO I=1,6
        TH(I)=1.0
      END DO
!
      CALL KCLEARV (SR,NDOFEL)
      CALL KCLEARV (VAUX,6)
!      CALL KCLEARIV(NODE,3)
!
!       Retrieves the element nodes according to the face
!       where the loads are applied
!
      GO TO(100,200,300,400) IDFACE
100   NODE(1)=1
      NODE(2)=2
      NODE(3)=5
      GO TO 600
200   NODE(1)=2
      NODE(2)=3
      NODE(3)=6
      GO TO 600
300   NODE(1)=3
      NODE(2)=4
      NODE(3)=7
      GO TO 600
400   NODE(1)=4
      NODE(2)=1
      NODE(3)=8
      GO TO 600
     
600   CONTINUE

!     Calculates, J-determinant,Gauss points
!     and weights, Shape matrix along the
!     element surface and assembles element
!     load vector
!
      XI=COORD(1,NODE(1))
      XM=COORD(1,NODE(3))
      XJ=COORD(1,NODE(2))
      YI=COORD(2,NODE(1))
      YM=COORD(2,NODE(3))
      YJ=COORD(2,NODE(2))
      DXS=(XI-XJ)**2
      DYS=(YI-YJ)**2
      DL=DSQRT(DXS+DYS)
      IF(XI.EQ.XJ) THEN
        IF(YI.GT.YJ) THEN
          F(1)=-PMAG
          F(2)=0.0
        ELSE
          F(1)=PMAG
          F(2)=0.0
        END IF
      ELSE
        IF(YI.EQ.YJ) THEN
          F(1)=0.0
          F(2)=PMAG
        ELSE
          SM=(YJ-YI)/(XJ-XI)
          CT=DXS/DL
          ST=DYS/DL
          IF(SM.LT.0) THEN
            F(1)=PMAG*CT
            F(2)=PMAG*ST
          ELSE
            IF(YI.GT.YJ) THEN
              F(1)=-PMAG*ST
              F(2)= PMAG*CT
            ELSE
              F(1)=PMAG*ST
              F(2)=-PMAG*CT
            END IF
          END IF
        END IF
      END IF
      IF(YJ.GT.YI) THEN
        F(1)= PMAG*(YI-YJ)/DL
        F(2)=-PMAG*(XI-XJ)/DL
      ELSE
        F(1)=-PMAG*(YI-YJ)/DL
        F(2)= PMAG*(XI-XJ)/DL
      END IF
      DETJAC=DL/2.0
      CALL KGAUSS1D(XGP,WGP)
      DO 10 N=1,6
        ETA =XGP(N)
        ALFA=WGP(N)
        CALL KSHAPE1D(AUX,ETA)
        T=TH(N)
        CALL KMMULT(AUX,6,2,F,2,1,FF)
        VAUX=VAUX+ALFA*T*FF*DETJAC
   10 CONTINUE
!
      SR(2*NODE(1)-1)=VAUX(1)
      SR(2*NODE(1))=VAUX(2)
      SR(2*NODE(2)-1)=VAUX(3)
      SR(2*NODE(2))=VAUX(4)
      SR(2*NODE(3)-1)=VAUX(5)
      SR(2*NODE(3))=VAUX(6)
!
      RETURN
!
      END
!
!
!23456789012345678901234567890123456789012345678901234567890123456789012
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
!   SUBROUTINE SHAPE1D                                                 C
!                                                                      C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!
      SUBROUTINE KSHAPE1D(SM,ETA)
!
      IMPLICIT REAL*8(A-H,O-Z)
!
      DIMENSION SM(6,2)
!
      CALL KCLEAR(SM,6,2)
!
      EP=1.0+ETA
      EM=1.0-ETA
      EMS=1.0-ETA**2
!
      SM(1,1)=0.5*EP-0.5*EMS
      SM(2,2)=SM(1,1)
      SM(3,1)=0.5*EM-0.5*EMS
      SM(4,2)=SM(3,1)
      SM(5,1)=EMS
      SM(6,2)=SM(5,1)
!
      RETURN
!
      END
!
!
!23456789012345678901234567890123456789012345678901234567890123456789012
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
!   SUBROUTINE GAUSS1D                                                 C
!                                                                      C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!23456789012345678901234567890123456789012345678901234567890123456789012
!
      SUBROUTINE KGAUSS1D(XP,XW)
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      DIMENSION XP(6),XW(6)
!
      XP(1)=-0.932469514203152
      XP(2)=-0.661209386466265
      XP(3)=-0.238619186083197
      XP(6)= 0.932469514203152
      XP(5)= 0.661209386466265
      XP(4)= 0.238619186083197
!
      XW(1)=0.171324492379170
      XW(2)=0.360761573048139
      XW(3)=0.467913934572691
      XW(6)=XW(1)
      XW(5)=XW(2)
      XW(4)=XW(3)
!
      RETURN
!
      END
