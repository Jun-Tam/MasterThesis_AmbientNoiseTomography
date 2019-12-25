      PROGRAM RAYLE2

C------------------------------- TYPE DECLARATION -----------------------------C
      DIMENSION  H(15),RHO(15),VP(15),VS(15),AP(1),AE(1),QP(15),
     1           QS(15),Y(4,40),P(5,40),WS(5,40),Y0(4),YIJ(42),
     2           PhSwKernel(40,40),PhaseVel(40),GroupVel(40)
      DATA  PI2/6.283185/
      EXTERNAL  RAYMRX
C------------------------------------------------------------------------------C

C--------------------------------- PARAMETERS ---------------------------------C
      FMN  = 0.15
      FMX  = 0.50
      DF   = 0.025
      CMN  = 2.0
      CMX  = 4.0
      DC   = 0.1
      DEPT = 0.0
      DEPY = 10.0
      TOL  = 0.001
      ITR  = 30
      NDIV = 1
      IY   = 1
      IP   = 1
      IA   = 1
      IQ   = 0
C------------------------------------------------------------------------------C

C--------------------------------- MODEL INPUT --------------------------------C
      I=0
      OPEN(100,FILE='Model.k')
      DO
      I=I+1
      READ(100,*,END=5) H(I),RHO(I),VP(I),VS(I)
      ENDDO
    5 CONTINUE
      CLOSE(100)
      IL=I-1
      L=I-1
      NumFrq = nint((FMX - FMN)/DF + 1)
      NumLyr = (L-1) * 3 +1
C------------------------------------------------------------------------------C

C-------------------------------- CALCULATION ---------------------------------C
      DO  15  J=0,NumFrq-1
      T = 1 / (FMN + J * DF)
      IF( T.LE.0 )  GO TO  15
      W = PI2/T
      CALL  RAYDSP(RAYMRX,H,RHO,VP,VS,AP,AE,L,W,CMN,CMX,DC,
     *             TOL,ITR,IA,C,UD,EKD,Y0,YIJ,IER1)
      IF( IER1.LT.0 .OR. IER1.GT.1 )  GO TO  15
      IF( IER1.EQ.1 )  DC = DC/2
      CALL  RAYEFX(H,RHO,VP,VS,QP,QS,L,LY,NDIV,W,C,UI,IQ,Q,
     *             EKD,EKI,ER,IY,Y,IP,P,NY,WS,IER2)

C--------------------------------------- PRINT ---------------------------------
C      WRITE(*,13)
C   13 FORMAT(/4X,'I',4X,'DEPTH',10X,'DW/DRHO',8X,'DW/DVP',9X,
C     *           'DW/DVS',9X,'DW/DXI',9X,'DW/DPHI',8X,'DW/DETA')
C      CALL  DSPRNX(P,5,NY,NDIV,H,0.)
C------------------------------------------------------------------------------C
      DO K=1,NumLyr
      PhSwKernel(J+1,K)=P(3,K)
      ENDDO
      PhaseVel(J+1)=C
      GroupVel(J+1)=UD
      IF( CMN.GE.CMX )  GO TO  90
   15 CONTINUE
C------------------------------------------------------------------------------C

C----------------------------------- OUTPUT -----------------------------------C
      OPEN(110,FILE='RwPhVel.k')
      OPEN(120,FILE='RwPhSwKernel.k')
      OPEN(130,FILE='RwGrVel.k')
      DO J=1,NumFrq
      WRITE(110,*) PhaseVel(J)
      WRITE(130,*) GroupVel(J)
      ENDDO
      DO K=1,NumLyr
      WRITE(120,31) (PhSwKernel(J,K),J=1,NumFrq)
   31 FORMAT(13E11.3)
      ENDDO
      CLOSE(110)
      CLOSE(120)
      CLOSE(130)
C------------------------------------------------------------------------------C
   90 STOP
      END


      SUBROUTINE  RAYDSP(DIFFEQ,H,RHO,VP,VS,AP,AE,L,                    RAY02770
     *                   W,CMN,CMX,DC,TOL,ITR,IA,C,U,EK,                RAY02780
     1                   Y0,YIJ,IER)                                    RAY02790
C                                                                       RAY02800
C     INPUT                                                             RAY02810
C       DIFFEQ : SUBROUTINE TO INTEGRATE EQ. OF MOTION                  RAY02820
C       H    : LAYER THICKNESS.  H(L) IS ARBITRARY                      RAY02830
C       RHO  : DENSITY                                                  RAY02840
C       VP   : COMPRESSIONAL WAVE VELOCITY                              RAY02850
C       VS   : SHEAR WAVE VELOCITY                                      RAY02860
C       AP   : ANISOTROPY FACTOR PHI                                    RAY02870
C       AE   : ANISOTROPY FACTOR ETA                                    RAY02880
C       L    : NUMBER OF LAYERS INCLUDING THE BOTTOM                    RAY02890
C              HALF-SPACE                                               RAY02900
C       W    : ANGULAR FREQUENCY                                        RAY02910
C       CMN  : LOWER LIMIT OF PHASE VELOCITY                            RAY02920
C       CMX  : UPPER LIMIT OF PHASE VELOCITY                            RAY02930
C       DC   : INCREMENT   OF PHASE VELOCITY                            RAY02940
C       TOL  : RELATIVE ACCURACY OF PHASE VELOCITY                      RAY02950
C       ITR  : MAXIMUM NUMBER OF ITERATIONS                             RAY02960
C       IA   : = 0 ;   ISOTROPIC MODEL                                  RAY02970
C              = 1 ; ANISOTROPIC MODEL                                  RAY02980
C     OUTPUT                                                            RAY02990
C       C    : PHASE VELOCITY                                           RAY03000
C       U    : GROUP VELOCITY BY DIFFERENTIATION                        RAY03010
C       EKD  : ENERGY INTEGRAL BY DIFFERENTIATION                       RAY03020
C              2*K**2*I3                                                RAY03030
C       Y0   : SURFACE VALUES OF EIGENFUNCTION                          RAY03040
C              Y0(1) ; Y1 (SCALE FACTOR)                                RAY03050
C              Y0(2) ; Y2/ABS(Y1), DISPERSION FUNCTION                  RAY03060
C              Y0(3) ; Y3/Y1                                            RAY03070
C       YIJ  : SURFACE VALUES OF YIJ (COMPOUNDS),                       RAY03080
C              C*DYIJ/DC, AND W*DYIJ/DW                                 RAY03090
C       IER  : RETURN CODE                                              RAY03100
C              < 0 ; INPUT ERROR                                        RAY03110
C              = 0 ; NO ERROR                                           RAY03120
C              = 1 ; SLOW CONVERGENCE                                   RAY03130
C              = 2 ; ROOT NOT FOUND                                     RAY03140
C                                                                       RAY03150
C     SUBROUTINE : DIFFEQ (RAYMRX OR RAYRKG)                            RAY03160
C                                                                       RAY03170
C     DISPER-80,  VER-86                                                RAY03180
C                                                                       RAY03190
C     M. SAITO  23/VI/79                                                RAY03200
C     REVISED   25/IX/85                                                RAY03210
C                                                                       RAY03220
      DIMENSION  H(L),RHO(L),VP(L),VS(L),AP(L),AE(L),                   RAY03230
     *           Y0(3),YIJ(15)                                          RAY03240
C                                                                       RAY03250
C     INITIALIZATION                                                    RAY03260
C                                                                       RAY03270
      IF( L.LE.2 .OR. W.LE.0 .OR. DC.EQ.0 )  GO TO  90                  RAY03280
      ONE = 1                                                           RAY03290
C                                                                       RAY03300
C     MACHINE EPSILON                                                   RAY03310
C                                                                       RAY03320
      EPS = 1                                                           RAY03330
    1 IF( (1+EPS).LE.1 )  GO TO  2                                      RAY03340
        EPS = EPS/2                                                     RAY03350
        GO TO  1                                                        RAY03360
    2 EPS = EPS/2                                                       RAY03370
C                                                                       RAY03380
      TOL1 = TOL                                                        RAY03390
      IF( TOL1.GT.0 )  TOL1 = MAX(TOL1,EPS)                             RAY03400
      IER1 = 0                                                          RAY03410
      WRITE(6,3)  W                                                     RAY03420
    3   FORMAT(/7X,'W',15X,'RAYLEIGH WAVE'/                             RAY03430
     *         1PE18.6//7X,'C',17X,'Y2',16X,'Y3')                       RAY03440
C    *         1PD18.6//7X,'C',17X,'Y2',16X,'Y3')                       RAY03450
      C3 = CMN                                                          RAY03460
      IF( C3.LE.0 )  GO TO  92                                          RAY03470
      CALL  DIFFEQ(H,RHO,VP,VS,AP,AE,L,W,C3,IA,1,U,                     RAY03480
     *             EK,Y0,YIJ,IER)                                       RAY03490
      IF( IER.NE.0 )  RETURN                                            RAY03500
      F3 = Y0(2)                                                        RAY03510
      WRITE(6,4)  C3,F3,Y0(3)                                           RAY03520
    4   FORMAT(1P3E18.6)                                                RAY03530
C   4   FORMAT(1P3D18.6)                                                RAY03540
      IF( F3.EQ.0 .AND. TOL1.GT.0 )  GO TO  9                           RAY03550
C                                                                       RAY03560
C     FIND A ZERO-CROSS                                                 RAY03570
C                                                                       RAY03580
      KX = (CMX - CMN)/DC + 0.5                                         RAY03590
      KX = MAX(1,KX)                                                    RAY03600
      DO  5  K=1,KX                                                     RAY03610
        CC = CMN + K*DC                                                 RAY03620
        C1 = C3                                                         RAY03630
        F1 = F3                                                         RAY03640
        C3 = CC                                                         RAY03650
        IF( C3.LE.0 )  GO TO  92                                        RAY03660
        CALL  DIFFEQ(H,RHO,VP,VS,AP,AE,L,W,C3,IA,1,U,                   RAY03670
     *               EK,Y0,YIJ,IER)                                     RAY03680
        IF( IER.NE.0 )  RETURN                                          RAY03690
        F3 = Y0(2)                                                      RAY03700
        WRITE(6,4)  C3,F3,Y0(3)                                         RAY03710
        IF(  TOL1.LE.0 )  GO TO  5                                      RAY03720
        IF( F3*SIGN(ONE,F1).LE.0 )  GO TO  6                            RAY03730
    5 CONTINUE                                                          RAY03740
C                                                                       RAY03750
      IER = 2                                                           RAY03760
      IF( TOL1.LE.0 )  RETURN                                           RAY03770
      GO TO  94                                                         RAY03780
C                                                                       RAY03790
C     INTERPOLATION                                                     RAY03800
C                                                                       RAY03810
    6 IF( F3.EQ.0 )  GO TO  9                                           RAY03820
      C2 = C3                                                           RAY03830
      F2 = F3                                                           RAY03840
      E  = C1 - C2                                                      RAY03850
      D  = E/2                                                          RAY03860
      C3 = C2 + D                                                       RAY03870
      KX = MAX(1,ITR)                                                   RAY03880
C                                                                       RAY03890
      DO  7  K=1,KX                                                     RAY03900
        CALL  DIFFEQ(H,RHO,VP,VS,AP,AE,L,W,C3,IA,1,U,                   RAY03910
     *               EK,Y0,YIJ,IER)                                     RAY03920
        F3 = Y0(2)                                                      RAY03930
        WRITE(6,4) C3,F3,Y0(3)                                          RAY03940
        IF( F3*SIGN(ONE,F2).GT.0 )  THEN                                RAY03950
          FF = C1                                                       RAY03960
          C1 = C2                                                       RAY03970
          C2 = FF                                                       RAY03980
          FF = F1                                                       RAY03990
          F1 = F2                                                       RAY04000
          F2 = FF                                                       RAY04010
        END IF                                                          RAY04020
        IF( ABS(F3).GT.ABS(F2) )  THEN                                  RAY04030
          FF = C2                                                       RAY04040
          C2 = C3                                                       RAY04050
          C3 = FF                                                       RAY04060
          FF = F2                                                       RAY04070
          F2 = F3                                                       RAY04080
          F3 = FF                                                       RAY04090
        END IF                                                          RAY04100
        E = C2 - C3                                                     RAY04110
        IF( F3.EQ.0 )  GO TO  9                                         RAY04120
        TOLC = C3*TOL1                                                  RAY04130
        DD = D                                                          RAY04140
        F32 = F3/F2                                                     RAY04150
        F31 = F3/F1                                                     RAY04160
        F21 = F2/F1                                                     RAY04170
        Q = F32*(E*(1 - F31) + F21*(F31 - F21)*(C1 - C3))               RAY04180
        S = (F21 - 1)*(F32 - 1)*(F31 - 1)                               RAY04190
C                                                                       RAY04200
C       TEST RANGE                                                      RAY04210
C                                                                       RAY04220
        IF( Q.LT.0 )  S =-S                                             RAY04230
        Q = ABS(Q)                                                      RAY04240
        IF( Q.GE.(E*S-ABS(TOLC*S)) )  THEN                              RAY04250
C                                                                       RAY04260
C       LINEAR INTERPOLATION                                            RAY04270
C                                                                       RAY04280
          D = E*F32/(F32 - 1)                                           RAY04290
C                                                                       RAY04300
C       INVERSE QUADRATIC INTERPOLATION                                 RAY04310
C                                                                       RAY04320
        ELSE                                                            RAY04330
          D = Q/S                                                       RAY04340
        END IF                                                          RAY04350
C                                                                       RAY04360
C       TEST CONVERGENCE                                                RAY04370
C                                                                       RAY04380
        C1 = C2                                                         RAY04390
        F1 = F2                                                         RAY04400
        C2 = C3                                                         RAY04410
        F2 = F3                                                         RAY04420
        C3 = C2 + D                                                     RAY04430
        IF( ABS(E).LE.TOLC )  GO TO  9                                  RAY04440
        IF( ABS(D).LE.TOLC )  THEN                                      RAY04450
C                                                                       RAY04460
C       BISECTION                                                       RAY04470
C                                                                       RAY04480
          IF( ABS(DD).LE.TOLC )  THEN                                   RAY04490
            D = E/2                                                     RAY04500
            C3 = C2 + D                                                 RAY04510
          ELSE                                                          RAY04520
            C3 = C2 + SIGN(TOLC,D)                                      RAY04530
          END IF                                                        RAY04540
        END IF                                                          RAY04550
C                                                                       RAY04560
    7 CONTINUE                                                          RAY04570
C                                                                       RAY04580
C     SLOW CONVERGENCE                                                  RAY04590
C                                                                       RAY04600
      WRITE(6,8)  KX                                                    RAY04610
    8   FORMAT(20X,5('?'),3X,'(RAYDSP)   SLOW CONV. ',                  RAY04620
     *         'AFTER',I5,' ITERATIONS',3X,5('?'))                      RAY04630
      IER = 1                                                           RAY04640
C                                                                       RAY04650
C     ROOT IS FOUND                                                     RAY04660
C                                                                       RAY04670
    9 CALL  DIFFEQ(H,RHO,VP,VS,AP,AE,L,W,C3,IA,3,U,                     RAY04680
     *             EK,Y0,YIJ,IER1)                                      RAY04690
      WRITE(6,4)  C3,Y0(2),Y0(3)                                        RAY04700
      C   = C3                                                          RAY04710
      RETURN                                                            RAY04720
C                                                                       RAY04730
C     INPUT ERROR                                                       RAY04740
C                                                                       RAY04750
   90 WRITE(6,91)  L,W,DC                                               RAY04760
   91   FORMAT(20X,5('?'),3X,'(RAYDSP)   INPUT ERROR',3X,               RAY04770
     *         'L  =',I5,3X,'W  =',1PE13.6,3X,                          RAY04780
     1         'DC =',E13.6,3X,5('?'))                                  RAY04790
      IER =-1                                                           RAY04800
      RETURN                                                            RAY04810
C                                                                       RAY04820
   92 WRITE(6,93)  C3                                                   RAY04830
   93   FORMAT(20X,5('?'),3X,'(RAYDSP)   INPUT ERROR',3X,               RAY04840
     *         'C  =',1PE13.6,3X,5('?'))                                RAY04850
      IER =-1                                                           RAY04860
      RETURN                                                            RAY04870
C                                                                       RAY04880
C     NO ROOT                                                           RAY04890
C                                                                       RAY04900
   94 WRITE(6,95)                                                       RAY04910
   95   FORMAT(20X,5('?'),3X,'(RAYDSP)',3X,                             RAY04920
     *         'ROOT NOT FOUND',3X,5('?'))                              RAY04930
      RETURN                                                            RAY04940
      END                                                               RAY04950
C                                                                       RAY04960
C                                                                       RAY08420
C     RAYLEIGH WAVE MATRIX METHOD INTEGRATION                           RAY08430
C                                                                       RAY08440
      SUBROUTINE  RAYMRX(H,RHO,VP,VS,AP,AE,L,W,C,                       RAY08450
     *                   IA,IG,U,EK,Y0,Y,IER)                           RAY08460
C                                                                       RAY08470
C     INPUT                                                             RAY08480
C       H    : LAYER THICKNESS.  H(L) IS ARBITRARY                      RAY08490
C       RHO  : DENSITY                                                  RAY08500
C       VP   : P WAVE VELOCITY                                          RAY08510
C       VS   : S WAVE VELOCITY                                          RAY08520
C       L    : NUMBER OF LAYERS INCLDING THE BOTTOM                     RAY08530
C              HALF-SPACE                                               RAY08540
C       W    : ANGULAR FREQUENCY                                        RAY08550
C       C    : PHASE VELOCITY                                           RAY08560
C       IG   : = 1 ; TO COMPUTE Y ONLY                                  RAY08570
C              = 2 ; TO COMPUTE Y AND C*CY/DC                           RAY08580
C              = 3 ; TO COMPUTE Y, C*DY/DC, AND W*DY/DW                 RAY08590
C     OUTPUT                                                            RAY08600
C       U    : GROUP VELOCITY BY DIFFERENTIATION (IG=3)                 RAY08610
C       EK   : ENERGY INTEGRAL 2*K**2*I3  (IG>=2)                       RAY08620
C       Y0   : SURFACE VALUES OF EIGENFUNCTIONS                         RAY08630
C              Y0(1) ; Y1 (SCALE FACTOR)                                RAY08640
C              Y0(2) ; Y2/ABS(Y1), DISPERSION FUNCTION                  RAY08650
C              Y0(3) ; Y3/Y1                                            RAY08660
C       Y    : SURFACE VALUES OF YIJ (COMPOUNDS),                       RAY08670
C              C*DYIJ/DC, AND W*DYIJ/DW                                 RAY08680
C              FOR SOLID LAYER                                          RAY08690
C                (IJ) = (12),(13),(14),(23),(24)                        RAY08700
C       IER  : RETURN CODE                                              RAY08710
C              < 0 ; INPUT ERROR                                        RAY08720
C              = 0 ; NO ERROR                                           RAY08730
C                                                                       RAY08740
C     ISOTROPIC MODEL ONLY.  AP, AE, IA ARE DUMMY.                      RAY08750
C                                                                       RAY08760
C     DISPER-80,  VER-86                                                RAY08770
C                                                                       RAY08780
C     M. SAITO  30/VII/79                                               RAY08790
C     REVISED   10/XII/86                                               RAY08800
C                                                                       RAY08810
      DIMENSION  H(L),RHO(L),VP(L),VS(L),Y0(3),Y(15),Z(15)              RAY08820
     *          ,AP(1),AE(1)                                            RAY08830
      DATA  EPS/1.E-10/,BIG/1.E+10/                                     RAY08840
C                                                                       RAY08850
C     DEFINE SINH(X)/X AND (COSH(X)-SINH(X)/X)/X**2                     RAY08860
C                                                                       RAY08870
      SH0(X) = 0.9999997 + X*(0.1666667 + X*(0.0083361                  RAY08880
     *       + X*0.0001984))                                            RAY08890
      SH1(X) = 0.3333333 + X*(0.0333333 + X*(0.0011907                  RAY08900
     *       + X*0.0000220))                                            RAY08910
C                                                                       RAY08920
C     FOR DOUBLE PRECISION USE                                          RAY08930
C                                                                       RAY08940
C     SH0(X) = 1.0D0 + X*(1.6666 66666 66666 7D-1                       RAY08950
C    *       + X*(8.33 33333 33334 0D-3                                 RAY08960
C    1       + X*(1.9 84126 98412 7D-4                                  RAY08970
C    2       + X*(2.7557 31918 9D-6 + X*(2.50 12108 4D-8                RAY08980
C    3       + X*(1.60596 1D-10 + X*7.64 7D-13))))))                    RAY08990
C     SH1(X) = 3.3333 33333 33333 3D-1                                  RAY09000
C    *       + X*(3.333 33333 33333 3D-2                                RAY09010
C    1       + X*(1.19 04761 90476 2D-3                                 RAY09020
C    2       + X*(2.20458 55379 2D-5                                    RAY09030
C    3       + X*(2.505 21083 7D-7 + X*(1.9 27085 3D-9                  RAY09040
C    4       + X*(1.0706 3D-11 + X*4.50D-14))))))                       RAY09050
C                                                                       RAY09060
C     INITIAL VALUE                                                     RAY09070
C                                                                       RAY09080
      IF( L.LE.0 .OR. W.LE.0 .OR. C.LE.0 )  GO TO  90                   RAY09090
      I   = L                                                           RAY09100
      RO  = RHO(I)                                                      RAY09110
      IF( RO.LE.0 .OR. C.GE.VP(I) )  GO TO  92                          RAY09120
      IER = 0                                                           RAY09130
      CC  = C*C                                                         RAY09140
      WN  = W/C                                                         RAY09150
      IGG = MAX(1,MIN(3,IG))                                            RAY09160
      ROC = RO*CC                                                       RAY09170
      SV  = VS(I)                                                       RAY09180
      CP  = C/VP(I)                                                     RAY09190
      RAA = (1 + CP)*(1 - CP)                                           RAY09200
      RA  = SQRT(RAA)                                                   RAY09210
      DO  1  J=1,15                                                     RAY09220
    1   Y(J) = 0                                                        RAY09230
C                                                                       RAY09240
C     LIQUID BOTTOM                                                     RAY09250
C                                                                       RAY09260
      IF( SV.LE.0 )  THEN                                               RAY09270
        Y(1) = RA*EPS                                                   RAY09280
        Y(2) =-ROC*EPS                                                  RAY09290
        JX = 2                                                          RAY09300
C                                                                       RAY09310
        IF( IGG.GE.2 )  THEN                                            RAY09320
          Y(3) =-CP**2*EPS/RA                                           RAY09330
          Y(4) =-2*ROC*EPS                                              RAY09340
          JX = 4                                                        RAY09350
C                                                                       RAY09360
          IF( IGG.GT.2 )  JX = 6                                        RAY09370
        ENDIF                                                           RAY09380
C                                                                       RAY09390
C     SOLID BOTTOM                                                      RAY09400
C                                                                       RAY09410
      ELSE                                                              RAY09420
        IF( C.GE.SV )  GO TO  92                                        RAY09430
        CS  = C/SV                                                      RAY09440
        RBB = (1 + CS)*(1 - CS)                                         RAY09450
        RB  = SQRT(RBB)                                                 RAY09460
        RG  = 2*RO*VS(I)**2                                             RAY09470
        Y(3) =-RA*EPS                                                   RAY09480
        Y(4) =-RB*EPS                                                   RAY09490
        Y(2) =-EPS*(CP**2*RBB + CS**2)/(ROC*(RA*RB + 1))                RAY09500
        Y(1) = RG*Y(2) + EPS                                            RAY09510
        Y(5) =-RG*(Y(1) + EPS) + ROC*EPS                                RAY09520
        JX = 5                                                          RAY09530
C                                                                       RAY09540
        IF( IGG.GE.2 )  THEN                                            RAY09550
          Y(8) = EPS*CP**2/RA                                           RAY09560
          Y(9) = EPS*CS**2/RB                                           RAY09570
          Y(7) =-(RB*Y(8) + RA*Y(9))/ROC - 2*Y(2)                       RAY09580
          Y(6) = RG*Y(7)                                                RAY09590
          Y(10)=-RG*Y(6) + EPS*ROC*2                                    RAY09600
          JX = 10                                                       RAY09610
C                                                                       RAY09620
          IF( IGG.GT.2 )  JX = 15                                       RAY09630
        ENDIF                                                           RAY09640
      ENDIF                                                             RAY09650
C                                                                       RAY09660
C     INTEGRATE UPWARD                                                  RAY09670
C                                                                       RAY09680
      IF( L.LE.1 )  GO TO  8                                            RAY09690
      DO  7  II=2,L                                                     RAY09700
        I  = I - 1                                                      RAY09710
        RO = RHO(I)                                                     RAY09720
        ROC= RO*CC                                                      RAY09730
        PV = VP(I)                                                      RAY09740
        SV = VS(I)                                                      RAY09750
        IF( PV.LE.0 .OR. RO.LE.0 .OR.                                   RAY09760
     *      H(I).LE.0 )  GO TO  92                                      RAY09770
        DO  2  J=1,15                                                   RAY09780
    2     Z(J) = Y(J)                                                   RAY09790
        IF( (SV.LE.0 .AND. VS(I+1).LE.0) .OR.                           RAY09800
     *      (SV.GT.0 .AND. VS(I+1).GT.0) )  GO TO  3                    RAY09810
C                                                                       RAY09820
C       SOLID  -----> LIQUID BOUNDARY                                   RAY09830
C                                                                       RAY09840
        IF( SV.LE.0 )  THEN                                             RAY09850
          Y0(3) =-Y(1)/Y(3)                                             RAY09860
          Z(1) = Y(3)                                                   RAY09870
          Z(2) = Y(5)                                                   RAY09880
          JX = 2                                                        RAY09890
C                                                                       RAY09900
          IF( IGG.GE.2 )  THEN                                          RAY09910
            Z(3) = Y(8)                                                 RAY09920
            Z(4) = Y(10)                                                RAY09930
            JX = 4                                                      RAY09940
C                                                                       RAY09950
            IF( IGG.GT.2 )  THEN                                        RAY09960
              Z(5) = Y(13)                                              RAY09970
              Z(6) = Y(15)                                              RAY09980
              JX = 6                                                    RAY09990
            ENDIF                                                       RAY10000
          ENDIF                                                         RAY10010
C                                                                       RAY10020
C       LIQUID -----> SOLID BOUNDARY                                    RAY10030
C                                                                       RAY10040
        ELSE                                                            RAY10050
          Z( 2) = Y(1)                                                  RAY10060
          Z( 4) = Y(2)                                                  RAY10070
          Z( 1) = 0                                                     RAY10080
          Z( 3) = 0                                                     RAY10090
          Z( 5) = 0                                                     RAY10100
          JX = 5                                                        RAY10110
C                                                                       RAY10120
          IF( IGG.GE.2 )  THEN                                          RAY10130
            Z( 7) = Y(3)                                                RAY10140
            Z( 9) = Y(4)                                                RAY10150
            Z( 6) = 0                                                   RAY10160
            Z( 8) = 0                                                   RAY10170
            Z(10) = 0                                                   RAY10180
            JX = 10                                                     RAY10190
C                                                                       RAY10200
            IF( IGG.GT.2 )  THEN                                        RAY10210
              Z(12) = Y(5)                                              RAY10220
              Z(14) = Y(6)                                              RAY10230
              Z(11) = 0                                                 RAY10240
              Z(13) = 0                                                 RAY10250
              Z(15) = 0                                                 RAY10260
              JX = 15                                                   RAY10270
            ENDIF                                                       RAY10280
          ENDIF                                                         RAY10290
        ENDIF                                                           RAY10300
C                                                                       RAY10310
    3   R2  = 1/ROC                                                     RAY10320
        CP  = C/PV                                                      RAY10330
        RAA = (1 + CP)*(1 - CP)                                         RAY10340
        CS  = 0                                                         RAY10350
        IF( SV.GT.0 )  CS = C/SV                                        RAY10360
        RBB = (1 + CS)*(1 - CS)                                         RAY10370
        HK  = H(I)*WN                                                   RAY10380
        HKK = HK**2                                                     RAY10390
        XX  = RAA*HKK                                                   RAY10400
        ONE = 1                                                         RAY10410
C                                                                       RAY10420
C       SINH(X)/X                                                       RAY10430
C                                                                       RAY10440
        DO  4  K=1,2                                                    RAY10450
          CHA = CHB                                                     RAY10460
          SHA = SHB                                                     RAY10470
          DHA = DHB                                                     RAY10480
          AA  = ABS(XX)                                                 RAY10490
          IF( AA.LE.1 )  THEN                                           RAY10500
            SHB = SH0(XX)                                               RAY10510
            CHB = 1 + XX*SH0(XX/4)**2/2                                 RAY10520
            IF( IGG.GE.2 )  DHB = SH1(XX)*HKK                           RAY10530
          ELSE                                                          RAY10540
            AA  = SQRT(AA)                                              RAY10550
            IF( XX.LE.0 )  THEN                                         RAY10560
              CHB = COS(AA)                                             RAY10570
              SHB = SIN(AA)/AA                                          RAY10580
            ELSE                                                        RAY10590
              IF( AA.GT.100 )  ONE = 0                                  RAY10600
              IF( AA.LE.100 )  ONE = ONE/COSH(AA)                       RAY10610
              CHB = 1                                                   RAY10620
              SHB = TANH(AA)/AA                                         RAY10630
            ENDIF                                                       RAY10640
            IF( IGG.GE.2 )  DHB = (HKK/XX)*(CHB - SHB)                  RAY10650
          ENDIF                                                         RAY10660
          XX  = HKK*RBB                                                 RAY10670
          SHB = HK*SHB                                                  RAY10680
    4   CONTINUE                                                        RAY10690
C                                                                       RAY10700
C     LAYER MATRICES                                                    RAY10710
C                                                                       RAY10720
C       LIQUID LAYER                                                    RAY10730
C                                                                       RAY10740
        IF( SV.LE.0 )  THEN                                             RAY10750
          B11 = CHA                                                     RAY10760
          B12 =-RAA*SHA*R2                                              RAY10770
          B21 =-ROC*SHA                                                 RAY10780
C                                                                       RAY10790
          Y(1) = B11*Z(1) + B12*Z(2)                                    RAY10800
          Y(2) = B21*Z(1) + B11*Z(2)                                    RAY10810
C                                                                       RAY10820
          IF( IGG.GE.2 )  THEN                                          RAY10830
            C11 =-HK*SHA                                                RAY10840
            C12 = (HK*CHA + (1 + RAA)*SHA)*R2                           RAY10850
            C21 = (HK*DHA - SHA)*ROC                                    RAY10860
C                                                                       RAY10870
            Y(3) = B11*Z(3) + B12*Z(4) + C11*Z(1)                       RAY10880
     *           + C12*Z(2)                                             RAY10890
            Y(4) = B21*Z(3) + B11*Z(4) + C21*Z(1)                       RAY10900
     *           + C11*Z(2)                                             RAY10910
C                                                                       RAY10920
            IF( IGG.GT.2 )  THEN                                        RAY10930
              W11 = HK*RAA*SHA                                          RAY10940
              W12 =-HK*RAA*CHA*R2                                       RAY10950
              W21 =-HK*ROC*CHA                                          RAY10960
C                                                                       RAY10970
              Y(5) = B11*Z(5) + B12*Z(6) + W11*Z(1)                     RAY10980
     *             + W12*Z(2)                                           RAY10990
              Y(6) = B21*Z(5) + B11*Z(6) + W21*Z(1)                     RAY11000
     *             + W11*Z(2)                                           RAY11010
            ENDIF                                                       RAY11020
          ENDIF                                                         RAY11030
C                                                                       RAY11040
C       SOLID LAYER                                                     RAY11050
C                                                                       RAY11060
        ELSE                                                            RAY11070
          G1  = 2/CS**2                                                 RAY11080
          RG  = G1*ROC                                                  RAY11090
          R4  = RG - ROC                                                RAY11100
          E1  = CHA*CHB                                                 RAY11110
          E2  = E1 - ONE                                                RAY11120
          E3  = SHA*SHB                                                 RAY11130
          E5  = SHA*CHB                                                 RAY11140
          E6  = SHB*CHA                                                 RAY11150
          F1  = E2 - E3                                                 RAY11160
          F2  = R2*F1                                                   RAY11170
          F3  = G1*F1 + E3                                              RAY11180
          B33 = E1                                                      RAY11190
          B34 = RAA*E3                                                  RAY11200
          B43 = RBB*E3                                                  RAY11210
          B25 =-R2*(F2 + R2*(E2 - RAA*B43))                             RAY11220
          B15 = RG*B25 + F2                                             RAY11230
          B16 =-RG*B15 - F3                                             RAY11240
          B22 = B16 + E1                                                RAY11250
          B12 = RG*B16 - R4*F3                                          RAY11260
          B52 =-RG*B12 + R4*(RG*F3 + R4*E3)                             RAY11270
          B23 = R2*(E5 - RBB*E6)                                        RAY11280
          B13 = RG*B23 - E5                                             RAY11290
          B42 =-RG*B13 + R4*E5                                          RAY11300
          B24 = R2*(E6 - RAA*E5)                                        RAY11310
          B14 = RG*B24 - E6                                             RAY11320
          B32 =-RG*B14 + R4*E6                                          RAY11330
          B11 = ONE - B16-B16                                           RAY11340
          B21 = B15+B15                                                 RAY11350
          B31 = B14+B14                                                 RAY11360
          B41 = B13+B13                                                 RAY11370
          B51 = B12+B12                                                 RAY11380
C                                                                       RAY11390
          Y(1) = B11*Z(1) + B12*Z(2) + B13*Z(3)                         RAY11400
     *         + B14*Z(4) + B15*Z(5)                                    RAY11410
          Y(2) = B21*Z(1) + B22*Z(2) + B23*Z(3)                         RAY11420
     *         + B24*Z(4) + B25*Z(5)                                    RAY11430
          Y(3) = B31*Z(1) + B32*Z(2) + B33*Z(3)                         RAY11440
     *         + B34*Z(4) + B24*Z(5)                                    RAY11450
          Y(4) = B41*Z(1) + B42*Z(2) + B43*Z(3)                         RAY11460
     *         + B33*Z(4) + B23*Z(5)                                    RAY11470
          Y(5) = B51*Z(1) + B52*Z(2) + B42*Z(3)                         RAY11480
     *         + B32*Z(4) + B22*Z(5)                                    RAY11490
C                                                                       RAY11500
          IF( IGG.GE.2 )  THEN                                          RAY11510
            RAAC =-2*CP*CP                                              RAY11520
            RBBC =-2*CS*CS                                              RAY11530
            R1C = ROC+ROC                                               RAY11540
            E1C =-HK*(E5 + E6)                                          RAY11550
            E3C =-E3-E3 - HK*(DHA*SHB + DHB*SHA)                        RAY11560
            E5C =-E5 - HK*(DHA*CHB + E3)                                RAY11570
            E6C =-E6 - HK*(DHB*CHA + E3)                                RAY11580
            F1C = E1C - E3C                                             RAY11590
            F2C = R2*(F1C - F1-F1)                                      RAY11600
            F3C = G1*(F1C - F1-F1) + E3C                                RAY11610
            C33 = E1C                                                   RAY11620
            C34 = RAA*E3C + RAAC*E3                                     RAY11630
            C43 = RBB*E3C + RBBC*E3                                     RAY11640
            C25 =-R2*(F2C + R2*(E1C - RAA*C43 - RAAC*B43))              RAY11650
     *          - 2*(B25+B25 + R2*F2)                                   RAY11660
            C15 = RG*C25 + F2C                                          RAY11670
            C16 =-RG*C15 - F3C                                          RAY11680
            C22 = C16 + E1C                                             RAY11690
            C12 = RG*C16 + R1C*F3 - R4*F3C                              RAY11700
            C52 =-RG*C12 + R4*(RG*F3C + R4*E3C)                         RAY11710
     *          - R1C*(RG*F3 + 2*R4*E3)                                 RAY11720
            C23 = R2*(E5C - RBB*E6C - RBBC*E6) - B23-B23                RAY11730
            C13 = RG*C23 - E5C                                          RAY11740
            C42 =-RG*C13 + R4*E5C - R1C*E5                              RAY11750
            C24 = R2*(E6C - RAA*E5C - RAAC*E5) - B24-B24                RAY11760
            C14 = RG*C24 - E6C                                          RAY11770
            C32 =-RG*C14 + R4*E6C - R1C*E6                              RAY11780
            C11 =-C16-C16                                               RAY11790
            C21 = C15+C15                                               RAY11800
            C31 = C14+C14                                               RAY11810
            C41 = C13+C13                                               RAY11820
            C51 = C12+C12                                               RAY11830
C                                                                       RAY11840
            Y( 6) = B11*Z(6) + B12*Z(7) + B13*Z(8)                      RAY11850
     *            + B14*Z(9) + B15*Z(10)                                RAY11860
     *            + C11*Z(1) + C12*Z(2) + C13*Z(3)                      RAY11870
     *            + C14*Z(4) + C15*Z(5)                                 RAY11880
            Y( 7) = B21*Z(6) + B22*Z(7) + B23*Z(8)                      RAY11890
     *            + B24*Z(9) + B25*Z(10)                                RAY11900
     *            + C21*Z(1) + C22*Z(2) + C23*Z(3)                      RAY11910
     *            + C24*Z(4) + C25*Z(5)                                 RAY11920
            Y( 8) = B31*Z(6) + B32*Z(7) + B33*Z(8)                      RAY11930
     *            + B34*Z(9) + B24*Z(10)                                RAY11940
     *            + C31*Z(1) + C32*Z(2) + C33*Z(3)                      RAY11950
     *            + C34*Z(4) + C24*Z(5)                                 RAY11960
            Y( 9) = B41*Z(6) + B42*Z(7) + B43*Z(8)                      RAY11970
     *            + B33*Z(9) + B23*Z(10)                                RAY11980
     *            + C41*Z(1) + C42*Z(2) + C43*Z(3)                      RAY11990
     *            + C33*Z(4) + C23*Z(5)                                 RAY12000
            Y(10) = B51*Z(6) + B52*Z(7) + B42*Z(8)                      RAY12010
     *            + B32*Z(9) + B22*Z(10)                                RAY12020
     *            + C51*Z(1) + C52*Z(2) + C42*Z(3)                      RAY12030
     *            + C32*Z(4) + C22*Z(5)                                 RAY12040
C                                                                       RAY12050
            IF( IGG.GT.2 )  THEN                                        RAY12060
              E1W = HK*(RAA*E5 + RBB*E6)                                RAY12070
              E3W = HK*(E5 + E6)                                        RAY12080
              E5W = HK*(E1 + B43)                                       RAY12090
              E6W = HK*(E1 + B34)                                       RAY12100
              F1W = E1W - E3W                                           RAY12110
              F2W = R2*F1W                                              RAY12120
              F3W = G1*F1W + E3W                                        RAY12130
              W33 = E1W                                                 RAY12140
              W34 = RAA*E3W                                             RAY12150
              W43 = RBB*E3W                                             RAY12160
              W25 =-R2*(F2W + R2*(E1W - RAA*W43))                       RAY12170
              W15 = RG*W25 + F2W                                        RAY12180
              W16 =-RG*W15 - F3W                                        RAY12190
              W22 = W16 + E1W                                           RAY12200
              W12 = RG*W16 - R4*F3W                                     RAY12210
              W52 =-RG*W12 + R4*(RG*F3W + R4*E3W)                       RAY12220
              W23 = R2*(E5W - RBB*E6W)                                  RAY12230
              W13 = RG*W23 - E5W                                        RAY12240
              W42 =-RG*W13 + R4*E5W                                     RAY12250
              W24 = R2*(E6W - RAA*E5W)                                  RAY12260
              W14 = RG*W24 - E6W                                        RAY12270
              W32 =-RG*W14 + R4*E6W                                     RAY12280
              W11 =-W16-W16                                             RAY12290
              W21 = W15+W15                                             RAY12300
              W31 = W14+W14                                             RAY12310
              W41 = W13+W13                                             RAY12320
              W51 = W12+W12                                             RAY12330
C                                                                       RAY12340
              Y(11) = B11*Z(11) + B12*Z(12) + B13*Z(13)                 RAY12350
     *              + B14*Z(14) + B15*Z(15)                             RAY12360
     1              + W11*Z( 1) + W12*Z( 2) + W13*Z( 3)                 RAY12370
     2              + W14*Z( 4) + W15*Z( 5)                             RAY12380
              Y(12) = B21*Z(11) + B22*Z(12) + B23*Z(13)                 RAY12390
     *              + B24*Z(14) + B25*Z(15)                             RAY12400
     1              + W21*Z( 1) + W22*Z( 2) + W23*Z( 3)                 RAY12410
     2              + W24*Z( 4) + W25*Z( 5)                             RAY12420
              Y(13) = B31*Z(11) + B32*Z(12) + B33*Z(13)                 RAY12430
     *              + B34*Z(14) + B24*Z(15)                             RAY12440
     1              + W31*Z( 1) + W32*Z( 2) + W33*Z( 3)                 RAY12450
     2              + W34*Z( 4) + W24*Z( 5)                             RAY12460
              Y(14) = B41*Z(11) + B42*Z(12) + B43*Z(13)                 RAY12470
     *              + B33*Z(14) + B23*Z(15)                             RAY12480
     1              + W41*Z( 1) + W42*Z( 2) + W43*Z( 3)                 RAY12490
     2              + W33*Z( 4) + W23*Z( 5)                             RAY12500
              Y(15) = B51*Z(11) + B52*Z(12) + B42*Z(13)                 RAY12510
     *              + B32*Z(14) + B22*Z(15)                             RAY12520
     1              + W51*Z( 1) + W52*Z( 2) + W42*Z( 3)                 RAY12530
     2              + W32*Z( 4) + W22*Z( 5)                             RAY12540
            ENDIF                                                       RAY12550
          ENDIF                                                         RAY12560
        ENDIF                                                           RAY12570
C                                                                       RAY12580
C     NORMALIZATION                                                     RAY12590
C                                                                       RAY12600
        Z1 = 0                                                          RAY12610
        DO  5  J=1,JX                                                   RAY12620
    5     Z1 = MAX(Z1,ABS(Y(J)))                                        RAY12630
        IF( Z1.GT.BIG )  THEN                                           RAY12640
          DO  6  J=1,JX                                                 RAY12650
    6       Y(J) = EPS*Y(J)                                             RAY12660
        ENDIF                                                           RAY12670
C                                                                       RAY12680
    7 CONTINUE                                                          RAY12690
C                                                                       RAY12700
C     NORMAL EXIT                                                       RAY12710
C                                                                       RAY12720
C     LIQUID SURFACE                                                    RAY12730
C                                                                       RAY12740
    8 IF( SV.LE.0 )  THEN                                               RAY12750
        Y0(1) = Y(1)                                                    RAY12760
        IF( ABS(Y(2))*EPS.LE.ABS(Y(1)) )  THEN                          RAY12770
          Y0(2) = Y(2)/ABS(Y(1))                                        RAY12780
        ELSE                                                            RAY12790
          Y0(2) = SIGN(BIG,Y(2))                                        RAY12800
        ENDIF                                                           RAY12810
C                                                                       RAY12820
        IF( IGG.GE.2 )  THEN                                            RAY12830
          EK =-WN*Y(4)/Y(1)                                             RAY12840
          IF( IGG.GT.2 )  U = C*Y(4)/(Y(4) + Y(6))                      RAY12850
        ENDIF                                                           RAY12860
C                                                                       RAY12870
C     SOLID SURFACE                                                     RAY12880
C                                                                       RAY12890
      ELSE                                                              RAY12900
        Y0(1) = Y(3)                                                    RAY12910
        IF( ABS(Y(5))*EPS.LE.ABS(Y(3)) )  THEN                          RAY12920
          Y0(2) = Y(5)/ABS(Y(3))                                        RAY12930
          Y0(3) =-Y(1)/Y(3)                                             RAY12940
        ELSE                                                            RAY12950
          Y0(2) = SIGN(BIG,Y(5))                                        RAY12960
        ENDIF                                                           RAY12970
C                                                                       RAY12980
        IF( IGG.GE.2 )  THEN                                            RAY12990
          EK =-WN*Y(10)/Y(3)                                            RAY13000
          IF( IGG.GT.2 )  U = C*Y(10)/(Y(10) + Y(15))                   RAY13010
        ENDIF                                                           RAY13020
      ENDIF                                                             RAY13030
      RETURN                                                            RAY13040
C                                                                       RAY13050
C     INPUT ERROR                                                       RAY13060
C                                                                       RAY13070
   90 WRITE(6,91)  L,W,C                                                RAY13080
   91   FORMAT(20X,5('?'),3X,'(RAYMRX)   INPUT ERROR',3X,               RAY13090
     *         'L  =',I5,3X,'W  =',1PE13.6,3X,                          RAY13100
     1         'C  =',E13.6,3X,5('?'))                                  RAY13110
      IER =-1                                                           RAY13120
C                                                                       RAY13130
C     DUMMY STATEMENTS                                                  RAY13140
C                                                                       RAY13150
      AP(1) = AP(1)                                                     RAY13160
      AE(1) = AE(1)                                                     RAY13170
      IA = IA                                                           RAY13180
      RETURN                                                            RAY13190
C                                                                       RAY13200
   92 WRITE(6,93)  I                                                    RAY13210
   93   FORMAT(20X,5('?'),3X,'(RAYMRX)   INPUT ERROR ',                 RAY13220
     *         'AT LAYER',I5,3X,5('?'))                                 RAY13230
      IER =-MAX(1,I)                                                    RAY13240
      RETURN                                                            RAY13250
      END                                                               RAY13260
C                                                                       RAY13270
C                                                                       RAY18380
C     RAYLEIGH WAVE EIGENFUNCTION AND PARTIAL DERIVATIVES               RAY18390
C                  BY MATRIX METHOD INTEGRATION                         RAY18400
C                                                                       RAY18410
      SUBROUTINE  RAYEFX(H,RHO,VP,VS,QP,QS,L,LY,LD,W,C,                 RAY18420
     *                   U,IQ,Q,EKD,EKI,ER,IY,Y,                        RAY18430
     1                   IP,P,NY,WS,IER)                                RAY18440
C                                                                       RAY18450
C     INPUT                                                             RAY18460
C       H    : LAYER THICKNESS.  H(L) IS ARBITRARY                      RAY18470
C       RHO  : DENSITY                                                  RAY18480
C       VP   : P WAVE VELOCITY                                          RAY18490
C       VS   : S WAVE VELOCITY                                          RAY18500
C       QP   : INVERSE OF P-WAVE Q                                      RAY18510
C       QS   : INVERSE OF S-WAVE Q                                      RAY18520
C       L    : NUMBER OF LAYERS INCLDING THE BOTTOM                     RAY18530
C              HALF-SPACE                                               RAY18540
C       LY   : DEEPEST LAYER WHERE EIGENFUNCTION HAS                    RAY18550
C              TO BE COMPUTED                                           RAY18560
C       LD   : NUMBER OF SUBLAYERS/LAYER                                RAY18570
C              ACTUALLY, 2*LD SUBLAYERS ARE INTRODUCED                  RAY18580
C       W    : ANGULAR FREQUENCY                                        RAY18590
C       C    : PHASE VELOCITY                                           RAY18600
C       EKD  : ENERGY INTEGRAL 2*K**2*I3                                RAY18610
C              COMPUTED BY RAYMRX                                       RAY18620
C       IY   : = 0 ; NOT TO COMPUTE EIGENFUNCTION                       RAY18630
C              = 1 ;     TO COMPUTE EIGENFUNCTION                       RAY18640
C       IP   : = 0 ; NOT TO COMPUTE PARTIALS                            RAY18650
C              = 1 ;     TO COMPUTE PARTIALS                            RAY18660
C       IQ   : = 0 ; NOT TO COMPUTE ATTENUATION                         RAY18670
C              = 1 ;     TO COMPUTE ATTENUATION                         RAY18680
C     OUTPUT                                                            RAY18690
C       U    : GROUP VELOCITY BY ENERGY INTEGRALS                       RAY18700
C       Q    : INVERSE OF RAYLEIGH WAVE Q                               RAY18710
C       EKI  : ENERGY INTEGRAL 2*K**2*I3 BY INTEGRATION                 RAY18720
C       ER   : RELATIVE ERROR IN ENERGY INTEGRALS                       RAY18730
C       Y    : EIGENFUNCTION                                            RAY18740
C       P    : PARTIALS                                                 RAY18750
C              P(1,J) ; DC/DRHO                                         RAY18760
C              P(2,J) ; DC/DVP                                          RAY18770
C              P(3,J) ; DC/DVS                                          RAY18780
C              P(4,J) ; DC/DPHI                                         RAY18790
C              P(5,J) ; DC/DETA                                         RAY18800
C       NY   : LENGTH OF Y AND/OR P                                     RAY18810
C       IER  : RETURN CODE                                              RAY18820
C              < 0 ; INPUT ERROR                                        RAY18830
C              = 0 ; NO ERROR                                           RAY18840
C              > 0 ; OVERFLOW                                           RAY18850
C     WORKING SPACE                                                     RAY18860
C       WS   : 5*L WORDS                                                RAY18870
C                                                                       RAY18880
C     ISOTROPIC MODEL ONLY                                              RAY18890
C                                                                       RAY18900
C     SUBPROGRAM : RAYENG                                               RAY18910
C                                                                       RAY18920
C     DISPER-80,  VER-86                                                RAY18930
C                                                                       RAY18940
C     M. SAITO  30/VII/79                                               RAY18950
C     REVISED   10/XII/86                                               RAY18960
C                                                                       RAY18970
      DIMENSION  H(L),RHO(L),VP(L),VS(L),QP(L),QS(L),                   RAY18980
     *           Y(4,*),P(5,*),WS(5,*),YY(5),Z(5),                      RAY18990
     1           F(4,3),S(4)                                            RAY19000
      DATA  EPS/1.E-10/,BIG/1.E+10/                                     RAY19010
C                                                                       RAY19020
C     DEFINE SINH(X)/X                                                  RAY19030
C                                                                       RAY19040
      SH0(X) = 0.9999997 + X*(0.1666667 + X*(0.0083361                  RAY19050
     *       + X*0.0001984))                                            RAY19060
C                                                                       RAY19070
C     FOR DOUBLE PRECISION USE                                          RAY19080
C                                                                       RAY19090
C     SH0(X) = 1.0D0 + X*(1.6666 66666 66666 7D-1                       RAY19100
C    *       + X*(8.33 33333 33334 0D-3                                 RAY19110
C    1       + X*(1.9 84126 96412 7D-4                                  RAY19120
C    2       + X*(2.7557 31918 9D-6 + X*(2.50 12108 4D-8                RAY19130
C    3       + X*(1.60596 1D-10 + X*7.64 7D-13))))))                    RAY19140
C                                                                       RAY19150
C     INITIAL VALUE                                                     RAY19160
C                                                                       RAY19170
      IF( L.LE.1 .OR. C.LE.0 .OR. W.LE.0 )  GO TO  90                   RAY19180
      I   = L                                                           RAY19190
      RO  = RHO(I)                                                      RAY19200
      IF( RO.LE.0 .OR. C.GE.VP(I) )  GO TO  92                          RAY19210
      IER = 0                                                           RAY19220
      NY  = 0                                                           RAY19230
      LLY = MAX(1,MIN0(L,LY))                                           RAY19240
      CC  = C*C                                                         RAY19250
      WN  = W/C                                                         RAY19260
      ROC = RO*CC                                                       RAY19270
      SV  = VS(I)                                                       RAY19280
      CP  = C/VP(I)                                                     RAY19290
      RAA = (1 + CP)*(1 - CP)                                           RAY19300
      RA  = SQRT(RAA)                                                   RAY19310
      DO  1  J=1,5                                                      RAY19320
    1   YY(J) = 0                                                       RAY19330
C                                                                       RAY19340
      IF( SV.LE.0 )  THEN                                               RAY19350
        YY(1) = RA                                                      RAY19360
        YY(2) =-ROC                                                     RAY19370
        JX = 2                                                          RAY19380
C                                                                       RAY19390
      ELSE                                                              RAY19400
        IF( C.GE.SV )  GO TO  92                                        RAY19410
        CS  = C/SV                                                      RAY19420
        RBB = (1 + CS)*(1 - CS)                                         RAY19430
        RB  = SQRT(RBB)                                                 RAY19440
        RG  = 2*ROC/CS**2                                               RAY19450
        YY(3) =-EPS*RA                                                  RAY19460
        YY(4) =-EPS*RB                                                  RAY19470
        YY(2) =-EPS*(CP**2*RBB + CS**2)/(ROC*(RA*RB + 1))               RAY19480
        YY(1) = RG*YY(2) + EPS                                          RAY19490
        YY(5) =-RG*(YY(1) + EPS) + EPS*ROC                              RAY19500
        JX = 5                                                          RAY19510
      ENDIF                                                             RAY19520
C                                                                       RAY19530
      DO  2  J=1,5                                                      RAY19540
    2   WS(J,I) = YY(J)                                                 RAY19550
C                                                                       RAY19560
C     INTEGRATE UPWARD                                                  RAY19570
C                                                                       RAY19580
      IF( L.LE.1 )  GO TO  17                                           RAY19590
      DO  9  II=2,L                                                     RAY19600
        I  = I - 1                                                      RAY19610
        RO = RHO(I)                                                     RAY19620
        ROC= RO*CC                                                      RAY19630
        PV = VP(I)                                                      RAY19640
        SV = VS(I)                                                      RAY19650
        IF( PV.LE.0 .OR. RO.LE.0 .OR.                                   RAY19660
     *      H(I).LE.0 )  GO TO  92                                      RAY19670
        DO  3  J=1,5                                                    RAY19680
    3     Z(J) = YY(J)                                                  RAY19690
        IF( (SV.LE.0 .AND. VS(I+1).LE.0) .OR.                           RAY19700
     *      (SV.GT.0 .AND. VS(I+1).GT.0) )  GO TO  4                    RAY19710
C                                                                       RAY19720
C     SOLID  -----> LIQUID                                              RAY19730
C                                                                       RAY19740
        IF( SV.LE.0 )  THEN                                             RAY19750
          Z(1) = YY(3)                                                  RAY19760
          Z(2) = YY(5)                                                  RAY19770
          JX = 2                                                        RAY19780
C                                                                       RAY19790
C     LIQUID -----> SOLID                                               RAY19800
C                                                                       RAY19810
        ELSE                                                            RAY19820
          Z(2) = YY(1)                                                  RAY19830
          Z(4) = YY(2)                                                  RAY19840
          Z(1) = 0                                                      RAY19850
          Z(3) = 0                                                      RAY19860
          Z(5) = 0                                                      RAY19870
          JX = 5                                                        RAY19880
        ENDIF                                                           RAY19890
C                                                                       RAY19900
    4   R2  = 1/ROC                                                     RAY19910
        CP  = C/PV                                                      RAY19920
        RAA = (1 + CP)*(1 - CP)                                         RAY19930
        CS  = 0                                                         RAY19940
        IF( SV.GT.0 )  CS = C/SV                                        RAY19950
        RBB = (1 + CS)*(1 - CS)                                         RAY19960
        HK  = H(I)*WN                                                   RAY19970
        HKK = HK**2                                                     RAY19980
        XX  = RAA*HKK                                                   RAY19990
        ONE = 1                                                         RAY20000
C                                                                       RAY20010
C       SINH(X)/X                                                       RAY20020
C                                                                       RAY20030
        DO  5  K=1,2                                                    RAY20040
          CHA = CHB                                                     RAY20050
          SHA = SHB                                                     RAY20060
          AX  = ABS(XX)                                                 RAY20070
          IF( AX.LE.1 )  THEN                                           RAY20080
            SHB = SH0(XX)                                               RAY20090
            CHB = 1 + XX*SH0(XX/4)**2/2                                 RAY20100
          ELSE                                                          RAY20110
            AX  = SQRT(AX)                                              RAY20120
            IF( XX.LE.0 )  THEN                                         RAY20130
              CHB = COS(AX)                                             RAY20140
              SHB = SIN(AX)/AX                                          RAY20150
            ELSE                                                        RAY20160
              IF( AX.GT.100 )  THEN                                     RAY20170
                ONE = 0                                                 RAY20180
              ELSE                                                      RAY20190
                ONE = ONE/COSH(AX)                                      RAY20200
              ENDIF                                                     RAY20210
              CHB = 1                                                   RAY20220
              SHB = TANH(AX)/AX                                         RAY20230
            ENDIF                                                       RAY20240
          ENDIF                                                         RAY20250
          XX  = HKK*RBB                                                 RAY20260
          SHB = HK*SHB                                                  RAY20270
    5   CONTINUE                                                        RAY20280
C                                                                       RAY20290
C     LIQUID LAYER                                                      RAY20300
C                                                                       RAY20310
        IF( SV.LE.0 )  THEN                                             RAY20320
          B11 = CHA                                                     RAY20330
          B12 =-RAA*SHA*R2                                              RAY20340
          B21 =-ROC*SHA                                                 RAY20350
C                                                                       RAY20360
          YY(1) = B11*Z(1) + B12*Z(2)                                   RAY20370
          YY(2) = B21*Z(1) + B11*Z(2)                                   RAY20380
C                                                                       RAY20390
C     SOLID LAYER                                                       RAY20400
C                                                                       RAY20410
        ELSE                                                            RAY20420
          G1  = 2/CS**2                                                 RAY20430
          RG  = G1*ROC                                                  RAY20440
          R4  = RG - ROC                                                RAY20450
          E1  = CHA*CHB                                                 RAY20460
          E2  = E1 - ONE                                                RAY20470
          E3  = SHA*SHB                                                 RAY20480
          E5  = SHA*CHB                                                 RAY20490
          E6  = SHB*CHA                                                 RAY20500
          F1  = E2 - E3                                                 RAY20510
          F2  = R2*F1                                                   RAY20520
          F3  = G1*F1 + E3                                              RAY20530
          B33 = E1                                                      RAY20540
          B34 = RAA*E3                                                  RAY20550
          B43 = RBB*E3                                                  RAY20560
          B25 =-R2*(F2 + R2*(E2 - RAA*B43))                             RAY20570
          B15 = RG*B25 + F2                                             RAY20580
          B16 =-RG*B15 - F3                                             RAY20590
          B22 = B16 + E1                                                RAY20600
          B12 = RG*B16 - R4*F3                                          RAY20610
          B52 =-RG*B12 + R4*(RG*F3 + R4*E3)                             RAY20620
          B23 = R2*(E5 - RBB*E6)                                        RAY20630
          B13 = RG*B23 - E5                                             RAY20640
          B42 =-RG*B13 + R4*E5                                          RAY20650
          B24 = R2*(E6 - RAA*E5)                                        RAY20660
          B14 = RG*B24 - E6                                             RAY20670
          B32 =-RG*B14 + R4*E6                                          RAY20680
          B11 = ONE - B16-B16                                           RAY20690
          B21 = B15+B15                                                 RAY20700
          B31 = B14+B14                                                 RAY20710
          B41 = B13+B13                                                 RAY20720
          B51 = B12+B12                                                 RAY20730
C                                                                       RAY20740
          YY(1) = B11*Z(1) + B12*Z(2) + B13*Z(3)                        RAY20750
     *          + B14*Z(4) + B15*Z(5)                                   RAY20760
          YY(2) = B21*Z(1) + B22*Z(2) + B23*Z(3)                        RAY20770
     *          + B24*Z(4) + B25*Z(5)                                   RAY20780
          YY(3) = B31*Z(1) + B32*Z(2) + B33*Z(3)                        RAY20790
     *          + B34*Z(4) + B24*Z(5)                                   RAY20800
          YY(4) = B41*Z(1) + B42*Z(2) + B43*Z(3)                        RAY20810
     *          + B33*Z(4) + B23*Z(5)                                   RAY20820
          YY(5) = B51*Z(1) + B52*Z(2) + B42*Z(3)                        RAY20830
     *          + B32*Z(4) + B22*Z(5)                                   RAY20840
        ENDIF                                                           RAY20850
C                                                                       RAY20860
        DO  6  J=1,5                                                    RAY20870
    6     WS(J,I) = YY(J)                                               RAY20880
C                                                                       RAY20890
C     NORMALIZATION                                                     RAY20900
C                                                                       RAY20910
        IF( EL.GT.0 )  THEN                                             RAY20920
          Z1 = 0                                                        RAY20930
          DO  7  J=1,JX                                                 RAY20940
    7       Z1 = MAX(Z1,ABS(YY(J)))                                     RAY20950
          IF( Z1.GT.BIG )  THEN                                         RAY20960
            DO  8  J=1,JX                                               RAY20970
    8         YY(J) = EPS*YY(J)                                         RAY20980
          ENDIF                                                         RAY20990
        ENDIF                                                           RAY21000
C                                                                       RAY21010
    9 CONTINUE                                                          RAY21020
C                                                                       RAY21030
C     INTEGRATE DOWNWARD                                                RAY21040
C                                                                       RAY21050
      DO  10  J=1,4                                                     RAY21060
          S(J) = 0                                                      RAY21070
        F(J,1) = 0                                                      RAY21080
        F(J,2) = 0                                                      RAY21090
        F(J,3) = 0                                                      RAY21100
   10 CONTINUE                                                          RAY21110
      LX  = MAX(1,LD)                                                   RAY21120
      DIV = 2*LX                                                        RAY21130
      IF( SV.GT.0 )  YY(3) =-YY(1)/YY(3)                                RAY21140
      YY(1) = 1                                                         RAY21150
      YY(2) = 0                                                         RAY21160
      YY(4) = 0                                                         RAY21170
C                                                                       RAY21180
      IF( LLY.LE.1 )  GO TO  17                                         RAY21190
      DO  16  II=2,LLY                                                  RAY21200
        RO  = RHO(I)                                                    RAY21210
        ROC = RO*CC                                                     RAY21220
        PV  = VP(I)                                                     RAY21230
        SV  = VS(I)                                                     RAY21240
        HH  = H(I)/DIV                                                  RAY21250
        EA  = RO*PV**2                                                  RAY21260
        EL  = RO*SV**2                                                  RAY21270
        EF  = EA - EL-EL                                                RAY21280
        R2  = 1/ROC                                                     RAY21290
        CP  = C/PV                                                      RAY21300
        RAA = (1 + CP)*(1 - CP)                                         RAY21310
        CS  = 0                                                         RAY21320
        IF( SV.GT.0 )  CS = C/SV                                        RAY21330
        RBB = (1 + CS)*(1 - CS)                                         RAY21340
        HK  = HH*WN                                                     RAY21350
        HKK = HK**2                                                     RAY21360
        XX  = RAA*HKK                                                   RAY21370
C                                                                       RAY21380
        DO  11  K=1,2                                                   RAY21390
          CHA = CHB                                                     RAY21400
          SHA = SHB                                                     RAY21410
          AX  = ABS(XX)                                                 RAY21420
          IF( AX.LE.1 )  THEN                                           RAY21430
            SHB = SH0(XX)                                               RAY21440
            CHB = 1 + XX*SH0(XX/4)**2/2                                 RAY21450
          ELSE                                                          RAY21460
            AX  = SQRT(AX)                                              RAY21470
            IF( XX.LE.0 )  THEN                                         RAY21480
              CHB = COS(AX)                                             RAY21490
              SHB = SIN(AX)/AX                                          RAY21500
            ELSE                                                        RAY21510
              IF( AX.GT.100 )  GO TO  94                                RAY21520
              EX  = EXP(AX)                                             RAY21530
              EXI = 1/EX                                                RAY21540
              CHB = (EX + EXI)/2                                        RAY21550
              SHB = (EX - EXI)/(2*AX)                                   RAY21560
            ENDIF                                                       RAY21570
          ENDIF                                                         RAY21580
          XX = HKK*RBB                                                  RAY21590
          SHB =-HK*SHB                                                  RAY21600
   11   CONTINUE                                                        RAY21610
C                                                                       RAY21620
C     LAYER MATRIX                                                      RAY21630
C                                                                       RAY21640
C       LIQUID LAYER                                                    RAY21650
C                                                                       RAY21660
        IF( SV.LE.0 )  THEN                                             RAY21670
          B11 = CHA                                                     RAY21680
          B12 =-RAA*SHA*R2                                              RAY21690
          B21 =-ROC*SHA                                                 RAY21700
C                                                                       RAY21710
C       SOLID LAYER                                                     RAY21720
C                                                                       RAY21730
        ELSE                                                            RAY21740
          RG  = EL+EL                                                   RAY21750
          R4  = RG - ROC                                                RAY21760
          B14 = (CHA - CHB)*R2                                          RAY21770
          B11 =-RG*B14 + CHA                                            RAY21780
          B33 = RG*B14 + CHB                                            RAY21790
          B23 = RG*(B33 - CHA)                                          RAY21800
          B34 = (SHA - RBB*SHB)*R2                                      RAY21810
          B24 = RG*B34 - SHA                                            RAY21820
          B21 =-RG*B24 + R4*SHA                                         RAY21830
          B12 = (SHB - RAA*SHA)*R2                                      RAY21840
          B13 =-RG*B12 + SHB                                            RAY21850
          B43 = RG*B13 + R4*SHB                                         RAY21860
        ENDIF                                                           RAY21870
C                                                                       RAY21880
        NY = NY + 1                                                     RAY21890
        CALL  RAYENG(RO,EA,EA,EL,EF,QP(I),QS(I),C,                      RAY21900
     *               YY,F,IY,Y(1,NY),IP,P(1,NY),IQ)                     RAY21910
C                                                                       RAY21920
        DO  15  LL=1,LX                                                 RAY21930
C                                                                       RAY21940
          DO  13  K=2,3                                                 RAY21950
            DO  12  J=1,4                                               RAY21960
   12         Z(J) = YY(J)                                              RAY21970
C                                                                       RAY21980
            IF( SV.LE.0 )  THEN                                         RAY21990
              YY(1) = B11*Z(1) + B12*Z(2)                               RAY22000
              YY(2) = B21*Z(1) + B11*Z(2)                               RAY22010
C                                                                       RAY22020
            ELSE                                                        RAY22030
              YY(1) = B11*Z(1) + B12*Z(2) + B13*Z(3)                    RAY22040
     *              + B14*Z(4)                                          RAY22050
              YY(2) = B21*Z(1) + B11*Z(2) + B23*Z(3)                    RAY22060
     *              + B24*Z(4)                                          RAY22070
              YY(3) =-B24*Z(1) - B14*Z(2) + B33*Z(3)                    RAY22080
     *              + B34*Z(4)                                          RAY22090
              YY(4) =-B23*Z(1) - B13*Z(2) + B43*Z(3)                    RAY22100
     *              + B33*Z(4)                                          RAY22110
            ENDIF                                                       RAY22120
C                                                                       RAY22130
            NY = NY + 1                                                 RAY22140
            CALL  RAYENG(RO,EA,EA,EL,EF,QP(I),QS(I),C,YY,               RAY22150
     *                   F(1,K),IY,Y(1,NY),IP,P(1,NY),IQ)               RAY22160
C                                                                       RAY22170
   13     CONTINUE                                                      RAY22180
C                                                                       RAY22190
          H6 = HH/3                                                     RAY22200
          DO  14  J=1,4                                                 RAY22210
              S(J) = S(J)                                               RAY22220
     *             + H6*(F(J,1) + 4*F(J,2) + F(J,3))                    RAY22230
            F(J,1) = F(J,3)                                             RAY22240
   14     CONTINUE                                                      RAY22250
C                                                                       RAY22260
   15   CONTINUE                                                        RAY22270
C                                                                       RAY22280
        I = I + 1                                                       RAY22290
C                                                                       RAY22300
        IF( VS(I).LE.0 )  THEN                                          RAY22310
          YY(2) = (WS(2,I)/WS(1,I))*YY(1)                               RAY22320
C                                                                       RAY22330
        ELSE                                                            RAY22340
          IF( SV.LE.0 )  YY(3) =-(WS(1,I)/WS(3,I))*YY(1)                RAY22350
          YY(2) = (WS(4,I)*YY(1) + WS(1,I)*YY(3))/WS(2,I)               RAY22360
          YY(4) = (WS(1,I)*YY(1) + WS(3,I)*YY(3))/WS(2,I)               RAY22370
        ENDIF                                                           RAY22380
C                                                                       RAY22390
   16 CONTINUE                                                          RAY22400
C                                                                       RAY22410
C     HALF SPACE INTEGRAL                                               RAY22420
C                                                                       RAY22430
   17 RO  = RHO(I)                                                      RAY22440
      ROC = RO*CC                                                       RAY22450
      EA  = RO*VP(I)**2                                                 RAY22460
      EL  = RO*VS(I)**2                                                 RAY22470
      EF  = EA - EL-EL                                                  RAY22480
      NY  = NY + 1                                                      RAY22490
      CALL  RAYENG(RO,EA,EA,EL,EF,QP(I),QS(I),C,                        RAY22500
     *             YY,F,IY,Y(1,NY),IP,P(1,NY),IQ)                       RAY22510
      CP  = C/VP(I)                                                     RAY22520
      RAA = (1 + CP)*(1 - CP)                                           RAY22530
      RA  = SQRT(RAA)                                                   RAY22540
C                                                                       RAY22550
C     LIQUID HALF SPACE                                                 RAY22560
C                                                                       RAY22570
      IF( VS(I).LE.0 )  THEN                                            RAY22580
        AA = (YY(1)/RA)**2/(2*WN*RA)                                    RAY22590
        S(1) = S(1) + RO*(1 + RAA)*AA                                   RAY22600
        S(2) = S(2) + ROC*CP**2*AA                                      RAY22610
        S(3) = S(3) + ROC*AA                                            RAY22620
        IF( IQ.NE.0 )  S(4) = S(4) + QP(I)*CP**2*AA                     RAY22630
C                                                                       RAY22640
C     SOLID HALF SPACE                                                  RAY22650
C                                                                       RAY22660
      ELSE                                                              RAY22670
        CS  = C/VS(I)                                                   RAY22680
        RBB = (1 + CS)*(1 - CS)                                         RAY22690
        RB  = SQRT(RBB)                                                 RAY22700
        RG  = EL+EL                                                     RAY22710
        CB  = (RA*RB + 1)/(CP**2 + CS**2*RAA)                           RAY22720
        CA  = (YY(3) - RB*YY(1))*CB                                     RAY22730
        CB  = (YY(1) - RA*YY(3))*CB                                     RAY22740
        AA  = CA**2/(2*WN*RA)                                           RAY22750
        AB  = CA*CB/(WN*(RA + RB))                                      RAY22760
        BB  = CB**2/(2*WN*RB)                                           RAY22770
        YA2 = RG - ROC                                                  RAY22780
        YA4 = RG*RA                                                     RAY22790
        YB2 = RG*RB                                                     RAY22800
        S11 = RAA*AA + 2*RA*AB + BB                                     RAY22810
        S22 = YA2**2*AA + 2*YA2*YB2*AB + YB2**2*BB                      RAY22820
        S33 = AA + 2*RB*AB + RBB*BB                                     RAY22830
        S44 = YA4**2*AA + 2*YA4*YA2*AB + YA2**2*BB                      RAY22840
        S14 = RA*YA4*AA + (RA*YA2 + YA4)*AB + YA2*BB                    RAY22850
        S23 = YA2*AA + (RB*YA2 + YB2)*AB + RB*YB2*BB                    RAY22860
        AF  = (EA - EF*EF/EA)*S33                                       RAY22870
        S(1) = S(1) + RO*(S11 + S33)                                    RAY22880
        S(2) = S(2) + S22/EA + S44/EL + AF                              RAY22890
        S(3) = S(3) + AF + S14 - EF*S23/EA                              RAY22900
        IF( IQ.NE.0 )  THEN                                             RAY22910
          DP = (S22 - 4*EL*(S23 - EL*S33))/EA                           RAY22920
          DS = S44/EL + 4*EL*(S23 + EF*S33)/EA                          RAY22930
          S(4) = S(4) + QP(I)*DP + QS(I)*DS                             RAY22940
        ENDIF                                                           RAY22950
      ENDIF                                                             RAY22960
C                                                                       RAY22970
      U   = S(3)/(C*S(1))                                               RAY22980
      ER  = 1 - S(2)/(CC*S(1))                                          RAY22990
      EKI = 2*WN*WN*S(3)                                                RAY23000
      IF( IQ.NE.0 )  Q = S(4)/S(3)                                      RAY23010
C                                                                       RAY23020
      IF( IP.NE.0 )  THEN                                               RAY23030
        IF( EKD.GT.0 )  THEN                                            RAY23040
          Z1 = WN**2/EKD                                                RAY23050
        ELSE                                                            RAY23060
          Z1 = 1/(2*S(3))                                               RAY23070
        ENDIF                                                           RAY23080
        DO  18  K=1,5                                                   RAY23090
        DO  18  J=1,NY                                                  RAY23100
          P(K,J) = Z1*P(K,J)                                            RAY23110
   18   CONTINUE                                                        RAY23120
      ENDIF                                                             RAY23130
C                                                                       RAY23140
      IF( IY.NE.0 )  THEN                                               RAY23150
        DO  19  J=1,NY                                                  RAY23160
          Y(2,J) = WN*Y(2,J)                                            RAY23170
          Y(4,J) = WN*Y(4,J)                                            RAY23180
   19   CONTINUE                                                        RAY23190
      ENDIF                                                             RAY23200
      RETURN                                                            RAY23210
C                                                                       RAY23220
C     INPUT ERROR                                                       RAY23230
C                                                                       RAY23240
   90 WRITE(6,91)  L,W,C                                                RAY23250
   91   FORMAT(20X,5('?'),3X,'(RAYEFX)   INPUT ERROR',3X,               RAY23260
     *         'L  =',I5,3X,'W  =',1PE13.6,3X,                          RAY23270
     1         'C  =',E13.6,3X,5('?'))                                  RAY23280
      IER =-1                                                           RAY23290
      RETURN                                                            RAY23300
C                                                                       RAY23310
   92 WRITE(6,93)  I                                                    RAY23320
   93   FORMAT(20X,5('?'),3X,'(RAYEFX)   INPUT ERROR ',                 RAY23330
     *         'AT LAYER',I5,3X,5('?'))                                 RAY23340
      IER =-MAX(1,I)                                                    RAY23350
      RETURN                                                            RAY23360
C                                                                       RAY23370
C     OVERFLOW                                                          RAY23380
C                                                                       RAY23390
   94 WRITE(6,95)  I                                                    RAY23400
   95   FORMAT(20X,5('?'),3X,'(RAYEFX)   OVERFLOW  ',                   RAY23410
     *         'AT LAYER',I5,3X,5('?'))                                 RAY23420
      IER = I                                                           RAY23430
      RETURN                                                            RAY23440
      END                                                               RAY23450
C                                                                       RAY23460
C     RAYLEIGH WAVE ENERGY INTEGRAL                                     RAY23470
C                                                                       RAY23480
      SUBROUTINE  RAYENG(RHO,EA,EC,EL,EF,QP,QS,C,YY,F,                  RAY23490
     *                   IY,Y,IP,P,IQ)                                  RAY23500
C                                                                       RAY23510
C     DISPER-80,  VER-86                                                RAY23520
C                                                                       RAY23530
C     M. SAITO  30/VII/79                                               RAY23540
C     REVISED   27/IX/85                                                RAY23550
C                                                                       RAY23560
      DIMENSION  YY(4),F(4),Y(4),P(5)                                   RAY23570
C                                                                       RAY23580
      CC = C*C                                                          RAY23590
      Y1 = YY(1)                                                        RAY23600
      Y2 = YY(2)                                                        RAY23610
      Y3 = YY(3)                                                        RAY23620
      Y4 = YY(4)                                                        RAY23630
C                                                                       RAY23640
      IF( EL.LE.0 )  THEN                                               RAY23650
        Y3 =-Y2/(RHO*CC)                                                RAY23660
        Y4 = 0                                                          RAY23670
      ENDIF                                                             RAY23680
C                                                                       RAY23690
      IF( IY.NE.0 )  THEN                                               RAY23700
        Y(1) = Y1                                                       RAY23710
        Y(2) = Y2                                                       RAY23720
        Y(3) = Y3                                                       RAY23730
        Y(4) = Y4                                                       RAY23740
      ENDIF                                                             RAY23750
C                                                                       RAY23760
      Y42  = 0                                                          RAY23770
      IF( EL.GT.0 )  Y42 = Y4**2/EL                                     RAY23780
      FC = EF/EC                                                        RAY23790
      F3 = (EA - EF*FC)*Y3**2                                           RAY23800
      F(1) = RHO*(Y1**2 + Y3**2)                                        RAY23810
      F(2) = Y2**2/EC + Y42 + F3                                        RAY23820
      F(3) = F3 + Y1*Y4 - FC*Y2*Y3                                      RAY23830
      IF( IP.EQ.0 .AND. IQ.EQ.0 )  RETURN                               RAY23840
C                                                                       RAY23850
      ET = EF/(EA - EL-EL)                                              RAY23860
      DPH= Y2 + EF*Y3                                                   RAY23870
      DP = (Y2 - 2*ET*EL*Y3)**2/EC                                      RAY23880
     *   + EA*(1 - (EA/EC)*ET**2)*Y3**2                                 RAY23890
      DS = Y42 + 4*ET*(EL/EC)*Y3*DPH                                    RAY23900
      IF( IQ.NE.0 )  F(4) = QP*DP + QS*DS                               RAY23910
C                                                                       RAY23920
      IF( IP.NE.0 )  THEN                                               RAY23930
        P(1) = F(2) - CC*F(1)                                           RAY23940
        P(2) = DP+DP                                                    RAY23950
        P(3) = DS+DS                                                    RAY23960
        P(4) = DPH**2/EC                                                RAY23970
        P(5) =-2*FC*Y3*DPH                                              RAY23980
      ENDIF                                                             RAY23990
C                                                                       RAY24000
      RETURN                                                            RAY24010
      END                                                               RAY24020
C                                                                       RAY24030
C     FIND GRID NUMBER FOR A GIVEN DEPTH                                RAY24040
C                                                                       RAY24050
      SUBROUTINE  DSPDEP(H,N,DEP,K,L)                                   RAY24060
C                                                                       RAY24070
C     INPUT                                                             RAY24080
C       H    : GRID INTERVAL                                            RAY24090
C       N    : NUMBER OF GRID POINTS                                    RAY24100
C       DEP  : DEPTH                                                    RAY24110
C       K    : SEARCH EVERY K GRID POINT                                RAY24120
C              = 2 ; TO FIND GRID POINT CONSISTENT                      RAY24130
C                    TO R-K-G SCHEME                                    RAY24140
C              = 4 ; TO FIND GRID POINT CONSISTENT                      RAY24150
C                    TO SIMPSON SCHEME                                  RAY24160
C                    (TO COMPUTE ENERGY INTEGRALS)                      RAY24170
C     OUTPUT                                                            RAY24180
C       L    : GRID NUMBER                                              RAY24190
C                                                                       RAY24200
C     DISPER-80,  VER-86                                                RAY24210
C                                                                       RAY24220
C     M. SAITO   28/I/78                                                RAY24230
C     REVISED     2/XII/86                                              RAY24240
C                                                                       RAY24250
      DIMENSION  H(N)                                                   RAY24260
C                                                                       RAY24270
      IF( K.LE.1 )  THEN                                                RAY24280
        KK = 1                                                          RAY24290
      ELSE IF( K.LE.3 )  THEN                                           RAY24300
        KK = 2                                                          RAY24310
      ELSE                                                              RAY24320
        KK = 4                                                          RAY24330
      ENDIF                                                             RAY24340
      DD = 0                                                            RAY24350
      I  = 1                                                            RAY24360
C                                                                       RAY24370
    1 IF( I.GE.N )  GO TO  2                                            RAY24380
      IF( H(I).EQ.0 )  I = I + 1                                        RAY24390
      IF( DEP-DD.LE.H(I)*0.1 )  GO TO  3                                RAY24400
      DD = DD + H(I)                                                    RAY24410
      IF( KK.GT.1 )  THEN                                               RAY24420
        DD = DD + H(I+1)                                                RAY24430
        IF( KK.EQ.4 )  DD = DD + H(I+2) + H(I+3)                        RAY24440
      ENDIF                                                             RAY24450
      I  = I + KK                                                       RAY24460
      GO TO  1                                                          RAY24470
C                                                                       RAY24480
    2 IF( I.GT.N )  I = I - KK                                          RAY24490
    3 L  = I                                                            RAY24500
      RETURN                                                            RAY24510
      END                                                               RAY24520
C                                                                       RAY24530
C     PRINT EIGENFUNCTION OR PARTIAL                                    RAY24840
C                                                                       RAY24850
      SUBROUTINE  DSPRNX(Y,K,NY,LD,H,D0)                                RAY24860
C                                                                       RAY24870
C     DISPER-80,  VER-86                                                RAY24880
C                                                                       RAY24890
C     M. SAITO  15/XI/81                                                RAY24900
C     REVISED    3/XII/86                                               RAY24910
C                                                                       RAY24920
      DIMENSION  Y(K,NY),H(1)                                           RAY24930
C                                                                       RAY24940
      J   = 0                                                           RAY24950
      I   = 0                                                           RAY24960
      D   = D0                                                          RAY24970
      LD2 = LD*2                                                        RAY24980
      DIV = LD2                                                         RAY24990
C                                                                       RAY25000
    1 J = J + 1                                                         RAY25010
      I = I + 1                                                         RAY25020
      WRITE(*,2)  I,D,(Y(M,J),M=1,K)                                    RAY25030
    2	FORMAT(I5,1P6E15.6)
      IF( J.GE.NY )  RETURN                                             RAY25060
      HH = H(I)/DIV                                                     RAY25070
      DO  3  L=1,LD2                                                    RAY25080
        D = D + HH                                                      RAY25090
        J = J + 1                                                       RAY25100
        WRITE(*,2)  I,D,(Y(M,J),M=1,K)                                  RAY25110
    3 CONTINUE                                                          RAY25120
      IF( J.LT.NY )  GO TO  1                                           RAY25130
C                                                                       RAY25140
      RETURN                                                            RAY25150
      END                                                               RAY25160
C                                                                       RAY25170


