C   IMSL ROUTINE NAME   - ZSPWA                                         ZSPR0010
C                                                                       ZSPR0020
C-----------------------------------------------------------------------ZSPR0030
C                                                                       ZSPR0040
C   COMPUTER            - VAX/SINGLE                                    ZSPR0050
C                                                                       ZSPR0060
C   LATEST REVISION     - JUNE 1, 1982                                  ZSPR0070
C                                                                       ZSPR0080
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSPOW  ZSPR0090
C                                                                       ZSPR0100
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         ZSPR0110
C                       - SINGLE/H36,H48,H60                            ZSPR0120
C                                                                       ZSPR0130
C   REQD. IMSL ROUTINES - SINGLE/VBLA=SNRM2,ZSPWB,ZSPWC,ZSPWD,ZSPWE,    ZSPR0140
C                           ZSPWF,ZSPWG                                 ZSPR0150
C                       - DOUBLE/VBLA=DNRM2,ZSPWB,ZSPWC,ZSPWD,ZSPWE,    ZSPR0160
C                           ZSPWF,ZSPWG                                 ZSPR0170
C                                                                       ZSPR0180
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           ZSPR0190
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      ZSPR0200
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  ZSPR0210
C                                                                       ZSPR0220
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.       ZSPR0230
C                                                                       ZSPR0240
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN ZSPR0250
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    ZSPR0260
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        ZSPR0270
C                                                                       ZSPR0280
C-----------------------------------------------------------------------ZSPR0290
C                                                                       ZSPR0300
      SUBROUTINE ZSPWA1(FCN1,N,X,FVEC,XTOL,MAXFEV,ML,MU,EPSFCN,DIAG,MODEZSPR0310
     *                  ,FACTOR,NPRINT,INFO,NFEV,FJAC,LDFJAC,R,LR,QTF,  ZSPR0320
     *                   WA1,WA2,WA3,WA4,PAR)                           ZSPR0330
C                                  SPECIFICATIONS FOR ARGUMENTS         ZSPR0340
      INTEGER            N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,LR ZSPR0350
      REAL               X(N),FVEC(N),XTOL,EPSFCN,DIAG(N),FACTOR,       ZSPR0360
     *                   FJAC(LDFJAC,N),R(LR),QTF(N),WA1(N),WA2(N),     ZSPR0370
     *                   WA3(N),WA4(N),PAR(N+1)                           ZSPR0380
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   ZSPR0390
      INTEGER            IFLAG,ITER,IWA(1),I,JM1,J,L,MSUM,NCFAIL,NCSUC, ZSPR0400
     *                   NSLOW1,NSLOW2                                  ZSPR0410
      REAL               ACTRED,DELTA,EPSMCH,FNORM1,FNORM,ONE,P0001,    ZSPR0420
     *                   P001,P1,P5,PNORM,PRERED,RATIO,SPMPAR,SUM,TEMP, ZSPR0430
     *                   XNORM,ZERO                                     ZSPR0440
      REAL               SNRM2                                          ZSPR0450
      LOGICAL            JEVAL,SING                                     ZSPR0460
      EXTERNAL           FCN1                                           ZSPR0470
c      DATA               SPMPAR /2.2204E-16/                            ZSPR0480
      DATA		         SPMPAR /1.1921E-07/
      DATA               ONE,P1,P5,P001,P0001,ZERO /1.0E0,1.0E-1,5.0E-1,ZSPR0490
     *                   1.0E-3,1.0E-4,0.0E0/                           ZSPR0500
C                                  EPSMCH IS THE MACHINE PRECISION.     ZSPR0510
C                                  FIRST EXECUTABLE STATEMENT           ZSPR0520
      EPSMCH = SPMPAR                                                   ZSPR0530
      INFO = 0                                                          ZSPR0540
      IFLAG = 0                                                         ZSPR0550
      NFEV = 0                                                          ZSPR0560
C                                  CHECK THE INPUT PARAMETERS FOR       ZSPR0570
C                                  ERRORS.                              ZSPR0580
      IF (N.LE.0 .OR. XTOL.LT.ZERO .OR. MAXFEV.LE.0 .OR. ML.LT.0 .OR.   ZSPR0590
     *MU.LT.0 .OR. FACTOR.LE.ZERO .OR. LDFJAC.LT.N .OR.                 ZSPR0600
     *LR.LT.(N*(N+1))/2) GO TO 150                                      ZSPR0610
      IF (MODE.NE.2) GO TO 10                                           ZSPR0620
      DO 5 J=1,N                                                        ZSPR0630
         IF (DIAG(J).LE.ZERO) GO TO 150                                 ZSPR0640
    5 CONTINUE                                                          ZSPR0650
   10 CONTINUE                                                          ZSPR0660
C                                  EVALUATE THE FUNCTION AT THE STARTINGZSPR0670
C                                  POINT AND CALCULATE ITS NORM.        ZSPR0680
      IFLAG = 1                                                         ZSPR0690
      CALL FCN1(X,FVEC,N,PAR)                                           ZSPR0700
      NFEV = 1                                                          ZSPR0710
      IF (IFLAG.LT.0) GO TO 150                                         ZSPR0720
      FNORM = SNRM2(N,FVEC,1)                                           ZSPR0730
C                                  DETERMINE THE NUMBER OF CALLS TO FCN ZSPR0740
C                                  NEEDED TO COMPUTE THE JACOBIAN       ZSPR0750
C                                  MATRIX.                              ZSPR0760
C                                                                       ZSPR0770
      MSUM = min0(ML+MU+1,N)                                            ZSPR0780
C                                                                       ZSPR0790
C                                  INITIALIZE ITERATION COUNTER AND     ZSPR0800
C                                  MONITORS.                            ZSPR0810
      ITER = 1                                                          ZSPR0820
      NCSUC = 0                                                         ZSPR0830
      NCFAIL = 0                                                        ZSPR0840
      NSLOW1 = 0                                                        ZSPR0850
      NSLOW2 = 0                                                        ZSPR0860
C                                  BEGINNING OF THE OUTER LOOP.         ZSPR0870
   15 CONTINUE                                                          ZSPR0880
      JEVAL = .TRUE.                                                    ZSPR0890
C                                  CALCULATE THE JACOBIAN MATRIX.       ZSPR0900
      IFLAG = 2                                                         ZSPR0910
      CALL ZSPWB1(FCN1,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,WA1,WA2, ZSPR0920
     *PAR)                                                              ZSPR0930
      NFEV = NFEV+MSUM                                                  ZSPR0940
      IF (IFLAG.LT.0) GO TO 150                                         ZSPR0950
C                                  COMPUTE THE QR FACTORIZATION OF THE  ZSPR0960
C                                  JACOBIAN.                            ZSPR0970
      CALL ZSPWG(N,N,FJAC,LDFJAC,.FALSE.,IWA,1,WA1,WA2,WA3)             ZSPR0980
C                                  ON THE FIRST ITERATION AND IF MODE ISZSPR0990
C                                  1, SCALE ACCORDING TO THE NORMS OF   ZSPR1000
C                                  THE COLUMNS OF THE INITIAL JACOBIAN. ZSPR1010
      IF (ITER.NE.1) GO TO 35                                           ZSPR1020
      IF (MODE.EQ.2) GO TO 25                                           ZSPR1030
      DO 20 J=1,N                                                       ZSPR1040
         DIAG(J) = WA2(J)                                               ZSPR1050
         IF (WA2(J).EQ.ZERO) DIAG(J) = ONE                              ZSPR1060
   20 CONTINUE                                                          ZSPR1070
   25 CONTINUE                                                          ZSPR1080
C                                  ON THE FIRST ITERATION, CALCULATE THEZSPR1090
C                                  NORM OF THE SCALED X AND INITIALIZE  ZSPR1100
C                                  THE STEP BOUND DELTA.                ZSPR1110
      DO 30 J=1,N                                                       ZSPR1120
         WA3(J) = DIAG(J)*X(J)                                          ZSPR1130
   30 CONTINUE                                                          ZSPR1140
      XNORM = SNRM2(N,WA3,1)                                            ZSPR1150
      DELTA = FACTOR*XNORM                                              ZSPR1160
      IF (DELTA.EQ.ZERO) DELTA = FACTOR                                 ZSPR1170
   35 CONTINUE                                                          ZSPR1180
C                                  FORM (Q TRANSPOSE)*FVEC AND STORE IN ZSPR1190
C                                  QTF.                                 ZSPR1200
      DO 40 I=1,N                                                       ZSPR1210
         QTF(I) = FVEC(I)                                               ZSPR1220
   40 CONTINUE                                                          ZSPR1230
      DO 60 J=1,N                                                       ZSPR1240
         IF (FJAC(J,J).EQ.ZERO) GO TO 55                                ZSPR1250
         SUM = ZERO                                                     ZSPR1260
         DO 45 I=J,N                                                    ZSPR1270
            SUM = SUM+FJAC(I,J)*QTF(I)                                  ZSPR1280
   45    CONTINUE                                                       ZSPR1290
         TEMP = -SUM/FJAC(J,J)                                          ZSPR1300
         DO 50 I=J,N                                                    ZSPR1310
            QTF(I) = QTF(I)+FJAC(I,J)*TEMP                              ZSPR1320
   50    CONTINUE                                                       ZSPR1330
   55    CONTINUE                                                       ZSPR1340
   60 CONTINUE                                                          ZSPR1350
C                                  COPY THE TRIANGULAR FACTOR OF THE QR ZSPR1360
C                                  FACTORIZATION INTO R.                ZSPR1370
      SING = .FALSE.                                                    ZSPR1380
      DO 75 J=1,N                                                       ZSPR1390
         L = J                                                          ZSPR1400
         JM1 = J-1                                                      ZSPR1410
         IF (JM1.LT.1) GO TO 70                                         ZSPR1420
         DO 65 I=1,JM1                                                  ZSPR1430
            R(L) = FJAC(I,J)                                            ZSPR1440
            L = L+N-I                                                   ZSPR1450
   65    CONTINUE                                                       ZSPR1460
   70    CONTINUE                                                       ZSPR1470
         R(L) = WA1(J)                                                  ZSPR1480
         IF (WA1(J).EQ.ZERO) SING = .TRUE.                              ZSPR1490
   75 CONTINUE                                                          ZSPR1500
C                                  ACCUMULATE THE ORTHOGONAL FACTOR IN  ZSPR1510
C                                  FJAC.                                ZSPR1520
      CALL ZSPWF(N,N,FJAC,LDFJAC,WA1)                                   ZSPR1530
C                                  RESCALE IF NECESSARY.                ZSPR1540
      IF (MODE.EQ.2) GO TO 85                                           ZSPR1550
      DO 80 J=1,N                                                       ZSPR1560
         DIAG(J) = amax1(DIAG(J),WA2(J))                                ZSPR1570
   80 CONTINUE                                                          ZSPR1580
   85 CONTINUE                                                          ZSPR1590
C                                  BEGINNING OF THE INNER LOOP.         ZSPR1600
   90 CONTINUE                                                          ZSPR1610
C                                  IF REQUESTED, CALL FCN TO ENABLE     ZSPR1620
C                                  PRINTING OF ITERATES.                ZSPR1630
      IF (NPRINT.LE.0) GO TO 95                                         ZSPR1640
      IFLAG = 0                                                         ZSPR1650
      IF (IFLAG.LT.0) GO TO 150                                         ZSPR1660
   95 CONTINUE                                                          ZSPR1670
C                                  DETERMINE THE DIRECTION P.           ZSPR1680
      CALL ZSPWC(N,R,LR,DIAG,QTF,DELTA,WA1,WA2,WA3)                     ZSPR1690
C                                  STORE THE DIRECTION P AND X + P.     ZSPR1700
C                                  CALCULATE THE NORM OF P.             ZSPR1710
      DO 100 J=1,N                                                      ZSPR1720
         WA1(J) = -WA1(J)                                               ZSPR1730
         WA2(J) = X(J)+WA1(J)                                           ZSPR1740
         WA3(J) = DIAG(J)*WA1(J)                                        ZSPR1750
  100 CONTINUE                                                          ZSPR1760
      PNORM = SNRM2(N,WA3,1)                                            ZSPR1770
C                                  ON THE FIRST ITERATION, ADJUST THE   ZSPR1780
C                                  INITIAL STEP BOUND.                  ZSPR1790
      IF (ITER.EQ.1) DELTA = amin1(DELTA,PNORM)                         ZSPR1800
C                                  EVALUATE THE FUNCTION AT X + P AND   ZSPR1810
C                                  CALCULATE ITS NORM.                  ZSPR1820
      IFLAG = 1                                                         ZSPR1830
      CALL FCN1(WA2,WA4,N,PAR)                                          ZSPR1840
      NFEV = NFEV+1                                                     ZSPR1850
      IF (IFLAG.LT.0) GO TO 150                                         ZSPR1860
      FNORM1 = SNRM2(N,WA4,1)                                           ZSPR1870
C                                  COMPUTE THE SCALED ACTUAL REDUCTION. ZSPR1880
      ACTRED = -ONE                                                     ZSPR1890
      IF (FNORM1.LT.FNORM) ACTRED = ONE-(FNORM1/FNORM)**2               ZSPR1900
C                                  COMPUTE THE SCALED PREDICTED         ZSPR1910
C                                  REDUCTION.                           ZSPR1920
      L = 1                                                             ZSPR1930
      DO 110 I=1,N                                                      ZSPR1940
         SUM = ZERO                                                     ZSPR1950
         DO 105 J=I,N                                                   ZSPR1960
            SUM = SUM+R(L)*WA1(J)                                       ZSPR1970
            L = L+1                                                     ZSPR1980
  105    CONTINUE                                                       ZSPR1990
         WA3(I) = QTF(I)+SUM                                            ZSPR2000
  110 CONTINUE                                                          ZSPR2010
      TEMP = SNRM2(N,WA3,1)                                             ZSPR2020
      PRERED = ONE                                                      ZSPR2030
      IF (TEMP.LT.FNORM) PRERED = ONE-(TEMP/FNORM)**2                   ZSPR2040
C                                  COMPUTE THE RATIO OF THE ACTUAL TO   ZSPR2050
C                                  THE PREDICTED REDUCTION.             ZSPR2060
      RATIO = ZERO                                                      ZSPR2070
      IF (PRERED.GT.ZERO) RATIO = ACTRED/PRERED                         ZSPR2080
C                                  UPDATE THE STEP BOUND.               ZSPR2090
      IF (RATIO.GE.P1) GO TO 115                                        ZSPR2100
      NCSUC = 0                                                         ZSPR2110
      NCFAIL = NCFAIL+1                                                 ZSPR2120
      DELTA = P5*DELTA                                                  ZSPR2130
      GO TO 120                                                         ZSPR2140
  115 CONTINUE                                                          ZSPR2150
      NCFAIL = 0                                                        ZSPR2160
      NCSUC = NCSUC+1                                                   ZSPR2170
      IF (RATIO.GE.P5 .OR. NCSUC.GT.1) DELTA = amax1(DELTA,PNORM/P5)    ZSPR2180
      IF (abs(RATIO-ONE).LE.P1) DELTA = PNORM/P5                        ZSPR2190
  120 CONTINUE                                                          ZSPR2200
C                                  TEST FOR SUCCESSFUL ITERATION.       ZSPR2210
      IF (RATIO.LT.P0001) GO TO 130                                     ZSPR2220
C                                  SUCCESSFUL ITERATION. UPDATE X, FVEC,ZSPR2230
C                                  AND THEIR NORMS.                     ZSPR2240
      DO 125 J=1,N                                                      ZSPR2250
         X(J) = WA2(J)                                                  ZSPR2260
         WA2(J) = DIAG(J)*X(J)                                          ZSPR2270
         FVEC(J) = WA4(J)                                               ZSPR2280
  125 CONTINUE                                                          ZSPR2290
      XNORM = SNRM2(N,WA2,1)                                            ZSPR2300
      FNORM = FNORM1                                                    ZSPR2310
      ITER = ITER+1                                                     ZSPR2320
  130 CONTINUE                                                          ZSPR2330
C                                  DETERMINE THE PROGRESS OF THE        ZSPR2340
C                                  ITERATION.                           ZSPR2350
      NSLOW1 = NSLOW1+1                                                 ZSPR2360
      IF (ACTRED.GE.P001) NSLOW1 = 0                                    ZSPR2370
      IF (JEVAL) NSLOW2 = NSLOW2+1                                      ZSPR2380
      IF (ACTRED.GE.P1) NSLOW2 = 0                                      ZSPR2390
C                                  TEST FOR CONVERGENCE.                ZSPR2400
      IF (DELTA.LE.XTOL*XNORM .OR. FNORM.EQ.ZERO) INFO = 1              ZSPR2410
      IF (INFO.NE.0) GO TO 150                                          ZSPR2420
C                                  TESTS FOR TERMINATION AND STRINGENT  ZSPR2430
C                                  TOLERANCES.                          ZSPR2440
      IF (NFEV.GE.MAXFEV) INFO = 2                                      ZSPR2450
      IF (P1*amax1(P1*DELTA,PNORM).LE.EPSMCH*XNORM) INFO = 3            ZSPR2460
      IF (NSLOW2.EQ.5) INFO = 4                                         ZSPR2470
      IF (NSLOW1.EQ.10) INFO = 5                                        ZSPR2480
      IF (INFO.NE.0) GO TO 150                                          ZSPR2490
C                                  CRITERION FOR RECALCULATING JACOBIAN ZSPR2500
C                                  APPROXIMATION BY FORWARD DIFFERENCES.ZSPR2510
      IF (NCFAIL.EQ.2) GO TO 145                                        ZSPR2520
C                                  CALCULATE THE RANK ONE MODIFICATION  ZSPR2530
C                                  TO THE JACOBIAN AND UPDATE QTF IF    ZSPR2540
C                                  NECESSARY.                           ZSPR2550
      DO 140 J=1,N                                                      ZSPR2560
         SUM = ZERO                                                     ZSPR2570
         DO 135 I=1,N                                                   ZSPR2580
            SUM = SUM+FJAC(I,J)*WA4(I)                                  ZSPR2590
  135    CONTINUE                                                       ZSPR2600
         WA2(J) = (SUM-WA3(J))/PNORM                                    ZSPR2610
         WA1(J) = DIAG(J)*((DIAG(J)*WA1(J))/PNORM)                      ZSPR2620
         IF (RATIO.GE.P0001) QTF(J) = SUM                               ZSPR2630
  140 CONTINUE                                                          ZSPR2640
C                                  COMPUTE THE QR FACTORIZATION OF THE  ZSPR2650
C                                  UPDATED JACOBIAN.                    ZSPR2660
      CALL ZSPWE(N,N,R,LR,WA1,WA2,WA3,SING)                             ZSPR2670
      CALL ZSPWD(N,N,FJAC,LDFJAC,WA2,WA3)                               ZSPR2680
      CALL ZSPWD(1,N,QTF,1,WA2,WA3)                                     ZSPR2690
C                                  END OF THE INNER LOOP.               ZSPR2700
      JEVAL = .FALSE.                                                   ZSPR2710
      GO TO 90                                                          ZSPR2720
  145 CONTINUE                                                          ZSPR2730
C                                  END OF THE OUTER LOOP.               ZSPR2740
      GO TO 15                                                          ZSPR2750
  150 CONTINUE                                                          ZSPR2760
C                                  TERMINATION, EITHER NORMAL OR USER   ZSPR2770
C                                  IMPOSED.                             ZSPR2780
      IF (IFLAG.LT.0) INFO = IFLAG                                      ZSPR2790
      IFLAG = 0                                                         ZSPR2800
      RETURN                                                            ZSPR2810
      END                                                               ZSPR2820
