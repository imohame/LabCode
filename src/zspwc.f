C   IMSL ROUTINE NAME   - ZSPWC                                         ZSPT0010
C                                                                       ZSPT0020
C-----------------------------------------------------------------------ZSPT0030
C                                                                       ZSPT0040
C   COMPUTER            - VAX/SINGLE                                    ZSPT0050
C                                                                       ZSPT0060
C   LATEST REVISION     - JUNE 1, 1982                                  ZSPT0070
C                                                                       ZSPT0080
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSPOW  ZSPT0090
C                                                                       ZSPT0100
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         ZSPT0110
C                       - SINGLE/H36,H48,H60                            ZSPT0120
C                                                                       ZSPT0130
C   REQD. IMSL ROUTINES - SINGLE/VBLA=SNRM2                             ZSPT0140
C                       - DOUBLE/VBLA=DNRM2                             ZSPT0150
C                                                                       ZSPT0160
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           ZSPT0170
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      ZSPT0180
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  ZSPT0190
C                                                                       ZSPT0200
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.       ZSPT0210
C                                                                       ZSPT0220
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN ZSPT0230
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    ZSPT0240
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        ZSPT0250
C                                                                       ZSPT0260
C-----------------------------------------------------------------------ZSPT0270
C                                                                       ZSPT0280
      SUBROUTINE ZSPWC (N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)                ZSPT0290
C                                  SPECIFICATIONS FOR ARGUMENTS         ZSPT0300
      INTEGER            N,LR                                           ZSPT0310
      REAL               R(LR),DIAG(N),QTB(N),DELTA,X(N),WA1(N),WA2(N)  ZSPT0320
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   ZSPT0330
      INTEGER            I,JJ,JP1,J,K,L                                 ZSPT0340
      REAL               ALPHA,BNORM,EPSMCH,GNORM,ONE,QNORM,SGNORM,     ZSPT0350
     *                   SPMPAR,SUM,TEMP,ZERO                           ZSPT0360
      REAL               SNRM2                                          ZSPT0370
      DATA               SPMPAR /1.1921E-07/                            ZSPT0380
c      DATA				 SPMPAR /2.2204E-16/
      DATA               ONE,ZERO /1.0E0,0.0E0/                         ZSPT0390
C                                  EPSMCH IS THE MACHINE PRECISION.     ZSPT0400
C                                  FIRST EXECUTABLE STATEMENT           ZSPT0410
      EPSMCH = SPMPAR                                                   ZSPT0420
C                                  FIRST, CALCULATE THE GAUSS-NEWTON    ZSPT0430
C                                  DIRECTION.                           ZSPT0440
      JJ = (N*(N+1))/2+1                                                ZSPT0450
      DO 25 K=1,N                                                       ZSPT0460
         J = N-K+1                                                      ZSPT0470
         JP1 = J+1                                                      ZSPT0480
         JJ = JJ-K                                                      ZSPT0490
         L = JJ+1                                                       ZSPT0500
         SUM = ZERO                                                     ZSPT0510
         IF (N.LT.JP1) GO TO 10                                         ZSPT0520
         DO 5 I=JP1,N                                                   ZSPT0530
            SUM = SUM+R(L)*X(I)                                         ZSPT0540
            L = L+1                                                     ZSPT0550
    5    CONTINUE                                                       ZSPT0560
   10    CONTINUE                                                       ZSPT0570
         TEMP = R(JJ)                                                   ZSPT0580
         IF (TEMP.NE.ZERO) GO TO 20                                     ZSPT0590
         L = J                                                          ZSPT0600
         DO 15 I=1,J                                                    ZSPT0610
            TEMP = AMAX1(TEMP,ABS(R(L)))                                ZSPT0620
            L = L+N-I                                                   ZSPT0630
   15    CONTINUE                                                       ZSPT0640
         TEMP = EPSMCH*TEMP                                             ZSPT0650
         IF (TEMP.EQ.ZERO) TEMP = EPSMCH                                ZSPT0660
   20    CONTINUE                                                       ZSPT0670
         X(J) = (QTB(J)-SUM)/TEMP                                       ZSPT0680
   25 CONTINUE                                                          ZSPT0690
C                                  TEST WHETHER THE GAUSS-NEWTON        ZSPT0700
C                                  DIRECTION IS ACCEPTABLE.             ZSPT0710
      DO 30 J=1,N                                                       ZSPT0720
         WA1(J) = ZERO                                                  ZSPT0730
         WA2(J) = DIAG(J)*X(J)                                          ZSPT0740
   30 CONTINUE                                                          ZSPT0750
      QNORM = SNRM2(N,WA2,1)                                            ZSPT0760
      IF (QNORM.LE.DELTA) GO TO 70                                      ZSPT0770
C                                  THE GAUSS-NEWTON DIRECTION IS NOT    ZSPT0780
C                                  ACCEPTABLE. NEXT, CALCULATE THE      ZSPT0790
C                                  SCALED GRADIENT DIRECTION.           ZSPT0800
      L = 1                                                             ZSPT0810
      DO 40 J=1,N                                                       ZSPT0820
         TEMP = QTB(J)                                                  ZSPT0830
         DO 35 I=J,N                                                    ZSPT0840
            WA1(I) = WA1(I)+R(L)*TEMP                                   ZSPT0850
            L = L+1                                                     ZSPT0860
   35    CONTINUE                                                       ZSPT0870
         WA1(J) = WA1(J)/DIAG(J)                                        ZSPT0880
   40 CONTINUE                                                          ZSPT0890
C                                  CALCULATE THE NORM OF THE SCALED     ZSPT0900
C                                  GRADIENT AND TEST FOR THE SPECIAL    ZSPT0910
C                                  CASE IN WHICH THE SCALED GRADIENT IS ZSPT0920
C                                  ZERO.                                ZSPT0930
      GNORM = SNRM2(N,WA1,1)                                            ZSPT0940
      SGNORM = ZERO                                                     ZSPT0950
      ALPHA = DELTA/QNORM                                               ZSPT0960
      IF (GNORM.EQ.ZERO) GO TO 60                                       ZSPT0970
C                                  CALCULATE THE POINT ALONG THE SCALED ZSPT0980
C                                  GRADIENT AT WHICH THE QUADRATIC IS   ZSPT0990
C                                  MINIMIZED.                           ZSPT1000
      DO 45 J=1,N                                                       ZSPT1010
         WA1(J) = (WA1(J)/GNORM)/DIAG(J)                                ZSPT1020
   45 CONTINUE                                                          ZSPT1030
      L = 1                                                             ZSPT1040
      DO 55 J=1,N                                                       ZSPT1050
         SUM = ZERO                                                     ZSPT1060
         DO 50 I=J,N                                                    ZSPT1070
            SUM = SUM+R(L)*WA1(I)                                       ZSPT1080
            L = L+1                                                     ZSPT1090
   50    CONTINUE                                                       ZSPT1100
         WA2(J) = SUM                                                   ZSPT1110
   55 CONTINUE                                                          ZSPT1120
      TEMP = SNRM2(N,WA2,1)                                             ZSPT1130
      SGNORM = (GNORM/TEMP)/TEMP                                        ZSPT1140
C                                  TEST WHETHER THE SCALED GRADIENT     ZSPT1150
C                                  DIRECTION IS ACCEPTABLE.             ZSPT1160
      ALPHA = ZERO                                                      ZSPT1170
      IF (SGNORM.GE.DELTA) GO TO 60                                     ZSPT1180
C                                  THE SCALED GRADIENT DIRECTION IS NOT ZSPT1190
C                                  ACCEPTABLE. FINALLY, CALCULATE THE   ZSPT1200
C                                  POINT ALONG THE DOGLEG AT WHICH THE  ZSPT1210
C                                  QUADRATIC IS MINIMIZED.              ZSPT1220
      BNORM = SNRM2(N,QTB,1)                                            ZSPT1230
      TEMP = (BNORM/GNORM)*(BNORM/QNORM)*(SGNORM/DELTA)                 ZSPT1240
      TEMP = TEMP-(DELTA/QNORM)*(SGNORM/DELTA)**2+SQRT((TEMP-(DELTA     ZSPT1250
     */QNORM))**2+(ONE-(DELTA/QNORM)**2)*(ONE-(SGNORM/DELTA)**2))       ZSPT1260
      ALPHA = ((DELTA/QNORM)*(ONE-(SGNORM/DELTA)**2))/TEMP              ZSPT1270
   60 CONTINUE                                                          ZSPT1280
C                                  FORM APPROPRIATE CONVEX COMBINATION  ZSPT1290
C                                  OF THE GAUSS-NEWTON DIRECTION AND THEZSPT1300
C                                  SCALED GRADIENT DIRECTION.           ZSPT1310
      TEMP = (ONE-ALPHA)*AMIN1(SGNORM,DELTA)                            ZSPT1320
      DO 65 J=1,N                                                       ZSPT1330
         X(J) = TEMP*WA1(J)+ALPHA*X(J)                                  ZSPT1340
   65 CONTINUE                                                          ZSPT1350
   70 CONTINUE                                                          ZSPT1360
      RETURN                                                            ZSPT1370
      END                                                               ZSPT1380
