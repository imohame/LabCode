C   IMSL ROUTINE NAME   - ZSPWB                                         ZSPS0010
C                                                                       ZSPS0020
C-----------------------------------------------------------------------ZSPS0030
C                                                                       ZSPS0040
C   COMPUTER            - VAX/SINGLE                                    ZSPS0050
C                                                                       ZSPS0060
C   LATEST REVISION     - JUNE 1, 1982                                  ZSPS0070
C                                                                       ZSPS0080
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSPOW  ZSPS0090
C                                                                       ZSPS0100
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         ZSPS0110
C                       - SINGLE/H36,H48,H60                            ZSPS0120
C                                                                       ZSPS0130
C   REQD. IMSL ROUTINES - NONE REQUIRED                                 ZSPS0140
C                                                                       ZSPS0150
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           ZSPS0160
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      ZSPS0170
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  ZSPS0180
C                                                                       ZSPS0190
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.       ZSPS0200
C                                                                       ZSPS0210
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN ZSPS0220
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    ZSPS0230
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        ZSPS0240
C                                                                       ZSPS0250
C-----------------------------------------------------------------------ZSPS0260
C                                                                       ZSPS0270
      SUBROUTINE ZSPWB (FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,WA1,ZSPS0280
     *                   WA2,PAR)                                       ZSPS0290
C                                  SPECIFICATIONS FOR ARGUMENTS         ZSPS0300
      INTEGER            N,LDFJAC,IFLAG,ML,MU                           ZSPS0310
      REAL               X(N),FVEC(N),FJAC(LDFJAC,N),EPSFCN,WA1(N),     ZSPS0320
     *                   WA2(N),PAR(1)                                  ZSPS0330
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   ZSPS0340
      INTEGER            I,J,K,MSUM                                     ZSPS0350
      REAL               EPSMCH,EPS,H,SPMPAR,TEMP,ZERO                  ZSPS0360
      DATA               SPMPAR /1.1921E-07/                            ZSPS0370
c      DATA				 SPMPAR /2.2204E-16/
      DATA               ZERO /0.0E0/                                   ZSPS0380
C                                  EPSMCH IS THE MACHINE PRECISION.     ZSPS0390
C                                  FIRST EXECUTABLE STATEMENT           ZSPS0400
      EPSMCH = SPMPAR                                                   ZSPS0410
      EPS = SQRT(AMAX1(EPSFCN,EPSMCH))                                  ZSPS0420
      MSUM = ML+MU+1                                                    ZSPS0430
      IF (MSUM.LT.N) GO TO 20                                           ZSPS0440
C                                  COMPUTATION OF DENSE APPROXIMATE     ZSPS0450
C                                  JACOBIAN.                            ZSPS0460
      DO 10 J=1,N                                                       ZSPS0470
         TEMP = X(J)                                                    ZSPS0480
         H = EPS*ABS(TEMP)                                              ZSPS0490
         IF (H.EQ.ZERO) H = EPS                                         ZSPS0500
         X(J) = TEMP+H                                                  ZSPS0510
         CALL FCN(X,WA1,N,PAR)                                          ZSPS0520
         IF (IFLAG.LT.0) GO TO 15                                       ZSPS0530
         X(J) = TEMP                                                    ZSPS0540
         DO 5 I=1,N                                                     ZSPS0550
            FJAC(I,J) = (WA1(I)-FVEC(I))/H                              ZSPS0560
    5    CONTINUE                                                       ZSPS0570
   10 CONTINUE                                                          ZSPS0580
   15 CONTINUE                                                          ZSPS0590
      GO TO 50                                                          ZSPS0600
   20 CONTINUE                                                          ZSPS0610
C                                  COMPUTATION OF BANDED APPROXIMATE    ZSPS0620
C                                  JACOBIAN.                            ZSPS0630
      DO 40 K=1,MSUM                                                    ZSPS0640
         DO 25 J=K,N,MSUM                                               ZSPS0650
            WA2(J) = X(J)                                               ZSPS0660
            H = EPS*ABS(WA2(J))                                         ZSPS0670
            IF (H.EQ.ZERO) H = EPS                                      ZSPS0680
            X(J) = WA2(J)+H                                             ZSPS0690
   25    CONTINUE                                                       ZSPS0700
         CALL FCN(X,WA1,N,PAR)                                          ZSPS0710
         IF (IFLAG.LT.0) GO TO 45                                       ZSPS0720
         DO 35 J=K,N,MSUM                                               ZSPS0730
            X(J) = WA2(J)                                               ZSPS0740
            H = EPS*ABS(WA2(J))                                         ZSPS0750
            IF (H.EQ.ZERO) H = EPS                                      ZSPS0760
            DO 30 I=1,N                                                 ZSPS0770
               FJAC(I,J) = ZERO                                         ZSPS0780
               IF (I.GE.J-MU .AND. I.LE.J+ML) FJAC(I,J) =               ZSPS0790
     *         (WA1(I)-FVEC(I))/H                                       ZSPS0800
   30       CONTINUE                                                    ZSPS0810
   35    CONTINUE                                                       ZSPS0820
   40 CONTINUE                                                          ZSPS0830
   45 CONTINUE                                                          ZSPS0840
   50 CONTINUE                                                          ZSPS0850
      RETURN                                                            ZSPS0860
      END                                                               ZSPS0870
