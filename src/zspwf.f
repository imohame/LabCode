C   IMSL ROUTINE NAME   - ZSPWF                                         ZSPW0010
C                                                                       ZSPW0020
C-----------------------------------------------------------------------ZSPW0030
C                                                                       ZSPW0040
C   COMPUTER            - VAX/SINGLE                                    ZSPW0050
C                                                                       ZSPW0060
C   LATEST REVISION     - JUNE 1, 1982                                  ZSPW0070
C                                                                       ZSPW0080
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSPOW  ZSPW0090
C                                                                       ZSPW0100
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         ZSPW0110
C                       - SINGLE/H36,H48,H60                            ZSPW0120
C                                                                       ZSPW0130
C   REQD. IMSL ROUTINES - NONE REQUIRED                                 ZSPW0140
C                                                                       ZSPW0150
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           ZSPW0160
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      ZSPW0170
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  ZSPW0180
C                                                                       ZSPW0190
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.       ZSPW0200
C                                                                       ZSPW0210
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN ZSPW0220
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    ZSPW0230
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        ZSPW0240
C                                                                       ZSPW0250
C-----------------------------------------------------------------------ZSPW0260
C                                                                       ZSPW0270
      SUBROUTINE ZSPWF (M,N,Q,LDQ,WA)                                   ZSPW0280
C                                  SPECIFICATIONS FOR ARGUMENTS         ZSPW0290
      INTEGER            M,N,LDQ                                        ZSPW0300
      REAL               Q(LDQ,M),WA(M)                                 ZSPW0310
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   ZSPW0320
      INTEGER            I,JM1,J,K,L,MINMN,NP1                          ZSPW0330
      REAL               ONE,SUM,TEMP,ZERO                              ZSPW0340
      DATA               ONE,ZERO /1.0E0,0.0E0/                         ZSPW0350
C                                  ZERO OUT UPPER TRIANGLE OF Q IN THE  ZSPW0360
C                                  FIRST MIN(M,N) COLUMNS.              ZSPW0370
C                                  FIRST EXECUTABLE STATEMENT           ZSPW0380
      MINMN = MIN0(M,N)                                                 ZSPW0390
      IF (MINMN.LT.2) GO TO 15                                          ZSPW0400
      DO 10 J=2,MINMN                                                   ZSPW0410
         JM1 = J-1                                                      ZSPW0420
         DO 5 I=1,JM1                                                   ZSPW0430
            Q(I,J) = ZERO                                               ZSPW0440
    5    CONTINUE                                                       ZSPW0450
   10 CONTINUE                                                          ZSPW0460
   15 CONTINUE                                                          ZSPW0470
C                                  INITIALIZE REMAINING COLUMNS TO THOSEZSPW0480
C                                  OF THE IDENTITY MATRIX.              ZSPW0490
      NP1 = N+1                                                         ZSPW0500
      IF (M.LT.NP1) GO TO 30                                            ZSPW0510
      DO 25 J=NP1,M                                                     ZSPW0520
         DO 20 I=1,M                                                    ZSPW0530
            Q(I,J) = ZERO                                               ZSPW0540
   20    CONTINUE                                                       ZSPW0550
         Q(J,J) = ONE                                                   ZSPW0560
   25 CONTINUE                                                          ZSPW0570
   30 CONTINUE                                                          ZSPW0580
C                                  ACCUMULATE Q FROM ITS FACTORED FORM. ZSPW0590
      DO 60 L=1,MINMN                                                   ZSPW0600
         K = MINMN-L+1                                                  ZSPW0610
         DO 35 I=K,M                                                    ZSPW0620
            WA(I) = Q(I,K)                                              ZSPW0630
            Q(I,K) = ZERO                                               ZSPW0640
   35    CONTINUE                                                       ZSPW0650
         Q(K,K) = ONE                                                   ZSPW0660
         IF (WA(K).EQ.ZERO) GO TO 55                                    ZSPW0670
         DO 50 J=K,M                                                    ZSPW0680
            SUM = ZERO                                                  ZSPW0690
            DO 40 I=K,M                                                 ZSPW0700
               SUM = SUM+Q(I,J)*WA(I)                                   ZSPW0710
   40       CONTINUE                                                    ZSPW0720
            TEMP = SUM/WA(K)                                            ZSPW0730
            DO 45 I=K,M                                                 ZSPW0740
               Q(I,J) = Q(I,J)-TEMP*WA(I)                               ZSPW0750
   45       CONTINUE                                                    ZSPW0760
   50    CONTINUE                                                       ZSPW0770
   55    CONTINUE                                                       ZSPW0780
   60 CONTINUE                                                          ZSPW0790
      RETURN                                                            ZSPW0800
      END                                                               ZSPW0810
