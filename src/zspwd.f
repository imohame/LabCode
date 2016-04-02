C   IMSL ROUTINE NAME   - ZSPWD                                         ZSPU0010
C                                                                       ZSPU0020
C-----------------------------------------------------------------------ZSPU0030
C                                                                       ZSPU0040
C   COMPUTER            - VAX/SINGLE                                    ZSPU0050
C                                                                       ZSPU0060
C   LATEST REVISION     - JUNE 1, 1982                                  ZSPU0070
C                                                                       ZSPU0080
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSPOW  ZSPU0090
C                                                                       ZSPU0100
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         ZSPU0110
C                       - SINGLE/H36,H48,H60                            ZSPU0120
C                                                                       ZSPU0130
C   REQD. IMSL ROUTINES - NONE REQUIRED                                 ZSPU0140
C                                                                       ZSPU0150
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           ZSPU0160
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      ZSPU0170
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  ZSPU0180
C                                                                       ZSPU0190
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.       ZSPU0200
C                                                                       ZSPU0210
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN ZSPU0220
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    ZSPU0230
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        ZSPU0240
C                                                                       ZSPU0250
C-----------------------------------------------------------------------ZSPU0260
C                                                                       ZSPU0270
      SUBROUTINE ZSPWD (M,N,A,LDA,V,W)                                  ZSPU0280
C                                  SPECIFICATIONS FOR ARGUMENTS         ZSPU0290
      INTEGER            M,N,LDA                                        ZSPU0300
      REAL               A(LDA,N),V(N),W(N)                             ZSPU0310
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   ZSPU0320
      INTEGER            I,J,NM1,NMJ                                    ZSPU0330
      REAL               TEMP1,ONE,TEMP2,TEMP                           ZSPU0340
      DATA               ONE /1.0E0/                                    ZSPU0350
C                                  APPLY THE FIRST SET OF GIVENS        ZSPU0360
C                                  ROTATIONS TO A.                      ZSPU0370
C                                  FIRST EXECUTABLE STATEMENT           ZSPU0380
      NM1 = N-1                                                         ZSPU0390
      IF (NM1.LT.1) GO TO 25                                            ZSPU0400
      DO 10 NMJ=1,NM1                                                   ZSPU0410
         J = N-NMJ                                                      ZSPU0420
         IF (ABS(V(J)).GT.ONE) TEMP1 = ONE/V(J)                         ZSPU0430
         IF (ABS(V(J)).GT.ONE) TEMP2 = SQRT(ONE-TEMP1**2)               ZSPU0440
         IF (ABS(V(J)).LE.ONE) TEMP2 = V(J)                             ZSPU0450
         IF (ABS(V(J)).LE.ONE) TEMP1 = SQRT(ONE-TEMP2**2)               ZSPU0460
         DO 5 I=1,M                                                     ZSPU0470
            TEMP = TEMP1*A(I,J)-TEMP2*A(I,N)                            ZSPU0480
            A(I,N) = TEMP2*A(I,J)+TEMP1*A(I,N)                          ZSPU0490
            A(I,J) = TEMP                                               ZSPU0500
    5    CONTINUE                                                       ZSPU0510
   10 CONTINUE                                                          ZSPU0520
C                                  APPLY THE SECOND SET OF GIVENS       ZSPU0530
C                                  ROTATIONS TO A.                      ZSPU0540
      DO 20 J=1,NM1                                                     ZSPU0550
         IF (ABS(W(J)).GT.ONE) TEMP1 = ONE/W(J)                         ZSPU0560
         IF (ABS(W(J)).GT.ONE) TEMP2 = SQRT(ONE-TEMP1**2)               ZSPU0570
         IF (ABS(W(J)).LE.ONE) TEMP2 = W(J)                             ZSPU0580
         IF (ABS(W(J)).LE.ONE) TEMP1 = SQRT(ONE-TEMP2**2)               ZSPU0590
         DO 15 I=1,M                                                    ZSPU0600
            TEMP = TEMP1*A(I,J)+TEMP2*A(I,N)                            ZSPU0610
            A(I,N) = -TEMP2*A(I,J)+TEMP1*A(I,N)                         ZSPU0620
            A(I,J) = TEMP                                               ZSPU0630
   15    CONTINUE                                                       ZSPU0640
   20 CONTINUE                                                          ZSPU0650
   25 CONTINUE                                                          ZSPU0660
      RETURN                                                            ZSPU0670
      END                                                               ZSPU0680
