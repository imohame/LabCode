C   IMSL ROUTINE NAME   - ZSPWE                                         ZSPV0010
C                                                                       ZSPV0020
C-----------------------------------------------------------------------ZSPV0030
C                                                                       ZSPV0040
C   COMPUTER            - VAX/SINGLE                                    ZSPV0050
C                                                                       ZSPV0060
C   LATEST REVISION     - JUNE 1, 1982                                  ZSPV0070
C                                                                       ZSPV0080
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSPOW  ZSPV0090
C                                                                       ZSPV0100
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         ZSPV0110
C                       - SINGLE/H36,H48,H60                            ZSPV0120
C                                                                       ZSPV0130
C   REQD. IMSL ROUTINES - NONE REQUIRED                                 ZSPV0140
C                                                                       ZSPV0150
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           ZSPV0160
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      ZSPV0170
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  ZSPV0180
C                                                                       ZSPV0190
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.       ZSPV0200
C                                                                       ZSPV0210
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN ZSPV0220
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    ZSPV0230
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        ZSPV0240
C                                                                       ZSPV0250
C-----------------------------------------------------------------------ZSPV0260
C                                                                       ZSPV0270
      SUBROUTINE ZSPWE (M,N,S,LS,U,V,W,SING)                            ZSPV0280
C                                  SPECIFICATIONS FOR ARGUMENTS         ZSPV0290
      INTEGER            M,N,LS                                         ZSPV0300
      REAL               S(LS),U(M),V(N),W(M)                           ZSPV0310
      LOGICAL            SING                                           ZSPV0320
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   ZSPV0330
      INTEGER            I,JJ,J,L,NM1,NMJ                               ZSPV0340
      REAL               TEMP1,TEMP2,GIANT,ONE,P25,P5,TEMP3,SPMPAR,     ZSPV0350
     *                   TEMP4,TAU,TEMP,ZERO                            ZSPV0360
      DATA               GIANT /1.7E+38/                                ZSPV0370
      DATA               ONE,P5,P25,ZERO /1.0E0,5.0E-1,2.5E-1,0.0E0/    ZSPV0380
C                                  INITIALIZE THE DIAGONAL ELEMENT      ZSPV0390
C                                  POINTER.                             ZSPV0400
C                                  FIRST EXECUTABLE STATEMENT           ZSPV0410
      JJ = (N*(2*M-N+1))/2-(M-N)                                        ZSPV0420
C                                  MOVE THE NONTRIVIAL PART OF THE LAST ZSPV0430
C                                  COLUMN OF S INTO W.                  ZSPV0440
      L = JJ                                                            ZSPV0450
      DO 5 I=N,M                                                        ZSPV0460
         W(I) = S(L)                                                    ZSPV0470
         L = L+1                                                        ZSPV0480
    5 CONTINUE                                                          ZSPV0490
C                                  ROTATE THE VECTOR V INTO A MULTIPLE  ZSPV0500
C                                  OF THE N-TH UNIT VECTOR IN SUCH A WAYZSPV0510
C                                  THAT A SPIKE IS INTRODUCED INTO W.   ZSPV0520
      NM1 = N-1                                                         ZSPV0530
      IF (NM1.LT.1) GO TO 35                                            ZSPV0540
      DO 30 NMJ=1,NM1                                                   ZSPV0550
         J = N-NMJ                                                      ZSPV0560
         JJ = JJ-(M-J+1)                                                ZSPV0570
         W(J) = ZERO                                                    ZSPV0580
         IF (V(J).EQ.ZERO) GO TO 25                                     ZSPV0590
C                                  DETERMINE A GIVENS ROTATION WHICH    ZSPV0600
C                                  ELIMINATES THE J-TH ELEMENT OF V.    ZSPV0610
         IF (ABS(V(N)).GE.ABS(V(J))) GO TO 10                           ZSPV0620
         TEMP2 = V(N)/V(J)                                              ZSPV0630
         TEMP3 = P5/SQRT(P25+P25*TEMP2**2)                              ZSPV0640
         TEMP1 = TEMP3*TEMP2                                            ZSPV0650
         TAU = ONE                                                      ZSPV0660
         IF (ABS(TEMP1)*GIANT.GT.ONE) TAU = ONE/TEMP1                   ZSPV0670
         GO TO 15                                                       ZSPV0680
   10    CONTINUE                                                       ZSPV0690
         TEMP4 = V(J)/V(N)                                              ZSPV0700
         TEMP1 = P5/SQRT(P25+P25*TEMP4**2)                              ZSPV0710
         TEMP3 = TEMP1*TEMP4                                            ZSPV0720
         TAU = TEMP3                                                    ZSPV0730
   15    CONTINUE                                                       ZSPV0740
C                                  APPLY THE TRANSFORMATION TO V AND    ZSPV0750
C                                  STORE THE INFORMATION NECESSARY TO   ZSPV0760
C                                  RECOVER THE GIVENS ROTATION.         ZSPV0770
         V(N) = TEMP3*V(J)+TEMP1*V(N)                                   ZSPV0780
         V(J) = TAU                                                     ZSPV0790
C                                  APPLY THE TRANSFORMATION TO S AND    ZSPV0800
C                                  EXTEND THE SPIKE IN W.               ZSPV0810
         L = JJ                                                         ZSPV0820
         DO 20 I=J,M                                                    ZSPV0830
            TEMP = TEMP1*S(L)-TEMP3*W(I)                                ZSPV0840
            W(I) = TEMP3*S(L)+TEMP1*W(I)                                ZSPV0850
            S(L) = TEMP                                                 ZSPV0860
            L = L+1                                                     ZSPV0870
   20    CONTINUE                                                       ZSPV0880
   25    CONTINUE                                                       ZSPV0890
   30 CONTINUE                                                          ZSPV0900
   35 CONTINUE                                                          ZSPV0910
C                                  ADD THE SPIKE FROM THE RANK 1 UPDATE ZSPV0920
C                                  TO W.                                ZSPV0930
      DO 40 I=1,M                                                       ZSPV0940
         W(I) = W(I)+V(N)*U(I)                                          ZSPV0950
   40 CONTINUE                                                          ZSPV0960
C                                  ELIMINATE THE SPIKE.                 ZSPV0970
      SING = .FALSE.                                                    ZSPV0980
      IF (NM1.LT.1) GO TO 70                                            ZSPV0990
      DO 65 J=1,NM1                                                     ZSPV1000
         IF (W(J).EQ.ZERO) GO TO 60                                     ZSPV1010
C                                  DETERMINE A GIVENS ROTATION WHICH    ZSPV1020
C                                  ELIMINATES THE J-TH ELEMENT OF THE   ZSPV1030
C                                  SPIKE.                               ZSPV1040
         IF (ABS(S(JJ)).GE.ABS(W(J))) GO TO 45                          ZSPV1050
         TEMP2 = S(JJ)/W(J)                                             ZSPV1060
         TEMP3 = P5/SQRT(P25+P25*TEMP2**2)                              ZSPV1070
         TEMP1 = TEMP3*TEMP2                                            ZSPV1080
         TAU = ONE                                                      ZSPV1090
         IF (ABS(TEMP1)*GIANT.GT.ONE) TAU = ONE/TEMP1                   ZSPV1100
         GO TO 50                                                       ZSPV1110
   45    CONTINUE                                                       ZSPV1120
         TEMP4 = W(J)/S(JJ)                                             ZSPV1130
         TEMP1 = P5/SQRT(P25+P25*TEMP4**2)                              ZSPV1140
         TEMP3 = TEMP1*TEMP4                                            ZSPV1150
         TAU = TEMP3                                                    ZSPV1160
   50    CONTINUE                                                       ZSPV1170
C                                  APPLY THE TRANSFORMATION TO S AND    ZSPV1180
C                                  REDUCE THE SPIKE IN W.               ZSPV1190
         L = JJ                                                         ZSPV1200
         DO 55 I=J,M                                                    ZSPV1210
            TEMP = TEMP1*S(L)+TEMP3*W(I)                                ZSPV1220
            W(I) = -TEMP3*S(L)+TEMP1*W(I)                               ZSPV1230
            S(L) = TEMP                                                 ZSPV1240
            L = L+1                                                     ZSPV1250
   55    CONTINUE                                                       ZSPV1260
C                                  STORE THE INFORMATION NECESSARY TO   ZSPV1270
C                                  RECOVER THE GIVENS ROTATION.         ZSPV1280
         W(J) = TAU                                                     ZSPV1290
   60    CONTINUE                                                       ZSPV1300
C                                  TEST FOR ZERO DIAGONAL ELEMENTS IN   ZSPV1310
C                                  THE OUTPUT S.                        ZSPV1320
         IF (S(JJ).EQ.ZERO) SING = .TRUE.                               ZSPV1330
         JJ = JJ+(M-J+1)                                                ZSPV1340
   65 CONTINUE                                                          ZSPV1350
   70 CONTINUE                                                          ZSPV1360
C                                  MOVE W BACK INTO THE LAST COLUMN OF  ZSPV1370
C                                  THE OUTPUT S.                        ZSPV1380
      L = JJ                                                            ZSPV1390
      DO 75 I=N,M                                                       ZSPV1400
         S(L) = W(I)                                                    ZSPV1410
         L = L+1                                                        ZSPV1420
   75 CONTINUE                                                          ZSPV1430
      IF (S(JJ).EQ.ZERO) SING = .TRUE.                                  ZSPV1440
      RETURN                                                            ZSPV1450
      END                                                               ZSPV1460
