C   IMSL ROUTINE NAME   - VBLA=SNRM2                                    VBSH0010
C                                                                       VBSH0020
C-----------------------------------------------------------------------VBSH0030
C                                                                       VBSH0040
C   COMPUTER            - VAX/SINGLE                                    VBSH0050
C                                                                       VBSH0060
C   LATEST REVISION     - JANUARY 1, 1978                               VBSH0070
C                                                                       VBSH0080
C   PURPOSE             - COMPUTE THE EUCLIDEAN LENGTH OR L2 NORM       VBSH0090
C                           OF A SINGLE PRECISION VECTOR                VBSH0100
C                                                                       VBSH0110
C   USAGE               - FUNCTION SNRM2 (N,SX,INCX)                    VBSH0120
C                                                                       VBSH0130
C   ARGUMENTS    SNRM2  - SQUARE ROOT OF THE SUM FROM I=1 TO N OF       VBSH0140
C                           X(I)**2. (OUTPUT)                           VBSH0150
C                           X(I) REFERS TO A SPECIFIC ELEMENT OF SX.    VBSH0160
C                           SEE INCX ARGUMENT DESCRIPTION.              VBSH0170
C                N      - LENGTH OF VECTOR X. (INPUT)                   VBSH0180
C                SX     - REAL VECTOR OF LENGTH N*INCX. (INPUT)         VBSH0190
C                INCX   - DISPLACEMENT BETWEEN ELEMENTS OF SX. (INPUT)  VBSH0200
C                           X(I) IS DEFINED TO BE SX(1+(I-1)*INCX).     VBSH0210
C                           INCX MUST BE GREATER THAN ZERO.             VBSH0220
C                                                                       VBSH0230
C   PRECISION/HARDWARE  - SINGLE/ALL                                    VBSH0240
C                                                                       VBSH0250
C   REQD. IMSL ROUTINES - NONE REQUIRED                                 VBSH0260
C                                                                       VBSH0270
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           VBSH0280
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      VBSH0290
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  VBSH0300
C                                                                       VBSH0310
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       VBSH0320
C                                                                       VBSH0330
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN VBSH0340
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    VBSH0350
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        VBSH0360
C                                                                       VBSH0370
C-----------------------------------------------------------------------VBSH0380
C                                                                       VBSH0390
      REAL FUNCTION SNRM2 (N,SX,INCX)                                   VBSH0400
C                                  SPECIFICATIONS FOR ARGUMENTS         VBSH0410
      INTEGER            N,INCX                                         VBSH0420
      REAL               SX(1)                                          VBSH0430
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   VBSH0440
      INTEGER            I,J,NEXT,NN                                    VBSH0450
      REAL               CUTLO,CUTHI,HITEST,SUM,XMAX,ZERO,ONE           VBSH0460
      DATA               ZERO, ONE /0.0E0, 1.0E0/                       VBSH0470
      DATA               CUTLO, CUTHI / 4.441E-16,  1.304E19/           VBSH0480
C                                  FIRST EXECUTABLE STATEMENT           VBSH0490
C                                                                       VBSH0500
      IF (N.GT.0) GO TO 5                                               VBSH0510
      SNRM2 = ZERO                                                      VBSH0520
      GO TO 70                                                          VBSH0530
C                                                                       VBSH0540
    5 ASSIGN 15 TO NEXT                                                 VBSH0550
      SUM = ZERO                                                        VBSH0560
      NN = N*INCX                                                       VBSH0570
C                                  BEGIN MAIN LOOP                      VBSH0580
      I = 1                                                             VBSH0590
   10 GO TO NEXT, (15,20,35,40)                                         VBSH0600
   15 IF (ABS(SX(I)).GT.CUTLO) GO TO 55                                 VBSH0610
      ASSIGN 20 TO NEXT                                                 VBSH0620
      XMAX = ZERO                                                       VBSH0630
C                                  PHASE 1. SUM IS ZERO                 VBSH0640
   20 IF (SX(I).EQ.ZERO) GO TO 65                                       VBSH0650
      IF (ABS(SX(I)).GT.CUTLO) GO TO 55                                 VBSH0660
C                                  PREPARE FOR PHASE 2.                 VBSH0670
      ASSIGN 35 TO NEXT                                                 VBSH0680
      GO TO 30                                                          VBSH0690
C                                  PREPARE FOR PHASE 4.                 VBSH0700
   25 I = J                                                             VBSH0710
      ASSIGN 40 TO NEXT                                                 VBSH0720
      SUM = (SUM/SX(I))/SX(I)                                           VBSH0730
   30 XMAX = ABS(SX(I))                                                 VBSH0740
      GO TO 45                                                          VBSH0750
C                                  PHASE 2. SUM IS SMALL. SCALE TO      VBSH0760
C                                    AVOID DESTRUCTIVE UNDERFLOW.       VBSH0770
   35 IF (ABS(SX(I)).GT.CUTLO) GO TO 50                                 VBSH0780
C                                  COMMON CODE FOR PHASES 2 AND 4. IN   VBSH0790
C                                    PHASE 4 SUM IS LARGE. SCALE TO     VBSH0800
C                                    AVOID OVERFLOW.                    VBSH0810
   40 IF (ABS(SX(I)).LE.XMAX) GO TO 45                                  VBSH0820
      SUM = ONE+SUM*(XMAX/SX(I))**2                                     VBSH0830
      XMAX = ABS(SX(I))                                                 VBSH0840
      GO TO 65                                                          VBSH0850
C                                                                       VBSH0860
   45 SUM = SUM+(SX(I)/XMAX)**2                                         VBSH0870
      GO TO 65                                                          VBSH0880
C                                  PREPARE FOR PHASE 3.                 VBSH0890
   50 SUM = (SUM*XMAX)*XMAX                                             VBSH0900
C                                  FOR REAL OR D.P. SET HITEST =        VBSH0910
C                                    CUTHI/N FOR COMPLEX SET HITEST =   VBSH0920
C                                    CUTHI/(2*N)                        VBSH0930
   55 HITEST = CUTHI/FLOAT(N)                                           VBSH0940
C                                  PHASE 3. SUM IS MID-RANGE. NO        VBSH0950
C                                    SCALING.                           VBSH0960
      DO 60 J=I,NN,INCX                                                 VBSH0970
         IF (ABS(SX(J)).GE.HITEST) GO TO 25                             VBSH0980
   60 SUM = SUM+SX(J)**2                                                VBSH0990
      SNRM2 = SQRT(SUM)                                                 VBSH1000
      GO TO 70                                                          VBSH1010
C                                                                       VBSH1020
   65 CONTINUE                                                          VBSH1030
      I = I+INCX                                                        VBSH1040
      IF (I.LE.NN) GO TO 10                                             VBSH1050
C                                  END OF MAIN LOOP. COMPUTE SQUARE     VBSH1060
C                                    ROOT AND ADJUST FOR SCALING.       VBSH1070
      SNRM2 = XMAX*SQRT(SUM)                                            VBSH1080
   70 CONTINUE                                                          VBSH1090
      RETURN                                                            VBSH1100
      END                                                               VBSH1110
