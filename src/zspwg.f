C   IMSL ROUTINE NAME   - ZSPWG                                         ZSPX0010
C                                                                       ZSPX0020
C-----------------------------------------------------------------------ZSPX0030
C                                                                       ZSPX0040
C   COMPUTER            - VAX/SINGLE                                    ZSPX0050
C                                                                       ZSPX0060
C   LATEST REVISION     - JUNE 1, 1982                                  ZSPX0070
C                                                                       ZSPX0080
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSPOW  ZSPX0090
C                                                                       ZSPX0100
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         ZSPX0110
C                       - SINGLE/H36,H48,H60                            ZSPX0120
C                                                                       ZSPX0130
C   REQD. IMSL ROUTINES - SINGLE/VBLA=SNRM2                             ZSPX0140
C                       - DOUBLE/VBLA=DNRM2                             ZSPX0150
C                                                                       ZSPX0160
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           ZSPX0170
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      ZSPX0180
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  ZSPX0190
C                                                                       ZSPX0200
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.       ZSPX0210
C                                                                       ZSPX0220
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN ZSPX0230
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    ZSPX0240
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        ZSPX0250
C                                                                       ZSPX0260
C-----------------------------------------------------------------------ZSPX0270
C                                                                       ZSPX0280
      SUBROUTINE ZSPWG (M,N,A,LDA,PIVOT,IPVT,LIPVT,RDIAG,ACNORM,WA)     ZSPX0290
C                                  SPECIFICATIONS FOR ARGUMENTS         ZSPX0300
      INTEGER            M,N,LDA,LIPVT,IPVT(LIPVT)                      ZSPX0310
      REAL               A(LDA,N),RDIAG(N),ACNORM(N),WA(N)              ZSPX0320
      LOGICAL            PIVOT                                          ZSPX0330
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   ZSPX0340
      INTEGER            I,JP1,J,KMAX,K,MINMN                           ZSPX0350
      REAL               AJNORM,EPSMCH,ONE,P05,SPMPAR,SUM,TEMP,ZERO     ZSPX0360
      REAL               SNRM2                                          ZSPX0370
      DATA               SPMPAR /1.1921E-07/                            ZSPX0380
c      DATA	             SPMPAR /2.2204E-16/
      DATA               ONE,P05,ZERO /1.0E0,5.0E-2,0.0E0/              ZSPX0390
C                                  EPSMCH IS THE MACHINE PRECISION.     ZSPX0400
C                                  FIRST EXECUTABLE STATEMENT           ZSPX0410
      EPSMCH = SPMPAR                                                   ZSPX0420
C                                  COMPUTE THE INITIAL COLUMN NORMS AND ZSPX0430
C                                  INITIALIZE SEVERAL ARRAYS.           ZSPX0440
      DO 5 J=1,N                                                        ZSPX0450
         ACNORM(J) = SNRM2(M,A(1,J),1)                                  ZSPX0460
         RDIAG(J) = ACNORM(J)                                           ZSPX0470
         WA(J) = RDIAG(J)                                               ZSPX0480
         IF (PIVOT) IPVT(J) = J                                         ZSPX0490
    5 CONTINUE                                                          ZSPX0500
C                                  REDUCE A TO R WITH HOUSEHOLDER       ZSPX0510
C                                  TRANSFORMATIONS.                     ZSPX0520
      MINMN = MIN0(M,N)                                                 ZSPX0530
      DO 55 J=1,MINMN                                                   ZSPX0540
         IF (.NOT.PIVOT) GO TO 20                                       ZSPX0550
C                                  BRING THE COLUMN OF LARGEST NORM INTOZSPX0560
C                                  THE PIVOT POSITION.                  ZSPX0570
         KMAX = J                                                       ZSPX0580
         DO 10 K=J,N                                                    ZSPX0590
            IF (RDIAG(K).GT.RDIAG(KMAX)) KMAX = K                       ZSPX0600
   10    CONTINUE                                                       ZSPX0610
         IF (KMAX.EQ.J) GO TO 20                                        ZSPX0620
         DO 15 I=1,M                                                    ZSPX0630
            TEMP = A(I,J)                                               ZSPX0640
            A(I,J) = A(I,KMAX)                                          ZSPX0650
            A(I,KMAX) = TEMP                                            ZSPX0660
   15    CONTINUE                                                       ZSPX0670
         RDIAG(KMAX) = RDIAG(J)                                         ZSPX0680
         WA(KMAX) = WA(J)                                               ZSPX0690
         K = IPVT(J)                                                    ZSPX0700
         IPVT(J) = IPVT(KMAX)                                           ZSPX0710
         IPVT(KMAX) = K                                                 ZSPX0720
   20    CONTINUE                                                       ZSPX0730
C                                  COMPUTE THE HOUSEHOLDER              ZSPX0740
C                                  TRANSFORMATION TO REDUCE THE J-TH    ZSPX0750
C                                  COLUMN OF A TO A MULTIPLE OF THE J-THZSPX0760
C                                  UNIT VECTOR.                         ZSPX0770
         AJNORM = SNRM2(M-J+1,A(J,J),1)                                 ZSPX0780
         IF (AJNORM.EQ.ZERO) GO TO 50                                   ZSPX0790
         IF (A(J,J).LT.ZERO) AJNORM = -AJNORM                           ZSPX0800
         DO 25 I=J,M                                                    ZSPX0810
            A(I,J) = A(I,J)/AJNORM                                      ZSPX0820
   25    CONTINUE                                                       ZSPX0830
         A(J,J) = A(J,J)+ONE                                            ZSPX0840
C                                  APPLY THE TRANSFORMATION TO THE      ZSPX0850
C                                  REMAINING COLUMNS AND UPDATE THE     ZSPX0860
C                                  NORMS.                               ZSPX0870
         JP1 = J+1                                                      ZSPX0880
         IF (N.LT.JP1) GO TO 50                                         ZSPX0890
         DO 45 K=JP1,N                                                  ZSPX0900
            SUM = ZERO                                                  ZSPX0910
            DO 30 I=J,M                                                 ZSPX0920
               SUM = SUM+A(I,J)*A(I,K)                                  ZSPX0930
   30       CONTINUE                                                    ZSPX0940
            TEMP = SUM/A(J,J)                                           ZSPX0950
            DO 35 I=J,M                                                 ZSPX0960
               A(I,K) = A(I,K)-TEMP*A(I,J)                              ZSPX0970
   35       CONTINUE                                                    ZSPX0980
            IF (.NOT.PIVOT .OR. RDIAG(K).EQ.ZERO) GO TO 40              ZSPX0990
            TEMP = A(J,K)/RDIAG(K)                                      ZSPX1000
            RDIAG(K) = RDIAG(K)*SQRT(AMAX1(ZERO,ONE-TEMP**2))           ZSPX1010
            IF (P05*(RDIAG(K)/WA(K))**2.GT.EPSMCH) GO TO 40             ZSPX1020
            RDIAG(K) = SNRM2(M-J,A(JP1,K),1)                            ZSPX1030
            WA(K) = RDIAG(K)                                            ZSPX1040
   40       CONTINUE                                                    ZSPX1050
   45    CONTINUE                                                       ZSPX1060
   50    CONTINUE                                                       ZSPX1070
         RDIAG(J) = -AJNORM                                             ZSPX1080
   55 CONTINUE                                                          ZSPX1090
      RETURN                                                            ZSPX1100
      END                                                               ZSPX1110
