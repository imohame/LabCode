C   IMSL ROUTINE NAME   - ZSPOW                                         ZSPO0010
C                                                                       ZSPO0020
C-----------------------------------------------------------------------ZSPO0030
C                                                                       ZSPO0040
C   COMPUTER            - VAX/SINGLE                                    ZSPO0050
C                                                                       ZSPO0060
C   LATEST REVISION     - JULY 1, 1983                                  ZSPO0070
C                                                                       ZSPO0080
C   PURPOSE             - SOLVE A SYSTEM OF NONLINEAR EQUATIONS         ZSPO0090
C                                                                       ZSPO0100
C   USAGE               - CALL ZSPOW (FCN,NSIG,N,ITMAX,PAR,X,FNORM,     ZSPO0110
C                           WK,IER)                                     ZSPO0120
C                                                                       ZSPO0130
C   ARGUMENTS    FCN    - THE NAME OF A USER-SUPPLIED SUBROUTINE WHICH  ZSPO0140
C                           EVALUATES THE SYSTEM OF EQUATIONS TO BE     ZSPO0150
C                           SOLVED. FCN MUST BE DECLARED EXTERNAL IN    ZSPO0160
C                           THE CALLING PROGRAM AND MUST HAVE THE       ZSPO0170
C                           FOLLOWING FORM,                             ZSPO0180
C                             SUBROUTINE FCN(X,F,N,PAR)                 ZSPO0190
C                             REAL X(N),F(N),PAR(1)                     ZSPO0200
C                             F(1)=                                     ZSPO0210
C                              .                                        ZSPO0220
C                             F(N)=                                     ZSPO0230
C                             RETURN                                    ZSPO0240
C                             END                                       ZSPO0250
C                           GIVEN X(1)...X(N), FCN MUST EVALUATE THE    ZSPO0260
C                           FUNCTIONS F(1)...F(N) WHICH ARE TO BE MADE  ZSPO0270
C                           ZERO. X SHOULD NOT BE ALTERED BY FCN. THE   ZSPO0280
C                           PARAMETERS IN VECTOR PAR (SEE ARGUMENT      ZSPO0290
C                           PAR BELOW) MAY ALSO BE USED IN THE          ZSPO0300
C                           CALCULATION OF F(1)...F(N).                 ZSPO0310
C                NSIG   - THE NUMBER OF DIGITS OF ACCURACY DESIRED      ZSPO0320
C                           IN THE COMPUTED ROOT. (INPUT)               ZSPO0330
C                N      - THE NUMBER OF EQUATIONS TO BE SOLVED AND      ZSPO0340
C                           THE NUMBER OF UNKNOWNS. (INPUT)             ZSPO0350
C                ITMAX  - THE MAXIMUM ALLOWABLE NUMBER OF ITERATIONS.   ZSPO0360
C                           (INPUT) THE MAXIMUM NUMBER OF CALLS TO FCN  ZSPO0370
C                           IS ITMAX*(N+1). SUGGESTED VALUE = 200.      ZSPO0380
C                PAR    - PAR CONTAINS A PARAMETER SET WHICH IS         ZSPO0390
C                           PASSED TO THE USER-SUPPLIED FUNCTION FCN.   ZSPO0400
C                           PAR MAY BE USED TO PASS ANY AUXILIARY       ZSPO0410
C                           PARAMETERS NECESSARY FOR COMPUTATION OF     ZSPO0420
C                           THE FUNCTION FCN. (INPUT)                   ZSPO0430
C                X      - A VECTOR OF LENGTH N. (INPUT/OUTPUT) ON INPUT,ZSPO0440
C                           X IS THE INITIAL APPROXIMATION TO THE ROOT. ZSPO0450
C                           ON OUTPUT, X IS THE BEST APPROXIMATION TO   ZSPO0460
C                           THE ROOT FOUND BY ZSPOW.                    ZSPO0470
C                FNORM  - ON OUTPUT, FNORM IS EQUAL TO                  ZSPO0480
C                           F(1)**2+...F(N)**2 AT THE POINT X.          ZSPO0490
C                WK     - WORK VECTOR OF LENGTH N*(3*N+15)/2            ZSPO0500
C                IER    - ERROR PARAMETER. (OUTPUT)                     ZSPO0510
C                         TERMINAL ERROR                                ZSPO0520
C                           IER = 129 INDICATES THAT THE NUMBER OF      ZSPO0530
C                             CALLS TO FCN HAS EXCEEDED ITMAX*(N+1).    ZSPO0540
C                             THE USER MAY TRY A NEW INITIAL GUESS.     ZSPO0550
C                           IER = 130 INDICATES THAT NSIG IS TOO        ZSPO0560
C                             LARGE.  NO FURTHER IMPROVEMENT IN THE     ZSPO0570
C                             APPROXIMATE SOLUTION X IS POSSIBLE.       ZSPO0580
C                             THE USER SHOULD DECREASE NSIG.            ZSPO0590
C                           IER = 131 INDICATES THAT THE ITERATION      ZSPO0600
C                             HAS NOT MADE GOOD PROGRESS.  THE USER     ZSPO0610
C                             MAY TRY A NEW INITIAL GUESS.              ZSPO0620
C                                                                       ZSPO0630
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         ZSPO0640
C                       - SINGLE/H36,H48,H60                            ZSPO0650
C                                                                       ZSPO0660
C   REQD. IMSL ROUTINES - SINGLE/uertst,UGETIO,VBLA=SNRM2,ZSPWA,        ZSPO0670
C                           ZSPWB,ZSPWC,ZSPWD,ZSPWE,ZSPWF,ZSPWG         ZSPO0680
C                       - DOUBLE/uertst,UGETIO,VBLA=DNRM2,ZSPWA,        ZSPO0690
C                           ZSPWB,ZSPWC,ZSPWD,ZSPWE,ZSPWF,ZSPWG         ZSPO0700
C                                                                       ZSPO0710
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           ZSPO0720
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      ZSPO0730
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  ZSPO0740
C                                                                       ZSPO0750
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.       ZSPO0760
C                                                                       ZSPO0770
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN ZSPO0780
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    ZSPO0790
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        ZSPO0800
C                                                                       ZSPO0810
C-----------------------------------------------------------------------ZSPO0820
C                                                                       ZSPO0830
      SUBROUTINE ZSPOW1 (FCN1,NSIG,N,ITMAX,PAR,X,FNORM,WK,IER)          ZSPO0840
C                                  SPECIFICATIONS FOR ARGUMENTS         ZSPO0850
      INTEGER            NSIG,N,ITMAX,IER                               ZSPO0860
      REAL               PAR(N+1),X(N),FNORM,WK(N*(3*N+15)/2)                       ZSPO0870
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   ZSPO0880
      INTEGER            INDEX2,INDEX,INFO,I,J,LR,MAXFEV,ML,MODE,MU,    ZSPO0890
     1                   NFEV,NPRINT                                    ZSPO0900
      REAL               EPSFCN,FACTOR,ONE,XTOL,ZERO                    ZSPO0910
      EXTERNAL           FCN1                                            ZSPO0920
      DATA               FACTOR,ONE,ZERO /1.0E2,1.0E0,0.0E0/            ZSPO0930
C                                  FIRST EXECUTABLE STATEMENT           ZSPO0940
      INFO = 0                                                          ZSPO0950
C                                  CALL ZSPWA                           ZSPO0960
      MAXFEV = ITMAX*(N + 1)                                            ZSPO0970
      XTOL = 0.1**NSIG                                                  ZSPO0980
      ML = N - 1                                                        ZSPO0990
      MU = N - 1                                                        ZSPO1000
      EPSFCN = ZERO                                                     ZSPO1010
      MODE = 2                                                          ZSPO1020
      DO 5 J = 1, N                                                     ZSPO1030
         WK(J) = ONE                                                    ZSPO1040
    5 CONTINUE                                                          ZSPO1050
      NPRINT = 0                                                        ZSPO1060
      LR = (N*(N + 1))/2                                                ZSPO1070
      INDEX = 7*N + LR                                                  ZSPO1080
      CALL ZSPWA1(FCN1,N,X,WK(6*N+1),XTOL,MAXFEV,ML,MU,EPSFCN,WK(1),    ZSPO1090
     * MODE,FACTOR,NPRINT,INFO,NFEV,WK(INDEX+1),N,WK(7*N+1),LR,         ZSPO1100
     * WK(N+1),WK(2*N+1),WK(3*N+1),WK(4*N+1),WK(5*N+1),PAR)             ZSPO1110
      IF (INFO .EQ. 5) INFO = 4                                         ZSPO1120
      FNORM = 0.0                                                       ZSPO1130
      DO 10 I=1,N                                                       ZSPO1140
         INDEX2 = 6*N+I                                                 ZSPO1150
         FNORM = FNORM+WK(INDEX2)*WK(INDEX2)                            ZSPO1160
   10 CONTINUE                                                          ZSPO1170
      IER = 0                                                           ZSPO1180
      IF (INFO .EQ. 2) IER = 129                                        ZSPO1190
      IF (INFO .EQ. 3) IER = 130                                        ZSPO1200
      IF (INFO .EQ. 4) IER = 131                                        ZSPO1210
      IF (IER .GT. 0) CALL uertst(IER,6HZSPOW )                         ZSPO1220
      write(*,15) ier
   15 format(5x,'IER1 =',i3/)
      RETURN                                                            ZSPO1230
      END                                                               ZSPO1240
