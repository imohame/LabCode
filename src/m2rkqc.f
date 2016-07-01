      SUBROUTINE M2RKQC(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS,
     $HMIN)
          PARAMETER (NMAX=2,PGROW=-.20,PSHRNK=-0.25,FCOR=1./15.,
     1ONE=1.,SAFETY=0.9,ERRCON=6.0E-04,BOUND=1.000E+08)
           EXTERNAL DERIVS
        DIMENSION Y(N),DYDX(N),YSCAL(N),YTEMP(NMAX),YSAV(NMAX),
     1DYSAV(NMAX)
C SAVE INITIAL VALUES
           XSAV=X
           DO 11 I=1,N
           YSAV(I)=Y(I)
            DYSAV(I)=DYDX(I)
11        CONTINUE
          H=HTRY
C   TWO HALF STEPS
1         HH=0.5*H
          CALL m2rk4(YSAV,DYSAV,N,XSAV,HH,YTEMP,DERIVS)
          X=XSAV+HH
          CALL m2derivs(X,YTEMP,DYDX)
          CALL m2rk4(YTEMP,DYDX,N,X,HH,Y,DERIVS)
          X=XSAV+H
C   TAKE THE LARGE STEP
          CALL m2rk4(YSAV,DYSAV,N,XSAV,H,YTEMP,DERIVS)
          ERRMAX=0.0
C   EVALUATE ACCURACY
          DO 12 I=1,N
C  YTEMP NOW CONTAINS THE ERROR ESTIMATE
           YTEMP(I)=Y(I)-YTEMP(I)
c      WRITE(6,*)YTEMP(I)
           ERRMAX=max(ERRMAX,abs(YTEMP(I)/YSCAL(I)))
12         CONTINUE
C SCALE RELATIVE TO THE REQUIRE EPS(ACCURACY)
           ERRMAX=ERRMAX/EPS
ck           write(27,121) errmax
121        format(5x,'errmax =',e17.9/)
           IF(ERRMAX .GT. ONE)THEN
C
C  TRUNCUATION ERROR TOO LARGE REDUCE STEPSIZE
c           WRITE(6,*)H
           H=SAFETY*H*(ERRMAX**PSHRNK)
c  CHECK  FORsSTIFFNESS AND POSSIBLE STABILITY PROBLEMS
ck           write(27,122) h,hmin
122        format(5x,'H =',e17.9,2x,'HMIN =',e17.9/)
       IF(H .LE. HMIN )THEN
       CALL M2IMPLICIT(YSAV,XSAV,HTRY)
          WRITE(27,31)
31       FORMAT('STIFFNESS DETECTED')
            X=XSAV+HTRY
ck            WRITE(27,*)H,ERRMAX,x
           DO 100 I=1,N
            Y(I)=YSAV(I)
100      CONTINUE
            HDID=HTRY
             HNEXT=HTRY
            RETURN
         ENDIF
           GO TO 1
C  GO FOR ANOTHER TRY
           ELSE
C   STEP SUCCEEDED
           HDID=H
           IF(ERRMAX .GT. ERRCON)THEN
            HNEXT=SAFETY*H*(ERRMAX**PGROW)
           ELSE
             HNEXT=4.00*H
C  INCREASE STEPSIZE
           ENDIF
           ENDIF
           DO 13 I=1,N
C  MOP UP FIFTH ORDER TRUNCUATION ERROR
           Y(I)=Y(I)+YTEMP(I)*FCOR
13         CONTINUE
99         CONTINUE
            RETURN
            END


























