       SUBROUTINE RKQC(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,HMIN)

       PARAMETER (NMAX   = 24,
     >            PGROW  = -.20,
     >            PSHRNK = -0.25,
     >            FCOR   = 1./15.,
     >            ONE    = 1.,
     >            SAFETY = 0.9,
     >            ERRCON = 6.0E-04,
     >            BOUND  = 1.000E+08 )

       DIMENSION DYDX(NMAX),YSCAL(NMAX),YTEMP(NMAX),
     >           DYSAV(NMAX),YSAV(NMAX)

       real Y(NMAX)

C SAVE INITIAL VALUES
           XSAV=X
           DO 11 I=1,N
            YSAV(I)=Y(I)
            DYSAV(I)=DYDX(I)
11        CONTINUE
          H = HTRY

C   TWO HALF STEPS
1         HH=0.5*H
	      !write(*,*) 'in rkqc before rk4'
		  !write(*,*) N
          CALL RK4(YSAV,DYSAV,N,XSAV,HH,YTEMP)
		  !write(*,*) 'in rkqc after rk4'
		  !write(*,*) N
          X=XSAV+HH
          CALL DERIVS(X,YTEMP,DYDX,N)
		  !write(*,*) 'in rkqc after derivs'
		  !write(*,*) N
          CALL RK4(YTEMP,DYDX,N,X,HH,Y)
          X=XSAV+H

C   TAKE THE LARGE STEP
          CALL RK4(YSAV,DYSAV,N,XSAV,H,YTEMP)
          ERRMAX = 0.0

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
121        format(5x,'errmax =',e15.9/)
           IF(ERRMAX .GT. ONE)THEN
C
C  TRUNCUATION ERROR TOO LARGE REDUCE STEPSIZE
            WRITE(*,*)'reducing stepsize', H
             H=SAFETY*H*(ERRMAX**PSHRNK)
c  CHECK  FORsSTIFFNESS AND POSSIBLE STABILITY PROBLEMS
ck           write(27,122) h,hmin
122          format(5x,'H =',e15.9,2x,'HMIN =',e15.9/) 


             IF (H .LE. HMIN)THEN
			   write(*,*) 'calling implicit',N
                CALL IMPLICIT(YSAV(1:N),XSAV,HTRY,N)
c                WRITE(27,31)
c31              FORMAT('STIFFNESS DETECTED')
                X = XSAV+HTRY
                DO 100 I=1,N
                   Y(I)=YSAV(I)
100             CONTINUE
                HDID = HTRY
                HNEXT= HTRY
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


























