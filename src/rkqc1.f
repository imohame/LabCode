       SUBROUTINE RKQC1(Y1,DYDX1,N1,X1,HTRY1,EPS1,YSCAL1,HDID1,HNEXT1,
     .hmin1)
          PARAMETER (PGROW1=-.20,PSHRNK1=-0.25,FCOR1=1./15.,
     1ONE1=1.,SAFETY1=0.9,ERRCON1=6.0E-04,BOUND1=1.000E+08)
ck           EXTERNAL derivs1
      common /test11/ ink
	  integer N1, ink
	  
      DIMENSION Y1(N1),DYDX1(N1),YSCAL1(N1),YTEMP1(N1),YSAV1(N1),
     1DYSAV1(N1)
C SAVE INITIAL VALUES
           XSAV1=X1
           DO 11 I=1,N1
           YSAV1(I)=Y1(I)
            DYSAV1(I)=DYDX1(I)
c			write(*,*) 'beginning',Y1(I)
11        CONTINUE
          H1=HTRY1
C   TWO HALF STEPS
1         HH1=0.5*H1
c          write(*,*) 'in rkqc1, calling rk41'
          CALL RK41(YSAV1,DYSAV1,N1,XSAV1,HH1,YTEMP1)
c		  write(*,*) 'in rkqc1, calling derivs1'
          X1=XSAV1+HH1
          CALL derivs1(X1,YTEMP1,DYDX1,N1)
c		  write(*,*) 'in rkqc1, calling rk41'
          CALL RK41(YTEMP1,DYDX1,N1,X1,HH1,Y1)
          X1=XSAV1+H1
C   TAKE THE LARGE STEP
c          write(*,*) 'in rkqc1, calling rk41'
          CALL RK41(YSAV1,DYSAV1,N1,XSAV1,H1,YTEMP1)
          ERRMAX1=0.0
C   EVALUATE ACCURACY
          DO 12 I=1,N1
C  YTEMP NOW CONTAINS THE ERROR ESTIMATE
           YTEMP1(I)=Y1(I)-YTEMP1(I)
c      WRITE(6,*)YTEMP(I)
           ERRMAX1=max(ERRMAX1,abs(YTEMP1(I)/YSCAL1(I)))
12         CONTINUE
C SCALE RELATIVE TO THE REQUIRE EPS(ACCURACY)
           ERRMAX1=ERRMAX1/EPS1
ck           write(27,121) errmax
121        format(5x,'errmax =',e17.9/)
           IF(ERRMAX1 .GT. ONE1)THEN
C
C  TRUNCUATION ERROR TOO LARGE REDUCE STEPSIZE
           WRITE(*,*) 'stiffness detected rkqc1',N1,ink
           H1=SAFETY1*H1*(ERRMAX1**PSHRNK1)
c  CHECK  FORsSTIFFNESS AND POSSIBLE STABILITY PROBLEMS
ck           write(27,122) h,hmin
122        format(5x,'H =',e17.9,2x,'HMIN =',e17.9/)
       IF(H1 .LE. HMIN1 )THEN
	   write(*,*) 'calling implicit1'
       CALL IMPLICIT1(YSAV1,XSAV1,HTRY1,N1)
c          WRITE(27,31)
31       FORMAT('stiffness detected 1')
            X1=XSAV1+HTRY1
ck            WRITE(27,*)H,ERRMAX,x
           DO 100 I=1,N1
            Y1(I)=YSAV1(I)
100      CONTINUE
            HDID1=HTRY1
             HNEXT1=HTRY1
            RETURN
         ENDIF
           GO TO 1
C  GO FOR ANOTHER TRY
           ELSE
C   STEP SUCCEEDED
           HDID1=H1
           IF(ERRMAX1 .GT. ERRCON1)THEN
            HNEXT1=SAFETY1*H1*(ERRMAX1**PGROW1)
           ELSE
             HNEXT1=4.00*H1
C  INCREASE STEPSIZE
           ENDIF
           ENDIF
           DO 13 I=1,N1
C  MOP UP FIFTH ORDER TRUNCUATION ERROR
           Y1(I)=Y1(I)+YTEMP1(I)*FCOR1
c		   write(*,*) 'end',Y1(I)
13         CONTINUE
99         CONTINUE
            RETURN
            END


























