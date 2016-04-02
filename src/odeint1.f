      SUBROUTINE ODEINT1 (YSTART1,NVAR1,X11,X21,EPS1,H11,HMIN1,NOK1,
     >                    NBAD1)

      PARAMETER (MAXSTP1=100,TWO1=2.0,ZERO1=0.00,TINY1=1.0E-20)
      integer NVAR1
      COMMON /PATH1/KMAX1,KOUNT1,DXSAV1,XP1(200),YP1(10,200)
      DIMENSION YSTART1(NVAR1),YSCAL1(NVAR1),Y1(NVAR1),DYDX1(NVAR1)

      X1    = X11
      H1    = sign(H11,X21-X11)
      NOK1  = 0
      NBAD1 = 0
      KMAX1 = 0
c      write(*,*) 'in odeint1,NVAR1',NVAR1
      DO 11 I=1,NVAR1
         Y1(I) = YSTART1(I)
c		 write(*,*) YSTART1(I),Y1(I)
11    CONTINUE

C  ASSURE STORAGE OF FIRST STEP
      XSAV1 = X1 - DXSAV1*TWO1
      
C  TAKE AT MOST MAXSTP STEPS
      DO 16 NSTP1 = 1, MAXSTP1
c       write(*,*) 'in odeint1, calling derivs1'
         CALL derivs1(X1,Y1,DYDX1,NVAR1)										!This Line: DYDX1 = dRho_m(alpha)/dt & dRho_im(alpha)/dt 
c       write(*,*) 'in odeint1, after derivs1'    
         DO 12 I = 1, NVAR1
C  SCALING USED TO MONITOR ACCURACY. THIS CAN BE MODIFIED AS NEEDED
            YSCAL1(I) = abs(Y1(I)) + abs(H1*DYDX1(I))+TINY1
12       CONTINUE

         IF((X1+H1-X21)*(X1+H1-X11) .GT. ZERO1) H1=X21-X1

C IF STEP CAN OVERSHOOT END,CUT DOWN STEPSIZE
c         write(*,*) 'in odeint1, calling rkqc1,NVAR1',NVAR1
         CALL RKQC1(Y1,DYDX1,NVAR1,X1,H1,EPS1,YSCAL1,HDID1,HNEXT1,
     >              HMIN1)
c	     write(*,*) 'in odeint1, after rkqc1'


         IF ((X1-X21)*(X21-X11) .GE. ZERO1) THEN

C ARE WE DONE?
            DO 14 I = 1, NVAR1
               YSTART1(I)=Y1(I)
c			write(*,*) YSTART1(I),Y1(I),'end',NVAR1
14          CONTINUE

            RETURN
         ENDIF
         
16    CONTINUE

      IF(abs(HNEXT1) .LT. HMIN1) H1 = HNEXT1
      WRITE(*,31)
31    FORMAT(' STEPSIZE SMALLER THAN MINIMUM')

      
      RETURN
      END

