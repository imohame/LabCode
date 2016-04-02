      SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT)

      PARAMETER (NMAX  = 24,
     >           BIG   = 1.00E+03,
     >           SMALL = 1.00E-03 )
      DIMENSION Y(NMAX),DYDX(NMAX),YOUT(NMAX),YT(NMAX),
     >          DYT(NMAX),DYM(NMAX)

      HH = H*0.5
      H6 = H/6.00
      XH = X + HH
	  !write(*,*) 'in rk4'
	  !write(*,*) N
      DO 11 I=1,N
         YT(I)=Y(I)+HH*DYDX(I)
11    CONTINUE
      CALL DERIVS(XH,YT,DYT,N)
      DO 12 I=1,N
         YT(I)=Y(I)+HH*DYT(I)
12    CONTINUE
      CALL DERIVS(XH,YT,DYM,N)
      DO 13 I=1,N
         YT(I) = Y(I) + H*DYM(I)
         DYM(I) = DYT(I) + DYM(I)
13    CONTINUE
      CALL DERIVS(X+H,YT,DYT,N)
      DO 14 I=1,N
         YOUT(I) = Y(I) + H6*(DYDX(I) + DYT(I) + 2.00*DYM(I))
14    CONTINUE


      RETURN
      END
