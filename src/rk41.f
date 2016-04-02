           SUBROUTINE RK41(Y1,DYDX1,N1,X1,H1,YOUT1)
           PARAMETER(BIG1=1.00E+03,SMALL1=1.00E-03)
           DIMENSION Y1(N1),DYDX1(N1),YOUT1(N1),YT1(N1),DYT1(N1),
     > DYM1(N1)
         HH1=H1*0.5
         H61=H1/6.00
         XH1=X1+HH1
c		 write(*,*) 'in rk41, N1',N1
          DO 11 I=1,N1
           YT1(I)=Y1(I)+HH1*DYDX1(I)
 11        CONTINUE
           CALL derivs1(XH1,YT1,DYT1,N1)
           DO 12 I=1,N1
           YT1(I)=Y1(I)+HH1*DYT1(I)
 12         CONTINUE
           CALL derivs1(XH1,YT1,DYM1,N1)
           DO 13 I=1,N1
           YT1(I)=Y1(I)+H1*DYM1(I)
           DYM1(I)=DYT1(I)+DYM1(I)
 13         CONTINUE
           CALL derivs1(X1+H1,YT1,DYT1,N1)
           DO 14 I=1,N1
           YOUT1(I)=Y1(I)+H61*(DYDX1(I)+DYT1(I)+2.00*DYM1(I))
 14         CONTINUE
           RETURN
             END
