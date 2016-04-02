        SUBROUTINE FCN2(X2,f2,n2,par2)
        common/density/denim(2),denm(2),gdot(2)
        DIMENSION x2(n2),f2(n2),par2(3),yprime2(2)
        b=3.0e-10
        ak1=0.0127/b
        ak2=6.67
        ak3=5.53
        ak4=0.0000276/(b*b)
        yprime2(1)=abs(gdot(2))*(ak1*(x2(1)**0.5)+ak3*x2(2)-ak2*x2(1))
        yprime2(2)=abs(gdot(2))*(ak4*x2(1)/x2(2)-ak1*(x2(1)**0.5)
     1-ak3*x2(2))
        f2(1)=x2(1)-par2(1)-par2(3)*yprime2(1)
        f2(2)=x2(2)-par2(2)-par2(3)*yprime2(2)
             RETURN
             END
