       subroutine m2derivs(x,y,yprime)
       common /myvar/ a,b,c,e,taur,twomu,xmhigh1,xmhigh2
       common /myvar1/ gdotr(2),tau(2)
       dimension y(2),yprime(2)
       sat=1.0/1.00
       ratio1=((abs(y(1)/taur)*sat)**(xmhigh1-1.))*(y(1)/taur)*sat
       ratio2=((abs(y(2)/taur)*sat)**(xmhigh2-1.))*(y(2)/taur)*sat
       yprime(1)=twomu*(gdotr(1)-a*ratio1-b*ratio2)
       yprime(2)=twomu*(gdotr(2)-c*ratio1-e*ratio2)
       return
       end
