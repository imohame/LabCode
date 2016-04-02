      subroutine m2rk4(y,dydx,n,x,h,yout,derivs)
      parameter(nmax=2,big=1.00e+03,small=1.00e-03)
      dimension y(n),dydx(n),yout(n),yt(nmax),dyt(nmax),dym(nmax)
      hh=h*0.5
      h6=h/6.00
      xh=x+hh
      do 11 i=1,n
      yt(i)=y(i)+hh*dydx(i)
11    continue
      call m2derivs(xh,yt,dyt)
      do 12 i=1,n
      yt(i)=y(i)+hh*dyt(i)
12    continue
      call m2derivs(xh,yt,dym)
      do 13 i=1,n
      yt(i)=y(i)+h*dym(i)
      dym(i)=dyt(i)+dym(i)
13    continue
      call m2derivs(x+h,yt,dyt)
      do 14 i=1,n
      yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.00*dym(i))
14    continue
      return
      end
