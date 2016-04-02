      function dfe (ix,r,z,x)
c     implicit double precision (a-h,o-z)                                    dp
c
      character*6 probid
      common /ident/ probid(12)
      common /cont/ had(4), ityp2d
      dimension ix(*), r(*), z(*), x(2,*), xdr(4), xdz(4), xur(4), xuz(4
     1 )
c
      do 10 i=1,4
      xur(i)=x(1,ix(i))
      xuz(i)=x(2,ix(i))
      xdr(i)=r(ix(i))
   10 xdz(i)=z(ix(i))
      ac1=+xdr(2)*(xdz(3)-xdz(4))+xdr(3)*(xdz(4)-xdz(2))+xdr(4)*(xdz(2)
     1 -xdz(3))
      ac2=+xdr(2)*(xdz(4)-xdz(1))+xdr(4)*(xdz(1)-xdz(2))+xdr(1)*(xdz(2)
     1 -xdz(4))
      au1=+xur(2)*(xuz(3)-xuz(4))+xur(3)*(xuz(4)-xuz(2))+xur(4)*(xuz(2)
     1 -xuz(3))
      au2=+xur(2)*(xuz(4)-xuz(1))+xur(4)*(xuz(1)-xuz(2))+xur(1)*(xuz(2)
     1 -xuz(4))
      if (ityp2d.ge.1) go to 20
      dfe=((xdr(2)+xdr(3)+xdr(4))*ac1+(xdr(1)+xdr(2)+xdr(4))*ac2)/((xur(
     1 2)+xur(3)+xur(4))*au1+(xur(1)+xur(2)+xur(4))*au2)
      return
   20 dfe=(ac1+ac2)/(au1+au2)
      return
      end
