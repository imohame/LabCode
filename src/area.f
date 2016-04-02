      function area(r,z,ix)
c     implicit double precision (a-h,o-z)                                    dp
      dimension r(*),z(*),ix(*),xdr(4),xdz(4)
      do 10 i=1,4
      xdr(i)=r(ix(i))
   10 xdz(i)=z(ix(i))
      zd42=xdz(4)-xdz(2)
      ac1=xdr(2)*(xdz(3)-xdz(4))+xdr(3)*zd42+xdr(4)*(xdz(2)-xdz(3))
      ac2=xdr(2)*(xdz(4)-xdz(1))+xdr(4)*(xdz(1)-xdz(2))-xdr(1)*zd42
      area=ac1+ac2
      return
      end
