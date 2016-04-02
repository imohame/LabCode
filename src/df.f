      function df (y,x)
c     implicit double precision (a-h,o-z)                                    dp
c
      common /cont/ had(4), ityp2d
      dimension y(2,4), x(2,4)
c
      ac1=+y(1,2)*(y(2,3)-y(2,4))+y(1,3)*(y(2,4)-y(2,2))+y(1,4)*(y(2,2)
     1 -y(2,3))
      ac2=+y(1,2)*(y(2,4)-y(2,1))+y(1,4)*(y(2,1)-y(2,2))+y(1,1)*(y(2,2)
     1 -y(2,4))
      au1=+x(1,2)*(x(2,3)-x(2,4))+x(1,3)*(x(2,4)-x(2,2))+x(1,4)*(x(2,2)
     1 -x(2,3))
      au2=+x(1,2)*(x(2,4)-x(2,1))+x(1,4)*(x(2,1)-x(2,2))+x(1,1)*(x(2,2)
     1 -x(2,4))
      if (ityp2d.ge.1) go to 20
      dfe=((y(1,2)+y(1,3)+y(1,4))*ac1+(y(1,1)+y(1,2)+y(1,4))*ac2)/((x(1,
     1 2)+x(1,3)+x(1,4))*au1+(x(1,1)+x(1,2)+x(1,4))*au2)
      return
   20 df=(ac1+ac2)/(au1+au2)
      return
      end
