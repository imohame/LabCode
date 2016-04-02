      subroutine basisp
c     implicit double precision (a-h,o-z)                                    dp
      common/bk30/numlp,numpc,h22(2,2),pl2(2,2),h33(3,2),pl3(3,2)
      dimension s(2)
      data s/-.5773502691895,.5773502691896/
      do 10 j=1,2
      h22(1,j)=.5*(1.0-s(j))
      h22(2,j)=.5*(1.0+s(j))
      h33(1,j)=.5*s(j)*(s(j)-1.0)
      h33(2,j)=.5*s(j)*(s(j)+1.0)
      h33(3,j)=1.0-s(j)*s(j)
      pl2(1,j)=-.50
      pl2(2,j)=.50
      pl3(1,j)=s(j)-.50
      pl3(2,j)=s(j)+.50
      pl3(3,j)=-2.*s(j)
   10 continue
      return
      end
