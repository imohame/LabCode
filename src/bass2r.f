      subroutine bass2r (h,det,xx,ipt)
c     implicit double precision (a-h,o-z)                                    dp
c
      common/bk01/h4(4,5),p14(4,5),p24(4,5)
      dimension h(*),p(2,4),xj(2,2),xx(2,*)
c
      do 10 i=1,4
      h(i)=h4(i,ipt)
      p(1,i)=p14(i,ipt)
   10 p(2,i)=p24(i,ipt)
      do 30 i=1,2
      do 30 j=1,2
      xj(i,j)=0.0
      do 20 k=1,4
   20 xj(i,j)=xj(i,j)+p(i,k)*xx(j,k)
   30 continue
      det=xj(1,1)*xj(2,2)-xj(2,1)*xj(1,2)
      return
      end
