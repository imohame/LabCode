      subroutine bssf (h,p)
c     implicit double precision (a-h,o-z)                                    dp
c
      common/bk01/h4(4,5),p14(4,5),p24(4,5)
      common/intgrt/nintg
      dimension h(*),p(2,1)
      do 10 i=1,4
      h(i)=h4(i,nintg)
      p(1,i)=p14(i,nintg)
   10 p(2,i)=p24(i,nintg)
      return
      end
