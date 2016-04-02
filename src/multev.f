      subroutine multev(d,p,nv)
c     implicit double precision (a-h,o-z)                                    dp
      dimension d(*),p(nv,1)
      do 10 i=1,nv
      do 10 j=1,nv
   10 p(j,i)=p(j,i)/d(i)
      return
      end
