      subroutine yield (eps,sige,ep,ak,qh)
c     implicit double precision (a-h,o-z)                                    dp
c
c     scalar subroutine for interpolating tabulated curve
c
      common/pnts/npoint
      dimension eps(*),sige(1)
c
      do 10 i=2,npoint
      m=i-1
      if (ep.gt.eps(i)) go to 10
      qh=(sige(i)-sige(i-1))/(eps(i)-eps(i-1))
      ak=sige(i-1)+qh*(ep-eps(i-1))
      go to 20
   10 continue
      m=m+1
      qh=(sige(m)-sige(m-1))/(eps(m)-eps(m-1))
      ak=sige(m)+qh*(ep-eps(m))
   20 continue
c
      return
c
      end
