      function arclen (r,z,crv,n)
c     implicit double precision (a-h,o-z)                                    dp
      dimension crv(2,*)
      tol=.001
      arclen=0.
      if (abs(r-crv(1,1)).le..001.and.abs(z-crv(2,1)).le..001) return
      do 10 i=2,n
      r1=crv(1,i-1)
      z1=crv(2,i-1)
      r2=crv(1,i)
      z2=crv(2,i)
      dl=sqrt((r2-r1)**2+(z2-z1)**2)
      if (r.lt.min(r1,r2)-tol) go to 10
      if (r.gt.max(r1,r2)+tol) go to 10
      if (z.lt.min(z1,z2)-tol) go to 10
      if (z.gt.max(z1,z2)+tol) go to 10
      go to 20
   10 arclen=arclen+dl
      return
   20 dl=sqrt((r-r1)**2+(z-z1)**2)
      arclen=arclen+dl
      return
      end
