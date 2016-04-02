      subroutine addmas(lm,ym,gm)
c     implicit double precision (a-h,o-z)                                    dp
c
      dimension lm(*),ym(*),gm(*)
      do 10 i=1,4
      ny=lm(2*i-1)
      nz=lm(2*i)
      if (ny.ne.0) gm(ny)=gm(ny)+ym(i)
      if (nz.ne.0) gm(nz)=gm(nz)+ym(i)
   10 continue
      return
      end
