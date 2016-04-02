      subroutine epflag(matype,ipfl,ipelem,lpb,numel)
c     implicit double precision (a-h,o-z)                                    dp
c
c     set element print flags
c
      dimension ipelem(2,*),matype(*),ipfl(*)
      do 20 i=1,numel
      do 10 j=1,lpb
      k=ipelem(1,j)
      l=ipelem(2,j)
      if (i.ge.k.and.i.le.l) go to 20
   10 continue
      ipfl(i)=1
   20 continue
c
      return
c
      end
