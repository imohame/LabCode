      subroutine psorte(d,v,n,eivec)
c     implicit double precision (a-h,o-z)                                    dp
      logical eivec
      dimension d(*),v(n,1)
      n1=n-1
      do 30 i=1,n1
      e=d(i)
      ii=i
      i1=i+1
      do 10 j=i1,n
      if (abs(e).gt.abs(d(j))) go to 10
      e=d(j)
      ii=j
   10 continue
      d(ii)=d(i)
      d(i)=e
      if (.not.eivec) go to 30
      do 20 j=1,n
      e=v(j,i)
      v(j,i)=v(j,ii)
   20 v(j,ii)=e
   30 continue
      return
      end
