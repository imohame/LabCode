      subroutine recol(j,is,il,ia,x,y,d)
c     implicit double precision (a-h,o-z)                                    dp
c
c     reduce the components of a column or load vector
c
c     note- an assembly language version of this routine
c     is available for the cray-1 and for the cdc-7600.
c
      dimension ia(*),x(*),y(*),d(*)
c
      k=0
      i2=ia(is-1)
      do 20 ii=is,il
      i1=i2+1
      i2=ia(ii)
      n=i2-i1
      if(j.ge.0) n=min(n,k+j)
      k=k+1
      if(j.lt.0) d(ii)=x(i2)
      if(n.eq.0) go to 20
      s=0.
      do 10 i=1,n
      s=s+x(i+i2-n-1)*y(i+k-n-1)
   10 continue
      y(k)=y(k)-s
   20 continue
      return
      end
