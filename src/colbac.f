      subroutine colbac(ja,a,b)
c     implicit double precision (a-h,o-z)                                    dp
c
c     backsubstitute equations, by columns, from block a
c
c     note- an assembly language version of this routine
c     is available for the cray-1 and for the cdc-7600.
c
      dimension ja(*),a(*),b(*)
      common/fissl4/je,jc
c
      j=jc
      j1=ja(j)
  100 continue
      j=j-1
      if(j.lt.0) return
      d=-b(j+1)
      j2=j1-1
      j1=ja(j)
      n=j2-j1
      if((d.eq.0.).or.(n.eq.0)) go to 100
      do 120 i=1,n
      b(i+j-n)=b(i+j-n)+d*a(i+j1)
  120 continue
      go to 100
      end
