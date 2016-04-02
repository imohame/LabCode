      subroutine intvec(bl,v,neq,nv)

c     implicit double precision (a-h,o-z)                                    dp
      common/main_block/ a(1)
      dimension bl(*),v(*)

      k16=igtpnt(16)
      j=1
      do 20 i=1,nv
      k=i*neq/nv
      m=neq*(i-1)
      call pzero (v,neq)
      do 10 n=j,k
   10 v(n)=bl(n)
      if (fdot(v,v,neq).eq.0.0) v(j)=1.0
      call blkcpy (v,a(k16+m),neq)
   20 j=k+1
      return
      end
