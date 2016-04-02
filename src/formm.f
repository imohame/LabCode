      subroutine formm(wv,v,s,xm,neq,nv)

c     implicit double precision (a-h,o-z)                                    dp
      common/main_block/ a(1)
      dimension wv(*),v(*),s(*),xm(*)

      k54=igtpnt(54)
      k16=igtpnt(16)
      do 30 j=1,nv
      l=neq*(j-1)
      call blkcpy (a(k16+l),v,neq)
      do 10 i=1,neq
   10 wv(i)=v(i)*xm(i)
      k=j*(j+1)/2
      do 20 i=j,nv
      m=neq*(i-1)
      call blkcpy (a(k16+m),v,neq)
      s(k)=fdot(wv,v,neq)
   20 k=k+i
   30 call blkcpy (wv,a(k16+l),neq)
      return
      end
