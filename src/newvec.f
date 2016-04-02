      subroutine newvec(v,p,t,nv,neq)
c     implicit double precision (a-h,o-z)                                    dp
      common/main_block/ a(1)
      dimension v(*),p(nv,1),t(*)

      k54=igtpnt(54)
      k16=igtpnt(16)
      do 20 i=1,nv
      call pzero (v,neq)
      do 10 j=1,nv
      l=neq*(j-1)
      call blkcpy (a(k54+l),t,neq)
      do 10 k=1,neq
   10 v(k)=v(k)+t(k)*p(j,i)
      l=neq*(i-1)
   20 call blkcpy (v,a(k16+l),neq)
      return
      end
