      subroutine forms(wv,v,s,neq,nv)

c     implicit double precision (a-h,o-z)                                    dp
      common/main_block/ a(1)
      dimension wv(*),v(*),s(*)

      k54=igtpnt(54)
      k16=igtpnt(16)
      do 20 i=1,nv
      j=neq*(i-1)
      call blkcpy (a(k54+j),wv,neq)
      do 10 j=i,nv
      m=j*(j-1)/2+i
      k=neq*(j-1)
      call blkcpy (a(k16+k),v,neq)
   10 s(m)=fdot(wv,v,neq)
   20 continue
      return
      end
