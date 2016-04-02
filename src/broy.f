      subroutine broy(d,w,v,step)
c     implicit double precision (a-h,o-z)                                    dp
      common/bk12/ntlen
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk35/nbroy,numupd
      dimension d(*),w(*),v(*)
      common/main_block/ a(1)

      k54=igtpnt(54)
      stepi=1./step
      do 10 i=1,neq
   10 d(i)=step*d(i)
      call blkcpi (0,d,a,k54+nbroy,neq)
      nbroy=nbroy+neq
      call blkcpi (0,w,a,k54+nbroy+neq,neq)
      call blkcpy (w,d,neq)
      call bsolvr (d,a(ntlen),3,8)
      if (numupd.eq.0) go to 30
      do 20 i=1,numupd
      call blkcpi (1,a,w,k54+2*(i-1)*neq,2*neq)
      dw=dotprd(d,w)
      do 20 j=1,neq
   20 d(j)=d(j)+dw*v(j)
   30 nbroy=2*neq*numupd
      numupd=numupd+1
      call blkcpi (1,a,w,k54+nbroy,neq)
      nbroy=nbroy+neq
      do 40 i=1,neq
   40 v(i)=stepi*w(i)-d(i)
      fac=1./dotprd(v,w)
      do 50 i=1,neq
   50 v(i)=(w(i)-v(i))*fac
      call blkcpi (0,v,a,k54+nbroy,neq)
      nbroy=nbroy+neq
      dw=dotprd(d,w)
      do 60 i=1,neq
   60 d(i)=d(i)+dw*v(i)
      return
      end
