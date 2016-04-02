      subroutine dfp(a,d,w,v,sp,neq,ntlen,nbfgs,numupd)
c     implicit double precision (a-h,o-z)                                    dp
      dimension a(*),d(*),w(*),v(*),fac2(201),alp(201)
      k13=igtpnt(13)
      k54=igtpnt(54)
      call blkcpi (1,a,v,k54+nbfgs,neq)
      call blkcpy (w,a(k13),neq)
      do 10 i=1,neq
      v(i)=v(i)-w(i)
   10 w(i)=sp*d(i)
      call blkcpy (a(k13),d,neq)
      call blkcpi (0,w,a,k54+nbfgs,3*neq)
      call bsolvr (d,a(ntlen),3,8)
      if (numupd.gt.0) then
      do 30 j=1,numupd
      locat=k54+2*neq*(j-1)
      call blkcpi (1,a,w,locat,2*neq)
      alfa0=fac2(j)*dotprd(v,a(k13))
      alfa1= alp(j)*dotprd(w,a(k13))
      do 20 i=1,neq
   20 d(i)=d(i)-alfa1*w(i)+alfa0*v(i)
   30 continue
      endif
      numupd=numupd+1
      call blkcpi (1,a,w,k54+nbfgs,2*neq)
      call blkcpy (w,a(k13),neq)
      do 40 i=1,neq
   40 w(i)=w(i)/sp-d(i)
      alp(numupd) =1./dotprd(w,v)
      fac2(numupd)=1./dotprd(a(k13),v)
      call blkcpi (0,w,a,k54+nbfgs,neq)
      call blkcpi (0,a(k13),a,k54+nbfgs+neq,neq)
      nbfgs=nbfgs+2*neq
      call blkcpi (1,a,v,k54+nbfgs,neq)
      alfa0=fac2(numupd)*dotprd(a(k13),v)
      alfa1=alp(numupd)*dotprd(w,v)
      do 50 i=1,neq
   50 d(i)=d(i)-alfa1*w(i)+alfa0*a(k13+i-1)
      return
      end
