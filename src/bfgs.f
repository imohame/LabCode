      subroutine bfgs(d,w,v,step,g,g0,dnorm)
c     implicit double precision (a-h,o-z)                                    dp
      common/bk12/ntlen
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk35/nbfgs,numupd
      dimension d(*),w(*),v(*)
      common/main_block/ a(1)
	  common/mbsize/numelt2, numnp2, neq2

      k13=igtpnt(13)
      k54=igtpnt(54)
      deltag=step*(g0-g)
      deltak=step*step*g0
      if (deltag.gt.0.0.and.deltak.gt.0.0) go to 10
      call blkcpy (w,d,neq)
      go to 30
   10 nupd=1
      fact1=1.0+step*sqrt(deltag/deltak)
      fact2=-step/deltag
      call blkcpi (1,a,v,k54+nbfgs,neq)
      call blkcpy (w,a(k13),neq)
      do 20 i=1,neq
      v(i)=fact1*v(i)-w(i)
   20 w(i)=fact2*d(i)
      call blkcpy (a(k13),d,neq)
      vw4=abs(4.*fact2*(fact1*g0-g)+4.0)
      vvf=dotprd(v,v)*((dnorm/deltag)**2)
      if (((sqrt(vvf)+sqrt(vvf+vw4))**2)/vw4.lt.1.e+05) go to 40
   30 nupd=0
      call blkcpi (0,d,a,k54+nbfgs,neq)
      go to 60
   40 call blkcpi (0,w,a,k54+nbfgs,3*neq2)
      dw=fact2*g
      do 50 i=1,neq
   50 d(i)=d(i)+dw*v(i)
   60 if (numupd.eq.0) go to 80
      do 70 j=1,numupd
      locat=k54+2*neq2*(numupd-j)
      call blkcpi (1,a,w,locat,2*neq2)
      dw=dotprd(d,w)
      do 70 i=1,neq
   70 d(i)=d(i)+dw*v(i)
   80 call bsolvr (d,a(ntlen),3,8)
      numupd=numupd+nupd
      if (numupd.eq.0) go to 100
      nbfgs=2*neq2*numupd
      do 90 j=1,numupd
      locat=k54+2*neq2*(j-1)
      call blkcpi (1,a,w,locat,2*neq2)
      dv=dotprd(d,v)
      do 90 i=1,neq
   90 d(i)=d(i)+dv*w(i)
  100 return
      end
