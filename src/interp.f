      subroutine interp(p,tau,numlp,f,xmag,ierr)
c     implicit double precision (a-h,o-z)                                    dp
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk26/delt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      dimension p(2,*)
      if (mthsol.lt.6.and.delt.ge.0.0) then
      if (tau-p(1,numlp).gt.1.0e-07) go to 30
      do 10 m=2,numlp
      n=m
      if (tau-p(1,m).le.0.0) go to 20
      if (abs(tau-p(1,m)).le.1.0e-07) go to 20
   10 continue
      go to 30
   20 dt=tau-p(1,n-1)
      d1=p(1,n)-p(1,n-1)
      d2=p(2,n)-p(2,n-1)
      f=p(2,n-1)+dt*d2/d1
      f=xmag*f
      return
   30 ierr=1
      f=0.0
      return
      else
      dt=tau-p(1,1)
      d1=p(1,2)-p(1,1)
      d2=p(2,2)-p(2,1)
      f=p(2,1)+dt*d2/d1
      f=xmag*f
      return
      endif
      end
