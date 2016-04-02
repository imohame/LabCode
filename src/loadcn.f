      subroutine loadcn(id,u,r,npc,p,cr,nod,idirn,ncur,fac,tt)
c     implicit double precision (a-h,o-z)                                    dp
c
c     apply concentrated loads and follower forces
c
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk27/nlcur,nptst,nthpy,nthpz,nthps,xmy,xmz,xms,nload,nptm
      common/bks17/fmult(5)
      dimension id(2,*),u(*),r(*),cr(4,*),nod(*),idirn(*),ncur(*),
     1          fac(*),npc(*),p(*)
c
      do 10 n=1,neq
   10 r(n)=0.0
c
      if (nload.eq.0) return
c
      do 40 n=1,nload
      ierr=0
      xmag=fac(n)
      iflwer=0
      lcc=ncur(n)
      if (lcc.gt.0) go to 20
      iflwer=1
      lcc=-lcc
   20 loc=npc(lcc)
      npoint=(npc(lcc+1)-loc)/2
      f=0.0
c.... island (loadset command)
      if(lcc.le.5)xmag=xmag*fmult(lcc)
      call interp (p(loc),tt,npoint,f,xmag,ierr)
      nod1=nod(n)
      nod2=idirn(n)
      if (iflwer.eq.1) go to 30
      idof=id(nod2,nod1)
      if (idof.eq.0) go to 40
      r(idof)=r(idof)+f
      go to 40
c
   30 id1y=id(1,nod1)
      id1z=id(2,nod1)
      id2y=id(1,nod2)
      id2z=id(2,nod2)
      y1=cr(1,n)
      z1=cr(2,n)
      y2=cr(3,n)
      z2=cr(4,n)
      if (id1y.ne.0) y1=y1+u(id1y)
      if (id1z.ne.0) z1=z1+u(id1z)
      if (id2y.ne.0) y2=y2+u(id2y)
      if (id2z.ne.0) z2=z2+u(id2z)
      y12=y1-y2
      z12=z1-z2
      xl=sqrt(y12**2+z12**2)
      if (id1y.ne.0) r(id1y)=r(id1y)+f*z12/xl
      if (id1z.ne.0) r(id1z)=r(id1z)-f*y12/xl
   40 continue
c
      return
      end
