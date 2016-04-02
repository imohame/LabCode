      subroutine ldmas(r,u,ui,udt,udtt,id,inod,xmass,yz,fval,
     1        jnod,xdamp,iter,roller,lm,sf,iprec)
c     implicit double precision (a-h,o-z)                                    dp
c
c     concentrated nodal dampers, masses, and roller boundary conditions
c
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk03/numdc,imassn,idampn,irller,penstf
      common/bk08/kprint,nstep,ite,ilimit,newstf
      common/bk12/ntlen
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk26/dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      common/bk27/nlcur,nptst,nthpy,nthpz,nthps,xmy,xmz,xms,nload,nptm
      common/fissn1/melemt,nnns,ntpe2,n2g,llls
      common/main_block/ a(1)
c
      dimension r(*),u(*),ui(*),udt(*),udtt(*),id(2,*),inod(*),
     1        xmass(2,*),yz(2,*),fval(*),jnod(*),xdamp(2,*),roller(*),
     2        lm(5*iprec,*),sf(5,*)
c
      equivalence (nnns,nd),(llls,nms)
c
      nnns=2
      llls=5
c
      if (imassn.eq.0) go to 110
c
      lcount=0
      do 100 i=1,imassn
      node=inod(i)
      sf(3,i)=0.
      sf(4,i)=0.
      sf(5,i)=0.
      lm(1,i)=id(1,node)
      lm(2,i)=id(2,node)
      ly=lm(1,i)
      lz=lm(2,i)
      xmassy=xmass(1,i)
      xmassz=xmass(2,i)
c
      if (imass.ne.1) go to 40
c
      if (ly.eq.0) go to 30
      wv=-a2*udt(ly)-a3*udtt(ly)
      if (iter.eq.1) wv=wv+a0*ui(ly)
      r(ly)=r(ly)-wv*xmassy
      sf(3,i)=sf(3,i)+a0*xmassy
   30 if (lz.eq.0) go to 40
      wv=-a2*udt(lz)-a3*udtt(lz)
      if (iter.eq.1) wv=wv+a0*ui(lz)
      r(lz)=r(lz)-wv*xmassz
      sf(5,i)=sf(5,i)+a0*xmassz
c
c     body force loads
c
   40 if (nthpy+nthpz+abs(nthps).eq.0) go to 90
c
      if (nthpy.eq.0.or.ly.eq.0) go to 50
      r(ly)=r(ly)-xmassy*fval(nthpy)*xmy
   50 if (nthpz.eq.0.or.lz.eq.0) go to 60
      r(lz)=r(lz)-xmassz*fval(nthpz)*xmz
   60 if (nthps.eq.0) go to 90
      ispny=0
      if (nthps.lt.0) ispny=1
      nthps=abs(nthps)
      facs=(fval(nthps)*xms)**2
      if (ly.eq.0) go to 70
      radius=yz(1,i)
      radius=radius+u(ly)
      r(ly)=r(ly)-xmassy*radius*facs
   70 if (ispny.eq.0.or.lz.eq.0) go to 80
      radius=yz(2,i)
      radius=radius+u(lz)
      r(lz)=r(lz)-xmassz*radius*facs
   80 if (ispny.eq.1) nthps=-nthps
   90 lcount=lcount+1
      if (lcount.lt.64) go to 100
      if (iphase.eq.3) go to 100
      if (newstf.ne.0) go to 100
      melemt=64
      lcount=0
      if (ioofc.eq.1) call fissln (-5,a(ntlen),a(ntlen),lm)
      if (ioofc.eq.0) call wdiskf (lm,melemt)
  100 continue
      if (lcount.eq.0) go to 110
      if (iphase.eq.3) go to 110
      if (newstf.ne.0) go to 110
      melemt=lcount
      if (ioofc.eq.1) call fissln (-5,a(ntlen),a(ntlen),lm)
      if (ioofc.eq.0) call wdiskf (lm,melemt)
c
  110 if (idampn.eq.0) go to 150
c
      lcount=0
      do 140 i=1,idampn
      node=jnod(i)
      sf(3,i)=0.
      sf(4,i)=0.
      sf(5,i)=0.
      lm(1,i)=id(1,node)
      lm(2,i)=id(2,node)
      ly=lm(1,i)
      lz=lm(2,i)
      xdampy=xdamp(1,i)
      xdampz=xdamp(2,i)
      if (ly.eq.0) go to 120
      wv=-a4*udt(ly)-a5*udtt(ly)
      if (iter.eq.1) wv=wv+a1*ui(ly)
      r(ly)=r(ly)-wv*xdampy
      sf(3,i)=sf(3,i)+a1*xdampy
  120 if (lz.eq.0) go to 130
      wv=-a4*udt(lz)-a5*udtt(lz)
      if (iter.eq.1) wv=wv+a1*ui(lz)
      r(lz)=r(lz)-wv*xdampz
      sf(5,i)=sf(5,i)+a1*xdampz
  130 lcount=lcount+1
      if (lcount.lt.64) go to 140
      if (iphase.eq.3) go to 140
      if (newstf.ne.0) go to 140
      melemt=64
      lcount=0
      if (ioofc.eq.1) call fissln (-5,a(ntlen),a(ntlen),lm)
      if (ioofc.eq.0) call wdiskf (lm,melemt)
  140 continue
      if (lcount.eq.0) go to 150
      if (iphase.eq.3) go to 150
      if (newstf.ne.0) go to 150
      melemt=lcount
      if (ioofc.eq.1) call fissln (-5,a(ntlen),a(ntlen),lm)
      if (ioofc.eq.0) call wdiskf (lm,melemt)
c
  150 if(irller.eq.0) return
c
      lcount=0
      do 160 i=1,numnp
      if(roller(i).eq.0.) go to 160
      angle=roller(i)/57.29577951
      sinthe=sin(angle)
      costhe=cos(angle)
      lcount=lcount+1
      lm(1,lcount)=id(1,i)
      lm(2,lcount)=id(2,i)
      sf(3,lcount)= penstf*sinthe*sinthe
      sf(4,lcount)=-penstf*sinthe*costhe
      sf(5,lcount)= penstf*costhe*costhe
      rtan=r(id(1,i))*costhe+r(id(2,i))*sinthe
      r(id(1,i))=rtan*costhe
      r(id(2,i))=rtan*sinthe
      if (lcount.lt.64) go to 160
      if (iphase.eq.3) go to 160
      if (newstf.ne.0) go to 160
      melemt=64
      lcount=0
      if (ioofc.eq.1) call fissln (-5,a(ntlen),a(ntlen),lm)
      if (ioofc.eq.0) call wdiskf (lm,melemt)
  160 continue
      if (lcount.eq.0) return
      if (iphase.eq.3) return
      if (newstf.ne.0) return
      melemt=lcount
      if (ioofc.eq.1) call fissln (-5,a(ntlen),a(ntlen),lm)
      if (ioofc.eq.0) call wdiskf (lm,melemt)
c
      return
c
      end
