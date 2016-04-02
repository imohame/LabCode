      subroutine ulbc2(iequit,lm,s,fval,idir,lc,xmag,iretn,ix,id,
     1 nodes,iprec)
c     implicit double precision (a-h,o-z)                                    dp
c
c     impose velocity boundary conditions at element level (only for plane strain/stress, doesnt work for axisymmetric)
c
       use mod_parameters
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk03/numdc,imassn,idampn,irller,penstf
      common/bk27/nlcur,nptst,nthpy,nthpz,nthps,xmy,xmz,xms,nload,nptm
      common/bk32/nsref,nequit,time,timep,lprint,nprint
      common/range/mft,mlt,lft,llt,nftm1
      common/vect1/ r(nelemg,10)
      common/vect5/
     1 wv1(nelemg),wv2(nelemg),wv3(nelemg),wv4(nelemg),
     2 wv5(nelemg),wv6(nelemg),wv7(nelemg),wv8(nelemg)
      common/vect9/scl(8,nelemg),yz(2,4,nelemg)
      common/bk26/dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      dimension fval(*),idir(*),lc(*),xmag(8,1),iretn(*),
     1 lm(44*iprec,*),s(44,1),ix(4,*),id(2,*),dsbr(2),nodes(1)
      if (numdc.eq.0) return
      do 190 l=lft,llt
      if (iretn(l).eq.1) go to 190
      iretn(l)=1
      do 170 i=1,numdc
      if (time.ge.xmag(2,i)) go to 170
      ndf=idir(i)
      if (ndf.ne.3000000) then
      dsb=0.
      if (iequit.eq.1) go to 10
      lcv=lc(i)
      dsb=(fval(lcv))*dt*xmag(1,i)
c	  write(*,*) dsb
   10 do 20 j=1,8
      k=j
      if (lm(j,l).eq.ndf) go to 30
   20 continue
      go to 170
   30 n=1
      iretn(l)=0
      do 40 m=1,k
      ml=(k*(k-1))/2+m
      r(l,n)=r(l,n)+s(ml,l)*dsb
      s(ml,l)=0.0
   40 n=n+1
      if (k.eq.8) go to 60
      k1=k+1
      do 50 m=k1,8
      kl=(m*(m-1))/2+k
      r(l,n)=r(l,n)+s(kl,l)*dsb
      s(kl,l)=0.0
   50 n=n+1
   60 s(ml,l)=1.0
      r(l,k)=-dsb
      else
      radian=45./atan(1.)
      dsbr(1)=0.0
      dsbr(2)=0.0
      nodnum=nodes(i)
      if     (nodnum.eq.ix(1,l)) then
      isave=1
      elseif (nodnum.eq.ix(2,l)) then
      isave=2
      elseif (nodnum.eq.ix(3,l)) then
      isave=3
      elseif (nodnum.eq.ix(4,l)) then
      isave=4
      else
      go to 170
      endif
      if (iequit.eq.1) go to 100
      lcv=lc(i)
      fvlnew=fval(lcv)*xmag(1,i)/radian
      fvlold=fval(lcv+nlcur)*xmag(1,i)/radian
      rinc=yz(1,isave,l)-xmag(3,i)
      zinc=yz(2,isave,l)-xmag(4,i)
      angle0=atan2(zinc,rinc)
      angln0=angle0+fvlold
      angln1=angle0+fvlnew
      radius=sqrt(rinc**2+zinc**2)
      dsbr(1)=radius*(cos(angln1)-cos(angln0))
      dsbr(2)=radius*(sin(angln1)-sin(angln0))
  100 do 165 itwice=1,2
      ndf=id(itwice,nodnum)
      dsb=dsbr(itwice)
      do 120 j=1,8
      k=j
      if (lm(j,l).eq.ndf) go to 130
  120 continue
      go to 170
  130 n=1
      iretn(l)=0
      do 140 m=1,k
      ml=(k*(k-1))/2+m
      r(l,n)=r(l,n)+s(ml,l)*dsb
      s(ml,l)=0.0
  140 n=n+1
      if (k.eq.8) go to 160
      k1=k+1
      do 150 m=k1,8
      kl=(m*(m-1))/2+k
      r(l,n)=r(l,n)+s(kl,l)*dsb
      s(kl,l)=0.0
  150 n=n+1
  160 s(ml,l)=1.0
      r(l,k)=-dsb
  165 continue
      endif
  170 continue
  190 continue
      return
      end
