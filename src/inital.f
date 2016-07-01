      subroutine inital(y,z,u,udt,udtt,id,v,icon)
c     implicit double precision (a-h,o-z)                                    dp
c
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk14/lfna(15),lfnt(6)
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      dimension y(*),z(*),u(*),udt(*),udtt(*),id(2,*),v(*)
      character*80 txts,mssg
c
      do 10 i=1,neq
   10 u(i)=0.
      if (imass.eq.0) go to 180
      do 20 i=1,neq
      udt(i)=0.
   20 udtt(i)=0.
      if (icon.eq.0) go to 180
      if (icon.lt.0) go to 200
      iflag=0
      if (icon.eq.2) iflag=1
   30 i=0
      ie=1
   40 call gttxsg (txts,lcount)
      read(unit=txts,fmt=220,err=214) n,v(2*n-1),v(2*n),nd
      if (nd) 50,60,50
   50 ie=nd
   60 if (i) 70,130,70
   70 nl=n-i
      if (nl-1) 130,120,80
   80 nl=nl/ie
      if (i+nl*ie-n) 190,90,190
   90 if (nl-1) 120,120,100
  100 anl=nl
      dr=(v(2*n-1)-v(2*i-1))/anl
      dz=(v(2*n)-v(2*i))/anl
      nl=n-2*ie
      do 110 j=i,nl,ie
      i1=j+ie
      v(2*i1-1)=v(2*j-1)+dr
  110 v(2*i1)=v(2*j)+dz
  120 if (numnp-n) 190,140,130
  130 i=n
      go to 40
  140 do 170 i=1,numnp
      ii=id(1,i)
      jj=id(2,i)
      if (iflag) 150,150,160
  150 if (ii.ne.0) udt(ii)=v(2*i-1)
      if (jj.ne.0) udt(jj)=v(2*i)
      go to 170
  160 if (ii.ne.0) udtt(ii)=v(2*i-1)
      if (jj.ne.0) udtt(jj)=v(2*i)
  170 continue
      iflag=iflag+1
      if (iflag.eq.1.and.icon.eq.3) go to 30
c
  180 return
c
  190 if (iflag.eq.0) write(lfnt(2),230)
      if (iflag.eq.1) write(lfnt(2),240)
      call bye (2)
c
c     angular velocity
c
  200 mssg=' error reading angular velocity card for velocity initializa
     1tion'
      call gttxsg (txts,lcount)
      read(unit=txts,fmt=250,err=216) omega,ynod,znod
      do 210 i=1,numnp
      idofy=id(1,i)
      idofz=id(2,i)
      dy=y(i)-ynod
      dz=z(i)-znod
      if (idofy.ne.0) udt(idofy)=-dz*omega
      if (idofz.ne.0) udt(idofz)=dy*omega
  210 continue
c
      return
c
  214 i=i+1
      write (unit=mssg,fmt=260) i
  216 call termin (txts,mssg,lcount,1)
c
  220 format(i5,2e10.1,i5)
  230 format(/' fatal error on initial velocity cards ')
  240 format(/' fatal error on initial acceleration cards ')
  250 format(3e10.1)
  260 format (' error read initial velocity cards, probably card #',i5)
      end
