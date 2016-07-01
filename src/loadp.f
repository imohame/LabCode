      subroutine loadp(y,z,yn,zn,pmult,nodes,strt,lc,nit)
c     implicit double precision (a-h,o-z)                                    dp
c
c     generates and stores pressure loads
c
      common/bk14/lfna(15),lfnt(6)
      common/bk30/numlp,numpc,h22(2,2),pl2(2,2),h33(3,2),pl3(3,2)
      common/cn0/iconv
      dimension y(*),z(*),yn(3,*),zn(3,*),pmult(4,*),nodes(3,*),strt(*),
     1          lc(1)
      character*80 txts,mssg
      logical nit
c
      if (numpc.eq.0) return
      call basisp
      i=0
      m=0
      inc0=0
      nds3=0
      pml3=0.
      mprnt=0
   10 inc=inc0
      mm=m
      if (.not.nit) then
      iconsv=iconv
      iconv=0
      call gttxsg (txts,lcount)
      read (unit=txts,fmt=110,err=102)
     1 m,lcm,nds1,nds2,nds3,pml1,pml2,pml3,strtm,inc0,ishear
      if(iconsv.eq.1)then
ck      write(lfnt(4),151)
ck     1 m,lcm,nds1,nds2,pml1,pml2,strtm,inc0,ishear,pml3,pml4
      endif
      iconv=iconsv
      endif
      if (nit) then
      call gttxsg (txts,lcount)
      read (unit=txts,fmt=150,err=102)
     1 m,lcm,nds1,nds2,pml1,pml2,strtm,inc0,ishear,pml3,pml4
      endif
      if (m.eq.0) m=i+1
      lc(m)=lcm
      if (pml3.eq.0.) pml3=1.e+20
      if (pml4.eq.0.) pml4=1.e+20
      nodes(1,m)=nds1
      nodes(2,m)=nds2
      nodes(3,m)=nds3
      pmult(1,m)=pml1
      pmult(2,m)=pml2
      pmult(3,m)=pml3
      pmult(4,m)=pml4
      strt(m)=strtm
      if (ishear.eq.1) lc(m)=-lc(m)
      if (inc0.eq.0) inc0=1
      if (pmult(1,m).eq.0.0) pmult(1,m)=1.0
      if (pmult(2,m).eq.0.0) pmult(2,m)=1.0
      if (pmult(3,m).eq.0.0) pmult(3,m)=1.0
   20 i=i+1
      if (m-i) 80,50,30
   30 do 40 j=1,2
      nodes(j,i)=nodes(j,i-1)+inc
   40 pmult(j,i)=pmult(j,mm)
      pmult(3,i)=pmult(3,mm)
      pmult(4,i)=pmult(4,mm)
      strt(i)=strt(mm)
      lc(i)=lc(mm)
      nodes(3,i)=0
      if (nodes(3,i-1).ne.0) nodes(3,i)=nodes(3,i-1)+inc
      ishear=0
   50 if (mprnt.gt.0) go to 60
      mprnt=50
      call header
      write(lfnt(2),120)
   60 mprnt=mprnt-1
      lcc=abs(lc(i))
      ishr=0
      if (lc(i).lt.0) ishr=1
      write(lfnt(2),130) i,lcc,(nodes(j,i),j=1,2),
     1 (pmult(j,i),j=1,4),strt(i),ishr
      if (m-i) 80,70,20
   70 if (numpc-i) 90,90,10
   80 write(lfnt(2),140)
      call bye (2)
   90 continue
c
      do 100 i=1,numpc
      do 100 j=1,3
      k=nodes(j,i)
      if (k.eq.0) go to 100
      yn(j,i)=y(k)
      zn(j,i)=z(k)
  100 continue
c
      return
c
  102 i=i+1
      write (unit=mssg,fmt=160) i
      call termin (txts,mssg,lcount,1)
c
  110 format(5i5,4e10.1,2i5)
  120 format(///' p r e s s u r e  c a r d s'//' card no.',2x,'loadc',
     1 2x,'node1',2x,'node2',6x,'pmult1',6x,'pmult2',6x,
     2 'y-max ',6x,'z-max ',2x,'start time',2x,'shear flag')
  130 format(2i7,2x,2i7,2x,5e12.3,i10)
  140 format(/'fatal error on pressure b.c. cards ')
  150 format(4i5,3e10.1,2i5,2e10.1)
  151 format(4i5,3e10.3,2i5,2e10.3)
  160 format(' error reading pressure cards, probably card #',i4)
      end
