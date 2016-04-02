      subroutine cnlods(y,z,cr,nod,idirn,ncur,fac)
c     implicit double precision (a-h,o-z)                                    dp
c
      common/bk14/lfna(15),lfnt(6)
      common/bk27/ndummy,nptst,nthpy,nthpz,nthps,xmy,xmz,xms,nload,nptm
c
      dimension y(*),z(*),cr(4,*),nod(*),idirn(*),ncur(*),fac(*)
      character*80 txts,mssg
c
      if (nload.eq.0) return
      call header
      write(lfnt(2),30)
      ifolow=0
      do 10 i=1,nload
      call gttxsg (txts,lcount)
      read(unit=txts,fmt=20,err=18) nod(i),idirn(i),ncur(i),fac(i),ifw
      if (fac(i).eq.0.0) fac(i)=1.0
      write(lfnt(2),40) nod(i),idirn(i),ncur(i),fac(i),ifw
      fact=fac(i)
      if (ifw.eq.0) go to 10
      ifolow=1
      ncur(i)=-ncur(i)
      nod1=nod(i)
      nod2=idirn(i)
      cr(1,i)=y(nod1)
      cr(2,i)=z(nod1)
      cr(3,i)=y(nod2)
      cr(4,i)=z(nod2)
   10 continue
c
      return
c
   18 write (unit=mssg,fmt=50) i
      call termin (txts,mssg,lcount,1)
c
   20 format(3i5,f10.0,i5)
   30 format(///' c o n c e n t r a t e d   l o a d s'/
     1/4x,' node1     node2/    load curve     load curve     follower'
     2/4x,'         direction                   multiplier       flag')
   40 format(3x,i6,4x,i7,9x,i5,8x,e12.4,i6)
   50 format(' error reading concentrated nodal load card #',i5)
      end
